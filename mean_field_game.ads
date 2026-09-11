--  Mean_Field_Game — Ada 2023 educational package for finite-state,
--  discrete-time mean-field games (MFG). Continuum limit of large symmetric
--  N-player games: a representative agent optimises against the population
--  distribution (mean field) m; consistency requires that the optimal
--  feedback induces a flow of m matching the assumed field.
--  Classroom scope: finite states / actions, congestion-style running costs
--  depending on (s,a,m), action-controlled transitions, best-response value
--  iteration vs fixed m, occupancy (Fokker–Planck) update, damped fixed-point
--  iteration. Not the continuous Lasry–Lions PDE system.
--  Reference: https://en.wikipedia.org/wiki/Mean-field_game_theory
--  Sibling sheets (README only — do not `with`): Potential Game, N-player
--  Nash — RobertBoettcherSF Ada algorithm series.

pragma Ada_2022;

package Mean_Field_Game
  with SPARK_Mode => Off
is

   ---------------------------------------------------------------------------
   -- Capacity (educational)
   ---------------------------------------------------------------------------

   Max_States  : constant Positive := 8;
   Max_Actions : constant Positive := 4;

   ---------------------------------------------------------------------------
   -- Identifiers and numeric types
   ---------------------------------------------------------------------------

   type State_Id is range 1 .. Max_States;
   type Action_Id is range 1 .. Max_Actions;

   --  Share 'Base so slices 1 .. S / 1 .. A type-check with discriminants.
   subtype State_Count is State_Id'Base range 0 .. State_Id'Base (Max_States);
   subtype Action_Count is Action_Id'Base range 0 .. Action_Id'Base (Max_Actions);

   type Real is digits 15;
   subtype Probability is Real;
   subtype Cost is Real;

   --  Population distribution m on states (probability simplex).
   type Distribution is array (State_Id range <>) of Probability;

   --  Deterministic Markov policy: π(s) ∈ Actions.
   type Policy is array (State_Id range <>) of Action_Id;

   --  Soft (mixed) policy: Mix(s,a) = Prob[action a | state s].
   type Soft_Policy is
     array (State_Id range <>, Action_Id range <>) of Probability;

   type Value_Vector is array (State_Id range <>) of Cost;

   --  Trans(s, a, s').
   type Transition_Tensor is
     array (State_Id range <>, Action_Id range <>, State_Id range <>)
       of Probability;

   --  Base running cost c0(s,a) before mean-field coupling.
   type Cost_Matrix is
     array (State_Id range <>, Action_Id range <>) of Cost;

   --  Per-state congestion weight: running cost adds Congestion(s) * m(s).
   type Congestion_Weights is array (State_Id range <>) of Cost;

   --  Markov kernel K(s, s' ) = P(s' | s) under a fixed policy.
   type Kernel_Matrix is
     array (State_Id range <>, State_Id range <>) of Probability;

   ---------------------------------------------------------------------------
   -- Finite-state discrete-time MFG
   ---------------------------------------------------------------------------

   --  Discounted infinite-horizon MFG on a finite state/action space.
   --  Running cost of representative agent in state s taking action a given
   --  mean field m:
   --    c(s,a,m) = Base_Cost(s,a) + Congestion(s) * m(s)
   --  Transitions depend on a only (not on m) — standard classroom toy.
   type MFG
     (N_States  : State_Count;
      N_Actions : Action_Count) is
   record
      Trans      : Transition_Tensor
                     (1 .. N_States, 1 .. N_Actions, 1 .. N_States);
      Base_Cost  : Cost_Matrix (1 .. N_States, 1 .. N_Actions);
      Congestion : Congestion_Weights (1 .. N_States);
      Discount   : Real := 0.9;  -- γ ∈ [0, 1)
   end record;

   --  Result of fixed-point / Solve iteration.
   type Solve_Result
     (N_States  : State_Count;
      N_Actions : Action_Count) is
   record
      Mean_Field : Distribution (1 .. N_States);
      Pi         : Policy (1 .. N_States);
      Values     : Value_Vector (1 .. N_States);
      Residual   : Real := 0.0;
      Iterations : Natural := 0;
      Converged  : Boolean := False;
   end record;

   ---------------------------------------------------------------------------
   -- Exceptions
   ---------------------------------------------------------------------------

   Invalid_Argument : exception;
   --  Raised for zero / oversized dimensions, non-simplex distributions,
   --  non-stochastic transition rows, discount ∉ [0,1), negative tolerances,
   --  damping ∉ (0,1], or indices out of range.

   ---------------------------------------------------------------------------
   -- Tolerances / Near
   ---------------------------------------------------------------------------

   Default_Tol : constant Real := 1.0E-9;

   function Near
     (X, Y : Real; Tol : Real := Default_Tol) return Boolean
     with Global => null;
   --  |X − Y| ≤ Tol. Tol must be ≥ 0 (else Invalid_Argument).

   ---------------------------------------------------------------------------
   -- Distribution / simplex helpers
   ---------------------------------------------------------------------------

   function Mass (M : Distribution) return Probability
     with Global => null;
   --  Σ_s M(s).

   function Is_Simplex
     (M : Distribution; Tol : Real := Default_Tol) return Boolean
     with Global => null;
   --  Nonempty, every entry ≥ −Tol, |Mass − 1| ≤ Tol. Tol ≥ 0.

   function Normalize (M : Distribution) return Distribution
     with Global => null;
   --  Scale so Mass = 1. Raises Invalid_Argument if Mass ≤ 0.

   function Uniform_Distribution (N : State_Count) return Distribution
     with Global => null;
   --  1/N on each state. Raises if N = 0.

   function Dirac_Distribution
     (N : State_Count; S : State_Id) return Distribution
     with Global => null;
   --  Unit mass on S. Raises if N = 0 or S out of 1 .. N.

   function L1_Distance (A, B : Distribution) return Real
     with Global => null;
   --  Σ_s |A(s) − B(s)|. Raises on length mismatch.

   function Total_Variation (A, B : Distribution) return Real
     with Global => null;
   --  (1/2) L1_Distance. Raises on length mismatch.

   ---------------------------------------------------------------------------
   -- Transition / cost helpers
   ---------------------------------------------------------------------------

   function Is_Stochastic_Row
     (G   : MFG;
      S   : State_Id;
      A   : Action_Id;
      Tol : Real := Default_Tol) return Boolean
     with Global => null;
   --  Trans(S,A,·) nonnegative within Tol and sums to 1 within Tol.

   function Is_Valid_Transitions
     (G : MFG; Tol : Real := Default_Tol) return Boolean
     with Global => null;
   --  Every (s,a) row is stochastic. Tol ≥ 0. N_States, N_Actions > 0.

   function Running_Cost
     (G : MFG; S : State_Id; A : Action_Id; M : Distribution) return Cost
     with Global => null;
   --  c(s,a,m) = Base_Cost(s,a) + Congestion(s)*M(s). Raises on bad dims /
   --  indices or if M'Length /= G.N_States.

   function Require_Compatible_Distribution
     (G : MFG; M : Distribution) return Boolean
     with Global => null;
   --  True when M'First = 1 and M'Last = G.N_States; else raises.

   ---------------------------------------------------------------------------
   -- Policies
   ---------------------------------------------------------------------------

   function Is_Valid_Policy (G : MFG; Pi : Policy) return Boolean
     with Global => null;
   --  Pi covers 1 .. N_States and each Pi(s) ∈ 1 .. N_Actions.

   function Pure_Soft_Policy
     (G : MFG; Pi : Policy) return Soft_Policy
     with Global => null;
   --  Dirac soft policy from deterministic Pi. Raises if Pi invalid.

   function Uniform_Soft_Policy (G : MFG) return Soft_Policy
     with Global => null;
   --  1/N_Actions in every state. Raises if N_Actions = 0.

   ---------------------------------------------------------------------------
   -- Best response vs fixed mean field (Bellman / value iteration)
   ---------------------------------------------------------------------------

   procedure Value_Iteration
     (G         : MFG;
      M         : Distribution;
      Values    : out Value_Vector;
      Pi        : out Policy;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000;
      Iters     : out Natural)
     with Global => null;
   --  Discounted value iteration for the representative agent's MDP with
   --  frozen mean field M. On exit Values ≈ Bellman fixed point and Pi is a
   --  greedy (argmin) policy. Raises on incompatible M, invalid G (discount
   --  ∉ [0,1), empty dims), or Tol < 0.

   function Best_Response_Policy
     (G         : MFG;
      M         : Distribution;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Policy
     with Global => null;
   --  Deterministic greedy policy from Value_Iteration vs M.

   function Best_Response_Values
     (G         : MFG;
      M         : Distribution;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Value_Vector
     with Global => null;
   --  Value vector from Value_Iteration vs M.

   function Soft_Best_Response
     (G         : MFG;
      M         : Distribution;
      Beta      : Real := 5.0;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Soft_Policy
     with Global => null;
   --  Logit / Boltzmann soft BR: Mix(s,a) ∝ exp(−Beta * Q(s,a)) with
   --  Q from Best_Response_Values. Beta > 0 (else Invalid_Argument).
   --  Useful when pure-policy consistency cannot hold exactly.

   function Q_Value
     (G : MFG;
      M : Distribution;
      V : Value_Vector;
      S : State_Id;
      A : Action_Id) return Cost
     with Global => null;
   --  Q(s,a) = c(s,a,m) + γ Σ_{s'} P(s'|s,a) V(s').

   function Policy_Evaluation
     (G         : MFG;
      M         : Distribution;
      Pi        : Policy;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Value_Vector
     with Global => null;
   --  Evaluate fixed Pi against frozen M (linear Bellman iteration).

   ---------------------------------------------------------------------------
   -- Occupancy / Fokker–Planck update
   ---------------------------------------------------------------------------

   function Evolve_Occupancy
     (G  : MFG;
      M  : Distribution;
      Pi : Policy) return Distribution
     with Global => null;
   --  One-step mean-field update:
   --    m'(s') = Σ_s m(s) P(s' | s, π(s)).
   --  Raises if M / Pi incompatible with G.

   function Evolve_Occupancy_Soft
     (G   : MFG;
      M   : Distribution;
      Mix : Soft_Policy) return Distribution
     with Global => null;
   --  m'(s') = Σ_{s,a} m(s) Mix(s,a) P(s'|s,a).

   function Evolve_Occupancy_N
     (G     : MFG;
      M     : Distribution;
      Pi    : Policy;
      Steps : Natural) return Distribution
     with Global => null;
   --  Apply Evolve_Occupancy Steps times (Steps = 0 returns a copy of M).

   function Induced_Kernel
     (G  : MFG;
      Pi : Policy) return Kernel_Matrix
     with Global => null;
   --  K(s, s') = P(s' | s, π(s)). Prefer Evolve_Occupancy for classroom use.

   ---------------------------------------------------------------------------
   -- Fixed-point iteration (consistency of m and π)
   ---------------------------------------------------------------------------

   procedure Fixed_Point_Iterate
     (G          : MFG;
      M0         : Distribution;
      Result     : out Solve_Result;
      Damping    : Real := 0.5;
      Tol        : Real := 1.0E-8;
      Max_Iters  : Natural := 500;
      Occupancy_Steps : Natural := 1)
     with Global => null;
   --  Alternate best-response policy vs current m and damped occupancy
   --  update until L1 residual ≤ Tol or Max_Iters reached.
   --    m ← (1 − α) m + α Evolve^k(m, π*(m))
   --  with α = Damping ∈ (0,1], k = Occupancy_Steps ≥ 1.
   --  Raises on bad args / incompatible M0.

   function Solve
     (G          : MFG;
      M0         : Distribution;
      Damping    : Real := 0.5;
      Tol        : Real := 1.0E-8;
      Max_Iters  : Natural := 500;
      Occupancy_Steps : Natural := 1) return Solve_Result
     with Global => null;
   --  Functional wrapper around Fixed_Point_Iterate.

   function Consistency_Residual
     (G  : MFG;
      M  : Distribution;
      Pi : Policy) return Real
     with Global => null;
   --  L1_Distance(M, Evolve_Occupancy(G, M, Pi)).

   ---------------------------------------------------------------------------
   -- Classic classroom constructors
   ---------------------------------------------------------------------------

   function Congestion_Toy return MFG;
   --  3 states (locations), 2 actions (Stay / Move). Moving jumps to a
   --  neighbour; staying keeps state with small leak. Congestion cost rises
   --  with local occupancy — agents prefer sparse locations.

   function Flocking_Toy return MFG;
   --  3 states on a line. Actions Left / Stay / Right. Mild preference to
   --  move toward the population mean (attraction / flocking) encoded as
   --  negative congestion on popular states plus move costs.

   function Entry_Toy return MFG;
   --  2 states: Out / In. Actions Stay / Switch. Being In is valuable unless
   --  the market is crowded (congestion on In); switching has a fixed cost.
   --  Classic sparse/entry MFG flavour.

   function Two_State_Symmetric return MFG;
   --  Minimal 2-state / 2-action symmetric congestion game for unit tests.

   ---------------------------------------------------------------------------
   -- Validation of a whole MFG record
   ---------------------------------------------------------------------------

   function Is_Valid_MFG
     (G : MFG; Tol : Real := Default_Tol) return Boolean
     with Global => null;
   --  N_States, N_Actions > 0, Discount ∈ [0,1), transitions stochastic,
   --  Tol ≥ 0.

end Mean_Field_Game;
