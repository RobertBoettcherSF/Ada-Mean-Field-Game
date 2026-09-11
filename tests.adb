--  Standalone test suite for Mean_Field_Game.

pragma Ada_2022;

with Ada.Command_Line;
with Ada.Text_IO; use Ada.Text_IO;
with Mean_Field_Game; use Mean_Field_Game;

procedure Tests is

   Pass_Count : Natural := 0;
   Fail_Count : Natural := 0;

   procedure Check (Condition : Boolean; Message : String) is
   begin
      if Condition then
         Pass_Count := Pass_Count + 1;
         Put_Line ("  PASS: " & Message);
      else
         Fail_Count := Fail_Count + 1;
         Put_Line ("  FAIL: " & Message);
      end if;
   end Check;

   procedure Section (Title : String) is
   begin
      New_Line;
      Put_Line ("=== " & Title & " ===");
   end Section;

   --  Non-static views (avoid -gnatwa constant-condition warnings).
   function Nat (X : Natural) return Natural is (X);
   function SC (X : Natural) return State_Count is (State_Count (X));
   function AC (X : Natural) return Action_Count is (Action_Count (X));
   function Sid (X : Positive) return State_Id is (State_Id (X));
   function Aid (X : Positive) return Action_Id is (Action_Id (X));
   function Rf (X : Real) return Real is (X);

   ---------------------------------------------------------------------------
   -- Exception helpers
   ---------------------------------------------------------------------------

   function Near_Raises (Tol : Real) return Boolean is
      Unused : Boolean;
   begin
      Unused := Near (0.0, 0.0, Tol);
      return False;
   exception
      when Invalid_Argument =>
         return True;
   end Near_Raises;

   function Normalize_Raises (M : Distribution) return Boolean is
   begin
      declare
         R : constant Distribution := Normalize (M);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Normalize_Raises;

   function Uniform_Raises (N : State_Count) return Boolean is
   begin
      declare
         R : constant Distribution := Uniform_Distribution (N);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Uniform_Raises;

   function Dirac_Raises (N : State_Count; S : State_Id) return Boolean is
   begin
      declare
         R : constant Distribution := Dirac_Distribution (N, S);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Dirac_Raises;

   function L1_Raises (A, B : Distribution) return Boolean is
   begin
      declare
         D : constant Real := L1_Distance (A, B);
         pragma Unreferenced (D);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end L1_Raises;

   function Running_Cost_Raises
     (G : MFG; S : State_Id; A : Action_Id; M : Distribution) return Boolean
   is
   begin
      declare
         C : constant Cost := Running_Cost (G, S, A, M);
         pragma Unreferenced (C);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Running_Cost_Raises;

   function BR_Raises (G : MFG; M : Distribution) return Boolean is
   begin
      declare
         P : constant Policy := Best_Response_Policy (G, M);
         pragma Unreferenced (P);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end BR_Raises;

   function Evolve_Raises
     (G : MFG; M : Distribution; Pi : Policy) return Boolean
   is
   begin
      declare
         R : constant Distribution := Evolve_Occupancy (G, M, Pi);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Evolve_Raises;

   function Solve_Damping_Raises
     (G : MFG; M : Distribution; Damp : Real) return Boolean
   is
   begin
      declare
         R : constant Solve_Result :=
           Solve (G, M, Damping => Damp, Max_Iters => 2);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Solve_Damping_Raises;

   function Solve_Occ_Raises
     (G : MFG; M : Distribution) return Boolean
   is
   begin
      declare
         R : constant Solve_Result :=
           Solve (G, M, Occupancy_Steps => 0, Max_Iters => 2);
         pragma Unreferenced (R);
      begin
         return False;
      end;
   exception
      when Invalid_Argument =>
         return True;
   end Solve_Occ_Raises;

   ---------------------------------------------------------------------------
   -- Near / Mass / Simplex
   ---------------------------------------------------------------------------

   procedure Test_Near_And_Simplex is
      U : constant Distribution := Uniform_Distribution (SC (3));
      D : constant Distribution := Dirac_Distribution (SC (3), Sid (2));
      Z : constant Distribution (1 .. 3) := [0.0, 0.0, 0.0];
      B : constant Distribution (1 .. 3) := [0.5, 0.5, 0.2];
      Neg : constant Distribution (1 .. 2) := [-0.1, 1.1];
   begin
      Section ("Near / Mass / Simplex");
      Check (Near (1.0, 1.0), "Near equal");
      Check (Near (1.0, 1.0 + 1.0E-12), "Near within default tol");
      Check (not Near (1.0, 2.0), "Near far");
      Check (Near_Raises (Rf (-1.0)), "Near negative tol raises");
      Check (Near (Mass (U), 1.0), "uniform mass 1");
      Check (Near (Mass (D), 1.0), "dirac mass 1");
      Check (Is_Simplex (U), "uniform is simplex");
      Check (Is_Simplex (D), "dirac is simplex");
      Check (not Is_Simplex (Z), "zero vector not simplex");
      Check (not Is_Simplex (B), "mass 1.2 not simplex");
      Check (not Is_Simplex (Neg), "negative entry not simplex");
      Check (Near (U (Sid (1)), 1.0 / 3.0), "uniform entry");
      Check (Near (D (Sid (2)), 1.0), "dirac peak");
      Check (Near (D (Sid (1)), 0.0), "dirac off-peak");
      Check (Normalize_Raises (Z), "normalize zero raises");
      Check (Uniform_Raises (SC (0)), "uniform N=0 raises");
      Check (Dirac_Raises (SC (0), Sid (1)), "dirac N=0 raises");
      Check (Dirac_Raises (SC (3), Sid (4)), "dirac bad state raises");
      declare
         N : constant Distribution := Normalize (B);
      begin
         Check (Is_Simplex (N), "normalize mass 1.2");
         Check (Near (Mass (N), 1.0), "normalized mass");
      end;
      Check (Near (L1_Distance (U, U), 0.0), "L1 self 0");
      Check (Near (Total_Variation (U, U), 0.0), "TV self 0");
      Check (L1_Distance (U, D) > 0.0, "L1 uniform vs dirac > 0");
      Check
        (Near (Total_Variation (U, D), 0.5 * L1_Distance (U, D)),
         "TV = half L1");
      declare
         Short : constant Distribution (1 .. 2) := [0.5, 0.5];
      begin
         Check (L1_Raises (U, Short), "L1 length mismatch raises");
      end;
      Check (Nat (Max_States) = 8, "Max_States = 8");
      Check (Nat (Max_Actions) = 4, "Max_Actions = 4");
   end Test_Near_And_Simplex;

   ---------------------------------------------------------------------------
   -- Constructors / validity
   ---------------------------------------------------------------------------

   procedure Test_Constructors is
      C : constant MFG := Congestion_Toy;
      F : constant MFG := Flocking_Toy;
      E : constant MFG := Entry_Toy;
      T : constant MFG := Two_State_Symmetric;
   begin
      Section ("Constructors / Is_Valid_MFG");
      Check (C.N_States = SC (3), "Congestion states");
      Check (C.N_Actions = AC (2), "Congestion actions");
      Check (Is_Valid_MFG (C), "Congestion valid");
      Check (Is_Valid_Transitions (C), "Congestion transitions");
      Check (F.N_States = SC (3), "Flocking states");
      Check (F.N_Actions = AC (3), "Flocking actions");
      Check (Is_Valid_MFG (F), "Flocking valid");
      Check (E.N_States = SC (2), "Entry states");
      Check (E.N_Actions = AC (2), "Entry actions");
      Check (Is_Valid_MFG (E), "Entry valid");
      Check (T.N_States = SC (2), "TwoState states");
      Check (Is_Valid_MFG (T), "TwoState valid");
      Check (C.Discount > 0.0 and then C.Discount < 1.0, "Congestion discount");
      Check (F.Congestion (Sid (1)) < 0.0, "Flocking attraction");
      Check (E.Congestion (Sid (2)) > 0.0, "Entry congestion on In");
      Check (E.Congestion (Sid (1)) = 0.0, "Entry no congestion Out");
      Check
        (Is_Stochastic_Row (C, Sid (1), Aid (1)), "Congestion row (1,1)");
      Check
        (Is_Stochastic_Row (C, Sid (2), Aid (2)), "Congestion row (2,2)");
      Check
        (Is_Stochastic_Row (F, Sid (3), Aid (3)), "Flocking row (3,3)");
      Check
        (Is_Stochastic_Row (E, Sid (1), Aid (2)), "Entry row (1,2)");
      --  Per-state congestion weights
      for S in State_Id range 1 .. 3 loop
         Check (C.Congestion (S) > 0.0, "Congestion weight positive");
      end loop;
      for S in State_Id range 1 .. 3 loop
         Check (F.Congestion (S) < 0.0, "Flocking weight negative");
      end loop;
   end Test_Constructors;

   ---------------------------------------------------------------------------
   -- Running cost / policies
   ---------------------------------------------------------------------------

   procedure Test_Cost_And_Policy is
      G : constant MFG := Two_State_Symmetric;
      M : constant Distribution := Uniform_Distribution (SC (2));
      D1 : constant Distribution := Dirac_Distribution (SC (2), Sid (1));
      Pi : Policy (1 .. 2);
      Mix : Soft_Policy (1 .. 2, 1 .. 2);
      Bad_M : constant Distribution (1 .. 3) := [0.3, 0.3, 0.4];
   begin
      Section ("Running_Cost / Policies");
      Check
        (Near
           (Running_Cost (G, Sid (1), Aid (1), M),
            G.Base_Cost (Sid (1), Aid (1)) + G.Congestion (Sid (1)) * M (1)),
         "running cost formula uniform");
      Check
        (Running_Cost (G, Sid (1), Aid (1), D1) >
           Running_Cost (G, Sid (1), Aid (1), M),
         "more mass raises congestion cost at s=1");
      Check
        (Running_Cost_Raises (G, Sid (3), Aid (1), M),
         "bad state raises");
      Check
        (Running_Cost_Raises (G, Sid (1), Aid (3), M),
         "bad action raises");
      Check
        (Running_Cost_Raises (G, Sid (1), Aid (1), Bad_M),
         "bad dist size raises");
      Check (Require_Compatible_Distribution (G, M), "compatible M");
      Pi := [Aid (1), Aid (2)];
      Check (Is_Valid_Policy (G, Pi), "valid policy");
      Mix := Pure_Soft_Policy (G, Pi);
      Check (Near (Mix (Sid (1), Aid (1)), 1.0), "pure soft (1,1)");
      Check (Near (Mix (Sid (1), Aid (2)), 0.0), "pure soft (1,2)");
      Check (Near (Mix (Sid (2), Aid (2)), 1.0), "pure soft (2,2)");
      Mix := Uniform_Soft_Policy (G);
      Check (Near (Mix (Sid (1), Aid (1)), 0.5), "uniform soft");
      Check (Near (Mix (Sid (2), Aid (2)), 0.5), "uniform soft 2");
      declare
         Bad_Pi : constant Policy (1 .. 3) := [1, 1, 1];
      begin
         Check (not Is_Valid_Policy (G, Bad_Pi), "wrong length policy");
      end;
   end Test_Cost_And_Policy;

   ---------------------------------------------------------------------------
   -- Best response / value iteration
   ---------------------------------------------------------------------------

   procedure Test_Best_Response is
      G : constant MFG := Congestion_Toy;
      M : constant Distribution := Uniform_Distribution (SC (3));
      Crowded : constant Distribution (1 .. 3) := [0.8, 0.1, 0.1];
      Pi : Policy (1 .. 3);
      V  : Value_Vector (1 .. 3);
      Iters : Natural;
      Bad : constant Distribution (1 .. 2) := [0.5, 0.5];
   begin
      Section ("Best_Response / Value_Iteration");
      Value_Iteration (G, M, V, Pi, Iters => Iters);
      Check (Iters > 0, "VI ran at least one iter");
      Check (Is_Valid_Policy (G, Pi), "VI policy valid");
      for S in State_Id range 1 .. 3 loop
         Check (Rf (V (S)) < 1.0E6, "value not exploding");
         Check (Rf (V (S)) > -1.0E6, "value not collapsing");
      end loop;
      Pi := Best_Response_Policy (G, M);
      Check (Is_Valid_Policy (G, Pi), "BR policy valid");
      V := Best_Response_Values (G, M);
      for S in 1 .. G.N_States loop
         declare
            Q_Best : constant Cost := Q_Value (G, M, V, S, Pi (S));
         begin
            Check (Near (Q_Best, V (S), 1.0E-6), "Bellman greediness");
            for A in 1 .. G.N_Actions loop
               Check
                 (Q_Value (G, M, V, S, A) >= V (S) - 1.0E-6,
                  "no improving deviation");
            end loop;
         end;
      end loop;
      --  Crowded state 1: same-state congestion is identical for both
      --  actions; incentive to move is through future value. BR from a
      --  heavily crowded site should prefer Move (action 2).
      Pi := Best_Response_Policy (G, Crowded);
      Check (Is_Valid_Policy (G, Pi), "crowded BR valid");
      Check (Pi (Sid (1)) = Aid (2), "crowded: BR moves away from state 1");
      declare
         Stay_C : constant Cost :=
           Running_Cost (G, Sid (1), Aid (1), Crowded);
         Uniform_C : constant Cost :=
           Running_Cost (G, Sid (1), Aid (1), M);
      begin
         Check (Stay_C > Uniform_C, "crowded: local cost exceeds uniform");
      end;
      V := Policy_Evaluation (G, M, Pi);
      Check (Rf (V (Sid (1))) < 1.0E6, "PE returns values");
      for S in 1 .. G.N_States loop
         Check (Near (Q_Value (G, M, V, S, Pi (S)), V (S), 1.0E-5),
                "PE Bellman for Pi");
      end loop;
      declare
         Mix : constant Soft_Policy := Soft_Best_Response (G, M, Beta => 5.0);
         Row_Sum : Probability;
      begin
         for S in 1 .. G.N_States loop
            Row_Sum := 0.0;
            for A in 1 .. G.N_Actions loop
               Check (Mix (S, A) >= -1.0E-12, "soft BR nonnegative");
               Row_Sum := Row_Sum + Mix (S, A);
            end loop;
            Check (Near (Row_Sum, 1.0, 1.0E-8), "soft BR row simplex");
         end loop;
      end;
      Check (BR_Raises (G, Bad), "BR bad M raises");
   end Test_Best_Response;

   ---------------------------------------------------------------------------
   -- Occupancy evolution
   ---------------------------------------------------------------------------

   procedure Test_Occupancy is
      G : constant MFG := Two_State_Symmetric;
      M : constant Distribution := Dirac_Distribution (SC (2), Sid (1));
      Pi_Stay : constant Policy (1 .. 2) := [Aid (1), Aid (1)];
      Pi_Switch : constant Policy (1 .. 2) := [Aid (2), Aid (2)];
      M2, M3 : Distribution (1 .. 2);
      Mix : Soft_Policy (1 .. 2, 1 .. 2);
      K : Kernel_Matrix (1 .. 2, 1 .. 2);
      Bad_Pi : constant Policy (1 .. 3) := [1, 1, 1];
   begin
      Section ("Evolve_Occupancy");
      M2 := Evolve_Occupancy (G, M, Pi_Stay);
      Check (Is_Simplex (M2), "evolved stay is simplex");
      Check (M2 (Sid (1)) > M2 (Sid (2)), "stay keeps mass mostly at 1");
      M2 := Evolve_Occupancy (G, M, Pi_Switch);
      Check (M2 (Sid (2)) > M2 (Sid (1)), "switch moves mass to 2");
      M3 := Evolve_Occupancy_N (G, M, Pi_Stay, 0);
      Check (Near (L1_Distance (M3, M), 0.0), "N=0 identity");
      M3 := Evolve_Occupancy_N (G, M, Pi_Stay, 5);
      Check (Is_Simplex (M3), "N=5 simplex");
      Mix := Uniform_Soft_Policy (G);
      M2 := Evolve_Occupancy_Soft (G, M, Mix);
      Check (Is_Simplex (M2), "soft evolve simplex");
      K := Induced_Kernel (G, Pi_Stay);
      Check (Near (K (Sid (1), Sid (1)), G.Trans (1, 1, 1)), "kernel match");
      Check
        (Near (K (Sid (1), Sid (1)) + K (Sid (1), Sid (2)), 1.0),
         "kernel row stochastic");
      Check (Evolve_Raises (G, M, Bad_Pi), "evolve bad policy raises");
      --  Mass conservation over many steps
      declare
         Cur : Distribution (1 .. 2) := Uniform_Distribution (SC (2));
      begin
         for I in 1 .. 20 loop
            Cur := Evolve_Occupancy (G, Cur, Pi_Stay);
         end loop;
         Check (Is_Simplex (Cur, 1.0E-8), "long-run simplex stay");
         Cur := Uniform_Distribution (SC (2));
         for I in 1 .. 20 loop
            Cur := Evolve_Occupancy (G, Cur, Pi_Switch);
         end loop;
         Check (Is_Simplex (Cur, 1.0E-8), "long-run simplex switch");
      end;
   end Test_Occupancy;

   ---------------------------------------------------------------------------
   -- Fixed point / Solve
   ---------------------------------------------------------------------------

   procedure Test_Solve is
      G : constant MFG := Congestion_Toy;
      E : constant MFG := Entry_Toy;
      F : constant MFG := Flocking_Toy;
      T : constant MFG := Two_State_Symmetric;
      M0 : constant Distribution := Uniform_Distribution (SC (3));
      M0e : constant Distribution := Uniform_Distribution (SC (2));
      R : Solve_Result (3, 2);
      Re : Solve_Result (2, 2);
      Rf3 : Solve_Result (3, 3);
   begin
      Section ("Fixed_Point_Iterate / Solve");
      R := Solve (G, M0, Damping => 0.5, Tol => 1.0E-6, Max_Iters => 200);
      Check (R.Iterations > 0, "congestion solve ran");
      Check (Is_Simplex (R.Mean_Field, 1.0E-6), "solve m simplex");
      Check (Is_Valid_Policy (G, R.Pi), "solve policy valid");
      Check (R.Residual >= 0.0, "residual nonnegative");
      Check
        (Consistency_Residual (G, R.Mean_Field, R.Pi) < 0.05
           or else R.Converged,
         "near-consistent or converged");
      --  Congestion equilibrium should not pile everything on one state.
      Check
        (R.Mean_Field (Sid (1)) < 0.95
           and then R.Mean_Field (Sid (2)) < 0.95
           and then R.Mean_Field (Sid (3)) < 0.95,
         "congestion not fully concentrated");
      Re := Solve (E, M0e, Damping => 0.4, Tol => 1.0E-6, Max_Iters => 300);
      Check (Is_Simplex (Re.Mean_Field, 1.0E-6), "entry m simplex");
      Check (Is_Valid_Policy (E, Re.Pi), "entry policy valid");
      --  With strong congestion on In, mass In should be bounded.
      Check (Re.Mean_Field (Sid (2)) < 0.9, "entry not fully In");
      Rf3 := Solve (F, M0, Damping => 0.35, Tol => 1.0E-5, Max_Iters => 400);
      Check (Is_Simplex (Rf3.Mean_Field, 1.0E-5), "flocking m simplex");
      Check (Is_Valid_Policy (F, Rf3.Pi), "flocking policy valid");
      declare
         Rt : constant Solve_Result :=
           Solve (T, M0e, Damping => 0.5, Tol => 1.0E-7, Max_Iters => 250);
      begin
         Check (Rt.Iterations > 0, "two-state solve ran");
         Check (Is_Simplex (Rt.Mean_Field, 1.0E-6), "two-state m");
         --  Symmetry: equal congestion → near-equal mass.
         Check
           (abs (Rt.Mean_Field (Sid (1)) - Rt.Mean_Field (Sid (2))) < 0.15,
            "symmetric congestion near-balanced");
      end;
      Check (Solve_Damping_Raises (G, M0, 0.0), "damping 0 raises");
      Check (Solve_Damping_Raises (G, M0, 1.5), "damping >1 raises");
      Check (Solve_Occ_Raises (G, M0), "occupancy_steps 0 raises");
      --  Damping = 1 allowed
      declare
         R1 : constant Solve_Result :=
           Solve (T, M0e, Damping => 1.0, Tol => 1.0E-5, Max_Iters => 100);
      begin
         Check (R1.Iterations > 0, "damping 1 runs");
         Check (Is_Simplex (R1.Mean_Field, 1.0E-5), "damping 1 simplex");
      end;
      --  Consistency residual identity
      declare
         Pi : constant Policy := Best_Response_Policy (G, M0);
         Resid : constant Real := Consistency_Residual (G, M0, Pi);
      begin
         Check (Resid >= 0.0, "consistency residual >= 0");
         Check
           (Near (Resid, L1_Distance (M0, Evolve_Occupancy (G, M0, Pi))),
            "consistency = L1(m, evolve)");
      end;
   end Test_Solve;

   ---------------------------------------------------------------------------
   -- Extra numerical / edge cases
   ---------------------------------------------------------------------------

   procedure Test_Edges is
      G : constant MFG := Entry_Toy;
      M : Distribution (1 .. 2);
      Pi : Policy (1 .. 2);
      V : Value_Vector (1 .. 2);
      Iters : Natural;
   begin
      Section ("Edge cases");
      M := Dirac_Distribution (SC (2), Sid (1));
      Value_Iteration (G, M, V, Pi, Tol => 1.0E-10, Max_Iters => 5000,
                       Iters => Iters);
      Check (Iters >= 1, "VI iters");
      Check (Is_Valid_Policy (G, Pi), "VI policy from dirac Out");
      M := Dirac_Distribution (SC (2), Sid (2));
      Value_Iteration (G, M, V, Pi, Iters => Iters);
      Check (Is_Valid_Policy (G, Pi), "VI policy from dirac In");
      --  Soft evolve with pure policy matches hard evolve
      declare
         Hard : constant Distribution := Evolve_Occupancy (G, M, Pi);
         Soft : constant Distribution :=
           Evolve_Occupancy_Soft (G, M, Pure_Soft_Policy (G, Pi));
      begin
         Check (Near (L1_Distance (Hard, Soft), 0.0, 1.0E-12),
                "soft pure = hard evolve");
      end;
      --  Invalid MFG: break a transition row
      declare
         Bad : MFG := G;
      begin
         Bad.Trans (1, 1, 1) := 0.0;
         Bad.Trans (1, 1, 2) := 0.0;
         Check (not Is_Valid_Transitions (Bad), "broken row invalid");
         Check (not Is_Valid_MFG (Bad), "broken MFG invalid");
         Bad := G;
         Bad.Discount := 1.0;
         Check (not Is_Valid_MFG (Bad), "discount 1 invalid");
         Bad := G;
         Bad.Discount := -0.1;
         Check (not Is_Valid_MFG (Bad), "negative discount invalid");
      end;
      Check (not Near_Raises (0.0), "tol 0 ok for Near");
      Check (Near (0.0, 0.0, 0.0), "exact with tol 0");
      --  Many BR calls under different m
      for K in 0 .. 10 loop
         declare
            T : constant Real := Real (K) / 10.0;
            Mk : constant Distribution :=
              [1.0 - T, T];
            Pk : constant Policy := Best_Response_Policy (G, Mk);
         begin
            Check (Is_Valid_Policy (G, Pk), "BR along path");
            Check (Is_Simplex (Mk), "path simplex");
         end;
      end loop;
   end Test_Edges;

   ---------------------------------------------------------------------------
   -- Batch property checks for coverage volume
   ---------------------------------------------------------------------------

   procedure Test_Batch_Properties is
      C : constant MFG := Congestion_Toy;
      E : constant MFG := Entry_Toy;
      F : constant MFG := Flocking_Toy;
      T : constant MFG := Two_State_Symmetric;
   begin
      Section ("Batch properties");
      --  Every transition row stochastic for all toys
      for S in 1 .. C.N_States loop
         for A in 1 .. C.N_Actions loop
            Check (Is_Stochastic_Row (C, S, A), "C stochastic row");
            declare
               Sum : Probability := 0.0;
            begin
               for Sp in 1 .. C.N_States loop
                  Sum := Sum + C.Trans (S, A, Sp);
                  Check (C.Trans (S, A, Sp) >= -1.0E-12, "C nonnegative");
               end loop;
               Check (Near (Sum, 1.0, 1.0E-9), "C row sum 1");
            end;
         end loop;
      end loop;
      for S in 1 .. E.N_States loop
         for A in 1 .. E.N_Actions loop
            Check (Is_Stochastic_Row (E, S, A), "E stochastic row");
         end loop;
      end loop;
      for S in 1 .. F.N_States loop
         for A in 1 .. F.N_Actions loop
            Check (Is_Stochastic_Row (F, S, A), "F stochastic row");
         end loop;
      end loop;
      for S in 1 .. T.N_States loop
         for A in 1 .. T.N_Actions loop
            Check (Is_Stochastic_Row (T, S, A), "T stochastic row");
         end loop;
      end loop;
      --  Normalize / Dirac / Uniform family
      for N in State_Count range 1 .. 5 loop
         declare
            U : constant Distribution := Uniform_Distribution (N);
         begin
            Check (Is_Simplex (U), "uniform simplex family");
            Check (Near (U (Sid (1)), 1.0 / Real (N)), "uniform value");
            for S in 1 .. N loop
               declare
                  D : constant Distribution := Dirac_Distribution (N, S);
               begin
                  Check (Is_Simplex (D), "dirac simplex family");
                  Check (Near (D (S), 1.0), "dirac peak family");
               end;
            end loop;
         end;
      end loop;
      --  L1 symmetry and triangle on small dists
      declare
         Dist_A : constant Distribution := [0.7, 0.3];
         Dist_B : constant Distribution := [0.4, 0.6];
         Dist_C : constant Distribution := [0.5, 0.5];
         D_AB : constant Real := L1_Distance (Dist_A, Dist_B);
         D_BA : constant Real := L1_Distance (Dist_B, Dist_A);
         D_AC : constant Real := L1_Distance (Dist_A, Dist_C);
         D_CB : constant Real := L1_Distance (Dist_C, Dist_B);
      begin
         Check (Near (X => D_AB, Y => D_BA), "L1 symmetric");
         Check (D_AC + D_CB >= D_AB - 1.0E-12, "L1 triangle");
         Check (Near (Total_Variation (Dist_A, Dist_A), 0.0), "TV zero");
      end;
      --  Q_Value linearity in V under zero discount
      declare
         G0 : MFG := T;
         M : constant Distribution := Uniform_Distribution (SC (2));
         V0 : constant Value_Vector (1 .. 2) := [1.0, 2.0];
         V1 : constant Value_Vector (1 .. 2) := [3.0, 4.0];
      begin
         G0.Discount := 0.0;
         for S in 1 .. G0.N_States loop
            for A in 1 .. G0.N_Actions loop
               Check
                 (Near
                    (Q_Value (G0, M, V0, S, A),
                     Running_Cost (G0, S, A, M)),
                  "Q = cost when gamma=0");
               Check
                 (Near
                    (Q_Value (G0, M, V1, S, A),
                     Running_Cost (G0, S, A, M)),
                  "Q ignores V when gamma=0");
            end loop;
         end loop;
      end;
   end Test_Batch_Properties;

   ---------------------------------------------------------------------------
   -- Solve from several initial conditions
   ---------------------------------------------------------------------------

   procedure Test_Solve_Inits is
      G : constant MFG := Two_State_Symmetric;
      Inits : constant array (1 .. 5) of Distribution (1 .. 2) :=
        [[0.5, 0.5],
         [0.9, 0.1],
         [0.1, 0.9],
         [0.75, 0.25],
         [0.25, 0.75]];
   begin
      Section ("Solve from varied initials");
      for I in Inits'Range loop
         declare
            R : constant Solve_Result :=
              Solve
                (G, Inits (I),
                 Damping => 0.3,
                 Tol => 1.0E-5,
                 Max_Iters => 500,
                 Occupancy_Steps => 3);
         begin
            Check (Is_Simplex (R.Mean_Field, 1.0E-6), "init solve simplex");
            Check (Is_Valid_Policy (G, R.Pi), "init solve policy");
            Check (R.Residual < 0.25, "init solve bounded residual");
            Check (R.Iterations >= 1, "init solve iterated");
            --  Symmetric costs: mass should not fully collapse.
            Check
              (R.Mean_Field (1) > 0.05 and then R.Mean_Field (2) > 0.05,
               "init solve not collapsed");
         end;
      end loop;
   end Test_Solve_Inits;

begin
   Put_Line ("Mean_Field_Game test suite (Ada 2023)");
   Test_Near_And_Simplex;
   Test_Constructors;
   Test_Cost_And_Policy;
   Test_Best_Response;
   Test_Occupancy;
   Test_Solve;
   Test_Edges;
   Test_Batch_Properties;
   Test_Solve_Inits;
   New_Line;
   Put_Line
     ("Result:" & Natural'Image (Pass_Count) & " PASS,"
      & Natural'Image (Fail_Count) & " FAIL");
   if Fail_Count > 0 then
      Ada.Command_Line.Set_Exit_Status (Ada.Command_Line.Failure);
   end if;
end Tests;
