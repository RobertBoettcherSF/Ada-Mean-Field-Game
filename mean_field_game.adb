--  Mean_Field_Game body — finite-state discrete-time MFG helpers.

pragma Ada_2022;

with Ada.Numerics.Long_Elementary_Functions;

package body Mean_Field_Game
  with SPARK_Mode => Off
is
   package Math renames Ada.Numerics.Long_Elementary_Functions;

   ---------------------------------------------------------------------------
   -- Local helpers
   ---------------------------------------------------------------------------

   procedure Check_Tol (Tol : Real) is
   begin
      if Tol < 0.0 then
         raise Invalid_Argument with "tolerance must be >= 0";
      end if;
   end Check_Tol;

   procedure Check_Nonempty (G : MFG) is
   begin
      if G.N_States = 0 or else G.N_Actions = 0 then
         raise Invalid_Argument with "MFG must have positive dimensions";
      end if;
   end Check_Nonempty;

   procedure Check_Discount (G : MFG) is
   begin
      if G.Discount < 0.0 or else G.Discount >= 1.0 then
         raise Invalid_Argument with "Discount must be in [0, 1)";
      end if;
   end Check_Discount;

   procedure Check_State (G : MFG; S : State_Id) is
   begin
      if S not in 1 .. G.N_States then
         raise Invalid_Argument with "state index out of range";
      end if;
   end Check_State;

   procedure Check_Action (G : MFG; A : Action_Id) is
   begin
      if A not in 1 .. G.N_Actions then
         raise Invalid_Argument with "action index out of range";
      end if;
   end Check_Action;

   procedure Check_Distribution (G : MFG; M : Distribution) is
   begin
      if M'First /= 1 or else M'Last /= G.N_States then
         raise Invalid_Argument with "distribution size mismatch";
      end if;
   end Check_Distribution;

   procedure Check_Policy (G : MFG; Pi : Policy) is
   begin
      if Pi'First /= 1 or else Pi'Last /= G.N_States then
         raise Invalid_Argument with "policy size mismatch";
      end if;
      for S in 1 .. G.N_States loop
         if Pi (S) not in 1 .. G.N_Actions then
            raise Invalid_Argument with "policy action out of range";
         end if;
      end loop;
   end Check_Policy;

   procedure Check_Soft (G : MFG; Mix : Soft_Policy) is
   begin
      if Mix'First (1) /= 1 or else Mix'Last (1) /= G.N_States
        or else Mix'First (2) /= 1 or else Mix'Last (2) /= G.N_Actions
      then
         raise Invalid_Argument with "soft policy size mismatch";
      end if;
   end Check_Soft;

   function Abs_R (X : Real) return Real is
   begin
      if X >= 0.0 then
         return X;
      else
         return -X;
      end if;
   end Abs_R;

   ---------------------------------------------------------------------------
   -- Near
   ---------------------------------------------------------------------------

   function Near
     (X, Y : Real; Tol : Real := Default_Tol) return Boolean
   is
   begin
      Check_Tol (Tol);
      return Abs_R (X - Y) <= Tol;
   end Near;

   ---------------------------------------------------------------------------
   -- Distribution helpers
   ---------------------------------------------------------------------------

   function Mass (M : Distribution) return Probability is
      S : Probability := 0.0;
   begin
      if M'Length = 0 then
         return 0.0;
      end if;
      for I in M'Range loop
         S := S + M (I);
      end loop;
      return S;
   end Mass;

   function Is_Simplex
     (M : Distribution; Tol : Real := Default_Tol) return Boolean
   is
   begin
      Check_Tol (Tol);
      if M'Length = 0 then
         return False;
      end if;
      for I in M'Range loop
         if M (I) < -Tol then
            return False;
         end if;
      end loop;
      return Abs_R (Mass (M) - 1.0) <= Tol;
   end Is_Simplex;

   function Normalize (M : Distribution) return Distribution is
      Total : constant Probability := Mass (M);
      R     : Distribution (M'Range);
   begin
      if Total <= 0.0 then
         raise Invalid_Argument with "cannot normalize non-positive mass";
      end if;
      for I in M'Range loop
         R (I) := M (I) / Total;
      end loop;
      return R;
   end Normalize;

   function Uniform_Distribution (N : State_Count) return Distribution is
      R : Distribution (1 .. N);
   begin
      if N = 0 then
         raise Invalid_Argument with "Uniform_Distribution: N = 0";
      end if;
      for I in R'Range loop
         R (I) := 1.0 / Real (N);
      end loop;
      return R;
   end Uniform_Distribution;

   function Dirac_Distribution
     (N : State_Count; S : State_Id) return Distribution
   is
      R : Distribution (1 .. N) := [others => 0.0];
   begin
      if N = 0 then
         raise Invalid_Argument with "Dirac_Distribution: N = 0";
      end if;
      if S not in 1 .. N then
         raise Invalid_Argument with "Dirac_Distribution: state out of range";
      end if;
      R (S) := 1.0;
      return R;
   end Dirac_Distribution;

   function L1_Distance (A, B : Distribution) return Real is
      D : Real := 0.0;
   begin
      if A'First /= B'First or else A'Last /= B'Last then
         raise Invalid_Argument with "L1_Distance: length mismatch";
      end if;
      for I in A'Range loop
         D := D + Abs_R (A (I) - B (I));
      end loop;
      return D;
   end L1_Distance;

   function Total_Variation (A, B : Distribution) return Real is
   begin
      return 0.5 * L1_Distance (A, B);
   end Total_Variation;

   ---------------------------------------------------------------------------
   -- Transition / cost
   ---------------------------------------------------------------------------

   function Is_Stochastic_Row
     (G   : MFG;
      S   : State_Id;
      A   : Action_Id;
      Tol : Real := Default_Tol) return Boolean
   is
      Sum : Probability := 0.0;
   begin
      Check_Tol (Tol);
      Check_Nonempty (G);
      Check_State (G, S);
      Check_Action (G, A);
      for Sp in 1 .. G.N_States loop
         if G.Trans (S, A, Sp) < -Tol then
            return False;
         end if;
         Sum := Sum + G.Trans (S, A, Sp);
      end loop;
      return Abs_R (Sum - 1.0) <= Tol;
   end Is_Stochastic_Row;

   function Is_Valid_Transitions
     (G : MFG; Tol : Real := Default_Tol) return Boolean
   is
   begin
      Check_Tol (Tol);
      if G.N_States = 0 or else G.N_Actions = 0 then
         return False;
      end if;
      for S in 1 .. G.N_States loop
         for A in 1 .. G.N_Actions loop
            if not Is_Stochastic_Row (G, S, A, Tol) then
               return False;
            end if;
         end loop;
      end loop;
      return True;
   end Is_Valid_Transitions;

   function Running_Cost
     (G : MFG; S : State_Id; A : Action_Id; M : Distribution) return Cost
   is
   begin
      Check_Nonempty (G);
      Check_Distribution (G, M);
      Check_State (G, S);
      Check_Action (G, A);
      return G.Base_Cost (S, A) + G.Congestion (S) * M (S);
   end Running_Cost;

   function Require_Compatible_Distribution
     (G : MFG; M : Distribution) return Boolean
   is
   begin
      Check_Distribution (G, M);
      return True;
   end Require_Compatible_Distribution;

   ---------------------------------------------------------------------------
   -- Policies
   ---------------------------------------------------------------------------

   function Is_Valid_Policy (G : MFG; Pi : Policy) return Boolean is
   begin
      if G.N_States = 0 or else G.N_Actions = 0 then
         return False;
      end if;
      if Pi'First /= 1 or else Pi'Last /= G.N_States then
         return False;
      end if;
      for S in 1 .. G.N_States loop
         if Pi (S) not in 1 .. G.N_Actions then
            return False;
         end if;
      end loop;
      return True;
   end Is_Valid_Policy;

   function Pure_Soft_Policy
     (G : MFG; Pi : Policy) return Soft_Policy
   is
      Mix : Soft_Policy (1 .. G.N_States, 1 .. G.N_Actions) :=
        [others => [others => 0.0]];
   begin
      Check_Nonempty (G);
      Check_Policy (G, Pi);
      for S in 1 .. G.N_States loop
         Mix (S, Pi (S)) := 1.0;
      end loop;
      return Mix;
   end Pure_Soft_Policy;

   function Uniform_Soft_Policy (G : MFG) return Soft_Policy is
      Mix : Soft_Policy (1 .. G.N_States, 1 .. G.N_Actions);
      P   : Probability;
   begin
      Check_Nonempty (G);
      P := 1.0 / Real (G.N_Actions);
      for S in 1 .. G.N_States loop
         for A in 1 .. G.N_Actions loop
            Mix (S, A) := P;
         end loop;
      end loop;
      return Mix;
   end Uniform_Soft_Policy;

   ---------------------------------------------------------------------------
   -- Bellman / value iteration
   ---------------------------------------------------------------------------

   function Q_Value
     (G : MFG;
      M : Distribution;
      V : Value_Vector;
      S : State_Id;
      A : Action_Id) return Cost
   is
      Exp : Cost := 0.0;
   begin
      Check_Nonempty (G);
      Check_Distribution (G, M);
      if V'First /= 1 or else V'Last /= G.N_States then
         raise Invalid_Argument with "value vector size mismatch";
      end if;
      Check_State (G, S);
      Check_Action (G, A);
      for Sp in 1 .. G.N_States loop
         Exp := Exp + G.Trans (S, A, Sp) * V (Sp);
      end loop;
      return Running_Cost (G, S, A, M) + G.Discount * Exp;
   end Q_Value;

   procedure Value_Iteration
     (G         : MFG;
      M         : Distribution;
      Values    : out Value_Vector;
      Pi        : out Policy;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000;
      Iters     : out Natural)
   is
      V      : Value_Vector (1 .. G.N_States) := [others => 0.0];
      V_Next : Value_Vector (1 .. G.N_States);
      Max_Diff : Real;
      Best_A : Action_Id;
      Best_Q : Cost;
      Q      : Cost;
   begin
      Check_Nonempty (G);
      Check_Discount (G);
      Check_Distribution (G, M);
      Check_Tol (Tol);
      if Values'First /= 1 or else Values'Last /= G.N_States
        or else Pi'First /= 1 or else Pi'Last /= G.N_States
      then
         raise Invalid_Argument with "Value_Iteration output size mismatch";
      end if;

      Iters := 0;
      loop
         Max_Diff := 0.0;
         for S in 1 .. G.N_States loop
            Best_A := 1;
            Best_Q := Q_Value (G, M, V, S, 1);
            for A in Action_Id range 2 .. G.N_Actions loop
               Q := Q_Value (G, M, V, S, A);
               if Q < Best_Q then
                  Best_Q := Q;
                  Best_A := A;
               end if;
            end loop;
            V_Next (S) := Best_Q;
            Pi (S) := Best_A;
            declare
               D : constant Real := Abs_R (V_Next (S) - V (S));
            begin
               if D > Max_Diff then
                  Max_Diff := D;
               end if;
            end;
         end loop;
         V := V_Next;
         Iters := Iters + 1;
         exit when Max_Diff <= Tol or else Iters >= Max_Iters;
      end loop;
      Values := V;
   end Value_Iteration;

   function Best_Response_Policy
     (G         : MFG;
      M         : Distribution;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Policy
   is
      Values : Value_Vector (1 .. G.N_States);
      Pi     : Policy (1 .. G.N_States);
      Iters  : Natural;
   begin
      Value_Iteration (G, M, Values, Pi, Tol, Max_Iters, Iters);
      return Pi;
   end Best_Response_Policy;

   function Best_Response_Values
     (G         : MFG;
      M         : Distribution;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Value_Vector
   is
      Values : Value_Vector (1 .. G.N_States);
      Pi     : Policy (1 .. G.N_States);
      Iters  : Natural;
   begin
      Value_Iteration (G, M, Values, Pi, Tol, Max_Iters, Iters);
      return Values;
   end Best_Response_Values;


   function Soft_Best_Response
     (G         : MFG;
      M         : Distribution;
      Beta      : Real := 5.0;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Soft_Policy
   is
      V   : Value_Vector (1 .. G.N_States);
      Pi  : Policy (1 .. G.N_States);
      It  : Natural;
      Mix : Soft_Policy (1 .. G.N_States, 1 .. G.N_Actions);
   begin
      if Beta <= 0.0 then
         raise Invalid_Argument with "Soft_Best_Response: Beta must be > 0";
      end if;
      Value_Iteration (G, M, V, Pi, Tol, Max_Iters, It);
      for S in 1 .. G.N_States loop
         declare
            Qs   : array (Action_Id range 1 .. G.N_Actions) of Cost;
            Mmin : Cost;
            Sum  : Probability := 0.0;
         begin
            for A in 1 .. G.N_Actions loop
               Qs (A) := Q_Value (G, M, V, S, A);
            end loop;
            Mmin := Qs (1);
            for A in Action_Id range 2 .. G.N_Actions loop
               if Qs (A) < Mmin then
                  Mmin := Qs (A);
               end if;
            end loop;
            for A in 1 .. G.N_Actions loop
               declare
                  Shift : Real := -Beta * (Qs (A) - Mmin);
                  E     : Real;
               begin
                  if Shift > 80.0 then
                     Shift := 80.0;
                  elsif Shift < -80.0 then
                     Shift := -80.0;
                  end if;
                  E := Real (Math.Exp (Long_Float (Shift)));
                  Mix (S, A) := E;
                  Sum := Sum + E;
               end;
            end loop;
            if Sum <= 0.0 then
               for A in 1 .. G.N_Actions loop
                  Mix (S, A) := 1.0 / Real (G.N_Actions);
               end loop;
            else
               for A in 1 .. G.N_Actions loop
                  Mix (S, A) := Mix (S, A) / Sum;
               end loop;
            end if;
         end;
      end loop;
      return Mix;
   end Soft_Best_Response;

   function Policy_Evaluation
     (G         : MFG;
      M         : Distribution;
      Pi        : Policy;
      Tol       : Real := Default_Tol;
      Max_Iters : Natural := 10_000) return Value_Vector
   is
      V      : Value_Vector (1 .. G.N_States) := [others => 0.0];
      V_Next : Value_Vector (1 .. G.N_States);
      Max_Diff : Real;
      Iters  : Natural := 0;
   begin
      Check_Nonempty (G);
      Check_Discount (G);
      Check_Distribution (G, M);
      Check_Policy (G, Pi);
      Check_Tol (Tol);
      loop
         Max_Diff := 0.0;
         for S in 1 .. G.N_States loop
            V_Next (S) := Q_Value (G, M, V, S, Pi (S));
            declare
               D : constant Real := Abs_R (V_Next (S) - V (S));
            begin
               if D > Max_Diff then
                  Max_Diff := D;
               end if;
            end;
         end loop;
         V := V_Next;
         Iters := Iters + 1;
         exit when Max_Diff <= Tol or else Iters >= Max_Iters;
      end loop;
      return V;
   end Policy_Evaluation;

   ---------------------------------------------------------------------------
   -- Occupancy updates
   ---------------------------------------------------------------------------

   function Evolve_Occupancy
     (G  : MFG;
      M  : Distribution;
      Pi : Policy) return Distribution
   is
      M_Next : Distribution (1 .. G.N_States) := [others => 0.0];
   begin
      Check_Nonempty (G);
      Check_Distribution (G, M);
      Check_Policy (G, Pi);
      for S in 1 .. G.N_States loop
         declare
            A : constant Action_Id := Pi (S);
         begin
            for Sp in 1 .. G.N_States loop
               M_Next (Sp) := M_Next (Sp) + M (S) * G.Trans (S, A, Sp);
            end loop;
         end;
      end loop;
      return M_Next;
   end Evolve_Occupancy;

   function Evolve_Occupancy_Soft
     (G   : MFG;
      M   : Distribution;
      Mix : Soft_Policy) return Distribution
   is
      M_Next : Distribution (1 .. G.N_States) := [others => 0.0];
   begin
      Check_Nonempty (G);
      Check_Distribution (G, M);
      Check_Soft (G, Mix);
      for S in 1 .. G.N_States loop
         for A in 1 .. G.N_Actions loop
            for Sp in 1 .. G.N_States loop
               M_Next (Sp) :=
                 M_Next (Sp) + M (S) * Mix (S, A) * G.Trans (S, A, Sp);
            end loop;
         end loop;
      end loop;
      return M_Next;
   end Evolve_Occupancy_Soft;

   function Evolve_Occupancy_N
     (G     : MFG;
      M     : Distribution;
      Pi    : Policy;
      Steps : Natural) return Distribution
   is
      Cur : Distribution (1 .. G.N_States);
   begin
      Check_Nonempty (G);
      Check_Distribution (G, M);
      Check_Policy (G, Pi);
      Cur := M;
      for K in 1 .. Steps loop
         Cur := Evolve_Occupancy (G, Cur, Pi);
      end loop;
      return Cur;
   end Evolve_Occupancy_N;

   function Induced_Kernel
     (G  : MFG;
      Pi : Policy) return Kernel_Matrix
   is
      K : Kernel_Matrix (1 .. G.N_States, 1 .. G.N_States);
   begin
      Check_Nonempty (G);
      Check_Policy (G, Pi);
      for S in 1 .. G.N_States loop
         for Sp in 1 .. G.N_States loop
            K (S, Sp) := G.Trans (S, Pi (S), Sp);
         end loop;
      end loop;
      return K;
   end Induced_Kernel;

   ---------------------------------------------------------------------------
   -- Fixed-point iteration
   ---------------------------------------------------------------------------

   function Consistency_Residual
     (G  : MFG;
      M  : Distribution;
      Pi : Policy) return Real
   is
   begin
      return L1_Distance (M, Evolve_Occupancy (G, M, Pi));
   end Consistency_Residual;

   procedure Fixed_Point_Iterate
     (G          : MFG;
      M0         : Distribution;
      Result     : out Solve_Result;
      Damping    : Real := 0.5;
      Tol        : Real := 1.0E-8;
      Max_Iters  : Natural := 500;
      Occupancy_Steps : Natural := 1)
   is
      M      : Distribution (1 .. G.N_States);
      M_New  : Distribution (1 .. G.N_States);
      Pi     : Policy (1 .. G.N_States);
      Values : Value_Vector (1 .. G.N_States);
      Iters  : Natural;
      Resid  : Real;
      Alpha  : Real;
   begin
      Check_Nonempty (G);
      Check_Discount (G);
      Check_Distribution (G, M0);
      Check_Tol (Tol);
      if Damping <= 0.0 or else Damping > 1.0 then
         raise Invalid_Argument with "Damping must be in (0, 1]";
      end if;
      if Occupancy_Steps = 0 then
         raise Invalid_Argument with "Occupancy_Steps must be >= 1";
      end if;
      if Result.N_States /= G.N_States
        or else Result.N_Actions /= G.N_Actions
      then
         raise Invalid_Argument with "Solve_Result discriminant mismatch";
      end if;

      Alpha := Damping;
      M := M0;
      Result.Converged := False;
      Result.Iterations := 0;
      Result.Residual := Real'Last;

      for It in 1 .. Max_Iters loop
         declare
            M_Prev : constant Distribution := M;
            Mix    : Soft_Policy (1 .. G.N_States, 1 .. G.N_Actions);
            Cur    : Distribution (1 .. G.N_States);
         begin
            --  Soft BR for occupancy flow; hard BR stored in Result.Pi.
            Value_Iteration (G, M, Values, Pi, Default_Tol, 10_000, Iters);
            Mix := Soft_Best_Response (G, M, Beta => 8.0);
            Cur := M;
            for Step in 1 .. Occupancy_Steps loop
               Cur := Evolve_Occupancy_Soft (G, Cur, Mix);
            end loop;
            M_New := Cur;
            for S in 1 .. G.N_States loop
               M (S) := (1.0 - Alpha) * M (S) + Alpha * M_New (S);
            end loop;
            --  Iteration residual: change in m (fixed point of the map).
            Resid := L1_Distance (M_Prev, M);
            Result.Iterations := It;
            Result.Residual := Resid;
            Result.Mean_Field := M;
            Result.Pi := Pi;
            Result.Values := Values;
            if Resid <= Tol then
               Result.Converged := True;
               exit;
            end if;
         end;
      end loop;
   end Fixed_Point_Iterate;

   function Solve
     (G          : MFG;
      M0         : Distribution;
      Damping    : Real := 0.5;
      Tol        : Real := 1.0E-8;
      Max_Iters  : Natural := 500;
      Occupancy_Steps : Natural := 1) return Solve_Result
   is
      R : Solve_Result (G.N_States, G.N_Actions);
   begin
      Fixed_Point_Iterate
        (G, M0, R, Damping, Tol, Max_Iters, Occupancy_Steps);
      return R;
   end Solve;

   ---------------------------------------------------------------------------
   -- Classic constructors
   ---------------------------------------------------------------------------

   --  Helper: zero-init then fill.
   function Blank (Ns : State_Count; Na : Action_Count) return MFG is
      G : MFG (Ns, Na);
   begin
      G.Trans := [others => [others => [others => 0.0]]];
      G.Base_Cost := [others => [others => 0.0]];
      G.Congestion := [others => 0.0];
      G.Discount := 0.9;
      return G;
   end Blank;

   procedure Set_Row
     (G : in out MFG; S : State_Id; A : Action_Id; Row : Distribution)
   is
      Sum : Probability := 0.0;
   begin
      if Row'Length /= Natural (G.N_States) then
         raise Invalid_Argument with "Set_Row length mismatch";
      end if;
      for Sp in 1 .. G.N_States loop
         G.Trans (S, A, Sp) := Row (Row'First + (Sp - 1));
         Sum := Sum + G.Trans (S, A, Sp);
      end loop;
      if Sum <= 0.0 then
         raise Invalid_Argument with "Set_Row zero mass";
      end if;
      for Sp in 1 .. G.N_States loop
         G.Trans (S, A, Sp) := G.Trans (S, A, Sp) / Sum;
      end loop;
   end Set_Row;

   function Congestion_Toy return MFG is
      --  States 1..3 = locations. Actions: 1=Stay, 2=Move (clockwise).
      G : MFG := Blank (3, 2);
   begin
      G.Discount := 0.95;
      --  Stay: mostly remain, small leak to neighbours.
      Set_Row (G, 1, 1, Distribution'(1 => 0.85, 2 => 0.10, 3 => 0.05));
      Set_Row (G, 2, 1, Distribution'(1 => 0.05, 2 => 0.85, 3 => 0.10));
      Set_Row (G, 3, 1, Distribution'(1 => 0.10, 2 => 0.05, 3 => 0.85));
      --  Move: rotate 1→2→3→1 with small noise.
      Set_Row (G, 1, 2, Distribution'(1 => 0.05, 2 => 0.90, 3 => 0.05));
      Set_Row (G, 2, 2, Distribution'(1 => 0.05, 2 => 0.05, 3 => 0.90));
      Set_Row (G, 3, 2, Distribution'(1 => 0.90, 2 => 0.05, 3 => 0.05));
      --  Base move cost; congestion penalises crowded locations.
      for S in State_Id range 1 .. 3 loop
         G.Base_Cost (S, 1) := 0.1;
         G.Base_Cost (S, 2) := 0.5;
         G.Congestion (S) := 2.0;
      end loop;
      return G;
   end Congestion_Toy;

   function Flocking_Toy return MFG is
      --  Line 1—2—3. Actions: 1=Left, 2=Stay, 3=Right.
      --  Attraction to crowd: negative congestion (reward for being where
      --  others are) plus movement friction.
      G : MFG := Blank (3, 3);
   begin
      G.Discount := 0.9;
      --  Left
      Set_Row (G, 1, 1, Distribution'(1 => 1.0, 2 => 0.0, 3 => 0.0));
      Set_Row (G, 2, 1, Distribution'(1 => 0.9, 2 => 0.1, 3 => 0.0));
      Set_Row (G, 3, 1, Distribution'(1 => 0.0, 2 => 0.9, 3 => 0.1));
      --  Stay
      Set_Row (G, 1, 2, Distribution'(1 => 0.9, 2 => 0.1, 3 => 0.0));
      Set_Row (G, 2, 2, Distribution'(1 => 0.05, 2 => 0.9, 3 => 0.05));
      Set_Row (G, 3, 2, Distribution'(1 => 0.0, 2 => 0.1, 3 => 0.9));
      --  Right
      Set_Row (G, 1, 3, Distribution'(1 => 0.1, 2 => 0.9, 3 => 0.0));
      Set_Row (G, 2, 3, Distribution'(1 => 0.0, 2 => 0.1, 3 => 0.9));
      Set_Row (G, 3, 3, Distribution'(1 => 0.0, 2 => 0.0, 3 => 1.0));
      for S in State_Id range 1 .. 3 loop
         G.Base_Cost (S, 1) := 0.3;  -- Left
         G.Base_Cost (S, 2) := 0.05; -- Stay cheap
         G.Base_Cost (S, 3) := 0.3;  -- Right
         G.Congestion (S) := -1.0;   -- flocking reward ∝ m(s)
      end loop;
      return G;
   end Flocking_Toy;

   function Entry_Toy return MFG is
      --  States: 1=Out, 2=In. Actions: 1=Stay, 2=Switch.
      G : MFG := Blank (2, 2);
   begin
      G.Discount := 0.92;
      --  Stay
      Set_Row (G, 1, 1, Distribution'(1 => 1.0, 2 => 0.0));
      Set_Row (G, 2, 1, Distribution'(1 => 0.0, 2 => 1.0));
      --  Switch
      Set_Row (G, 1, 2, Distribution'(1 => 0.05, 2 => 0.95));
      Set_Row (G, 2, 2, Distribution'(1 => 0.95, 2 => 0.05));
      --  Out is mildly costly (missed opportunity); In has negative base
      --  cost (profit) but strong congestion when crowded.
      G.Base_Cost (1, 1) := 0.2;
      G.Base_Cost (1, 2) := 0.4;   -- switch cost from Out
      G.Base_Cost (2, 1) := -0.5;  -- profit when In
      G.Base_Cost (2, 2) := -0.2;  -- switch away (lose some profit + fee)
      G.Congestion (1) := 0.0;
      G.Congestion (2) := 3.0;     -- market congestion
      return G;
   end Entry_Toy;

   function Two_State_Symmetric return MFG is
      G : MFG := Blank (2, 2);
   begin
      G.Discount := 0.8;
      --  Stay
      Set_Row (G, 1, 1, Distribution'(1 => 0.9, 2 => 0.1));
      Set_Row (G, 2, 1, Distribution'(1 => 0.1, 2 => 0.9));
      --  Switch
      Set_Row (G, 1, 2, Distribution'(1 => 0.1, 2 => 0.9));
      Set_Row (G, 2, 2, Distribution'(1 => 0.9, 2 => 0.1));
      G.Base_Cost (1, 1) := 0.0;
      G.Base_Cost (1, 2) := 0.2;
      G.Base_Cost (2, 1) := 0.0;
      G.Base_Cost (2, 2) := 0.2;
      G.Congestion (1) := 1.5;
      G.Congestion (2) := 1.5;
      return G;
   end Two_State_Symmetric;

   ---------------------------------------------------------------------------
   -- Validation
   ---------------------------------------------------------------------------

   function Is_Valid_MFG
     (G : MFG; Tol : Real := Default_Tol) return Boolean
   is
   begin
      Check_Tol (Tol);
      if G.N_States = 0 or else G.N_Actions = 0 then
         return False;
      end if;
      if G.Discount < 0.0 or else G.Discount >= 1.0 then
         return False;
      end if;
      return Is_Valid_Transitions (G, Tol);
   end Is_Valid_MFG;

end Mean_Field_Game;
