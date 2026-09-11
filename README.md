# Mean-Field Games in Ada 2023

## Project Overview

A **mean-field game (MFG)** studies strategic decision-making by a continuum of
identical agents. It is the continuum limit of a large symmetric $N$-player
game as $N\to\infty$: each infinitesimal **representative agent** optimises
against the **population distribution** (mean field) $m$, and **consistency**
requires that the optimal feedback induces a flow of $m$ matching the assumed
field. The name echoes mean-field theory in physics (individual particles have
negligible impact on the aggregate).

Foundational continuous-time work is due independently to **Lasry–Lions** and
to **Huang–Caines–Malhamé**; related discrete ideas appear earlier in
economics (Jovanovic–Rosenthal). In continuous state/time an MFG couples a
Hamilton–Jacobi–Bellman equation for the representative value with a
Fokker–Planck equation for $m$. This package is deliberately
**numerical/educational** and **not** PDE-heavy: it implements **finite-state,
discrete-time** classroom toys.

States $s\in\{1,\ldots,S\}$, actions $a\in\{1,\ldots,A\}$, discount
$\gamma\in[0,1)$, running cost depending on own state/action and the mean
field (congestion: cost rises with $m(s)$), and transitions $P(s'\mid s,a)$
depending on $a$. Against a frozen $m$, the agent solves a standard MDP by
value iteration. The occupancy (discrete Fokker–Planck) update advances $m$
under the best-response policy. Alternating the two with damping yields a
numerical fixed point $(\pi^\star,m^\star)$.

Primary source:
[Wikipedia — Mean-field game theory](https://en.wikipedia.org/wiki/Mean-field_game_theory).

Part of the **RobertBoettcherSF** Ada algorithm series.

## Classroom model

Running cost of the representative agent:

$$
c(s,a,m)=c_0(s,a)+\kappa(s)\,m(s).
$$

Bellman optimality against frozen $m$:

$$
V(s)=\min_a\Bigl\{c(s,a,m)+\gamma\sum_{s'}P(s'\mid s,a)\,V(s')\Bigr\},
\qquad
\pi^\star(s)\in\arg\min_a Q_m(s,a).
$$

One-step mean-field (occupancy) update under a deterministic policy $\pi$:

$$
m'(s')=\sum_s m(s)\,P(s'\mid s,\pi(s)).
$$

Damped fixed-point iteration:

$$
m\leftarrow(1-\alpha)\,m+\alpha\,\mathrm{Evolve}^{k}(m,\pi^\star(m)),
$$

until the $L^1$ residual is small (consistency of $m$ with $\pi$).

## Contrast with Potential Games and $N$-player Nash (README only)

| Concept | Role | Notes |
| --- | --- | --- |
| **This package** (`Ada-Mean-Field-Game`) | Continuum representative agent vs $m$ | Consistency fixed point; not an $N\times N$ payoff tensor |
| **Potential Game** (sibling sheet) | Exact/weighted potential $P$ aligns incentives | Finite $N$; BR dynamics on potential; no mean field |
| **$N$-player Nash** (sibling sheet) | Mutual best responses in normal form | Exact finite-player NE; MFG is the $N\to\infty$ limit |

README links only — **no** package `with` of siblings. A potential game
organises finite-player best-reply structure via a scalar potential. Classical
Nash works with a fixed, usually small, player set. An MFG replaces pairwise
coupling by interaction through the distribution $m$ and asks for a policy
that regenerates that $m$.

## Classroom examples

### Congestion toy (3 states, 2 actions)

Locations on a cycle; Stay vs Move. Congestion weight $\kappa(s)>0$ makes
crowded sites expensive, so equilibrium mass spreads rather than collapsing
onto one location.

### Flocking toy (3 states, 3 actions)

Line $1$—$2$—$3$ with Left / Stay / Right. Negative $\kappa$ rewards co-location
(attraction / flocking) while movement carries friction.

### Entry toy (2 states, 2 actions)

Out / In with Stay / Switch. Being In is profitable until market congestion
on In bites — a sparse/entry MFG flavour.

### Two-state symmetric

Minimal symmetric congestion game for unit tests; equilibrium mass near
balanced under equal $\kappa$.

## Build

```bash
make        # gnatmake -gnatwa -gnat2022 -Pmean_field_game.gpr
make test   # run bin/tests
make clean
```

Requires GNAT with Ada 2022 support (`-gnat2022`). The project file
`mean_field_game.gpr` builds the standalone `tests` main into `bin/`.

## API summary

| Entity | Role |
| --- | --- |
| `Max_States`, `Max_Actions` | Caps ($8$, $4$) |
| `Distribution`, `Policy`, `Soft_Policy`, `Value_Vector` | Mean field, controls, values |
| `Transition_Tensor`, `Cost_Matrix`, `Congestion_Weights` | Dynamics and costs |
| `MFG` | Finite-state discounted MFG record |
| `Solve_Result` | Fixed-point output $(m,\pi,V,\mathrm{residual})$ |
| `Near`, `Default_Tol` | Numeric comparison |
| `Mass`, `Is_Simplex`, `Normalize`, `Uniform_Distribution`, `Dirac_Distribution` | Simplex helpers |
| `L1_Distance`, `Total_Variation` | Distribution metrics |
| `Running_Cost`, `Is_Valid_Transitions`, `Is_Stochastic_Row` | Cost / kernel checks |
| `Best_Response_Policy`, `Best_Response_Values`, `Soft_Best_Response`, `Value_Iteration`, `Q_Value` | BR / logit soft BR vs frozen $m$ |
| `Policy_Evaluation` | Evaluate fixed $\pi$ vs $m$ |
| `Evolve_Occupancy`, `Evolve_Occupancy_Soft`, `Evolve_Occupancy_N` | Occupancy / FP update |
| `Induced_Kernel`, `Consistency_Residual` | Kernel and residual |
| `Fixed_Point_Iterate`, `Solve` | Damped BR $\leftrightarrow$ $m$ iteration |
| `Congestion_Toy`, `Flocking_Toy`, `Entry_Toy`, `Two_State_Symmetric` | Classroom constructors |
| `Invalid_Argument` | Bad simplex / dims / discount / damping |

## License / series note

Educational reference code in the **RobertBoettcherSF** Ada 2023 algorithm
series. Not a continuous Lasry–Lions PDE solver, not a general
mean-field-type control (social planner) package, and not an $N$-player
exact Nash engine.
