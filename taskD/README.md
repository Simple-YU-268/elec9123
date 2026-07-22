# ELEC9123 Design Task D - Discrete Event Simulation of Queues

Author: z5615351 YU, DTB  
Course: ELEC9123 Wireless Communications  
Year: 2026

## Folder Structure

This submission follows the Tiered Learning Taxonomy (TLT) Level 5 structure. Each folder contains a Main, Analysis, and Simulation `.m` file triplet.

| Folder    | System                         | TLT Level |
|-----------|--------------------------------|-----------|
| M_M_1     | M/M/1 single-server queue      | Level 2   |
| M_M_1_K   | M/M/1/K finite-capacity queue  | Level 3   |
| M_G_1     | M/G/1 with general service     | Level 4   |
| M_M_m     | Two-class M/M/m airport model  | Level 5   |

## Files per Folder

- `*_main.m`        - driver script, plots theory vs simulation
- `*_analysis.m`    - closed-form theoretical metrics
- `*_simulation.m`  - hand-written discrete-event simulation engine

## Running the Code

Open MATLAB in the `taskD` folder and run the top-level driver:

```matlab
z5615351_YU_DTB_DTD_2026
```

To run individual models, uncomment the corresponding section in the driver or run the `*_main.m` file directly from its folder.

## Parameters Summary

- M/M/1: λ = 10/h, service = 5 min
- M/M/1/K: K = 10, λ = 10/h, service = 5 and 6 min
- M/G/1: λ = 10/h, mean service = 5 min (exponential, uniform 2.5-7.5, deterministic)
- M/M/m: Business λ = 20/h, m = 2, s = 2.5 min; Economy λ = 100/h, m = 5, s = 2.5 min

## AI Tool Declaration

AI tools were used during the design and coding of this task, including code generation and explanatory text drafting. The final verification and design journal interpretation remain the author's own work.
