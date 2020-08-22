# BendersDecomposition

This is a benders decomposition implementation in Python using Gurobi for solving the Uncapacitated Facility Location (UFL) Problem.

## Uncapacitated Facility Location Problem (UFL)
Given: n facilities and m customers, profit matrix of assigning a customer to facility and cost matrix of opening a new facility.

Variables: Which facilities to open (binary) which facility to assign to each customer (continuous)

Constraints: Assign one facility to each customer and assign a facility to customer if and only if it is open



