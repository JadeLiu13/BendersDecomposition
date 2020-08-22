# BendersDecomposition

This is a benders decomposition implementation in Python using Gurobi for solving the Uncapacitated Facility Location (UFL) Problem. The method solves two problems namely, masters problem and sub problem to get the UB and LB respectively (for the maximization problem). The sub problem provides cut to strengthen the lower bound in each iteration.

## Uncapacitated Facility Location Problem (UFL)
*Given*: n facilities and m customers, profit matrix of assigning a customer to facility and cost matrix of opening a new facility.

*Variables*: Which facilities to open (binary) which facility to assign to each customer (continuous).

*Constraints*: Assign one facility to each customer and assign a facility to a customer if and only if it is open.

*Objective*: Maximize the profit of assigning and minimize the cost of opening new facilities.

## Uncapacitated Facility Location Problem (UFL)
The script generates random data for profit matrix and cost vector. 

```
solveUFL(tol, x_initial, maxIter, verbose)
```
 - *tol* tolerance between the upper bound and lower bound
 
 - *x_initial* is the initial value of x variables (master problem variables)
 
 - *maxIter* Maximum no. of iterations allowed before the method fail to converge
 
 - *verbose* the value 1 prints the details of the convergence
 
