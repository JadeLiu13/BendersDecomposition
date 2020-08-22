# BendersDecomposition

This is a benders decomposition implementation in Python using Gurobi for solving the Uncapacitated Facility Location (UFL) Problem. The classic Benders partitioning method solves two problems namely, masters problem and sub problem to get the UB and LB respectively (for the maximization problem). The sub problem provides cut to strengthen the lower bound in each iteration. To read more about this method, please refer to Prof. Rubin's [blog]{https://orinanobworld.blogspot.com/2011/10/benders-decomposition-then-and-now.html}.

## Uncapacitated Facility Location Problem (UFL)
*Given*: n facilities and m customers, profit matrix of assigning a customer to facility and cost matrix of opening a new facility.

*Variables*: Which facilities to open (x binary) and which facility to assign to each customer (y continuous).

*Constraints*: Assign one facility to each customer and assign a facility to a customer if and only if it is open.

*Objective*: Maximize the profit of assigning and minimize the cost of opening new facilities.
```
MAXIMIZE      \sum_{i \in C} \sum_{j \in F} p_{ij} y_{ij} - \sum_{j \in F} x_j
s.t.          \sum_{j \in F} y_{ij} = 1, \forall i \in C
              0 <= y_{ij} <= x_j, \forall i \in C, \forall j in F
              x binary, y continuous
```
## Information about Different Scripts
### bendersClassic.py
The script generates random data for profit matrix and cost vector and program the classic Benders Decomposition.

```
solveUFL(tol, x_initial, maxIter, verbose)
```
 - *tol* tolerance between the upper bound and lower bound
 
 - *x_initial* is the initial value of x variables (master problem variables)
 
 - *maxIter* Maximum no. of iterations allowed before the method fail to converge
 
 - *verbose* the value 1 prints the details of the convergence
 
