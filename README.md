# BendersDecomposition

This is a benders decomposition implementation in Python using Gurobi for solving the Uncapacitated Facility Location (UFL) Problem. The classic Benders partitioning method solves two problems namely, master problem and sub problem to get the UB and LB respectively (for the maximization problem). The sub problem provides cuts to strengthen the lower bound in each iteration. To read more about this method, please refer to [Prof. Rubin's blog](https://orinanobworld.blogspot.com/2011/10/benders-decomposition-then-and-now.html).

## Uncapacitated Facility Location Problem (UFL)
*Given*: n facilities and m customers, profit matrix of assigning a customer to facility and cost matrix of opening a new facility.

*Variables*: Which facilities to open (x binary) and which facility to assign to each customer (y continuous).

*Constraints*: Assign one facility to each customer and assign a facility to a customer if and only if it is open.

*Objective*: Maximize the profit of assigning and minimize the cost of opening new facilities.
```
MAXIMIZE      \sum_{i \in C} \sum_{j \in F} p_{ij} y_{ij} - \sum_{j \in F} c_j x_j
s.t.          \sum_{j \in F} y_{ij} = 1, \forall i \in C
              0 <= y_{ij} <= x_j, \forall i \in C, \forall j in F
              x binary, y continuous
```
## Details about various Scripts
### 1. bendersClassic.py
The script generates random data for profit matrix and cost vector and program the classic Benders Decomposition. The explaination of some of the methods is given below:

```
1. solveUFLBenders(tol, x_initial, maxIter, verbose)
```
for solving the problem using Benders Decomposition. It returns the optimal value of the decision variables as well as the optimal objective value.

 - *tol* tolerance between the upper bound and lower bound
 
 - *x_initial* is the initial value of x variables (master problem variables)
 
 - *maxIter* Maximum no. of iterations allowed before the method fail to converge
 
 - *verbose* the value 1 prints the details of the convergence
 
 
```
2. solveModelGurobi()
```
 You can also solve the model using Gurobi using this function. It returns the optimal value of the decision variables as well as the optimal objective value.
 
 ### 2. bendersClassicWithCallBackLazy.py
 This script adds another function over on top of functions in bendersClassic.py. This new function is called 'runCallBackBenders'. The classic Benders Decomposition requires solving the Master Problem from scratch every time a new cut is added. The callback procedure adds a lazy cut and continue the branch and bound tree. Once a feasible solution is found then the bounds are checked. To read more about this procedure, please refer to [Prof. Rubin's blog](https://orinanobworld.blogspot.com/2011/10/benders-decomposition-then-and-now.html). 
 
 
```
1. runCallBackBenders()  
```
 
 It also returns the optimal value of the decision variables as well as the optimal objective value. 
 Although the callback procedure performs better than the classic approach, you will observe that both these methods cannot beat the Gurobi Solver in CPU time.
 
 ### 3. bendersDisaggregatedCuts.py
 By fixing the solution of the variables related to which fcacilities to open, we obtain the Benders subproblem. This benders subproblem can be further decomposed for each customer. Therefore, multiple disaggregated cuts (for each customer) can be added to the master problem after solving the subproblem. This script codes above two classic  procedures with disaggregated cuts.  Read more [here](https://orinanobworld.blogspot.com/2012/09/separable-benders-decomposition.html)
 
 
 ```
 1. bendersDisaggCuts(tol, x_initial, maxIter, verbose)
 2. runCallBackBendersDisagg()
```
 These methods computes the solution in lower computational time than the Gurobi solver for this problem. 
