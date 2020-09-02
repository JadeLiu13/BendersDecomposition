'''
This script is written by Pramesh Kumar for testing various acceleration techniques for Benders Decomposition

Problem:

    
MAXIMIZE      \sum_{i \in C} \sum_{j \in F} p_{ij} y_{ij} - \sum_{j \in F} x_j
s.t.          \sum_{j \in F} y_{ij} = 1, \forall i \in C
              0 <= y_{ij} <= x_j, \forall i \in C, \forall j in F
              x binary, y binary
'''


import numpy as np
from gurobipy import *
import time






def generateFacilityLocationData(C, F):
    # Unbounded ray instance seed 159
    np.random.seed(15645)
    p =  np.random.randint(1000, size=(C, F))
    f = np.random.randint(1000, size=(F))
    for j in range(F):
        for i in range(C):
            f[j] += round(0.05*p[i,j])

    return C, F, p, f


# Step 1: Initialize variables
C = 100
F = 50



# Step 2: Start clock
ts = time.time()
# Step 3: Generate instance
C, F, p, f = generateFacilityLocationData(C,F)

############################################################################################################################
##############################################################################################################################

    


def solveModelGurobi():
    m2 = Model()
    x = {j: m2.addVar(lb=0, vtype=GRB.BINARY) for j in range(F)}
    y = {(i, j): m2.addVar(lb=0, vtype=GRB.BINARY) for i in range(C) for j in range(F)}
    for i in range(C):
        m2.addConstr(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            m2.addConstr(y[i, j] <= bigM*x[j])
    obj = 0
    for j in range(F):
        obj = obj -f[j] * x[j]
        for i in range(C):
            obj += p[i, j] * y[i, j]
    m2.setObjective(obj, sense=GRB.MAXIMIZE)
    m2.update()
    m2.Params.OutputFlag = 0
    m2.setParam('OutputFlag', False)
    m2.optimize()
    xVal= [x[j].x for j in range(F)]
    yVal =[y[i, j].x for i in range(C) for j in range(F)]
    return m2.objVal, xVal, yVal


def subProblem(x):
    m1 = Model()
    y = {(i, j):  m1.addVar(lb=0, vtype=GRB.CONTINUOUS) for i in range(C) for j in range(F)}
    constrMu = {}
    constrNu = {}
    for i in range(C):
        constrMu[i] = m1.addConstr(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            constrNu[i, j] = m1.addConstr(y[i, j] <= bigM*x[j])

    obj = sum([p[i, j] * y[i, j]  for j in range(F) for i in range(C)]) -sum([f[j] * x[j] for j in range(F)])

    m1.setObjective(obj, sense=GRB.MAXIMIZE)

    m1.update()
    m1.Params.OutputFlag = 0
    m1.Params.InfUnbdInfo = 1
    m1.Params.DualReductions = 0

    m1.optimize()

    if m1.status == GRB.OPTIMAL:
        mu = {}
        nu = {}
        for i in range(C):
            mu[i] = constrMu[i].pi
        for j in range(F):
            for i in range(C):
                nu[i, j] = constrNu[i, j].pi
        return m1.objVal, mu, nu, [y[i, j].x for i in range(C) for j in range(F)], m1.status
    else:
        mu = {}
        nu = {}
        for i in range(C):
            mu[i] = constrMu[i].FarkasDual
        for j in range(F):
            for i in range(C):
                nu[i, j] = constrNu[i, j].FarkasDual
        return -float("inf"), mu, nu, [], m1.status
    
    

def setupMasterProblemModel():
    m = Model()
    eta = {i: m.addVar(vtype=GRB.CONTINUOUS, ub = 1e5, name ='eta_' + str(i)) for i in range(C)}  
    x = {j: m.addVar(lb=0, vtype=GRB.BINARY, name = str(j)) for j in range(F)}
    
    m.setObjective(sum([eta[i] for i in range(C)]) -sum([f[j] * x[j] for j in range(F)]), sense=GRB.MAXIMIZE)
    m.update()
    m.Params.OutputFlag = 0
    m.Params.lazyConstraints = 1
    
    # Limit how many solutions to collect
    m.setParam(GRB.Param.PoolSolutions, 500)
    # Limit the search space by setting a gap for the worst possible solution
    # that will be accepted
    m.setParam(GRB.Param.PoolGap, 0.8)
    # do a systematic search for the k-best solutions
    m.setParam(GRB.Param.PoolSearchMode, 0)

    return m

def solveMaster(m,  optCuts_mu, optCuts_nu, fesCuts_mu, fesCuts_nu):

    '''
    Adding optimality cut
    '''
    for mu in range(len(optCuts_nu)):            
        if len(optCuts_mu[mu]) != 0:   
            for i in range(C):
                tot = optCuts_mu[mu][i]
                #tot = sum(optCuts_mu[mu].values())    
                for j in range(F):
                    tot += optCuts_nu[mu][i, j]*m.getVarByName(str(j))*bigM        
                m.addConstr(m.getVarByName('eta_' + str(i)) <= tot)
        '''
    Adding feasibility cut
    '''
    for mu in range(len(fesCuts_mu)):    
        if len(fesCuts_mu[mu]) != 0:     
            for i in range(C):
                tot = fesCuts_mu[mu][i]
                for j in range(F):
                    tot += fesCuts_nu[mu][i, j] * m.getVarByName(str(j))* bigM
                m.addConstr (tot >= 0)
    

    
    m.optimize()
      
    
    
    if m.status == GRB.OPTIMAL:
        nSolutions = m.SolCount     # Print number of solutions stored
        optSol = [round(m.getVarByName(str(k)).x) for k in range(F)]
        otherSols = [optSol]
        if (nSolutions >= 2):
            for solNum in range(nSolutions - 1):                    
                m.setParam(GRB.Param.SolutionNumber, solNum+1)
                sol = [round(m.getVarByName(str(k)).Xn) for k in range(F)]
                if sol not in otherSols:
                    otherSols.append(sol)   
        return m.objVal, optSol, m.getVarByName('eta'), m, otherSols
    else:
        print("Sth went wrong in the master problem and it is ", m.status)





def solveUFLBendersMultipleDisaggCuts(eps, x_initial, maxit, verbose=0):
    UB = float("inf")
    LB = -float("inf")
    optCuts_mu = []
    optCuts_nu = []
    fesCuts_mu = []
    fesCuts_nu = []
    tol = float("inf")
    it = 0
    x = x_initial
    m = setupMasterProblemModel()
    solutions = [x]
    solutionsEverFound = [list(x)]
    while eps <= tol  and it < maxit :
        muCutsOpt = []
        nuCutsOpt = []
        muCutsFeas = []
        nuCutsFeas = []
        for s in solutions:
            ob, mu, nu, y, status = subProblem(s)
            LB = max(LB, ob)
            if status == 2:                    
                muCutsOpt.append(mu)
                nuCutsOpt.append(nu)
            else:
                muCutsFeas.append(mu)
                nuCutsFeas.append(nu)
                
        obj, x, eta, m, otherSol = solveMaster(m, muCutsOpt, nuCutsOpt, muCutsFeas, nuCutsFeas)
        UB = min(UB, obj)
        solutions = []
        for sol in otherSol:
            if sol not in solutionsEverFound:
                solutions.append(sol)
                solutionsEverFound.append(sol)     
                
        
                            
        print(len(solutions))
        
        tol = UB - LB
        it += 1
        if verbose == 1:            
            print('----------------------iteration '  + str(it) +'-------------------' )
            print ('LB = ', LB, ', UB = ', UB, ', tol = ', tol)
            if len([k for k in range(F) if round(x[k]) != 0]) != 0:
                print('Opened Facilities: \t ', [k for k in range(F) if round(x[k]) != 0])
            else:
                print('No open facilities')
            
            '''
            print('Assignment....')
            for i in range(C):
                print('Customer \t | \t Facility')
                print(str(i)+ '\t | \t ' + str([k for k in for j in range(F) for i in range(C) if y[0] == i and y[k] == 1][0]))
            print(LB, UB, tol, x, it)
            '''
        else:
            continue
    return x, y, UB


def checkGurobiBendersSimilarity(xb, yb, xg, yg):

    ind = 0
    for j in range(F):
        if xb[j] != xg[j]:
            ind += 1
    k = 0
    for j in range(F):
        for i in range(C):
            if yb[k] != yg[k]:
                ind += 1
            k += 1
    if ind == 0:
        print('Solution obtained from both methods are same')
    else:
        print('Solution obtained from both methods are different!!')


bigM = 1
x_initial = np.zeros(F)
x_initial[1] = 1
x_initial[2] = 0
start = time.time()
xb, yb, obb = solveUFLBendersMultipleDisaggCuts(1, x_initial, 1000, 0)
print("Benders took...", round(time.time() - start, 2), "seconds")
start = time.time()
obg, xg, yg = solveModelGurobi()
print("Gurobi took...", round(time.time() - start, 2), "seconds")
checkGurobiBendersSimilarity(xb, yb, xg, yg)
