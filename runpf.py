# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 15:52:16 2020

@author: ignac
"""
import pf as pf
import numpy as np
import loaddata as ld 


B,G = pf.createBG()
v,theta = pf.init_variables()

filename='data_3_bus.xlsx'
GENDATA,DEMDATA,LINDATA,STODATA = ld.loadsystemdata(filename)

# Pesp = GENDATA[:,1]
# Qesp = GENDATA[:,2]



# Pcalc = 3*[0]
# Qcalc = 3*[0]

Pesp = np.array([200,250])
Qesp = np.array([0])

Yesp = np.array([200,250,0])
Ycalc = np.array([0,0,0])

Xk = np.array([0,0,1])
Xk1 = np.array([0,0,0])

xP = [1,2]
xdelta = [2]

nodes_slack = [0]
nodes_PV = [1]
nodes_PQ = [2]

def updateV(v, nodes_slack,nodes_PV,nodes_PQ,Xk,xP,xdelta):
    numnodespv = len(xP)
    # numnodespq = len(xdelta)
    for i in nodes_slack:
        v[i] = 1
    for i in nodes_PV:
        v[i] = 1.02
    for i in range(len(nodes_PQ)):
        v[nodes_PQ[i]] = Xk[numnodespv+i]
    return v
    

def updateTheta(theta, nodes_slack,nodes_PV,nodes_PQ,Xk,xP,xdelta):
    # numnodespv = len(xP)
    numnodespq = len(xdelta)
    for i in nodes_slack:
        theta[i] = 0
    for i in range(len(nodes_PV)):
        theta[nodes_PV[i]] = Xk[i]
    for i in range(len(nodes_PQ)):
        theta[nodes_PQ[i]] = Xk[numnodespq+i]
    return theta


PQesp = np.array([200,250,0])
PQcalc = np.array([0,0,0])


J = np.ndarray((3,3))

MAXITER = 5
it=0
epsilon = 0.01
maxdif = 10000


while((it<MAXITER or maxdif<epsilon)):
    J = pf.updateJacobian(J, xP, xdelta, v,theta)
    dx = J.dot(Xk)
    # print('Iteración',it)
    # print('Pesp - Pcalc = ',dx)
    print('Pcalc = ',PQcalc)
    maxdif = max(abs(dx))
    PQcalc = PQesp - dx
    Xk1 = Xk - np.linalg.inv(J).dot(PQcalc)
    Xk = Xk1
    v = updateV(v,nodes_slack,nodes_PV,nodes_PQ,Xk,xP,xdelta)
    theta = updateTheta(theta,nodes_slack,nodes_PV,nodes_PQ,Xk,xP,xdelta)
    it=it+1
    
    
    