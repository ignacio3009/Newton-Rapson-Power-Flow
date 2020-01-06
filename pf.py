# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 14:55:03 2020

@author: ignac
"""
import loaddata as ld
import math
filename='data_3_bus.xlsx'
GENDATA,DEMDATA,LINDATA,STODATA = ld.loadsystemdata(filename)

X = LINDATA[:,4]
R = LINDATA[:,5]
BO = LINDATA[:,2]
BD = LINDATA[:,3]
LARGE = LINDATA[:,6]


"""
x= [delta, \v\]

"""
B = [[0,0,0],[0,0,0],[0,0,0]]
G = [[0.5,0.5,0.5],[0.5,0.5,0.5],[0.5,0.5,0.5]]
numNodes = 3

def init_variables():

    v = numNodes*[1]
    theta = numNodes*[0]
    return v, theta

v, theta = init_variables()
# v_k_1, theta_k_1 = init_variables()


npq=1 #not known: V, theta
npv=1 #not known: theta 
m = 2*npq + npv  
x = m*[0] #delta, V

ntheta = npq+npv
nv = npq


for i in range(ntheta): x[i] = 0
for i in range(ntheta,nv+ntheta): x[i] = 1

J=[[0 for i in range(m)] for j in range(m)]
nodes = [0,1,2]

Pesp = [10,10,10]
Pcalc = [0,0,0]
Qesp = []

vb=220
sb=100

z_base = vb**2/sb


def createBG():
    for i in range(len(LINDATA)):
        ni = int(BO[i])
        nj = int(BD[i])
        large_km = LARGE[i]
        x_ohm_km = X[i]
        x_ohm = x_ohm_km*large_km
        x_pu = x_ohm/z_base
        
        r_ohm_km = R[i]
        r_ohm = r_ohm_km*large_km
        r_pu = r_ohm/z_base
        
        den = (x_pu**2 + r_pu**2)
        B[ni][nj] = x_pu/den
        B[nj][ni] = x_pu/den
        G[ni][nj] = r_pu/den
        G[nj][ni] = r_pu/den
        
    for i in range(len(B)):
        B[i][i] = -1*sum(B[i][:])
        G[i][i] = sum(G[i][:])
    return B,G
        
        

def dP_ddelta(i,j,V,theta):
    if i==j:
        s=0
        for k in range(numNodes):
            if k!=i:
                deltatheta =  theta[i] - theta[k]
                cos_deltatheta = math.cos(deltatheta)
                sin_deltatheta = math.sin(deltatheta) 
                Bik = B[i][k]
                Gik = G[i][k]
                Vik = v[i]*v[k]
                s = s + Vik*(Bik*cos_deltatheta - Gik*sin_deltatheta)
        return s
    else:
        Vij = v[i]*v[j]
        Bij = B[i][j]
        Gij = G[i][j]
        deltatheta =  theta[i] - theta[j]
        cos_deltatheta = math.cos(deltatheta)
        sin_deltatheta = math.sin(deltatheta) 
        return Vij*(Gij*sin_deltatheta - Bij*cos_deltatheta)
    

def dQ_ddelta(i,j,V,theta):
    if i==j:
        s=0
        for k in range(numNodes):
            if k!=i:
                deltatheta =  theta[i] - theta[k]
                cos_deltatheta = math.cos(deltatheta)
                sin_deltatheta = math.sin(deltatheta) 
                Bik = B[i][k]
                Gik = G[i][k]
                Vik = v[i]*v[k]
                s = s + Vik*(Gik*cos_deltatheta + Bik*sin_deltatheta)
        return s
    else:
        Vij = v[i]*v[j]
        Bij = B[i][j]
        Gij = G[i][j]
        deltatheta =  theta[i] - theta[j]
        cos_deltatheta = math.cos(deltatheta)
        sin_deltatheta = math.sin(deltatheta) 
        return -1*Vij*(Gij*cos_deltatheta - Bij*sin_deltatheta)
    
 
def dP_dV(i,j,V,theta):
    s=0
    if i==j:
        for k in range(numNodes):
            if k!=i:
                deltatheta =  theta[i] - theta[k]
                cos_deltatheta = math.cos(deltatheta)
                sin_deltatheta = math.sin(deltatheta) 
                Bik = B[i][k]
                Gik = G[i][k]
                Vik = v[i]*v[k]
                s = s + Vik*(Gik*cos_deltatheta + Bik*sin_deltatheta) 
        Gii = G[i][i]
        Vi2 = v[i]**2
        return s + 2*Gii*Vi2
    else:
        Vij = v[i]*v[j]
        Bij = B[i][j]
        Gij = G[i][j]
        deltatheta =  theta[i] - theta[j]
        cos_deltatheta = math.cos(deltatheta)
        sin_deltatheta = math.sin(deltatheta) 
        return Vij*(Gij*cos_deltatheta + Bij*sin_deltatheta)
    
    
    
            

def dQ_dV(i,j,V,theta):
    if i==j:
        s=0
        for k in range(numNodes):
            if k!=i:
                deltatheta =  theta[i] - theta[k]
                cos_deltatheta = math.cos(deltatheta)
                sin_deltatheta = math.sin(deltatheta) 
                Bik = B[i][k]
                Gik = G[i][k]
                Vik = v[i]*v[k]
                s = s + Vik*(Gik*sin_deltatheta - Bik*cos_deltatheta)
        Bii = B[i][i]
        Vi2 = v[i]**2
        return s - 2*Bii*Vi2
    else:
        Vij = v[i]*v[j]
        Bij = B[i][j]
        Gij = G[i][j]
        deltatheta =  theta[i] - theta[j]
        cos_deltatheta = math.cos(deltatheta)
        sin_deltatheta = math.sin(deltatheta) 
        return Vij*(Gij*sin_deltatheta - Bij*cos_deltatheta)
    
    
def updateJacobian(J,xP,xdelta, v, theta):
    nxp = len(xP)
    nxd = len(xdelta)
    numP = range(nxp)
    numQ = range(nxd)
    numPQ = range(nxp,nxp+nxd)
    
    for i in numP:
        for j in numP:
            J[i][j] = dP_ddelta(i,j,v,theta)
        for j in numPQ:
            J[i][j] = dP_dV(i,j,v,theta)
    for i in numPQ:
        for j in numP:
            J[i][j] = dQ_ddelta(i,j,v,theta)
        for j in numPQ:
            J[i][j] = dQ_dV(i,j,v,theta)
    return J


