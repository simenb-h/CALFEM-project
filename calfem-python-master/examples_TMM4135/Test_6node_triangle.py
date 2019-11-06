# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""

# example exs5
# ----------------------------------------------------------------
# PURPOSE 
#    Analysis of a simply supported beam.
# ----------------------------------------------------------------

# REFERENCES
#     G"oran Sandberg 94-03-08 
#     Karl-Gunnar Olsson 95-09-28
#     Ola Dahlblom 2004-09-21
# ----------------------------------------------------------------

import numpy as np
import calfem.core as cfc
import triangles as tri

cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

# ----- Topology -------------------------------------------------
ex = np.array([0.0,1.0,0.0,0,0,0])
ey = np.array([0.0,0.0,1.0,0,0,0])

for i in range(3):
    j = cyclic_ijk[i+1]
    k = cyclic_ijk[i+2]
    ex[i+3] = (ex[i] + ex[j])/2
    ey[i+3] = (ey[i] + ey[j])/2

th = 0.1
ep = [1,th]

E  = 2.1e11
nu = 0.3

D = np.mat([
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

eq = [1.0, 3.0]

#Ke, fe = tri.plante(ex,ey,ep,D,eq)

Ke = np.mat(np.zeros((12,12)))
fe = np.mat(np.zeros((12,1)))



rigX = np.mat(np.zeros((12,1)))
rigY = np.mat(np.zeros((12,1)))
rigR = np.mat(np.zeros((12,1)))



for i in range(6):
    rigX[i*2  ,0] = 1.0
    rigY[i*2+1,0] = 1.0
    rigR[i*2  ,0] = ey[i]
    rigR[i*2+1,0] = -ex[i]

Ke, fe = tri.tri6e(ex,ey,D,th,eq)

print('Stiffness matrix:\n', Ke)
print('Consistent forces:\n', fe)

fx = Ke * rigX
fy = Ke * rigY
fr = Ke * rigR

print('Force from rigX translation:\n',fx)
print('Force from rigY translation:\n',fy)
print('Force from rigR rotation:\n',fr)


constEx = np.mat(np.zeros((12,1)))
constEy = np.mat(np.zeros((12,1)))
constGamma1 = np.mat(np.zeros((12,1)))
constGamma2 = np.mat(np.zeros((12,1)))


for i in range(6):
    constEx[i*2  ,0] = ex[i]
    constEy[i*2+1,0] = ey[i]
    constGamma1[i*2  ,0] = ey[i]
    constGamma2[i*2+1,0] = ex[i]
    
zetaInt = np.array([[0.5,0.5,0.0],
                    [0.0,0.5,0.5],
                    [0.5,0.0,0.5]]) 
    
for i in range(3):
    Be = tri.tri6_Bmatrix(zetaInt[i],ex,ey)
    
    print('Be\n', Be)
    
    Ex = Be * constEx
    Ey = Be * constEy
    G1  = Be * constGamma1
    G2  = Be * constGamma2
    print('Ex:\n',Ex)
    print('Ey:\n',Ey)
    print('G:\n',G1)
    print('G:\n',G2)

    