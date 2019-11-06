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

# ----- Topology -------------------------------------------------
ex = np.array([0.,1.,0.])
ey = np.array([0.,0.,1.])

th = 0.1
ep = [1,th]

E  = 2.1e11
nu = 0.3

D = np.mat([
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

eq = [1.0, 3.0]

Ke, fe = tri.plante(ex,ey,ep,D,eq)

print('Stiffness matrix:\n', Ke)
print('Consistent forces:\n', fe)

rigX = np.mat([[1.0],[0.0],[1.0],[0.0],[1.0],[0.0]])
rigY = np.mat([[0.0],[1.0],[0.0],[1.0],[0.0],[1.0]])
rigR = np.mat([[-ey[0]],[ex[0]],[-ey[1]],[ex[1]],[-ey[2]],[ex[2]]])

fx = Ke * rigX
fy = Ke * rigY
fr = Ke * rigR

print('Force from rigX translation:\n',fx)
print('Force from rigY translation:\n',fy)
print('Force from rigR rotation:\n',fr)
