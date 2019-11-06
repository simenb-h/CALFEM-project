# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""
import numpy as np
import calfem.core as cfc
import quads as quads
import calfem.vis as cfv

nElx = 40
nEly = 20

bDrawMesh = True

# Cantilever with dimensions H x L x thickness
H         = 2.5
L         = 12.5
thickness = 0.1

# Distributed load in x and y
eq = np.array([3.0,0.])
eqTotal = eq * L * H * thickness

# Material properties and thickness

ep = [1,thickness]
E  = 2.1e11
nu = 0.3
Dmat = np.mat([
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

nEls = nElx * nEly
nNodx = nElx +1
nNody = nEly +1
nNods = nNodx * nNody

L_elx = L / nElx
L_ely = H / nEly

nelnod = 4

coords = np.zeros((nNods,2))
dofs   = np.zeros((nNods,2),int)       #Dofs is starting on 1 on first dof
edofs  = np.zeros((nEls,nelnod*2),int) #edofs are also starting on 1 based dof


inod = 0 # The coords table starts numbering on 0
idof = 1 # The values of the dofs start on 1
ndofs = nNods *2

# Set the node coordinates and node dofs

for i in range(nNodx):
    for j in range(nNody):
        coords[inod,0] = L_elx * i
        coords[inod,1] = L_ely * j
        dofs[inod,0] = idof
        dofs[inod,1] = idof+1
        idof += 2
        inod += 1

# Set the element connectivites and element dofs

elnods = np.zeros((nEls,nelnod),int)
eldofs = np.zeros((nEls,nelnod*2),int)

iel = 0
for i in range(nElx):
    for j in range(nEly):
        # 0 based node numbers
        i1 =     i*nNody + j
        i2 = (i+1)*nNody + j
        i3 = (i+1)*nNody + (j+1)
        i4 =     i*nNody + (j+1)
        # 1 based dof numbers
        idof1 = i1*2 +1  # 1 based dof number for the element
        idof2 = i2*2 +1
        idof3 = i3*2 +1
        idof4 = i4*2 +1

        elnods[iel,:] = [i1, i2, i3, i4]
        eldofs[iel,:] = [idof1,idof1+1, idof2,idof2+1, idof3,idof3+1, idof4,idof4+1]
        iel += 1


# Draw the mesh.
if bDrawMesh:
    cfv.drawMesh(
        coords=coords,
        edof=eldofs,
        dofsPerNode=2,
        elType=3,
        filled=True,
        title="4 node quad elements")
    cfv.showAndWait()

# Set fixed boundary condition on left side, i.e. nodes 0-nNody
bc = np.array(np.zeros(nNody*2),'i')
idof = 1
for i in range(nNody):
    idx = i*2
    bc[idx]   = idof
    bc[idx+1] = idof+1
    idof += 2

# Assemble stiffness matrix

# Extract element coordinates
ex, ey = cfc.coordxtr(eldofs,coords,dofs)

K = np.zeros((ndofs,ndofs))
R = np.zeros((ndofs,1))
#r = np.zeros( ndofs )

for iel in range(nEls):
    K_el, f_el = quads.quad4e(ex[iel],ey[iel],Dmat,thickness,eq)
    cfc.assem(eldofs[iel],K,K_el,R,f_el)

r, R0 = cfc.solveq(K,R,bc)

xTR = r[-2,0]
yTR = r[-1,0]
xBR = r[-(nNody)*2  ,0]
yBR = r[-(nNody)*2+1,0]
print("Displacement upper-right, x:{:12.3e}       y:{:12.3e}".format(xTR, yTR))
print("Displacement lower-right, x:{:12.3e}       y:{:12.3e}".format(xBR, yBR))

R0Top = R0[nNody*2 -2,0]
R0Bot = R0[0,0]
print("\nReaction forces in X,   top:{:12.3e}  bottom:{:12.3e}".format(R0Top, R0Bot))

RySum = 0
RxSum = 0
for i in range(0,(nNody*2),2):
    RxSum += R0[i,0]
    RySum += R0[i+1,0]

print("Total reaction force in x:{:12.3e}   y:{:12.3e} (Total x{:12.3e}  y:{:12.3e})".format(RxSum,RySum,eqTotal[0],eqTotal[1]))

# Draw the displacements

if bDrawMesh:
    disp = np.array(np.zeros((nNods, 2)), 'f')
    rMax = max(abs(max(r)), abs(min(r)))
    scale = 0.15 * L / rMax

    for i in range(np.size(disp, 0)):
        disp[i, 0] = r[i * 2, 0] * scale
        disp[i, 1] = r[i * 2 + 1, 0] * scale

    cfv.drawDisplacements(displacements=disp,
        coords=coords,
        edof=eldofs,
        dofsPerNode=2,
        elType=3,
        title="4 node quad elements")
    cfv.showAndWait()