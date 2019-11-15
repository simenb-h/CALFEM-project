# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""
import numpy as np

def plante(ex,ey,ep,D,eq=None):
    
    Dshape = D.shape
    if Dshape[0] != 3:
        raise NameError('Wrong constitutive dimension in plante')
        
    if ep[0] == 1 :
        return tri3e(ex,ey,D,ep[1],eq)
    else:
        Dinv = np.inv(D)
        return tri3e(ex,ey,Dinv,ep[1],eq)


def tri3e(ex,ey,D,th,eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """
    
    tmp = np.matrix([[1,ex[0],ey[0]],
                     [1,ex[1],ey[1]],
                     [1,ex[2],ey[2]]])
    
    A2 = np.linalg.det(tmp)  # Double of triangle area
    A  = A2 / 2.0
       
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

    # TODO: fill out missing parts (or reformulate completely)

    ''' MY CODE '''
 
    zeta_px, zeta_py = zeta_partials_x_and_y(ex,ey)
    
    B  = (1/A2)*np.array([[zeta_px[0], 0, zeta_px[1], 0, zeta_px[2], 0],
                            [0, zeta_py[0], 0, zeta_py[1], 0, zeta_py[2]],
                            [zeta_py[0], zeta_px[0], zeta_py[1], zeta_px[1], zeta_py[2], zeta_px[2]]])
    

    Ke = np.mat(np.zeros((6, 6)))

    Ke = (B.T * D * B) * A * th


    if eq is None:
        return Ke
    else:
        fe = np.mat(np.zeros((6,1)))
        return Ke, fe
    
def zeta_partials_x_and_y(ex,ey):
    """
    Compute partials of area coordinates with respect to x and y.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    """
    
    tmp = np.matrix([[1,ex[0],ey[0]],
                     [1,ex[1],ey[1]],
                     [1,ex[2],ey[2]]])
    
    A2 = np.linalg.det(tmp)  # Double of triangle area
       
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k
    
    zeta_px = np.zeros(3)           # Partial derivative with respect to x
    zeta_py = np.zeros(3)           # Partial derivative with respect to y

    # TODO: fill out missing parts (or reformulate completely)
   
    a_i = []
    v_i = []
    c_i = []
    for i in range(0,2): #Hvilken range? 
        a_i.append(ex[cyclic_ijk[i+1]]*ey[cyclic_ijk[i+2]] - ex[cyclic_ijk[i+2]]*ey[cyclic_ijk[i+1]])
        v_i.append(ey[cyclic_ijk[i+1]]-ey[cyclic_ijk[i+2]])
        c_i.append(ex[cyclic_ijk[i+2]]-ex[cyclic_ijk[i+1]])
        #A_i?

    for i in range(0,2):
        zeta_px[i] = v_i[i]/A2
        zeta_py[i] = c_i[i]/A2
     
    return zeta_px, zeta_py


# Functions for 6 node triangle
    
def tri6_area(ex,ey):
        
    tmp = np.matrix([[1,ex[0],ey[0]],
                     [1,ex[1],ey[1]],
                     [1,ex[2],ey[2]]])
    
    A = np.linalg.det(tmp) / 2
    
    return A


def tri6_shape_functions(zeta):
    
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

    N6 = np.zeros(6)

    # TODO: fill out missing parts (or reformulate completely) 
    '''DONE?'''

    for i in range(0,2):
       z = zeta[i]
       N6.append(z*2(z-1))
    
    for i in range(0,2):
        j = cyclic_ijk[i+1]
        N6.append(4*zeta[i]*zeta[j])


    return N6


def tri6_shape_function_partials_x_and_y(zeta,ex,ey):
    
    zeta_px, zeta_py = zeta_partials_x_and_y(ex,ey)
    
    N6_px = np.zeros(6)
    N6_py = np.zeros(6)
    
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

    # TODO: fill out missing parts (or reformulate completely)
    ''' DONE? '''
    N6_px = tri6_shape_functions(zeta_px)
    N6_py = tri6_shape_functions(zeta_py)


    return N6_px, N6_py


def tri6_Bmatrix(zeta,ex,ey):
    
    nx,ny = tri6_shape_function_partials_x_and_y(zeta, ex, ey)

    Bmatrix = np.matrix(np.zeros((3,12)))

    # TODO: fill out missing parts (or reformulate completely)
    ''' DONE ? '''
    A = tri6_area(ex,ey)
    iA2 = 1/(2*A) #i before A2 because it's the inverse of A2

    B  = (iA2)*np.array([[nx[0], 0, nx[1], 0, nx[2], 0],
                  [0, ny[0], 0, ny[1], 0, ny[2]],
                  [ny[0],nx[0],ny[1],nx[1],ny[2],nx[2]]])

    return Bmatrix


def tri6_Kmatrix(ex,ey,D,th,eq=None):
    
    zetaInt = np.array([[0.5,0.5,0.0],
                        [0.0,0.5,0.5],
                        [0.5,0.0,0.5]])
    
    '''HVA ER DENNE VERDIEN OG SKAL DEN INN I B? '''
    wInt = np.array([1.0/3.0,1.0/3.0,1.0/3.0])

    A    = tri6_area(ex,ey)
    
    Ke = np.matrix(np.zeros((12,12)))

    # TODO: fill out missing parts (or reformulate completely)
    B = tri6_Bmatrix(zetaInt, ex, ey)


    Ke = (B.T * D * B)* A * th

    if eq is None:
        return Ke
    else:
        fe = np.matrix(np.zeros((12,1)))
        ''' HVA SKJER OM DET ER EN FORDELT LAST PÅ? ''' 
        # TODO: fill out missing parts (or reformulate completely)

        return Ke, fe

def tri6e(ex,ey,D,th,eq=None):
    return tri6_Kmatrix(ex,ey,D,th,eq)

