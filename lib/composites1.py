# coding: utf-8

'''
Module for composite material analysis

Hyer-Stress Analysis of Fiber-Reinforced Composite Materials
Herakovich-Mechanics of Fibrous Composites
Daniel-Engineering Mechanics of Composite Materials
Kollar-Mechanics of COmposite Structures
NASA- Basic Mechancis of Lamianted Composites
    https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950009349.pdf
    
TODO:
    
 * transverse shear stress reddy pg 136 or daniel pg 139
 * include line loads (Qx,Qy) for combined loading 
 * calculate capability of panel based on margin
 
'''
#==============================================================================
# Import Modules
#==============================================================================
from __future__ import print_function, division

__author__ = 'Neal Gordon <nealagordon@gmail.com>'
__date__ =   '2016-12-02'
__version__ = 0.1

from copy import copy
from numpy import pi, zeros, ones, linspace, arange, array, sin, cos, sqrt, pi
from numpy.linalg import solve, inv

#from scipy import linalg
import numpy as np
#np.set_printoptions(suppress=False,precision=2)   # suppress scientific notation
np.set_printoptions(precision=3, linewidth=200)#, threshold=np.inf)

import scipy
from scipy.spatial import ConvexHull
#np.set_printoptions(formatter={'float': lambda x: "{:.2f}".format(x)})

import pandas as pd

import sympy as sp
from sympy import Function, dsolve, Eq, Derivative, symbols, pprint
from sympy.plotting import plot3d

#from sympy import cos, sin
#sp.init_printing(use_latex='mathjax')
#sp.init_printing(wrap_line=False, pretty_print=True)
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = (8,5)
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 14
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,figure,xlim,ylim,title,legend, \
grid, show, xlabel,ylabel, tight_layout
from mpl_toolkits.mplot3d import axes3d
    
# if using ipython console, turn off inline plotting
#mpl.use('Qt5Agg')
# inline plotting
from IPython import get_ipython
#get_ipython().magic('matplotlib inline')

###disable inline plotting
try:
    get_ipython().magic('matplotlib')
except:
    pass

from IPython.display import display

import os

plt.close('all')
#==============================================================================
# Functions
#==============================================================================


def import_matprops(mymaterial=['T300_5208','AL_7075']):
    '''
    import material properties
    '''

    matprops = pd.read_csv(os.path.join(os.path.dirname(__file__), "compositematerials.csv"), index_col=0)
    

    if mymaterial==[] or mymaterial=='':
        print(matprops.columns.tolist())

    mat = matprops[mymaterial]
    #mat.applymap(lambda x:np.float(x))
    mat = mat.applymap(lambda x:pd.to_numeric(x, errors='ignore'))
    return mat


def Sf(E1,E2,nu12,G12):
    '''transversely isptropic compliance matrix. pg 58 herakovich'''
    nu21 = E2*nu12/E1
    S = array([[1/E1,    -nu21/E2, 0],
               [-nu12/E1, 1/E2,    0],
               [0,        0,       1/G12]])
    return S

def S6f(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23):
    '''
    daniel pg 74
    transversely isotropic compliance matrix.
    For transversly isotropic
    E2=E3, nu12=nu13,G12=G13,G23=E2/(2(1+nu23))
    '''
    S6 = array( [[    1/E1, -nu12/E1, -nu12/E1,     0,     0,       0],
                 [-nu12/E1,     1/E2, -nu23/E2,     0,     0,       0],
                 [-nu12/E1, -nu23/E2,     1/E2,     0,     0,       0],
                 [     0,        0,        0,      1/G23,  0,       0],
                 [     0,        0,        0,       0,    1/G13,    0],
                 [     0,        0,        0,       0,     0,   1/G12]])
    return S6

def C6f(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23):
    '''
    daniel pg 74
    transversely isotropic stiffness matrix.
    '''
    C6 = inv(S6f(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23))
    return C6

def Qf(E1,E2,nu12,G12):
    '''transversly isptropic compliance matrix. pg 58 herakovich
    G12 = E1/(2*(1+nu12))  if isotropic'''
    nu21 = E2*nu12/E1
    Q = array([[E1/(1-nu12*nu21),    E2*nu12/(1-nu12*nu21), 0],
               [ E2*nu12/(1-nu12*nu21), E2/(1-nu12*nu21),    0],
               [0,        0,       G12]])
    return Q

def T61(th):
    '''Stress
    th=ply angle in degrees
    voight notation for stress tranform. sigma1 = T1 @ sigmax
    reddy pg 91'''
    n = sin(th*pi/180)
    m = cos(th*pi/180)
    T1 = array( [[m**2, n**2, 0, 0, 0, 2*m*n],
                 [n**2, m**2, 0, 0, 0,-2*m*n],
                 [0,    0,    1, 0, 0, 0],
                 [0,    0,    0, m,-n, 0],
                 [0,    0,    0, n, m, 0],
                 [-m*n, m*n,  0, 0, 0,(m**2-n**2)]])
    return T1

def T62(th):
    '''Strain
    voight notation for strain transform. epsilon1 = T2 @ epsilonx
    th=ply angle in degrees
    reddy pg 91
    '''
    n = sin(th*pi/180)
    m = cos(th*pi/180)
    T2 = array( [[m**2, n**2, 0, 0, 0, m*n],
                 [n**2, m**2, 0, 0, 0,-m*n],
                 [0,    0,    1, 0, 0, 0],
                 [0,    0,    0, m,-n, 0],
                 [0,    0,    0, n, m, 0],
                 [-2*m*n, 2*m*n,  0, 0, 0,(m**2-n**2)]])
    return T2


def T1(th):
    '''Stress Transform for Plane Stress
    th=ply angle in degrees
    voight notation for stress tranform. sigma1 = T1 @ sigmax
    recall T1(th)**-1 == T1(-th)'''
    n = sin(th*pi/180)
    m = cos(th*pi/180)
    T1 = array( [[m**2, n**2, 2*m*n],
                 [n**2, m**2,-2*m*n],
                 [-m*n, m*n,(m**2-n**2)]])
    return T1

def T2(th):
    '''Strain Transform for Plane Stress
    th=ply angle in degrees
    voight notation for strain transform. epsilon1 = T2 @ epsilonx'''
    n = sin(th*pi/180)
    m = cos(th*pi/180)
    T2 = array( [[m**2, n**2, m*n],
                 [n**2, m**2,-m*n],
                 [-2*m*n, 2*m*n,  (m**2-n**2)]])
    return T2

def T1s(th):
    '''Symbolic Stress Transform for Plane Stress
    th=ply angle in degrees
    voight notation for stress tranform. sigma1 = T1 @ sigmax
    recall T1(th)**-1 == T1(-th)'''
    n = sp.sin(th*sp.pi/180)
    m = sp.cos(th*sp.pi/180)
    T1 = sp.Matrix( [[m**2, n**2, 2*m*n],
                 [n**2, m**2,-2*m*n],
                 [-m*n, m*n,(m**2-n**2)]])
    return T1

def T2s(th):
    '''Symbolic Strain Transform for Plane Stress
    th=ply angle in degrees
    voight notation for strain transform. epsilon1 = T2 @ epsilonx'''
    n = sp.sin(th*sp.pi/180)
    m = sp.cos(th*sp.pi/180)
    T2 = sp.Matrix( [[m**2, n**2, m*n],
                 [n**2, m**2,-m*n],
                 [-2*m*n, 2*m*n,  (m**2-n**2)]])
    return T2



def failure_envelope():
    # failure envelopes

    # max stress criteria
    # 1 direction in first row
    # 2 direction in second row

    # failure strength in compression
    #Fc = matrix([[-1250.0, -600.0],
    #            [-200.0,  -120.0]]) # ksi
    #
    ##failure strength in tension
    #Ft =  matrix([[1500, 1000]
    #              [50,     30]]) # ksi
    #
    ##Failure strength in shear
    #Fs = matrix( [100,    70] ) # Shear

    Fc1 = [-1250, -600] # Compression 1 direction
    Fc2 = [-200,  -120] # Compression 2 direction
    Ft1 = [1500, 1000]  # Tension 1 direction
    Ft2 = [50,     30]  # Tension 2 direction
    Fs =  [100,    70]   # Shear

    # F1 = Ft(1);
    # F2 = Ft(1);
    # F6 = Fs(1);

    for c in range(2):# mattype
        factor = 1.25
        # right
        plot( [Ft1[c], Ft1[c]], [Fc2[c], Ft2[c]])

        # left
        plot( [Fc1[c], Fc1[c]] , [Fc2[c], Ft2[c]])
        # top
        plot( [Fc1[c], Ft1[c]] , [Ft2[c], Ft2[c]])
        # bottom
        plot( [Fc1[c], Ft1[c]]  , [Fc2[c], Fc2[c]])
        # center horizontal
        plot( [Fc1[c], Ft1[c]]  , [0, 0])
        # center vertical
        plot( [0, 0]            , [Fc2[c], Ft2[c]])

        #xlim([min(Fc1) max(Ft1)]*factor)
        #ylim([min(Fc2) max(Ft2)]*factor)
        xlabel('$\sigma_1,ksi$')
        ylabel('$\sigma_2,ksi$')
        title('failure envelope with Max-Stress Criteria')

def material_plots(materials = ['Carbon_cloth_AGP3705H']):
    '''
    plotting composite properties
    
    Sf(E1,E2,nu12,G12)
    
    '''
#    plt.rcParams['figure.figsize'] = (10, 8)
#    plt.rcParams['font.size'] = 14
#    plt.rcParams['legend.fontsize'] = 14

                 
    plyangle = arange(-45, 45.1, 0.1)
    h = 1      # lamina thickness
    
    layupname='[0]'
    mat = import_matprops(materials)
    Ex = mat[materials[0]].E1
    Ey = mat[materials[0]].E2
    nuxy = mat[materials[0]].nu12
    Gxy = mat[materials[0]].G12    
             

#    layupname = '[0, 45, 45, 0]'
#    Ex=   2890983.38
#    Ey=   2844063.06
#    nuxy= 0.27
#    Gxy=  1129326.25
#    h = 0.0600   
    
    plt.close('all')
    
    S = Sf(Ex,Ey,nuxy,Gxy)

    C = inv(S)
    
    C11 = [(inv(T1(th)) @ C @ T2(th))[0,0] for th in plyangle]
    C22 = [(inv(T1(th)) @ C @ T2(th))[1,1] for th in plyangle]
    C33 = [(inv(T1(th)) @ C @ T2(th))[2,2] for th in plyangle]
    C12 = [(inv(T1(th)) @ C @ T2(th))[0,1] for th in plyangle]

    Exbar = zeros(len(plyangle))
    Eybar = zeros(len(plyangle))
    Gxybar = zeros(len(plyangle))

    Q = Qf(Ex,Ey,nuxy,Gxy)

    Qbar = zeros((len(plyangle),3,3))
    for i,th in enumerate(plyangle):
        Qbar[i] = solve(T1(th), Q) @ T2(th)
    #Qbar = [solve(T1(th),Q) @ T2(th) for th in plyangle]

    Qbar11 = Qbar[:,0,0]
    Qbar22 = Qbar[:,1,1]
    Qbar66 = Qbar[:,2,2]
    Qbar12 = Qbar[:,0,1]
    Qbar16 = Qbar[:,0,2]
    Qbar26 = Qbar[:,1,2]

    Aij = Qbar*h

    # laminate Stiffness
    #     | Exbar    Eybar    Gxybar   |
    # A = | vxybar   vyxbar   etasxbar |
    #     | etaxsbar etaysbar etasybar |

    # laminate Comnpliance
    aij = zeros((len(plyangle),3,3))
    for i, _Aij in enumerate(Aij):
        aij[i] = inv(_Aij)

    # material properties for whole laminate (Daniel, pg183)
    Exbar  = [1/(h*_aij[0,0]) for _aij in aij]
    Eybar  = [1/(h*_aij[1,1]) for _aij in aij]
    Gxybar = [1/(h*_aij[2,2]) for _aij in aij]

    # Global Stress
    s_xy = array([[100],
                  [10],
                  [5]])

    # local ply stress
    s_12 = np.zeros((3,len(plyangle)))
    for i,th in enumerate(plyangle):
        #s_12[:,i] = np.transpose(T1(th) @ s_xy)[0]   # local stresses
        s_12[:,[i]] = T1(th) @ s_xy

    # Plotting
    figure()#, figsize=(10,8))
    plot(plyangle, C11, plyangle, C22, plyangle, C33, plyangle, C12)
    legend(['$\overline{C}_{11}$','$\overline{C}_{22}$', '$\overline{C}_{44}$', '$\overline{C}_{66}$'])
    title('Transversly Isotropic Stiffness properties of carbon fiber T300_5208')
    xlabel("$\Theta$")
    ylabel('$\overline{C}_{ii}$, ksi')
    grid()

    figure()#, figsize=(10,8))
    plot(plyangle, Exbar, label = r"Modulus: $E_x$")
    plot(plyangle, Eybar, label = r"Modulus: $E_y$")
    plot(plyangle, Gxybar, label = r"Modulus: $G_{xy}$")
    title("Constitutive Properties in various angles")
    xlabel("$\Theta$")
    ylabel("modulus, psi")
    legend()
    grid()

    figure()#,figsize=(10,8))
    plot(plyangle, s_12[0,:], label = '$\sigma_{11},ksi$' )
    plot(plyangle, s_12[1,:], label = '$\sigma_{22},ksi$' )
    plot(plyangle, s_12[2,:], label = '$\sigma_{12},ksi$' )
    legend(loc='lower left')
    xlabel("$\Theta$")
    ylabel("Stress, ksi")
    grid()

    # plot plyangle as a function of time
    figure()#,figsize=(10,8))
    plot(plyangle,Qbar11, label = "Qbar11")
    plot(plyangle,Qbar22, label = "Qbar22")
    plot(plyangle,Qbar66, label = "Qbar66")
    legend(loc='lower left')
    xlabel("$\Theta$")
    ylabel('Q')
    grid()

    # plot plyangle as a function of time
    figure()#,figsize=(10,8))
    plot(plyangle,Qbar12, label = "Qbar12")
    plot(plyangle,Qbar16, label = "Qbar16")
    plot(plyangle,Qbar26, label = "Qbar26")
    legend(loc='lower left')
    xlabel("$\Theta$")
    ylabel('Q')
    grid()

    titlename = 'Laminate Properties varying angle for {} {}'.format(materials[0], layupname)
    #df = pd.DataFrame({'plyangle':plyangle, 'Exbar':Exbar, 'Eybar':Eybar,'Gxybar':Gxybar})
    #print(df)
    #df.to_csv(titlename+'.csv')
    
    plt.figure(figsize=(9,6))
    plot(plyangle, Exbar, label = r"Modulus: $E_x$")
    plot(plyangle, Eybar, label = r"Modulus: $E_y$")
    plot(plyangle, Gxybar, label = r"Modulus: $G_{xy}$")
    title(titlename)
    xlabel("$\Theta$")
    ylabel("modulus, psi")
    legend(loc='best')
    grid()
    
    #plt.savefig(titlename+'.png')

    show()


def laminate_gen(lamthk=1.5, symang=[45,0,90], plyratio=2.0, matrixlayers=False, balancedsymmetric=True):
    '''
    ## function created to quickly create laminates based on given parameters
    lamthk=1.5    # total #thickness of laminate
    symang = [45,0,90, 30]  #symmertic ply angle
    plyratio=2.0  # lamina/matrix ratio
    matrixlayers=False  # add matrix layers between lamina plys
    nonsym=False    # symmetric
    mat = material type, as in different plies, matrix layer, uni tapes, etc
    #ply ratio can be used to vary the ratio of thickness between a matrix ply
         and lamina ply. if the same thickness is desired, plyratio = 1,
         if lamina is 2x as thick as matrix plyratio = 2
    '''
    if matrixlayers:
        nply = (len(symang)*2+1)*2
        nm = nply-len(symang)*2
        nf = len(symang)*2
        tm = lamthk / (plyratio*nf + nm)
        tf = tm*plyratio
        plyangle = zeros(nply//2)
        mat = 2*ones(nply//2)  #  orthotropic fiber and matrix = 1, isotropic matrix=2,
        mat[1:-1:2] = 1   #  [2 if x%2 else 1 for x in range(nply//2) ]
        plyangle[1:-1:2] = symang[:]  # make a copy
        thk = tm*ones(nply//2)
        thk[2:2:-1] = tf
        lamang = list(symang) + list(symang[::-1])
        plyangle = list(plyangle) + list(plyangle[::-1])
        mat = list(mat) + list(mat[::-1])
        thk = list(thk) + list(thk[::-1])
    else: # no matrix layers, ignore ratio
        if balancedsymmetric:
            nply = len(symang)*2
            mat = list(3*np.ones(nply))
            thk = list(lamthk/nply*np.ones(nply))
            lamang = list(symang) + list(symang[::-1])
            plyangle = list(symang) + list(symang[::-1])
        else:
            nply = len(symang)
            mat =[1]*nply
            thk = list(lamthk/nply*np.ones(nply))
            lamang = symang[:]
            plyangle = symang[:]

    return thk,plyangle,mat,lamang

def make_quasi(n0=4,n45=4):
    #n0 = 4
    #n45 = 13
    #
    #ply0 = [0]*n0
    #ply45 = [45]*n45
    #plyangle = []
    #from itertools import zip_longest
    #for x,y in zip_longest(ply0,ply45):
    #    if len(plyangle)<min(len(ply0),len(ply45))*2:
    #        plyangle.append(x)
    #        plyangle.append(y)
    #    else:
    #        plyangle.append(x)
    #        plyangle.reverse()
    #        plyangle.append(y)
    #plyangle = [x for x in plyangle if x is not None]
    #plyangle

    ntot = n45+n0
    plyangle = [45]*int(n45)
    for p in [0]*int(n0):
        plyangle.append(p)
        plyangle.reverse()
    return plyangle

#@xw.func
def laminate_calcs(NM,ek,q0,plyangle,plymatindex,materials,platedim, zoffset,SF,plots,prints):
    '''
    code to compute composite properties, applied mechanical and thermal loads
    and stress and strain

    inputs
    NM # force/moments lbs/in
    ek # strain, curvature  in/in
    q0 = pressure
    plyangle # angle for each ply
    plymatindex  # material for each ply
    materials   # list materials used,

    general outline for computing elastic properties of composites

    1) Determine engineering properties of unidirectional laminate. E1, E2, nu12, G12
    2) Calculate ply stiffnesses Q11, Q22, Q12, Q66 in the principal/local coordinate system
    3) Determine Fiber orientation of each ply
    4) Calculate the transformed stiffness Qxy in the global coordinate system
    5) Determine the through-thicknesses of each ply
    6) Determine the laminate stiffness Matrix (ABD)
    7) Calculate the laminate compliance matrix by inverting the ABD matrix
    8) Calculate the laminate engineering properties

    # Stress Strain Relationship for a laminate, with Q=reduced stiffness matrix
    |sx | |Qbar11 Qbar12 Qbar16| |ex +z*kx |
    |sy |=|Qbar12 Qbar22 Qbar26|=|ey +z*ky |
    |sxy| |Qbar16 Qbar26 Qbar66| |exy+z*kxy|

    # Herakovich pg 84
    Qbar =  inv(T1) @ Q @ T2 == solve(T1, Q) @ T2

    transformation reminders - see Herakovich for details
    sig1 = T1*sigx
    sigx = inv(T1)*sig1
    eps1 = T2*epsx
    epsx = inv(T2)*epsx
    sigx = inv(T1)*Q*T2*epsx
    Qbar = inv(T1)*Q*T2
    Sbar = inv(T2)*inv(Q)*T2


    Notes, core transverse direction is G13, ribbon direction is G23
    
    a_width =  50  # plate width (inches or meters)
    b_length =  50  # laminate length, inches or meters    
    '''

    #==========================================================================
    # Initialize python settings
    #==========================================================================
    #get_ipython().magic('matplotlib')
    plt.close('all')
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 13
    #plt.rcParams['legend.fontsize'] = 14


    #==========================================================================
    # Define composite properties
    #==========================================================================

    assert(len(plyangle)==len(plymatindex))
    a_width, b_length = platedim

    # either apply strains or loads , lb/in
    Nx_, Ny_, Nxy_, Mx_, My_, Mxy_ = NM
    NMbarapp =      array([[Nx_],[Ny_],[Nxy_],[Mx_],[My_],[Mxy_]])

    ex_, ey_, exy_, kx_, ky_, kxy_ = ek
    epsilonbarapp = array([[ex_],[ey_],[exy_],[kx_],[ky_],[kxy_]])
    
    Ti = 0   # initial temperature (C)
    Tf = 0 # final temperature (C)

    #SF = 1.0 # safety factor

    #==========================================================================
    # Import Material Properties
    #==========================================================================
    mat  = import_matprops(materials)
    #mat  = import_matprops(['E-Glass Epoxy cloth','rohacell2lb'])  # Herakovich
    alphaf = lambda mat: array([[mat.alpha1], [mat.alpha2], [0]])

    ''' to get ply material info, use as follows
    alpha = alphaf(mat[materials[plymatindex[i]]])

    mat[materials[1]].E2

    '''

    laminatethk = array([mat[materials[i]].plythk for i in plymatindex ])

    nply = len(laminatethk) # number of plies
    H =   np.sum(laminatethk) # plate thickness
    #    area = a_width*H
    z = zeros(nply+1)
    zmid = zeros(nply)
    z[0] = -H/2
    for i in range(nply):
        z[i+1] = z[i] + laminatethk[i]
        zmid[i] = z[i] + laminatethk[i]/2

    #==========================================================================
    # ABD Matrix Compute
    #==========================================================================
    # Reduced stiffness matrix for a plane stress ply in principal coordinates
    # calcluating Q from the Compliance matrix may cause cancE1ation errors

    A = zeros((3,3)); B = zeros((3,3)); D = zeros((3,3))
    for i in range(nply):  # = nply
        Q = Qf(mat[materials[plymatindex[i]]].E1, mat[materials[plymatindex[i]]].E2, mat[materials[plymatindex[i]]].nu12, mat[materials[plymatindex[i]]].G12 )

        Qbar = solve(T1(plyangle[i]), Q) @ T2(plyangle[i])  # inv(T1(plyangle[i])) @ Q   @ T2(plyangle[i])
        A += Qbar*(z[i+1]-z[i])
        # coupling  stiffness
        B += (1/2)*Qbar*(z[i+1]**2-z[i]**2)
        # bending or flexural laminate stiffness relating moments to curvatures
        D += (1/3)*Qbar*(z[i+1]**3-z[i]**3)

    #Cbar6 = T61 @ C6 @ np.transpose(T61)

    # laminate stiffness matrix
    ABD = zeros((6,6))
    ABD[0:3,0:3] = A
    ABD[0:3,3:6] = B + zoffset*A
    ABD[3:6,0:3] = B + zoffset*A
    ABD[3:6,3:6] = D + 2*zoffset*B + zoffset**2*A

    # laminatee compliance
    abcd = inv(ABD)
    a = abcd[0:3,0:3]

    #==========================================================================
    # Laminate Properties
    #==========================================================================

    # effective laminate shear coupling coefficients
    etasxbar = a[0,2]/a[2,2]
    etasybar = a[1,2]/a[2,2]
    etaxsbar = a[2,0]/a[0,0]
    etaysbar = a[2,1]/a[1,1]

    # laminate engineer properties
    Exbar  = 1 / (H*a[0,0])
    Eybar  = 1 / (H*a[1,1])
    Gxybar = 1 / (H*a[2,2])
    nuxybar = -a[0,1]/a[0,0]
    nuyxbar = -a[0,1]/a[1,1]

    # TODO: validate results, does not appear to be correct
    # strain centers, pg 72, NASA-Basic mechanics of lamianted composites
    # added divide by zero epsilon
    z_eps0_x  = -B[0,0] / (D[0,0] + 1e-16)
    z_eps0_y  = -B[0,1] / (D[0,1] + 1e-16)
    z_eps0_xy = -B[0,2] / (D[0,2] + 1e-16)
    z_sc = -B[2,2] / (D[2,2] +1e-16) # shear center
    
    # --------------------- Double Check ---------------------
#    # Laminate compliance matrix
#    LamComp = array([ [1/Exbar,       -nuyxbar/Eybar,  etasxbar/Gxybar],
#                      [-nuxybar/Exbar,  1/Eybar ,       etasybar/Gxybar],
#                      [etaxsbar/Exbar, etaysbar/Eybar, 1/Gxybar]] )
#    # Daniel pg 183
#    # combines applied loads and applied strains
#    strain_laminate = LamComp @ Nxyzapplied[:3]/H + strainxyzapplied[:3]
#    Nxyz = A @ strain_laminate
#    stress_laminate = Nxyz/H
    # --------------------------------------------------------

    #==========================================================================
    # Pressure Load
    #==========================================================================
    #==========================================================================
    # pressure displacement and moments
    #==========================================================================
    D11,D12,D22,D66 = D[0,0], D[0,1], D[1,1], D[2,2]
    B11 = B[0,0]
    A11, A12 = A[0,0], A[0,1]
    
    # reddy pg 247 Navier  displacement solution for a simply supported plate
    s = b_length/a_width
    x = a_width/2
    y = b_length/2

    # 5.2.8, reddy, or hyer 13.123
    terms = 5
    w0 = 0
    for m in range(1,terms,2):
        for n in range(1,terms,2):
            dmn = pi**4/b_length**4 * (D11*m**4*s**4 + 2*(D12 + 2*D66)*m**2*n**2*s**2 + D22*n**4)
            alpha = m*pi/a_width
            beta = n*pi/b_length
            # for uniformly distributed loads, m,n = 1,3,5,...
            Qmn = 16*q0/(pi**2*m*n)
            Wmn = Qmn/dmn
            w0 += Wmn * sin(alpha*x) * sin(beta*y)
    w0_simplesupport = w0

    # 5.2.12a, reddy
    # mid span moments
    Mxq=Myq=Mxyq=0
    for m in range(1,terms,2):
        for n in range(1,terms,2):
            dmn = pi**4/b_length**4 * (D11*m**4*s**4 + 2*(D12 + 2*D66)*m**2*n**2*s**2 + D22*n**4)
            alpha = m*pi/a_width
            beta = n*pi/b_length
            # for uniformly distributed loads, m,n = 1,3,5,...
            Qmn = 16*q0/(pi**2*m*n)
            Wmn = Qmn/dmn
            Mxq += (D11*alpha**2 + D12*beta**2 ) * Wmn * sin(m*pi*x/a_width) * sin(n*pi*y/b_length)
            Myq += (D12*alpha**2 + D22*beta**2 ) * Wmn * sin(m*pi*x/a_width) * sin(n*pi*y/b_length)
            Mxyq += alpha*beta*D66 * Wmn * cos(m*pi*x/a_width) * cos(n*pi*y/b_length)
    Mxyq = -2*Mxyq
    NMq = [[0],[0],[0],[Mxq],[Myq],[Mxyq]]
    # hyer, x-pin-pin, y-free-free plate reaction forces, pg 619
    # Forces and Moments across the width of the plate
    A11R = A11*(1-B11**2/(A11*D11))
    D11R = D11*(1-B11**2/(A11*D11))
    Nxq0 = lambda x: B11/D11 * q0 * a_width**2 /12
    Nyq0 = lambda x: B11 * A12*q0 * a_width**2 / (D11*A11R*12) * (6*(x/a_width)**2-1/2)
    Nxyq0 = lambda x: 0
    Mxq0 = lambda x: q0 * a_width**2/8 * (1-4*(x/a_width)**2)
    Myq0 = lambda x: D12 * q0 * a_width**2 / (D11R*8) * ((1-2*B11**2/(3*A11*D11))-(4*(x/a_width)**2))
    Mxyq0 = lambda x: 0
    
    
    
    # clamped plate 5.4.11, reddy 
    #w0_clamped = ( 49 * q0*a_width**4 * (x/a_width - (x/a_width)**2 )**2 * (y/b_length - (y/b_length)**2)**2) / (8 * (7*D11+4*(D12 + 2*D66)*s**2 + 7*D22*s**4) )

    # reddy, 5.4.12
    w0_clamped = 0.00342 * (q0*a_width**4) / (D11+0.5714*(D12+2*D66)*s**2+D22*s**4)
    
    # reddy, 5.4.15
    #w0_clamped = 0.00348 * (q0*a_width**4) / (D11*b_length**4+0.6047*(D12+2*D66)*s**2+D22*s**4)
    
    # reddy 5.4.15, for isotropic D11=D
    w0_clamped_isotropic = 0.00134*q0*a_width**4/D11    


    #==========================================================================
    #  Applied Loads and pressure loads
    #==========================================================================
    NMbarapptotal = NMbarapp + NMq + ABD @ epsilonbarapp

    #==========================================================================
    # Thermal Loads
    #==========================================================================
    '''
    if the material is isotropic and unconstrained, then no thermal stresses
        will be experienced. If there are constraints, then the material will experience
        thermally induced stresses. As with orthotropic materials, various directions will have
        different stresses, and when stacked in various orientations, stresses can be
        unintuitive and complicated. Global Thermal strains are subtracted from applied strains
    # 1) determine the free unrestrained thermal strains in each layer, alphabar
    '''
    dT = Tf-Ti

    Nhatth= zeros((3,1))  # unit thermal force in global CS
    Mhatth = zeros((3,1)) # unit thermal moment in global CS
    alphabar = zeros((3,nply))    # global ply CTE
    for i in range(nply):  # = nply
        Q = Qf(mat[materials[plymatindex[i]]].E1, mat[materials[plymatindex[i]]].E2, mat[materials[plymatindex[i]]].nu12, mat[materials[plymatindex[i]]].G12 )
        alpha = alphaf(mat[materials[plymatindex[i]]])
        Qbar = inv(T1(plyangle[i])) @ Q   @ T2(plyangle[i])
        alphabar[:,[i]] = solve(T2(plyangle[i]),  alpha)
        #alphabar[:,[i]] = inv(T2(plyangle[i])) @ alpha # Convert to global CS
        Nhatth += Qbar @ (alphabar[:,[i]])*(z[i+1] - z[i]) # Hyer method for calculating thermal unit loads
        Mhatth += 0.5*Qbar@(alphabar[:,[i]])*(z[i+1]**2-z[i]**2)

    NMhatth = np.vstack((Nhatth,Mhatth))
    NMbarth = NMhatth*dT # resultant thermal loads

    # Laminate CTE
    epsilonhatth = abcd@NMhatth # laminate CTE

    # applied loads and thermal loads
    epsilonbarapp = abcd @ NMbarapptotal
    epsilonbarth  = abcd @ NMbarth  # resultant thermal strains
    epsilonbartotal = epsilonbarapp + epsilonbarth

    # Composite respone from applied mechanical loads and strains. Average
    # properties only. Used to compare results from tensile test.
    #epsilon_laminate = abcd@NMbarapptotal
    #sigma_laminate = ABD@epsilon_laminate/H

    epsilon_laminate = epsilonbartotal[:]
    sigma_laminate = ABD@epsilonbartotal/H
    alpha_laminate = a@Nhatth

    # determine thermal load and applied loads or strains Hyer pg 435,452
    Nx = NMbarapptotal[0,0]*a_width # units kiloNewtons, total load as would be applied in a tensile test
    Ny = NMbarapptotal[1,0]*b_length # units kN

    #==========================================================================
    # Thermal and mechanical local and global stresses at the ply interface
    #==========================================================================
    # Declare variables for plotting
    epsilon_app         = zeros((3,2*nply))
    sigma_app           = zeros((3,2*nply))
    epsilonbar_app      = zeros((3,2*nply))
    sigmabar_app        = zeros((3,2*nply))
    epsilon_th          = zeros((3,2*nply))
    sigma_th            = zeros((3,2*nply))
    epsilonbar_th       = zeros((3,2*nply))
    sigmabar_th         = zeros((3,2*nply))
    epsilon             = zeros((3,2*nply))
    epsilonbar          = zeros((3,2*nply))
    sigma               = zeros((3,2*nply))
    sigmabar            = zeros((3,2*nply))

    for i,k in enumerate(range(0,2*nply,2)):
        # stress is calcuated at top and bottom of each ply
        Q = Qf(mat[materials[plymatindex[i]]].E1, mat[materials[plymatindex[i]]].E2, mat[materials[plymatindex[i]]].nu12, mat[materials[plymatindex[i]]].G12 )
        Qbar = inv(T1(plyangle[i])) @ Q   @ T2(plyangle[i])

        ### transverse shear, herakovich pg 254
        #Q44 = mat[materials[plymatindex[i]]].G23
        #Q55 = mat[materials[plymatindex[i]]].G13
        #Qbar44 = Q44*cos(plyangle[i])**2+Q55*sin(plyangle[i])**2
        #Qbar55 = Q55*cos(plyangle[i])**2 + Q44*sin(plyangle[i])**2
        #Qbar45 = (Q55-Q44)*cos(plyangle[i])*sin(plyangle[i])
        #epsilontransverse = array([[gammayz],[gammaxz]])
        #sigmatransverse = array([[Qbar44, Qbar45],[Qbar45, Qbar55]]) @ epsilontransverse

         # Global stresses and strains, applied load only
        epsbarapp1 = epsilonbarapp[0:3] + z[i]*epsilonbarapp[3:7]
        epsbarapp2 = epsilonbarapp[0:3] + z[i+1]*epsilonbarapp[3:7]
        sigbarapp1 = Qbar @ epsbarapp1
        sigbarapp2 = Qbar @ epsbarapp2
        # Local stresses and strains, appplied load only
        epsapp1 = T2(plyangle[i]) @ epsbarapp1
        epsapp2 = T2(plyangle[i]) @ epsbarapp2
        sigapp1 = Q @ epsapp1
        sigapp2 = Q @ epsapp2
        # Interface Stresses and Strains
        epsilon_app[:,k:k+2]    = np.column_stack((epsapp1,epsapp2))
        epsilonbar_app[:,k:k+2] = np.column_stack((epsbarapp1,epsbarapp2))
        sigma_app[:,k:k+2]      = np.column_stack((sigapp1,sigapp2))
        sigmabar_app[:,k:k+2]   = np.column_stack((sigbarapp1,sigbarapp2))

        # Global stress and strains, thermal loading only
        epsbarth1 = epsilonbarth[0:3] + z[i]*epsilonbarth[3:7]   - dT*alphabar[:,[i]]
        epsbarth2 = epsilonbarth[0:3] + z[i+1]*epsilonbarth[3:7] - dT*alphabar[:,[i]]
        sigbarth1 = Qbar @ epsbarth1
        sigbarth2 = Qbar @ epsbarth2

        # Local stress and strains, thermal loading only
        epsth1 = T2(plyangle[i]) @ epsbarth1
        epsth2 = T2(plyangle[i]) @ epsbarth2
        sigth1 = Q @ epsth1
        sigth2 = Q @ epsth2

        # Interface Stresses and Strains
        epsilon_th[:,k:k+2]    = np.column_stack((epsth1,epsth2))
        epsilonbar_th[:,k:k+2] = np.column_stack((epsbarth1+dT*alphabar[:,[i]],epsbarth2+dT*alphabar[:,[i]])) # remove the local thermal loads for plotting. only use local thermal strains for calculating stress
        sigma_th[:,k:k+2]      = np.column_stack((sigth1,sigth2))
        sigmabar_th[:,k:k+2]   = np.column_stack((sigbarth1,sigbarth2))

        # TOTAL global stresses and strains, applied and thermal
        epsbar1 = epsbarapp1 + epsbarth1
        epsbar2 = epsbarapp2 + epsbarth2
        sigbar1 = Qbar @ epsbar1
        sigbar2 = Qbar @ epsbar2
        # TOTAL local stresses and strains , applied and thermal
        eps1 = T2(plyangle[i]) @ epsbar1
        eps2 = T2(plyangle[i]) @ epsbar2
        sig1 = Q @ eps1
        sig2 = Q @ eps2
        # Interface Stresses and Strains
        epsilon[:,k:k+2]     = np.column_stack((eps1,eps2))
        epsilonbar[:,k:k+2]  = np.column_stack((epsbar1+dT*alphabar[:,[i]],epsbar2+dT*alphabar[:,[i]])) # remove the local thermal loads for plotting. only use local thermal strains for calculating stress
        sigma[:,k:k+2]       = np.column_stack((sig1,sig2))
        sigmabar[:,k:k+2]    = np.column_stack((sigbar1,sigbar2))


    #==========================================================================
    # Strength Failure Calculations
    #==========================================================================
    # Strength Ratio
    STRENGTHRATIO_MAXSTRESS = zeros((3,2*nply))
    # Failure Index
    FAILUREINDEX_MAXSTRESS = zeros((3,2*nply))
    STRENGTHRATIO_TSAIWU = zeros((nply))
    for i,k in enumerate(range(0,2*nply,2)):

        # stress
        s1 = sigma[0,k]
        s2 = sigma[1,k]
        s12 = np.abs(sigma[2,k])

        # strength
        F1 =  mat[materials[plymatindex[i]]].F1t  if s1 > 0 else  mat[materials[plymatindex[i]]].F1c
        F2 =  mat[materials[plymatindex[i]]].F2t  if s2 > 0 else  mat[materials[plymatindex[i]]].F2c
        F12 = mat[materials[plymatindex[i]]].F12

        # Max Stress failure index ,failure if > 1, then fail, FI = 1/SR
        FAILUREINDEX_MAXSTRESS[0,k:k+2] = s1  / F1
        FAILUREINDEX_MAXSTRESS[1,k:k+2] = s2  / F2
        FAILUREINDEX_MAXSTRESS[2,k:k+2] = s12 / F12


        # Tsai Wu, failure occures when > 1
        F1t = mat[materials[plymatindex[i]]].F1t
        F1c = mat[materials[plymatindex[i]]].F1c
        F2t = mat[materials[plymatindex[i]]].F2t
        F2c = mat[materials[plymatindex[i]]].F2c
        F12 = mat[materials[plymatindex[i]]].F12

        # inhomogeneous Tsai-Wu criterion # from Daniel
        # http://www2.mae.ufl.edu/haftka/composites/mcdaniel-nonhomogenous.pdf
        f1 =  1/F1t + 1/F1c 
        f2 =  1/F2t + 1/F2c
        f11 = -1/(F1t*F1c)
        f22 = -1/(F2t*F2c)
        f66 = 1/F12**2
        f12 = -0.5*sqrt(f11*f22)
        #TW = f1*s1 + f2*s2 + f11*s1**2 + f22*s2**2 + f66*s12**2 + 2*f12*s1*s2
        # polynomial to solve. Added a machine epsilon to avoid divide by zero errors
        lam1 = f11*s1**2 + f22*s2**2 + f66*s12**2 + 2*f12*s1*s2 + 1e-16
        lam2 = f1*s1 + f2*s2 + 1e-16
        lam3 = -1 
        # smallest positive root
        roots = array([(-lam2+sqrt(lam2**2-4*lam1*lam3)) / (2*lam1) ,
                       (-lam2-sqrt(lam2**2-4*lam1*lam3)) / (2*lam1)] )
        STRENGTHRATIO_TSAIWU[i] = roots[roots>=0].min()  # strength ratio

#        f1 =  1/F1t - 1/F1c
#        f2 =  1/F2t - 1/F2c
#        f11 = 1/(F1t*F1c)
#        f22 = 1/(F2t*F2c)
#        f66 = 1/F12**2
#        STRENGTHRATIO_TSAIWU[i] =  2 / (f1*s2 + f2*s2 + sqrt((f1*s1+f2*s2)**2+4*(f11*s1**2+f22*s2**2+f66*s12**2)))

    ### Apply safety factors
    FAILUREINDEX_MAXSTRESS = FAILUREINDEX_MAXSTRESS * SF
    STRENGTHRATIO_TSAIWU = STRENGTHRATIO_TSAIWU / SF
    
    ### 
    MARGINSAFETY_TSAIWU = STRENGTHRATIO_TSAIWU-1   # margin of safety

    # strength ratio for max stress, if < 1, then fail, SR = 1/FI
    STRENGTHRATIO_MAXSTRESS = 1/(FAILUREINDEX_MAXSTRESS+1e-16)
    # margin of safety based on max stress criteria
    MARGINSAFETY_MAXSTRESS = STRENGTHRATIO_MAXSTRESS-1
    
    # minimum margin of safety for Max stress failure
    MARGINSAFETY_MAXSTRESS_min = MARGINSAFETY_MAXSTRESS.min().min()
    FAILUREINDEX_MAXSTRESS_max = FAILUREINDEX_MAXSTRESS.max().max()
    
    # minimum margin of safety of both Tsai-Wu and Max Stress
    #MARGINSAFETY_MAXSTRESS_min = np.minimum(MARGINSAFETY_MAXSTRESS.min().min(), MARGINSAFETY_TSAIWU.min() )
    
    # find critial values for all failure criteria
    #MARGINSAFETY_MAXSTRESS = MARGINSAFETY_MAXSTRESS[~np.isinf(MARGINSAFETY_MAXSTRESS)] # remove inf
    #MARGINSAFETY_TSAIWU = MARGINSAFETY_TSAIWU[~np.isinf(MARGINSAFETY_TSAIWU)] # remove inf


    #==========================================================================
    # Buckling Failure Calculations
    #==========================================================================
    ''' Buckling of Clamped plates under shear load, reddy, 5.6.17'''

    
    k11 = 537.181*D11/a_width**4 + 324.829*(D12+2*D66)/(a_width**2*b_length**2) + 537.181*D22/b_length**4
    k12 = 23.107/(a_width*b_length)
    k22 = 3791.532*D11/a_width**4 + 4227.255*(D12+2*D66)/(a_width**2*b_length**2) + 3791.532*D22/b_length**4
    Nxycrit0 = 1/k12*np.sqrt(k11*k22)
    FI_clamped_shear_buckling = (abs(Nxy_)*SF) / Nxycrit0  # failure if > 1

    MS_clamped_shear_buckling = 1/(FI_clamped_shear_buckling+1e-16)-1
    '''Kassapoglous pg 126,137

    simply supported plate buckling, assumes Nx>0 is compression

    Nxcrit0 is the axial load that causes buckling
    Nxycrit0 is the shear load that cause buckling

    Nxcrit is the axial load part of a combined load that causes buckling
    Nxycrit is the shear load part of a combined load that causes buckling
    '''

    # no buckling issues if Nx is positive
    # buckling calcuations assumes Nx compression is positive.
    Nx__ = abs(Nx_) if Nx_ < 0 else np.float64(0)
    Nxy__ = np.float64(0) if Nxy_ == 0 else abs(Nxy_) # assume shear in 1 direction although both directions are ok

    # Nxy=0
    Nxcrit0 = pi**2/a_width**2 * (D11 + 2*(D12 + 2*D66)*a_width**2/b_length**2 + D22*a_width**4/b_length**4)

    # Nx=0
    Nxycrit0 =  9*pi**4*b_length / (32*a_width**3) * (D11 + 2*(D12 + 2*D66)*a_width**2/b_length**2 + D22*a_width**4/b_length**4)

    FI_Nxy0_buckling, FI_Nx0_buckling, FI_Nx_buckling, FI_Nxy_buckling = 0,0,0,0

    if Nx__ == 0 or Nxy__ == 0:
        FI_Nxy0_buckling =  (Nxy__*SF)/Nxycrit0
        FI_Nx0_buckling =  (Nx__*SF)/Nxcrit0
    else:
        # interaction term
        k = Nxy__ / Nx__
        Nxcrit = min( abs((pi**2/a_width**2) * (D11 + 2*(D12 + 2*D66)*a_width**2/b_length**2 +D22*a_width**4/b_length**4 ) / (2-8192*a_width**2*k**2/(81*b_length**2*pi**4)) * (5 + sqrt(9 + 65536*a_width**2*k**2/(81*pi**4*b_length**2)))) ,
                      abs((pi**2/a_width**2) * (D11 + 2*(D12 + 2*D66)*a_width**2/b_length**2 +D22*a_width**4/b_length**4 ) / (2-8192*a_width**2*k**2/(81*b_length**2*pi**4)) * (5 - sqrt(9 + 65536*a_width**2*k**2/(81*pi**4*b_length**2)))) )

        Nxycrit = Nxycrit0*sqrt(1-Nxcrit/Nxcrit0)
        # interactive calc
        FI_Nx_buckling =  (Nx__ *SF)/Nxcrit
        FI_Nxy_buckling =  (Nxy__*SF)/Nxycrit

    FI_combinedload_simplesupport_buckle = max([FI_Nxy0_buckling,
                                                FI_Nx0_buckling,
                                                FI_Nx_buckling,
                                                FI_Nxy_buckling] )

    MS_min_buckling = 1/(FI_combinedload_simplesupport_buckle+1e-16)-1
    
    #==========================================================================
    # Facesheet Wrinkling
    #==========================================================================


    #==========================================================================
    # principal lamainte stresses
    #==========================================================================
    sigma_principal_laminate = np.linalg.eig(array([[sigma_laminate[0,0],sigma_laminate[2,0],0],
                                                   [sigma_laminate[2,0],sigma_laminate[1,0],0],
                                                   [0,0,0]]))[0]
    tauxy_p = sigma_laminate[2,0]
    sigmax_p = sigma_laminate[0,0]
    sigmay_p = sigma_laminate[1,0]
    thetap = 0.5 * np.arctan( 2*tauxy_p / ((sigmax_p-sigmay_p+1e-16))) * 180/np.pi
                  
    #==========================================================================
    # Printing Results
    #==========================================================================
    if prints:
        print('--------------- laminate1 Stress analysis of fibers----------')
        print('(z-) plyangles (z+)'); print(plyangle)
        print('(z-) plymatindex (z+)'); print(plymatindex)
        print('ply layers') ; print(z)
        print('lamiante thickness, H = {:.4f}'.format(H))
        
        #print('x- zero strain laminate center, z_eps0_x  = {:.4f}'.format(z_eps0_x))
        #print('y- zero strain laminate center, z_eps0_y  = {:.4f}'.format(z_eps0_y))
        #print('xy-zero strain laminate center, z_eps0_xy = {:.4f}'.format(z_eps0_xy))       
        #print('shear center laminate center, z_sc = {:.4f}'.format(z_sc))  
        
        print('Applied Loads'); print(NM)
        print('ABD=');print(ABD)
        print('Ex=   {:.2f}'.format(Exbar) )
        print('Ey=   {:.2f}'.format(Eybar) )
        print('nuxy= {:.2f}'.format(nuxybar) )
        print('Gxy=  {:.2f}'.format(Gxybar) )
        print('epsilon_laminate') ; print(epsilon_laminate)
        print('sigma_laminate') ; print(sigma_laminate)
        print('sigma_principal_laminate') ; print(sigma_principal_laminate)
        print('principal_angle = {:.2f} deg'.format(thetap))
        print('NMbarapp') ; print(NMbarapp)
        print('sigma') ; print(sigma)
    
        print('\nMax Stress Percent Margin of Safety, failure < 0, minimum = {:.4f}'.format( MARGINSAFETY_MAXSTRESS_min ) )
        print(MARGINSAFETY_MAXSTRESS)
        
        print('\nTsai-Wu Percent Margin of Safety, failure < 0, minimum = {:.4f}'.format(MARGINSAFETY_TSAIWU.min()))
        print(MARGINSAFETY_TSAIWU)
        
        print('\nmaximum failure index = {:.4f}'.format(  FAILUREINDEX_MAXSTRESS_max ))      
        print(FAILUREINDEX_MAXSTRESS)
        
        print('\nBuckling MS for Nxy only for clamped edges = {:.4f}\n'.format(MS_clamped_shear_buckling))

    #    print('---- Individual Buckling Failure Index (fail>1) combined loads and simple support -----')
    #    print('FI_Nxy0 = {:.2f}'.format(FI_Nxy0_buckling) )
    #    print('FI_Nx0  = {:.2f}'.format(FI_Nx0_buckling) )
    #    print('---- Interactive Buckling Failure Index (fail>1) combined loads and simple support -----')
    #    print('FI_Nx   = {:.2f}'.format(FI_Nx_buckling) )
    #    print('FI_Nxy  = {:.2f}'.format(FI_Nxy_buckling) )
    #    print('---- Buckling Failure Index (fail>1) combined loads and simple support -----')
    #    print(FI_combinedload_simplesupport_buckle)
        print('buckling combined loads and simple support MS = {:.4f}\n'.format((MS_min_buckling)))
        print('Mx_midspan = {:.2f}'.format(Mxq) )
        print('My_midspan = {:.2f}'.format(Myq) ) 
        print('Mxy_midspan = {:.2f}'.format(Mxyq) )    
        print('w0_simplesupport =    {:.6f}'.format(w0_simplesupport) )    
        print('w0_clamped =          {:.6f}'.format(w0_clamped) )
        print('w0_clamped_isotropic= {:.6f}'.format(w0_clamped_isotropic) )
    
    #display(sp.Matrix(sigmabar))

    #==========================================================================
    # Plotting
    #==========================================================================
    if plots:

        windowwidth = 800
        windowheight = 450
        zplot = zeros(2*nply)
        for i,k in enumerate(range(0,2*nply,2)):  # = nply
            zplot[k:k+2] = z[i:i+2]

        #legendlab = ['total','thermal','applied','laminate']
        # global stresses and strains
        mylw = 1.5 #linewidth
        # Global Stresses and Strains
        f1, ((ax1,ax2,ax3), (ax4,ax5,ax6)) = plt.subplots(2,3, sharex='row', sharey=True)
        f1.canvas.set_window_title('Global Stress and Strain of %s laminate' % (plyangle))
        stresslabel = ['$\sigma_x$','$\sigma_y$','$\\tau_{xy}$']
        strainlabel = ['$\epsilon_x$','$\epsilon_y$','$\gamma_{xy}$']

        for i,ax in enumerate([ax1,ax2,ax3]):
            ## the top axes
            ax.set_ylabel('thickness,z')
            ax.set_xlabel(strainlabel[i])
            ax.set_title(' Ply Strain '+strainlabel[i])
            ax.ticklabel_format(axis='x', style='sci', scilimits=(1,4))  # scilimits=(-2,2))
            ax.plot(epsilonbar[i,:],     zplot, color='blue', lw=mylw, label='total')
            ax.plot(epsilonbar_th[i,:],  zplot, color='red', lw=mylw, alpha=0.75, linestyle='--',  label='thermal')
            ax.plot(epsilonbar_app[i,:], zplot, color='green', lw=mylw, alpha=0.75,linestyle='-.', label='applied')
            ax.plot([epsilon_laminate[i], epsilon_laminate[i]],[np.min(z) , np.max(z)], color='black', lw=mylw, label='laminate')
            ax.grid(True)
            #ax.set_xticks(linspace( min(ax.get_xticks()) , max(ax.get_xticks()) ,6))

        for i,ax in enumerate([ax4,ax5,ax6]):
            ax.set_ylabel('thickness,z')
            ax.set_xlabel(stresslabel[i])
            ax.set_title(' Ply Stress '+stresslabel[i])
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3)) # scilimits=(-2,2))
            ax.plot(sigmabar[i,:],     zplot, color='blue', lw=mylw, label='total')
            ax.plot(sigmabar_th[i,:], zplot, color='red', lw=mylw, alpha=0.75,linestyle='--', label='thermal')
            ax.plot(sigmabar_app[i,:], zplot, color='green', lw=mylw, alpha=0.75,linestyle='-.', label='applied')
            ax.plot([sigma_laminate[i], sigma_laminate[i]],[np.min(z) , np.max(z)], color='black', lw=mylw, label='laminate')
            ax.grid(True)

        leg = legend(fancybox=True) ; leg.get_frame().set_alpha(0.3)
        tight_layout()
        
        try:
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(25,50,windowwidth,windowheight)
        except:
            pass
        f1.show()
        #plt.savefig('global-stresses-strains.png')

        ### Local Stresses and Strains
        f2, ((ax1,ax2,ax3), (ax4,ax5,ax6)) = plt.subplots(2,3, sharex='row', sharey=True)
        f2.canvas.set_window_title('Local Stress and Strain of %s laminate' % (plyangle))
        stresslabel = ['$\sigma_1$','$\sigma_2$','$\\tau_{12}$']
        strainlabel = ['$\epsilon_1$','$\epsilon_2$','$\gamma_{12}$']
        strengthplot = [ [ [F1t,F1t],[zplot.min(), zplot.max()], [F1c,  F1c],[zplot.min(), zplot.max()] ] ,
                         [ [F2t,F2t],[zplot.min(), zplot.max()], [F2c,  F2c],[zplot.min(), zplot.max()] ] ,
                         [ [F12,F12],[zplot.min(), zplot.max()], [-F12,-F12],[zplot.min(), zplot.max()] ] ]

        for i,ax in enumerate([ax1,ax2,ax3]):
            ## the top axes
            ax.set_ylabel('thickness,z')
            ax.set_xlabel(strainlabel[i])
            ax.set_title(' Ply Strain '+strainlabel[i])
            ax.ticklabel_format(axis='x', style='sci', scilimits=(1,4))  # scilimits=(-2,2))
            ax.plot(epsilon[i,:],     zplot, color='blue', lw=mylw, label='total')
            ax.plot(epsilon_th[i,:], zplot, color='red', lw=mylw, alpha=0.75,linestyle='--', label='thermal')
            ax.plot(epsilon_app[i,:], zplot, color='green', lw=mylw, alpha=0.75,linestyle='-.', label='applied')
            ax.plot([epsilon_laminate[i], epsilon_laminate[i]],[np.min(z) , np.max(z)], color='black', lw=mylw, label='laminate')
            ax.grid(True)

        for i,ax in enumerate([ax4,ax5,ax6]):
            ax.set_ylabel('thickness,z')
            ax.set_xlabel(stresslabel[i])
            ax.set_title(' Ply Stress '+stresslabel[i])
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3,3)) # scilimits=(-2,2))
            ax.plot(sigma[i,:],     zplot, color='blue', lw=mylw, label='total')
            ax.plot(sigma_th[i,:], zplot, color='red', lw=mylw, alpha=0.75,linestyle='--', label='thermal')
            ax.plot(sigma_app[i,:], zplot, color='green', lw=mylw, alpha=0.75,linestyle='-.', label='applied')
            ax.plot([sigma_laminate[i], sigma_laminate[i]],[np.min(z) , np.max(z)], color='black', lw=mylw, label='laminate')
            ### plots strengths
            #ax.plot(strengthplot[i][0],strengthplot[i][1], color='yellow', lw=mylw)
            ax.grid(True)

        leg = legend(fancybox=True) ; leg.get_frame().set_alpha(0.3)
        tight_layout()
        try:
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(windowwidth+50,50,windowwidth,windowheight)
        except:
            pass
        f2.show()
        #plt.savefig('local-stresses-strains.png')

        ### Failure
        f3, ((ax1,ax2,ax3)) = plt.subplots(1,3, sharex=True, sharey=True)
        f3.canvas.set_window_title('Failure Index(failure if > 1),  %s laminate' % (plyangle))
        stresslabel = ['$\sigma_1/F_1$','$\sigma_2/F_2$','$\\tau_{12}/F_{12}$']
        for i,ax in enumerate([ax1,ax2,ax3]):
            ## the top axes
            ax.set_ylabel('thickness,z')
            ax.set_xlabel(stresslabel[i])
            #ax.set_title(' Ply Strain at $\epsilon=%f$' % (epsxapp*100))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(1,4))  # scilimits=(-2,2))
            ax.plot(FAILUREINDEX_MAXSTRESS[i,:],     zplot, color='blue', lw=mylw, label='total')
            ax.grid(True)
            ax.set_title('Failure Index, fail if > 1')
        #leg = legend(fancybox=True) ; leg.get_frame().set_alpha(0.3)
        tight_layout()
        try:
            mngr = plt.get_current_fig_manager() 
            mngr.window.setGeometry(25,windowheight+100,windowwidth,windowheight)
        except:
            pass
        f2.show()
        #plt.savefig('local-stresses-strains.png')

        ### warpage
        res = 100
        Xplt,Yplt = np.meshgrid(np.linspace(-a_width/2,a_width/2,res), np.linspace(-b_length/2,b_length/2,res))
        epsx = epsilon_laminate[0,0]
        epsy = epsilon_laminate[1,0]
        epsxy = epsilon_laminate[2,0]
        kapx = epsilon_laminate[3,0]
        kapy = epsilon_laminate[4,0]
        kapxy = epsilon_laminate[5,0]
        ### dispalcement
        w = -0.5*(kapx*Xplt**2 + kapy*Yplt**2 + kapxy*Xplt*Yplt)
        u = epsx*Xplt  # pg 451 hyer
        fig = plt.figure('plate-warpage')
        ax = fig.gca(projection='3d')
        ax.plot_surface(Xplt, Yplt, w+zmid[0], cmap=mpl.cm.jet, alpha=0.3)
        ###ax.auto_scale_xyz([-(a_width/2)*1.1, (a_width/2)*1.1], [(b_length/2)*1.1, (b_length/2)*1.1], [-1e10, 1e10])
        ax.set_xlabel('plate width,y-direction,in')
        ax.set_ylabel('plate length,x-direction, in')
        ax.set_zlabel('warpage,in')
        #ax.set_zlim(-0.01, 0.04)
        #mngr = plt.get_current_fig_manager() ; mngr.window.setGeometry(450,550,600, 450)
        try:
            mngr = plt.get_current_fig_manager() 
            mngr.window.setGeometry(windowwidth+50,windowheight+100,windowwidth,windowheight)
        except:
            pass
        plt.show()
        #plt.savefig('plate-warpage')

    return MARGINSAFETY_MAXSTRESS_min, FAILUREINDEX_MAXSTRESS_max



def plate():
    '''
    composite plate mechanics

    TODO - results need vetted
    '''


    #==========================================================================
    # Initialize
    #==========================================================================
    get_ipython().magic('matplotlib')
    plt.close('all')
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 13
    #plt.rcParams['legend.fontsize'] = 14

    #==========================================================================
    # Import Material Properties
    #==========================================================================

    plythk = 0.0025
    plyangle = array([0,90,-45,45,0]) * np.pi/180 # angle for each ply
    nply = len(plyangle) # number of plies
    laminatethk = np.zeros(nply) + plythk
    H =   sum(laminatethk) # plate thickness
    # Create z dimensions of laminate
    z_ = np.linspace(-H/2, H/2, nply+1)

    a =   20  # plate width;
    b =   10  # plate height
    q0_ = 5.7 # plate load;
    # Transversly isotropic material properties
    E1 = 150e9
    E2 = 12.1e9
    nu12 = 0.248
    G12 = 4.4e9
    nu23 = 0.458
    G23 = E2 / (2*(1+nu23))
    # Failure Strengths
    F1t =  1500e6
    F1c = -1250e6
    F2t =  50e6
    F2c = -200e6
    F12t =  100e6
    F12c =  -100e6
    Strength = np.array([[F1t, F1c],
                            [F2t, F2c],
                            [F12t, F12c]])


    th = sp.symbols('th')

    # Stiffnes matrix in material coordinates
    Cijm6 = inv(Sij6)


    # reduced stiffness in structural
    Cij = sp.Matrix([[Cij6[0,0], Cij6[0,1], 0],
                     [Cij6[0,1], Cij6[1,1], 0],
                     [0, 0, Cij6[5,5] ]] )

    Tij = sp.Matrix([[cos(th)**2, sin(th)**2, 2*sin(th)*cos(th)],
                     [sin(th)**2, cos(th)**2, -2*sin(th)*cos(th)],
                     [-cos(th)*sin(th), sin(th)*cos(th), (cos(th)**2-sin(th)**2)]])



    ## Cylindrical Bending of a laminated plate

    # displacement in w (z direction)
    from sympy.abc import x
    f = Function('f')
    eq = dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x), hint = '1st_homogeneous_coeff_best', simplify=False)
    pprint(eq)
    #==============================================================================

    th,x,y,z,q0,C1,C2,C3,C4,C5,C6,C7,A11,B11,D11,A16,B16 = symbols('th x y z q0 C1 C2 C3 C4 C5 C6 C7 A11 B11 D11 A16 B16')

    wfun = Function('wfun')
    ufun = Function('ufun')

    ## EQ 4.4.1a
    eq1 = A11*ufun(x).diff(x,2) - B11*wfun(x).diff(x,3)
    #eq1   = A11*diff(ufun,x,2) - B11*diff(wfun,x,3); # C5 C1

    ## EQ 4.4.1b
    #eq2   = A16*diff(ufun,x,2) - B16*diff(wfun,x,3); # C5 C1
    eq2 = A16*ufun(x).diff(x,2) - B16*wfun(x).diff(x,3)

    ## EQ 4.4.1c
    #eq3 = B11*diff(ufun,x,3) - D11*diff(wfun,x,4) + q0;
    eq3 = B11*ufun(x).diff(x,3) - D11*wfun(x).diff(x,4) + q0

    ################## python conversion eded here ################################

    # solve eq1 eq2 and eq3 to get the w and u functions

    # displacement in w (z direction) from eq1,eq2,eq3
    wfun = A11*q0*x**4 / (4*(6*B11**2-6*A11*D11)) + C1 + C2*x + C3*x**2 + C4*x**3 #  C1 C2 C3 C4

    # displacement in u (x direction) from eq1,eq2,eq3
    ufun = B11*q0*x**3 / (6*(B11**2-A11*D11)) + C7 + x*C6 + 3*B11*x**2*C5/A11 # C5 C6 C7

    # Cij6.evalf(subs={th:plyangle[i]}) * (z_[i+1]**3-z_[i]**3)

    # cond1 -> w(0)=0 at x(0), roller
    C1sol = sp.solve(wfun.subs(x,0), C1)[0] # = 0
    # cond2 -> angle at dw/dx at x(0) is 0, cantilever
    C2sol = sp.solve(wfun.diff(x).subs(x,0),C2)[0]  # =  0
    # cond3 -> w(z) = 0 at x(a), roller
    C4sol1 =  sp.solve(wfun.subs({x:a,C1:C1sol,C2:C2sol}),C4)[0] # C3
    # cond4 u = 0 at x = 0
    C7sol = sp.solve(ufun.subs(x,0),C7)[0] #=0
    # u=0 at x = a
    C5sol1 = sp.solve(ufun.subs({x:a, C7:C7sol}),C5)[0] #C6
    # cond 5 EQ 4.4.14a Myy = 0 @ x(a) (Mxx , B11 D11) (Myy, B12 D12) roller no moment
    C6sol1 = sp.solve( ( ((B11*ufun.diff(x)+0.5*wfun.diff(x)**2 ) - D11*wfun.diff(x,2)).subs({x:a, C1:C1sol, C2:C2sol, C4:C4sol1, C5:C5sol1, C7:C7sol})), C6)[0] # C6 C3
    # EQ 4.4.13a, Nxx = 0 @ x(0) roller has no Nxx
    C6sol2 = sp.solve( ((A11* ufun.diff(x) + 0.5*wfun.diff(x)**2)-B11*wfun.diff(x,2)).subs({x:a, C1:C1sol, C2:C2sol, C4:C4sol1, C5:C5sol1, C7:C7sol}),C6)[0] # C6 C3
    C3sol = sp.solve(C6sol1 - C6sol2,C3)[0]
    C4sol = C4sol1.subs(C3,C3sol)
    C6sol = sp.simplify(C6sol2.subs(C3,C3sol))
    C5sol = sp.simplify(C5sol1.subs(C6,C6sol))
    # substitute integration constants with actual values( _ is actual number)
    C1_ = copy(C1sol)
    C2_ = copy(C2sol)
    C7_ = copy(C7sol)
    C3_ = C3sol.subs({q0:q0_, A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})
    C4_ = C4sol.subs({q0:q0_, A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})
    C5_ = C5sol.subs({q0:q0_, A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})
    C6_ = C6sol.subs({q0:q0_, A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})

    # function w(x) vertical displacement w along z with actual vaules
    wsol = wfun.subs({q0:q0_, C1:C1_, C2:C2_, C3:C3_, C4:C4_,  A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})
    # function u(x) horizontal displacement u along x with actual vaules
    usol = ufun.subs({q0:q0_, C5:C5_, C6:C6_, C7:C7_,  A11:Aij[0,0], B11:Bij[0,0], D11:Dij[0,0]})

    # 3d plots
    plot3d(wsol,(x,0,a), (y,0,b))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Cylindrical Bending -Displacement of a plate With CLPT')

    ## Strain calculation
    # eq 3.3.8 (pg 116 reddy (pdf = 138))
    epstotal = array([[usol.diff(x) + 0.5* wsol.diff(x)**5 - z*wsol.diff(x,2)],[0],[0]])
    epsx = epstotal[0,0]
    ## Calculating and plotting Stress in each layer
    res = 8 # accuracy of finding max and min stress
    xplot = linspace(0,a,res)
    yplot = linspace(0,b,res)
    G0 = sp.symbols('G0')
    Globalminstress = np.zeros((3, nply))
    Globalmaxstress = np.zeros((3, nply))

    for kstress in range(3): # stress state s_x, s_y, s_xz
        plt.figure(kstress+1)

        for klay in range(nply): # loop through all layers
            thplot = plyangle[klay]
            zplot = linspace(z_[klay],z_[klay+1],res)
            stressplot = np.zeros((len(zplot),len(xplot)))
            ## Calc Stresses
            if kstress == 2:
                # Shear stresses

                G0_ = -sp.integrate(s_stress[0].diff(x),z)+G0
                # solve for shear stresses from s_1
                s_xz = sp.solve(G0_,G0)[0]
                # out of plane shear S_xz does not need to be transformed ??
                plot3d(s_xz, (x,0, a), (z, z_[klay], z_[klay+1]) )
            else:
                # normal stresses
                # Cij = reduced structural stiffness in strictural coordinates 3x3
                # stress in structural coordinates
                s_stress = Cij.subs(th,thplot) @ epstotal
                # stressin material coordinates
                m_stress = Tij.subs(th,thplot) @ s_stress

                #ezsurf(m_stress(kstress),[0,a,z_(klay),z_(klay+1)])

            ## find max stress in each layer
            ii=0
            for i in xplot:
                jj=0
                for j in zplot:
                    if kstress == 2:
                        stressplot[ii,jj] = s_xz.subs({x:i, z:j})
                    else:
                        stressplot[ii,jj] = m_stress[kstress].subs({x:i, z:j})
                    jj+=jj
                ii+=ii

            Globalminstress[kstress,klay] = np.min(stressplot)
            Globalmaxstress[kstress,klay] = np.max(stressplot)
            #

        plt.title('\sigma_%i' % kstress)

    ## Plot max stress and failure strength
    plt.figure()
    for i in range(3):

        plt.subplot(1, 3, i+1)
        plt.bar(range(nply), Globalmaxstress[i,:])

        plt.bar(range(nply), Globalminstress[i,:])
        plt.scatter(range(nply),np.ones(nply) * Strength[i,0])
        plt.scatter(range(nply),np.ones(nply) * Strength[i,1])

        plt.xlabel('layer')
        plt.title('\sigma%i' % i)


def plate_navier():
    '''
    composite plate bending with navier solution

    TODO - code needs to be converted from matlab
    '''

    ## Plate a*b*h simply supported under q = q0 CLPT

    pass
    '''
    q0,a,b,m,n,x,y = sp.symbols('q0 a b m n x y')

    Qmn = 4/(a*b)*sp.integrate( sp.integrate( q0*sp.sin(m*pi*x/a)*sp.sin(n*pi*y/b),(x,0,a)) ,(y,0,b))

    dmn = pi**4 / b**4 * (DTij(1,1)*m**4*(b/a)**4 + 2* (DTij(1,2)+2*DTij(6,6)) *m**2*n**2*(b/a)**2 + DTij(2,2)*n**4)

    Wmn = Qmn/dmn;

    w0 = Wmn * sin(m*pi*x/a) * sin(n*pi*y/b);

    w0_ = subs(w0,[q0 a b],[-q0_ a_ b_] );

    figure
    w0sum = 0;
    for n_ = 1:10
        for m_ = 1:10
            w0sum = w0sum + subs(w0_,[n m],[n_ m_]);
        end
    end
    w0sum;

    % xplot = linspace(0,a_,res);
    % yplot = linspace(0,b_,res);

    ii=1;
    for i = xplot
        jj=1;
        for j = yplot
            w0plot(ii,jj) = subs(w0sum,[x y],[i j]);
            jj=jj+1;
        end
        ii=ii+1;
    end

    surf(xplot,yplot,w0plot)
    colorbar
    set(gca,'PlotBoxAspectRatio',[2 1 1]);
    xlabel('length a, u(x)')
    ylabel('length b, v(y)')
    zlabel('w(z)')
    '''


class laminate(object):
    """
    IN-WORK - laminate object for composite material analysis
    """

    # constructor
    def __init__(self, plyangle, matindex, matname):
        # run when laminate is instantiated

        # loads materials used
        self.plyangle = plyangle
        self.matindex = matindex
        self.matname = matname


        self.__mat = self.__import_matprops(matname)

        # create a simple function to handle CTE properties
        def __alphaf(self, mat):
            return array([[mat.alpha1], [mat.alpha2], [0]])

        self.laminatethk = array([self.__mat[matname[i]].plythk for i in matindex ])

        self.nply = len(self.laminatethk) # number of plies
        self.H =   np.sum(self.laminatethk) # plate thickness
        #    area = a_width*H
        z = zeros(self.nply+1)
        zmid = zeros(self.nply)
        z[0] = -self.H/2
        for i in range(self.nply):
            z[i+1] = z[i] + self.laminatethk[i]
            zmid[i] = z[i] + self.laminatethk[i]/2
        self.z = z
        self.zmid = zmid

        self.__abdmatrix()

    def __Qf(self, E1,E2,nu12,G12):
        '''transversly isptropic compliance matrix. pg 58 herakovich
        G12 = E1/(2*(1+nu12))  if isotropic'''
        nu21 = E2*nu12/E1
        Q = array([[E1/(1-nu12*nu21),    E2*nu12/(1-nu12*nu21), 0],
                   [ E2*nu12/(1-nu12*nu21), E2/(1-nu12*nu21),    0],
                   [0,        0,       G12]])
        return Q


    def __T1(self, th):
        '''Stress Transform for Plane Stress
        th=ply angle in degrees
        voight notation for stress tranform. sigma1 = T1 @ sigmax
        recall T1(th)**-1 == T1(-th)'''
        n = sin(th*pi/180)
        m = cos(th*pi/180)
        T1 = array( [[m**2, n**2, 2*m*n],
                     [n**2, m**2,-2*m*n],
                     [-m*n, m*n,(m**2-n**2)]])
        return T1

    def __T2(self, th):
        '''Strain Transform for Plane Stress
        th=ply angle in degrees
        voight notation for strain transform. epsilon1 = T2 @ epsilonx'''
        n = sin(th*pi/180)
        m = cos(th*pi/180)
        T2 = array( [[m**2, n**2, m*n],
                     [n**2, m**2,-m*n],
                     [-2*m*n, 2*m*n,  (m**2-n**2)]])
        return T2

    # private method
    def __abdmatrix(self):
        '''used within the object but not accessible outside'''
        #==========================================================================
        # ABD Matrix Compute
        #==========================================================================
        # Reduced stiffness matrix for a plane stress ply in principal coordinates
        # calcluating Q from the Compliance matrix may cause cancE1ation errors

        A = zeros((3,3)); B = zeros((3,3)); D = zeros((3,3))
        for i in range(self.nply):  # = nply
            Q = self.__Qf(self.__mat[self.matname[self.matindex[i]]].E1,
                   self.__mat[self.matname[self.matindex[i]]].E2,
                   self.__mat[self.matname[self.matindex[i]]].nu12,
                   self.__mat[self.matname[self.matindex[i]]].G12 )

            Qbar = inv(self.__T1(self.plyangle[i])) @ Q @ self.__T2(self.plyangle[i]) # solve(T1(plyangle[i]), Q) @ T2(plyangle[i])
            A += Qbar*(self.z[i+1]-self.z[i])
            # coupling  stiffness
            B += (1/2)*Qbar*(self.z[i+1]**2-self.z[i]**2)
            # bending or flexural laminate stiffness relating moments to curvatures
            D += (1/3)*Qbar*(self.z[i+1]**3-self.z[i]**3)

        # laminate stiffness matrix
        ABD = zeros((6,6))
        ABD[0:3,0:3] = A
        ABD[0:3,3:6] = B
        ABD[3:6,0:3] = B
        ABD[3:6,3:6] = D
        self.ABD = ABD


    # method
    def available_materials(self):
        '''show the materials available in the library'''
        matprops = pd.read_csv(os.path.join(os.path.dirname(__file__), "compositematerials.csv"), index_col=0)
        print('---available materials---')
        for k in matprops.columns.tolist():
            print(k)
        print('-------------------------')

    # private method to be used internally
    def __import_matprops(self, mymaterial=['T300_5208','AL_7075']):
        '''
        import material properties
        '''
        matprops = pd.read_csv(os.path.join(os.path.dirname(__file__), "compositematerials.csv"), index_col=0)

        if mymaterial==[] or mymaterial=='':
            print(matprops.columns.tolist())

        mat = matprops[mymaterial]
        #mat.applymap(lambda x:np.float(x))
        mat = mat.applymap(lambda x:pd.to_numeric(x, errors='ignore'))
        return mat




def failure_envelope_laminate(Nx,Ny,Nxy,Mx,My,Mxy,q0,mymat,layup):
    '''
    find the miniumu margin give load conditions
    '''
    # create a 45 carbon cloth panel with a 0.5 inch rohacell core
    
    _, FAILUREINDEX_MAXSTRESS_max = laminate_calcs(NM=[Nx,Ny,Nxy,Mx,My,Mxy],
                             ek=[0,0,0,0,0,0],
                             q0=q0,
                             plyangle=   layup,
                             plymatindex=[0,0,0,0],
                             materials = [mymat],
                             platedim=[10,10],
                             zoffset=0,
                             SF=1.0,
                             plots=0,
                             prints=0)
    return FAILUREINDEX_MAXSTRESS_max


def plot_single_max_failure_loads(mymat='E-Glass Epoxy fabric M10E-3783', mylayup=[0,45,45,0] ):
    '''
    loops through and tries to find a load that is close to 0 and then 
    attempts to find the root (ie margin=0)
    
    older version used newton method for root finding
    scipy.optimize.newton(laminate_min, guess)
    
    TODO: Current calculation is stupid using random points to plot. fix it 
            by use FI, failure index instead of margin to generate a 
            linear relationship and envelope    
    '''
    #laminate_min = lambda N: failure_envelope_laminate(N,0,0,0,0,0,0)
    
    loadnamelist = ['Nx','Ny','Nxy','Mx','My','Mxy','q0']
    laminate_min_list = []
    
    laminate_min_list.append(lambda N: failure_envelope_laminate(N,0,0,0,0,0,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,N,0,0,0,0,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,0,N,0,0,0,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,0,0,N,0,0,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,0,0,0,N,0,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,0,0,0,0,N,0,mymat,mylayup))
    laminate_min_list.append(lambda N: failure_envelope_laminate(0,0,0,0,0,0,N,mymat,mylayup))
    
    envelope_loads = []
    N_t = array([0,1])
    N_c = array([0,-1])
    
    for loadname,laminate_min in zip(loadnamelist,laminate_min_list):

        # tension
        FI = [laminate_min(N) for N in N_t]
        m = (FI[1]-FI[0]) / (N_t[1] - N_t[0])
        b = FI[1]-m*N_t[1]
        N_crit_t = (1-b) / m
                 
        # compression
        FI = [laminate_min(N) for N in N_c]
        m = (FI[1]-FI[0]) / (N_c[1] - N_c[0])
        b = FI[1]-m*N_c[1]
        N_crit_c = (1-b) / m                    
                  
        envelope_loads.append('{} = {:.1f} , {:.1f}'.format(loadname,N_crit_t, N_crit_c))

    print('------------- enveloped loads for {} {} -----------------'.format(mylayup, mymat))
    for k in envelope_loads:
        print(k)
    
    # plot envelope
    Nx_env = []
    Nxy_env = []
    laminate_min = lambda N: failure_envelope_laminate(N,0,0,0,0,0,0,mymat,mylayup)
    # compression
    FI = [laminate_min(N) for N in N_c]
    m = (FI[1]-FI[0]) / (N_c[1] - N_c[0])
    b = FI[1]-m*N_c[1]
    Nx_env.append( (1-b) / m  )
    Nxy_env.append( 0  )
    # tension
    FI = [laminate_min(N) for N in N_t]
    m = (FI[1]-FI[0]) / (N_t[1] - N_t[0])
    b = FI[1]-m*N_t[1]
    Nx_env.append( (1-b) / m )
    Nxy_env.append( 0  )
    
    
    laminate_min = lambda N: failure_envelope_laminate(0,0,N,0,0,0,0,mymat,mylayup)
    # compression
    FI = [laminate_min(N) for N in N_c]
    m = (FI[1]-FI[0]) / (N_c[1] - N_c[0])
    b = FI[1]-m*N_c[1]
    Nxy_env.append( (1-b) / m  )   
    Nx_env.append( 0  )
    # tension
    FI = [laminate_min(N) for N in N_t]
    m = (FI[1]-FI[0]) / (N_t[1] - N_t[0])
    b = FI[1]-m*N_t[1]
    Nxy_env.append( (1-b) / m )   
    Nx_env.append( 0  )

    laminate_min_Nx_Nxy_func = lambda Nx,Nxy: failure_envelope_laminate(Nx,0,Nxy,0,0,0,0,mymat,mylayup)

    n = 500
    f = 1.25   # < 1
#    arr1 = np.random.randint(Nx_env[0]-abs(Nx_env[0]*f),Nx_env[0]+abs(Nx_env[0])*f,n)
#    arr2 = np.random.randint(Nx_env[1]-abs(Nx_env[1]*f),Nx_env[1]+abs(Nx_env[1])*f,n)
#    Nx_r = np.concatenate((arr1, arr2))
#    
#    arr1 = np.random.randint(Nxy_env[2]-abs(Nxy_env[2])*f,Nxy_env[2]+abs(Nxy_env[2])*f,n)
#    arr2 = np.random.randint(Nxy_env[3]-abs(Nxy_env[3])*f,Nxy_env[3]+abs(Nxy_env[3])*f,n)
#    Nxy_r = np.concatenate((arr1, arr2))    
    
    Nx_r = np.random.randint(Nx_env[0]*f,Nx_env[1]*f, n)
    Nxy_r = np.random.randint(Nxy_env[2]*f,Nxy_env[3]*f, n)
    for Nx_ri, Nxy_ri in zip(Nx_r, Nxy_r):
        FI = laminate_min_Nx_Nxy_func(Nx_ri, Nxy_ri)
        if FI < 1:
            Nx_env.append(Nx_ri)
            Nxy_env.append(Nxy_ri)
    
    points = array([ [x,xy] for x,xy in zip(Nx_env, Nxy_env)])
    
    hull = scipy.spatial.ConvexHull(points)
    plot(points[:,0], points[:,1], 'bo')
    for simplex in hull.simplices:
        plot(points[simplex, 0], points[simplex, 1], 'k-') 
    xlabel('Nx, lb/in')
    ylabel('Nxy, lb/in')
    title('Failure envelope')
    
    return envelope_loads


def my_laminate_with_loading():
    # loads lbs/in
    Nx  = 50
    Ny  = 0
    Nxy = 0
    Mx  = 0
    My  = 0
    Mxy = 0
    q0 =  0 # pressure
    # Qx = 0
    # Qy = 0
    a_width = 50
    b_length = 3.14*6.75
    
    ## sandwich laminate
    # plyangle=   [45,45,0, 45,45],
    # plymatindex=[0, 0, 1, 0, 0],    
    
    # create a 45 carbon cloth panel with a 0.5 inch rohacell core
    laminate_calcs(NM=[Nx,Ny,Nxy,Mx,My,Mxy],
             ek=[0,0,0,0,0,0],
             q0=q0,
             plyangle=   [0,60,-60,-60,60,0],
             plymatindex=[0,0,0,0,0,0],
             materials = ['E-Glass Epoxy Uni'],
             platedim=[a_width,b_length],
             zoffset=0,
             SF=2.0,
             plots=0,
             prints=1)

if __name__=='__main__':

    
    #plot_single_max_failure_loads()
    #plot_failure_index()
    my_laminate_with_loading()
    #material_plots(['E-Glass Epoxy fabric M10E-3783'])
    #plate()

    #plot_Nx_Nxy_failure_envelope(['Carbon_cloth_AGP3705H'])
    #plot_single_max_failure_loads()
    
#    # reload modules
#    import importlib ; importlib.reload
#    from composites import laminate
#    plyangle = [0,45]
#    matindex = [0,0]
#    matname = ['graphite-polymer_SI']
#    lam1 = laminate(plyangle, matindex, matname)
#    lam1.ABD
