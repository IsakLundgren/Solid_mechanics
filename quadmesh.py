import numpy as np
import matplotlib.pyplot as plt
import calfem.core as cco
#import calfem.vis as vis

def ex_ey_quadmesh(p1,p2,nelx,nely,ndofs):
    xv=np.linspace(p1[0],p2[0],nelx+1)
    yv=np.linspace(p1[1],p2[1],nely+1)

    nel  = nelx*nely;
    Ex=np.zeros((nel,4))
    Ey=np.zeros((nel,4))
    for m in range(0,nely):
        for n in range(0,nelx):
            Ex[n+m*nelx,0]=xv[n]
            Ex[n+m*nelx,1]=xv[n+1]
            Ex[n+m*nelx,2]=xv[n+1]
            Ex[n+m*nelx,3]=xv[n]
        #
            Ey[n+m*nelx,0]=yv[m]
            Ey[n+m*nelx,1]=yv[m]
            Ey[n+m*nelx,2]=yv[m+1]
            Ey[n+m*nelx,3]=yv[m+1]
            
    return Ex, Ey

def edof_quadmesh(nelx,nely,ndofs):
    Edof=np.zeros((nelx*nely,4*ndofs),'i')
    for m in range(0,nely):
        for n in range(0,nelx):
            Edof[n+m*nelx,0]=n*ndofs+1+m*(nelx+1)*ndofs
            Edof[n+m*nelx,1]=n*ndofs+2+m*(nelx+1)*ndofs
            Edof[n+m*nelx,2]=(n+1)*ndofs+1+m*(nelx+1)*ndofs
            Edof[n+m*nelx,3]=(n+1)*ndofs+2+m*(nelx+1)*ndofs
        #
            Edof[n+m*nelx,4]=(n+1)*ndofs+1+(m+1)*(nelx+1)*ndofs
            Edof[n+m*nelx,5]=(n+1)*ndofs+2+(m+1)*(nelx+1)*ndofs
            Edof[n+m*nelx,6]=n*ndofs+1+(m+1)*(nelx+1)*ndofs      
            Edof[n+m*nelx,7]=n*ndofs+2+(m+1)*(nelx+1)*ndofs
    return Edof

def B1B2B3B4_quadmesh(nelx,nely,ndofs):
    #lower boundary, dofs
    B1=np.linspace(1,(nelx+1)*ndofs,(nelx+1)*ndofs)
    B1=B1.astype(int)
    B2=np.zeros(((nely+1)*ndofs),'i')
    nn=0
    for n in range(0,nely+1):
        B2[nn]=(nelx+1)*ndofs*(n+1)-1
        if ndofs>1:
            B2[nn+1]=(nelx+1)*ndofs*(n+1)+0
        nn=nn+ndofs

    B3=np.linspace(1,(nelx+1)*ndofs,(nelx+1)*ndofs)+(nelx+1)*ndofs*nely
    B3=B3.astype(int)

    B4=np.zeros(((nely+1)*ndofs),'i')
    nn=0
    for n in range(0,nely+1):
        B4[nn]=(nelx+1)*ndofs*n+1
        if ndofs>1:
            B4[nn+1]=(nelx+1)*ndofs*n+2
        nn=nn+ndofs
   
    P1=np.zeros((2),'i'); P2=np.zeros((2),'i'); P3=np.zeros((2),'i'); P4=np.zeros((2),'i')
    for m in range(0,2):
        P1[m]=B1[m]
        P2[m]=B2[m]
        P4[m]=B3[m]
    P3[0]=B3[-1]-1
    P3[1]=B3[-1]
    return B1,B2,B3,B4,P1,P2,P3,P4

def quadmesh(p1,p2,nelx,nely,ndofs):
    Ex, Ey=ex_ey_quadmesh(p1,p2,nelx,nely,ndofs)
    Edof=edof_quadmesh(nelx,nely,ndofs)
    B1,B2,B3,B4,P1,P2,P3,P4=B1B2B3B4_quadmesh(nelx,nely,ndofs)
    return Ex,Ey,Edof,B1,B2,B3,B4,P1,P2,P3,P4

sigma_yield = 400e+9 #Pa
E = 200e+9 #Pa
nu = 0.3
rho = 7800 #kg/m^3
g = 9.81 #m/s^2
m = 100 #kg
h = 0.05 #m
b = 0.05 #m
A = b * h
D = cco.hooke(1, E, nu) #Stiffness matrix

p1 = [0,0] #m
p2 = [2,0.025] #m
thickness = 1 #m
nelx = int(2 / 0.025 * 1)
nely = 1
nel = nelx * nely
node_dofs = 2 #2 dimensional

Ex,Ey,Edof,B1,B2,B3,B4,P1,P2,P3,P4 = quadmesh(p1,p2,nelx,nely,node_dofs)

ndofs = np.max(Edof)
bc = np.array([B4, np.zeros(np.size(B4))]) #Clamp left hand side

K = np.zeros((ndofs, ndofs))
f = np.zeros(ndofs)
eq = np.zeros(ndofs)

#Set boundary forces
#for i in B3:
#    if (i - 1 % 2 == 0):
#        eq[i - 1] = -rho * g * A


print(Ex[0, :])
for el in range(nel):
    Ke, fe = cco.platre(Ex[el, :], Ey[el, :], [thickness], D, eq[el])
    K, f = cco.assem(Edof[el,:], K, Ke, f, fe)


    #done hehe