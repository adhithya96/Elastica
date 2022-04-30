# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 18:42:05 2022

@author: Utilisateur
"""

# Construction of a beam fem code
import numpy as np
import matplotlib.pyplot as plt

# Geometric parameters
L = 1  # beam length (m)
Aire = 0.0001  # beam section (m^2)
# material parameters
E = 22000000000  # Young modulus (Pa, N/m^2)
I=50000


#Select the dimensionv of the geometric space

#dimension: int = 1
dimension: int = 2
#dimension: int = 3

# EF parameters
ne= 5  # nb of elements
  

        #dof=2 # degree of freedom
if (dimension==2):
    ddlpnode=3
elif (dimension==1):
    ddlpnode=1
else:   
    ddlpnode=6
N=ne*ddlpnode
#nodes = N + 1  # nb of nodes
Le = L / N  # Element length
if (dimension==2):
    # elementary stiffness matrix (E*I / L**3) 
    K_el =  np.array([[12, 6*Le, -12, 6*Le], [6*Le, 4*Le**2, -6*Le, 2*Le**2], [-12, -6*Le, 12, -6*Le], [6*Le, 2*Le**2, -6*Le, 4*Le**2]])
    # Global stiffness matrix
    K = np.zeros((N, N))
    for i in range(4):
        for j in range(4):
            K[i, j] += K_el[i, j]
    for p in range(N - 5):
        for k in range(N - 5):
            K[p + 4, k + 4] += K[p, k]
    K[N - 1, N - 1] = K_el[1, 1]
    #Solution directe
    # U_ef = np.zeros(N)
    # K_inv = np.linalg.inv(K[1:N, 1:N])# U_ef[1:N] = np.dot(K_inv, F[1:N]) #K*U=F==> U=inv(K)*F  
    #Incremental oads
    #Now i will try to apply load incrementally 
    loads=np.arange(0,1000,50)
    Y_last_node=[]

    for f in loads: 
        #f = 1  # applied force at the beam's end ( N)
        F = np.zeros(N)
        F[N - 1] = f
        m=0
        F[N - 2] = m
        #k=5*N
        #def conjugate_gradient(): iterative solution with residual norm
        crit=1e-12 #critère de convergeance sous forme d'erreur
        U_iter=np.zeros(N)  #( Y1, Theta1,Y2, Theta2,....)
        r = np.dot(K, U_iter) - F
        p = - r #direction de recherche
        r_k_norm = np.dot(r, r) #erreur
        iter=0
        max_iter=20
        while r_k_norm > crit and iter<max_iter:
            iter+=1
            Ap = np.dot(K, p)
            alpha = r_k_norm / np.dot(p, Ap) # calcul du pas optimal
            U_iter += alpha * p
            r += alpha * Ap
            r_kplus1_norm = np.dot(r, r)
            beta = r_kplus1_norm / r_k_norm # reorthoganalisation
            r_k_norm = r_kplus1_norm
            if r_kplus1_norm < crit:
                print ("The solution converged after", iter,"iterations")  #("The solution converged after", i,"iterations")
                break
            if r_k_norm > crit and iter>max_iter:
                print("No convergence and max iter reached ")
            p = beta * p - r
        print("\n The iterative solution is U_iter = \n",U_iter)
        Y_last_node.append(U_iter[-1])
        #return U_iter
    # Plot the solution
    x = np.arange(0, L + L / (N - 1), L / (N - 1))  # Length of the beam
    Y=[]
    Theta=[]
    Y=[U_iter[i] for i in range(len(U_iter)) if i%2==0]
    Theta=[U_iter[i] for i in range(len(U_iter)) if i%2!=0]

    print("Vertical displacement and rotation are")
    if len(Y)==len(Theta):
        for i in range(len(Y) or len(Theta)):
            print(Y[i],Theta[i])

    plt.plot(Y_last_node,loads)
    plt.grid()
    plt.title(" Displacement vs loads") 
    plt.xlabel("Displacement at the last node")
    plt.ylabel("Force")
    plt.show()
            



