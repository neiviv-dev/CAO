# CHAMBE Vivien - CHEMIN Melvyn

from cmath import pi
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def sieve_of_eratosthene ():
    nmax = int(input("Entrez un Nmax: "))
    # Crée une liste d'entier de 1 à nmax en supprimant 0 et 1
    l = []
    for i in range (nmax):
        if i == 0 or i == 1:
            l.append("x")
        else:
            l.append(i)
    i = 0
    # Remplace les nombres non premiers par des "x"
    while i != nmax:
        if l[i] == "x":
            i+=1
        else:
            mult = l[i]
            j= i+1
            while j != nmax:
                if l[j] == "x":
                    j += 1
                elif l[j]% mult == 0:
                        l[j] = "x"
                        j+=1
                else: j+=1
            i+=1
    # Retire les x de la liste 
    primed = []
    for i in range (nmax):
        if l[i] != "x":
            primed.append(l[i])

    
    return primed

#print (sieve_of_eratosthene())

# Exercice 2

# M1
A = np.ones((2,4))
B = np.zeros((2,3))
M1 = np.hstack((A,B))

print("M1:")
print(M1)

# M2
A = np.linspace(1,9,5)
# Il est également possible d'utiliser np.arange 
A = np.arange(1,10,2)
B = np.linspace(8,0,5)
M2 = np.vstack((A,B,B))

print("M2:")
print(M2)

# M3
A = np.zeros((3,2))
B = np.ones((2,5))
C = np.linspace(1,5,5)

D = np.vstack((B,C))
M3 = np.hstack((A,D))

print("M3:")
print(M3)

# M4
A = np.ones(7)
B = np.ones(8)*4
B[0] = 2
B[-1] = 2
C = np.diag(A,1)
D = np.diag(A,-1)
E = np.diag(B)
M4 = C + D + E

print("M4:")
print(M4)


#### Exercice 3

fig  = plt.figure()
ax = plt.axes(projection='3d')

nbpoints = 240

u = np.linspace(0,2,nbpoints)
v = np.linspace(0,2*np.pi,nbpoints)

U,V = np.meshgrid(u,v)

def F1(U,V):
    return U*np.cos(V)
def F2(U,V):
    return U*np.sin(V)

def F3(U,V):
    return pow(U,3)*np.cos(3*V)

x = F1(U,V)
y = F2(U,V)
z = F3(U,V)

ax.plot_wireframe(x,y,z)


