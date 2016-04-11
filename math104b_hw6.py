#Madeline Whitlow
#Math 104B 
#Homework 6

import numpy as np
import matplotlib.pyplot as plt

#this function was originally written to generate A in a different homework
# and has been revised to generate A, f, u based on n

def FiniteDiffMethod(n):
    A = np.zeros((n,n))
    h = float(1)/(n+1)
    u = np.array([np.sin((np.pi*i)/(n+1)) for i in range(1,n+1)])
    for i in range(0,n-1):
        A[i,i] = 2
        A[i+1,i] = -1
        A[i,i+1] = -1
    A[n-1,n-1] = 2
    A*=1/h**2
    f = np.dot(A,u)
    return A,f,u

#NOTE: it is important to note that both the JacobiItMethod and GaussSeidelMethod
# only work for tridiagnolalized systems

# problem (1) part (b) asks to write the jacobi iteration for system (3) 
# problem (1) part (f) askes to solve system using Jacobi method


def JacobiItMethod(A,f,x,u,n,eps):
    p = np.zeros((n))
    norm = eps*np.linalg.norm(u)
    j = 0
    while np.linalg.norm(x-u) > norm:
    #for j in range(0,M):
        j+=1
        p[0] = (f[0] - A[0][1]*x[1])/A[0][0]
        p[n-1] = (f[n-1] - A[n-1][n-2]*x[n-2])/A[n-1][n-1]
        for i in range(1,n-1):
            p[i] = (f[i] - (A[i][i-1]*x[i-1] + A[i][i+1]*x[i+1]))/A[i][i]
        for i in range(0,n):
            x[i] = p[i]
    return j

eps = 10**(-3)
#jlist = []
#for k in range(6):
#    A,f,u= FiniteDiffMethod(10*2**k)
#    x = np.ones((10*2**k))
#    j = JacobiItMethod(A,f,x,u,10*2**k,eps)
#    jlist.append(j)
#    print("k = ",10*2**k, ":iteration =", j)
#
#
n = [10*2**k for k in range(0,6)]
L = list(map(np.log,n))
#listvals = list(map(np.log,jlist)) 
#slope = np.polyfit(L,listvals,1)
#print(slope)  
#plt.plot(L,listvals)
#plt.show()


# problem (1) part (g) askes to solve system using Gauss-Seidel method

def GaussSeidelMethod(A,f,x,u,n,eps):
    norm = eps*np.linalg.norm(u)
    j = 0
    while np.linalg.norm(x-u) > norm:
        j+=1
        x[0] = (f[0] - A[0][1]*x[1])/A[0][0]
        for i in range(1,n-1):
            x[i] = (f[i] - (A[i][i-1]*x[i-1] + A[i][i+1]*x[i+1]))/A[i][i] 
        x[n-1] = (f[n-1] - A[n-1][n-2]*x[n-2])/A[n-1][n-1]
    return j

klist = []
for k in range(6):
    A,f,u= FiniteDiffMethod(10*2**k)
    x = np.ones((10*2**k))
    kk = GaussSeidelMethod(A,f,x,u,10*2**k,eps)
    klist.append(kk)
    print("k = ",10*2**k, ":iteration =", kk)
 
listvals2 = list(map(np.log,klist))   
slope2 = np.polyfit(L,listvals2,1)
print(slope2) 
plt.plot(L,listvals2)
plt.show()