#Madeline Whitlow
#math 104b
#take home final

import numpy as np
import matplotlib.pyplot as plt
import copy

def f(x):
    return float(-np.cos(x)/np.sin(1))*(3 - np.sin(np.pi*x))

def u(x):
    return float(-np.cos(x)/np.sin(1))


def FiniteDiffMethod(n,f,u):
    A = np.zeros((n+1,n+1))
    x0 = np.ones((n+1))
    h = float(1)/n
    x = [i*h for i in range(0,n+1)]
    for i in range(0,n):
        A[i,i] = 2
        A[i+1,i] = -1
        A[i,i+1] = -1
    A[n,n] = 2
    A[0,1] = -2
    A[n,n-1] = - 2
    A*= 1/h**2
    for i in range(1,n+1):
       A[i,i] += float(2) - np.sin(np.pi*x[i])
    A[0,0] += float(2)
    U = np.array(map(u,x))
    F = np.array(map(f,x))
    F[n] += 2/h
    return A,F,x,x0,U

def GaussianElim(A,b):
    n = len(A)
    newA = copy.deepcopy(A)
    newb = copy.deepcopy(b)
    for i in range(0,n-1):
        for j in range(i+1,n):
            C = (newA[j,i])/float(newA[i,i])
            newA[j,i:] = newA[j,i:] - C*newA[i,i:]
            newA[j,i] = C
            newb[j] = newb[j] - C*newb[i]
    return newA, newb

def BackwardSub(A,b):
    n = len(b) 
    x = np.zeros(n)
    x[n-1] = b[n-1]/A[n-1,n-1]
    for j in range(n-2,-1,-1):
        x[j]=b[j]
        for k in range(j+1,n):
            x[j] = x[j] - A[j,k]*x[k]
        x[j] = x[j]/A[j,j]
    return x


hlist = []
errorlistGE = []
for k in range(6):
    AA,FF,xx,xxstart,ux = FiniteDiffMethod(10*2**k,f,u)
    CC,dd = GaussianElim(AA,FF)
    E = BackwardSub(CC,dd)
    norm = max(abs(E-ux))/max(abs(ux))
    errorlistGE.append(norm)
    hlist.append(10*2**k)
    print("k = ",10*2**k, ":error =", norm)

errorlistGE2 = np.array(map(np.log,errorlistGE))
hlist2 = np.array(map(np.log,hlist))
plt.plot(hlist2,errorlistGE2)
plt.show()

alpha = [(np.log(errorlistGE[i]/errorlistGE[i+1]))/np.log(2) for i in range(0,len(errorlistGE)-1)]
print alpha

def JacobiItMethod(A,f,x):
    n = len(x)
    p = np.zeros((n))
    q = np.zeros((n))
    j = 0
    while max(abs(x-q)) > 10**(-5):
        j+=1
        q = copy.deepcopy(x)
        p[0] = (f[0] - A[0][1]*x[1])/A[0][0]
        p[n-1] = (f[n-1] - A[n-1][n-2]*x[n-2])/A[n-1][n-1]
        for i in range(1,n-1):
            p[i] = (f[i] - (A[i][i-1]*x[i-1] + A[i][i+1]*x[i+1]))/A[i][i]
        for i in range(0,n):
            x[i] = p[i]
    return j,x


jlist = []
errorlistJ = []
for k in range(6):
    A,F,x,x0,U= FiniteDiffMethod(10*2**k,f,u)
    j,xnew1 = JacobiItMethod(A,F,x0)
    jlist.append(j)
    print("k = ",10*2**k, ":iteration =", j)

N = [10*2**k for k in range(0,6)]
plt.plot(N,jlist)
plt.show()

def GaussSeidelMethod(A,f,x):
    j = 0
    n = len(x)
    p = np.zeros(n)
    while max(abs(x-p)) > 10**(-5):
        j+=1
        p = copy.deepcopy(x)
        x[0] = (f[0] - A[0][1]*x[1])/A[0][0]
        for i in range(1,n-1):
            x[i] = (f[i] - (A[i][i-1]*x[i-1] + A[i][i+1]*x[i+1]))/A[i][i] 
        x[n-1] = (f[n-1] - A[n-1][n-2]*x[n-2])/A[n-1][n-1]
    return j,x




klistGS = []
errorlistGS = []
for k in range(6):
    AGS,fGS,xGS,xstartGS,uGS= FiniteDiffMethod(10*2**k,f,u)
    kGS , xnew2 = GaussSeidelMethod(AGS,fGS,xstartGS)
    klistGS.append(kGS)
    print("k = ",10*2**k, ":iteration =", kGS)

plt.plot(N,klistGS)
plt.show()

