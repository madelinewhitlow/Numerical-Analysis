#Madeline Whitlow
#Homework #4
#MATH 104B
#6344816

import time
import numpy as np
import copy
import matplotlib.pyplot as plt

# problem (1) part (a)
# f(x) as defined in part (a)
def f(x):
    return float((np.pi**2))*np.sin(np.pi*x)*np.cosh(np.sin(np.pi*x))

# u(x) as defined in part (a)
def u(x):
    return np.sinh(np.sin(np.pi*x))

# problem (1) part (b) the function FiniteDiffMethod uses discretization and
# the second order central finite differences to approximate the second derivative 
# which outputs a tridiagonal system 

def FiniteDiffMethod(n,f):
    A = np.zeros((n,n))
    F = np.zeros((n+1,1))
    h = float(1)/(n+1)
    x = [i*h for i in range(0,n+2)]
    for i in range(0,n-1):
        A[i,i] = 2
        A[i+1,i] = -1
        A[i,i+1] = -1
    A[n-1,n-1] = 2
    A*=1/h**2
    for i in range(0,n):
        A[i,i] += np.pi**2*(np.cos(np.pi*x[i+1]))**2
    F = np.array(map(f,x[1:n+1]))
    return A,F,x

# problem (1) part (c)
# gaussian elimination and backward substitution functions from homework (3)

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
    x = np.zeros((n,1))
    x[n-1] = b[n-1]/A[n-1,n-1]
    for j in range(n-2,-1,-1):
        x[j]=b[j]
        for k in range(j+1,n):
            x[j] = x[j] - A[j,k]*x[k]
        x[j] = x[j]/A[j,j]
    return x

# the TriDiag functions takes lists which represent the 3 diagonals of the 
# tridiagonal matrix and then solves the system of equations

def TriDiag(a,b,c,v):
    # a and c should have same length and b should be 1 more
    #length of v shoul be length of b
    a.insert(0,0)
    c.append(0)
    
    m = len(b) # result should be a list of length m
    
    cp = []
    vp = []
    cp.append(c[0]/b[0])
    vp.append(v[0]/b[0])

    for j in range(1,m-1):
        cp.append(c[j]/(b[j]-a[j]*cp[j-1]))
    
    for j in range(1,m):
        vp.append((v[j]-a[j]*vp[j-1])/(b[j]-a[j]*cp[j-1]))
    
    z = [0]*m
    z[m-1]=vp[m-1]
    for j in range(m-2,-1,-1):
        z[j] = vp[j]-cp[j]*z[j+1]
    return z

# timefunction uses GaussianElim and Backward sub to calculate how long the 
# algorithm takes for each value of n

def timefunction():
    timeSolve=[]
    for k in range(7):
        A,F,x = FiniteDiffMethod(10*2**k,f)
        tic = time.time()
        C,d = GaussianElim(A,F)
        E = BackwardSub(C,d)
        toc = time.time() 
        timeSolve.append(abs(tic-toc))
    return timeSolve

# timefunction2 uses TriDiag to calculate how long the 
# algorithm takes for each value of n

def timefunction2():
    timeSolve=[]
    for k in range(7):
        A,F,x = FiniteDiffMethod(10*2**k,f)
        a = [A[i+1,i] for i in range(0,10*2**k-1)]
        b = [A[i,i] for i in range(0,10*2**k)]
        c = [A[i,i+1] for i in range(0,10*2**k-1)]
        tic = time.time()
        L = TriDiag(a,b,c,F)
        toc = time.time() 
        timeSolve.append(abs(tic-toc))
    return timeSolve

#creating a loglog plot to compare the numerical complexity of each algorithm

n = [10*2**i for i in range(0,7)]
L = list(map(np.log,n))
timevalsSolve= timefunction()
timevalsSolve2 = timefunction2()
SolveListvals = list(map(np.log,timevalsSolve))
SolveListvals2 = list(map(np.log,timevalsSolve2))

# green is using tridiag and blue is regular gaussian
# green is faster for all values and blue is slowest for all values
# the slope of green is approximately .88 
# the slope of blue is approximately 1.777

plt.plot(L,SolveListvals)
plt.plot(L,SolveListvals2)
plt.show()



errorlist = []
h = [float(1)/(10*2**m+1) for m in range(6)]

for k in range(6):
    A,F,x= FiniteDiffMethod(10*2**k,f)
    a = [A[i+1,i] for i in range(0,10*2**k-1)]
    b = [A[i,i] for i in range(0,10*2**k)]
    c = [A[i,i+1] for i in range(0,10*2**k-1)]
    L = TriDiag(a,b,c,F)
    newx = x[1:10*2**k+1]
    U = np.array(map(u,newx))
    error = max(abs(U-L))
    errorlist.append(error)
    #print("k = ",10*2**k, ":solution error =", error)

Lerror = list(map(np.log,h))
ErrorListvals = list(map(np.log,errorlist))   
plt.plot(Lerror,ErrorListvals)
plt.show()


alpha = [(np.log(errorlist[i]/errorlist[i+1]))/np.log(2) for i in range(0,len(ErrorListvals)-1)]
print('alpha = order of accuracy for Finite diff method and TriDiag solve approximation',alpha)