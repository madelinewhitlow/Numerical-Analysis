#Madeline Whitlow
#Math 104B
#Homework 3

import time
import numpy as np
import scipy.linalg 
import matplotlib.pyplot as plt
import copy
import math

#problem (1) 

#RandAxb will generate A an n x n matrix, and x and b two n x 1 matrices which
# will verify that the functions GaussianElim and BackwardSub 
#generate the correct solution to the system of equations

def RandAxb(n):
    A = 4*np.identity(n) + np.random.rand(n,n)
    x = np.random.rand(n,1)
    b = np.dot(A,x)
    return A,x,b

#Gaussian Elimination function which takes inputs A and b and returns an
#upper triangular matrix with the multipliers in the lower triangular part

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

#accepts A and b with assumption that A is already upper traingular (from 
# GaussianElim)  and returns x the solution to the equation Ax = b

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

#problem (1) part (a) 
#the linalg.norm procedure verifies that the difference is zero
#which means that the computed values are accurate
#print('Problem (1) part (a)')
#for k in range(7):
#    A,x,b = RandAxb(10*2**k)
#    C,d = GaussianElim(A,b)
#    Solution = BackwardSub(C,d)
#    print("k = ",k, ":solution error =", np.linalg.norm(x-Solution))

#problem (1) part (b)
#constructA builds an nxn matrix but constructs x using np.linalg.solve function
def constructA(n):
    A = 4*np.identity(n) + np.random.rand(n,n)
    b = np.random.rand(n,1)
    x = np.linalg.solve(A,b)
    return A,x,b

#we know the code is working properly by applying the linalg norm again
# to the difference between the constructed value of x and the generated one
# and getting a result of 0 

#print('Problem (1) part (b)')

#for k in range(7):
#    Ab,xb,bb = constructA(10*2**k)
#    Cb, db = GaussianElim(Ab,bb)
#    newxb = BackwardSub(Cb,db)
#    print("k = ",k, ":solution error =", np.linalg.norm(xb-newxb))


#problem (1) part (c)
# we know the code is working properly because the difference between the generated
#value xc and the solution using the built in LU decomp U is zero

#print('problem (1) part (c)')
#
#for k in range(7):
#    Ac,xc,bc = constructA(10*2**k)
#    Cc, dc = GaussianElim(Ac,bc)
#    newxc = BackwardSub(Cc,dc)
#    F = scipy.linalg.lu_factor(Ac)
#    U = scipy.linalg.lu_solve(F,bc)
#    print("k = ",k, ":solution error =", np.linalg.norm(U-xc))

#problem (1) part (d)
# the function timepartD tracks the time for each computation of GaussianElim
# and for BackSub

def timepartD():
    time1=[]
    time2=[]
    for k in range(7):
        Ab,xb,bb = constructA(10*2**k)
        tic1=time.time()
        Cb, db = GaussianElim(Ab,bb)
        toc1=time.time()
        time1.append(abs(toc1-tic1))
        tic2=time.time()
        newxb = BackwardSub(Cb,db)
        toc2=time.time()
        time2.append(abs(toc2-tic2))
    return time1,time2

#print('problem (1) part (d)')
##creating lists of values for n and corresponding run times of functions
#timeGaussD,timeBackSubD=timepartD()
n = [10*2**i for i in range(0,7)]
L = list(map(math.log,n))
#GaussList = list(map(math.log,timeGaussD))
#BackSubList = list(map(math.log,timeBackSubD))
#print('GaussList of run times:',GaussList)
#print('BackSubList of run times:',BackSubList)
#print('Blue is GaussList and Green is BackSubList')
#f = plt.figure(1)
#plt.plot(L,GaussList)
#plt.plot(L,BackSubList)
#f.show()
#the scale of the plot is linear with a slope of 2

#problem (1) part (e)

def timefunction():
    timeSolve=[]
    for k in range(7):
        Ac,xc,bc = constructA(10*2**k)
        tic = time.time()
        Cc = np.linalg.solve(Ac,bc)
        toc = time.time() 
        timeSolve.append(abs(tic-toc))
    return timeSolve

print('problem (1) part (e)')
timevalsSolve= timefunction()
SolveListvals = list(map(math.log,timevalsSolve))
print('SolveListvals of run times:',SolveListvals)
t = plt.figure(6)
plt.plot(L,SolveListvals)
t.show()

#the curve of the timing is relatively linear with slope 1.5



def f(x):
   return 3*x**2 + 4 - 2*np.e**x - np.e**(1-x)
    
def K(x,t):
    return np.e**(abs(x-t))

def U(x):
    return x**2


def TrapezoidalRuleMatrix(f,K,m):
    h = float(1)/m
    x = [j*h for j in range(m+1)]
    b = -np.array(map(f,x))
    matrix = np.zeros((m+1,m+1))
    for i in range(0,m+1):
        for j in range(0,m+1):
            matrix[i,j] = K(x[i],x[j])
    matrix = h*matrix
    matrix[:,0]*=.5
    matrix[:,m]*=.5
    matrix = matrix - np.identity(m+1)
    return matrix,b,x
    
#print('Problem (2) part (b)')

    
#CompTrapezoidalRuleMatrix generates the matrix of K() values using the
#trapezoidal rule and then by using the Gaussian Elimination code and 
#Backward sub we can solve to approximate the integral in (1)

def CompTrapezoidalRuleMatrix(f,K,m):
    h = float(1)/m
    x = [j*h for j in range(m+1)]
    b = -np.array(map(f,x))
    matrix = np.zeros((m+1,m+1))
    for i in range(0,m+1):
        for j in range(0,m+1):
            matrix[i,j] = K(x[i],x[j])
    matrix2 = copy.deepcopy(matrix)
    matrix2[:,0]*=3
    matrix2[:,1]*=(-4)
    matrix2[:,m-1]*=(-4)
    matrix2[:,m]*=3
    matrix2[:,3:m-2]*=0
    matrix2*=h/24
    matrix = h*matrix
    matrix[:,0]*=.5
    matrix[:,m]*=.5
    matrix = matrix + matrix2 - np.identity(m+1)
    return matrix,b,x

errorlist = []
h = [float(1)/(10*2**m) for m in range(5)]
for k in range(5):
    matrix,b,x=TrapezoidalRuleMatrix(f,K,10*2**k)
    matrixA,mb = GaussianElim(matrix,b)
    g = BackwardSub(matrixA,mb) 
    u = np.array(map(U,x))
    error = max(abs(u-g.flatten()))
    errorlist.append(error)
    print("k = ",10*2**k, ":solution error =", error)

#p=plt.figure(3)
#Lerror = list(map(math.log,h))
#ErrorListvals = list(map(math.log,errorlist))   
#plt.plot(Lerror,ErrorListvals)
#p.show()
#
###alpha is the order of accuracy of the TrapezoidalApprox
#alpha = [(math.log(errorlist[i]/errorlist[i+1]))/math.log(2) for i in range(0,len(ErrorListvals)-1)]
#print('alpha = order of accuracy for trapezoidal approximation',alpha)

#the order of accuracy is 2


#print('Problem (2) part (c)')
#
#errorlist2 = []
#for k in range(5):
#    matrix2,b2,x2=CompTrapezoidalRuleMatrix(f,K,10*2**k)
#    matrixA2,mb2 = GaussianElim(matrix2,b2)
#    g2 = BackwardSub(matrixA2,mb2) 
#    u2 = np.array(map(U,x2))
#    error2 = max(abs(u2-g2.flatten()))
#    errorlist2.append(error2)
#    #print("k = ",10*2**k, ":solution error =", error2)
#y = plt.figure(2)
#Lerror2 = list(map(math.log,h))
#ErrorListvals2 = list(map(math.log,errorlist2))   
#plt.plot(Lerror2,ErrorListvals2)
#y.show()
#
#alpha2 = [(math.log(errorlist2[i]/errorlist2[i+1]))/math.log(2) for i in range(0,len(ErrorListvals2)-1)]
#print('alpha2 = order of accuracy for composite trapezoidal approximation',alpha2)

# the order of accuracy is again 2
