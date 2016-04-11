#Math 104B
#Homework 8
#madeline Whitlow
#6344816


import numpy as np
import copy
import matplotlib.pyplot as plt


# problem (1) part (a) f_test is the constructed function for which we know 
# the exact analytic solution
def f_test(x):
    return (np.pi**2)*np.sin(np.pi*x) + np.cosh(np.pi*x)**2*np.sin(np.pi*x)

#this is the exact analytic solution to the problem with f_test

def u_test(x):
    return np.sin(np.pi*x)
    
#this is the actual f(x) that we want to find u(x) for

def f_actual(x):
    return np.sin(np.pi*x)
# We construct our matrix A using the discretization technique as well as F which
# is the right hand side of Ax = b, x which are the x_i points and x0 which is the 
# intital guess 

def FiniteDiffMethod(n,f):
    A = np.zeros((n,n))
    F = np.zeros((n))
    x0 = np.ones(n)
    h = float(1)/(n+1)
    x = np.array([i*h for i in range(1,n+1)])
    for i in range(0,n-1):
        A[i,i] = 2
        A[i+1,i] = -1
        A[i,i+1] = -1
    A[n-1,n-1] = 2
    A*=1/h**2
    for i in range(0,n):
        A[i,i] += (np.cosh(np.pi*x[i]))**2
    F = np.array(map(f,x))
    return A,F,x,x0

def steepestDescentIt(A,x,b):
    p = b - np.dot(A,x)
    k = 0
    while np.linalg.norm(p) > 10**(-10):
        p = b - np.dot(A,x)
        alpha = float(np.dot(p,p)/np.dot(np.dot(p,A),p))
        x += alpha*p
        k+=1
    return k,x


def GradientMethodIt(A,x0,b):
    x = x0
    r = b - np.dot(A,x)
    p = r
    k = 0
    while np.linalg.norm(r) > 10**(-10):
        alpha = float(np.dot(r.T,r)/np.dot(np.dot(p.T,A),p))
        x += alpha*p
        ri  = r - np.dot(alpha*A,p)
        beta = float(np.dot(ri.T,ri)/np.dot(r.T,r))
        p = ri + beta*p
        r = ri
        k+=1
    return k,x
    

kvals = [10*2**i for i in range(5)]
errorList = []
klist = []
#for k in range(5):
#    A,F,xval,x0= FiniteDiffMethod(10*2**k,f_test)
#    kk,y = steepestDescentIt(A,x0,F)
#    newU = np.array(map(u_test,xval))
#    SDerror = np.linalg.norm(y - newU)
#    errorList.append(SDerror)
#    klist.append(kk)
#    print("k = ",10*2**k, ":iteration =", kk, SDerror)

errorList2 = []
klist2 = []
for j in range(6):
    AA,FF,xval2,xx0= FiniteDiffMethod(10*2**j,f_test)
    kk2,y2 = GradientMethodIt(AA,xx0,FF)
    newU2 = np.array(map(u_test,xval2))
    GMerror = np.linalg.norm(y2 - newU2)
    errorList2.append(GMerror)
    klist2.append(kk2)
    print("k = ",10*2**j, ":iteration =", kk2, GMerror)

newError = np.array(map(np.log,errorList))
newError2 = np.array(map(np.log,errorList2))
newK = np.array(map(np.log,klist))
newK2 = np.array(map(np.log,klist2))

slope = np.polyfit(newK,newError,1)
print(slope)  
plt.plot(newK,newError)
plt.plot(newK2,newError2)
plt.show()




kvalsB = [10*2**i for i in range(6)]
errorListB = []
klistB = []
for k in range(6):
    AB,FB,xvalB,x0B= FiniteDiffMethod(10*2**k,f_actual)
    kkB,yB = steepestDescentIt(AB,x0B,FB)
    klist.append(kkB)
    norm = np.linalg.norm(yB)
    print("k = ",10*2**k, ":iteration =", kkB)

klist2C = []
for j in range(6):
    AAC,FFC,xval2C,xx0C= FiniteDiffMethod(10*2**j,f_actual)
    kk2C,y2C = GradientMethodIt(AAC,xx0C,FFC)
    klist2C.append(kk2C)
    print("k = ",10*2**j, ":iteration =", kk2C)


X = np.array([i/10.0 for i in range(0,11)])
Y = np.array([0.95161,1.12879,1.28257,1.58077,1.77259,2.15300,2.17509,2.29357,2.76034,2.92346,2.97703])

def LeastSquaresMatrices(X,Y,Deg):
    A = np.zeros((Deg+1,Deg+1))
    F = np.zeros((Deg+1))
    for m in range(0,Deg+1):
        for j in range(0, Deg+1):
            A[m,j] = sum(X**(j+m))
        F[m] = np.dot(X**m,Y)
    return A,F

def LeastSquaresSolver(X,Y,Deg):
    A,F = LeastSquaresMatrices(X,Y,Deg)
    return np.linalg.solve(A,F)

X1 = LeastSquaresSolver(X,Y,1)
X2 = LeastSquaresSolver(X,Y,2)

X1List = [X1[0] + X[i]*X1[1] for i in range(0,11)]
X2List = [X2[0] + X[i]*X2[1] + X[i]*X2[2]**2 for i in range(0,11)]

error1 = sum((X1List - Y)**2)
error2 = sum((X2List - Y)**2)

plt.plot(X,X1List)
plt.plot(X,X2List)
plt.scatter(X,Y)
plt.show()