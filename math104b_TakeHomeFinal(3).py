import numpy as np
import matplotlib.pyplot as plt

# problem (3)

def f(x):
    return float(-np.cos(x)/np.sin(1))*(3 - np.sin(np.pi*x))

def u(x):
    return float(-np.cos(x)/np.sin(1))

def FiniteDiffMethod(n,f,u):
    A = np.zeros((n+1,n+1))
    x0 = np.zeros((n+1))
    h = float(1)/n
    x = [i*h for i in range(0,n+1)]
    U = np.array(map(u,x))
    for i in range(0,n):
        A[i,i] = 2
        A[i+1,i] = -1
        A[i,i+1] = -1
    A[n,n] = 1
    A[0,0] = 1
    A*= 1/h**2
    for i in range(0,n+1):
       A[i,i] += float(2) - np.sin(np.pi*x[i])
    F = np.array(map(f,x))
    F[n] += 1/h
    return A,F,x,x0,U
   
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

errorList2 = []
klist2 = []
for j in range(6):
    AA,FF,xval2,xx0,newU2= FiniteDiffMethod(10*2**j,f,u)
    kk2,y2 = GradientMethodIt(AA,xx0,FF)
    GMerror = max(abs(y2 - newU2))/max(abs(newU2))
    errorList2.append(GMerror)
    klist2.append(kk2)
    print("k = ",10*2**j, ":iteration =", kk2, GMerror)

newError2 = np.array(map(np.log,errorList2))
newK2 = np.array(map(np.log,klist2))

plt.plot(newK2,newError2)
plt.show()
alpha = [(np.log(errorList2[i]/errorList2[i+1]))/np.log(2) for i in range(0,len(errorList2)-1)]
