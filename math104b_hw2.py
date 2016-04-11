#math104b 
#homework 2
#madeline whitlow

import math
import numpy as np

#problem (3) part (a) Hermite function is a quadrature rule which
#returns a function p(x) the hermite interpolation polynomial and 
#takes values fa = f(a) fb = f(b) fpa = f'(a) fpb = f'(b) and endpoints a,b


# function constructed for problem (3) part (d) to verify quadrature rule
def g(x):
    return math.cosh(x)

def gprime(x):
    return math.sinh(x)

def Hermite(fa,fpa,fb,fpb,a,b):
    h = float(b-a)
    alpha0 = fa
    alpha1 = fpa
    alpha2 = (fb - alpha0 - alpha1*h)/h**2
    alpha3 = (fpb -alpha1 - 2*alpha2*h)/h**2
    def P(x):
        return alpha0 + alpha1*(x-a) + alpha2*(x-a)**2 + alpha3*(x-b)*(x-a)**2
    return P

#constructed HermiteCoeffs to return a list of the coeffecient for
# the function p(x) 

def HermiteCoeffs(fa,fpa,fb,fpb,a,b):
    h = float(b-a)
    alpha0 = fa
    alpha1 = fpa
    alpha2 = (fb - alpha0 - alpha1*h)/h**2
    alpha3 = (fpb -alpha1 - 2*alpha2*h)/h**2
    return [alpha0,alpha1,alpha2,alpha3]

#constructed HermiteInt to integrate p(x) which calls the list of
#coefficients from HermiteCoeffs

def HermiteInt(fa,fpa,fb,fpb,a,b):
    p = HermiteCoeffs(fa,fpa,fb,fpb,a,b)
    intp = p[0]*(b-a) + p[1]/float(2)*(b-a)**2 + p[2]/float(3)*(b-a)**3 - p[3]/float(12)*(b-a)**4
    return intp

# constructed HermiteCompInt to cut the intervals into smaller chunks
#in order to verify the order of accuracy of the quadrature

def HermiteCompInt(f,fp,a,b,n):
    h = float(b-a)/n
    x = [a+i*h for i in range(0,n+1)]
    fvals = list(map(f,x))
    fpvals = list(map(fp,x))
    sum1 = 0
    for j in range(0,n):
        sum1+= HermiteInt(fvals[j],fpvals[j],fvals[j+1],fpvals[j+1],x[j],x[j+1])       
    return sum1 

a = 0
b = 1

#actual value that quadrature should converge towards
RealVal = math.sinh(1)

#values of hermite quadrature
IvalsHermite = [HermiteCompInt(g,gprime,a,b,2**i) for i in range (0,10)]

#error between approximated values and real values
ErrorHermite = [abs(RealVal - IvalsHermite[i]) for i in range(0,len(IvalsHermite))]

#alpha value that shoes the order of the approximation of error
alpha = [(math.log(ErrorHermite[i]/ErrorHermite[i+1]))/math.log(2) for i in range(0,len(IvalsHermite)-1)]



#problem (4) Gaussian quadrature for integrals from [-1,1]

def f(x):
    if x == 0:
        return 0
    else:
        return math.sin(x)/x

def Gaussian(f):
    R = float(5)/9*f(-math.sqrt(3/5)) + float(8)/9*f(0) + float(5)/9*f(math.sqrt(3/5))
    return R

def CompGaussian(f):
    h = float(1)/2
    list1=[]
    R = float(5)/9*f(-math.sqrt(float(3)/5)) + float(5)/9*f(math.sqrt(float(3)/5))
    T = float(5)/9*f(-math.sqrt(3/5)*h + h) + float(8)/9*f(h) + float(5)/9*f(math.sqrt(3/5)*h + h)
    list1.append(R)
    list1.append(T)
    return list1

SineReal= 0.9460830703671830149413533138231796578123379547381117

#values of gaussian quadrature for sin(x)/x part (c)
IvalsGaussian = CompGaussian(f)

#error between approximated values and real values
ErrorGaussian = [abs(SineReal - IvalsGaussian[i]) for i in range(0,len(IvalsGaussian))]

#alpha value that shows the order of accuracy in quadrature approx
alphaGaussian = [(math.log(ErrorGaussian[i]/ErrorGaussian[i+1]))/math.log(2) for i in range(0,len(IvalsGaussian)-1)]

#function defined to output nth taylor term for sin(x)/(x)
def TaylorTerm(n):
    return float(1)/((2*n+3)*math.factorial(2*n+3))    

#tells what n we need to get 10 digits of accuracy (n=5)
def TaylorAccuracy(tolerance):
    n = 0
    while TaylorTerm(n) > tolerance:
        n+=1
    return n

def TaylorForSinIntegral(n):
    Q = 0
    for i in range(0,n+1):
        Q+= float((-1)**i)/((2*i+1)*math.factorial(2*i+1)) 
    return Q