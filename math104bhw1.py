#Madeline whitlow
#Math 104B Homework 1

import math
import numpy as np

a = 0
b = 1
n = 1000
h = float(b-a)/n

def few(c):
    return c + np.e**(-c**2)

ureal = 1.635257257949994654568013599782201031869552084055950138982

# trapezoidal rule approximation built for problem (3) part (a) but also referenced in problem (1)
def TrapezoidalRule(f,a,b,n,usefirst):
    h = float(b-a)/n
    #if usefirst:
    s = f(a) + f(b)
    #else:
        #s = f(b)
    for i in range(1,n):
        s += 2*f(a+i*h)
    return(s*h/2)





#function using to test corrected trapezoidal rule in (1)
def g(u):
    return math.cos(u)
    
#corrected trapezoidal rule approximation using taylor's expansion
#to get approximations for function derivatives
def CorrectedTrapRule(f,a,b,n):
    h = float(b-a)/n
    FprimeofA = (4*f(a+h) - f(a + 2*h) -3*f(a))/(2*h)
    FprimeofB = (-4*f(b-h) + f(b - 2*h) + 3*f(b))/(2*h)
    return TrapezoidalRule(f,a,b,n,True)+(h**2/12)*(FprimeofB - FprimeofA)

#real value of interval trying to approximate in (1)
Inew=math.sin(1)

#values of corrected trapezoidal rule 
IvalsTrapCorrected = [CorrectedTrapRule(g,a,b,2**i) for i in range (0,20)]

#error between approximated values and real values
ErrorTrapCorrected = [abs(Inew - IvalsTrapCorrected[i]) for i in range(0,len(IvalsTrapCorrected))]

#alpha value that shoes the order of the approximation of error
alphaWithCorrected = [math.log((ErrorTrapCorrected[i]/ErrorTrapCorrected[i+1])*4)/math.log(2) for i in range(0,len(IvalsTrapCorrected)-1)]

#problem (2)
#function trying to approximate in part (c)
def newg(u):
    return 1/(u+1)

#function in part (d) trying to approximate
def partdFunction(u):
    return math.sqrt(u)

#approximation for (2) part (d) using quadrature formula 
#derived in problem (2) part (b) on paper

#unsure of what problem asks, so approximating function newg() by corrected trap
def Ifunction(f,a,b,n):
   # h = float(b-a)/n
    x = float(b-a)/(90)*(7*f(a)+32*f((b+3*a)/float(4))+12*f((a+b)/float(2))+32*f((a+3*b)/float(4))+7*f(b))
    return x

#real value of function trying to approximate (ln(2)) in part (c)
realI = math.log(2)

#run of numbers approximating function
IvalsApprox = [CorrectedTrapRule(newg,a,b,i**2) for i in range (2,22)]

#error between approx and real values
ErrorApprox = [realI - IvalsApprox[i] for i in range(0,len(IvalsApprox))]

#alpha term that indicates order of error term 
alphaApprox = [math.log(ErrorApprox[i]/ErrorApprox[i+1])/math.log(2) for i in range(0,len(IvalsApprox)-1)]


# problem (2) part (d)

def sqrtfunct(u):
    return math.sqrt(u)

actualsqrt = float(2)/3

#run of numbers approximating function
IvalsApproxSqrt = [CorrectedTrapRule(sqrtfunct,a,b,i**2) for i in range (2,22)]

#error between approx and real values
ErrorApproxSqrt = [actualsqrt - IvalsApproxSqrt[i] for i in range(0,len(IvalsApproxSqrt))]

#alpha term that indicates order of error term 
alphaApproxSqrt = [math.log(ErrorApproxSqrt[i]/ErrorApproxSqrt[i+1])/math.log(2) for i in range(0,len(IvalsApproxSqrt)-1)]

#function trying to approximate in problem (3)
def f(u):
   return math.cos(u)/math.sqrt(u)


# simpson's rule approximation problem (3) part (b)
def SimpsonsRule(f,a,b,n,usefirst):
    if usefirst:
        sum1 = f(a) + f(b)
    else:
        sum1=f(b)
    h = float(b-a)/n
    sum2=0
    sum3=0
    for i in range(2,n/2+1):
        sum2 += 2*f(a+(2*i-2)*h)
    for i in range(1,n/2+1):
        sum3 += 4*f(a+(2*i-1)*h)
    return (h/3)*(sum1+sum2+sum3)

#actual value of integral
#part (a) number (3) showing order of accuracy  

#real value trying to approximate in (a)
I = 1.809048475800308

#approximated values by trapezoidal rule
IvalsTrap = [TrapezoidalRule(f,a,b,2**i,False) for i in range (0,20)]

#error between approximated values and real
ErrorTrap = [I - IvalsTrap[i] for i in range(0,len(IvalsTrap))]

#constant created by the error term
TrapC = [ErrorTrap[i]*2**(float(i)/2) for i in range(0, len(IvalsTrap))]

alphaTrap = [math.log(ErrorTrap[i]/ErrorTrap[i+1])/math.log(2) for i in range(0,len(IvalsTrap)-1)]


# part (b) number (3) showing order of accuracy

#approximated values by simpsons rule
IvalsSim = [SimpsonsRule(f,a,b,2**i,False) for i in range (0,20)]

#error values of approximated from simpsons
ErrorSim = [I - IvalsSim[i] for i in range(0,len(IvalsSim))]

#constant error term
SimC = [ErrorSim[i]*2**(float(i)/2) for i in range(0, len(IvalsSim))]

#alpha that indicates order of accuracy of error for simpsons rule
alpha = [math.log(ErrorSim[i]/ErrorSim[i+1])/math.log(2) for i in range(0,len(IvalsSim)-1)]

#new function for problem (3) part (c)
def newf(u):
   return (math.cos(u)-1)/math.sqrt(u)

#new values for trapezoidal rule
IvalsTrapNew = [2 + TrapezoidalRule(newf,a,b,2**i,False) for i in range (0,20)]

#error for new trap approx
ErrorTrapNew = [abs(I - IvalsTrapNew[i]) for i in range(0,len(IvalsTrapNew))]

#error constant
TrapNewC = [ErrorTrapNew[i]*2**(float(i)/2) for i in range(0, len(IvalsTrapNew))]

#alpha indicating order of accuracy for new trap approx
alphaTrapNew = [math.log(ErrorTrapNew[i]/ErrorTrapNew[i+1])/math.log(2) for i in range(0,len(IvalsTrapNew)-1)]

#part (d) approximation following interpolation polynomial
def LinearInterp(f, a, b, n,usefirst):
    sum1 = 0
    h = float(b-a)/n
    if usefirst:
        sum1 = f(a) + f(b)
    else:
        sum1=f(b)
    for i in range(n):
        sum1 += float(2)*((a+h*(i+1))**(0.5)-(a+h*i)**(0.5))*(f(a+h*i) - (a+h*i)*(f(a+h*(i+1))-f(a+h*i))/(a+h*(i+1)-(a+h*i))) 
        + float(2)/float(3)*((a+h*i)**(1.5)-(a+h*i)**(1.5))*(f(a+h*(i+1))-f(a+h*i))/(a+h*(i+1) - (a+h*i))
    return(sum1)

IvalsInterp = [LinearInterp(f,a,b,2**i,False) for i in range (2,20)]

#error for new trap approx
ErrorInterp = [abs(I - IvalsInterp[i]) for i in range(0,len(IvalsInterp))]

#alpha indicating order of accuracy for new trap approx
alphaInterp = [math.log(ErrorInterp[i]/ErrorInterp[i+1])/math.log(2) for i in range(0,len(IvalsInterp)-1)]


