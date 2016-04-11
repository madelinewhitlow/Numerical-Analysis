
import numpy as np

#problem (1) 

actualVal = np.log(3)

def f(x):
    return float(1)/x

def Q(f,a,b,k):
    h = .1/2**k
    n = 10*(b-a)*2**k
    s = h*(f(a)/2.0 + f(b)/2.0)
    for i in range(1,n):
        s += h*f(a+i*h) 
    return s


def CompTrapezoidalRule(f,a,b,k):
    slist = []
    for i in range(0,5):
        Q1 = Q(f,a,b,i)
        Q2 = Q(f,a,b,i+1)
        #Q3 = Q(f,a,b,i+2)
        #check1 = Q2 - Q1
        #check2 = Q3 - Q2
        new = 4*Q2 - Q1
        #check = check1/check2
        #print check
        slist.append(new/float(3))
    return slist
    
slist = CompTrapezoidalRule(f,1,3,5)  
errorlist = abs(actualVal-slist)
 
alpha = [(np.log(errorlist[i]/errorlist[i+1]))/np.log(2) for i in range(0,len(errorlist)-1)]
