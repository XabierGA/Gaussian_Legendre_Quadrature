import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from home_exercise_3 import pl_1 , pl_0, function, gaussian_quadrature , generate_legendre , bisection_legendre
import math
sns.set()
sns.set_style("whitegrid")

N = 2000

def integrate(f,a,b, method = 'simpsons' , N = 4*150):
    h = (b-a)/N
#we may use a vectorized version
#    x = np.linspace(a,b , N+1)
#    y = f(x)
    integral = 0
    n = np.linspace(a,b,N)
    #print(n)
    if method == 'trapezoid':
        for i in n:
            value = (f(i-h) + 2*f(i) +f(i+h))/2.0
            if (not math.isnan(value)):
                integral+=value
        return h*integral/2.0
    elif method =='simpsons':
        for i in n:
            value = (f(i-h) + 4*f(i) +f(i+h))/2.0
            if (not math.isnan(value)):
                integral+=value
        return h*integral/3.0
    elif method =='bode':
        if N%4!=0:
            print("N must be multiple of 4 for this method")
        else:
            for i in n:
                value = (7*f(i) + 32*f(i+h)+12*f(i+2*h)+32*f(i+3*h)+7*f(i+4*h))/4.0
                if (not math.isnan(value)):
                    integral+=value
        return 2*h/45.0*integral


N = np.linspace(100 , 1000 ,10)

result_bode=[]
result_simpsons=[]
result_trapezoid = []
H = [2/n for n in N]
error_trazezoid = []
error_bode = []
error_simp = []
for i, numbers in enumerate(N):
    result_trapezoid.append(integrate(function , -0.999 , 0.999 , method ="trapezoid" ,N = numbers))
    result_simpsons.append(integrate(function , -1 , 1 , method = "simpsons" ,N = numbers))
    result_bode.append(integrate(function , -1 , 1 , method = "bode" , N = numbers))
    error_trazezoid.append(np.abs(result_trapezoid[i] - np.pi/2.0))
    error_simp.append(np.abs(result_simpsons[i] - np.pi/2.0))
    error_bode.append(np.abs(result_bode[i] - np.pi/2.0))

result_quadrature = [1.5709096321424787 for n in N]


fig_new = plt.figure(num=None, figsize=(10, 12), dpi=80, facecolor='w', edgecolor='k')
plt.plot(H , result_trapezoid , label = "Trapezoid Rule ")
plt.plot(H , result_simpsons , label = "Simpsons Rule")
plt.plot(H , result_bode , label = "Bode Rule")
plt.plot([0.0025] , [1.5709096321424787] , 'or', label = "Gaussian Quadrature, degree 15")
plt.xlabel("h")
plt.legend(loc = 3)
plt.ylabel("Integral Result")
plt.title("Comparison between different methods")
fig_new.savefig("comparison_methods.pdf")
plt.show(True)


fig_new_2 = plt.figure(num=None, figsize=(10, 12), dpi=80, facecolor='w', edgecolor='k')
plt.plot(H , error_trazezoid , label = "Trapezoid Rule ")
plt.plot(H , error_simp , label = "Simpsons Rule")
plt.plot(H , error_bode , label = "Bode Rule")
plt.plot([0.0025] , [0.00011330534758213773] , 'or', label = "Gaussian Quadrature, degree 15")
plt.xlabel("h")
plt.legend(loc = 2)
plt.ylabel("Integral Result")
plt.title("Comparison between errors for different methods")
fig_new_2.savefig("comparison_methods_errors.pdf")
plt.show(True)
