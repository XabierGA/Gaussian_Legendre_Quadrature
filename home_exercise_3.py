import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_style("whitegrid")

def pl_1(x):
    return x

def pl_0(x):
    return 1


#we can use recursive functions C-style with the recurrence formula for legendre polynomials
def generate_legendre(x , l , pl_0 , pl_1 ):
    if (l == 0):

        return pl_0(x)
    elif (l == 1):

        return pl_1(x)
    else:

        return ((2*l-1)*x*generate_legendre(x , l-1 , pl_0 , pl_1)-(l-1)*generate_legendre(x , l-2 , pl_0 , pl_1))/l

def bisection_legendre(n , f , a ,b , epsilon , maxsteps):
    if (f(a, n , pl_0 , pl_1)*f(b,n , pl_0 , pl_1)<=0):
        it = 1
        while ((np.abs(a-b)>=epsilon) and (it<=maxsteps)):
            c = (a+b)/2.0
            if (f(a, n , pl_0 , pl_1)*f(c, n, pl_0 , pl_1)>0):
                a = c
            elif (f(a, n , pl_0 , pl_1)*f(c, n, pl_0 , pl_1)<0):
                b = c
            elif(f(c,n, pl_0 , pl_1) == 0):
                return c
            it+=1

        return c
    else:
        #print("No roots found here")
        return None


def function(x):
    return (1-x**2)**(1/2)

def gaussian_quadrature(bisection_legendre , generate_legendre , func , polynomial_degree):
    integral = 0
    grid = np.arange(-1 , 1 , 0.001)
    r_found = []
    for i in range(len(grid)-1):
        root = bisection_legendre( polynomial_degree , generate_legendre , grid[i] , grid[i+1] , 0.0001 , 100000)
        #print(root)
        if (root !=None):
            r_found.append(root)
    N = len(r_found)
    for i in range(N):
        weight = (2 * (1- r_found[i] **2))/((N+1)**2*(generate_legendre(r_found[i] , N+1 , pl_0 , pl_1))**2)
        integral += weight*func(r_found[i])

    return integral



degrees = np.arange(2 , 15 , 1)
integrals = []
errors = []
for degree in degrees:
    result = gaussian_quadrature(bisection_legendre , generate_legendre , function , int(degree))
    errors.append(np.abs(np.pi/2 - result))
    integrals.append(result)
print(integrals[-1])
print(errors[-1])
fig = plt.figure(num=None, figsize=(10, 12), dpi=80, facecolor='w', edgecolor='k')
plt.plot(degrees , integrals , 'b-' , label = "Result of the integral" , linewidth= 7)
plt.xlabel("Degree of the polynomial")
plt.ylabel("Result of the integral")
plt.title("Numerical Integration using Gaussian-Legendre quadrature")
plt.legend(loc = 1)
fig.savefig("gaussian-legendre.pdf")
#plt.show(True)

fig_2 =plt.figure(num=None, figsize=(10, 12), dpi=80, facecolor='w', edgecolor='k')
plt.plot(degrees, errors , 'b-', label = "Absolute error" , linewidth = 7)
plt.xlabel("Degree of the polynomial")
plt.ylabel("Absolute error")
plt.title("Absolute error between analytical and numerical solution")
plt.legend(loc = 1)
fig_2.savefig("error-gaussian.pdf")
#plt.show(True)
