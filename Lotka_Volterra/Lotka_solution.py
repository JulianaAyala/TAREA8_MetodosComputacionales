
# In[50]:

import numpy
import scipy
import matplotlib
#from scipy.integrate import odeint


# In[51]:

data = loadtxt("lotka_volterra_obs.dat")

time = data.T[0]
n_prey = data.T[1]
n_predator = data.T[1]

print data.T[0]


# In[52]:

def likelihood(x_obs, x_model, y_obs, y_model):
    chi_x = 0.5*( ( sum(x_obs-x_model)**2 ) )
    chi_y = 0.5*( ( sum(y_obs-y_model)**2 ) )
    chi = chi_x+chi_y 
    return -chi


# In[53]:

def analitic_model(values, t=0):
    dxdt = values[0]*(alpha-(beta*values[1]) )  #NOOOOOOOOOOOOOOOOOOOOOOOOOOO
    dydt = -values[1]*(gamma-(delta*values[0]) )
    return dxdt, dydt


# In[54]:

def model(x, y, alpha, beta, gamma, delta):
    dxdt = x*(alpha-(beta*y) )
    dydt = -y*(gamma-(delta*x) )
    return dxdt, dydt


# In[55]:

# Solucion del sistema de ecuaciones diferenciales:
# Valores tentativos para los parametros libres:

alpha = 29
beta = 4
gamma = 50
delta = 2


# In[56]:

# Condiciones de equilibro de la solucion

equilib_1 = numpy.zeros(2)

# Teniendo las ecuaciones igualadas a cero, obtenemos para x y para y:

equilib_2 = [gamma/delta , alpha/beta]
equilib_2 = numpy.array(equilib_2)
# El orden de las variables alpha, beta, etc parece invertirse debido a que estas son las constantes de las cuales se
# acompañan las variables x & y respectivamente al igualarse a cero el sistema de ecuaciones

# Debe cumplirse que nuestro modelo sea igual a cero para estos dos casos:

all(model(equilib_1[0], equilib_1[1], alpha, beta, gamma, delta) == np.zeros(2))  and all(model(equilib_2[0], equilib_2[1], alpha, beta, gamma, delta)) == numpy.zeros(2)


# In[57]:

# Definimos una matriz jacobiana para poder solucionar la ecuacion por medio de scipy.integrate.odeint

def Jacob(values, alpha, beta, gamma, delta):
    J = [  [ alpha-(beta*values[1]) , -beta*values[0] ] , [ delta*values[1] , -gamma+(delta*values[0]) ]  ]
    J = numpy.array(J)
    return J


# In[58]:

# Soluciones del equilibrio:

zero_value1 = Jacob(equilib_1, alpha, beta, gamma, delta) # Solucion trivial: Extincion de ambas especies
zero_value2 = Jacob(equilib_2, alpha, beta, gamma, delta) # Solucion estable


# In[59]:

# Valores propios: matriz Jacobiana

eigen = linalg.eigvals(zero_value2) # Este valor nos da la frecuencia de la oscialcion de las problaciones, el perido es:

T = 2*numpy.pi/abs(eigen[1])


# In[60]:

# Integracion y solucion final del sistema de ecuaciones:

from scipy import integrate

t = linspace(time[0], time[-1], 1000)
x_init = [n_prey[0], n_predator[0]]
x_init = numpy.array(x_init)

values = integrate.odeint(analitic_model, x_init, t)


# In[61]:

prey, pred = values.T

fig1 = pylab.figure()
pylab.plot(t, prey, c="r", label="Presa" )
pylab.plot(t, pred, label="Cazador")
pylab.legend(loc="best")
pylab.xlabel("Tiempo")
pylab.ylabel("Poblacion")
pylab.title("Evolucion Temporal: Cazador-Presa")

fig1.savefig("analitic_sol.png")


# In[62]:

# Diagrama de fase

fig2 = pylab.figure()
pylab.plot(prey, pred, c="r", label="Presa-Cazador" )
pylab.xlabel("Presa")
pylab.ylabel("Predador")
pylab.title("Evolucion Temporal: Presa vs. Cazador")

fig2.savefig("fase.png")


# In[63]:

alpha_walk = []
beta_walk = []
gamma_walk = []
delta_walk = []

l_walk = []

alpha_walk = numpy.append(alpha_walk, numpy.random.random())
beta_walk = numpy.append(beta_walk, numpy.random.random())
gamma_walk = numpy.append(gamma_walk, numpy.random.random())
delta_walk = numpy.append(delta_walk, numpy.random.random())

# Resultados iniciales presa-cazador

init_values = model(n_prey, n_predator, alpha_walk[0], beta_walk[0], gamma_walk[0], delta_walk[0])


# In[64]:

n_iterations = 50000

for i in range(n_iterations):
    alpha_prime = numpy.random.normal(alpha_walk[i], 0.1)
    beta_prime = numpy.random.normal(beta_walk[i], 0.1)
    gamma_prime = numpy.random.normal(gamma_walk[i], 0.1)
    delta_prime = numpy.random.normal(delta_walk[i], 0.1)
    
    init_values = model(n_prey, n_predator, alpha_walk[i], beta_walk[i], gamma_walk[i], delta_walk[i])
    init_values_prime = model(n_prey, n_predator, alpha_prime, beta_prime, gamma_prime, delta_prime)
    
    l_prime = likelihood(n_prey, init_values_prime[0], n_predator, init_values_prime[1])
    l_init = likelihood(n_prey, init_values[0], n_predator, init_values[1])
    
    alpha_check = l_prime/l_init
    
    if(alpha_check<=1.0):
        alpha_walk = numpy.append(alpha_walk, alpha_prime)
        beta_walk = numpy.append(beta_walk, beta_prime)
        gamma_walk = numpy.append(gamma_walk, gamma_prime)
        delta_walk = numpy.append(delta_walk, delta_prime)
        
        l_walk = numpy.append(l_walk, l_prime)
    else: 
        beta_check = numpy.random.random()
        if(beta_check<=exp(-alpha_check)):
            alpha_walk = numpy.append(alpha_walk, alpha_prime)
            beta_walk = numpy.append(beta_walk, beta_prime)
            gamma_walk = numpy.append(gamma_walk, gamma_prime)
            delta_walk = numpy.append(delta_walk, delta_prime)
                
            l_walk = append(l_walk, l_prime)
        else:
            alpha_walk = numpy.append(alpha_walk, alpha_walk[i])
            beta_walk = numpy.append(beta_walk, beta_walk[i])
            gamma_walk = numpy.append(gamma_walk, gamma_walk[i])
            delta_walk = numpy.append(delta_walk, delta_walk[i])
            
            l_walk = append(l_walk, l_init)
                


# In[77]:

pyplot.scatter(alpha_walk, beta_walk)


# In[66]:

pyplot.scatter(gamma_walk, delta_walk)


# In[67]:

r = len(alpha_walk), len(l_walk)
print r
pyplot.scatter(alpha_walk[:-1], (l_walk))


# In[68]:

min_alpha = min(alpha_walk)
min_beta = min(beta_walk)
min_gamma = min(gamma_walk)
min_delta = min(delta_walk)

max_alpha = max(alpha_walk)
max_beta = max(beta_walk)
max_gamma = max(gamma_walk)
max_delta = max(delta_walk)

grid_alpha, grid_beta = mgrid[min_alpha:max_alpha:200j, min_beta:max_beta:200j]


# In[69]:

from scipy.interpolate import griddata

n_points = size(alpha_walk)
points = numpy.ones((n_points, 2))
print shape(points)
points[:,0] = alpha_walk
points[:,1] = beta_walk
grid_l = griddata(points, l_walk, (grid_alpha, grid_beta), method='cubic')
imshow(grid_l.T, extent=(min_alpha,max_alpha,min_beta,max_beta), aspect='auto',origin='lower')


# In[78]:

hist_alpha = pylab.figure()
count, bins, ignored =plt.hist(alpha_walk, 60, normed=True)
hist_alpha.savefig("hist_alpha.png")


# In[79]:

hist_beta = pylab.figure()
count, bins, ignored =plt.hist(beta_walk, 100, normed=True)
hist_beta.savefig("hist_beta.png")


# In[80]:

hist_gamma = pyplot.figure()
count, bins, ignored =plt.hist(gamma_walk, 50, normed=True)
hist_gamma.savefig("hist_gamma.png")


# In[81]:

hist_delta = pyplot.figure()
count, bins, ignored =plt.hist(delta_walk, 150, normed=True)
hist_delta.savefig("hist_delta.png")


# In[82]:

max_likelihood_id = argmax(l_walk)
best_alpha = alpha_walk[max_likelihood_id]
best_beta = beta_walk[max_likelihood_id]
best_gamma = gamma_walk[max_likelihood_id]
best_delta = delta_walk[max_likelihood_id]

print l_walk[max_likelihood_id]
print best_alpha
print best_beta
print best_gamma
print best_delta


# In[83]:

f_alpha = model(n_prey, n_predator, best_alpha, best_beta, best_gamma, best_delta)

plt.plot(time, f_alpha[0])
plt.plot(time, f_alpha[1])
scatter(time, n_prey, c="y")
scatter(time, n_predator, c="r", s=0.5)
plot(t, values)


# In[ ]:



