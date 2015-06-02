import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;

def f(F, rho, T, l):
	return 2.0*F*l*l*T/(l*l + rho*rho)**2

def EAigrand(F, rho, T, l):
	return rho*np.exp(-m2*T)/T**3*(f(F, rho, T, l)/np.tanh(f(F, rho, T, l))-1-1/3*f(F, rho, T, l)**2)
	

def EAconst(F, l):
	return 1.0/(4.0*np.pi)*numerical_integral(lambda rho: numerical_integral(lambda T: EAigrand(F, rho, T, l),0,Infinity)[0], 0,Infinity)[0]

print(EAconst(100.0, 10.0))


