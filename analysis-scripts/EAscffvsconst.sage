import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

e=0.3028221
q1=0.443991
q2=0.0742478
q3=0.0187671
F=1.0
a=1.0
td=1.0
lmin=0.1*td
An=(F/e)*(2.0/(q1*np.sqrt(2.0)))

def bump(x):
    if abs(x)<1:
        return exp(-1.0/(1.0-x**2))
    else:
        return 0.0
        
def B(rho,l,a):
    Bbg=6.0*F/(e*a**2)*(l-lmin)/(a-lmin)
    n=floor(rho/a+0.5)
    if n == 0:
       A0=4.0*F/(l**2*e*q2)*(1.0-0.75*(l-lmin)/(a-lmin))
       B=A0*bump(2.0*rho/l) + Bbg
    else:
       B=Bbg + 12.0*F/(q1*e*a*l)*((a-l)/(a-lmin))*bump(2.0*(rho-n*a)/l)
    return B


def DELO(rho, l, a):
        var('T')
        f=numerical_integral(exp(-T)/T**3*(e*B(rho,l, a)*T/tanh(e*B(rho,l, a)*T)
                -1.0+1.0/6.0*(e*B(rho,l, a)*T)**2),(0,100))
        return f[0]/(2.0*np.pi)


f=open("/home/mazur/worldline/EAscff.F1l2p09a21.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I=list()
xcm=list()
eI=list()
for row in varList:
	if row[0]=='Leffvsrho:':
		xcm.append(eval(row[1])/a)
		I.append(eval(row[2]))
		eI.append(eval(row[3]))

x=np.arange(0.0,0.3,5e-3)
#DE=[DELO(xi,0.3, a) for xi in x]


#DE = [-3209.71725336,-2936.73013681,-2490.31314204,-1893.55401841,-1201.84413114,-535.815812501,
#-98.6231019947,-0.550748224290,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
#-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
#-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
#-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228]

DE = [-3293.28555160,-3268.32687001,-3227.08914982,-3169.55070950,-3095.55024139,
-3005.60011341,-2899.58309952,-2778.08130014,-2641.45440100,-2490.31314204,
-2325.39478220,-2147.84432326,-1958.80989950,-1760.12013070,-1553.99342531,
-1343.23209559,-1131.33125294,-922.575232475,-722.092884924,-535.815812501,
-370.244365454,-231.872091647,-126.085122294,-55.3961420686,-17.2828780744,
-2.99638285276,-0.195028060671,-0.0101282156775,-0.00610544135068,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,
-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228,-0.00610357265228]

data = zip(xcm,I,eI)
data.sort()
[xcm,I,eI] = zip(*data)

interp = spline([(xcm[i], I[i]) for i in range(len(xcm))])
splineI = [interp(xi) for xi in x]

def fill_between_intersection(x, S, Z):
      #find the intersection point
      ind = 25
      # compute a new curve which we will fill below
      plt.fill_between(x[:ind], S[:ind], Z[:ind], facecolor='green', alpha=0.31)
      plt.fill_between(x[ind-1:], S[ind-1:], Z[ind-1:], facecolor='red', alpha=0.31)

plt.clf()
plt.cla()
#plt.errorbar(xcm,I,eI,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.plot(x,splineI,'k--', label='WLN result')
plt.plot(x,DE,'k', label='LCF approximation')
fill_between_intersection(x,splineI, DE)

#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}/a$")
plt.ylabel('Effective Action Density')
plt.legend(loc=0)
ax=plt.gca()
#ax.set_ylim(-0.1,0.1)
ax.set_xlim(0,0.20)
#ax.fill_between(x, splineI, DE, facecolor='blue', alpha=0.5)

plt.savefig('EAscffvsconst.png')
plt.savefig('EAscffvsconst.eps')

call('scp EAscffvsconst.eps mazur-desktop:~/Documents/Thesis/images/.', shell = True)
