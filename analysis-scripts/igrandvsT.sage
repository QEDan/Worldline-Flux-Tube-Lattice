import numpy as np
import matplotlib.pyplot as plt



f=open("/home/mazur/worldline/EA.smooth3.F10l21.out", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

T1=list()
exact1=list()
igmean1=list()
igmeanerr1=list()

T2=list()
exact2=list()
igmean2=list()
igmeanerr2=list()

T3=list()
exact3=list()
igmean3=list()
igmeanerr3=list()
for row in varList:
	if row[0]=='WLvconstExpression:' and eval(row[1])==0.111111:
		T1.append(eval(row[2]))
		exact1.append(eval(row[3]))
		igmean1.append(eval(row[5]))
		igmeanerr1.append(eval(row[6]))
	if row[0]=='WLvconstExpression:' and eval(row[1])==0.444444:
		T2.append(eval(row[2]))
		exact2.append(eval(row[3]))
		igmean2.append(eval(row[5]))
		igmeanerr2.append(eval(row[6]))
	if row[0]=='WLvconstExpression:' and eval(row[1])==1.0:
		T3.append(eval(row[2]))
		exact3.append(eval(row[3]))
		igmean3.append(eval(row[5]))
		igmeanerr3.append(eval(row[6]))

data = zip(T1,exact1,igmean1, igmeanerr1)
data.sort()
[T1,exact1,igmean1,igmeanerr1] = zip(*data)

data = zip(T2,exact2,igmean2, igmeanerr2)
data.sort()
[T2,exact2,igmean2,igmeanerr2] = zip(*data)

data = zip(T3,exact3,igmean3, igmeanerr3)
data.sort()
[T3,exact3,igmean3,igmeanerr3] = zip(*data)


plt.clf()
plt.errorbar(T1,igmean1,igmeanerr1,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$\rho = 0.1111$')
plt.plot(T1,exact1,'k')
plt.errorbar(T2,igmean2,igmeanerr2,marker='x', color='k', ms=4, linewidth=0.0, elinewidth=1, label=r'$\rho = 0.4444$')
plt.plot(T2,exact2,'k')
plt.errorbar(T3,igmean3,igmeanerr3,marker='^', color='k', ms=4, linewidth=0.0, elinewidth=1, label=r'$\rho = 1.0000$')
plt.plot(T3,exact3,'k')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$T$")
plt.ylabel('Integrand')
plt.legend()
ax=plt.gca()
#ax.set_ylim(-1.0,0.02)
ax.set_xlim(0.0,2.0)
plt.savefig('igrandvsT4.F10l21.png')
plt.savefig('igrandvsT4.F10l21.eps')
