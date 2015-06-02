import numpy as np
import matplotlib.pyplot as plt



f=open("/home/mazur/worldline/EAff.F1l2p09a2.out", 'r')

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
	if row[0]=='WLvconstExpression:' and eval(row[1])==0.00:
		T1.append(eval(row[2]))
		exact1.append(eval(row[3]))
		igmean1.append(eval(row[5]))
		igmeanerr1.append(eval(row[6]))
	if row[0]=='WLvconstExpression:' and eval(row[1])==0.00125:
		T2.append(eval(row[2]))
		exact2.append(eval(row[3]))
		igmean2.append(eval(row[5]))
		igmeanerr2.append(eval(row[6]))
	if row[0]=='WLvconstExpression:' and eval(row[1])==2.0:
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

#data = zip(T3,exact3,igmean3, igmeanerr3)
#data.sort()
#[T3,exact3,igmean3,igmeanerr3] = zip(*data)


plt.clf()
plt.errorbar(T1,igmean1,igmeanerr1,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
plt.plot(T1,exact1,'k')
plt.errorbar(T2,igmean2,igmeanerr2,marker='.', color='g', ms=8, linewidth=0.0, elinewidth=1)
plt.plot(T2,exact2,'g')
#plt.errorbar(T3,igmean3,igmeanerr3,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1)
#plt.plot(T3,exact3,'k')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="-#>",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$T$")
plt.ylabel('Integrand')
ax=plt.gca()
ax.set_ylim(-300.0,350.0)
ax.set_xlim(0.0,2.0)
plt.savefig('igrandvsTff.F1l2p64a2.png')
plt.savefig('igrandvsTff.F1l2p64a2.eps')
