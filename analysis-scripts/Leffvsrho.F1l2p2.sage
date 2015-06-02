import numpy as np
import matplotlib.pyplot as plt

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;


f=open("/home/mazur/worldline/EA.F1l2p2.2048.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I1=list()
xcm1=list()
eI1=list()
for row in varList:
	xcm1.append(eval(row[1]))
	I1.append(eval(row[2]))
	eI1.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1l2p214.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I2=list()
xcm2=list()
eI2=list()
for row in varList:
	xcm2.append(eval(row[1]))
	I2.append(eval(row[2]))
	eI2.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1l2p23.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I3=list()
xcm3=list()
eI3=list()
for row in varList:
	xcm3.append(eval(row[1]))
	I3.append(eval(row[2]))
	eI3.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1l2p267.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I4=list()
xcm4=list()
eI4=list()
for row in varList:
	xcm4.append(eval(row[1]))
	I4.append(eval(row[2]))
	eI4.append(eval(row[3]))

f=open("/home/mazur/worldline/EA.F1l2p3.2048.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

I5=list()
xcm5=list()
eI5=list()
for row in varList:
	xcm5.append(eval(row[1]))
	I5.append(eval(row[2]))
	eI5.append(eval(row[3]))


data=zip(xcm1,I1,eI1)
data.sort()
[xcm1,I1,eI1]=zip(*data)

data=zip(xcm2,I2,eI2)
data.sort()
[xcm2,I2,eI2]=zip(*data)

data=zip(xcm3,I3,eI3)
data.sort()
[xcm3,I3,eI3]=zip(*data)

data=zip(xcm4,I4,eI4)
data.sort()
[xcm4,I4,eI4]=zip(*data)

data=zip(xcm5,I5,eI5)
data.sort()
[xcm5,I5,eI5]=zip(*data)


plt.clf()
p1 =plt.errorbar(xcm1,I1,eI1,marker='^', color='k', label=r'$\lambda=0.447$',ms=6, 
	linewidth=1.0,elinewidth=1)
p2 =plt.errorbar(xcm2,I2,eI2,marker='o', color='k', label=r'$\lambda=0.463$', ms=6, 
	linewidth=1.0, elinewidth=1)
p3 =plt.errorbar(xcm3,I3,eI3,marker='x', color='k', label=r'$\lambda=0.480$',ms=6, 
	linewidth=1.0, elinewidth=1)
p4 =plt.errorbar(xcm4,I4,eI4,marker='.', color='k', label=r'$\lambda=0.517$', ms=8, 
	linewidth=1.0, linestyle='dashed', elinewidth=1)
p4 =plt.errorbar(xcm5,I5,eI5,marker='v', color='k', label=r'$\lambda=0.548$', ms=6, 
	linewidth=1.0, linestyle='dashed', elinewidth=1)



#plt.plot(x,DE,'k',x,DE2,'g')
#plt.annotate(r'LO Derivative Expansion',(x[250],DE[250]),xytext=(7.0,-0.01),arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"),bbox=dict(boxstyle="round", fc="1.0"))
plt.xlabel(r"$\rho_{cm}$")
plt.ylabel('Effective Action Density')
plt.legend()
ax=plt.gca()
#ax.set_ylim(-0.5,0.5)
ax.set_xlim(0.0,0.4)
plt.savefig('EA.F1l2p2.eps')
plt.savefig('EA.F1l2p2.png')
