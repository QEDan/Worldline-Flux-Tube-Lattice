import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from sage.gsl.all import spline

m2=1.0;
e=0.30282212;
#m2=1.0;
#e=1.0;
a = np.sqrt(8.0)


call('grep EAvsLambda EAscff2*a2.out > EAscffvsLambda.dat', shell=True)

f=open("/home/mazur/worldline/EAscffvsLambda.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

l=list()
EAscal=list()
eEAscal=list()
ratio = list()
eratio = list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	l.append(eval(row[1])/a)
	EAscal.append(-eval(row[3]))
	eEAscal.append(eval(row[4]))
	ratio.append(-EAscal[-1]/eval(row[2]))
	eratio.append(eEAscal[-1]/eval(row[2]))

call('grep EAvsLambda EAscff.F1l2p*a21.out > EAscffvsLambda.a21.dat', shell=True)

f=open("/home/mazur/worldline/EAscffvsLambda.a21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

la21=list()
EAscala21=list()
eEAscala21=list()
ratioa21 = list()
eratioa21 = list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	la21.append(eval(row[1]))
	EAscala21.append(-eval(row[3]))
	eEAscala21.append(eval(row[4]))
	ratioa21.append(-EAscala21[-1]/eval(row[2]))
	eratioa21.append(eEAscala21[-1]/eval(row[2]))

call('grep EAvsLambda EAscff*a2p125.out > EAscffvsLambda.a2p125.dat', shell=True)

f=open("/home/mazur/worldline/EAscffvsLambda.a2p125.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

la2p125=list()
EAscala2p125=list()
eEAscala2p125=list()
ratioa2p125 = list()
eratioa2p125 = list()
for row in varList:
  if eval(row[4]) < abs(eval(row[3])):
	la2p125.append(eval(row[1])*np.sqrt(8.0))
	EAscala2p125.append(-eval(row[3]))
	eEAscala2p125.append(eval(row[4]))
	ratioa2p125.append(-EAscala2p125[-1]/eval(row[2]))
	eratioa2p125.append(eEAscala2p125[-1]/eval(row[2]))

f=open("/home/mazur/worldline/EAscffvsl.lcf.a28.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

la28lcf=list()
EAscala28lcf=list()
for row in varList:
	la28lcf.append(eval(row[0]))
	EAscala28lcf.append(-eval(row[1]))

f=open("/home/mazur/worldline/EAscffvsl.lcf.a21.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

la21lcf=list()
EAscala21lcf=list()
for row in varList:
	la21lcf.append(eval(row[0]))
	EAscala21lcf.append(-eval(row[1]))

f=open("/home/mazur/worldline/EAscffvsl.lcf.a2p125.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

la2p125lcf=list()
EAscala2p125lcf=list()
for row in varList:
	la2p125lcf.append(eval(row[0]))
	EAscala2p125lcf.append(-eval(row[1]))

plt.clf()
plt.errorbar(l,EAscal,eEAscal,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a = \sqrt{8}\lambda_c$')
plt.errorbar(la21,EAscala21,eEAscala21,marker='x', color='k', ms=8, linewidth=0.0, elinewidth=1, label = r'$a=\lambda_c$')
plt.errorbar(la2p125,EAscala2p125,eEAscala2p125,marker='^', color='k', ms=8, linewidth=0.0, elinewidth=1, label = r'$a=\frac{1}{\sqrt{8}}\lambda_c$')
plt.plot(la21lcf,EAscala21lcf,'k--')
plt.plot(la28lcf,EAscala28lcf,'k--')
plt.plot(la2p125lcf,EAscala2p125lcf, 'k--')
plt.xlabel(r"$\lambda/a$")
plt.ylabel(r"$-\Gamma^{(1)}_{\rm Scal}$")
ax=plt.gca()
plt.yscale('log')
#plt.xscale('log')
#ax.set_ylim(1e-10,1e20)
#ax.set_xlim(0.1,1.0)
plt.legend(loc=0)
plt.savefig('EAscffvsl.eps')
plt.savefig('EAscffvsl.png')

plt.clf()
plt.errorbar(l,EAscal,eEAscal,marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a = \sqrt{8}\lambda_c$')
plt.errorbar(la21,EAscala21,eEAscala21,marker='x', color='k', ms=8, linewidth=0.0, elinewidth=1, label = r'$a=\lambda_c$')
plt.errorbar(la2p125,EAscala2p125,eEAscala2p125,marker='^', color='k', ms=8, linewidth=0.0, elinewidth=1, label = r'$a=\frac{1}{\sqrt{8}}\lambda_c$')
plt.plot(la21lcf,EAscala21lcf,'k--')
plt.plot(la28lcf,EAscala28lcf,'k--')
plt.plot(la2p125lcf, EAscala2p125lcf, 'k--')
plt.xlabel(r"$\lambda/a$")
plt.ylabel(r"$E^{(1)}_{\rm Scal}$")
ax=plt.gca()
#plt.yscale('log')
#plt.xscale('log')
#ax.set_ylim(1e-10,1e20)
#ax.set_xlim(0.1,1.0)
plt.legend(loc=0)
plt.savefig('EAscffvsllinear.eps')
plt.savefig('EAscffvsllinear.png')

plt.clf()
plt.errorbar(l,ratio, eratio, marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a = \sqrt{8}\lambda_c$')
plt.errorbar(la21,ratioa21, eratioa21, marker='x', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a=\lambda_c$')
plt.errorbar(la2p125,ratioa2p125, eratioa2p125, marker='^', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a=\frac{1}{\sqrt{8}}\lambda_c$')
plt.xlabel(r"$\lambda/a$")
plt.ylabel(r"$\Gamma^{(1)}_{\rm Scal}/\Gamma_{\rm class}$")
plt.yscale('log')
plt.legend(loc=0)
plt.savefig('ActRatscff.eps')
plt.savefig('ActRatscff.png')

interpolationa21 = spline([(la21lcf[i],EAscala21lcf[i]) for i in range(len(la21lcf))])


resida21 = list()
eresida21 = list()
for i in range(len(la21)):
	resida21.append(EAscala21[i]/interpolationa21(la21[i]) - 1.0)
	eresida21.append(eEAscala21[i]/interpolationa21(la21[i]))

interpolationa28 = spline([(la28lcf[i],EAscala28lcf[i]) for i in range(len(la28lcf))])

resida28 = list()
eresida28 = list()
for i in range(len(l)):
	resida28.append(EAscal[i]/interpolationa28(l[i]) - 1.0)
	eresida28.append(eEAscal[i]/interpolationa28(l[i]))

interpolationa2p125 = spline([(la2p125lcf[i],EAscala2p125lcf[i]) for i in range(len(la2p125lcf))])

resida2p125 = list()
eresida2p125 = list()
for i in range(len(la2p125)):
	resida2p125.append(EAscal[i]/interpolationa2p125(l[i]) - 1.0)
	eresida2p125.append(eEAscal[i]/interpolationa2p125(l[i]))

plt.clf()
plt.errorbar(l, resida28, eresida28, marker='.', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a = \sqrt{8}\lambda_c$')
plt.errorbar(la21, resida21, eresida21, marker='x', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a=\lambda_c$')
plt.errorbar(la2p125, resida2p125, eresida2p125, marker='^', color='k', ms=8, linewidth=0.0, elinewidth=1, label=r'$a=\frac{1}{\sqrt{8}}\lambda_c$')
plt.xlabel(r"$\lambda/a$")
plt.ylabel(r"$\Gamma^{(1)}_{\rm Scal}/\Gamma^{(1)}_{\rm lcf}-1$")
#plt.yscale('log')
plt.axhline(y=0, c='k')
plt.legend(loc=0)
plt.savefig('EAscffresid.eps')
plt.savefig('EAscffresid.png')

call('scp EAscffvsl.eps ActRatscff.eps EAscffresid.eps mazur-desktop:~/Documents/Thesis/images/.',shell=True)

