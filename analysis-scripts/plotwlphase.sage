import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

m2=1.0;
e=0.30282212;

def phases(N, x, y):
  phi = list()
  phi.append(0.0)
  for i in range(N-1):
    phi.append((phi[-1] + 0.5*e*(x[i]*y[i+1] - y[i]*x[i+1])))
  return phi


f=open("/home/mazur/worldline/onebig.worldline.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()


wlx = list()
wly = list()

for row in varList:
  wlx.append(eval(row[0]))
  wly.append(eval(row[1]))



wlx.append(wlx[0])
wly.append(wly[0])

T = list()
for i in range(len(wlx)):
	T.append(1.0*i/len(wlx))
phase = phases(len(wlx),wlx,wly)


plt.figure(figsize=(5,8))
plt.cla()
plt.subplot(211)
plt.plot(wlx,wly,'k')
plt.annotate(r'$\tau = 0,1$',(wlx[0],wly[0]),xytext=(0.2,0.8),
	arrowprops=dict(arrowstyle="->", color='red'), color='black', size='x-large')
plt.annotate(r'$\tau = \frac{1}{4}$',(wlx[(len(wlx)-1)/4],wly[(len(wly)-1)/4]),xytext=(-0.75,-0.5),
	arrowprops=dict(arrowstyle="->", color='red'), color='black', size='x-large')
plt.annotate(r'$\tau = \frac{1}{2}$',(wlx[(len(wlx)-1)/2],wly[(len(wly)-1)/2]),xytext=(0.7,0.6),
	arrowprops=dict(arrowstyle="->", color='red'), color='black', size='x-large')
plt.annotate(r'$\tau = \frac{3}{4}$',(wlx[3*(len(wlx)-1)/4],wly[3*(len(wly)-1)/4]),xytext=(-0.85,0.9),
	arrowprops=dict(arrowstyle="->", color='red'), color='black', size='x-large')
plt.xlabel('x',size='large')
plt.ylabel('y',size='large')

plt.subplot(212)
plt.plot(T, phase, 'k')
plt.xlabel(r'$\tau$', size='x-large')
#plt.ylabel(r'$e \int_0^\tau d\tau A_\mu(x^\nu(\tau')) \dot{x}^\mu(\tau')$', size='x-large')
plt.ylabel('Phase', size='large')
ax = plt.gca()
ax.yaxis.set_ticks([1.*i*(np.pi/50)-np.pi/50 for i in range(7)])
ax.yaxis.set_ticklabels([r'$-\pi/50$',r'$0$',r'$\pi/50$',r'$2\pi/50$',r'$3\pi/50$',
	r'$4\pi/50$',r'$5\pi/50$'])
#plt.gcf().subplots_adjust(left=0.20)
plt.savefig('phiplot.eps')
plt.savefig('phiplot.png')

call('scp phiplot.eps mazur-desktop:~/Documents/Thesis/images/.', shell=True)

fig=plt.figure(frameon=False)

plt.cla()
plt.plot(wlx,wly, 'k')
ax = plt.gca()
for a in (ax.xaxis, ax.yaxis):
    for t in a.get_ticklines()+a.get_ticklabels():
        t.set_visible(False)
ax.set_frame_on(False)
plt.savefig('wlclean.eps')
plt.savefig('wlclean.png')

call('scp wlclean.eps mazur-desktop:~/Documents/Thesis/images/.', shell=True)
call('scp wlclean.png mazur-desktop:~/Documents/Thesis/images/.', shell=True)
