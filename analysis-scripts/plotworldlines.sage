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

def cmphase(N, x, y):
  phase = phases(N, x, y)
  cmphi = list()
  for phi in phase:
    cmphi.append(cm((phi - min(phase))/(max(phase)-min(phase))))
  return cmphi

MAP='jet' 


f=open("/home/mazur/worldline/onebig.worldline.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()


wlx = list()
wly = list()
wlx2 = list()
wly2 = list()
wlx3 = list()
wly3 = list()
wlx4 = list()
wly4 = list()
wlx5 = list()
wly5 = list()
for row in varList:
  wlx.append(eval(row[0]))
  wly.append(eval(row[1]))
  if len(wlx)%64 == 0:
    wlx2.append(eval(row[0]))
    wly2.append(eval(row[1]))
  if len(wlx)%1024 == 0:
    wlx3.append(eval(row[0]))
    wly3.append(eval(row[1]))
  if len(wlx)%8192 == 0:
    wlx4.append(eval(row[0]))
    wly4.append(eval(row[1]))
  if len(wlx)%8 == 0:
    wlx5.append(eval(row[0]))
    wly5.append(eval(row[1]))

print(len(wlx))
print(len(wlx2))
print(len(wlx3))
print(len(wlx4))

wlx.append(wlx[0])
wly.append(wly[0])
wlx2.append(wlx2[0])
wly2.append(wly2[0])
wlx3.append(wlx3[0])
wly3.append(wly3[0])
wlx4.append(wlx4[0])
wly4.append(wly4[0])
wlx5.append(wlx5[0])
wly5.append(wly5[0])

fig = plt.figure()
ax1 = fig.add_subplot(221) # regular resolution color map
ax2 = fig.add_subplot(222) # regular resolution alpha
ax3 = fig.add_subplot(223) # high resolution color map
ax4 = fig.add_subplot(224) # high resolution alpha

# Choose a color map, loop through the colors, and assign them to the color 
# cycle. You need NPOINTS-1 colors, because you'll plot that many lines 
# between pairs. In other words, your line is not cyclic, so there's 
# no line from end to beginning
cm = plt.get_cmap(MAP)
plt.subplot(221)
ax = plt.gca()
ax.set_color_cycle(cmphase(len(wlx4),wlx4,wly4))
for i in range(len(wlx4)-1):
    ax1.plot(wlx4[i:i+2],wly4[i:i+2])
plt.title(r'$N_{\rm ppl} = 16$')

for a in (ax.xaxis, ax.yaxis):
    for t in a.get_ticklines()+a.get_ticklabels():
        t.set_visible(False)
plt.subplot(222)
ax = plt.gca()
ax.set_color_cycle(cmphase(len(wlx3),wlx3,wly3))
for i in range(len(wlx3)-1):
    ax2.plot(wlx3[i:i+2],wly3[i:i+2])
plt.title(r'$N_{\rm ppl} = 128$')
#plt.plot(wlx3,wly3, 'k')
ax = plt.gca()
for a in (ax.xaxis, ax.yaxis):
    for t in a.get_ticklines()+a.get_ticklabels():
        t.set_visible(False)
plt.subplot(223)
ax = plt.gca()
ax.set_color_cycle(cmphase(len(wlx2),wlx2,wly2))
for i in range(len(wlx2)-1):
    ax3.plot(wlx2[i:i+2],wly2[i:i+2])
plt.title(r'$N_{\rm ppl} = 2048$')
ax = plt.gca()
for a in (ax.xaxis, ax.yaxis):
    for t in a.get_ticklines()+a.get_ticklabels():
        t.set_visible(False)
plt.subplot(224)
ax = plt.gca()
ax.set_color_cycle(cmphase(len(wlx5),wlx5,wly5))
for i in range(len(wlx5)-1):
    ax4.plot(wlx5[i:i+2],wly5[i:i+2])
plt.title(r'$N_{\rm ppl} = 16384$')
ax = plt.gca()
for a in (ax.xaxis, ax.yaxis):
    for t in a.get_ticklines()+a.get_ticklabels():
        t.set_visible(False)
plt.savefig('worldlineplot.png')
plt.savefig('worldlineplot.eps')

call('scp worldlineplot.eps mazur-desktop:~/Documents/Thesis/images/.', shell=True)
