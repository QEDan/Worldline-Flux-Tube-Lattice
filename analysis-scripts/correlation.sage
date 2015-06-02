import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

f=open("/home/mazur/worldline/intoutratio.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

x = list()
y = list()
for row in varList:
  if 'intoutsc' in row[0]:
	x.append(eval(row[1]))
  else:
	y.append(eval(row[1]))

meanx = np.mean(x)
meany = np.mean(y)
sdx = np.std(x)
sdy = np.std(y)

r=0.0
for i in range(len(x)):
	r += (x[i]-meanx)*(y[i]-meany)/(sdx*sdy)

r /= (len(x)-1)
print(r)
