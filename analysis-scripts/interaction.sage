import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

f=open("/home/mazur/worldline/interaction.dat", 'r')

varList = [line.split() for line in f.readlines()]
f.close()

Ea210000 = list()
Ea21p4 = list()
for row in varList:
  if 'a210000' in row[0]:
	Ea210000.append(-eval(row[1]))
  else:
	Ea21p4.append(-eval(row[1]))

diff = list()
ediff = list()
diff = [Ea21p4[i] - Ea210000[i] for i in range(len(Ea21p4))]
meand = mean(diff)
toMeV = 510.9989
ediff = sqrt(sageobj(r.var(Ea21p4)) + sageobj(r.var(Ea210000)) 
	- 2.0*sageobj(r.cov(Ea21p4,Ea210000)))
print(str(mean(diff)*toMeV) + ' +/- ' +str(ediff*toMeV))



