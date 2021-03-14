import numpy as np
from matplotlib import pyplot
import matplotlib.backends.backend_pdf
from mpl_toolkits.mplot3d import Axes3D
import math

dirObservables = 'Data/Observables/'

Nsamples = 80
Temps = 46
Ti = 5.0
Tf = 0.5
Ts = np.linspace(Ti,Tf,Temps)

Ls = [10,15,20,25]


#Apertura de archivo de observables

pyplot.figure()

for L in Ls:
	filenameO = '%sOT_L%d.dat' %(dirObservables,L) 
	file1 = open(filenameO, 'r')
	Observables = np.loadtxt(file1)
	Observables = np.array(Observables)

	Enes = list()

	for TT in range(Temps):

		Enes.append(np.mean(Observables[0 + TT*Nsamples: (TT + 1)*Nsamples,1]))

	pyplot.plot(Ts, Enes, '-o', label='L=%i'%L)
	
pyplot.xlabel("Temperatura")
pyplot.ylabel("Energía")
pyplot.grid()
pyplot.legend(loc='best')
pyplot.savefig("EvsT.pdf")

pyplot.figure()

for L in Ls:
	filenameO = '%sOT_L%d.dat' %(dirObservables,L) 
	file1 = open(filenameO, 'r')
	Observables = np.loadtxt(file1)
	Observables = np.array(Observables)

	Mags = list()

	for TT in range(Temps):

		Mags.append(np.mean(Observables[0 + TT*Nsamples: (TT + 1)*Nsamples,2]**2))

	pyplot.plot(Ts, Mags,'-o', label='L=%i'%L)

pyplot.xlabel("Temperatura")
pyplot.ylabel("Magnetización")
pyplot.grid()
pyplot.legend(loc='best')
pyplot.savefig("MvsT.pdf")

# PARAMETOS DE ORDEN DESDE PRIMERAS COMPONENTES


filenameP0 = 'Data/PCA/Parametros0.txt'
file = open(filenameP0, 'r')
X = np.loadtxt(file)
X = np.array(X)
print(X)

pyplot.figure()

for L in range(len(Ls)):
	print(Ts)
	print(X[L])
	pyplot.plot(Ts,X[L],'-o', label='L=%i'%Ls[L])

pyplot.title('Grafica parametro de orden 0')	
pyplot.xlabel('$Temperatura$')
pyplot.ylabel('$<|P_{0}|>/L^{3/2}$')
pyplot.grid()
pyplot.legend(loc='best')
pyplot.savefig("Parametro P0")



filenameP1 = 'Data/PCA/Parametros1.txt'
file = open(filenameP1, 'r')
X = np.loadtxt(file)
X = np.array(X)

pyplot.figure()

for L in range(len(Ls)):
	pyplot.plot(Ts,X[L],'-o', label='L=%i'%Ls[L])

pyplot.title('Grafica parametro de orden 1')	
pyplot.xlabel('$Temperatura$')
pyplot.ylabel('$<|P_{1}|>/L^{3/2}$')
pyplot.grid()
pyplot.legend(loc='best')
pyplot.savefig("Parametro P1")


