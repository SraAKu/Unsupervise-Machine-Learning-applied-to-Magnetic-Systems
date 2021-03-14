import numpy as np
from matplotlib import pyplot
import matplotlib.backends.backend_pdf
from mpl_toolkits.mplot3d import Axes3D
import math

L = 25
Size = L*L*L
Temps = 46
Nsamples = 80
Ti = 5.0
Tf = 0.5
Ts = np.linspace(Ti,Tf,Temps)

##Calculo RGB##

Dt = Ti - Tf
m = 1/Dt
Tmax = Ti
Tmin = Tf

####------####
dirComponents = 'Data/PCA/'


#Extracción parametros de orden: P0

filenameC = '%sComps.txt' %(dirComponents)
file = open(filenameC, 'r')
X = np.loadtxt(file)
X = np.array(X)

PProms = list()
PProm = 0

g = 0
for i in range(Temps*Nsamples):
	if( i and not(i%Nsamples)):
		PProm = PProm / Nsamples / math.sqrt(Size)
		print(PProm)
		PProms.append(PProm)
		PProm = 0
		g += 1

	PProm += abs(X[0][i])

PProms.append(PProm / Nsamples / math.sqrt(Size))

Parametros0 = 'Data/PCA/Parametros0.txt'
PParametros0 = open(Parametros0, 'a')

for pp in range(len(PProms)):
	PParametros0.write('%.4f\t' %(PProms[pp]))
PParametros0.write('\n')


#Extracción parametros de orden: P1

#filenameC = '%sComps.txt' %(dirComponents)
#file = open(filenameC, 'r')
#X = np.loadtxt(file)
#X = np.array(X)

PProms = list()
PProm = 0

g = 0
for i in range(Temps*Nsamples):
	if( i and not(i%Nsamples)):
		PProm = PProm / Nsamples / math.sqrt(Size)
		PProms.append(PProm)
		PProm = 0
		g += 1

	PProm += abs(X[1][i])

PProms.append(PProm / Nsamples / math.sqrt(Size))

Parametros1 = 'Data/PCA/Parametros1.txt'
PParametros1 = open(Parametros1, 'a')

for pp in range(len(PProms)):
	PParametros1.write('%.4f\t' %(PProms[pp]))
PParametros1.write('\n')


figs = list()

#Graficador componente vs componente

for k in range(5):
	for j in range(5):
		fig = pyplot.figure()

		filenameC = '%sComps.txt' %(dirComponents)
		file = open(filenameC, 'r')
		X = np.loadtxt(file)
		X = np.array(X)

		g = 0
		for i in range(Temps*Nsamples):
			if( i and not(i%Nsamples)):
				g += 1

			pyplot.scatter(X[k][i], X[j][i], c = [Ts[g]/Dt - (Tmin/Dt), (Ts[g]*(1/Ti))*0.001, -Ts[g]/Dt + (Tmax/Dt)], alpha=0.8)

		pyplot.title('Grafica primeras componentes: muestra $L^{3}$ = %d ' %(Size))	
		pyplot.xlabel('Componente 0 $P_{%d}$' %k)
		pyplot.ylabel('Componente %d $P_{%d}$' %(j,j))
		figs.append(fig)


#Graficador 3D primeras 3 componentes

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')

filenameC = '%sComps.txt' %(dirComponents)
file = open(filenameC, 'r')
X = np.loadtxt(file)
X = np.array(X)
		
g = 0	
for i in range(Temps*Nsamples):
	if( i and not(i%Nsamples)):
		g += 1

	ax.scatter(X[0][i],X[1][i],X[2][i], c = [Ts[g]/Dt - (Tmin/Dt), (Ts[g]*(1/Ti))*0.001, -Ts[g]/Dt + (Tmax/Dt)], alpha=0.8)

pyplot.title('Grafica componentes 0, 1 y 2: muestra $L^{3}$ = %d ' %(Size))	
figs.append(fig)

pdf = matplotlib.backends.backend_pdf.PdfPages("Graficas/Graficas de componentes L=%d.pdf" %(L))
for fig in figs :
    pdf.savefig( fig )
pdf.close()

filenameL = 'Data/PCA/Lambs.txt'
fileL = open(filenameL, 'r')
lambAcum = np.loadtxt(fileL)
lambAcum = np.array(lambAcum)

pyplot.figure()
pyplot.plot(range(20), lambAcum[0:20], "*")
pyplot.grid()
pyplot.savefig("lambda.pdf")
