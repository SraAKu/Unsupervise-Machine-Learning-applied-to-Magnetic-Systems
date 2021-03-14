import numpy as np
from matplotlib import pyplot

dirObservables = 'Data/Observables/'
dirSamples = 'Data/Samples/'
dirComponents = 'Data/PCA/'

L = 25;  #Dimensión de la muestra
Size = L*L*L  #Tamaño de la muestra
Ti = 5.0
Tf = 0.5
Temps = 46
MCRsteps = 2000
MCsteps  = 500
Nsamples = 80

lambAcum = np.zeros(Size)

Ts = np.linspace(Ti,Tf,Temps)

#Análisis de Samples

M = np.zeros((Nsamples * Temps, Size))

sMean = np.zeros(Size)

for s in range(Nsamples):

	for (iT,T) in enumerate(Ts):

		#Apertura de archivo de samples
		filenameS = '%sS_L%dT_%.4f.dat' %(dirSamples,L,T)
		file2 = open(filenameS, 'r')
		Samples = np.loadtxt(file2)
		Samples = np.array(Samples)

		M[s + (iT*Nsamples)][:] = Samples[s,:]
		sMean = sMean + Samples[s,:]

sMean = sMean / Temps /Nsamples

for i in range(Temps * Nsamples):
	M[i] = M[i] - sMean

(lamb, P) = np.linalg.eig(np.dot(M.T,M))

lamb = lamb.real
P = P.real

X = np.matmul(M,P)

Comps        = '%sComps.txt' %(dirComponents)
Comps = open(Comps, 'w')

for c in range(Size):
	Xc = X[:,c]
	for tt in range(Temps*Nsamples):
		Comps.write('%.4f\t' %Xc[tt])

	Comps.write('\n')

Lambs         = 'Data/PCA/Lambs.txt'
Lamb = open(Lambs, 'w')

for l in range(len(lamb)):
	Lamb.write('%.4f\t' %(lamb[l]))
