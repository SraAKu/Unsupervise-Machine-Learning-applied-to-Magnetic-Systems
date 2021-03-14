#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

double const J = 1.0;
int const Temps = 10;
double const Ti = 1.4;
double const Tf = 0.5;
int const L = 25;
unsigned int base = L*L;
unsigned int sites = L*L*L;       //Geometría
int samples = 80;
unsigned int MCRsteps = 2000;
unsigned int MCsteps  = 500; 

int ciclo(double, FILE*, FILE*, int[][6], int);
double getEnergy(float [], int[][6]);
double localEnergy(float [], int [][6], int);
double getMag(float []);
void sweep(float[], int[][6], double);

int main()
{   

	int neighbours[sites][6];

	//Definición de vecinos con condiciones de frontera periodicas

	unsigned int z = 0;

	for (unsigned int j = 0; j < sites; j++) {
		//printf("%.4f\t", spins[j]);

		//Vecino a derecha  OK
		neighbours[j][0] = j + 1;
		if ((j+1)%(L) == 0) {
			neighbours[j][0] = j + 1 - L;
		}

		//Vecino abajo  OK
		neighbours[j][1] = j - L;
		if (j< (L + (base * z))) {
			neighbours[j][1] = j - L + base;
		}

		//Vecino a izquierda  OK
		neighbours[j][2] = j - 1;
		if (j%(L) == 0) {
			neighbours[j][2] = j - 1 + L;
		}

		//Vecino arriba  OK
		neighbours[j][3] = j + L;
		if ((j + 1) > ( (base * (z + 1)) - L)) {
			neighbours[j][3] = j + L - base;
		}

		//Vecino atras  OK
		neighbours[j][4] = j - base;
		if (j < (base)) {
			neighbours[j][4] = j - base + sites;
		}

		//Vecino adelante OK
		neighbours[j][5] = j + base;
		if ((j + 1) > (sites - base)) {
			neighbours[j][5] = j + base - sites;
		}

		if ((j+1)%(base) == 0) {
			z += 1;
		}

		//printf("spin %d:%f\n",j, spins[j]);
		//printf("Indice de vecinos: %d %d %d %d %d %d\n",neighbours[j][0],neighbours[j][1],neighbours[j][2],neighbours[j][3],neighbours[j][4],neighbours[j][5]);
		//printf("ESPINES:           %f %f %f %f %f %f\n",spins[neighbours[j][0]],spins[neighbours[j][1]],spins[neighbours[j][2]],spins[neighbours[j][3]],spins[neighbours[j][4]],spins[neighbours[j][5]]);

	}

	FILE *ptrOT;
	char ObservablesT[30];
	sprintf(ObservablesT, "Data/Observables/OT_L%d.dat", L);
	ptrOT = fopen(ObservablesT,"w");

	double T = Ti;		
	double dt = (Tf - Ti) / (Temps - 1);

	for (unsigned int t = 0; t < Temps; t++) {
		printf("Temperatura %f\n", T);

		FILE *ptrS;
		char Samples[30];
		sprintf(Samples, "Data/Samples/S_L%dT_%.4f.dat", L, T);
		ptrS = fopen(Samples, "w");

		for (unsigned int sam = 0; sam < samples; sam++){
			printf("Sample %d\n", sam);
			ciclo(T,ptrOT, ptrS, neighbours, sam);
		}
		T = T + dt;
	}
}

int ciclo(double T, FILE *ptrOT, FILE *ptrS, int neighbours[][6], int sam)
{

	float spins[sites];

	/*Inicialización de spines*/

	for(unsigned int i = 0; i < sites; i++) { 
		srand(sam);
		spins[i] = 2.0 * ((float) rand() / (RAND_MAX)) - 1.0;
		if(spins[i]>=0.0) {
			spins[i] = 1.0;
		}
		else {
			spins[i] = -1.0;
		}
		//printf("%f\n",spins[i]);
	}

	//Pasos de relajación
	for(unsigned int q = 0; q < MCRsteps; q++) {
		sweep(spins, neighbours, T);
		srand(q);
	}

	//Archivo observables
	FILE *ptrO;
	char Observables[30];
	sprintf(Observables, "Data/Observables/O_L%dT_%.4f.dat", L, T);
	ptrO = fopen(Observables,"w");

	//Almacenamiento de observables

	double Eacum = 0.0;
	double Macum = 0.0;

	for(unsigned int p = 0; p < MCsteps; p++) {
		sweep(spins, neighbours, T);
		srand(p);
		double Energy = getEnergy(spins, neighbours);
		double Magnetization = getMag(spins);

		//printf("%f\n", Energy);

		Eacum += Energy;
		Macum += Magnetization;

		if(p >= (MCsteps-10)){
			fprintf(ptrO, "%d\t%.4f\t%.4f\n", p, Energy, Magnetization/sites);
		}
	}

	Eacum = Eacum/MCsteps/sites;
	Macum = Macum/MCsteps/sites;
	printf("%f\n", Eacum);

	fprintf(ptrOT, "%.4f\t%.4f\t%.4f\n", T, Eacum, Macum);

	for(unsigned int site = 0; site < sites; site++) {
		if(spins[site] < 0.0){
			fprintf(ptrS, "-1.0\t");
		}
		else{
			fprintf(ptrS, "1.0\t");
		}
	}

	fprintf(ptrS, "\n");



}

double getEnergy(float spins[], int neighbours[][6])
{
	double currEnergy = 0.0;

	for (unsigned int k = 0; k < sites; k++) {
		//printf("%f\n",spins[k]);
		//printf("*****VECINOS****\n");
		//printf("%f\t%f\t%f\n", spins[neighbours[k][0]], spins[neighbours[k][3]], spins[neighbours[k][5]]);
		double calc = spins[k] * (spins[neighbours[k][0]] + spins[neighbours[k][3]] + spins[neighbours[k][5]]) / 3;
		currEnergy += -J * calc;
		//printf("pos:%d  E=%f Ecum=%f\n", k,calc, currEnergy);
	}
	//printf("-----------------------\n");

	return currEnergy;
}

double localEnergy(float spins[], int neighbours[][6], int i)
{
	double local = 0.0;

	local = -J * spins[i] * (spins[neighbours[i][0]] + spins[neighbours[i][1]] + spins[neighbours[i][2]] + spins[neighbours[i][3]] + spins[neighbours[i][4]] + spins[neighbours[i][5]]);
	return local;
}

double getMag(float spins[])
{
	float M = 0;
	for (unsigned int l = 0; l < sites; l++) {
		M += spins[l];
	}

	return M;
}

void sweep(float spins[], int neighbours[][6], double T)
{
	for (unsigned int n = 0; n < sites; n++) {
		int pos = rand() % sites + 0;
		double E1 = localEnergy(spins, neighbours, pos);
		spins[pos] = -spins[pos];
		double E2 = localEnergy(spins, neighbours, pos);
		spins[pos] = -spins[pos];

		//printf("%.4f\t%.4f\t", E1, E2);

		double deltaE = E2 - E1;

		double r = ((double) rand() / (RAND_MAX));

		//printf("%.4f\n",deltaE);
		if ((deltaE <= 0.0) || (r < exp(-deltaE/T))) {
			spins[pos] = -spins[pos];
		}
	}
}
