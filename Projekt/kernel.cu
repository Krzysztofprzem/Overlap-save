#include <device_functions.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <Windows.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <cstdlib>

using namespace std;

struct complex {
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						    POLA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	float Re;
	float Im;

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					    KONSTRUKTORY
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ complex(float Re = 0, float Im = 0) {
		this->Re = Re;
		this->Im = Im;
	}
	__host__ __device__ complex(complex &arg) {
		this->Re = arg.Re;
		this->Im = arg.Im;
	}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					     DESTRUKTOR
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ ~complex() {
	}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  GETTERY I SETTERY
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ void  set_Re(float Re) { this->Re = Re; }
	__host__ __device__ void  set_Im(float Im) { this->Im = Im; }
	__host__ __device__ float get_Re(		 ) { return Re;		}
	__host__ __device__ float get_Im(		 ) { return Im;		}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  PRZYDATNE METODY
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ string str_complex() {
		if (this->Im < 0) return to_string(this->Re) + " "  + to_string(this->Im) + "j";
		else			  return to_string(this->Re) + " +" + to_string(this->Im) + "j";
	}
	__host__ __device__ float magnitude() {
		return sqrt(Re * Re + Im * Im);
	}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						  OPERATORY
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ complex operator+(complex &arg) {
		complex result(this->get_Re() + arg.get_Re(),
					   this->get_Im() + arg.get_Im());
		return result;
	}
	__host__ __device__ complex operator-(complex &arg) {
		complex result(this->get_Re() - arg.get_Re(),
					   this->get_Im() - arg.get_Im());
		return result;
	}
	__host__ __device__ complex operator*(complex &arg) {
		complex result(this->get_Re() * arg.get_Re() - this->get_Im() * arg.get_Im(),
					   this->get_Re() * arg.get_Im() + this->get_Im() * arg.get_Re());
		return result;
	}
	__host__ __device__ complex operator*(float multiplier) {
		complex result(this->get_Re() * multiplier,
					   this->get_Im() * multiplier);
		return result;
	}
	__host__ __device__ complex operator*(int multiplier) {
		complex result(this->get_Re() * multiplier,
					   this->get_Im() * multiplier);
		return result;
	}
	__host__ __device__ complex operator/(int N) {
		complex result(this->get_Re() / N,
					   this->get_Im() / N);
		return result;
	}
	__host__ __device__ void operator=(complex &arg) {
		this->set_Re(arg.get_Re());
		this->set_Im(arg.get_Im());
	}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					METODY DLA OBIEKTÓW 
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	__host__ __device__ static complex twiddle_factor(int k, int N) { // <-zmienia obiekt w współczynnik obrotu
		float x = 0;
		float y = 0;
		// Szczególne przypadki dla których tangens się psuje :>
		if		(k	   == 0)	 { x =  1;		y =  0; }
		else if (k * 4 == N)	 { x =  0;		y = -1; }
		else if (k * 2 == N)	 { x = -1;		y =  0; }
		else if (k * 3 == N * 4) { x =  0;		y =  1; }
		// Reszta przypadków
		else {
			float PI = 3.141592653589793238;
			float fi = -(2 * PI * k / N);
			x = 1 / sqrt(tan(fi)*tan(fi) + 1);
			// Z wyprowadzonych wzorów wynika, że x jest zawsze dodatni, 
			// a to nie do końca prawda bo tangens nie uwzględnia ćwiartek
			if ((fi < -(PI / 2)) && (fi > -(3 * PI / 2)))
				x = -x;
			y = x * tan(fi);
		}
		complex result(x, y);
		return result;
	}
	__host__ __device__ static complex inv_twiddle_factor(int k, int N) { // <-zmienia obiekt w współczynnik obrotu
		float x = 0;
		float y = 0;
		// Szczególne przypadki dla których tangens się psuje :>
		if		(k	   == 0	   ) { x =  1;		y =  0; }
		else if (k * 4 == N	   ) { x =  0;		y =  1; }
		else if (k * 2 == N	   ) { x = -1;		y =  0; }
		else if (k * 3 == N * 4) { x =  0;		y = -1; }
		// Reszta przypadków
		else {
			float PI = 3.141592653589793238;
			float fi = 2 * PI * k / N;
			x = 1 / sqrt(tan(fi)*tan(fi) + 1);
			// Z wyprowadzonych wzorów wynika, że x jest zawsze dodatni, 
			// a to nie do końca prawda bo tangens nie uwzględnia ćwiartek
			if ((fi > (PI / 2)) && (fi < (3 * PI / 2)))
				x = -x;
			y = x * tan(fi);
		}
		complex result(x, y);
		return result;
	}
};


// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
//												 SPIS TREŚCI
//													 GPU
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
cudaError_t FIR(complex* h, float Fs, float Fc, int M, int N);
cudaError_t Blocks(complex* x, complex* b, int Nx, int M, int B, int N);
cudaError_t BlockFFT(complex* b, int B, int N, bool dir_inv);
cudaError_t Convolution(complex* b, complex* h, int B, int N);
cudaError_t Merge(complex* b, complex* y, int Ny, int M, int B, int N);


__global__ void FIRresponseGPU(complex *h, float Fs, float Fc, int M, int N) {		// M- długość odpowiedzi impulsowej Filtra (nieparzysta)
	int index = threadIdx.x;
	float PI = 3.141592653589793238;

	if (index == (M - 1) / 2)  h[index].set_Re(2 * Fc / Fs);
	else if (index < M)				h[index].set_Re((2 * Fc / Fs)*sin(2 * PI*Fc*(index - ((M - 1) / 2)) / Fs) / (2 * PI*Fc*(index - ((M - 1) / 2)) / Fs));
	else							h[index].set_Re(0);

	float OvH = 0.54 - 0.46*cos(2 * PI*index / (M - 1));	// Okno von Hanna
	h[index].set_Re(h[index].get_Re()*OvH);
}

__global__ void BlocksGPU(complex *x, complex *b, int Nx, int M) {
	int Bindex = blockIdx.x*blockDim.x + threadIdx.x;								// Jeden wątek - jedna próbka bloku
	int Xindex = blockIdx.x*blockDim.x + threadIdx.x - (blockIdx.x + 1)*(M - 1);	// Każdemu kolejnemu indeksowi bloku odpowiada odpowiednio przesunięty indeks z sygnału
	if (Xindex >= 0)						// W pierwszym bloku pierwsze M-1 próbek to zera
		if (Xindex < Nx)			// Pilnujemy, żeby nie wyjść poza zakres (w ostatnim bloku jest szansa, że trzeba dopisać zera (do 1024 punktowego FFT))
			b[Bindex] = x[Xindex];
}

__global__ void BlockFFTGPU(complex *b, int N, bool dir_inv) {
	// Każdy wątek na podstawie swojego indeksu określa 
	// ostateczną wartość próbki FFT o tym samym indeksie.
	int Bindex = threadIdx.x;						//indeks w pojedyńczym bloku
	int Oindex = Bindex + blockIdx.x*blockDim.x;	//indeks w całym "blokowisku"
	int Tindex = 0;									//nowy indeks
	extern __shared__ complex temp[];
	//complex* W = (complex*)&temp[N];
	// Trzeba przyjrzeć się diagramowi separującemu próbki parzyste i nieparzyste.
	// Na jego podstawie wyznaczone zostały wzory (męka niebywała) xD
	// SEPARACJA
	for (int i = 1; N / i >= 4; i *= 2) {
		if (Bindex % 2 == 0)
			Tindex = Bindex / 2 + (N / (2 * i) *(Bindex / (N / i)));
		else
			Tindex = (Bindex - 1) / 2 + (N / (2 * i) *((Bindex / (N / i)) + 1));
		__syncthreads();
		temp[Tindex] = b[Oindex];
		__syncthreads();
		b[Oindex] = temp[Bindex];
	}
	temp[Bindex] = b[Oindex];
	__syncthreads();

	// PROBLEM Z TWIDDLE FACTOR (YYY... Psuło się przy tworzeniu nowych obiektów. Może zapychały pamięć?)
	// Trzeba przyjrzeć się diagramowi motylkowemu próbki parzyste.
	// Na jego podstawie wyznaczone zostały wzory (niezbyt trudne).
	// DIAGRAM MOTYLKOWY
	for (int i = 2; i <= N; i *= 2) {
		__syncthreads();
		int butterflyindex = Bindex / i;
		if (Bindex - butterflyindex * i < i / 2) {
			int k = Bindex - butterflyindex * i;
			if (dir_inv) temp[Bindex] = b[Oindex] + (b[Oindex + i / 2] * complex::twiddle_factor(k, i));
			else		 temp[Bindex] = b[Oindex] + (b[Oindex + i / 2] * complex::inv_twiddle_factor(k, i));
		}
		else {
			int k = Bindex - butterflyindex * i - i / 2;
			if (dir_inv) temp[Bindex] = b[Oindex - i / 2] - (b[Oindex] * complex::twiddle_factor(k, i));
			else		 temp[Bindex] = b[Oindex - i / 2] - (b[Oindex] * complex::inv_twiddle_factor(k, i));
		}
		__syncthreads();
		b[Oindex] = temp[Bindex];
	}
	if (!dir_inv)
		b[Oindex] = b[Oindex] / N;
}

__global__ void ConvolutionGPU(complex*B, complex*H) {
	int Bindex = blockIdx.x*blockDim.x + threadIdx.x;
	int Findex = threadIdx.x;
	B[Bindex] = (B[Bindex] * H[Findex]);
}

__global__ void MergeGPU(complex *b, complex *y, int Ny, int M, int N) {
	int Yindex = blockIdx.x*blockDim.x + threadIdx.x;								// Jeden wątek - jedna próbka bloku
	int Bindex = blockIdx.x*blockDim.x + threadIdx.x + (blockIdx.x + 1)*(M - 1);	// Każdemu kolejnemu indeksowi wyjścia odpowiada odpowiednio przesunięty indeks z bloku
	y[Yindex].set_Re(b[Bindex].get_Re());
}
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
//												 SPIS TREŚCI
//													 GPU
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O


// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
//												 SPIS TREŚCI
//													 CPU
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
void FIRresponseCPU(complex* h, float Fs, float Fc, int M, int N) {
	float PI = 3.141592653589793238;
	float OvH = 0;
	for (int i = 0; i < (M - 1) / 2; i++) {
		OvH = 0.54 - 0.46*cos(2 * PI*i / (M - 1));	// Okno von Hanna
		h[i].set_Re(((2 * Fc / Fs)*sin(2 * PI*Fc*(i - ((M - 1) / 2)) / Fs) / (2 * PI*Fc*(i - ((M - 1) / 2)) / Fs))*OvH);
	}
	OvH = 0.54 - 0.46*cos((2 * PI*(M - 1) / 2) / (M - 1));
	h[(M - 1) / 2].set_Re((2 * Fc / Fs)*OvH);
	for (int i = ((M - 1) / 2 + 1); i < M; i++) {
		OvH = 0.54 - 0.46*cos(2 * PI*i / (M - 1));	// Okno von Hanna
		h[i].set_Re(((2 * Fc / Fs)*sin(2 * PI*Fc*(i - ((M - 1) / 2)) / Fs) / (2 * PI*Fc*(i - ((M - 1) / 2)) / Fs))*OvH);
	}
}

void BlocksCPU(complex* x, complex* b, int Nx, int M, int B, int N) {
	for (int i = 0; i < B; i++) {					// <-Bloki
		for (int j = 0; j < N; j++) {				// <-Próbki w bloku
			int Bindex = i * N + j;
			int Xindex = i * N + j - (i + 1)*(M - 1);
			if (Xindex >= 0)
				if (Xindex < Nx)
					b[Bindex] = x[Xindex];
		}
	}
}

void SeparateCPU(complex* X, complex* temp, int N) {
	for (int i = 0; i < N / 2; i++)
		temp[i] = X[i * 2 + 1];
	for (int i = 0; i < N / 2; i++)
		X[i] = X[i * 2];
	for (int i = 0; i < N / 2; i++)
		X[i + N / 2] = temp[i];
}

void FFTCPU(complex* X, int N, complex* temp, complex* e, complex* o) {
	if (N < 2) {

	}
	else {
		SeparateCPU(X, temp, N);
		FFTCPU(X		, N / 2, temp, e, o);
		FFTCPU(X + N / 2, N / 2, temp, e, o);
		for (int k = 0; k < N / 2; k++) {
			*e = X[k		];
			*o = X[k + N / 2];
			X[k		   ] = *e + (*o*complex::twiddle_factor(k, N));
			X[k + N / 2] = *e - (*o*complex::twiddle_factor(k, N));
		}
	}
}

void IFFTCPU(complex* X, int N, complex* temp, complex* e, complex* o) {
	if (N < 2) {

	}
	else {
		SeparateCPU(X, temp, N);
		IFFTCPU(X		 , N / 2, temp, e, o);
		IFFTCPU(X + N / 2, N / 2, temp, e, o);
		for (int k = 0; k < N / 2; k++) {
			*e = X[k		];
			*o = X[k + N / 2];
			X[k		   ] = *e + (*o*complex::inv_twiddle_factor(k, N));
			X[k + N / 2] = *e - (*o*complex::inv_twiddle_factor(k, N));
		}
	}
}

void BlockFFTCPU(complex* b, int B, int N, bool dir_inv) {
	complex* temp = new complex[N / 2];
	complex* e	  = new complex();
	complex* o	  = new complex();
	if (dir_inv) {
		for (int i = 0; i < B; i++)
			FFTCPU(b + i * N, N, temp, e, o);
	}
	else {
		for (int i = 0; i < B; i++)
			IFFTCPU(b + i * N, N, temp, e, o);
		for (int i = 0; i < B*N; i++)
			b[i] = b[i] / N;
	}
	delete[] temp;
	delete e, o;
}

void ConvolutionCPU(complex* b, complex* h, int B, int N) {
	for (int i = 0; i < B; i++)
		for (int j = 0; j < N; j++)
			b[i*N + j] = (b[i*N + j] * h[j]);

}

void MergeCPU(complex* b, complex*y, int M, int B, int N) {
	for (int i = 0; i < B; i++)
		for (int j = 0; j < N; j++) {
			int Yindex = i * N + j;
			int Bindex = i * N + j + (i + 1)*(M - 1);
			y[Yindex] = b[Bindex];
		}
}
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
//												 SPIS TREŚCI
//													 CPU
// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O


void GetSamples(complex** x, int channels) {
	fstream file;
	file.open("Samples.csv", ios::in);
	if (!file) {
		cerr << "NIEUDANE OTWARCIE PLIKU" << endl;
		exit(-1);
	}
	int n = 0;
	float sample;
	char c;
	while (file >> sample) {
		x[0][n].set_Re(sample);
		for (int i = 1; i < channels; i++) {
			if (i == channels - 1)
				file >> c;
			file >> sample;
			x[i][n].set_Re(sample);
		}
		n++;
		if (n % 100000 == 0)
			cout << "WCZYTANO " << n << " SAMPLI NA KANAŁ" << endl;
	}
	cout << "WCZYTANO " << n << " SAMPLI NA KANAŁ" << endl;
	file.close();
}

void GiveSamples(complex** y, int Ny, int channels) {
	fstream file;
	file.open("Output.csv", ios::out);
	if (!file) {
		cerr << "NIEUDANE OTWARCIE PLIKU" << endl;
		exit(-1);
	}
	int n = 0;
	while (n < Ny) {
		for (int j = 0; j < channels; j++) {
			file << y[j][n].get_Re();
			if (j != channels - 1) file << ",";
			else				   file << endl;
		}
		n++;
		if (n % 100000 == 0)
			cout << "ZAPISANO " << n << " SAMPLI NA KANAŁ" << endl;
	}
	cout << "ZAPISANO " << n << " SAMPLI NA KANAŁ" << endl;
	file.close();
}

int main() {
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						DANE DO EDYCJI
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	// DANE O SAMPLACH
	int   Nx	   = 5996784; //11685360;//5996784;				// Liczba próbek na kanał
	int	  channels = 2;						// Liczba kanałów
	float Fs	   = 44100;					// Częstotliwość próbkowania
	// DANE O PRZETWARZANIU
	int   N   = 1024;						// Rozmiar bloku									<- Zalecane 1024 (wykorzystanie pełnych bloków)
	int   M   = 61;							// Długość odpowiedzi impulsowej filtra FIR			<- Arbitralnie przyjęta
	float Fc  = 5000;						// Częstotliwość graniczna filtra					<- Arbitralnie przyjęta
	int	  NB  = 30000;						// Ilość bloków przetwarzanych w jednej turze		<- Edytować tylko w ekstremalnych przypadkach
	bool  CPU = false;						// Uruchomić przetwarzanie na CPU dla porównania?




	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						DANE DO PORÓWNANIA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	time_t startGPU;
	time_t stopGPU;
	time_t startCPU;
	time_t stopCPU;


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//							 PROLOG 
	//					 "WPROWADZENIE PRÓBEK"
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WPROWADZENIE PRÓBEK" << endl;
	complex** x = new complex*[channels];
	for (int i = 0; i < channels; i++)
		x[i] = new complex[Nx];
	GetSamples(x, channels);


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//													TOM I
	//											PRZETWARZANIE NA GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << endl << "ROZPOCZĘCIE PRZETWARZANIA NA GPU" << endl;
	startGPU = clock();


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//							 AKT 1 
	//					       ROZDZIAŁ 1
	//		 "TWORZENIE ODPOWIEDZI IMPULSOWEJ FILTRA FIR"										GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "TWORZENIE ODPOWIEDZI IMPULSOWEJ FILTRA FIR" << endl;
	complex* h = new complex[N];
	cudaError_t cudaStatus = FIR(h, Fs, Fc, M, N);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FIRresponse failed!");
		return 1;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 1 
	//						   ROZDZIAŁ 2
	//	  "WYZNACZANIE FFT ODPOWIEDZI IMPULSOWEJ FILTRA FIR"									GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WYZNACZANIE FFT ODPOWIEDZI IMPULSOWEJ FILTRA FIR" << endl;
	cudaStatus = BlockFFT(h, 1, N, true);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BlockFFT failed!");
		cout << endl << cudaGetErrorString(cudaStatus) << endl;
		return 1;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 2 
	//						   ROZDZIAŁ 1
	//			   "WYDZIELANIE BLOKÓW Z SYGNAŁU X"												GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WYDZIELANIE BLOKÓW Z SYGNAŁU X" << endl;
	int B = Nx % (N - (M - 1)) == 0 ? Nx / (N - (M - 1)) : Nx / (N - (M - 1)) + 1;	//ilość bloków
	complex** b = new complex*[channels];
	for (int i = 0; i < channels; i++)
		b[i] = new complex[B * N];

	for (int i = 0; i < channels; i++) {
		cudaStatus = Blocks(x[i], b[i], Nx, M, B, N);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Blocks failed!");
			return 1;
		}
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 2 
	//						   ROZDZIAŁ 2
	//			"WYZNACZANIE FFT POSZCZEGÓLNYCH BLOKÓW"											GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WYZNACZANIE FFT POSZCZEGÓLNYCH BLOKÓW" << endl;
	int rounds = (B%NB == 0) ? B / NB : B / NB + 1;
	int NBrest = B - ((rounds - 1)*NB);

	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < rounds; j++) {
			if (j == rounds - 1)	cudaStatus = BlockFFT(b[i] + j * NB*N, NBrest, N, true);
			else					cudaStatus = BlockFFT(b[i] + j * NB*N, NB,	   N, true);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "BlockFFT failed!");
				cout << endl << cudaGetErrorString(cudaStatus) << endl;
				return 1;
			}
		}
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 2 
	//						   ROZDZIAŁ 3
	//							"SPLOT"															GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "SPLOT" << endl;
	for (int i = 0; i < channels; i++) {
		cudaStatus = Convolution(b[i], h, B, N);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Convolution failed!");
			return 1;
		}
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 3 
	//						   ROZDZIAŁ 1
	//			"WYZNACZANIE IFFT POSZCZEGÓLNYCH BLOKÓW"										GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WYZNACZANIE IFFT POSZCZEGÓLNYCH BLOKÓW" << endl;
	for (int i = 0; i < channels; i++) {
		for (int j = 0; j < rounds; j++) {
			if (j == rounds - 1)	cudaStatus = BlockFFT(b[i] + j * NB*N, NBrest, N, false);
			else					cudaStatus = BlockFFT(b[i] + j * NB*N, NB,	   N, false);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "BlockIFFT failed!");
				cout << endl << cudaGetErrorString(cudaStatus) << endl;
				return 1;
			}
		}
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						     AKT 3 
	//						   ROZDZIAŁ 2
	//			  "ŁĄCZENIE BLOKÓW W SYGNAŁ WYNIKOWY"											GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "ŁĄCZENIE BLOKÓW W SYGNAŁ WYNIKOWY" << endl;
	int Ny = B * N - B * (M - 1);

	complex **y = new complex*[channels];
	for (int i = 0; i < channels; i++)
		y[i] = new complex[Ny];
	for (int i = 0; i < channels; i++) {
		cudaStatus = Merge(b[i], y[i], Ny, M, B, N);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Merge failed!");
			return 1;
		}
	}


	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}


	stopGPU = clock();
	double timeGPU = (double)(stopGPU - startGPU) / CLOCKS_PER_SEC;
	cout << "POMYŚLNIE ZAKOŃCZONO PRZETWARZANIE NA GPU" << endl;
	cout << "CZAS PRZETWARZANIA NA GPU WYNOSI: " << timeGPU << "s" << endl << endl;
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//												KONIEC TOMU I
	//											PRZETWARZANIE NA GPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O



	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//													TOM II
	//											PRZETWARZANIE NA CPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	if (CPU) {
		for (int i = 0; i < channels; i++) {
			delete[]b[i];
			delete[]y[i];
		}
		delete[]h;
		delete[]b;
		delete[]y;


		cout << "ROZPOCZĘCIE PRZETWARZANIA NA CPU" << endl;
		startCPU = clock();


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//							 AKT 1 
		//					       ROZDZIAŁ 1
		//		 "TWORZENIE ODPOWIEDZI IMPULSOWEJ FILTRA FIR"									CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "TWORZENIE ODPOWIEDZI IMPULSOWEJ FILTRA FIR" << endl;
		h = new complex[N];
		FIRresponseCPU(h, Fs, Fc, M, N);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 1 
		//						   ROZDZIAŁ 2
		//	  "WYZNACZANIE FFT ODPOWIEDZI IMPULSOWEJ FILTRA FIR"								CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "WYZNACZANIE FFT ODPOWIEDZI IMPULSOWEJ FILTRA FIR" << endl;
		BlockFFTCPU(h, 1, N, true);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 2 
		//						   ROZDZIAŁ 1
		//			   "WYDZIELANIE BLOKÓW Z SYGNAŁU X"											CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "WYDZIELANIE BLOKÓW Z SYGNAŁU X" << endl;
		int B = Nx % (N - (M - 1)) == 0 ? Nx / (N - (M - 1)) : Nx / (N - (M - 1)) + 1;	//ilość bloków
		b = new complex*[channels];
		for (int i = 0; i < channels; i++)
			b[i] = new complex[B * N];
		for (int i = 0; i < channels; i++)
			BlocksCPU(x[i], b[i], Nx, M, B, N);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 2 
		//						   ROZDZIAŁ 2
		//			"WYZNACZANIE FFT POSZCZEGÓLNYCH BLOKÓW"										CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "WYZNACZANIE FFT POSZCZEGÓLNYCH BLOKÓW" << endl;
		for (int i = 0; i < channels; i++)
			BlockFFTCPU(b[i], B, N, true);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 2 
		//						   ROZDZIAŁ 3
		//						   "SPLOT"														CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "SPLOT" << endl;
		for (int i = 0; i < channels; i++)
			ConvolutionCPU(b[i], h, B, N);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 3 
		//						   ROZDZIAŁ 1
		//			"WYZNACZANIE IFFT POSZCZEGÓLNYCH BLOKÓW"									CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "WYZNACZANIE IFFT POSZCZEGÓLNYCH BLOKÓW" << endl;
		for (int i = 0; i < channels; i++)
			BlockFFTCPU(b[i], B, N, false);


		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		//						     AKT 3 
		//						   ROZDZIAŁ 2
		//			  "ŁĄCZENIE BLOKÓW W SYGNAŁ WYNIKOWY"										CPU
		// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
		cout << "ŁĄCZENIE BLOKÓW W SYGNAŁ WYNIKOWY" << endl;
		int Ny = B * N - B * (M - 1);

		y = new complex*[channels];
		for (int i = 0; i < channels; i++)
			y[i] = new complex[Ny];
		for (int i = 0; i < channels; i++)
			MergeCPU(b[i], y[i], M, B, N - (M - 1));



		stopCPU = clock();
		double timeCPU = (double)(stopCPU - startCPU) / CLOCKS_PER_SEC;
		cout << "POMYŚLNIE ZAKOŃCZONO PRZETWARZANIE NA CPU" << endl;
		cout << "CZAS PRZETWARZANIA NA CPU WYNOSI: " << timeCPU << "s" << endl;
		cout << "GPU JEST " << timeCPU / timeGPU << " RAZY SZYBSZY OD CPU!" << endl << endl;
	}

	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//											   KONIEC TOMU II
	//											PRZETWARZANIE NA CPU
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O



	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						    EPILOG 
	//					"WYPROWADZENIE PRÓBEK"
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	cout << "WYPROWADZENIE PRÓBEK" << endl;
	GiveSamples(y, Ny, channels);


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//						  POSŁOWIE
	//					"CZYSZCZENIE PAMIĘCI"
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	for (int i = 0; i < channels; i++) {
		delete[]x[i];
		delete[]b[i];
		delete[]y[i];
	}
	delete[]h;
	delete[]x;
	delete[]b;
	delete[]y;


	//-----<>-----<>-----<>-----<>-----
	//	   TESTOWANIE TWORZENIA FIR
	//-----<>-----<>-----<>-----<>-----
	/*
	// TWORZENIE ODPOWIEDZI IMPULSOWEJ FILTRA FIR
	N = 8;
	Fs = 44100;
	Fc = 8500;
	M = 7;
	complex *h = new complex[N];

	cudaError_t cudaStatus = FIR(h, Fs, Fc, M, N);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FIRresponse failed!");
		return 1;
	}

	cout << "FIR" << endl;
	for (int i = 0; i < N; i++) {
		cout << h[i].str_complex() << endl;
	}

	cudaStatus = BlocksFFT(h, N, 1, N, true);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BlockFFT failed!");
		return 1;
	}
	cout << "FIR FFT" << endl;
	for (int i = 0; i < N; i++)
		cout << h[i].magnitude() << endl;
		*/

	//-----<>-----<>-----<>-----<>-----
	//	   TESTOWANIE BLOKOWEGO FFT
	//-----<>-----<>-----<>-----<>-----
	/*
	// x FFT
	int N = 8;
	int B = 50000;
	int size = B*N;
	complex*test = new complex[size];
	for (int i = 0; i < B * N; i+=N) {
		test[i].set_Re(1);
		test[i+2].set_Re(-1);
	}

	cout << "FFT" << endl;
	///for (int i = 0; i < size; i++)
	///	cout << "n= " << i << " " << test[i].str_complex() << endl;

	cudaError_t cudaStatus = BlocksFFT(test, B, N, true);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BlockFFT failed!");
		return 1;
	}
	for (int i = 0; i < size; i++)
		cout << "n= " << i << " " << test[i].str_complex() << endl;
	//cout << "IFFT" << endl;
	cudaStatus = BlocksFFT(test, B, N, false);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BlockIFFT failed!");
		return 1;
	}
	///for (int i = 0; i < size; i++)
	///	cout << "n= " << i << " " << test[i].str_complex() << endl;
	*/

	//-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----
	//			   TESTOWANIE TWORZENIA BLOKÓW, SPLOTU I ŁĄCZENIA W SYGNAŁ WYNIKOWY
	//-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----<>-----
	/*
	// TWORZENIE BLOKÓW
	int Nx = 23;	//długość sygnału
	int Nb = 0;
	int B = 0;
	N = 7;			//długość bloku
	M = 4;			//długość filtra

	B = Nx % (N - (M - 1)) == 0 ? Nx / (N - (M - 1)) : Nx / (N - (M - 1)) + 1;		//ilość bloków
	cout << "B=" << B << endl;
	Nb = B * N;		//łączna długość wszystkich bloków
	cout << "Nb=" << Nb << endl;


	complex *x = new complex[Nx];
	cout << "Sygnał x" << endl;
	for (int i = 0; i < Nx; i++) {
		x[i].set_Re(i);
		cout << x[i].str_complex() << endl;
	}

	complex *b = new complex[Nb];
	cout << "Bloki" << endl;
	for (int i = 0; i < Nb; i++) {
		cout << b[i].str_complex() << endl;
	}
	cout << "ładujemy x do bloków" << endl;
	cudaError_t cudaStatus2 = Blocks(x, b, Nx, M, B, N);

	if (cudaStatus2 != cudaSuccess) {
		fprintf(stderr, "Blocks failed!");
		return 1;
	}
	cout << "Załadowane Bloki" << endl;
	for (int i = 0; i < Nb; i++)
		cout << b[i].str_complex() << endl;

	complex *h = new complex[N];
	for (int i = 0; i < M; i++)
		h[i].set_Re(i + 1);
	cout << "Filtr" << endl;
	for (int i = 0; i < N; i++)
		cout << h[i].str_complex() << endl;

	// SPLOT
	cudaError_t cudaStatus3 = Convolution(b, h, B, N);
	if (cudaStatus3 != cudaSuccess) {
		fprintf(stderr, "Blocks failed!");
		return 1;
	}
	cout << "SPLOT" << endl;
	for (int i = 0; i < Nb; i++)
		cout << b[i].str_complex() << endl;

	// SCALANIE
	int Ny = Nb - B * (M - 1);
	complex *y = new complex[Ny];
	cout << "Sygnał wyjściowy" << endl;
	for (int i = 0; i < Ny; i++) {
		cout << y[i].str_complex() << endl;
	}
	cudaError_t cudaStatus4 = Merge(b, y, Ny, M, B, N);
	if (cudaStatus4 != cudaSuccess) {
		fprintf(stderr, "Merge failed!");
		return 1;
	}
	cout << "Scalanie" << endl;
	for (int i = 0; i < Ny; i++)
		cout << y[i].str_complex() << endl;
	*/


	return 0;
}


cudaError_t FIR(complex* h, float Fs, float Fc, int M, int N) {
	complex* dev_h = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	// Allocate GPU buffers.
	cudaStatus = cudaMalloc((void**)&dev_h, N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_h, h, N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  WYWOŁANIE KERNELA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	FIRresponseGPU << <1, N >> > (dev_h, Fs, Fc, M, N);


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "FIRresponse launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching FIRresponse!\n", cudaStatus);
		goto Error;
	}
	// Copy output from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(h, dev_h, N * sizeof(complex), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
Error:
	cudaFree(dev_h);
	return cudaStatus;
}

cudaError_t Blocks(complex* x, complex* b, int Nx, int M, int B, int N) {
	complex* dev_x = 0;
	complex* dev_b = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	// Allocate GPU buffers.
	cudaStatus = cudaMalloc((void**)&dev_x, Nx * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_b, B*N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_x, x, Nx * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_b, b, B*N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  WYWOŁANIE KERNELA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	BlocksGPU << <B, N >> > (dev_x, dev_b, Nx, M);


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "SplitIntoBlocks launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching SplitIntoBlocks!\n", cudaStatus);
		goto Error;
	}
	// Copy output from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(b, dev_b, B*N * sizeof(complex), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_x);
	cudaFree(dev_b);
	return cudaStatus;
}

cudaError_t BlockFFT(complex* b, int B, int N, bool dir_inv) {
	complex* dev_b = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	// Allocate GPU buffers.
	cudaStatus = cudaMalloc((void**)&dev_b, B*N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, B*N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  WYWOŁANIE KERNELA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	BlockFFTGPU << <B, N, N * sizeof(complex) >> > (dev_b, N, dir_inv);


	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "BlockFFT launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching BlockFFT!\n", cudaStatus);
		goto Error;
	}
	// Copy output from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(b, dev_b, B*N * sizeof(complex), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_b);
	return cudaStatus;
}

cudaError_t Convolution(complex* b, complex* h, int B, int N) {
	complex* dev_b = 0;
	complex* dev_h = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	// Allocate GPU buffers.
	cudaStatus = cudaMalloc((void**)&dev_b, B*N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_h, N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, B*N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_h, h, N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  WYWOŁANIE KERNELA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	ConvolutionGPU << <B, N >> > (dev_b, dev_h);


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Conv launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching ConvolutionGPU!\n", cudaStatus);
		goto Error;
	}
	// Copy output from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(b, dev_b, B*N * sizeof(complex), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_b);
	cudaFree(dev_h);
	return cudaStatus;
}

cudaError_t Merge(complex* b, complex* y, int Ny, int M, int B, int N) {
	complex* dev_b = 0;
	complex* dev_y = 0;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	// Allocate GPU buffers.
	cudaStatus = cudaMalloc((void**)&dev_b, B*N * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_y, Ny * sizeof(complex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_b, b, B*N * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_y, y, Ny * sizeof(complex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}


	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	//					  WYWOŁANIE KERNELA
	// O=====<+>=====<+>=====<+>=====<+>=====<+>=====<+>=====O
	MergeGPU << <B, N - (M - 1) >> > (dev_b, dev_y, Ny, M, N);


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Merg launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching MergeGPU!\n", cudaStatus);
		goto Error;
	}
	// Copy output from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(y, dev_y, Ny * sizeof(complex), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_b);
	cudaFree(dev_y);
	return cudaStatus;
}