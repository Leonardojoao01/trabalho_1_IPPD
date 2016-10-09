#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

#define threads 2
double minR= -1.5;
double maxR= 0.7;
double minI= -1.0;
double maxI= 1.0;

int findMandelbrot(double cr, double ci, int max_iterations);

double mapToReal(int x, int imageWidth, double minR, double maxR);

double mapToImaginary(int y, int imageHeight, double minI, double maxI);

void teste(int imageHeight,int imageWidth, int maxN, int** vec);

void teste_1(int imageHeight,int imageWidth, int maxN, int** vec);

void blocos(int itmax, int largImg, int altImg, int largIni, int largFim, int altIni, int altFim, int** vec);

void parallel_mandelbrot_blocos_task_omp(int itmax, int largImg, int altImg, int largIni, int largFim, int altIni, int altFim, int** vec);

void geraImagem(int larg, int alt, int** vec);


int main(int argc, char ** argv){

	clock_t Time[2];
  Time[0] = clock();

	// Get the required input values from file...
	ifstream fin("input.txt");
	int imageWidth, imageHeight, maxN;
	//double minR, maxR, minI, maxI;

	if (!fin){
		cout << "Could not open file!" << endl;
		return 1;
	}

	fin >> imageWidth >> imageHeight >> maxN;
	//fin >> minR >> maxR >> minI >> maxI;
	fin.close(); // Not necessary, good practice :D

	int **vec;

	//------------ALOCA A MATRIZ P/ SALVAR O CÁLCULO DA PROFUNDIDADE--------------
	// Primeiro momento é alocado as linhas e após as colunas
	vec = (int **)malloc((imageWidth) * sizeof(int **));
    if(vec == NULL){
        exit(1);
    }
  for(int i = 0; i<imageWidth; i++){
      vec[i] = (int *)malloc(imageHeight * sizeof(int));
      if(vec[i] == NULL) {
          exit(1);
      }
  }

	//============================================================================
	//----------------------------------------------------------------------------

	teste(imageHeight, imageWidth, maxN, vec);

	//teste_1(imageHeight, imageWidth, maxN, vec);

	//blocos(maxN, imageWidth, imageHeight, 0, imageWidth, 0, imageHeight, vec);

	//----------------------------------------------------------------------------
	//============================================================================

	geraImagem(imageWidth, imageHeight, vec);
	 //Desaloca a espaço de memória, utilizado para salva os dados calculados em relação a profundidade
	free(vec);

	Time[1] = clock();
  double tempo_decorrido = (Time[1] - Time[0]) * 1000.0 / CLOCKS_PER_SEC;

	cout << tempo_decorrido << endl;
	return 0;
}

//==============================================================================
//----------------------------------Funções-------------------------------------

double mapToReal(int x, int imageWidth, double minR, double maxR){
	double range = maxR - minR;
	return x * (range / imageWidth) + minR;
}
				//----------------------------------------------------------------------
double mapToImaginary(int y, int imageHeight, double minI, double maxI){
	double range = maxI - minI;
	return y * (range / imageHeight) + minI;
}
				//----------------------------------------------------------------------
int findMandelbrot(double cr, double ci, int max_iterations){
	int i = 0;
	double zr = 0.0, zi = 0.0;
	//Formula retirada da literatura para gerar a profundidade
	while (i < max_iterations && zr * zr + zi * zi < 4.0){
		double temp = zr * zr - zi * zi + cr;
		zi = 2.0 * zr * zi + ci;
		zr = temp;
		i++;
	}
	return i;
}
				//----------------------------------------------------------------------
// Paraleliza a LINHA
//void teste(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec){

void teste(int imageHeight,int imageWidth, int maxN, int** vec){
	int x,y,n;
	double cr=0, ci=0;

	#pragma omp parallel for num_threads(threads) shared(vec, y) private(x,cr,ci,n) firstprivate(maxR, minR, maxI, minI)
	for (y = 0; y < imageHeight; y++){
		for (x = 0; x < imageWidth; x++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos

			cr = mapToReal(x, imageWidth, minR, maxR);
			ci = mapToImaginary(y, imageHeight, minI, maxI);

			n = findMandelbrot(cr, ci, maxN);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			vec[y][x]=n;
		}
	}
}
				//----------------------------------------------------------------------
// Paraleliza a COLUNA
//void teste_1(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec){
void teste_1(int imageHeight,int imageWidth, int maxN, int** vec){

	int x,y,n;
	double cr=0, ci=0;

	#pragma omp parallel for num_threads(threads) shared(vec,x) private(y,cr,ci,n) firstprivate(maxR, minR, maxI, minI)
	for (x = 0; x < imageHeight; x++){

		for (y = 0; y < imageWidth; y++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos

			cr = mapToReal(x, imageWidth, minR, maxR);
			ci = mapToImaginary(y, imageHeight, minI, maxI);

			n = findMandelbrot(cr, ci, maxN);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			vec[y][x]=n;
		}
	}

}

void blocos(int itmax, int largImg, int altImg, int largIni, int largFim, int altIni, int altFim, int** vec){
	#pragma omp parallel num_threads(threads)
        {
            #pragma omp single
            {
                //Time[0] = clock();
                parallel_mandelbrot_blocos_task_omp(itmax,largImg,altImg,0,largImg,0,altImg, vec);
                //Time[1] = clock();
            }
        }
}

void parallel_mandelbrot_blocos_task_omp(int itmax, int largImg, int altImg, int largIni, int largFim, int altIni, int altFim,  int **vec){
    int diffLarg, diffAlt;
		//double cr=0, ci=0;

		diffLarg = largFim - largIni;
    diffAlt = altFim - altIni;

		// maxR

		if(diffAlt <= 64 || diffLarg <= 64){
        for(int x = largIni; x < largFim; x++) {
            for(int y = altIni; y < altFim; y++) {

							double cr = mapToReal(y, largImg, minR, maxR);
							double ci = mapToImaginary(x, altImg, minI, maxI);

							int n = findMandelbrot(cr, ci, itmax);

							vec[x][y] = n;
            }
        }
    }
    else{
            #pragma omp task untied shared(vec) firstprivate(itmax,largImg, altImg,largIni,largFim,altIni,altFim)
            {
                parallel_mandelbrot_blocos_task_omp(itmax, largImg, altImg, largIni, ceil((largFim+largIni)/2), altIni, ceil((altFim+altIni)/2), vec);
            }
            #pragma omp task untied shared(vec) firstprivate(itmax,largImg, altImg,largIni,largFim,altIni,altFim)
            {
                parallel_mandelbrot_blocos_task_omp(itmax, largImg, altImg, ceil((largFim+largIni)/2), largFim, altIni, ceil((altFim+altIni)/2), vec);
            }
            #pragma omp task untied shared(vec) firstprivate(itmax,largImg, altImg,largIni,largFim,altIni,altFim)
            {
                parallel_mandelbrot_blocos_task_omp(itmax, largImg, altImg, largIni, ceil((largFim+largIni)/2), ceil((altFim+altIni)/2), altFim, vec);
            }
            #pragma omp task untied shared(vec) firstprivate(itmax,largImg, altImg,largIni,largFim,altIni,altFim)
            {
                parallel_mandelbrot_blocos_task_omp(itmax, largImg, altImg, ceil((largFim+largIni)/2), largFim, ceil((altFim+altIni)/2), altFim, vec);
            }
            #pragma omp taskwait
    }

}


void geraImagem(int larg, int alt, int** vec){
  //ofstream fout(nome);
	ofstream fout("output_image.ppm");
	fout << "P3" << endl; // "Magic Number" - PPM file
	fout << larg << " " << alt << endl; // Dimensions
	fout << "255" << endl;

	for(int i = 0; i < larg; i++){
	    for(int j = 0; j < alt; j++){
	        int prof = vec[i][j];  //descobre a profundidade do ponto
	        int red,green,blue;
            red = (prof % 16) *50;
						//green = prof % 256;
						//blue = prof % 256;
            green = (prof % 256) * 16;
            blue = (prof % 32) * 8;
            fout << red << " " << green << " " << blue << " ";
	    }
	    fout << endl;
	}
    fout.close();
}
