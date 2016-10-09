#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

#define threads 8

int findMandelbrot(double cr, double ci, int max_iterations);

double mapToReal(int x, int imageWidth, double minR, double maxR);

double mapToImaginary(int y, int imageHeight, double minI, double maxI);

void teste(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec);

void teste_1(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec);

void geraImagem(int larg, int alt, int** vec);


int main(int argc, char ** argv){

	clock_t Time[2];
  Time[0] = clock();

	// Get the required input values from file...
	ifstream fin("input.txt");
	int imageWidth, imageHeight, maxN;
	double minR, maxR, minI, maxI;

	if (!fin){
		cout << "Could not open file!" << endl;
		return 1;
	}

	fin >> imageWidth >> imageHeight >> maxN;
	fin >> minR >> maxR >> minI >> maxI;
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
	//--------------------------------------------------------------------------

	//teste(imageHeight, imageWidth, maxN, minR, maxR, minI, maxI, vec);

	teste_1(imageHeight, imageWidth, maxN, minR, maxR, minI, maxI, vec);

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
void teste(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec){

	int x,y;
	double cr=0, ci=0;

	#pragma omp parallel for num_threads(threads) shared(vec, y) private(x,cr,ci)
	for (y = 0; y < imageHeight; y++){
		for (x = 0; x < imageWidth; x++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos
			// double cr = mapToReal(x, imageWidth, minR, maxR);
			// double ci = mapToImaginary(y, imageHeight, minI, maxI);

			cr = mapToReal(x, imageWidth, minR, maxR);
			ci = mapToImaginary(y, imageHeight, minI, maxI);

			int n = findMandelbrot(cr, ci, maxN);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			vec[y][x]=n;
		}
	}
}
				//----------------------------------------------------------------------
// Paraleliza a COLUNA
void teste_1(int imageHeight,int imageWidth, int maxN, double minR, double maxR, double minI, double maxI, int** vec){

	int x,y;
	double cr=0, ci=0;

	for (y = 0; y < imageHeight; y++){
		#pragma omp parallel for num_threads(threads) shared(vec,y) private(x)
		for (x = 0; x < imageWidth; x++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos
			// double cr = mapToReal(x, imageWidth, minR, maxR);
			// double ci = mapToImaginary(y, imageHeight, minI, maxI);

			cr = mapToReal(x, imageWidth, minR, maxR);
			ci = mapToImaginary(y, imageHeight, minI, maxI);

			int n = findMandelbrot(cr, ci, maxN);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			vec[y][x]=n;
		}
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
						green = prof % 256;
						blue = prof % 256;
            //green = (prof % 256) * 16;
            //blue = (prof % 32) * 8;
            fout << red << " " << green << " " << blue << " ";
	    }
	    fout << endl;
	}
    fout.close();
}
