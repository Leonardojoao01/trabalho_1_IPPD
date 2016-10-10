#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

#define threads 2
double min_R= -1.5;
double max_R= 0.7;
double min_I= -1.0;
double max_I= 1.0;

int calculo_profundidade(double cr, double ci, int max_iterations);

double valor_real(int x, int largura_imagem, double min_R, double max_R);

double valor_imaginario(int y, int altura_imagem, double min_I, double max_I);

void paralelismo_linha(int altura_imagem,int largura_imagem, int iter_max, int** matriz);

void paralelismo_coluna(int altura_imagem,int largura_imagem, int iter_max, int** matriz);

void paralelismo_bloco(int iter_max, int largImg, int altImg, int** matriz);

void paralelismo_blocos_task_omp_recursivo(int iter_max, int largImg, int altImg, int larg_ini, int larg_fim, int alt_ini, int alt_fim, int** matriz);

void calcula_imagem(int larg, int alt, int** matriz);


int main(int argc, char ** argv){

	clock_t Time[2];
  Time[0] = clock();

	// Valores de entrada no arquivo
	ifstream fin("entrada.txt");
	int largura_imagem, altura_imagem, iter_max;

	if (!fin){
		return 1;
	}

	fin >> largura_imagem >> altura_imagem >> iter_max;
	fin.close(); // Fecha o arquivo, para não ficar ocupando memória

	int **matriz;		// utilizado para gerar uma matriz e passar como argumento aos outros modulos
	int i;
	//------------ALOCA A MATRIZ P/ SALVAR O CÁLCULO DA PROFUNDIDADE--------------
	// Primeiro momento é alocado as linhas e após as colunas
	matriz = (int **)malloc((largura_imagem) * sizeof(int **));
    if(matriz == NULL){
        exit(1);
    }
	// Quando é alocado as colunas
  for(i = 0; i<largura_imagem; i++){
      matriz[i] = (int *)malloc(altura_imagem * sizeof(int));
      if(matriz[i] == NULL) {
          exit(1);
      }
  }

	//============================================================================
	//----------------------------------------------------------------------------

	//paralelismo_linha(altura_imagem, largura_imagem, iter_max, matriz);

	//paralelismo_coluna(altura_imagem, largura_imagem, iter_max, matriz);

	paralelismo_bloco(iter_max, largura_imagem, altura_imagem, matriz);

	//----------------------------------------------------------------------------
	//============================================================================

	calcula_imagem(largura_imagem, altura_imagem, matriz);
	 //Desaloca a espaço de memória, utilizado para salva os dados calculados em relação a profundidade
	free(matriz);

	Time[1] = clock();
  double tempo_decorrido = (Time[1] - Time[0]) * 1000.0 / CLOCKS_PER_SEC;

	cout << tempo_decorrido << endl;
	return 0;
}

//==============================================================================
//----------------------------------Funções-------------------------------------

double valor_real(int x, int largura_imagem, double min_R, double max_R){
	double range = max_R - min_R;
	return x * (range / largura_imagem) + min_R;
}
				//----------------------------------------------------------------------
double valor_imaginario(int y, int altura_imagem, double min_I, double max_I){
	double range = max_I - min_I;
	return y * (range / altura_imagem) + min_I;
}
				//----------------------------------------------------------------------
int calculo_profundidade(double cr, double ci, int max_iterations){
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
//void paralelismo_linha(int altura_imagem,int largura_imagem, int iter_max, double min_R, double max_R, double min_I, double max_I, int** matriz){

void paralelismo_linha(int altura_imagem,int largura_imagem, int iter_max, int** matriz){
	int x,y,n;
	double cr=0, ci=0;

	#pragma omp parallel for num_threads(threads) shared(matriz, y) private(x,cr,ci,n) firstprivate(max_R, min_R, max_I, min_I)
	for (y = 0; y < altura_imagem; y++){
		for (x = 0; x < largura_imagem; x++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos

			cr = valor_real(x, largura_imagem, min_R, max_R);
			ci = valor_imaginario(y, altura_imagem, min_I, max_I);

			n = calculo_profundidade(cr, ci, iter_max);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			matriz[y][x]=n;
		}
	}
}
				//----------------------------------------------------------------------
// Paraleliza a COLUNA
//void paralelismo_coluna(int altura_imagem,int largura_imagem, int iter_max, double min_R, double max_R, double min_I, double max_I, int** matriz){
void paralelismo_coluna(int altura_imagem,int largura_imagem, int iter_max, int** matriz){

	int x,y,n;
	double cr=0, ci=0;

	#pragma omp parallel for num_threads(threads) shared(matriz,x) private(y,cr,ci,n) firstprivate(max_R, min_R, max_I, min_I)
	for (x = 0; x < altura_imagem; x++){

		for (y = 0; y < largura_imagem; y++){

			// Cálcula o valor COMPLEXO REAL e IMAGINARIO
			// // utilizando a linha e coluna para esses cálculos

			cr = valor_real(x, largura_imagem, min_R, max_R);
			ci = valor_imaginario(y, altura_imagem, min_I, max_I);

			n = calculo_profundidade(cr, ci, iter_max);	//Gera a profundidade

			//Salva o valor da profundidade em um vetor para posteriormente efetual o
			// cálculo da cor e salvar em um arquivo, gerando a imagem ".ppm"
			matriz[y][x]=n;
		}
	}

}

void paralelismo_bloco(int iter_max, int largImg, int altImg, int** matriz){

	#pragma omp parallel num_threads(threads)
        {
            #pragma omp single
            {
                paralelismo_blocos_task_omp_recursivo(iter_max,largImg,altImg,0,largImg,0,altImg, matriz);
            }
        }
}

void paralelismo_blocos_task_omp_recursivo(int iter_max, int largImg, int altImg, int larg_ini, int larg_fim, int alt_ini, int alt_fim,  int **matriz){
    int dif_larg, dif_alt;

		dif_larg = larg_fim - larg_ini;
    dif_alt = alt_fim - alt_ini;

		if(dif_alt <= 64 || dif_larg <= 64){
        for(int x = larg_ini; x < larg_fim; x++) {
            for(int y = alt_ini; y < alt_fim; y++) {

							double cr = valor_real(y, largImg, min_R, max_R);
							double ci = valor_imaginario(x, altImg, min_I, max_I);

							int n = calculo_profundidade(cr, ci, iter_max);

							matriz[x][y] = n;
            }
        }
    }
    else{
            #pragma omp task untied shared(matriz) firstprivate(iter_max,largImg, altImg,larg_ini,larg_fim,alt_ini,alt_fim)
            {
                paralelismo_blocos_task_omp_recursivo(iter_max, largImg, altImg, larg_ini, ceil((larg_fim+larg_ini)/2), alt_ini, ceil((alt_fim+alt_ini)/2), matriz);
            }
            #pragma omp task untied shared(matriz) firstprivate(iter_max,largImg, altImg,larg_ini,larg_fim,alt_ini,alt_fim)
            {
                paralelismo_blocos_task_omp_recursivo(iter_max, largImg, altImg, ceil((larg_fim+larg_ini)/2), larg_fim, alt_ini, ceil((alt_fim+alt_ini)/2), matriz);
            }
            #pragma omp task untied shared(matriz) firstprivate(iter_max,largImg, altImg,larg_ini,larg_fim,alt_ini,alt_fim)
            {
                paralelismo_blocos_task_omp_recursivo(iter_max, largImg, altImg, larg_ini, ceil((larg_fim+larg_ini)/2), ceil((alt_fim+alt_ini)/2), alt_fim, matriz);
            }
            #pragma omp task untied shared(matriz) firstprivate(iter_max,largImg, altImg,larg_ini,larg_fim,alt_ini,alt_fim)
            {
                paralelismo_blocos_task_omp_recursivo(iter_max, largImg, altImg, ceil((larg_fim+larg_ini)/2), larg_fim, ceil((alt_fim+alt_ini)/2), alt_fim, matriz);
            }
            #pragma omp taskwait
    }

}


void calcula_imagem(int larg, int alt, int** matriz){

	int i, j, red, green, blue, val_profundidade;

	ofstream criando_imagem("output_image.ppm");
	criando_imagem << "P3" << endl; 							// utilizado para gerar .ppm
	criando_imagem << larg << " " << alt << endl;
	criando_imagem << "255" << endl;							// Dimensão da imagem

	for(i = 0; i < larg; i++){
	    for(j = 0; j < alt; j++){
						//Valor da profundidade utilizado para calcula a cor no ponto
	        	val_profundidade = matriz[i][j];
            red = (val_profundidade % 16) *50;
						//green = prof % 256;
						//blue = prof % 256;
            green = (val_profundidade % 256) * 16;
            blue = (val_profundidade % 32) * 8;
            criando_imagem << red << " " << green << " " << blue << " "; // Escreve no arquivo
	    }
	    criando_imagem << endl;
	}
    criando_imagem.close();
}
