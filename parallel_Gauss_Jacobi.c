/* Paulo Yamagishi DRE:121072893 */
/* Gustavo Mariz DRE:121073784 */
/* Disciplina: Programacao Concorrente */
/* Prof.: Silvana Rossetto */
/* Trabalho Final*/
/* Codigo: Implementação do algoritmo de Gauss-Jacobi usando programação concorrente. */

// Bibliotecas
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h> 
#include <math.h>

// Variáveis globais
double *new_values, *A, *b, *old_values, *aux_global, epsilon; // Vetor que armazena os valores novos
int N, nthreads, seed; // Variáveis para controle de dimensão
int it = 0; // Variável para contagem de iterações
int conv = 1; //Variável que define quando convergiu
pthread_mutex_t mutex; // Variável para controle de exclusão mútua
pthread_cond_t cond; // Variável para controle de sincronização
double elapsed_time_seq; // Variável para controle de tempo


// Função que gera a matriz A e o vetor b do sistema linear Ax = b (A é uma matriz quadrada diagonal dominante) (b é um vetor aleatório)
void gera_Ab(){
    double sum;
    srand(seed);

    for(int i = 0; i < N; i++){
        sum = 0;
        for(int j = 0; j < N; j++){
            if(i != j){
                A[i*N + j] = (double)rand()/(double)rand();
                sum += A[i*N + j];
            }
        }
        A[i*N + i] = sum * 5; // Garante a diagonal dominante
        b[i] = (double)rand()/(double)rand(); // Randomiza o vetor b    
    }
}

// Função que imprime a matriz A e o vetor b
void imprime_Ab(){
    int i, j;

    printf("Matriz A:\n");
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            printf("%.2f ", A[i*N + j]);
        }
        printf("\n");
    }
    printf("\nVetor b:\n");
    for(i = 0; i < N; i++){
        printf("%.2f ", b[i]);
    }
    printf("\n");
}

//Função que inicia os vetores de valores antigos e novos
void inicia_vetores(){
    int i;

    for(i = 0; i < N; i++){
        new_values[i] = 1;  
        old_values[i] = 0;  
        aux_global[i] = 0;
    }
}

// Função que multiplica a matriz A por um vetor (Ax = b)
void multiplica_matriz_vetor(double *A, double *x, int inicio, int jump){
    int i, j;
    double soma = 0;

    for(i = inicio; i < N; i += jump){
        soma = 0;
        for(j = 0; j < N; j++){
            soma += A[i*N + j] * x[j];
        }
        aux_global[i] = soma;
    }

}

// Função que calcula a distância entre dois vetores
double distancia(double *x, double *y){
    double soma = 0;
    int i;

    for(i = 0; i < N; i++){
        soma += (x[i] - y[i])*(x[i] - y[i]);
    }
    return (double)sqrt(soma);
}

// Função que contém o calculo do Gauss-Jacobi sequencial
void gauss_jacobi_calculo() {
    // Atualiza os valores antigos
    double soma = 0;
    // Calcula os novos valores
    for(int i = 0; i < N; i++){
        old_values[i] = new_values[i];
    }
    for(int i = 0; i < N; i++){
        soma = 0;
        for(int j = 0; j < N; j++){
            if(i != j){
                soma += A[i*N + j] * old_values[j];
            }
        }
        new_values[i] = (b[i] - soma) / A[i*N + i];
    }
}

// Função que contém o calculo do Gauss-Jacobi concorrente
void gauss_jacobi_calculo_concorrente(int inicio) {
    // Atualiza os valores antigos
    double soma = 0;
    // Calcula os novos valores
    for(int i = inicio; i < N; i += nthreads){
        old_values[i] = new_values[i];
    }
    for(int i = inicio; i < N; i += nthreads){
        soma = 0;
        for(int j = 0; j < N; j++){
            if(i != j){
                soma += A[i*N + j] * old_values[j];
            }
        }
        new_values[i] = (b[i] - soma) / A[i*N + i];
    }
}

// Gauss-Jacobi sequencial
void gaussJacobi_Sequencial(int version){
    int i, j;
    int iter = 0;

    // Inicializa os vetores de valores antigos e novos
    inicia_vetores();

    if (version == 1){
        // Enquanto não convergir
        while(distancia(new_values, old_values) > epsilon){
            iter++;
            gauss_jacobi_calculo();
        }
    }
    else {
        // Enquanto não convergir
        while(distancia(aux_global, b) > epsilon){
            iter++;
            gauss_jacobi_calculo();
            multiplica_matriz_vetor(A, new_values, 0, 1);
        }
    }

    // Imprime o número de iterações
    printf("\n---------------------------------------------------------\n");
    printf("Iterações: %d\n", iter);
}

// Função barreira para sincronização das threads e checagem de convergência
void barreira(int nthreads, int version){
    static int count = 0;
    // Soma as iterações de cada thread
    pthread_mutex_lock(&mutex);
    it += 1;
    pthread_mutex_unlock(&mutex);

    pthread_mutex_lock(&mutex);
    count++;
    if(count < nthreads){
        pthread_cond_wait(&cond, &mutex);
    }else{
        count = 0;
        if (version == 1) {
            if (distancia(new_values, old_values) < epsilon) {conv = 0;}
        }
        else {
            if (distancia(aux_global, b) < epsilon) {conv = 0;}
        }
        pthread_cond_broadcast(&cond);
    }
    pthread_mutex_unlock(&mutex);
}

// Função que imprime os resultados 
void imprime_resultados_sequenciais(int version){
    // Variáveis para controle de tempo
    struct timespec start, end;
    // vetor auxiliar para armazenar o resultado da multiplicação de A por new_values

    inicia_vetores();
    elapsed_time_seq = 0;

    // Executa o algoritmo sequencial com a versão passada no parâmetro
    clock_gettime(CLOCK_MONOTONIC, &start);
    gaussJacobi_Sequencial(version);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_seq = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Tempo sequencial versão %d: %.25f\n", version, elapsed_time_seq);
    multiplica_matriz_vetor(A, new_values, 0, 1);
    printf("Distancia: %.25f\n", distancia(aux_global, b));
    printf("---------------------------------------------------------\n\n");
    
}

// Função Gauss-Jacobi concorrente
void *gaussJacobi_Concorrente(void *arg){
    int inicio = *((int *) arg);
    int iter = 0;
    double soma = 0;

    // Enquanto não convergir
    while(conv){
        iter++;
        // Atualiza os valores antigos
        gauss_jacobi_calculo_concorrente(inicio);
        barreira(nthreads, 1);
    }
    pthread_exit(NULL);
}

// Função Gauss-Jacobi concorrente
void *gaussJacobi_Concorrente2(void *arg){
    int inicio = *((int *) arg);
    int iter = 0;
    double soma = 0;

    // Enquanto não convergir
    while(conv){
        iter++;
        // Atualiza os valores antigos
        gauss_jacobi_calculo_concorrente(inicio);
        multiplica_matriz_vetor(A, new_values, inicio, nthreads);
        barreira(nthreads, 2);
    }
    pthread_exit(NULL);
}

// Função que cria as threads e imprime os resultados
void cria_threads(void* (*f)(void*), int version) {
    // Variáveis para controle de tempo
    struct timespec start, end;
    // Inicializa as variáveis de controle de tempo
    int i, j, k;
    // Variáveis para controle de threads
    pthread_t *tid;
    tid = (pthread_t *) malloc(sizeof(pthread_t)*nthreads);

    inicia_vetores();
    it = 0;
    conv = 1;
    elapsed_time_seq = 0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    // Cria as threads
    for (i = 0; i < nthreads; i++) {
        int *arg = malloc(sizeof(int)); *arg = i;
        if (pthread_create(tid + i, NULL, (*f), arg)) {
            printf("--ERRO: pthread_create()\n"); exit(-1);
        }
    }
    // Espera as threads terminarem
    for (i = 0; i < nthreads; i++) {
        if (pthread_join(*(tid+i), NULL)) {
            printf("--ERRO: pthread_join()\n"); exit(-1);
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_seq = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("\n---------------------------------------------------------\n");
    printf("Iterações: %d\n", it/nthreads);
    printf("Tempo concorrente Versão %d: %.25f\n", version, elapsed_time_seq);
    multiplica_matriz_vetor(A, new_values, 0, 1);
    printf("Distancia: %.25f\n", distancia(aux_global, b));
    printf("---------------------------------------------------------\n\n");

    //Libera a memória alocada  
    free(tid);
}

// Função que cria uma cópia de um vetor
double *copia_vetor(double *vet){
    double *result = (double *) malloc(sizeof(double)*N);
    for (int i = 0; i < N; i++) {
        result[i] = vet[i];
    }
    return result;
}

// Função que calcula a aceleção do algoritmo concorrente
double aceleracao(double t1, double t2){
    return t1/t2;
}

// Função que calcula a eficiência do algoritmo concorrente
double eficiencia(double t1, double t2){
    return aceleracao(t1, t2)/nthreads;
}

// Função principal
int main(int argc, char *argv[]) {

    // Variáveis para controle de iteração
    int i, j;

    // Variáveis para controle de distância
    double *vet1, *vet2, *vet3, *vet4;

    // Variáveis para controle de métricas
    double t1, t2, t3, t4;

    // Inicializa as variáveis de controle de dimensão
    if (argc < 5){
        printf("Digite: %s <dimensao da matriz> <numero de threads> <epsilon> <seed para a criação de matrizes>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    N = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    epsilon = atof(argv[3]);
    seed = atoi(argv[4]);
    if (N < 1 || nthreads < 1 || epsilon < 0 || seed < 0){
        printf("N > 1; nthreads >= 1; epsilon > 0; seed > 0\n");
        exit(EXIT_FAILURE);
    }

    // Inicializa as variáveis de controle de matrizes e vetores
    A = (double *) malloc(sizeof(double)*N*N);
    b = (double *) malloc(sizeof(double)*N);
    new_values = (double *) malloc(sizeof(double)*N);
    old_values = (double *) malloc(sizeof(double)*N);
    aux_global = (double *) malloc(sizeof(double)*N);

    // Gera a matriz A e o vetor b do sistema linear Ax = b
    gera_Ab();
    //imprime_Ab();

    // Executa o algoritmo sequencial 

    imprime_resultados_sequenciais(1);
    vet1 = copia_vetor(new_values);
    t1 = elapsed_time_seq;

    imprime_resultados_sequenciais(2);
    vet2 = copia_vetor(new_values);
    t2 = elapsed_time_seq;

    // Executa o algoritmo concorrente

    cria_threads(gaussJacobi_Concorrente, 1);
    vet3 = copia_vetor(new_values);
    t3 = elapsed_time_seq;

    cria_threads(gaussJacobi_Concorrente2, 2);
    vet4 = copia_vetor(new_values);
    t4 = elapsed_time_seq;

    // Imprime os resultados
    printf("---------------------------------------------------------\n");
    printf("Distancia entre a solução sequencial e concorrente VERSÃO 1: %.25f\n", distancia(vet1, vet3));
    printf("Aceleração VERSÃO 1: %.25f\n", aceleracao(t1, t3));
    printf("Eficiência VERSÃO 1: %.25f\n", eficiencia(t1, t3));
    printf("---------------------------------------------------------\n");
    printf("Distancia entre a solução sequencial e concorrente VERSÃO 2: %.25f\n", distancia(vet2, vet4));
    printf("Aceleração VERSÃO 2: %.25f\n", aceleracao(t2, t4));
    printf("Eficiência VERSÃO 2: %.25f\n", eficiencia(t2, t4));
    printf("---------------------------------------------------------\n");
    //Libera a memória alocada
    free(A);
    free(b);
    free(new_values);

    return 0;
}