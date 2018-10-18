#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <string.h>

#define N 1024 // number of row and col of matrix
#define BLOCK_SIZE 32 // block size for tiling method

int chunk = 0; // number of computations handled by each thread using OpenMP

//struct definition for  formatting & display resutls
typedef struct method {
        char name[50];
        double exec_time;
}method;


//Usuful functions
int** allocate_matrix(size_t size);
void free_matrix(int** matrix, size_t size);
int** transpose_matrix(int** A, size_t size);
int compare_matrix(int** A, int** B, size_t size);
int** zeros(size_t size);
int** ones(size_t size);
void display_results(method methods[10], double non_optim_exec_time);

// different matrix multiplication functions implementation
void multiply_matrix(int** A, int** B, int**C, size_t size);
void multiply_matrix_tmp(int** A, int** B, int** C, size_t size);
void multiply_matrix_omp(int** A, int** B, int** C, size_t size);
void multiply_matrix_tiling_2d(int** A, int** B, int** C,  size_t matrix_size, size_t block_size);
void multiply_matrix_tiling_3d(int** A, int** B, int** C,  size_t matrix_size, size_t block_size);
void multiply_matrix_transpose(int** A, int** B, int** C, size_t size);
void multiply_matrix_transpose_omp(int** A, int** B, int** C, size_t size);
void multiply_matrix_tiling_3d_with_transposition(int** A, int** B, int** C, size_t matrix_size, size_t block_size);
void multiply_matrix_tiling_2d_with_transposition(int** A, int** B, int** C, size_t matrix_size, size_t block_size );
void multiply_matrix_tiling_2d_with_transposition_omp(int** A, int** B, int** C, size_t matrix_size, size_t block_size );

int main(int argc, char** argv) {

    //Operand matrix
    int** A;
    int** B;
    // resutl matrix
    int** C;
    int** D;
    int** E;
    int** F;

    //time measurement vars
    struct timeval  tv1, tv2;
    double time_eclapse = 0.0;
    // Result of comparing mathods with a naive matrix multiplication
    int notEqual = 0;


    

    method methods[10];
    //execution time of the naive matrix multiplication
    int non_optim_exec_time = 0;
    
    //Init matrix
    A = ones(N);
    B = ones(N);
    // init the result matrix with zeros
    D = zeros(N);
    C = zeros(N);
    E = zeros(N); 
    F = zeros(N);

    
    printf("General infos:\nBlock size: %d\nMatrix size: %d x %d\n",BLOCK_SIZE,N,N);
    printf("_______________________________________________________\n");
    
    //Naive matrix multiplication
    printf("Multiply without optimization\n");

    gettimeofday(&tv1, NULL);
    
    multiply_matrix(A,B,C,N);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);
    non_optim_exec_time = time_eclapse;
    strcpy(methods[0].name, "Non-optimized");
    methods[0].exec_time = time_eclapse;
    
    //Second method
    printf("Multiply with temporary var optimization\n");

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tmp(A,B,D,N);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,D,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[1].name, "tmp optimization");
        methods[1].exec_time = time_eclapse;
    
    } 
    

    //Third
    printf("Multiply with OpenMP optimization\n");

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_omp(A,B,E,N);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,E,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[2].name, "OpenMP optimization");
        methods[2].exec_time = time_eclapse;
    
    }

    //Fourth
    printf("Multiply with matrix transpose optimization\n");

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_transpose(A,B,F,N);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,F,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[3].name, "transpose optimization");
        methods[3].exec_time = time_eclapse;
    
    }

    //fifth
    printf("Multiply with matrix transpose and omp optimization\n");

    free_matrix(D,N);
    D = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_transpose_omp(A,B,D,N);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,D,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[4].name, "transpose_OpenMP optimization");
        methods[4].exec_time = time_eclapse;
    
    }

    //Sixth
    printf("Multiply with 2D tiling optimization\n");

    free_matrix(E,N);
    E = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tiling_2d(A,B,E,N,BLOCK_SIZE);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,E,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[5].name, "2D_tiling optimization");
        methods[5].exec_time = time_eclapse;
    
    }


    //Seventh
    printf("Multiply with 3D tiling optimization\n");

    free_matrix(F,N);
    F = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tiling_3d(A,B,F,N,BLOCK_SIZE);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,F,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[6].name, "3D_tiling optimization");
        methods[6].exec_time = time_eclapse;
    
    }


    //8th
    printf("Multiply with 2D tiling and matrix transposition optimization\n");

    free_matrix(D,N);
    D = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tiling_2d_with_transposition(A,B,D,N,BLOCK_SIZE);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,D,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[7].name, "2D_tiling_transpose optimization");
        methods[7].exec_time = time_eclapse;
    
    }

    //9th
    printf("Multiply with 3D tiling and matrix transposition optimization\n");

    free_matrix(E,N);
    E = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tiling_2d_with_transposition(A,B,E,N,BLOCK_SIZE);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,E,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[8].name, "3D_tiling_transpose optimization");
        methods[8].exec_time = time_eclapse;
    
    }

    //10th
    printf("Multiply with 2D tiling and matrix transposition optimization with OpenMP\n");

    free_matrix(F,N);
    F = zeros(N);

    gettimeofday(&tv1, NULL);
    
    multiply_matrix_tiling_2d_with_transposition_omp(A,B,F,N,BLOCK_SIZE);
    
    gettimeofday(&tv2, NULL);

    time_eclapse = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf("Operation took: %f s\n",time_eclapse);

    notEqual = compare_matrix(C,F,N);
    if(notEqual) printf("The two multiply resutls are not equal\n");
    else{
        printf("The multiplication gave the same result for the two matrix\n");
        strcpy(methods[9].name, "2D_tiling_transpose_OpenMP optimization");
        methods[9].exec_time = time_eclapse;
    
    }

    //free all matrix
    free_matrix(A,N);
    free_matrix(B,N);
    free_matrix(C,N);
    free_matrix(D,N);
    free_matrix(E,N);
    free_matrix(F,N);

    //display execution time results and poucentage of imprevement
    display_results(methods, non_optim_exec_time);
    
    return 0;

}

int** allocate_matrix(size_t size)
{
    int **arr = malloc(size * sizeof(int *));
    for (int i=0; i<size; i++)
        arr[i] = malloc(size * sizeof(int));

    return arr;
}

void free_matrix(int** matrix, size_t size)
{
    for(int i = 0; i< size; i++)
        free(matrix[i]);
    free(matrix);
}

void multiply_matrix(int** A, int** B, int** C, size_t size)
{
    int tmp = 0;
    for(int i = 0; i<size; i++){
        for(int j = 0; j< size; j++){
            C[i][j] = 0;

            for(int k = 0; k < size; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
            
        }
    }
}

void multiply_matrix_tmp(int** A, int** B, int** C, size_t size)
{
    int tmp = 0;
    for(int i = 0; i<size; i++){
        for(int j = 0; j< size; j++){
            tmp = 0;

            for(int k = 0; k < size; k++){
                tmp += A[i][k]*B[k][j];
            }
            C[i][j] = tmp;
        }
    }
}

void multiply_matrix_omp(int** A, int** B, int** C, size_t size)
{
    int tmp = 0;
    int nthread = omp_get_num_threads();
    chunk = size/nthread;
    int tid ;
    #pragma omp parallel shared(A,B,C,nthread, chunk) private(tid,tmp)
        //tid =   omp_get_thread_num();
        #pragma omp for schedule(static,chunk) 
            for(int i = 0; i<size; i++){
            
                for(int j = 0; j< size; j++){
                    tmp = 0;
                    
                    for(int k = 0; k < size; k++){
                        tmp += A[i][k]*B[k][j];
                    }
                    C[i][j] = tmp;
                }
            }
        
        
    
    
}

int** transpose_matrix(int** A, size_t size){
    int** TrA = allocate_matrix(size);
    for(int i = 0 ; i < size ; i++){
        for(int j = 0; j < size; j++){
            TrA[i][j] = A[j][i];
        }
    }
    return TrA;
}

void multiply_matrix_transpose(int** A, int** B, int** C, size_t size)
{
    int** TrB = transpose_matrix(B,size);
    int tmp = 0;
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size ; j++)
        {
            tmp = 0;
            for(int k = 0; k < size; k++)
                tmp += A[i][k]*TrB[j][k];
            C[i][j] = tmp;
        }
    }
}

void multiply_matrix_transpose_omp(int** A, int** B, int** C, size_t size)
{
    int** TrB = transpose_matrix(B,size);
    
    int tmp = 0;
    int nthread = omp_get_num_threads();
    chunk = size/nthread;
    int tid ;
        #pragma omp parallel shared(A,B,C,nthread, chunk) private(tid,tmp)
        //tid =   omp_get_thread_num();
            #pragma omp for schedule(dynamic,chunk) 
            for(int i = 0; i < size; i++){
                for(int j = 0; j < size ; j++)
                {
                    tmp = 0;
                    for(int k = 0; k < size; k++)
                        tmp += A[i][k]*TrB[j][k];
                    C[i][j] = tmp;
                }
            }
}

int compare_matrix(int** A, int** B, size_t size)
{
    for(int i = 0; i<size; i++){
        for(int j = 0; j < size; j++){
            if(A[i][j] != B[i][j]){
                printf("A[%d][%d] = %d, B[%d][%d] = %d\n",i,j,A[i][j],i,j,B[i][j]);
                return 1;
            } 
        }
    }
    return 0;
}

void multiply_matrix_tiling_2d(int** A, int** B, int** C, size_t matrix_size, size_t block_size )
{
    /*
    for(int i =0 ; i< size; i++)
        for(int j = 0 ; j < size; j++)
            C[i][j] = 0; 
    */
    int tmp = 0;
    for(int i = 0; i < matrix_size/block_size; i++){
        for(int j = 0; j < matrix_size / block_size; j++){
            for(int k = block_size*i; k< block_size*(i+1); k++){
                for(int l = block_size*j; l< block_size*(j+1); l++){
                    tmp = 0;    
                    for(int n = 0; n < matrix_size; n++)
                            tmp += A[k][n]*B[n][l];
                    C[k][l] = tmp;
                }
            }
        }
    }
     
}


void multiply_matrix_tiling_3d(int** A, int** B, int** C, size_t matrix_size, size_t block_size)
{
    /*
    for(int i =0 ; i< size; i++)
        for(int j = 0 ; j < size; j++)
            C[i][j] = 0; 
    
    */
    for(int i = 0; i < matrix_size / block_size; i++){
        for(int j = 0; j < matrix_size / block_size; j++){
            for(int k = 0; k< matrix_size / block_size; k++){
                for(int l = block_size*i; l< block_size*(i+1); l++)
                    for(int m = block_size*j; m < block_size*(j+1); m++ )
                    {
                        for(int n = block_size*k; n< block_size*(k+1); n++)
                            C[l][m] += A[l][n]*B[n][m];
                    }    
                        
                
            }
        }
    }
}


void multiply_matrix_tiling_2d_with_transposition(int** A, int** B, int** C, size_t matrix_size, size_t block_size )
{
    int** TrB = transpose_matrix(B,N);
    
    int tmp = 0;
    for(int i = 0; i < matrix_size/block_size; i++){
        for(int j = 0; j < matrix_size / block_size; j++){
            for(int k = block_size*i; k< block_size*(i+1); k++){
                for(int l = block_size*j; l< block_size*(j+1); l++){
                    tmp = 0;    
                    for(int n = 0; n < matrix_size; n++)
                            tmp += A[k][n]*TrB[l][n];
                    C[k][l] = tmp;
                }
            }
        }
    }
}




void multiply_matrix_tiling_3d_with_transposition(int** A, int** B, int** C, size_t matrix_size, size_t block_size)
{
    int** TrB = transpose_matrix(B,matrix_size);

    for(int i = 0; i < matrix_size / block_size; i++){
        for(int j = 0; j < matrix_size / block_size; j++){
            for(int k = 0; k< matrix_size / block_size; k++){
                for(int l = block_size*i; l< block_size*(i+1); l++)
                    for(int m = block_size*j; m < block_size*(j+1); m++ )
                    {
                        for(int n = block_size*k; n< block_size*(k+1); n++)
                            C[l][m] += A[l][n]*TrB[m][n];
                    }    
                        
                
            }
        }
    }
}


void multiply_matrix_tiling_2d_with_transposition_omp(int** A, int** B, int** C, size_t matrix_size, size_t block_size )
{
    int** TrB = transpose_matrix(B,N);
    int tmp = 0;
    int nthread = omp_get_num_threads();
    chunk = matrix_size/nthread;
    //int tid ;
    
        #pragma omp parallel shared(A,B,C,nthread, chunk) private(tmp)
            #pragma omp for schedule(dynamic,chunk) 
            for(int i = 0; i < matrix_size/block_size; i++){
                for(int j = 0; j < matrix_size / block_size; j++){
                    for(int k = block_size*i; k< block_size*(i+1); k++){
                        for(int l = block_size*j; l< block_size*(j+1); l++){
                            tmp = 0;    
                            for(int n = 0; n < matrix_size; n++)
                                    tmp += A[k][n]*TrB[l][n];
                            C[k][l] = tmp;
                        }
                    }
                }
            }
}


// Init a matrix with zeros
int** zeros(size_t size)
{
    int** A = allocate_matrix(size);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            A[i][j] = 0;
    return A;
}

// fill matrix with half top with ones and the rest with 2
// NB: A general init matrix function can be considered 
//to replace this and the one above 
int** ones(size_t size)
{
    int** A = allocate_matrix(size);
    int half_size = size/2;

    for(int i = 0; i < half_size; i++)
        for(int j = 0; j < size; j++)
            A[i][j] = 1;
    
    
    for(int i = half_size; i < size; i++)
    {
        /* code */
        
        for(int j = 0; j < size; j++)
        {
            A[i][j] = 2;
        }
        
    }
    
    return A;
}


void display_results(method methods[10], double non_optim_exec_time)
{
    double pourcentage = 0.0;
    char vars[3][15] = {"Method", "Taken time", "Improvement"};
    
    for(size_t i = 0; i < 71; i++)
    {
        printf("-");
    }
    printf("\n");
    printf("| %40s | %10s | %10s |\n",vars[0],vars[1],vars[2]);

    for(size_t i = 0; i < 71; i++)
    {
        printf("-");
    }
    printf("\n");
    for(size_t i = 0; i < 10; i++)
    {
        pourcentage = ((double)(non_optim_exec_time - methods[i].exec_time)/(double)non_optim_exec_time)*100;
        printf("| %40s | %10f | %11.2f |\n",methods[i].name,methods[i].exec_time,pourcentage);
        for(size_t j = 0; j < 71; j++)
        {
            printf("-");
        }   
        printf("\n");
    }

}