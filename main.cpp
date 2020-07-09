#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#define N 5     //Number of rows/cols
#define NUMT 4  //Number of Threads, same as in run.bat
#define PARTMATRIX_SIZE (N*N)/(NUMT-1)
#define REST (N*N)%(NUMT-1)

const double pi = acos(-1);

using namespace std;

void initm(double* matrix) {
    matrix[0] = 1;
    matrix[N*N-1] = 1;

    double h = 1/(double)(N-1);
    for (int i = 1; i < N; i++) {
        matrix[i] = matrix[i - 1] - h;
        matrix[i*N] = matrix[(i-1)*N] - h;
        matrix[(N - 1) * N + i] = matrix[(N - 1) * N + i - 1] + h;
        matrix[i * N + (N - 1)] = matrix[(i - 1) * N + (N - 1)] + h;
    }
}


void printm(double* matrix) {
    cout << setprecision(4);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << matrix[i*N+j] << "\t\t";
        }
        cout << endl;
    }
}


double getv(double* matrix,int i, int j) {

    return (4 * matrix[i * N + j] - matrix[(i - 1) * N + j] - matrix[(i + 1) * N + j] - matrix[i * N + j - 1] - matrix[i * N + j + 1]);
}

double getf(double* matrix, int i, int j) {
    //return 2 * pi * pi * sin(pi*i) * sin(pi*j);   //Fall2
    return 0;   //Fall 1
}

void makeres(double* matrix, double* res) {
    double h = 1 / (double)(N - 1);
    for (int i = 1; i < N-1; i++) {
        for (int j = 1; j < N-1; j++) {
            res[i * N + j] = matrix[i * N + j] *h*h-getv(matrix, i, j);
        }
    }
}

void makenew(double* matrix, double* res) {
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            matrix[i * N + j] = matrix[i * N + j] + 0.25*res[i * N + j];
        }
    }
}

double getmaxres(double* res) {
    double erg = res[1 * N + 1];
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            if (res[i * N + j] > erg) erg = res[i * N + j];
        }
    }
    return erg;
}

int main(int argc, char **argv)
{
    int rank=0, size=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if(rank == 0){ //main thread
        double matrix[N][N] = { 0 }; //initalize matrix
        double res[N][N] = { 0 };

        int command = 1;

        cout << "Teilmatrix_Groeße1: " << PARTMATRIX_SIZE << endl;
        cout << "Teilmatrix_Rest1: " << REST << endl;

        initm(&matrix[0][0]);
        printm(&matrix[0][0]);
        cout << "\n\n";

        for(int i = 1; i < NUMT-1; i++){
            MPI_Send(&matrix[(PARTMATRIX_SIZE*(i-1))/N][(PARTMATRIX_SIZE*(i-1))%N], PARTMATRIX_SIZE, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

        }
            MPI_Send(&matrix[(PARTMATRIX_SIZE*(NUMT-1-1))/N][(PARTMATRIX_SIZE*(NUMT-1-1))%N], PARTMATRIX_SIZE+REST, MPI_DOUBLE, NUMT-1, 0, MPI_COMM_WORLD);

        for(int i = 1; i < NUMT; i++){
            MPI_Send(&command, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

    }else{
        int my_matrix_size = 0; //saves number of values for this specific thread (last one might have more than the others)
        double* matrix_ptr;     //pointer to matrix so I can allocate the memory in the if cases and still use it afterwards
        int command;
        double matrix[N][N];    //create auxiliary matrix for each thread
        int start_row = (PARTMATRIX_SIZE*(rank-1))/N;   //calculate starting point for values in matrix
        int start_column = (PARTMATRIX_SIZE*(rank-1))%N;

        for (int i = 0; i < N; i++) {   //fill matrix with -1
            for (int j = 0; j < N; j++) {
                matrix[i][j] = -1;
            }
        }


        if(rank<NUMT-1){
            my_matrix_size = PARTMATRIX_SIZE;
        }
        else{       //last thread may have more values
            my_matrix_size = PARTMATRIX_SIZE+REST;
        }


        MPI_Recv(&matrix[start_row][start_column], my_matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive values from master thread and save them at the same indices as the used to be



        MPI_Recv(&command, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(command){    //execute one iteration step
            cout << rank << endl;
        }


    }

    MPI_Finalize();
    return 0;
}
