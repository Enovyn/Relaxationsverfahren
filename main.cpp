#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#define N 9     //Number of rows/cols
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
    return 0;
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
    int rank=0, size=0, world_size=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int matrix_size = N*N;
    int partmatrix_size = matrix_size/(world_size-1);
    int rest = matrix_size%(world_size-1);

    if(rank == 0){ //main thread
        double matrix[N][N] = { 0 }; //initalize matrix
        double res[N][N] = { 0 };

        initm(&matrix[0][0]);
        printm(&matrix[0][0]);
        cout << "\n\n";

        int command = 1;

        for(int i = 1; i < world_size-1; i++){  //distribute values
            MPI_Send(&matrix[(partmatrix_size*(i-1))/N][(partmatrix_size*(i-1))%N], partmatrix_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        MPI_Send(&matrix[(partmatrix_size*(world_size-1-1))/N][(partmatrix_size*(world_size-1-1))%N], partmatrix_size+rest, MPI_DOUBLE, world_size-1, 0, MPI_COMM_WORLD);

        for(int i = 1; i < world_size; i++){  //send start command
            MPI_Send(&command, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

}

//important
 /*   do {
        makeres(&matrix[0][0], &res[0][0]);     //Bestimme Residuum
        makenew(&matrix[0][0], &res[0][0]);     //Bestimme neue ErgebnisMatrix
        cout << "\n\n";
        printm(&matrix[0][0]);
    } while (getmaxres(&res[0][0]) > 0.0001);   //max Residuum prüfen
    printm(&matrix[0][0]);
    cout << "number of threads: " << world_size << endl;
*/

    else{
        int command = 0;
        else if((0 < rank)&&(rank<world_size-1)){   //receive matrix for standard threads
            double* matrix = (double*) malloc(partmatrix_size);
            MPI_Recv(matrix, partmatrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        else if(rank == world_size-1){  //receive matrix for last thread
            double* matrix = (double*) malloc(partmatrix_size+rest);
            MPI_Recv(matrix, partmatrix_size+rest, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Recv(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    //wait for command

        if(command == 1){//execute one iteration step
            if(partmatrix_size >= N){    //the thread calculates at least one complete line
                if(rank == 1){  //first thread which also includes first line

                }
            }else if(partmatrix_size < N){

            }
        }
    }

    MPI_Finalize();
    return 0;
}
