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

void initm(double* matrix){
    matrix[0] = 1;
    matrix[N*N-1] = 1;

    double h = 1/(double)(N-1);
    for (int i = 1; i < N; i++){
        matrix[i] = matrix[i - 1] - h;
        matrix[i*N] = matrix[(i-1)*N] - h;
        matrix[(N - 1) * N + i] = matrix[(N - 1) * N + i - 1] + h;
        matrix[i * N + (N - 1)] = matrix[(i - 1) * N + (N - 1)] + h;
    }
}


void printm(double* matrix){
    cout << setprecision(4);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << matrix[i*N+j] << "\t\t";
        }
        cout << endl;
    }
}


double getv(double* matrix, int i, int j){
    return (4 * matrix[i * N + j] - matrix[(i - 1) * N + j] - matrix[(i + 1) * N + j] - matrix[i * N + j - 1] - matrix[i * N + j + 1]);
}

double getf(double* matrix, int i, int j){
    //return 2 * pi * pi * sin(pi*i) * sin(pi*j);   //Fall2
    return 0;   //Fall 1
}

void makeres(double* matrix, double* res, int thread_rank, int matrix_size){ //geändert
    double h = 1 / (double)(N - 1);
    for (int i = 0; i < matrix_size; i++){
            int current_row = (PARTMATRIX_SIZE*(thread_rank-1)+i)/N;   //row and column where the value we are currently working on is
            int current_column = (PARTMATRIX_SIZE*(thread_rank-1)+i)%N;
            if((current_row != 0) && (current_row != N-1) && (current_column != 0) && (current_column != N-1)){ //don't calculate for first and last line and column
                res[current_row * N + current_column] = matrix[current_row * N + current_column] *h*h-getv(matrix, current_row, current_column);
            }
    }
}

void makenew(double* matrix, double* res, int thread_rank, int matrix_size){
    for (int i = 0; i < matrix_size; i++){
            int current_row = (PARTMATRIX_SIZE*(thread_rank-1)+i)/N;   //row and column where the value we are currently working on is
            int current_column = (PARTMATRIX_SIZE*(thread_rank-1)+i)%N;
            matrix[current_row * N + current_column] = matrix[current_row * N + current_column] + 0.25*res[current_row * N + current_column];
    }
}

double getmaxres(double* res){
    double erg = res[1 * N + 1];
    for (int i = 1; i < N - 1; i++){
        for (int j = 1; j < N - 1; j++){
            if (res[i * N + j] > erg)
                erg = res[i * N + j];
        }
    }
    return erg;
}

int main(int argc, char **argv){
    int rank=0, size=0;
    double res[N][N] = { 0 };
    double max_res = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if(rank == 0){  //main thread
        double matrix[N][N] = { 0 }; //initalize matrix

        int command = 1;

        cout << "Teilmatrix_Groesse: " << PARTMATRIX_SIZE << endl;
        cout << "Teilmatrix_Rest: " << REST << endl;

        initm(&matrix[0][0]);
        printm(&matrix[0][0]);
        cout << "\n\n";

        for(int i = 1; i < NUMT-1; i++){
            MPI_Send(&matrix[(PARTMATRIX_SIZE*(i-1))/N][(PARTMATRIX_SIZE*(i-1))%N], PARTMATRIX_SIZE, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);//spread values among threads
        }
        MPI_Send(&matrix[(PARTMATRIX_SIZE*(NUMT-1-1))/N][(PARTMATRIX_SIZE*(NUMT-1-1))%N], PARTMATRIX_SIZE+REST, MPI_DOUBLE, NUMT-1, 0, MPI_COMM_WORLD);

        do{//repeat iterations until max_res < 0.0001
            for(int i = 1; i < NUMT; i++){//send start command
                MPI_Send(&command, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            for(int i = 1; i < NUMT; i++){ //receive max_res of each thread
                max_res = 0;
                double temp = 0;
                MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(temp > max_res){//save max max_res value from slaves in max_res from master
                    max_res = temp;
                }
            }
            cout << "max " << max_res << endl;
        }while(max_res > 0.0001);
        command = 0;
        for(int i = 1; i < NUMT; i++){//send end command
                MPI_Send(&command, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for(int i = 1; i < NUMT-1; i++){
            MPI_Recv(&matrix[(PARTMATRIX_SIZE*(i-1))/N][(PARTMATRIX_SIZE*(i-1))%N], PARTMATRIX_SIZE, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//spread values among threads
        }
        MPI_Recv(&matrix[(PARTMATRIX_SIZE*(NUMT-2))/N][(PARTMATRIX_SIZE*(NUMT-2))%N], PARTMATRIX_SIZE+REST, MPI_DOUBLE, NUMT-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        cout << "ENDE" << endl;
        printm(&matrix[0][0]);

    }
    else{  //slaves
        int my_matrix_size = 0; //saves number of values for this specific thread (last one might have more than the others)
        int command;
        double matrix[N][N];    //create auxiliary matrix for each thread
        int start_row = (PARTMATRIX_SIZE*(rank-1))/N;   //calculate starting point for values in matrix
        int start_column = (PARTMATRIX_SIZE*(rank-1))%N;

        for (int i = 0; i < N; i++){     //fill matrix with -1
            for (int j = 0; j < N; j++){
                matrix[i][j] = -1;
            }
        }

        if(rank<NUMT-1){
            my_matrix_size = PARTMATRIX_SIZE;
        }
        else{        //last thread may have more values
            my_matrix_size = PARTMATRIX_SIZE+REST;
        }

        MPI_Recv(&matrix[start_row][start_column], my_matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive values from master thread and save them at the same indices as the used to be

        //////get values from other threads
        int end_row = ((PARTMATRIX_SIZE*(rank-1)+my_matrix_size)-1)/N;   //calculate starting point for values in matrix
        int end_column = ((PARTMATRIX_SIZE*(rank-1)+ my_matrix_size)-1)%N;

        MPI_Recv(&command, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        while(command){    //execute one iteration step (exchanges edge values (send and receive), calculate new matrix one time)
            /////////////send values
            for (int i = 0; i < my_matrix_size; i++){
                int current_row = (PARTMATRIX_SIZE*(rank-1)+i)/N;   //row and column where the value we are currently working on is
                int current_column = (PARTMATRIX_SIZE*(rank-1)+i)%N;

                int prev_rank = 0;
                int next_rank = 0;

                if(current_row < N-1){  //if we are not in the last line, we have to send our value to the field in the next line (i+1)
                    //send this value to field in next line (i+1)
                    int t = my_matrix_size - i;     //how many fields from this thread are left
                    if(t > N); //do nothing because the rank that needs this value is this rank
                    else{
                        t = N - t + 1; //how many fields are between the end of this thread and the desired field
                        int r= PARTMATRIX_SIZE;
                        next_rank = rank + ceil((double)t/(double)r);   //next rank which contains the field (i+1)
                        if(next_rank > NUMT-1){
                            next_rank = NUMT-1;   //if next_rank is bigger than number of threads then the last thread has the remaining values assigned to
                        }
                        //send this value to rank with field in next line (i+1)
                        MPI_Send(&matrix[current_row][current_column], 1, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
                    }
                }

                if(current_row > 0){
                    if(i >= N); //do nothing because the rank that has this value is this rank
                    else{
                        int t = N - i; //how many fields are between the end of this thread and the desired field
                        int r= PARTMATRIX_SIZE;
                        prev_rank = rank - ceil((double)t/(double)r);   //previous rank which contains the field (i-1)

                        //send this value to rank with field in previous line (i-1)
                        MPI_Send(&matrix[current_row][current_column], 1, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD);
                    }
                }
            }

            if(1 < rank){  //the first value also has to be send to the field left to it (j-1) but only if there is another rank before
                MPI_Send(&matrix[start_row][start_column], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            }

            if(rank < NUMT-1){  //the last value also has to be send to the field right to it (j-1) but only if there is another rank behind
                MPI_Send(&matrix[end_row][end_column], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            }
            /////////////end send values


            /////////////receive values
            for (int i = 0; i < my_matrix_size; i++){
                int current_row = (PARTMATRIX_SIZE*(rank-1)+i)/N;   //row and column where the value is saved we are currently working on
                int current_column = (PARTMATRIX_SIZE*(rank-1)+i)%N;

                int prev_rank = 0;
                int next_rank = 0;

                if(current_row < N-1){  //if we are not in the last line, we have to receive the value from the field in the next line
                    //send this value to field in next line (i+1)
                    int t = my_matrix_size - i;     //how many fields from this thread are left
                    if(t > N); //do nothing because the rank that needs this value is this rank
                    else{
                        t = N - t + 1; //how many fields are between the end of this thread and the desired field
                        int r= PARTMATRIX_SIZE;
                        next_rank = rank + ceil((double)t/(double)r);   //next rank which contains the field (i+1)
                        if(next_rank > NUMT-1)
                        {
                            next_rank = NUMT-1;   //if next_rank is bigger than number of threads then the last thread has the remaining values assigned to
                        }
                        MPI_Recv(&matrix[current_row+1][current_column], 1, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }

                if(current_row > 0){
                    if(i >= N); //do nothing because the rank that has this value is this rank
                    else{
                        int t = N - i; //how many fields are between the end of this thread and the desired field
                        int r= PARTMATRIX_SIZE;
                        prev_rank = rank - ceil((double)t/(double)r);   //previous rank which contains the field (i-1)

                        //receive value from field in previous line (i-1)
                        MPI_Recv(&matrix[current_row-1][current_column], 1, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

            if(1 < rank){  //the first value also needs the field left to it (j-1) but only if there is another rank before
                MPI_Recv(&matrix[start_row][start_column-1], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if(rank < NUMT-1){  //the last value also needs the field right to it (j-1) but only if there is another rank behind
                MPI_Recv(&matrix[end_row][end_column+1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            /////////////end receive values

            //calculate residuum
            makeres(&matrix[0][0], &res[0][0], rank, my_matrix_size);

            //make new matrix
			makenew(&matrix[0][0], &res[0][0], rank, my_matrix_size);

            //Get maximal residuum
			max_res = getmaxres(&res[0][0]);

            //send max_res back to master
            MPI_Send(&max_res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            MPI_Recv(&command, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//check if we should continue
        }//end command
        MPI_Send(&matrix[start_row][start_column], my_matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);//when the calculations finished send our values back to master
    }//end slave description

    MPI_Finalize();
    return 0;
}
