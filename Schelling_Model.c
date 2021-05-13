#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

//////////////////////////////////
// Computation settings
//////////////////////////////////
#define SIZE 5
#define O_PERCENTAGE 33
#define X_PERCENTAGE 33
#define SAT_THRESHOLD 33.3
#define N_ITERACTION 100
#define ASSIGN_SEED 117
//////////////////////////////////
// Costants
//////////////////////////////////
#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_RESET "\x1b[0m"
#define MASTER 0
#define MAX_RAND_VALUE 99
///////////////////////////////////
// Definition of structures
///////////////////////////////////
typedef struct {
    //Assigned start and finish submatrix
    int r_start;
    int r_finish;

    //Assigned submatrix
    char *sub_mat;

    //Number of assigned empty locations and locations
    int n_my_empty;
    int *my_emp_slots;

    //Information for the subdivision scatter
    int *sec_size;
    int *sec_disp;

    //Information for the subdivision gather
    int *sec_gt_size;
    int *sec_gt_disp;
} Data;
typedef struct {
    //Arrival index of the moving agent
    int id_agent;

    //Starting index of the moving agent
    int id_reset;

    //Type of agent on the move
    char vl_agent;
} Move;
///////////////////////////////////

//? Debug functions //
//Fills the matrix with a constant set of values
void sampleMat(char *mat) {
    mat[0] = 'X';
    mat[1] = 'X';
    mat[2] = 'O';
    mat[3] = ' ';
    mat[4] = 'X';

    mat[5] = 'O';
    mat[6] = 'O';
    mat[7] = 'O';
    mat[8] = 'O';
    mat[9] = ' ';

    mat[10] = 'X';
    mat[11] = 'X';
    mat[12] = 'O';
    mat[13] = 'O';
    mat[14] = ' ';

    mat[15] = 'X';
    mat[16] = 'X';
    mat[17] = 'O';
    mat[18] = 'X';
    mat[19] = 'X';

    mat[20] = ' ';
    mat[21] = 'X';
    mat[22] = ' ';
    mat[23] = ' ';
    mat[24] = 'X';
}
//Inserts delays in computation to synchronize the processes
void syncProcess(unsigned int rank) {
    usleep(rank * 500);
}

//? Matrix generation functions //
//Converts a random integer value to an agent type
char randomValue() {
    int value = rand() % MAX_RAND_VALUE + 1;
    if (value > 0 && value <= O_PERCENTAGE)
        return 'O';
    else if (value > O_PERCENTAGE && value <= (X_PERCENTAGE + O_PERCENTAGE))
        return 'X';
    else
        return ' ';
}
//Randomly fills the matrix
void generateMat(char *mat) {
    srand(time(NULL) + MASTER);
    for (int i = 0; i < (SIZE * SIZE); i++)
        mat[i] = randomValue();
}

//? Functions for printing //
//Print a single agent
void printChar(char x) {
    if (x == 'O')
        printf(COLOR_GREEN " %c " COLOR_RESET, x);
    if (x == 'X')
        printf(COLOR_RED " %c " COLOR_RESET, x);
    if (x == ' ')
        printf("   ");
}
//Print the entire matrix
void printMat(char *mat) {
    printf("\n");
    for (int i = 0; i < (SIZE * SIZE); i++) {
        printf("|");
        printChar(mat[i]);
        if (((i + 1) % (SIZE) == 0) && (i != 0))
            printf("|\n");
    }
    printf("\n");
}

//? Functions for dividing the matrix between processes //
//Calculate section sizes and displacements for gather and scatter operations
void calcSizes(int wd_size, Data data) {
    int section = SIZE / (wd_size);
    int difference = SIZE % (wd_size);

    for (int i = 0; i < wd_size; i++) {
        data.sec_size[i] = i < difference ? (section + 1) * SIZE : section * SIZE;
        data.sec_gt_size[i] = data.sec_size[i];
        data.sec_gt_disp[i] = i == 0 ? 0 : data.sec_gt_disp[i - 1] + data.sec_gt_size[i - 1];
        data.sec_size[i] += ((i == 0) || (i == wd_size - 1)) ? SIZE : SIZE * 2;

        if (i == 0) {
            data.sec_disp[i] = 0;
        } else if (i == wd_size - 1)
            data.sec_disp[i] = (SIZE * SIZE) - data.sec_size[i];
        else {
            data.sec_disp[i] = data.sec_disp[i - 1] + data.sec_size[i - 1] - (SIZE * 2);
        }
    }
}
//Calculate start position of real submatrix
int calcStart(int rank, int wd_size) {
    if (rank == 0)
        return 0;
    else
        return SIZE;
}
//Calculate finish position of real submatrix
int calcFinish(Data data, int rank, int wd_size) {
    if (rank == 0)
        return data.sec_size[rank] - SIZE - 1;
    else if (rank == wd_size - 1)
        return data.sec_size[rank] - 1;
    else
        return data.sec_size[rank] - SIZE - 1;
}

//? Functions for calculating the degree of satisfaction
//Converts the index of a vector to its matrix counterpart
void convertIndex(int id, int *row_index, int *col_index) {
    *row_index = id / SIZE;
    *col_index = id % SIZE;
}
//Check if two agents are of the same type
int isMyKind(int my_index, int x_index, Data data) {
    return data.sub_mat[my_index] == data.sub_mat[x_index];
}
//Check the satisfaction of an agent in the corner
void satCorner(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 3;
    switch (pos) {
        //Top left corner
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            break;
        //Top right corner
        case 1:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            break;
        //Bottom left corner
        case 2:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            break;
        //Bottom right corner
        case 3:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
    }
}
//Check the satisfaction of an agent in the edge
void satEdge(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 5;
    switch (pos) {
        //Top edge
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            break;
        //Left edge
        case 1:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            break;
        //Right edge
        case 2:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
        //Bottom edge
        case 3:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
    }
}
//Check the satisfaction of a central agent
void satCenter(int id, Data data, int *neigh, int *my_kynd) {
    *neigh += 8;
    *my_kynd += isMyKind(id, id + 1, data);
    *my_kynd += isMyKind(id, id - 1, data);
    *my_kynd += isMyKind(id, (id - SIZE), data);
    *my_kynd += isMyKind(id, (id + SIZE), data);
    *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
    *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
    *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
    *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
}
//Calculate the satisfaction of an agent
int calcSat(int id, Data data, int rank) {
    int neigh = 0, my_kynd = 0;
    int row = data.sec_size[rank] / SIZE;
    int row_index, col_index;

    convertIndex(id, &row_index, &col_index);

    if (row_index == 0)
        if (col_index == 0)
            satCorner(id, data, &neigh, &my_kynd, 0);  //Upper left corner
        else if (col_index == SIZE - 1)
            satCorner(id, data, &neigh, &my_kynd, 1);  //Upper right corner
        else
            satEdge(id, data, &neigh, &my_kynd, 0);  //Upper edge

    else if (row_index == row - 1)
        if (col_index == 0)
            satCorner(id, data, &neigh, &my_kynd, 2);  //Lower left corner
        else if (col_index == SIZE - 1)
            satCorner(id, data, &neigh, &my_kynd, 3);  //Lower right corner
        else
            satEdge(id, data, &neigh, &my_kynd, 3);  //Lower edge

    else if (col_index == 0 && row_index > 0 && row_index < row - 1)
        satEdge(id, data, &neigh, &my_kynd, 1);  //Left edge

    else if (col_index == SIZE - 1 && row_index > 0 && row_index < row)
        satEdge(id, data, &neigh, &my_kynd, 2);  //Right edge
    else
        satCenter(id, data, &neigh, &my_kynd);  //Center

    float perc = (100 / (float)neigh) * my_kynd;
    return perc >= SAT_THRESHOLD;
}

//? Free locations management functions
//Randomizes the position of values in a vector
void shuffle(int *vet, int length) {
    for (int i = 0; i < length - 1; i++) {
        size_t j = i + rand() / (RAND_MAX / (length - i) + 1);
        int t = vet[j];
        vet[j] = vet[i];
        vet[i] = t;
    }
}
//Locate all empty locations in submatrices and associates them
int calcEmptySlots(Data data, int *my_emp_slots, int rank, int wd_size) {
    int empty_tot = 0;
    int vet_siz[wd_size];
    int vet_disp[wd_size];

    // Find free locations in the local sub-matrix
    int vet_emp[data.r_finish - data.r_start];
    int n_my_emp = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] == ' ') {
            vet_emp[n_my_emp] = data.sec_disp[rank] + i;
            n_my_emp++;
        }
    }

    //Gather of global free locations
    MPI_Allgather(&n_my_emp, 1, MPI_INT, vet_siz, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < wd_size; i++) {
        empty_tot += vet_siz[i];
        vet_disp[i] = i == 0 ? 0 : vet_disp[i - 1] + vet_siz[i - 1];
    }
    int emp_slots[empty_tot];
    MPI_Allgatherv(vet_emp, n_my_emp, MPI_INT, emp_slots, vet_siz, vet_disp, MPI_INT, MPI_COMM_WORLD);

    //Randomizes the array of free locations
    srand(ASSIGN_SEED);
    for (int i = 0; i < rand() % 10; i++) {
        shuffle(emp_slots, empty_tot);
    }

    //Assignment of free locations to processes
    data.n_my_empty = empty_tot / wd_size;
    int k = rank * data.n_my_empty;
    for (int i = 0; i < data.n_my_empty; i++) {
        my_emp_slots[i] = emp_slots[k];
        k++;
    }

    return data.n_my_empty;
}

//? Functions for moving agents
//Locates the agents to be moved in the submatrix
int findMoves(Data data, Move *my_moves, int rank, int wd_size) {
    int k = 0, is_over = 0;
    int n_moves = data.n_my_empty;

    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] != ' ')
            if (!calcSat(i, data, rank)) {
                //if (rank == 0) printf("ins: %d\n", i);
                my_moves[k].id_agent = data.my_emp_slots[k];
                if(rank == 0) printf("|%d|", data.my_emp_slots[n_moves]);
                my_moves[k].id_reset = data.sec_disp[rank] + i;
                my_moves[k].vl_agent = data.sub_mat[i];
                is_over++;
                if (n_moves-- == 0) break;
                k++;
            }
    }
    if (rank == 0) printf("\n");
    


    for (int i = 0; i < n_moves; i++) {
        my_moves[k].id_agent = -1;
        my_moves[k].id_reset = -1;
        my_moves[k].vl_agent = 'n';
        k++;
    }

    //Check if all processes have finished their movements
    int tot_over[wd_size];
    MPI_Allgather(&is_over, 1, MPI_INT, tot_over, 1, MPI_INT, MPI_COMM_WORLD);
    is_over = 0;

    if (rank == 0){
        for (int i = 0; i < wd_size; i++) {
            printf("%d\n", tot_over[i]);
        }
        printf("\n");
    }

    for (int i = 0; i < wd_size; i++) is_over += tot_over[i];

    return (is_over == 0);
}
//Given an index and a process checks if it is associated with it
int calcMembership(Data data, int rank, int index) {
    int start = data.sec_disp[rank];
    int finish = data.sec_disp[rank] + data.sec_size[rank];

    return (index >= start && index < finish);
}
//Carries out the movement operations of the agents
void move(Data data, Move *my_moves, MPI_Datatype move_data_type, int wd_size, int rank) {
    int sec_size[wd_size];
    int sec_disp[wd_size];

    for (int i = 0; i < wd_size; i++) {
        sec_size[i] = data.n_my_empty;
        sec_disp[i] = i == 0 ? 0 : sec_disp[i - 1] + sec_size[i - 1];
    }

    Move *moves = malloc(sizeof(Move) * data.n_my_empty * wd_size);
    MPI_Allgatherv(my_moves, data.n_my_empty, move_data_type, moves, sec_size, sec_disp, move_data_type, MPI_COMM_WORLD);

    for (int i = 0; i < data.n_my_empty * wd_size; i++) {
        if (moves[i].id_agent != -1) {
            if (calcMembership(data, rank, moves[i].id_agent)) {
                data.sub_mat[moves[i].id_agent - data.sec_disp[rank]] = moves[i].vl_agent;
            }
            if (calcMembership(data, rank, moves[i].id_reset)) {
                data.sub_mat[moves[i].id_reset - data.sec_disp[rank]] = ' ';
            }
        }
    }
    free(moves);
}
//Ends execution if no process has dissatisfied agents
void gatherResult(Data data, int rank, char *mat) {
    char section[data.sec_gt_size[rank]];
    int k = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        section[k] = data.sub_mat[i];
        k++;
    }
    MPI_Gatherv(section, data.sec_gt_size[rank], MPI_CHAR, mat, data.sec_gt_size, data.sec_gt_disp, MPI_CHAR, MASTER, MPI_COMM_WORLD);
}

void main() {
    //Variable definitions
    int n_itc = N_ITERACTION;
    char *mat;
    Data data;

    //MPI initialization
    int wd_size, rank;
    double start, end;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wd_size);

    //MPI_Struct definition
    MPI_Datatype move_data_type;
    int blockcounts[2] = {2, 1};
    MPI_Aint offsets[2] = {0, sizeof(int) * 2};
    MPI_Datatype oldtypes[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &move_data_type);
    MPI_Type_commit(&move_data_type);

    if (wd_size <= SIZE) {
        //Matrix generation
        if (rank == MASTER) {
            mat = malloc(SIZE * SIZE * sizeof(char));
            //generateMat(mat);
            sampleMat(mat);
            printf("Qui Master 🧑‍🎓, la matrice generata è: \n");
            printMat(mat);
        }

        //Matrix division
        data.sec_size = malloc(sizeof(int) * wd_size);
        data.sec_disp = malloc(sizeof(int) * wd_size);
        data.sec_gt_size = malloc(sizeof(int) * wd_size);
        data.sec_gt_disp = malloc(sizeof(int) * wd_size);
        calcSizes(wd_size, data);
        data.sub_mat = malloc(sizeof(char) * data.sec_size[rank]);
        MPI_Scatterv(mat, data.sec_size, data.sec_disp, MPI_CHAR, data.sub_mat, data.sec_size[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);

        //Start and finish calculation
        data.r_start = calcStart(rank, wd_size);
        data.r_finish = calcFinish(data, rank, wd_size);

        data.my_emp_slots = malloc(sizeof(int) * SIZE * SIZE);

        //Start computation
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while (n_itc != 0) {
            //Empty slot calculation
            data.n_my_empty = calcEmptySlots(data, data.my_emp_slots, rank, wd_size);

            if (n_itc == N_ITERACTION) {
                //End of computation test
                data.my_emp_slots = realloc(data.my_emp_slots, sizeof(int) * data.n_my_empty);
            }

            //Identification of the agents to be moved
            Move *my_moves = malloc(sizeof(Move) * data.n_my_empty);
            if (findMoves(data, my_moves, rank, wd_size)) break;

            //Mobilization of agents
            move(data, my_moves, move_data_type, wd_size, rank);

            free(my_moves);
            n_itc--;
        }

        //Creation of the results matrix
        gatherResult(data, rank, mat);

        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        //Finish computation
    }

    MPI_Finalize();

    //Results display
    if (wd_size <= SIZE) {
        if (rank == MASTER) {
            printf("Qui Master 🧑‍🎓, la matrice elaborata è: \n");
            printMat(mat);
            printf("Iterazioni effettuate: %d\n", N_ITERACTION - n_itc);
            printf("\n\n🕒 Time in ms = %f\n", end - start);

            free(mat);
        }
    } else if (rank == MASTER)
        printf("Qui Master 🧑‍🎓, computazione impossibile, numero di slave eccessivo. \n");
}