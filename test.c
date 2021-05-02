#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

///////////////////////////////////////
// Impostazioni della Matrice
///////////////////////////////////////
#define ROW 95
#define COLUMN 95
#define O_PERCENTAGE 33
#define X_PERCENTAGE 33
#define SAT_THRESHOLD 33.3
#define N_ITERACTION 10000
#define ASSIGN_SEED 117
///////////////////////////////////////
// Costants
///////////////////////////////////////
#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_RESET "\x1b[0m"
#define MASTER 0
#define MAX_RAND_VALUE 99
///////////////////////////////////////

//Structs definition
typedef struct {
    //Start and finish id of submat
    int r_start;
    int r_finish;

    //Submat of a specific process
    char *sub_mat;

    //Number of assigned empty slots and id
    int n_my_empty;
    int *my_emp_slots;

    //Generale section and displacement for al processes
    int *sec_size;
    int *sec_disp;

    //Generale section and displacement for gather operation
    int *sec_gt_size;
    int *sec_gt_disp;
} Data;
typedef struct {
    int id_agent;
    int id_reset;
    char vl_agent;
} Move;

// Debug functions
void sampleMat(char *mat) {
    mat[0] = 'O';
    mat[1] = 'O';
    mat[2] = 'O';
    mat[3] = 'O';
    mat[4] = 'O';

    mat[5] = ' ';
    mat[6] = 'O';
    mat[7] = 'O';
    mat[8] = ' ';
    mat[9] = 'O';

    mat[10] = 'X';
    mat[11] = ' ';
    mat[12] = 'O';
    mat[13] = 'X';
    mat[14] = 'O';

    mat[15] = 'O';
    mat[16] = ' ';
    mat[17] = 'X';
    mat[18] = ' ';
    mat[19] = 'X';

    mat[20] = 'X';
    mat[21] = ' ';
    mat[22] = ' ';
    mat[23] = ' ';
    mat[24] = 'O';
}
void syncProcess(unsigned int rank) {
    usleep(rank * 500);
}
void saveMatrixToFile(char *mat) {
    FILE *file = fopen("Matrix.txt", "w");
    int results = fputs(mat, file);
    fclose(file);
}

//Matrix generation
char randomValue() {
    int value = rand() % MAX_RAND_VALUE + 1;
    if (value > 0 && value <= O_PERCENTAGE)
        return 'O';
    else if (value > O_PERCENTAGE && value <= (X_PERCENTAGE + O_PERCENTAGE))
        return 'X';
    else
        return ' ';
}
void generateMat(char *mat) {
    srand(time(NULL) + MASTER);
    for (int i = 0; i < (ROW * COLUMN); i++)
        mat[i] = randomValue();
}

//Print funtions for debugging
void printChar(char x) {
    if (x == 'O')
        printf(COLOR_GREEN " %c " COLOR_RESET, x);
    if (x == 'X')
        printf(COLOR_RED " %c " COLOR_RESET, x);
    if (x == ' ')
        printf("   ");
}
void printCharYellow(char x) {
    if (x == 'O')
        printf(COLOR_YELLOW " %c " COLOR_RESET, x);
    if (x == 'X')
        printf(COLOR_YELLOW " %c " COLOR_RESET, x);
    if (x == ' ')
        printf(COLOR_YELLOW " _ " COLOR_RESET);
}
void printVetInt(int vet[], int len) {
    if (len > 0) {
        for (int i = 0; i < len; i++) {
            printf("| %d ", vet[i]);
        }
        printf("|\n");
    } else
        printf("\n");
}
void printVetChar(char vet[], int len) {
    if (len > 0) {
        for (int i = 0; i < len; i++) {
            printf("| %c ", vet[i]);
        }
        printf("|\n");
    } else
        printf("\n");
}
void printRows(char vet[], int start, int len, int is_Ylw) {
    if (is_Ylw)
        for (int i = start; i < len; i++) {
            printf("|");
            printCharYellow(vet[i]);
            if (((i + 1) % (ROW) == 0) && (i != 0))
                printf("|\n");
        }
    else
        for (int i = start; i < len; i++) {
            printf("|");
            printChar(vet[i]);
            if (((i + 1) % (ROW) == 0) && (i != 0))
                printf("|\n");
        }
}
void printMat(char *mat) {
    printf("\n");
    for (int i = 0; i < (ROW * COLUMN); i++) {
        printf("|");
        printChar(mat[i]);
        if (((i + 1) % (ROW) == 0) && (i != 0))
            printf("|\n");
    }
    printf("\n");
}
void printSubMat(int wd_size, char *mat, int row, int rank) {
    printf("\n");
    if (rank == 0) {
        printRows(mat, 0, (row - 1) * COLUMN, 0);
        printRows(mat, (row * COLUMN) - COLUMN, row * COLUMN, 1);
        printf("\n");
    } else if (rank == wd_size - 1) {
        printRows(mat, 0, COLUMN, 1);
        printRows(mat, COLUMN, row * COLUMN, 0);
        printf("\n");
    } else {
        printRows(mat, 0, COLUMN, 1);
        printRows(mat, COLUMN, (row * COLUMN) - COLUMN, 0);
        printRows(mat, (row * COLUMN) - COLUMN, row * COLUMN, 1);
        printf("\n");
    }
}

//Divide matrix in n row per process
void calcSizes(int wd_size, Data data) {
    int section = ROW / (wd_size);
    int difference = ROW % (wd_size);

    for (int i = 0; i < wd_size; i++) {
        data.sec_size[i] = i < difference ? (section + 1) * ROW : section * ROW;
        data.sec_gt_size[i] = data.sec_size[i];
        data.sec_gt_disp[i] = i == 0 ? 0 : data.sec_gt_disp[i - 1] + data.sec_gt_size[i - 1];
        data.sec_size[i] += ((i == 0) || (i == wd_size - 1)) ? COLUMN : COLUMN * 2;

        if (i == 0) {
            data.sec_disp[i] = 0;
        } else if (i == wd_size - 1)
            data.sec_disp[i] = (ROW * COLUMN) - data.sec_size[i];
        else {
            data.sec_disp[i] = data.sec_disp[i - 1] + data.sec_size[i - 1] - (COLUMN * 2);
        }
    }
    //! Debug //
    /*tf("Sezioni generate: ");
    printVetInt(sec_size, wd_size);
    printf("Displacements generati: ");
    printVetInt(displacement, wd_size);*/
    //! ----- //
}
int calcStart(int rank, int wd_size) {
    if (rank == 0)
        return 0;
    else
        return COLUMN;
}
int calcFinish(Data data, int rank, int wd_size) {
    if (rank == 0)
        return data.sec_size[rank] - COLUMN - 1;
    else if (rank == wd_size - 1)
        return data.sec_size[rank] - 1;
    else
        return data.sec_size[rank] - COLUMN - 1;
}

//Calculate the degree of satisfaction of a single agent
void convertIndex(int id, int *row_index, int *col_index) {
    *row_index = id / COLUMN;
    *col_index = id % COLUMN;
}
int isMyKind(int my_index, int x_index, Data data) {
    return data.sub_mat[my_index] == data.sub_mat[x_index];
}
void calcSatCorner(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 3;
    switch (pos) {
        //Upper left corner
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN), data);
            *my_kynd += isMyKind(id, (id + COLUMN) + 1, data);
            break;
        //Upper right corner
        case 1:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN), data);
            *my_kynd += isMyKind(id, (id + COLUMN) - 1, data);
            break;
        //Lower left corner
        case 2:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN), data);
            *my_kynd += isMyKind(id, (id - COLUMN) + 1, data);
            break;
        //Lower right corner
        case 3:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN), data);
            *my_kynd += isMyKind(id, (id - COLUMN) - 1, data);
            break;
    }
}
void calcSatEdge(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 5;
    switch (pos) {
        //Upper edge
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN), data);
            *my_kynd += isMyKind(id, (id + COLUMN) + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN) - 1, data);
            break;
        //Left edge
        case 1:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN), data);
            *my_kynd += isMyKind(id, (id - COLUMN), data);
            *my_kynd += isMyKind(id, (id + COLUMN) + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN) + 1, data);
            break;
        //Right edge
        case 2:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMN), data);
            *my_kynd += isMyKind(id, (id - COLUMN), data);
            *my_kynd += isMyKind(id, (id + COLUMN) - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN) - 1, data);
            break;
        //Lower edge
        case 3:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN), data);
            *my_kynd += isMyKind(id, (id - COLUMN) + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMN) - 1, data);
            break;
    }
}
void calcSatCenter(int id, Data data, int *neigh, int *my_kynd) {
    *neigh += 8;
    *my_kynd += isMyKind(id, id + 1, data);
    *my_kynd += isMyKind(id, id - 1, data);
    *my_kynd += isMyKind(id, (id - COLUMN), data);
    *my_kynd += isMyKind(id, (id + COLUMN), data);
    *my_kynd += isMyKind(id, (id - COLUMN) + 1, data);
    *my_kynd += isMyKind(id, (id - COLUMN) - 1, data);
    *my_kynd += isMyKind(id, (id + COLUMN) + 1, data);
    *my_kynd += isMyKind(id, (id + COLUMN) - 1, data);
}
int calcSat(int id, Data data, int rank) {
    int neigh = 0;
    int my_kynd = 0;
    int row = data.sec_size[rank] / COLUMN;

    int row_index;
    int col_index;

    convertIndex(id, &row_index, &col_index);

    if (row_index == 0)
        if (col_index == 0)
            calcSatCorner(id, data, &neigh, &my_kynd, 0);  //Upper left corner
        else if (col_index == COLUMN - 1)
            calcSatCorner(id, data, &neigh, &my_kynd, 1);  //Upper right corner
        else
            calcSatEdge(id, data, &neigh, &my_kynd, 0);  //Upper edge

    else if (row_index == row - 1)
        if (col_index == 0)
            calcSatCorner(id, data, &neigh, &my_kynd, 2);  //Lower left corner
        else if (col_index == COLUMN - 1)
            calcSatCorner(id, data, &neigh, &my_kynd, 3);  //Lower right corner
        else
            calcSatEdge(id, data, &neigh, &my_kynd, 3);  //Lower edge

    else if (col_index == 0 && row_index > 0 && row_index < row - 1)
        calcSatEdge(id, data, &neigh, &my_kynd, 1);  //Left edge

    else if (col_index == COLUMN - 1 && row_index > 0 && row_index < row)
        calcSatEdge(id, data, &neigh, &my_kynd, 2);  //Right edge
    else
        calcSatCenter(id, data, &neigh, &my_kynd);  //Center

    float perc = (100 / (float)neigh) * my_kynd;

    //! Debug //
    /*printf("Rank: %d\n", rank);
    printf("Value: %c\n", data.sub_mat[id]);
    printf("Index: %d\n", id);
    printf("Rel_Row: %d\n", row);
    printf("Row_index: %d\n", row_index);
    printf("Col_index: %d\n", col_index);
    printf("Neighborhood: %d\n", neigh);
    printf("Kind: %d\n", my_kynd);
    printf("Perc: %.1f\n\n", perc);*/
    //! ----- //

    return perc >= SAT_THRESHOLD;
}

//Empty slots management
void shuffle(int *vet, int length) {
    srand(ASSIGN_SEED);

    for (int i = 0; i < length - 1; i++) {
        size_t j = i + rand() / (RAND_MAX / (length - i) + 1);
        int t = vet[j];
        vet[j] = vet[i];
        vet[i] = t;
    }
}
int calcEmptySlots(Data data, int *my_emp_slots, int rank, int wd_size) {
    int empty_tot = 0;
    int vet_siz[wd_size];
    int vet_disp[wd_size];

    //Individuazione degli slot vuoti
    int *vet_emp = malloc(sizeof(int) * data.r_finish - data.r_start);
    int n_my_emp = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] == ' ') {
            vet_emp[n_my_emp] = data.sec_disp[rank] + i;
            n_my_emp++;
        }
    }

    MPI_Allgather(&n_my_emp, 1, MPI_INT, vet_siz, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < wd_size; i++) {
        empty_tot += vet_siz[i];
        vet_disp[i] = i == 0 ? 0 : vet_disp[i - 1] + vet_siz[i - 1];
    }

    int *emp_slots = malloc(sizeof(int) * empty_tot);

    MPI_Allgatherv(vet_emp, n_my_emp, MPI_INT, emp_slots, vet_siz, vet_disp, MPI_INT, MPI_COMM_WORLD);
    free(vet_emp);

    //Randomizzazione dell'array
    for (int i = 0; i < wd_size; i++) {
        shuffle(emp_slots, empty_tot);
    }

    //Assegnazione degli slot vuoti
    data.n_my_empty = empty_tot / wd_size;
    int k = rank * data.n_my_empty;

    for (int i = 0; i < data.n_my_empty; i++) {
        my_emp_slots[i] = emp_slots[k];
        k++;
    }

    return data.n_my_empty;
}

//Move operation
void findMoves(Data data, Move *my_moves, int rank) {
    int k = 0;
    int n_moves = data.n_my_empty;

    //Identificazione agenti da spostare
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] != ' ')
            if (!calcSat(i, data, rank)) {
                my_moves[k].id_agent = data.my_emp_slots[k];
                my_moves[k].id_reset = data.sec_disp[rank] + i;
                my_moves[k].vl_agent = data.sub_mat[i];
                if (n_moves-- == 0) break;
                k++;
            } else {
                my_moves[k].id_agent = -1;
                my_moves[k].id_reset = -1;
                my_moves[k].vl_agent = 'n';
                if (n_moves-- == 0) break;
                k++;
            }
    }
}
int calcMembership(Data data, int rank, int index) {
    int start = data.sec_disp[rank];
    int finish = data.sec_disp[rank] + data.sec_size[rank];

    return (index >= start && index < finish);
}
void move(Data data, Move *my_moves, Move *moves, MPI_Datatype move_data_type, int wd_size, int rank, int n_itc) {
    int sec_size[wd_size];
    int sec_disp[wd_size];

    for (int i = 0; i < wd_size; i++) {
        sec_size[i] = data.n_my_empty;
        sec_disp[i] = i == 0 ? 0 : sec_disp[i - 1] + sec_size[i - 1];
    }

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
}

//Stop operation
int isOver(int n_moves, int wd_size, int rank) {
    int all_n_moves[wd_size];

    MPI_Allgather(&n_moves, 1, MPI_INT, all_n_moves, 1, MPI_INT, MPI_COMM_WORLD);

    int sum = 0;
    for (int i = 0; i < wd_size; i++) sum += all_n_moves[i];

    return sum == 0;
}
void gatherResult(Data data, int rank, char *mat) {
    char *test = malloc(sizeof(char) * data.sec_gt_size[rank]);
    int k = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        test[k] = data.sub_mat[i];
        k++;
    }
    MPI_Gatherv(test, data.sec_gt_size[rank], MPI_CHAR, mat, data.sec_gt_size, data.sec_gt_disp, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    free(test);
}

void main() {
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

    if (wd_size <= ROW) {
        //Matrix generation
        if (rank == MASTER) {
            mat = malloc(ROW * COLUMN * sizeof(char));
            generateMat(mat);
            //sampleMat(mat);
            saveMatrixToFile(mat);
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

        //Calcolo start and finish
        data.r_start = calcStart(rank, wd_size);
        data.r_finish = calcFinish(data, rank, wd_size);

        data.my_emp_slots = malloc(sizeof(int) * ROW * COLUMN);
        Move *moves = malloc(sizeof(Move));

        //Start computation
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while (n_itc != 0) {
            //Empty slot calculation
            data.n_my_empty = calcEmptySlots(data, data.my_emp_slots, rank, wd_size);

            if (n_itc == N_ITERACTION) {
                data.my_emp_slots = realloc(data.my_emp_slots, sizeof(int) * data.n_my_empty);
                moves = realloc(moves, sizeof(Move) * data.n_my_empty * wd_size);
            }

            //Identificazione degli agenti da spostare
            Move my_moves[data.n_my_empty];
            findMoves(data, my_moves, rank);

            //Move agent
            move(data, my_moves, moves, move_data_type, wd_size, rank, n_itc);

            n_itc--;
        }

        gatherResult(data, rank, mat);

        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        //Finish computation
    }

    MPI_Finalize();

    if (wd_size <= ROW) {
        if (rank == MASTER) {
            printf("Qui Master 🧑‍🎓, la matrice elaborata è: \n");
            printMat(mat);
            printf("Iterazioni effettuate: %d\n", N_ITERACTION - n_itc);
            printf("\n\n🕒 Time in ms = %f\n", end - start);
        }

        free(data.sec_size);
        free(data.sec_disp);
        free(data.sec_gt_size);
        free(data.sec_gt_disp);
        free(data.sub_mat);
        free(data.my_emp_slots);
    } else if (rank == MASTER)
        printf("Qui Master 🧑‍🎓, computazione impossibile, numero di slave eccessivo. \n");
}