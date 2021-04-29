#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

// Matrix settings
#define ROW 5
#define COLUMN 5
#define O_PERCENTAGE 33
#define X_PERCENTAGE 33
//Satisfaction percentage setting
#define SAT_THRESHOLD 33.3

//Color costants
#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_RESET "\x1b[0m"
//Other costants
#define MASTER 0
#define MAX_RAND_VALUE 99

//Structs definition
typedef struct {
    int real_start;
    int real_finish;
    int *section_size;
    int *section_disp;
    char *sub_mat;
} Data;

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
int generateMat(char *mat, int *empty_slots) {
    int j = 0;
    srand(time(NULL) + MASTER);
    for (int i = 0; i < (ROW * COLUMN); i++) {
        mat[i] = randomValue();
        if (mat[i] == ' ') {
            printf("index: %d\n", j);
            empty_slots[j] = i;
            j++;
        }
    }
    empty_slots = realloc(empty_slots, sizeof(int) * j);
    return j;
}

//Print funtions
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
void printVetInt(int vet[], int length) {
    for (int i = 0; i < length; i++) {
        printf("| %d ", vet[i]);
    }
    printf("|\n");
}
void printRows(char vet[], int start, int length, int isYellow) {
    if (isYellow)
        for (int i = start; i < length; i++) {
            printf("|");
            printCharYellow(vet[i]);
            if (((i + 1) % (ROW) == 0) && (i != 0))
                printf("|\n");
        }
    else
        for (int i = start; i < length; i++) {
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
void printSubMat(int world_size, char *mat, int row, int rank) {
    printf("\n");
    if (rank == 0) {
        printRows(mat, 0, (row - 1) * COLUMN, 0);
        printRows(mat, (row * COLUMN) - COLUMN, row * COLUMN, 1);
        printf("\n");
    } else if (rank == world_size - 1) {
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
void calcSizes(int world_size, Data data) {
    int section = ROW / (world_size);
    int difference = ROW % (world_size);

    for (int i = 0; i < world_size; i++) {
        data.section_size[i] = i < difference ? (section + 1) * ROW : section * ROW;
        //printf("Real SS: %d\n", data[i].section_size / ROW);
        data.section_size[i] += ((i == 0) || (i == world_size - 1)) ? COLUMN : COLUMN * 2;

        if (i == 0) {
            data.section_disp[i] = 0;
        } else if (i == world_size - 1)
            data.section_disp[i] = (ROW * COLUMN) - data.section_size[i];
        else {
            data.section_disp[i] = data.section_disp[i - 1] + data.section_size[i - 1] - (COLUMN * 2);
        }
    }

    //* printf("Sezioni generate: ");
    //* printVetInt(section_size, world_size);
    //* printf("Displacements generati: ");
    //* printVetInt(displacement, world_size);
}
int calcStart(Data data, int rank, int world_size) {
    if (rank == 0)
        return 0;
    else
        return COLUMN;
}
int calcFinish(Data data, int rank, int world_size) {
    if (rank == 0)
        return data.section_size[rank] - COLUMN - 1;
    else if (rank == world_size - 1)
        return data.section_size[rank] - 1;
    else
        return data.section_size[rank] - COLUMN - 1;
}

//Calculate the degree of satisfaction of a single agent
void convertIndex(int index, int *row_index, int *col_index) {
    *row_index = index / COLUMN;
    *col_index = index % COLUMN;
}
int isMyKind(int my_index, int x_index, Data data) {
    return data.sub_mat[my_index] == data.sub_mat[x_index];
}
void calcSatCorner(int index, Data data, int *neighborhood, int *my_kynd, int position) {
    *neighborhood += 3;
    switch (position) {
        //Upper left corner
        case 0:
            *my_kynd += isMyKind(index, index + 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN), data);
            *my_kynd += isMyKind(index, (index + COLUMN) + 1, data);
            break;
        //Upper right corner
        case 1:
            *my_kynd += isMyKind(index, index - 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN), data);
            *my_kynd += isMyKind(index, (index + COLUMN) - 1, data);
            break;
        //Lower left corner
        case 2:
            *my_kynd += isMyKind(index, index + 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN), data);
            *my_kynd += isMyKind(index, (index - COLUMN) + 1, data);
            break;
        //Lower right corner
        case 3:
            *my_kynd += isMyKind(index, index - 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN), data);
            *my_kynd += isMyKind(index, (index - COLUMN) - 1, data);
            break;
    }
}
void calcSatEdge(int index, Data data, int *neighborhood, int *my_kynd, int position) {
    *neighborhood += 5;
    switch (position) {
        //Upper edge
        case 0:
            *my_kynd += isMyKind(index, index + 1, data);
            *my_kynd += isMyKind(index, index - 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN), data);
            *my_kynd += isMyKind(index, (index + COLUMN) + 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN) - 1, data);
            break;
        //Left edge
        case 1:
            *my_kynd += isMyKind(index, index + 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN), data);
            *my_kynd += isMyKind(index, (index - COLUMN), data);
            *my_kynd += isMyKind(index, (index + COLUMN) + 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN) + 1, data);
            break;
        //Right edge
        case 2:
            *my_kynd += isMyKind(index, index - 1, data);
            *my_kynd += isMyKind(index, (index + COLUMN), data);
            *my_kynd += isMyKind(index, (index - COLUMN), data);
            *my_kynd += isMyKind(index, (index + COLUMN) - 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN) - 1, data);
            break;
        //Lower edge
        case 3:
            *my_kynd += isMyKind(index, index + 1, data);
            *my_kynd += isMyKind(index, index - 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN), data);
            *my_kynd += isMyKind(index, (index - COLUMN) + 1, data);
            *my_kynd += isMyKind(index, (index - COLUMN) - 1, data);
            break;
    }
}
void calcSatCenter(int index, Data data, int *neighborhood, int *my_kynd) {
    *neighborhood += 8;
    *my_kynd += isMyKind(index, index + 1, data);
    *my_kynd += isMyKind(index, index - 1, data);
    *my_kynd += isMyKind(index, (index - COLUMN), data);
    *my_kynd += isMyKind(index, (index + COLUMN), data);
    *my_kynd += isMyKind(index, (index - COLUMN) + 1, data);
    *my_kynd += isMyKind(index, (index - COLUMN) - 1, data);
    *my_kynd += isMyKind(index, (index + COLUMN) + 1, data);
    *my_kynd += isMyKind(index, (index + COLUMN) - 1, data);
}
int calcSat(int index, Data data, int rank) {
    int neighborhood = 0;
    int my_kynd = 0;
    int row = data.section_size[rank] / COLUMN;

    int row_index;
    int col_index;

    convertIndex(index, &row_index, &col_index);

    if (row_index == 0)
        if (col_index == 0)
            calcSatCorner(index, data, &neighborhood, &my_kynd, 0);  //Upper left corner
        else if (col_index == COLUMN - 1)
            calcSatCorner(index, data, &neighborhood, &my_kynd, 1);  //Upper right corner
        else
            calcSatEdge(index, data, &neighborhood, &my_kynd, 0);  //Upper edge

    else if (row_index == row - 1)
        if (col_index == 0)
            calcSatCorner(index, data, &neighborhood, &my_kynd, 2);  //Lower left corner
        else if (col_index == COLUMN - 1)
            calcSatCorner(index, data, &neighborhood, &my_kynd, 3);  //Lower right corner
        else
            calcSatEdge(index, data, &neighborhood, &my_kynd, 3);  //Lower edge

    else if (col_index == 0 && row_index > 0 && row_index < row - 1)
        calcSatEdge(index, data, &neighborhood, &my_kynd, 1);  //Left edge

    else if (col_index == COLUMN - 1 && row_index > 0 && row_index < row)
        calcSatEdge(index, data, &neighborhood, &my_kynd, 2);  //Right edge
    else
        calcSatCenter(index, data, &neighborhood, &my_kynd);  //Center

    float perc = (100 / (float)neighborhood) * my_kynd;

    //* printf("Rank: %d\n", rank);
    //* printf("Value: %c\n", data.sub_mat[index]);
    //* printf("Index: %d\n", index);
    //* printf("Rel_Row: %d\n", row);
    //* printf("Row_index: %d\n", row_index);
    //* printf("Col_index: %d\n", col_index);
    //* printf("Neighborhood: %d\n", neighborhood);
    //* printf("Kind: %d\n", my_kynd);
    //* printf("Perc: %.1f\n\n", perc);

    return perc >= SAT_THRESHOLD;
}

void main() {
    char *mat, *sub_mat;
    int *empty_slots;
    Data data;

    //MPI initialization
    int world_size, rank;
    double start, end;

    MPI_Status status;
    MPI_Request request;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Matrix generation
    if (rank == MASTER) {
        mat = malloc(ROW * COLUMN * sizeof(char));
        empty_slots = malloc(ROW * COLUMN * sizeof(int));
        int n_empty = generateMat(mat, empty_slots);
        printf("Qui Master 🧑‍🎓, la matrice generata è: \n");
        printMat(mat);

        for (int i = 0; i < n_empty; i++) {
            printf("%d ", empty_slots[i]);
        }
    }

    data.section_size = malloc(sizeof(int) * world_size);
    data.section_disp = malloc(sizeof(int) * world_size);
    calcSizes(world_size, data);

    data.sub_mat = malloc(sizeof(char) * data.section_size[rank]);
    MPI_Scatterv(mat, data.section_size, data.section_disp, MPI_CHAR, data.sub_mat, data.section_size[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);
    if (rank == MASTER) free(mat);

    data.real_start = calcStart(data, rank, world_size);
    data.real_finish = calcFinish(data, rank, world_size);

    //Start computation
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //* for (int i = data.real_start; i <= data.real_finish; i++) {
    //*     if (data.sub_mat[i] != ' ')
    //*         if (calcSat(i, data, rank))
    //*             printf("In posizione %d se sta" COLOR_GREEN " alla grande bro ! \n" COLOR_RESET, i + data.section_disp[rank]);
    //*         else
    //*             printf("In posizione %d se sta" COLOR_RED " na merda !\n" COLOR_RESET, i + data.section_disp[rank]);
    //* }

    //* printf("Qui processo %d 👨‍🔧, ho ricevuto: \n", rank);
    //* printSubMat(world_size, data.sub_mat, data.section_size[rank] / COLUMN, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    //Finish computation

    MPI_Finalize();

    if (rank == MASTER)
        printf("\n\n🕒 Time in ms = %f\n", end - start);
}

//TODO Gestione pochi processi