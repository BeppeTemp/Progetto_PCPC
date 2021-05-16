#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

//Matrix settings
#define ROW 10
#define COLUMN 10
#define O_PERCENT 33
#define X_PERCENT 33
//Satisfaction setting
#define SAT_THRESHOLD 33.3

//Costants
#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_RESET "\x1b[0m"
#define MASTER 0
#define MAX_RAND_VALUE 99

//Struct definition
typedef struct {
    int obs_sec_size, obs_sec_disp, rel_sec_size, start_rel_sec;
    char *sub_mat;
} Node_Data;

//Matrix generation
char randomValue() {
    int value = rand() % MAX_RAND_VALUE + 1;
    if (value > 0 && value <= O_PERCENT)
        return 'O';
    else if (value > O_PERCENT && value <= (X_PERCENT + O_PERCENT))
        return 'X';
    else
        return ' ';
}
void generateMat(char *mat) {
    srand(time(NULL) + MASTER);
    for (int i = 0; i < (ROW * COLUMN); i++) {
        *(mat + i) = randomValue();
    }
}

//Print funtions
void printVetInt(int vet[], int length) {
    for (int i = 0; i < length; i++) {
        printf("| %d ", vet[i]);
    }
    printf("|\n");
}
void printVetChar(char vet[], int length) {
    for (int i = 0; i < length; i++) {
        if (vet[i] == 'O')
            printf(COLOR_GREEN " %c " COLOR_RESET, vet[i]);
        if (vet[i] == 'X')
            printf(COLOR_RED " %c " COLOR_RESET, vet[i]);
        if (vet[i] == ' ')
            printf("   ");
        printf("|");
    }
    printf("\n");
}
void printMat(char *mat) {
    printf("\n");
    for (int i = 0; i < (ROW * COLUMN); i++) {
        printf("|");
        if (*(mat + i) == 'O')
            printf(COLOR_GREEN " %c " COLOR_RESET, *(mat + i));
        if (*(mat + i) == 'X')
            printf(COLOR_RED " %c " COLOR_RESET, *(mat + i));
        if (*(mat + i) == ' ')
            printf("   ");
        if (((i + 1) % (ROW) == 0) && (i != 0))
            printf("|\n");
    }
    printf("\n");
}
void printSubMat(Node_Data data) {
    int j;
    for (j = 0; j < (data.obs_sec_disp % ROW); j++) printf("    ");
    for (int i = 0; i < data.obs_sec_size; i++) {
        if ((((i + j) % ROW) == 0) && (i != 0)) printf("|\n");
        printf("|");
        if ((i >= data.start_rel_sec) && (i < (data.start_rel_sec + data.rel_sec_size))) {
            if (data.sub_mat[i] == 'O')
                printf(COLOR_GREEN " %c " COLOR_RESET, data.sub_mat[i]);
            if (data.sub_mat[i] == 'X')
                printf(COLOR_RED " %c " COLOR_RESET, data.sub_mat[i]);
            if (data.sub_mat[i] == ' ')
                printf("   ");
        } else if (data.sub_mat[i] == ' ')
            printf(COLOR_YELLOW " _ " COLOR_RESET);
        else
            printf(COLOR_YELLOW " %c " COLOR_RESET, data.sub_mat[i]);
        if (i == data.obs_sec_size - 1) printf("|");
    }
    printf("\n\n");
}

//Divide matrix into worldsize - 1 parts
void calcSize(int world_size, Node_Data node_data[], char mat[], MPI_Request request) {
    int end_row, extra_space, j;
    int position = 0;

    int section = (ROW * COLUMN) / (world_size - 1);
    int difference = (ROW * COLUMN) % (world_size - 1);

    node_data[0].obs_sec_size = 0;
    node_data[0].obs_sec_disp = 0;
    node_data[0].rel_sec_size = 0;
    node_data[0].start_rel_sec = 0;
    node_data[0].sub_mat = NULL;

    for (int i = 1; i < world_size; i++) {
        node_data[i].rel_sec_size = i <= difference ? section + 1 : section;

        end_row = ROW - (node_data[i].rel_sec_size % ROW);
        extra_space = (end_row + (ROW - end_row) + 1);

        if (i == 1) {
            node_data[i].obs_sec_size = node_data[i].rel_sec_size + extra_space;
            node_data[i].start_rel_sec = 0;
            node_data[i].obs_sec_disp = 0;
        } else if (i == (world_size - 1)) {
            node_data[i].obs_sec_size = node_data[i].rel_sec_size + extra_space;
            node_data[i].start_rel_sec = node_data[i].obs_sec_size - node_data[i].rel_sec_size;
            node_data[i].obs_sec_disp = (ROW * COLUMN) - extra_space - node_data[i].rel_sec_size;
        } else {
            node_data[i].obs_sec_size = node_data[i].rel_sec_size + (extra_space);
            node_data[i].obs_sec_disp = position - extra_space;
            if (node_data[i].obs_sec_disp < 0) {
                node_data[i].obs_sec_size += node_data[i].obs_sec_disp + extra_space;
                node_data[i].start_rel_sec = node_data[i].obs_sec_disp + extra_space;
                node_data[i].obs_sec_disp = 0;
            } else {
                node_data[i].obs_sec_size += extra_space;
                node_data[i].start_rel_sec = extra_space;
            }
            //Cut overflow
            if ((node_data[i].obs_sec_disp + node_data[i].obs_sec_size) > (ROW * COLUMN))
                node_data[i].obs_sec_size -= ((node_data[i].obs_sec_disp + node_data[i].obs_sec_size) - (ROW * COLUMN));
        }
        position += node_data[i].rel_sec_size;

        j = 0;
        node_data[i].sub_mat = malloc(sizeof(char) * node_data[i].obs_sec_size);
        for (int k = node_data[i].obs_sec_disp; k < (node_data[i].obs_sec_disp + (node_data[i].obs_sec_size - 1)); k++) {
            *(node_data[i].sub_mat + j) = mat[k];
            j++;
        }

        MPI_Isend(node_data[i].sub_mat, node_data[i].obs_sec_size, MPI_CHAR, i, MASTER, MPI_COMM_WORLD, &request);
    }
    //* printf("Sezioni relative: ");
    //* printVetInt obs_sec_size, world_size);
    //* printf("Displacements relativi: ");
    //* printVetInt(obs_sec_disp, world_size);
    //* printf("Sezioni reali: ");
    //* printVetInt(rel_sec_size, world_size);
    //* printf("Punti di inizio della seguenza reale: ");
    //* printVetInt(start_rel_sec, world_size);
}

//Calculate the degree of satisfaction of a single agent
int isMyKind(int x_index, int y_index, char sub_mat[]) {
    return sub_mat[x_index] == sub_mat[y_index];
}
int calcSat(int obs_index, int rel_index, char sub_mat[]) {
    int neighborhood = 0;
    int my_kynd = 0;
    float perc = 0;

    //TODO prova a ragionare con le singole sotto matrici
    if ((obs_index % ROW) != 0 && ((obs_index + 1) % ROW) != 0) {
        neighborhood += 2;
        my_kynd += isMyKind(rel_index, rel_index - 1, sub_mat);
        my_kynd += isMyKind(rel_index, rel_index + 1, sub_mat);

        if ((obs_index - ROW) >= 0) {
            neighborhood += 3;
            my_kynd += isMyKind(rel_index, rel_index - ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index - ROW) + 1, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index - ROW) - 1, sub_mat);
        }

        if ((obs_index - ROW) <= (ROW * COLUMN)) {
            neighborhood += 3;
            my_kynd += isMyKind(rel_index, rel_index + ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index + ROW) + 1, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index + ROW) - 1, sub_mat);
        }
    } else if ((obs_index % ROW) == 0) {
        neighborhood++;
        my_kynd += isMyKind(rel_index, rel_index + 1, sub_mat);

        if ((obs_index - ROW) >= 0) {
            neighborhood += 2;
            my_kynd += isMyKind(rel_index, rel_index - ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index - ROW) + 1, sub_mat);
        }

        if ((obs_index - ROW) <= (ROW * COLUMN)) {
            neighborhood += 2;
            my_kynd += isMyKind(rel_index, rel_index + ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index + ROW) + 1, sub_mat);
        }
    } else {
        neighborhood++;
        my_kynd += isMyKind(rel_index, rel_index + 1, sub_mat);

        if ((obs_index - ROW) >= 0) {
            neighborhood += 2;
            my_kynd += isMyKind(rel_index, rel_index - ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index - ROW) - 1, sub_mat);
        }

        if ((obs_index - ROW) <= (ROW * COLUMN)) {
            neighborhood += 2;
            my_kynd += isMyKind(rel_index, rel_index + ROW, sub_mat);
            my_kynd += isMyKind(rel_index, (rel_index + ROW) - 1, sub_mat);
        }
    }
    perc = (float)(100 / neighborhood) * my_kynd;

    //* printf("Neighborhood: %d\n", neighborhood);
    //* printf("Kind: %d\n", my_kynd);
    //* printf("Perc: %.1f\n", perc);

    return perc >= SAT_THRESHOLD;
}

void main() {
    char *mat, *sub_mat;
    Node_Data *node_data;

    //MPI initialization
    int world_size, rank;
    double start, end;

    MPI_Status status;
    MPI_Request request;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //MPI_Struct definition
    MPI_Datatype node_data_type;
    int blockcounts[2] = {4, 1};
    MPI_Aint offsets[2] = {0, sizeof(int) * 4};
    MPI_Datatype oldtypes[2] = {MPI_INT, MPI_LONG};
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &node_data_type);
    MPI_Type_commit(&node_data_type);

    //Start computation
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    if (rank == MASTER) {
        //Matrix generation
        mat = malloc(sizeof(char) * ROW * COLUMN);
        generateMat(mat);
        printf("Qui Master 🧑‍🎓, la matrice generata è: \n");
        printMat(mat);

        //Matrix division
        node_data = malloc(sizeof(Node_Data) * world_size);
        calcSize(world_size, node_data, mat, request);
    }

    Node_Data section_data;
    MPI_Scatter(node_data, 1, node_data_type, &section_data, 1, node_data_type, MASTER, MPI_COMM_WORLD);

    if (rank != MASTER) {
        section_data.sub_mat = malloc(sizeof(char) * section_data.obs_sec_size);
        MPI_Recv(section_data.sub_mat, section_data.obs_sec_size, MPI_CHAR, MASTER, MASTER, MPI_COMM_WORLD, &status);

        printf("Qui Slave %d 👨‍🔧, ho ricevuto: \n\n", rank);
        printSubMat(section_data);

        // TODO non lo fa vede cosi hahahah
        /*for (int i = section_data.start_rel_sec; i < (section_data.start_rel_sec + section_data.rel_sec_size); i++) {
            if (section_data.sub_mat[i] != ' ')
                if (calcSat(section_data.obs_sec_disp + i, i, section_data.sub_mat) == 0)
                    printf("In posizione %d se sta" COLOR_RED " na merda !\n" COLOR_RESET, i);
                else
                    printf("In posizione %d se sta" COLOR_GREEN " alla grande bro ! \n" COLOR_RESET, i);
        }*/

        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //Finish computation

    end = MPI_Wtime();
    MPI_Finalize();

    if (rank == MASTER)
        printf("\n\n🕒 Time in ms = %f\n", end - start);
}