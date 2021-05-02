#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

///////////////////////////////////////
// Impostazioni della Matrice
///////////////////////////////////////
#define ROW 150
#define COLUMN 150
#define O_PERCENTAGE 33
#define X_PERCENTAGE 33
#define SAT_THRESHOLD 33.3
#define N_ITERACTION 500
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

    //Number and empty slots in general matrix
    int n_empty;
    int *emp_slots;

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
    //Numero di celle da aggiornare
    int n_agents;

    //Index delle celle da aggiornare e valore
    int *id_agents;
    char *vl_agents;

    //index delle celle da azzerare
    int *id_reset;
} Message;

// Debug functions
void sampleMat55(char *mat) {
    mat[0] = 'O';
    mat[1] = ' ';
    mat[2] = ' ';
    mat[3] = ' ';
    mat[4] = 'X';

    mat[5] = 'O';
    mat[6] = ' ';
    mat[7] = ' ';
    mat[8] = ' ';
    mat[9] = 'X';

    mat[10] = ' ';
    mat[11] = ' ';
    mat[12] = 'X';
    mat[13] = 'X';
    mat[14] = ' ';

    mat[15] = 'O';
    mat[16] = ' ';
    mat[17] = 'O';
    mat[18] = ' ';
    mat[19] = 'X';

    mat[20] = ' ';
    mat[21] = 'O';
    mat[22] = 'O';
    mat[23] = ' ';
    mat[24] = ' ';
}
void sampleMat1010(char *mat) {
    mat[0] = 'X';
    mat[1] = ' ';
    mat[2] = 'X';
    mat[3] = 'O';
    mat[4] = 'X';
    mat[5] = 'X';
    mat[6] = ' ';
    mat[7] = ' ';
    mat[8] = ' ';
    mat[9] = ' ';

    mat[10] = ' ';
    mat[11] = ' ';
    mat[12] = ' ';
    mat[13] = 'X';
    mat[14] = 'O';
    mat[15] = 'X';
    mat[16] = 'O';
    mat[17] = 'X';
    mat[18] = 'X';
    mat[19] = ' ';

    mat[20] = ' ';
    mat[21] = 'O';
    mat[22] = ' ';
    mat[23] = 'X';
    mat[24] = 'O';
    mat[25] = 'X';
    mat[26] = 'X';
    mat[27] = 'X';
    mat[28] = 'O';
    mat[29] = 'X';

    mat[30] = ' ';
    mat[31] = 'X';
    mat[32] = 'O';
    mat[33] = 'O';
    mat[34] = 'X';
    mat[35] = ' ';
    mat[36] = 'X';
    mat[37] = 'O';
    mat[38] = ' ';
    mat[39] = 'X';

    mat[40] = 'O';
    mat[41] = 'X';
    mat[42] = 'O';
    mat[43] = 'O';
    mat[44] = 'O';
    mat[45] = 'X';
    mat[46] = 'X';
    mat[47] = 'O';
    mat[48] = 'O';
    mat[49] = ' ';

    mat[50] = 'O';
    mat[51] = ' ';
    mat[52] = 'O';
    mat[53] = ' ';
    mat[54] = 'X';
    mat[55] = 'O';
    mat[56] = 'X';
    mat[57] = ' ';
    mat[58] = 'X';
    mat[59] = 'X';

    mat[60] = 'O';
    mat[61] = 'X';
    mat[62] = 'O';
    mat[63] = 'X';
    mat[64] = 'X';
    mat[65] = 'X';
    mat[66] = 'X';
    mat[67] = 'O';
    mat[68] = ' ';
    mat[69] = 'O';

    mat[70] = 'X';
    mat[71] = 'O';
    mat[72] = ' ';
    mat[73] = ' ';
    mat[74] = 'O';
    mat[75] = 'O';
    mat[76] = 'O';
    mat[77] = 'X';
    mat[78] = 'O';
    mat[79] = 'O';

    mat[80] = 'O';
    mat[81] = 'O';
    mat[82] = 'O';
    mat[83] = 'X';
    mat[84] = 'O';
    mat[85] = 'X';
    mat[86] = 'X';
    mat[87] = 'X';
    mat[88] = 'X';
    mat[89] = 'O';

    mat[90] = ' ';
    mat[91] = ' ';
    mat[92] = 'X';
    mat[93] = 'O';
    mat[94] = 'O';
    mat[95] = 'O';
    mat[96] = 'X';
    mat[97] = ' ';
    mat[98] = 'O';
    mat[99] = 'X';
}
void syncProcess(unsigned int rank) {
    usleep(rank * 100);
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
int calcEmptySlots(Data data, int *e_slots, int rank, int wd_size) {
    int empty_tot = 0;
    int vet_siz[wd_size];
    int vet_disp[wd_size];

    int vet_emp[ROW * COLUMN];
    int my_emp = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] == ' ') {
            vet_emp[my_emp] = data.sec_disp[rank] + i;
            my_emp++;
        }
    }

    MPI_Allgather(&my_emp, 1, MPI_INT, vet_siz, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < wd_size; i++) {
        empty_tot += vet_siz[i];
        vet_disp[i] = i == 0 ? 0 : vet_disp[i - 1] + vet_siz[i - 1];
    }

    MPI_Allgatherv(vet_emp, my_emp, MPI_INT, e_slots, vet_siz, vet_disp, MPI_INT, MPI_COMM_WORLD);

    return empty_tot;
}
int assignEmptySlots(Data data, int *m_e_slots, int rank, int wd_size) {
    int n_my_slots[wd_size];
    for (int i = 0; i < wd_size; i++) {
        n_my_slots[i] = data.n_empty / wd_size;
        shuffle(data.emp_slots, data.n_empty);
    }

    int k = rank * n_my_slots[rank];

    for (int i = 0; i < n_my_slots[rank]; i++) {
        m_e_slots[i] = data.emp_slots[k];
        k++;
    }

    return n_my_slots[rank];
}

//Move operation
int findMoves(Data data, int rank, int *moves) {
    int n_moves = data.n_my_empty;
    int k = 0;

    //Identificazione agenti da spostare
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] != ' ')
            if (!calcSat(i, data, rank)) {
                moves[k] = data.sec_disp[rank] + i;
                if (n_moves-- == 0) break;
                k++;
            }
    }

    return k;
}
int calcMembership(Data data, int rank, int emp_pos, int wd_size) {
    int start = data.sec_disp[rank];
    int finish = data.sec_disp[rank] + data.sec_size[rank];

    return (emp_pos >= start && emp_pos < finish);
}
void move(Data data, Message msg, int rank) {
    for (int i = 0; i < msg.n_agents; i++) {
        data.sub_mat[msg.id_agents[i] - data.sec_disp[rank]] = msg.vl_agents[i];

        //! Debug //
        /*printf("rank: %d\n", rank);
        printf("Indice reset: %d\n", msg.id_reset[i]);
        printf("Indice convertito: %d\n", msg.id_reset[i] - data.sec_disp[rank]);
        printf("\n");*/
        //! ----- //

        if (msg.id_reset[i] >= data.sec_disp[rank] && msg.id_reset[i] < (data.sec_disp[rank] + data.sec_size[rank]))
            data.sub_mat[msg.id_reset[i] - data.sec_disp[rank]] = ' ';
    }
}
void resetLocalIndex(Data data, Message msg, int rank) {
    for (int i = 0; i < msg.n_agents; i++)
        if (msg.id_reset[i] >= data.sec_disp[rank] && msg.id_reset[i] < (data.sec_disp[rank] + data.sec_size[rank]))
            data.sub_mat[msg.id_reset[i] - data.sec_disp[rank]] = ' ';
}
void moveAgents(Data data, int *moves, int n_moves, int rank, int wd_size, MPI_Request request, MPI_Status status) {
    Message all_msg[wd_size];

    //Messages inizialization
    for (int i = 0; i < wd_size; i++) {
        all_msg[i].n_agents = 0;
        all_msg[i].id_agents = malloc(sizeof(int) * n_moves);
        all_msg[i].vl_agents = malloc(sizeof(char) * n_moves);
        all_msg[i].id_reset = malloc(sizeof(int) * n_moves);
    }

    //Messages creation
    for (int i = 0; i < wd_size; i++) {
        for (int j = 0; j < n_moves; j++) {
            if (calcMembership(data, i, data.my_emp_slots[j], wd_size)) {
                all_msg[i].id_agents[all_msg[i].n_agents] = data.my_emp_slots[j];
                all_msg[i].vl_agents[all_msg[i].n_agents] = data.sub_mat[moves[j] - data.sec_disp[rank]];
                all_msg[i].id_reset[all_msg[i].n_agents] = moves[j];
                all_msg[i].n_agents++;
            }
        }
    }

    //! Debug //
    /*printf("I messaggi che manderò: \n\n");
    for (int i = 0; i < wd_size; i++) {
        if (i != rank) {
            if (all_msg[i].n_agents != 0) {
                printf("Al rank: %d\n", i);
                printf("N_element: %d\n", all_msg[i].n_agents);
                printf("Indici di arrivo: ");
                printVetInt(all_msg[i].id_agents, all_msg[i].n_agents);
                printf("Valori da scrivere: ");
                printVetChar(all_msg[i].vl_agents, all_msg[i].n_agents);
                printf("Indici da azzerare: ");
                printVetInt(all_msg[i].id_reset, all_msg[i].n_agents);
            } else {
                printf("Al rank %d niente.\n", i);
            }
        } else {
            if (all_msg[i].n_agents != 0) {
                printf("A me stesso\n");
                printf("Movimenti: %d\n", all_msg[i].n_agents);
                printf("Indici di arrivo: ");
                printVetInt(all_msg[i].id_agents, all_msg[i].n_agents);
                printf("Valori da spostare: ");
                printVetChar(all_msg[i].vl_agents, all_msg[i].n_agents);
                printf("Indici da azzerare: ");
                printVetInt(all_msg[i].id_reset, all_msg[i].n_agents);
            } else {
                printf("A me stesso niente.\n");
            }
        }
        printf("\n");
    }*/
    //! ----- //

    //Messages send
    char message[BUFSIZ];
    int position;
    for (int i = 0; i < wd_size; i++) {
        if (i != rank) {
            resetLocalIndex(data, all_msg[i], rank);

            //External moves
            position = 0;
            MPI_Pack(&all_msg[i].n_agents, 1, MPI_INT, message, BUFSIZ, &position, MPI_COMM_WORLD);
            MPI_Pack(all_msg[i].id_agents, all_msg[i].n_agents, MPI_INT, message, BUFSIZ, &position, MPI_COMM_WORLD);
            MPI_Pack(all_msg[i].vl_agents, all_msg[i].n_agents, MPI_CHAR, message, BUFSIZ, &position, MPI_COMM_WORLD);
            MPI_Pack(all_msg[i].id_reset, all_msg[i].n_agents, MPI_INT, message, BUFSIZ, &position, MPI_COMM_WORLD);

            MPI_Isend(message, BUFSIZ, MPI_PACKED, i, 0, MPI_COMM_WORLD, &request);
        } else {
            //Internal moves
            if (all_msg[i].n_agents != 0) {
                move(data, all_msg[i], rank);
            }
        }
    }

    free(all_msg[rank].id_agents);
    free(all_msg[rank].vl_agents);
    free(all_msg[rank].id_reset);

    //Messages receive
    for (int i = 0; i < wd_size; i++) {
        if (i != rank) {
            MPI_Recv(message, BUFSIZ, MPI_PACKED, i, 0, MPI_COMM_WORLD, &status);

            position = 0;
            MPI_Unpack(message, BUFSIZ, &position, &all_msg[rank].n_agents, 1, MPI_INT, MPI_COMM_WORLD);
            all_msg[rank].id_agents = malloc(sizeof(int) * n_moves);
            all_msg[rank].vl_agents = malloc(sizeof(char) * n_moves);
            all_msg[rank].id_reset = malloc(sizeof(int) * n_moves);

            if (all_msg[rank].n_agents != 0) {
                MPI_Unpack(message, BUFSIZ, &position, all_msg[rank].id_agents, all_msg[rank].n_agents, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(message, BUFSIZ, &position, all_msg[rank].vl_agents, all_msg[rank].n_agents, MPI_CHAR, MPI_COMM_WORLD);
                MPI_Unpack(message, BUFSIZ, &position, all_msg[rank].id_reset, all_msg[rank].n_agents, MPI_INT, MPI_COMM_WORLD);

                //! Debug //
                /*printf("%d: Numero di spostamenti: %d\n", rank, all_msg[rank].n_agents);
                printf("%d: Indirizzi spostamenti: ", rank);
                printVetInt(all_msg[rank].id_agents, all_msg[rank].n_agents);
                printf("%d: Valori spostati: ", rank);
                printVetChar(all_msg[rank].vl_agents, all_msg[rank].n_agents);
                printf("%d: Valori da resettare: ", rank);
                printVetInt(all_msg[rank].id_reset, all_msg[rank].n_agents);
                printf("\n");*/
                //! ----- //

                move(data, all_msg[rank], rank);
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

    if (wd_size <= ROW) {
        //Matrix generation
        if (rank == MASTER) {
            mat = malloc(ROW * COLUMN * sizeof(char));
            generateMat(mat);
            saveMatrixToFile(mat);
            //sampleMat55(mat);
            //sampleMat1010(mat);
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

        //Start computation
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while (n_itc != 0) {
            //Empty slot calculation
            data.emp_slots = malloc(sizeof(int) * ROW * COLUMN);
            data.n_empty = calcEmptySlots(data, data.emp_slots, rank, wd_size);
            data.emp_slots = realloc(data.emp_slots, sizeof(int) * data.n_empty);

            //Empty slot assignment
            data.my_emp_slots = malloc(sizeof(int) * data.n_empty);
            data.n_my_empty = assignEmptySlots(data, data.my_emp_slots, rank, wd_size);
            data.my_emp_slots = realloc(data.my_emp_slots, sizeof(int) * data.n_my_empty);

            //Identificazione degli agenti da spostare
            int moves[data.n_my_empty];
            int n_moves = findMoves(data, rank, moves);
            if (isOver(n_moves, wd_size, rank)) break;

            //! Debug //
            /*if (rank == 0) {
            printf("Displacement: ");
            printVetInt(data.sec_disp, wd_size);
            printf("Section size: ");
            printVetInt(data.sec_size, wd_size);
            printf("Disp_gt: ");
            printVetInt(data.sec_gt_disp, wd_size);
            printf("Sec_gt: ");
            printVetInt(data.sec_gt_size, wd_size);
            printf("\n");
        }

        syncProcess(rank);

        printf("***********************************************\n");
        printf("Rank: %d\n", rank);
        printf("Numero posizioni vuote assegnate: %d\n", data.n_my_empty);
        printf("Le posizioni vuote assegnate sono: ");
        printVetInt(data.my_emp_slots, data.n_my_empty);
        printf("I miei valori da spostare sono: ");
        printVetInt(moves, n_moves);
        printf("\n");*/
            //! ----- //

            //Spostamento degli agenti
            moveAgents(data, moves, n_moves, rank, wd_size, request, status);

            //! Debug //
            /*printf("Qui processo %d 👨‍🔧, ho elaborato: \n", rank);
        printSubMat(wd_size, data.sub_mat, data.sec_size[rank] / COLUMN, rank);*/
            //! ----- //

            free(data.emp_slots);
            free(data.my_emp_slots);

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
    } else if (rank == MASTER)
        printf("Qui Master 🧑‍🎓, computazione impossibile, numero di slave eccessivo. \n");
}