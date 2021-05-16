#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

///////////////////////////////////////
// Matrix settings
///////////////////////////////////////
#define ROW 5
#define COLUMN 5
#define O_PERCENTAGE 33
#define X_PERCENTAGE 33
#define SAT_THRESHOLD 33.3
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
    int n_agents;
    int *id_agents;
    char *vl_agents;
} Message;

// Debug functions
void sampleMat(char *mat) {
    mat[0] = 'X';
    mat[1] = 'O';
    mat[2] = 'O';
    mat[3] = 'O';
    mat[4] = 'X';

    mat[5] = 'O';
    mat[6] = 'O';
    mat[7] = 'O';
    mat[8] = 'X';
    mat[9] = 'O';

    mat[10] = 'X';
    mat[11] = ' ';
    mat[12] = 'X';
    mat[13] = 'X';
    mat[14] = ' ';

    mat[15] = 'O';
    mat[16] = 'O';
    mat[17] = 'O';
    mat[18] = 'O';
    mat[19] = 'O';

    mat[20] = 'X';
    mat[21] = 'X';
    mat[22] = 'X';
    mat[23] = 'O';
    mat[24] = 'O';
}
void syncProcess(unsigned int tms) {
    usleep(tms * 1000);
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
        //printf("Real SS: %d\n", data[i].sec_size / ROW);
        data.sec_size[i] += ((i == 0) || (i == wd_size - 1)) ? COLUMN : COLUMN * 2;

        if (i == 0) {
            data.sec_disp[i] = 0;
        } else if (i == wd_size - 1)
            data.sec_disp[i] = (ROW * COLUMN) - data.sec_size[i];
        else {
            data.sec_disp[i] = data.sec_disp[i - 1] + data.sec_size[i - 1] - (COLUMN * 2);
        }
    }

    //* printf("Sezioni generate: ");
    //* printVetInt(sec_size, wd_size);
    //* printf("Displacements generati: ");
    //* printVetInt(displacement, wd_size);
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

    //* printf("Rank: %d\n", rank);
    //* printf("Value: %c\n", data.sub_mat[id]);
    //* printf("Index: %d\n", id);
    //* printf("Rel_Row: %d\n", row);
    //* printf("Row_index: %d\n", row_index);
    //* printf("Col_index: %d\n", col_index);
    //* printf("Neighborhood: %d\n", neigh);
    //* printf("Kind: %d\n", my_kynd);
    //* printf("Perc: %.1f\n\n", perc);

    return perc >= SAT_THRESHOLD;
}

//Empty slots management
int calcEmptySlots(Data data, int *e_slots, int rank, int wd_size) {
    int empty_tot = 0;
    int vet_siz[wd_size];
    int vet_disp[wd_size];

    int *vet_emp = malloc(sizeof(int) * ROW * COLUMN);
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
    free(vet_emp);

    return empty_tot;
}
int assignEmptySlots(Data data, int *m_e_slots, int rank, int wd_size) {
    int section = data.n_empty / wd_size;
    int difference = data.n_empty % wd_size;

    int n_my_slots[wd_size];
    for (int i = 0; i < wd_size; i++) {
        n_my_slots[i] = i < difference ? section + 1 : section;
    }

    m_e_slots = realloc(m_e_slots, sizeof(int) * n_my_slots[rank]);

    int k = data.n_empty - 1 - rank;
    for (int i = 0; i < n_my_slots[rank]; i++) {
        m_e_slots[i] = data.emp_slots[k];
        k -= wd_size;
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
void offlineMove(Data data, Message msg, int *lc_reset, int rank) {
    for (int i = 0; i < msg.n_agents; i++) {
        data.sub_mat[msg.id_agents[i] - data.sec_disp[rank]] = msg.vl_agents[i];
        data.sub_mat[lc_reset[i]] = ' ';
    }
}
void onlineMove(Data data, Message msg, int rank) {
    for (int i = 0; i < msg.n_agents; i++) {
        data.sub_mat[msg.id_agents[i] - data.sec_disp[rank]] = msg.vl_agents[i];
    }
}
void moveAgents(Data data, int *moves, int n_moves, int rank, int wd_size, MPI_Request request, MPI_Status status) {
    Message *all_msg = malloc(sizeof(Message) * wd_size);
    int *lc_moves = malloc(sizeof(int) * n_moves);

    //Messages inizialization
    for (int i = 0; i < wd_size; i++) {
        all_msg[i].n_agents = 0;
        all_msg[i].id_agents = malloc(sizeof(int) * n_moves);
        all_msg[i].vl_agents = malloc(sizeof(char) * n_moves);
    }

    //Messages creation
    for (int i = 0; i < wd_size; i++) {
        for (int j = 0; j < n_moves; j++) {
            if (calcMembership(data, i, data.my_emp_slots[j], wd_size)) {
                all_msg[i].id_agents[all_msg[i].n_agents] = data.my_emp_slots[j];
                all_msg[i].vl_agents[all_msg[i].n_agents] = data.sub_mat[moves[j] - data.sec_disp[rank]];
                if (i == rank) lc_moves[all_msg[i].n_agents] = moves[j] - data.sec_disp[rank];
                all_msg[i].n_agents++;
            }
        }
    }

    //Free extra space
    if (all_msg[rank].n_agents == 0) free(lc_moves);

    //* Debug log
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
                printVetInt(lc_moves, all_msg[i].n_agents);
            } else {
                printf("A me stesso niente.\n");
            }
        }
        printf("\n");
    }*/

    //Messages send
    char message[BUFSIZ];
    int position;
    for (int i = 0; i < wd_size; i++) {
        if (i != rank) {
            //External moves
            position = 0;
            MPI_Pack(&all_msg[i].n_agents, 1, MPI_INT, message, BUFSIZ, &position, MPI_COMM_WORLD);
            MPI_Pack(all_msg[i].id_agents, all_msg[i].n_agents, MPI_INT, message, BUFSIZ, &position, MPI_COMM_WORLD);
            MPI_Pack(all_msg[i].vl_agents, all_msg[i].n_agents, MPI_CHAR, message, BUFSIZ, &position, MPI_COMM_WORLD);

            MPI_Isend(message, BUFSIZ, MPI_PACKED, i, 0, MPI_COMM_WORLD, &request);

            free(all_msg[i].id_agents);
            free(all_msg[i].vl_agents);
        } else {
            //Internal moves
            if (all_msg[i].n_agents != 0) {
                offlineMove(data, all_msg[i], lc_moves, rank);

                free(all_msg[i].id_agents);
                free(all_msg[i].vl_agents);
                free(lc_moves);
            }
        }
    }

    for (int i = 0; i < wd_size; i++) {
        if (i != rank) {
            MPI_Recv(message, BUFSIZ, MPI_PACKED, i, 0, MPI_COMM_WORLD, &status);

            position = 0;
            MPI_Unpack(message, BUFSIZ, &position, &all_msg[rank].n_agents, 1, MPI_INT, MPI_COMM_WORLD);
            all_msg[rank].id_agents = malloc(sizeof(int) * n_moves);
            all_msg[rank].vl_agents = malloc(sizeof(char) * n_moves);

            if (all_msg[rank].n_agents != 0) {
                MPI_Unpack(message, BUFSIZ, &position, all_msg[rank].id_agents, all_msg[rank].n_agents, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(message, BUFSIZ, &position, all_msg[rank].vl_agents, all_msg[rank].n_agents, MPI_INT, MPI_COMM_WORLD);

                //* printf("%d: Numero di spostamenti: %d\n", rank, all_msg[rank].n_agents);
                //* printf("%d: Indirizzi spostamenti: ", rank);
                //* printVetInt(all_msg[rank].id_agents, all_msg[rank].n_agents);
                //* printf("%d: Valori spostati: ", rank);
                //* printVetChar(all_msg[rank].vl_agents, all_msg[rank].n_agents);
                //* printf("\n");

                onlineMove(data, all_msg[rank], rank);

                free(all_msg[rank].id_agents);
                free(all_msg[rank].vl_agents);
            }
        }
    }
}

void main() {
    char *mat, *sub_mat;
    Data data;

    //MPI initialization
    int wd_size, rank;
    double start, end;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wd_size);

    //Matrix generation
    if (rank == MASTER) {
        mat = malloc(ROW * COLUMN * sizeof(char));
        //generateMat(mat);
        sampleMat(mat);
        printf("Qui Master 🧑‍🎓, la matrice generata è: \n");
        printMat(mat);
    }

    //Matrix division
    data.sec_size = malloc(sizeof(int) * wd_size);
    data.sec_disp = malloc(sizeof(int) * wd_size);

    //TODO vanno qui

    calcSizes(wd_size, data);
    data.sub_mat = malloc(sizeof(char) * data.sec_size[rank]);
    MPI_Scatterv(mat, data.sec_size, data.sec_disp, MPI_CHAR, data.sub_mat, data.sec_size[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);
    data.r_start = calcStart(rank, wd_size);
    data.r_finish = calcFinish(data, rank, wd_size);

    printf("start: %d\n", data.r_start);
    printf("finish: %d\n", data.r_finish);

    //Empty slot calculation
    data.emp_slots = malloc(sizeof(int) * ROW * COLUMN);
    data.n_empty = calcEmptySlots(data, data.emp_slots, rank, wd_size);
    data.emp_slots = realloc(data.emp_slots, sizeof(int) * data.n_empty);

    //Empty slot assignment
    data.my_emp_slots = malloc(sizeof(int) * data.n_empty);
    data.n_my_empty = assignEmptySlots(data, data.my_emp_slots, rank, wd_size);

    //Start computation
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //! Debug //
    if (rank == 0) {
        printf("Displacement: ");
        printVetInt(data.sec_disp, wd_size);
        printf("Section size: ");
        printVetInt(data.sec_size, wd_size);
        printf("\n");
    }
    syncProcess(rank);

    printf("***********************************************\n");
    printf("Rank: %d\n", rank);
    printf("Numero posizioni vuote assegnate: %d\n", data.n_my_empty);
    printf("Le posizioni vuote assegnate sono: ");
    printVetInt(data.my_emp_slots, data.n_my_empty);
    //! ----- //

    //Identificazione degli agenti da spostare
    int *moves = malloc(sizeof(int) * data.n_my_empty);
    int n_moves = findMoves(data, rank, moves);

    //! Debug //
    //* printf("I miei valori da spostare sono: ");
    //* printVetInt(moves, n_moves);
    //* printf("\n");
    //! ----- //

    //Spostamento degli agenti
    moveAgents(data, moves, n_moves, rank, wd_size, request, status);

    //! Debug //
    //* printf("Qui processo %d 👨‍🔧, ho ricevuto: \n", rank);
    //* printSubMat(wd_size, data.sub_mat, data.sec_size[rank] / COLUMN, rank);
    //! ----- //

    MPI_Gatherv(data.sub_mat, data.sec_size, MPI_CHAR, mat, data.sec_size, data.sec_disp, MPI_CHAR, MASTER, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    //Finish computation

    MPI_Finalize();

    free(data.sec_size);
    free(data.sec_disp);
    free(data.sub_mat);
    free(data.emp_slots);

    if (rank == MASTER)
        printf("\n\n🕒 Time in ms = %f\n", end - start);
}
//TODO Gestione pochi processi