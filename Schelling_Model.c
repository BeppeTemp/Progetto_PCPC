#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

//////////////////////////////////
// Impostazioni della Matrice
//////////////////////////////////
#define SIZE 20              //Dimensione della matrice
#define O_PERCENTAGE 33     //Percentuale di agenti di tipo O
#define X_PERCENTAGE 33     //Percentuale di agenti di tipo X
#define SAT_THRESHOLD 33.3  //Soglia di soddisfazione in percentuale
#define N_ITERACTION 10      //Numero di interazioni massime consentite
#define ASSIGN_SEED 117     //Seme per l'assegnazione degli slot vuoti
//////////////////////////////////
// Costanti
//////////////////////////////////
#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\033[0;33m"
#define COLOR_RESET "\x1b[0m"
#define MASTER 0
#define MAX_RAND_VALUE 99
///////////////////////////////////

//Definizione delle strutture
typedef struct {
    //Valori relativi di inizio e fine sottomatrice
    int r_start;
    int r_finish;

    //Sottomatrice assegnata al processo
    char *sub_mat;

    //Numero di posizioni vuote assegnate
    int n_my_empty;
    int *my_emp_slots;

    //Sezioni per la gather di divisione del workload
    int *sec_size;
    int *sec_disp;

    //Sezione per la gather di chiusura
    int *sec_gt_size;
    int *sec_gt_disp;
} Data;
typedef struct {
    //Indice di arrivo dell'agente in spostamento
    int id_agent;

    //Indice di partenza dell'agente in spostamento
    int id_reset;

    //Tipologia di agente in spostamente
    char vl_agent;
} Move;

//? Funzioni per il debug ?//
void sampleMat(char *mat) {
    mat[0] = 'X';
    mat[1] = 'O';
    mat[2] = 'O';
    mat[3] = 'O';
    mat[4] = 'O';

    mat[5] = 'O';
    mat[6] = 'O';
    mat[7] = ' ';
    mat[8] = 'X';
    mat[9] = 'X';

    mat[10] = ' ';
    mat[11] = 'O';
    mat[12] = ' ';
    mat[13] = 'O';
    mat[14] = 'O';

    mat[15] = 'O';
    mat[16] = 'O';
    mat[17] = 'O';
    mat[18] = 'O';
    mat[19] = ' ';

    mat[20] = 'O';
    mat[21] = 'X';
    mat[22] = 'X';
    mat[23] = 'X';
    mat[24] = 'O';
}
void syncProcess(unsigned int rank) {
    usleep(rank * 500);
}
/* Funzioni inutili
void saveMatrixToFile(char *mat) {
    FILE *file = fopen("Matrix.txt", "w");
    int results = fputs(mat, file);
    fclose(file);
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
            if (((i + 1) % (SIZE) == 0) && (i != 0))
                printf("|\n");
        }
    else
        for (int i = start; i < len; i++) {
            printf("|");
            printChar(vet[i]);
            if (((i + 1) % (SIZE) == 0) && (i != 0))
                printf("|\n");
        }
}
void printSubMat(int wd_size, char *mat, int row, int rank) {
    printf("\n");
    if (rank == 0) {
        printRows(mat, 0, (row - 1) * SIZE, 0);
        printRows(mat, (row * SIZE) - SIZE, row * SIZE, 1);
        printf("\n");
    } else if (rank == wd_size - 1) {
        printRows(mat, 0, SIZE, 1);
        printRows(mat, SIZE, row * SIZE, 0);
        printf("\n");
    } else {
        printRows(mat, 0, SIZE, 1);
        printRows(mat, SIZE, (row * SIZE) - SIZE, 0);
        printRows(mat, (row * SIZE) - SIZE, row * SIZE, 1);
        printf("\n");
    }
}*/
//?-----------------?//

//Funzioni di generazione della matrice
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
    for (int i = 0; i < (SIZE * SIZE); i++)
        mat[i] = randomValue();
}

//Funzioni di stampa
void printChar(char x) {
    if (x == 'O')
        printf(COLOR_GREEN " %c " COLOR_RESET, x);
    if (x == 'X')
        printf(COLOR_RED " %c " COLOR_RESET, x);
    if (x == ' ')
        printf("   ");
}
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

//Funzioni per la divisione della matrice
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
int calcStart(int rank, int wd_size) {
    if (rank == 0)
        return 0;
    else
        return SIZE;
}
int calcFinish(Data data, int rank, int wd_size) {
    if (rank == 0)
        return data.sec_size[rank] - SIZE - 1;
    else if (rank == wd_size - 1)
        return data.sec_size[rank] - 1;
    else
        return data.sec_size[rank] - SIZE - 1;
}

//Funzioni per il calcolo del grado di soddisfazione
void convertIndex(int id, int *row_index, int *col_index) {
    *row_index = id / SIZE;
    *col_index = id % SIZE;
}
int isMyKind(int my_index, int x_index, Data data) {
    return data.sub_mat[my_index] == data.sub_mat[x_index];
}
void satCorner(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 3;
    switch (pos) {
        //Angolo superiore sinistro
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            break;
        //Angolo superiore destro
        case 1:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            break;
        //Angolo inferiore sinistro
        case 2:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            break;
        //Angolo inferiore destro
        case 3:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
    }
}
void satEdge(int id, Data data, int *neigh, int *my_kynd, int pos) {
    *neigh += 5;
    switch (pos) {
        //Bordo superiore
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            break;
        //Bordo sinistro
        case 1:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            break;
        //Bordo destro
        case 2:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id + SIZE) - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
        //Bordo inferiore
        case 3:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - SIZE), data);
            *my_kynd += isMyKind(id, (id - SIZE) + 1, data);
            *my_kynd += isMyKind(id, (id - SIZE) - 1, data);
            break;
    }
}
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

//Funzioni di gestione degli slot liberi
void shuffle(int *vet, int length) {
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

    //Individuazione degli slot liberi nella sottomatrice
    int vet_emp[data.r_finish - data.r_start];
    int n_my_emp = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] == ' ') {
            vet_emp[n_my_emp] = data.sec_disp[rank] + i;
            n_my_emp++;
        }
    }

    //Aggregazione degli slot liberi globali
    MPI_Allgather(&n_my_emp, 1, MPI_INT, vet_siz, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < wd_size; i++) {
        empty_tot += vet_siz[i];
        vet_disp[i] = i == 0 ? 0 : vet_disp[i - 1] + vet_siz[i - 1];
    }
    int emp_slots[empty_tot];
    MPI_Allgatherv(vet_emp, n_my_emp, MPI_INT, emp_slots, vet_siz, vet_disp, MPI_INT, MPI_COMM_WORLD);

    //Randomizzazione dell'array degli slot liberi
    srand(ASSIGN_SEED);
    for (int i = 0; i < rand() % 10; i++) {
        shuffle(emp_slots, empty_tot);
    }

    //Assegnazione degli slot liberi ai processi
    data.n_my_empty = empty_tot / wd_size;
    int k = rank * data.n_my_empty;
    for (int i = 0; i < data.n_my_empty; i++) {
        my_emp_slots[i] = emp_slots[k];
        k++;
    }

    return data.n_my_empty;
}

//Funzioni per lo spostamento degli agenti
int findMoves(Data data, Move *my_moves, int rank, int wd_size) {
    int k = 0, is_over = 0;
    int n_moves = data.n_my_empty;

    //Identificazione agenti da spostare
    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] != ' ')
            if (!calcSat(i, data, rank)) {
                my_moves[k].id_agent = data.my_emp_slots[k];
                my_moves[k].id_reset = data.sec_disp[rank] + i;
                my_moves[k].vl_agent = data.sub_mat[i];
                is_over++;
                if (n_moves-- == 0) break;
                k++;
            }
    }

    for (int i = 0; i < n_moves; i++) {
        my_moves[k].id_agent = -1;
        my_moves[k].id_reset = -1;
        my_moves[k].vl_agent = 'n';
        k++;
    }

    int tot_over[wd_size];
    MPI_Allgather(&is_over, 1, MPI_INT, tot_over, 1, MPI_INT, MPI_COMM_WORLD);
    is_over = 0;
    for (int i = 0; i < wd_size; i++) is_over += tot_over[i];

    return (is_over == 0);
}
int calcMembership(Data data, int rank, int index) {
    int start = data.sec_disp[rank];
    int finish = data.sec_disp[rank] + data.sec_size[rank];

    return (index >= start && index < finish);
}
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
        printf("%d: id %d, vl %c, rs %d\n", rank, moves[i].id_agent, moves[i].vl_agent, moves[i].id_reset);
        if (moves[i].id_agent != -1) {
            if (calcMembership(data, rank, moves[i].id_agent)) {
                data.sub_mat[moves[i].id_agent - data.sec_disp[rank]] = moves[i].vl_agent;
            }
            if (calcMembership(data, rank, moves[i].id_reset)) {
                data.sub_mat[moves[i].id_reset - data.sec_disp[rank]] = ' ';
            }
        }
    }
    printf("\n");
    free(moves);
}

//Stop operation
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
    //Definizione delle variabili
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
            generateMat(mat);
            //sampleMat(mat);
            printf("Qui Master 🧑‍🎓, la matrice generata è: \n");
            //printMat(mat);
        }

        printf("qui ci sto");

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

        data.my_emp_slots = malloc(sizeof(int) * SIZE * SIZE);

        //Start computation
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        

        while (n_itc != 0) {
            //Empty slot calculation
            data.n_my_empty = calcEmptySlots(data, data.my_emp_slots, rank, wd_size);

            if (n_itc == N_ITERACTION) {
                data.my_emp_slots = realloc(data.my_emp_slots, sizeof(int) * data.n_my_empty);
            }

            

            //Identificazione degli agenti da spostare
            Move *my_moves = malloc(sizeof(Move) * data.n_my_empty);
            if (findMoves(data, my_moves, rank, wd_size)) break;

            //Move agent
            move(data, my_moves, move_data_type, wd_size, rank);

            free(my_moves);
            n_itc--;
        }

        //Creazione della matrice dei risultati
        gatherResult(data, rank, mat);

        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        //Finish computation
    }

    MPI_Finalize();

    //Visualizzazione risultati
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