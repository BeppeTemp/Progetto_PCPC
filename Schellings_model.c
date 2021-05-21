#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mpi.h"

//* Outuput type
#define OUTPUT_TYPE 0  //? 0 HTML output, 1 CLI output. 2 CLI reduced

//*#region Computation settings
#define ROWS 40           //? Number of rows
#define COLUMNS 60        //? Number of columns
#define O_PERCENTAGE 33   //? Percentage of O agents
#define X_PERCENTAGE 33   //? Percentage of X agents
#define SAT_THRESHOLD 35  //? Percentage of satisfaction required
#define N_ITERACTION 100  //? Number of iteration calculated
#define ASSIGN_SEED 10    //? Seed for free location assignment
//#endregion

//*#region Utils and Costants
#define PRINT_GREEN(str) printf("\x1b[32m %c \x1b[0m", str);
#define PRINT_RED(str) printf("\x1b[31m %c \x1b[0m", str);
#define MASTER 0
#define MAX_RAND_VALUE 99
//#endregion

//*#region Definition of structures
typedef struct {
    //Assigned start and finish submatrix
    int r_start;
    int r_finish;

    //Assigned submatrix
    char *sub_mat;

    //Number of assigned empty locations and locations
    int n_my_empty;
    int *my_emp_loc;

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
//#endregion

//!#region Debug functions
void sampleMat(char *i_mat) {
    //Fills the matrix with a constant set of values
    i_mat[0] = 'O';
    i_mat[1] = 'O';
    i_mat[2] = 'X';
    i_mat[3] = ' ';
}
void syncProcess(unsigned int rank) {
    //Inserts delays in computation to synchronize the processes
    usleep(rank * 500);
}
void printVetInt(int vet[], int len) {
    //Print integer array
    if (len > 0) {
        for (int i = 0; i < len; i++) {
            printf("| %d ", vet[i]);
        }
        printf("|\n");
    } else
        printf("\n");
}
//!#endregion

//*#region Matrix generation functions
char randomValue() {
    //Converts a random integer value to an agent type
    int value = rand() % MAX_RAND_VALUE + 1;
    if (value > 0 && value <= O_PERCENTAGE)
        return 'O';
    else if (value > O_PERCENTAGE && value <= (X_PERCENTAGE + O_PERCENTAGE))
        return 'X';
    else
        return ' ';
}
void generateMat(char *i_mat) {
    //Randomly fills the matrix
    srand(time(NULL) + MASTER);
    for (int i = 0; i < (ROWS * COLUMNS); i++)
        i_mat[i] = randomValue();
}
//#endregion

//*#region Functions for printing and saving
void printChar(char x) {
    //Print a single agent
    if (x == 'O')
        PRINT_RED(x);
    if (x == 'X')
        PRINT_GREEN(x);
    if (x == ' ')
        printf("   ");
}
void printMat(char *i_mat) {
    //Print the entire matrix
    printf("\n");
    for (int i = 0; i < (ROWS * COLUMNS); i++) {
        printf("|");
        printChar(i_mat[i]);
        if (((i + 1) % (COLUMNS) == 0) && (i != 0))
            printf("|\n");
        if ((ROWS * COLUMNS) == 1)
            printf("|\n");
        if (COLUMNS == 1 && ((i == 0)))
            printf("|\n");
    }
    printf("\n");
}
void printCharHTML(FILE *fp, char x) {
    if (x == 'O')
        fprintf(fp, "<th class=\"red\">%c</th>\n", x);
    if (x == 'X')
        fprintf(fp, "<th class=\"green\">%c</th>\n", x);
    if (x == ' ')
        fprintf(fp, "<th class=\"null\">N</th>\n");
}
void printResultHTML(char *i_mat, char *r_mat, int n_itc, double t, int wd_size) {
    //Print the entire matrix on file
    FILE *fp;
    time_t curtime;
    time(&curtime);
    fp = fopen("result.html", "w+");
    fprintf(fp,
            "<html>\n<head>\n<style>\n"
            "body {\n"
            "text-align: center;\n"
            "font-family: monospace;\n"
            "font-size: 18px;\n"
            "color: rgb(190, 190, 190);\n"
            "background-color: #1c1f24;\n"
            "}\n"
            "table {\n"
            "border-style: solid;\n"
            "border-color: rgb(83, 83, 83);\n"
            "border-width: 2px;\n"
            "margin-left: auto;\n"
            "margin-right: auto;\n"
            "}\n"
            "th {\n"
            "font-size: 8px;\n"
            "width: max-content;\n"
            "border-collapse: collapse;\n"
            "padding-left: 3px;\n"
            "padding-right: 3px;\n"
            "}\n"
            ".green { color: rgb(0, 255, 0); }\n"
            ".red { color: red; }\n"
            ".null{ color:rgba(255, 255, 255, 0)}\n"
            "</style>\n</head>\n<body>\n"
            "<h1>Schelling's model of segregation</h1>\n");
    fprintf(fp, "<h3>📆 %s</h3>\n", ctime(&curtime));
    fprintf(fp, "<h3>⚙ Computation settings:</h3>\n");
    fprintf(fp, "📏 Rows: <b>%d.</b><br>\n", ROWS);
    fprintf(fp, "📐 Collumns: <b>%d.</b><br>\n", COLUMNS);
    fprintf(fp, "<p></p>\n");
    fprintf(fp, "❎ Population: <b>%d%%.</b><br>\n", X_PERCENTAGE);
    fprintf(fp, "⭕ Population: <b>%d%%.</b><br>\n", O_PERCENTAGE);
    fprintf(fp, "🤷‍♂️ Empty slots: <b>%d%%.</b><br>\n", 100 - X_PERCENTAGE - O_PERCENTAGE);
    fprintf(fp, "👨‍👦‍👦 Number of processes: <b>%d.</b><br>\n", wd_size);
    fprintf(fp, "<p></p>\n");
    fprintf(fp, "😎 Satisfaction threshold: <b>%d%%.</b><br>\n", SAT_THRESHOLD);
    fprintf(fp, "📟 Empty slot assigment seed: <b>%d.</b><br>\n", ASSIGN_SEED);
    fprintf(fp, "🚀 Max iterations: <b>%d.</b> <br>\n", N_ITERACTION);
    fprintf(fp, "<h3>📈📉 Computations result:</h3>\n");
    fprintf(fp, "🔬 Number of iterations: <b>%d.</b> <br>\n", n_itc);
    fprintf(fp, "⏲ Time: <b>%fs.</b><br>\n", t);
    fprintf(fp,
            "<h3>Here the master 🧑‍🎓! The resulting matrix is:</h3>\n"
            "<table>\n<tr>\n");
    for (int i = 0; i < (ROWS * COLUMNS); i++) {
        printCharHTML(fp, r_mat[i]);
        if (i == (ROWS * COLUMNS) - 1)
            fprintf(fp, "</tr>\n</table>\n");
        else if (((i + 1) % (COLUMNS) == 0) && (i != 0))
            fprintf(fp, "</tr>\n<tr>");
        else if (COLUMNS == 1 && ((i == 0)))
            fprintf(fp, "</tr>\n<tr>");
    }
    fprintf(fp,
            "<h3>And this is the initial matrix, bye-bye 👋❤️:</h3>\n"
            "<table>\n<tr>\n");
    for (int i = 0; i < (ROWS * COLUMNS); i++) {
        printCharHTML(fp, i_mat[i]);
        if (i == (ROWS * COLUMNS) - 1)
            fprintf(fp, "</tr>\n</table>\n");
        else if (((i + 1) % (COLUMNS) == 0) && (i != 0))
            fprintf(fp, "</tr>\n<tr>");
        else if (COLUMNS == 1 && ((i == 0)))
            fprintf(fp, "</tr>\n<tr>");
    }
    if (ROWS * COLUMNS >= 20000) fprintf(fp, "<h3>Why are you so evil? 😭</h3>\n");
    fprintf(fp, "</body>\n</html>");
    fclose(fp);
}
//#endregion

//*#region Functions for dividing the matrix between processes
void calcSizes(int wd_size, Data data) {
    //Calculate section sizes and displacements for gather and scatter operations
    int section = ROWS / (wd_size);
    int difference = ROWS % (wd_size);

    if (wd_size >= 2)
        for (int i = 0; i < wd_size; i++) {
            data.sec_size[i] = i < difference ? (section + 1) * COLUMNS : section * COLUMNS;
            data.sec_gt_size[i] = data.sec_size[i];
            data.sec_gt_disp[i] = i == 0 ? 0 : data.sec_gt_disp[i - 1] + data.sec_gt_size[i - 1];
            data.sec_size[i] += ((i == 0) || (i == wd_size - 1)) ? COLUMNS : COLUMNS * 2;

            if (i == 0) {
                data.sec_disp[i] = 0;
            } else if (i == wd_size - 1)
                data.sec_disp[i] = (ROWS * COLUMNS) - data.sec_size[i];
            else {
                data.sec_disp[i] = data.sec_disp[i - 1] + data.sec_size[i - 1] - (COLUMNS * 2);
            }
        }
    else {
        data.sec_size[0] = ROWS * COLUMNS;
        data.sec_disp[0] = 0;
    }
}
int calcStart(int rank, int wd_size) {
    //Calculate start position of real submatrix
    if (rank == 0)
        return 0;
    else
        return COLUMNS;
}
int calcFinish(Data data, int rank, int wd_size) {
    //Calculate finish position of real submatrix
    if (rank == wd_size - 1)
        return data.sec_size[rank] - 1;
    else
        return data.sec_size[rank] - COLUMNS - 1;
}
//#endregion

//*#region Functions for calculating the degree of satisfaction
void convertIndex(int id, int *rowS_index, int *col_index) {
    //Converts the index of a vector to its matrix counterpart
    *rowS_index = id / COLUMNS;
    *col_index = id % COLUMNS;
}
int isMyKind(int my_index, int x_index, Data data) {
    //Check if two agents are of the same type
    return data.sub_mat[my_index] == data.sub_mat[x_index];
}
void satCorner(int id, Data data, int *neigh, int *my_kynd, int pos) {
    //Check the satisfaction of an agent in the corner
    *neigh += 3;
    switch (pos) {
        //Top left corner
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS), data);
            *my_kynd += isMyKind(id, (id + COLUMNS) + 1, data);
            break;
        //Top right corner
        case 1:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS), data);
            *my_kynd += isMyKind(id, (id + COLUMNS) - 1, data);
            break;
        //Bottom left corner
        case 2:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS), data);
            *my_kynd += isMyKind(id, (id - COLUMNS) + 1, data);
            break;
        //Bottom right corner
        case 3:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS), data);
            *my_kynd += isMyKind(id, (id - COLUMNS) - 1, data);
            break;
    }
}
void satEdge(int id, Data data, int *neigh, int *my_kynd, int pos) {
    //Check the satisfaction of an agent in the edge
    *neigh += 5;
    switch (pos) {
        //Top edge
        case 0:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS), data);
            *my_kynd += isMyKind(id, (id + COLUMNS) + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS) - 1, data);
            break;
        //Left edge
        case 1:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS), data);
            *my_kynd += isMyKind(id, (id - COLUMNS), data);
            *my_kynd += isMyKind(id, (id + COLUMNS) + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS) + 1, data);
            break;
        //Right edge
        case 2:
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id + COLUMNS), data);
            *my_kynd += isMyKind(id, (id - COLUMNS), data);
            *my_kynd += isMyKind(id, (id + COLUMNS) - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS) - 1, data);
            break;
        //Bottom edge
        case 3:
            *my_kynd += isMyKind(id, id + 1, data);
            *my_kynd += isMyKind(id, id - 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS), data);
            *my_kynd += isMyKind(id, (id - COLUMNS) + 1, data);
            *my_kynd += isMyKind(id, (id - COLUMNS) - 1, data);
            break;
    }
}
void satCenter(int id, Data data, int *neigh, int *my_kynd) {
    //Check the satisfaction of a central agent
    *neigh += 8;
    *my_kynd += isMyKind(id, id + 1, data);
    *my_kynd += isMyKind(id, id - 1, data);
    *my_kynd += isMyKind(id, (id - COLUMNS), data);
    *my_kynd += isMyKind(id, (id + COLUMNS), data);
    *my_kynd += isMyKind(id, (id - COLUMNS) + 1, data);
    *my_kynd += isMyKind(id, (id - COLUMNS) - 1, data);
    *my_kynd += isMyKind(id, (id + COLUMNS) + 1, data);
    *my_kynd += isMyKind(id, (id + COLUMNS) - 1, data);
}
int calcSat(int id, Data data, int rank) {
    //Calculate the satisfaction of an agent
    int neigh = 0, my_kynd = 0;
    int rowS = data.sec_size[rank] / COLUMNS;
    int rowS_index, col_index;

    convertIndex(id, &rowS_index, &col_index);

    if (rowS_index == 0)
        if (col_index == 0)
            satCorner(id, data, &neigh, &my_kynd, 0);  //Upper left corner
        else if (col_index == COLUMNS - 1)
            satCorner(id, data, &neigh, &my_kynd, 1);  //Upper right corner
        else
            satEdge(id, data, &neigh, &my_kynd, 0);  //Upper edge

    else if (rowS_index == rowS - 1)
        if (col_index == 0)
            satCorner(id, data, &neigh, &my_kynd, 2);  //Lower left corner
        else if (col_index == COLUMNS - 1)
            satCorner(id, data, &neigh, &my_kynd, 3);  //Lower right corner
        else
            satEdge(id, data, &neigh, &my_kynd, 3);  //Lower edge

    else if (col_index == 0 && rowS_index > 0 && rowS_index < rowS - 1)
        satEdge(id, data, &neigh, &my_kynd, 1);  //Left edge

    else if (col_index == COLUMNS - 1 && rowS_index > 0 && rowS_index < rowS)
        satEdge(id, data, &neigh, &my_kynd, 2);  //Right edge
    else
        satCenter(id, data, &neigh, &my_kynd);  //Center

    float perc = (100 / (float)neigh) * my_kynd;
    return perc >= SAT_THRESHOLD;
}
//#endregion

//*#region Free locations management functions
void swap(int *a, int *b) {
    //Swap two value
    int temp = *a;
    *a = *b;
    *b = temp;
}
void shuffle(int *vet, int length) {
    //Randomizes the position of values in a vector
    for (int i = length - 1; i > 0; i--) {
        // Pick a random index from 0 to i
        int j = rand() % (i + 1);

        // Swap arr[i] with the element at random index
        swap(&vet[i], &vet[j]);
    }
}
int calcEmptySlots(Data data, int *my_emp_loc, int rank, int wd_size, int n_itc) {
    //Locate all empty locations in submatrices and associates them
    int empty_tot = 0;
    int vet_siz[wd_size];
    int vet_disp[wd_size];

    // Find free locations in the local sub-matrix
    int *vet_emp = malloc(sizeof(int) * (data.r_finish - data.r_start));
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
    int *emp_slots = malloc(sizeof(int) * empty_tot);
    MPI_Allgatherv(vet_emp, n_my_emp, MPI_INT, emp_slots, vet_siz, vet_disp, MPI_INT, MPI_COMM_WORLD);

    //Randomizes the array of free locations
    srand(ASSIGN_SEED + n_itc);
    shuffle(emp_slots, empty_tot);

    //Assignment of free locations to processes
    data.n_my_empty = empty_tot / wd_size;
    int k = rank * data.n_my_empty;
    for (int i = 0; i < data.n_my_empty; i++) {
        my_emp_loc[i] = emp_slots[k];
        k++;
    }

    free(emp_slots);
    free(vet_emp);
    return data.n_my_empty;
}
//#endregion

//*#region Functions for moving agents
int findMoves(Data data, Move *my_moves, int rank, int wd_size) {
    //Locates the agents to be moved in the submatrix
    int k = 0, is_over = 0;
    int n_moves = data.n_my_empty;

    for (int i = data.r_start; i <= data.r_finish; i++) {
        if (data.sub_mat[i] != ' ')
            if (!calcSat(i, data, rank)) {
                if (n_moves == 0) break;
                my_moves[k].id_agent = data.my_emp_loc[k];
                my_moves[k].id_reset = data.sec_disp[rank] + i;
                my_moves[k].vl_agent = data.sub_mat[i];
                is_over++;
                n_moves--;
                k++;
            }
    }

    for (int i = 0; i < n_moves; i++) {
        my_moves[k].id_agent = -1;
        my_moves[k].id_reset = -1;
        my_moves[k].vl_agent = 'n';
        k++;
    }

    //Check if all processes have finished their movements
    int tot_over;
    MPI_Allreduce(&is_over, &tot_over, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return (tot_over == 0);
}
int calcMembership(Data data, int rank, int index) {
    //Given an index and a process checks if it is associated with it
    int start = data.sec_disp[rank];
    int finish = data.sec_disp[rank] + data.sec_size[rank];

    return (index >= start && index < finish);
}
void move(Data data, Move *my_moves, MPI_Datatype move_data_type, int wd_size, int rank) {
    //Carries out the movement operations of the agents
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
//#endregion

//*#region Results management
void gatherResult(Data data, int rank, char *i_mat) {
    //Gather results and save them to a file
    char *section = malloc(sizeof(char) * data.sec_gt_size[rank]);
    int k = 0;
    for (int i = data.r_start; i <= data.r_finish; i++) {
        section[k] = data.sub_mat[i];
        k++;
    }
    MPI_Gatherv(section, data.sec_gt_size[rank], MPI_CHAR, i_mat, data.sec_gt_size, data.sec_gt_disp, MPI_CHAR, MASTER, MPI_COMM_WORLD);

    free(section);
}
//#endregion

//* Main funtion
void main() {
    //Variable definitions
    int n_itc = N_ITERACTION;
    char *i_mat, *r_mat;
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

    if (wd_size <= ROWS) {
        //Matrix generation
        if (rank == MASTER) {
            printf("Computation \033[0;33mstarted...\x1b[0m \n");
            i_mat = malloc(ROWS * COLUMNS * sizeof(char));
            generateMat(i_mat);
            //sampleMat(i_mat);
        }

        //Matrix division
        data.sec_size = malloc(sizeof(int) * wd_size);
        data.sec_disp = malloc(sizeof(int) * wd_size);
        data.sec_gt_size = malloc(sizeof(int) * wd_size);
        data.sec_gt_disp = malloc(sizeof(int) * wd_size);
        calcSizes(wd_size, data);
        data.sub_mat = malloc(sizeof(char) * data.sec_size[rank]);
        MPI_Scatterv(i_mat, data.sec_size, data.sec_disp, MPI_CHAR, data.sub_mat, data.sec_size[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);

        //Start and finish calculation
        data.r_start = calcStart(rank, wd_size);
        data.r_finish = calcFinish(data, rank, wd_size);

        data.my_emp_loc = malloc(sizeof(int) * ROWS * COLUMNS);

        //Start computation
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        while (n_itc != 0) {
            //Empty slot calculation
            data.n_my_empty = calcEmptySlots(data, data.my_emp_loc, rank, wd_size, n_itc);

            if (n_itc == N_ITERACTION) {
                //End of computation test
                data.my_emp_loc = realloc(data.my_emp_loc, sizeof(int) * data.n_my_empty);
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
        if (rank == 0) r_mat = malloc(sizeof(char) * ROWS * COLUMNS);
        if (wd_size >= 2)
            gatherResult(data, rank, r_mat);
        else
            r_mat = data.sub_mat;

        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        //Finish computation
    }

    MPI_Finalize();

    //Results display
    if (wd_size <= ROWS) {
        free(data.sec_size);
        free(data.sec_disp);
        free(data.sec_gt_size);
        free(data.sec_gt_disp);
        free(data.sub_mat);
        free(data.my_emp_loc);

        if (rank == MASTER) {
            printf("Computation \x1b[32mfinished.\x1b[0m \n\n");
            if (OUTPUT_TYPE == 0) {
                printResultHTML(i_mat, r_mat, N_ITERACTION - n_itc, end - start, wd_size);
            } else {
                printf("🔬 Number of iterations: %d.\n", N_ITERACTION - n_itc);
                printf("⏲  Time: %fs.\n\n", end - start);

                if (OUTPUT_TYPE == 1) {
                    printf("Here the master 🧑‍🎓! The resulting matrix is: \n");
                    printMat(r_mat);
                    printf("And this is the initial matrix, bye-bye 👋: \n");
                    printMat(i_mat);
                }
                if (ROWS * COLUMNS >= 20000) printf("Why are you so evil? 😭\n");
            }
            free(i_mat);
            free(r_mat);
        }
    } else if (rank == MASTER)
        printf("Here the master 🧑‍🎓! \x1b[31mcomputation impossible\x1b[0m, excessive number of slaves. \n");
}