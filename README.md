# **[Schelling's model of segregation](https://en.wikipedia.org/wiki/Schelling's_model_of_segregation)**

Progetto di corso per l'esame di **Programmazione Concorrente e Parallela su Cloud 2020/21**.

- Studente: **Giuseppe Arienzo**
- Matricola: **0522501062**
- Data di consegna: **01/06/2021**

<img src="docs/Images/sample.png"/>

## **Descrizione dell'implementazione**

La seguente implementazione si basa su **quattro passi di computazione** ben definiti

- [Divisione della matrice tra i processi coinvolti](#divisione-della-matrice-tra-i-processi-coinvolti)
- [Calcolo della soddisfazione degli agenti](#calcolo-della-soddisfazione-degli-agenti)
- [Spostamento degli agenti insoddisfatti nella matrice](#spostamento-degli-agenti-insoddisfatti-nella-matrice)
- [Aggregazione dei risultati e presentazione](#aggregazione-dei-risultati-e-presentazione)

### **Divisione della matrice tra i processi coinvolti**

Per dividere il carico di lavoro tra i vari processi è stato inizialmente valutata la possibilità di dividere la matrice in modo perfettaente equo tra questi ultimi, andando cioè a valutare divisioni composte da sotto sezioni di righe. Successivamente si è però notato che tale soluzione creava un notevo dispendio di tempo di computazione, che rendeva vano ogni miglioramento dovuto all distribuzione più equa deglia genti, si è per questo motivo optato per una soluzione più classica dividendo la matrice per righe, in questo modo nel caso pessimo uno o più agenti gestiscono **COLLUMS** agenti in più rispetto agli altri.

La divisione effettiva viene infine realizzata tramite l'utilizzo di una **ScatterV**:

```c  
MPI_Scatterv(i_mat, data.sec_size, data.sec_disp, MPI_CHAR, data.sub_mat, data.sec_size[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);
```

ed il codice che si occupa del calcolo delle grandezze necessarie alla divisione è il seguente:

```c
void calcSizes(int wd_size, Data data) {
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
```

é importante notare come la seguente funzione si occupa anche del calcolo delle grandezze necessarie alle successive operationi di gather con lo scopo riunire tutte le sotto matrici e presentare il risultato finale. Tutti questi valori vengono inoltre conservati in una struttura dati **Data** che contiene tutti i dati necessari alla computazione per ogni processo, in particolare per quanto riguarda i dati relativi ai displacement e le section size, tutti i processi dispongono di tutti i dati, questo per evitare inutili e costose operazioni di comunicazione:

```c
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
```

### **Calcolo della soddisfazione degli agenti**

### **Spostamento degli agenti insoddisfatti nella matrice**

### **Aggregazione dei risultati e presentazione**

## **Note sull'implementazione**

### **Compilazione**

### **Esecuzione**

## **Risultati**

I={some value}
P={some value}

### K = {some value}

#### Scalabilità debole N={some value}

#### Scalabilità forte N={some value}

### K/2 = {some value}

#### Scalabilità debole N={some value}

#### Scalabilità forte N={some value} 

### 2K = {some value}

#### Scalabilità debole N={some value}

#### Scalabilità forte N={some value} 

## **Descrizione dei risultati**

## **Conclusioni**