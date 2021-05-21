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

Per dividere il carico di lavoro tra i vari processi è stato inizialmente valutata la possibilità di dividere la matrice in modo perfettaente equo tra i processi, andando cioè a valutare anche divisioni composte da sezioni di righe. Successivamente si è però notato che tale soluzione creava un notevo dispendio di tempo, rispetto alla soluzione più classica di dividere la matrice in righe, che è la soluzione adottata.

La divisione effettiva viene realizzata tramite l'utilizzo di una **ScatterV**:

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