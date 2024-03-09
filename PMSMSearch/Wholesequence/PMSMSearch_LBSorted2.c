//
// Created by Chris on 11.11.2023.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define INF 1e20       //Pseudo Infitinte number for this code

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define dist(x,y) (fabs(x-y))

#define C_COST 0.5 // cost for merge and split
#define INF 1e20   // pseudo infinite number for this code



typedef struct Index
{   double value;
    int    index;
} Index;

int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    double diff = x->value - y->value;
    if (diff > 0) return 1;
    else if (diff < 0) return -1;
    else return 0;
}

void error(int id)
{
    if(id==1)
    {
        printf("ERROR: Invalid Number of Arguments!");
        printf("Command Usage:   PMSMSearch.exe  data_file  query_file");
        printf("For example  :   PMSMSearch.exe  data.tsv   query.tsv");
    }
    else if ( id == 2 )
        printf("Error while open file!");
    exit(1);
}

double lb_sorted(Index *t, Index *q, int len, double bsf)
{

    double lb, d, tdif, qdif;

    for(int l = 1; l<len; l++){
        d = dist(t[l].value, q[l].value);
        if(d>2*C_COST){
            tdif = dist(t[l].value, t[l-1].value);
            tdif = min(tdif,dist(t[l].value, t[l+1].value));
            qdif = dist(q[l].value, q[l-1].value);
            qdif = min(tdif,dist(q[l].value, q[l+1].value));
            d=2*C_COST+tdif+qdif;
        }
        lb += d;
        if (lb >= bsf)   {
            return lb;
        }
    }

    return lb;
}

double* calculateMsmGreedyArray(double *X, double *Y, int m)
{

    double* greedyArray = malloc((m+1) * sizeof(double));
    for(int i = 0; i < m+1; i++) {
        greedyArray[i] = 0;
    }

    // compute upper Bounds for every diagonal entry
    // the upper bound is computed from right to left: Possibility to update the upper Bound when a diagonal entry is computed
    // upper Bound for the last entry is 0, because there is nothing left to compute
    greedyArray[m] = 0.;

    // assume that the time series are aligned at the end
    double distCurrent = fabs(X[m - 1] - Y[m - 1]);
    double distTmp = distCurrent;
    double xCurrent;
    double yCurrent;
    double xTmp = X[m - 1];
    double yTmp = Y[m - 1];
    int rel = (xTmp > yTmp) ? 1 : 2;

    greedyArray[m - 1] = distCurrent;

    for (int i = 2; i <= m; i++)
    {
        xCurrent = X[m - i];
        yCurrent = Y[m - i];
        distCurrent = xCurrent - yCurrent;

        if ((rel == 1) && distCurrent > 0)
        {
            if (distCurrent > 2 * C_COST && distTmp > 2 * C_COST)
            {
                greedyArray[m - i] = 2 * C_COST + fabs(xCurrent - xTmp) + fabs(yCurrent - yTmp) + greedyArray[m - i + 1];
            }
            else
            {
                greedyArray[m - i] = distCurrent + greedyArray[m - i + 1];
            }
        }
        else if ((rel == 1) && distCurrent <= 0)
        {
            greedyArray[m - i] = -1 * distCurrent + greedyArray[m - i + 1];
            rel = 2;
        }
        else if ((rel == 2) && distCurrent <= 0)
        {
            if (fabs(distCurrent) > 2 * C_COST && fabs(distTmp) > 2 * C_COST)
            {
                greedyArray[m - i] = 2 * C_COST + fabs(xCurrent - xTmp) + fabs(yCurrent - yTmp) + greedyArray[m - i + 1];
            }
            else
            {
                greedyArray[m - i] = -1 * distCurrent + greedyArray[m - i + 1];
            }
        }
        else
        {
            greedyArray[m - i] = distCurrent + greedyArray[m - i + 1];
            rel = 1;
        }

        distTmp = distCurrent;
        xTmp = xCurrent;
        yTmp = yCurrent;
    }

    return greedyArray;
}

unsigned int computeBandwidth(double upperBound)
{

    return (unsigned int)ceil(upperBound / C_COST);
}

/**
 * cost of Split/Merge operation
 *
 * @param new_point point to merge/ split to
 * @param x         xcoord
 * @param y         ycoord
 * @return cost for merge/split
 */
double C(double new_point, double x, double y)
{

    // c - cost of Split/Merge operation. Change this value to what is more
    // appropriate for your data.
    if (new_point < min(x, y) || new_point > max(x, y))
    {
        return C_COST + min(fabs(new_point - x), fabs(new_point - y));
    }

    return C_COST;
}

double getLowerBound(int xCoord, int yCoord)
{
    return fabs(xCoord - yCoord) * C_COST;
}

/**
 * compute the Msm distance table with pruned entries
 *
 * @return pair: first Double: Distance, second Double: relative amount of pruned cells
 * msmDistPruned by Jana Holznigenkemper
 */
double msmDistPruned(double *X, double *Y, int m, double bsf)
{

    double* upperBoundArray = calculateMsmGreedyArray(X, Y, m);

    double* ts1 = malloc((m+2) * sizeof(double));
    double* ts2 = malloc((m+2) * sizeof(double));
    double upperBound = upperBoundArray[0] + 0.0000001;
    ts1[0] = INF;
    ts2[0] = INF;

    for (int i=1; i <=m+1; i++) {
        ts1[i] = X[i-1]; // Copy elements from X to ts1
        ts2[i] = Y[i-1]; // Copy elements from X to ts1
    }

    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules

    int i, j, k;
    double* tmpArray = malloc((m+2) * sizeof(double));
    for(k=0; k<m+2; k++)    tmpArray[k]=INF;

    // value storing the first "real value" of the array before overwriting it
    //  the first value of the first row has to be 0
    double tmp = 0;

    // index for pruning start and end of row
    unsigned int sc = 1;
    unsigned int ec = 1;

    // remember if an entry smaller than UB was found -> cannot cut
    bool smallerFound, smaller_as_bsf;
    int ecNext;

    // initialize first row
    //  the first entry of this row is inf, this inf is used for computing the second row
    //  Arrays.fill(tmpArray, Double.POSITIVE_INFINITY); --> this was done at the initialize step

    //  int counterBandwidth =0;
    // row index
    for (i = 1; i < m+1; i++)
    {

        // compute bandwidth regarding the upper bound
        unsigned int bandwidth = computeBandwidth(upperBound);
        unsigned int start = (bandwidth > i) ? sc : max(sc, i - bandwidth);

        //unsigned int end = min(i + bandwidth + 1, tmpArray.size());
        unsigned int end = min(i + bandwidth + 1, m+1);
        double xi = ts1[i];
        // the index for the pruned end cannot be lower than the diagonal
        // All entries on the diagonal have to be equal or smaller than
        // the upper bound (Euclidean distance = diagonal path)
        ecNext = i;
        smallerFound = false;
        smaller_as_bsf = false;
        // column index
        for (j = start; j < end; j++)
        {

            double yj = ts2[j];
            double d1, d2, d3;
            d1 = tmp + fabs(xi - yj);
            // merge
            d2 = tmpArray[j] + C(xi, ts1[i - 1], yj);
            // split
            d3 = tmpArray[j - 1] + C(yj, xi, ts2[j - 1]);
            // store old entry before overwriting
            tmp = tmpArray[j];

            tmpArray[j] = min(d1, min(d2, d3));

            if (tmpArray[j] < bsf) {
                smaller_as_bsf = true;
            }
            // PruningExperiments strategy
            double lb = getLowerBound(i, j);

            if ((tmpArray[j] + lb) > upperBound)
            {
                if (!smallerFound)
                    sc = j + 1;
                if (j > ec)
                {
                    for(k=j+1; k<m+1; k++)    tmpArray[k]=INF;
                    break;
                }
            }
            else
            {
                smallerFound = true;
                ecNext = j + 1;
            }

            if (i == j)
            {
                upperBound = tmpArray[j] + upperBoundArray[j] + 0.00001;
            }
        }
        if (!smaller_as_bsf) return INF;

        for(k=1; k<sc; k++)    tmpArray[k]=INF;
        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }
    free(upperBoundArray);
    free(ts1);
    free(ts2);

    double result = tmpArray[m];
    free(tmpArray);
    return result;
}


int main(  int argc , char *argv[] )
{
    FILE *sp;
    FILE *qp;
    int m, query_size, sequence_size, i, j, tp = 0,lb_count=0;
    char dataset[50];
    char querypath[200];
    char sequencepath[200];
    double d, t1, t2, bsf, distance, bclass, acc, glb;
    Index *Q_tmp, *T_tmp;

    //read args
    if (argc<=4)
        error(4);
    // Copy dataset name from argv[1]
    strncpy(dataset, argv[1], 99);
    dataset[99] = '\0'; // Ensuring null-termination

    // Construct querypath
    snprintf(querypath, sizeof(querypath), "data/%s/%s_TEST.tsv", dataset, dataset);

    // Construct sequencepath
    snprintf(sequencepath, sizeof(sequencepath), "data/%s/%s_TRAIN.tsv", dataset, dataset);

    m = atol(argv[2]);
    query_size = atol(argv[3]);
    sequence_size = atol(argv[4]);

    double** q_file = (double**)malloc(query_size * sizeof(double*));
    double** s_file = (double**)malloc(sequence_size * sizeof(double*));
    double* qclass = (double*)malloc(query_size * sizeof(double));
    double* sclass = (double*)malloc(sequence_size * sizeof(double));
    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);
    T_tmp = (Index *)malloc(sizeof(Index)*m);
    if( T_tmp == NULL )
        error(1);

    for (i = 0; i < query_size; i++) {
        // Allocate a memory block of size m+1 for each row
        q_file[i] = (double*)malloc((m+1) * sizeof(double));
    }
    for (i = 0; i < sequence_size; i++) {
        // Allocate a memory block of size m+1 for each row
        s_file[i] = (double*)malloc((m+1) * sizeof(double));
    }

    qp = fopen(querypath,"r");
    i= 0;
    j=-1;
    while(fscanf(qp,"%lf",&d) != EOF && i < (m+1)*query_size)
    {
        if(i%(m+1)==0){
            j++;
            i=0;
            qclass[j]=d;
        }
        else {
            q_file[j][i-1] = d;

        }
        i++;
    }

    fclose(qp);
    sp = fopen(sequencepath,"r");
    i=0;
    j=-1;
    while(fscanf(sp,"%lf",&d) != EOF && i < (m+1)*sequence_size)
    {
        if(i%(m+1)==0){
            j++;
            i=0;
            sclass[j]=d;
        }
        else {
            s_file[j][i-1] = d;
        }
        i++;
    }
    fclose(sp);

    tp=0;
    t1 = clock();
    for (int i = 0; i < query_size; i++){
        bsf = INF;
        for (int j = 0; j < sequence_size; j++){


            for( int k = 0; k<m; k++)
            {
                Q_tmp[k].value = q_file[i][k];
                Q_tmp[k].index = k;
                T_tmp[k].value = s_file[j][k];
                T_tmp[k].index = k;
            }
            qsort(Q_tmp, m, sizeof(Index),comp);
            qsort(T_tmp, m, sizeof(Index),comp);

            glb = lb_sorted(T_tmp, Q_tmp, m, bsf);
            if(glb < bsf){
                distance = msmDistPruned(q_file[i], s_file[j], m, bsf);
                if(distance < bsf)
                {
                    bsf = distance;
                    bclass = sclass[j];
                }
            }else{
                lb_count++;
            }

        }
        if(qclass[i] == bclass)   tp++;
    }
    t2 = clock();


    for(i = 0; i < query_size; i++) {
        free(q_file[i]);
    }

    for(i = 0; i < sequence_size; i++) {
        free(s_file[i]);
    }
    free(q_file);
    free(s_file);
    free(qclass);
    free(sclass);

    acc = (double)tp / (double)query_size;
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");
    fprintf(rd,"%s,%s,%f,%f,%i\n", "PMSMSearch with LB_Sorted2",dataset,acc, (t2-t1)/CLOCKS_PER_SEC, lb_count);
    fclose(rd);
    return 0;
}

