#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) (fabs(x-y))
#define INF 1e20       //Pseudo Infitinte number for this code
#define C_COST 0.5 // cost for merge and split

int mergesplit_counter = 0; // global merge and split counter
int move_counter = 0; // global move counter
int mil_mergesplit_counter = 0; // global merge and split counter
int mil_move_counter = 0; // global move counter

/// Data structure for sorting the query
typedef struct Index
{   double value;
    int    index;
} Index;

/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    double diff = x->value - y->value;
    if (diff > 0) return 1;
    else if (diff < 0) return -1;
    else return 0;
}


double lb_sorted(Index *t, Index *q, int len, double bsf)
{
    double lb = 0, d;

    for(int l = 1; l<len; l++){
        d = dist(t[l].value, q[l].value);
        if(d>C_COST)d=C_COST;
        lb += d;
        if (lb >= bsf)   {
            return lb;
        }
    }
    return lb;
}

double calculateMsmGreedyArray(double *X, double *Y, int m, double *greedyArray)
{
    int i, k;
    for(k=0; k<m+1; k++)    greedyArray[k]=0;
    greedyArray[m] = 0.;
    double distCurrent = fabs(X[m - 1] - Y[m - 1]);
    double distTmp = distCurrent;
    double xCurrent;
    double yCurrent;
    double xTmp = X[m - 1];
    double yTmp = Y[m - 1];
    int rel = (xTmp > yTmp) ? 1 : 2;
    greedyArray[m - 1] = distCurrent;
    for (i = 2; i <= m; i++)
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
    return *greedyArray;
}

unsigned int computeBandwidth(double upperBound)
{
    return (unsigned int)ceil(upperBound / C_COST);
}

double C(double new_point, double x, double y)
{
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

double msmDistPruned(double *X, double *Y, int m, double bsf, double sakoe_bandwidth, double *tmpArray, double *upperBoundArray)
{
    *upperBoundArray = calculateMsmGreedyArray(X, Y, m, upperBoundArray);
    double upperBound = upperBoundArray[0] + 0.0000001;
    int i, j, k;
    for(k=0; k<m+1; k++)    tmpArray[k]=INF;
    double tmp = 0;
    unsigned int sc = 1;
    unsigned int ec = 1;
    bool smallerFound, smaller_as_bsf;
    int ecNext;
    for (i = 0; i < m+1; i++)
    {
        int start = max(1, i - sakoe_bandwidth);
        int end = min(m, i + sakoe_bandwidth);
        double xi = X[i];
        ecNext = i;
        smallerFound = false;
        smaller_as_bsf = false;
        for (j = start; j <= end; j++)
        {
            double yj = Y[j];
            double d1, d2, d3;
            d1 = tmp + fabs(xi - yj);
            // merge
            d2 = tmpArray[j] + C(xi, X[i - 1], yj);
            // split
            d3 = tmpArray[j - 1] + C(yj, xi, Y[j - 1]);
            // store old entry before overwriting
            tmp = tmpArray[j];
            if (d1 <= min(d2, d3)){
                move_counter++; // move
                tmpArray[j] = d1;
            }else{
                mergesplit_counter++; // merge or split
                tmpArray[j] = min(d2, d3);
            }
            if (move_counter == 1000000){
                mil_move_counter++;
                move_counter=0;
            }
            if (mergesplit_counter == 1000000){
                mil_mergesplit_counter++;
                mergesplit_counter=0;
            }
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
        tmp = INF;
        ec = ecNext;
    }
    return tmpArray[m];
}

void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
    }
    exit(1);
}

int main(  int argc , char *argv[] )
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;
    double d;
    double ex , ex2 , mean, std;
    double distCalc=0, global_lb=0, bandwidth;
    double *buffer, *u_buff, *l_buff;
    double t1,t2,t3;
    double *time_result;
    double *tmpArray, *upperBoundArray;
    int tr_count = 0;
    long long loc = 0;
    long long i , j;
    int m=-1, r=-1;
    int lbsorted = 0;
    Index *Q_tmp, *T_tmp;


    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;
    /// If not enough input, display an error.
    if (argc<=4)
        error(4);
    /// read args
    if (argc>3)
        m = atol(argv[3]);
    /// read Sakoe
    if (argc>4)
    {
        bandwidth = atol(argv[4]);
        bandwidth = ((double)bandwidth/100.0)*m;
    }
    fp = fopen(argv[1],"r");
    if( fp == NULL )
        error(2);
    qp = fopen(argv[2],"r");
    if( qp == NULL )
        error(2);
    /// start the clock
    t1 = clock();
    /// malloc everything here
    time_result = (double *)malloc(sizeof(double)*500);
    if( time_result == NULL )
        error(1);
    q = (double *)malloc(sizeof(double)*m);
    if( q == NULL )
        error(1);
    qo = (double *)malloc(sizeof(double)*m);
    if( qo == NULL )
        error(1);
    uo = (double *)malloc(sizeof(double)*m);
    if( uo == NULL )
        error(1);
    lo = (double *)malloc(sizeof(double)*m);
    if( lo == NULL )
        error(1);

    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);

    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);
    T_tmp = (Index *)malloc(sizeof(Index)*m);
    if( T_tmp == NULL )
        error(1);

    u = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);

    l = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);

    cb = (double *)malloc(sizeof(double)*m);
    if( cb == NULL )
        error(1);

    cb1 = (double *)malloc(sizeof(double)*m);
    if( cb1 == NULL )
        error(1);

    cb2 = (double *)malloc(sizeof(double)*m);
    if( cb2 == NULL )
        error(1);

    u_d = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);

    l_d = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);

    t = (double *)malloc(sizeof(double)*m*2);
    if( t == NULL )
        error(1);
    tz = (double *)malloc(sizeof(double)*m);
    if( tz == NULL )
        error(1);
    buffer = (double *)malloc(sizeof(double)*EPOCH);
    if( buffer == NULL )
        error(1);
    tmpArray = (double*)malloc(sizeof(double)*(m+2));
    if( tmpArray == NULL )
        error(1);
    upperBoundArray = (double*)malloc(sizeof(double)*(m+2));
    if( upperBoundArray == NULL )
        error(1);
    u_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( u_buff == NULL )
        error(1);
    l_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( l_buff == NULL )
        error(1);


    /// Read query file
    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;
    while(fscanf(qp,"%lf",&d) != EOF && i < m)
    {
        ex += d;
        ex2 += d*d;
        q[i] = d;
        i++;
    }
    fclose(qp);
    /// Do z-normalize the query, keep in same array, q
    mean = ex/m;
    std = ex2/m;
    std = sqrt(std-mean*mean);
    for( i = 0 ; i < m ; i++ )
        q[i] = (q[i] - mean)/std;

    /// Sort the query one time by abs(z-norm(q[i]))
    for( i = 0; i<m; i++)
    {
        Q_tmp[i].value = q[i];
        Q_tmp[i].index = i;
    }
    qsort(Q_tmp, m, sizeof(Index),comp);

    /// also create another arrays for keeping sorted envelop
    for( i=0; i<m; i++)
    {   int o = Q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }

    /// Initial the cummulative lower bound
    for( i=0; i<m; i++)
    {   cb[i]=0;
        cb1[i]=0;
        cb2[i]=0;
    }
    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t
    ex = ex2 = 0;
    bool done = false;
    int it=0, ep=0, k=0;
    long long I;    /// the starting index of the data in current chunk of size EPOCH
    while(!done)
    {
        /// Read first m-1 points
        ep=0;
        if (it==0)
        {   for(k=0; k<m-1; k++)
                if (fscanf(fp,"%lf",&d) != EOF)
                    buffer[k] = d;
        }
        else
        {   for(k=0; k<m-1; k++)
                buffer[k] = buffer[EPOCH-m+1+k];
        }
        /// Read buffer of size EPOCH or when all data has been read.
        ep=m-1;
        while(ep<EPOCH)
        {   if (fscanf(fp,"%lf",&d) == EOF)
                break;
            buffer[ep] = d;
            ep++;
        }
        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.
        if (ep<=m-1)
        {   done = true;
        } else
        {
            /// Get time for epochs.
            if (it%(100000/(EPOCH-m+1))==0){
                fprintf(stderr,".");
                t3 = clock();
                time_result[tr_count] = (t3-t1)/CLOCKS_PER_SEC;
                tr_count = tr_count + 1;
            }
            /// run main task for each epoch
            ex=0; ex2=0;
            for(i=0; i<ep; i++)
            {
                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];
                /// Calcualte sum and sum square
                ex += d;
                ex2 += d*d;
                /// t is a circular array for keeping current data
                t[i%m] = d;
                /// Double the size for avoiding using modulo "%" operator
                t[(i%m)+m] = d;
                /// Start the task when there are more than m-1 points in the current chunk
                if( i >= m-1 )
                {
                    mean = ex/m;
                    std = ex2/m;
                    std = sqrt(std-mean*mean);
                    /// compute the start location of the data in the current circular array, t
                    j = (i+1)%m;
                    /// the start location of the data in the current chunk
                    I = i-(m-1);
                    /// Use a constant lower bound to prune the obvious subsequence
                    for(k=0;k<m;k++)
                    {
                        tz[k] = (t[(k+j)] - mean)/std;
                        T_tmp[k].value = tz[k];
                        T_tmp[k].index = k;
                    }
                    qsort(T_tmp, m, sizeof(Index),comp);

                    global_lb = lb_sorted(T_tmp, Q_tmp,  m, bsf);
                    if (global_lb < bsf)
                    {

                        distCalc = msmDistPruned(tz,q,m,bsf,bandwidth, tmpArray, upperBoundArray);

                        if( distCalc < bsf )
                        {   /// Update bsf
                            bsf = distCalc;
                            loc = (it)*(EPOCH-m+1) + i-m+1;
                        }
                    } else{
                        lbsorted++;
                    }
                    /// Reduce obsolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j]*t[j];
                }
            }
            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep<EPOCH)
                done=true;
            else
                it++;
        }
    }
    i = (it)*(EPOCH-m+1) + ep;
    fclose(fp);
    free(q);
    free(u);
    free(l);
    free(uo);
    free(lo);
    free(qo);
    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(tmpArray);
    free(upperBoundArray);
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);
    free(Q_tmp);
    t2 = clock();
    /// Output
    FILE *rd = NULL;
    rd = fopen("subsequence_results.csv", "a");
    fprintf(rd,"%s %s,%i,%lli,%f,%lld,%f,%d\n", "PMSMSearch with LB_Sorted with Sakoe", argv[4], m,i,bsf,loc, (t2-t1)/CLOCKS_PER_SEC, lbsorted);
    fprintf(rd,"Times for every 100000: [");
    for (int i = 0; i < tr_count; i++){
        fprintf(rd,"%f, ", time_result[i]);
    }
    fprintf(rd,"]\n");
    fprintf(rd,"[%i,%i]\n", mil_move_counter, mil_mergesplit_counter);
    fclose(rd);
    free(time_result);

    return 0;
}
