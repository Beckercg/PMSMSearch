/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


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


/// Data structure for sorting the query
typedef struct Index
{   double value;
    int    index;
} Index;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct deque
{   int *dq;
    int size,capacity;
    int f,r;
};


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return abs(y->value) - abs(x->value);   // high to low
}

/// Initial the queue at the begining step of envelop calculation
void init(struct deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int)*d->capacity);
    d->f = 0;
    d->r = d->capacity-1;
}

/// Destroy the queue
void destroy(struct deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity-1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity-1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r+1)%d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity-1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r+1)%d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Calculate quick lower bound
/// Die Punkte zwischen zwei Moves die kleiner als C_COST sind haben immer mindestens Kosten von C_COST
double lb_edcost_h(double *t, double *q, int len, double bsf)
{
    double lb;
    lb = dist(t[(len-1)],q[len-1]);
    if (lb >= bsf)   return lb;
    for(int l = 2; l<len; l++){
        if(dist(t[(len-l)],q[len-l]) >= C_COST) {
            lb += C_COST;
        }else{
            lb += dist(t[(len-l)],q[len-l]);
        }
        if (lb >= bsf)   return lb;
    }
    return lb;
}

//vector<double> calculateMsmGreedyArray(const vector<double> &X, const vector<double> &Y)
double *calculateMsmGreedyArray(double *X, double *Y, int m)
{
    //const vector<double>::size_type m = X.size();
    int i, k;
    double *greedyArray;
    greedyArray = (double*)malloc(sizeof(double)*(m+1));
    for(k=0; k<m+1; k++)    greedyArray[k]=0;

    //vector<double> greedyArray;
    // greedyArray.reserve(m + 1);
    //greedyArray = vector<double>(m + 1, 0);

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

    //for (vector<double>::size_type i = 2; i <= m; i++)
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
//double msmDistPruned(const vector<double> &X, const vector<double> &Y)
double msmDistPruned(double *X, double *Y, int m, double bsf)
{

    // const vector<double>::size_type m = X.size();

    //vector<double> upperBoundArray = calculateMsmGreedyArray(X, Y);
    //double upperBound = upperBoundArray[0] + 0.0000001;
    double *upperBoundArray = calculateMsmGreedyArray(X, Y, m);
    double upperBound = upperBoundArray[0] + 0.0000001;

    /*
    vector<double> ts1 = vector<double>(1, INF);
    vector<double> ts2 = vector<double>(1, INF);

    ts1.reserve(m + 1);
    ts2.reserve(m + 1);

    ts1.insert(ts1.end(), X.begin(), X.end());
    ts2.insert(ts2.end(), Y.begin(), Y.end());
    */

    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules

    double *tmpArray;
    int i, j, k;

    tmpArray = (double*)malloc(sizeof(double)*(m+1));
    for(k=0; k<m+1; k++)    tmpArray[k]=INF;
    //vector<double> tmpArray = vector<double>(m + 1, INF);

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
    for (int i = 0; i < m+1; i++)
    {

        // compute bandwidth regarding the upper bound
        unsigned int bandwidth = computeBandwidth(upperBound);
        unsigned int start = (bandwidth > i) ? sc : max(sc, i - bandwidth);

        //unsigned int end = min(i + bandwidth + 1, tmpArray.size());
        unsigned int end = min(i + bandwidth + 1, m+1);

        double xi = X[i];
        // the index for the pruned end cannot be lower than the diagonal
        // All entries on the diagonal have to be equal or smaller than
        // the upper bound (Euclidean distance = diagonal path)
        ecNext = i;
        smallerFound = false;
        smaller_as_bsf = false;

        // column index
        for (j = start; j < end; j++)
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

        if(!smaller_as_bsf) return INF;
        // tmpArray = this.fillWithInf(1, sc, tmpArray);
        for(k=1; k<sc; k++)    tmpArray[k]=INF;
        //fill(tmpArray.begin() + 1, tmpArray.begin() + sc, INF);

        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }

    return tmpArray[m];
}

/// If expected error happens, teminated the program.
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
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
        printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
    }
    exit(1);
}

/// Main Function
int main(  int argc , char *argv[] )
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;


    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    int m=-1, r=-1;
    long long loc = 0;
    double t1,t2;
    int edcost = 0,keogh = 0, keogh2 = 0;
    double distCalc=0, lb_edcost=0, lb_k=0, lb_k2=0;
    double *buffer, *u_buff, *l_buff;
    Index *Q_tmp;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    /// If not enough input, display an error.
    if (argc<=3)
        error(4);

    /// read size of the query
    if (argc>3)
        m = atol(argv[3]);

    /// read warping windows
    if (argc>4)
    {   double R = atof(argv[4]);
        if (R<=1)
            r = floor(R*m);
        else
            r = floor(R);
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
    free(Q_tmp);

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
            /// Just for printing a dot for approximate a million point. Not much accurate.
            if (it%(1000000/(EPOCH-m+1))==0)
                fprintf(stderr,".");

            /// Do main task here..
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

                    for(k=0;k<m;k++)
                    {
                        tz[k] = (t[(k+j)] - mean)/std;
                    }
                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_edcost = lb_edcost_h(tz, q, m, bsf);
                    if (lb_edcost < bsf)
                    {

                        /// Compute MSM and early abandoning if possible

                        distCalc = msmDistPruned(tz,q,m,bsf);
                        if( distCalc < bsf )
                        {   /// Update bsf
                            bsf = distCalc;
                            loc = (it)*(EPOCH-m+1) + i-m+1;
                        }
                    } else
                        edcost++;
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
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);

    t2 = clock();

    printf("\n");
    printf("Pruned by LB_EDCostH    : %6.2f%%\n", ((double) edcost / i)*100);
    printf("PMSM Calculation     : %6.2f%%\n", 100-(((double)edcost+keogh+keogh2)/i*100));
    FILE *rd = NULL;    //result data
    rd = fopen("subsequence_results.csv", "a");
    fprintf(rd,"%s,%i,%lli,%f,%lld,%d,%f\n", "PMSMSearch with LB_EDCostH", m,i,bsf,loc,edcost, (t2-t1)/CLOCKS_PER_SEC);
    fclose(rd);

    return 0;
}
