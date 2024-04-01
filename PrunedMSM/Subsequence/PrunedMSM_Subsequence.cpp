#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) (fabs(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code
#define C_COST 0.5 // cost for merge and split

int mergesplit_counter = 0; // global merge and split counter
int move_counter = 0; // global move counter
int mil_mergesplit_counter = 0; // global merge and split counter
int mil_move_counter = 0; // global move counter

vector<double> calculateMsmGreedyArray(const vector<double> &X, const vector<double> &Y,int m)
{

    vector<double> greedyArray;
    // greedyArray.reserve(m + 1);
    greedyArray = std::vector<double>(m + 1, 0);

    // compute upper Bounds for every diagonal entry
    // the upper bound is computed from right to left: Possibility to update the upper Bound when a diagonal entry is computed
    // upper Bound for the last entry is 0, because there is nothing left to compute
    greedyArray[m] = 0.;

    // assume that the time series are aligned at the end
    double distCurrent = abs(X[m - 1] - Y[m - 1]);
    double distTmp = distCurrent;
    double xCurrent;
    double yCurrent;
    double xTmp = X[m - 1];
    double yTmp = Y[m - 1];
    int rel = (xTmp > yTmp) ? 1 : 2;

    greedyArray[m - 1] = distCurrent;

    for (std::vector<double>::size_type i = 2; i <= m; i++)
    {
        xCurrent = X[m - i];
        yCurrent = Y[m - i];
        distCurrent = xCurrent - yCurrent;

        if ((rel == 1) && distCurrent > 0)
        {
            if (distCurrent > 2 * C_COST && distTmp > 2 * C_COST)
            {
                greedyArray[m - i] = 2 * C_COST + abs(xCurrent - xTmp) + abs(yCurrent - yTmp) + greedyArray[m - i + 1];
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
            if (abs(distCurrent) > 2 * C_COST && abs(distTmp) > 2 * C_COST)
            {
                greedyArray[m - i] = 2 * C_COST + abs(xCurrent - xTmp) + abs(yCurrent - yTmp) + greedyArray[m - i + 1];
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
        return C_COST + min(abs(new_point - x), abs(new_point - y));
    }

    return C_COST;
}

double getLowerBound(int xCoord, int yCoord)
{

    return abs(xCoord - yCoord) * C_COST;
}

/**
 * compute the Msm distance table with pruned entries
 *
 * @return pair: first Double: Distance, second Double: relative amount of pruned cells
 * msmDistPruned by Jana Holznigenkemper
 */
double msmDistPruned(const vector<double> &X, const vector<double> &Y, int m)
{

    vector<double> upperBoundArray = calculateMsmGreedyArray(X, Y,m);
    double upperBound = upperBoundArray[0] + 0.0000001;


    int i,j,k;
    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules
    vector<double> tmpArray = std::vector<double>(m + 2);

    for(k=0; k<m+1; k++)    tmpArray[k]=INF;
    // value storing the first "real value" of the array before overwriting it
    //  the first value of the first row has to be 0
    double tmp = 0;

    // index for pruning start and end of row
    unsigned int sc = 1;
    unsigned int ec = 1;

    // remember if an entry smaller than UB was found -> cannot cut
    bool smallerFound;
    int ecNext;

    // initialize first row
    //  the first entry of this row is inf, this inf is used for computing the second row
    //  Arrays.fill(tmpArray, Double.POSITIVE_INFINITY); --> this was done at the initialize step

    //  int counterBandwidth =0;
    // row index
    for (i = 0; i < m+1; i++)
    {

        // compute bandwidth regarding the upper bound
        unsigned int bandwidth = computeBandwidth(upperBound);
        unsigned int start = (bandwidth > i) ? sc : max(sc, i - bandwidth);

        unsigned int end = min(i + bandwidth + 1, m+1);

        double xi = X[i];

        // the index for the pruned end cannot be lower than the diagonal
        // All entries on the diagonal have to be equal or smaller than
        // the upper bound (Euclidean distance = diagonal path)
        ecNext = i;
        smallerFound = false;

        // column index
        for (std::vector<double>::size_type j = start; j < end; j++)
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

        // tmpArray = this.fillWithInf(1, sc, tmpArray);

        for(k=1; k<sc; k++)    tmpArray[k]=INF;
        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }

    fprintf(stderr,"%f\n", tmpArray[m]);
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
    int tr_count = 0;


    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    int m=-1;
    long long loc = 0;
    double t1,t2,t3;
    double distCalc=0;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    /// If not enough input, display an error.
    if (argc<=3)
        error(4);

    /// read size of the query
    if (argc>3)
        m = atol(argv[3]);

    vector<double> tz;
    tz = std::vector<double>(m, 0);
    vector<double> t;
    t = std::vector<double>(2*m, 0);
    vector<double> q;
    q = std::vector<double>(m, 0);
    vector<double> buffer;
    buffer = std::vector<double>(EPOCH, 0);
    vector<double> time_result;
    time_result = std::vector<double>(m*20, 0);

    fp = fopen(argv[1],"r");
    if( fp == NULL )
        error(2);

    qp = fopen(argv[2],"r");
    if( qp == NULL )
        error(2);

    /// start the clock
    t1 = clock();




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
            if (it%(100000/(EPOCH-m+1))==0){
                fprintf(stderr,".");
                t3 = clock();
                time_result[tr_count] = (t3-t1)/CLOCKS_PER_SEC;
                tr_count = tr_count + 1;
            }

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

                    /// Use a constant lower bound to prune the obvious subsequence
                    for(k=0;k<m;k++)
                    {
                        tz[k] = (t[(k+j)] - mean)/std;
                    }
                    distCalc = msmDistPruned(tz,q,m);

                    if( distCalc < bsf )
                    {   /// Update bsf
                        bsf = distCalc;
                        loc = (it)*(EPOCH-m+1) + i-m+1;

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


    t2 = clock();

    FILE *rd = NULL;    //result data
    rd = fopen("subsequence_results.csv", "a");
    fprintf(rd,"%s,%i,%lli,%f,%lld,%f\n", "PrunedMSM", m,i,bsf,loc, (t2-t1)/CLOCKS_PER_SEC);
    fprintf(rd,"Times for every 100000: [");
    for (i = 0; i < tr_count; i++){
        fprintf(rd,"%f, ", time_result[i]);
    }
    fprintf(rd,"[%i,%i]\n", mil_move_counter, mil_mergesplit_counter);
    fprintf(rd,"]\n");
    fclose(rd);

    return 0;
}
