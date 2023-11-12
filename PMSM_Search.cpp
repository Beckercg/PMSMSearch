//
// Created by Chris on 11.11.2023.
//


#include <iostream>
#include <fstream>

#include <ctime>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define INF 1e20       //Pseudo Infitinte number for this code


using namespace std;

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define C_COST 0.5 // cost for merge and split
#define INF 1e20   // Pseudo Infitinte number for this code

void error(int id)
{
    if(id==1)
    {
        cout << "ERROR: Invalid Number of Arguments!" << endl;
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file query_length" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv  20" << endl;
    }
    else if ( id == 2 )
        cout << "Error while open file!" << endl;
    exit(1);
}
/*
double C(double new_point, double x, double y){
    double c = 0.1; // Change this value if needed.
    double dist;
    if ( ( (x <= new_point) && (new_point <= y) ) ||
         ( (y <= new_point) && (new_point <= x) ) ){
        dist = c;
    }
    else{
        dist = c + min( fabs(new_point - x), fabs(new_point - y) );
    }
    return dist;
}

double MSM_Distance(double* X, int m, double * Y, int n){
    int i,j;
    double Cost[m][n];
    double d1, d2, d3;
    // Initialize
    Cost[0][0] = fabs( X[0] - Y[0] );
    for (i = 1; i < m; i++){
        Cost[i][0] = Cost[i-1][0] + C(X[i], X[i-1], Y[0]);
    }
    for ( j = 1; j < n; j++){
        Cost[0][j] = Cost[0][j-1] + C(Y[j], X[0], Y[j-1]);
    }
    // Main Loop
    for( i = 1; i < m; i++){
        for (j = 1; j < n; j++){
            d1 = Cost[i-1][j-1] + fabs(X[i]-Y[j]);
            d2 = Cost[i-1][j] + C(X[i], X[i-1], Y[j]);
            d3 = Cost[i][j-1] + C(Y[j], X[i], Y[j-1]);
            Cost[i][j] = min( d1, min(d2,d3) );
        }
    }
    // Output
    return Cost[m-1][n-1];
}
*/

vector<double> calculateMsmGreedyArray(const vector<double> &X, const vector<double> &Y)
{
    const std::vector<double>::size_type m = X.size();

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
 */
double msmDistPruned(const vector<double> &X, const vector<double> &Y)
{
    const std::vector<double>::size_type m = X.size();

    vector<double> upperBoundArray = calculateMsmGreedyArray(X, Y);
    double upperBound = upperBoundArray[0] + 0.0000001;

    vector<double> ts1 = std::vector<double>(1, INF);
    vector<double> ts2 = std::vector<double>(1, INF);

    ts1.reserve(m + 1);
    ts2.reserve(m + 1);

    ts1.insert(ts1.end(), X.begin(), X.end());
    ts2.insert(ts2.end(), Y.begin(), Y.end());

    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules
    vector<double> tmpArray = std::vector<double>(m + 1, INF);

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
    for (std::vector<double>::size_type i = 1; i < tmpArray.size(); i++)
    {

        // compute bandwidth regarding the upper bound
        unsigned int bandwidth = computeBandwidth(upperBound);
        unsigned int start = (bandwidth > i) ? sc : max(sc, i - bandwidth);

        unsigned int end = min(i + bandwidth + 1, tmpArray.size());

        double xi = ts1[i];
        // the index for the pruned end cannot be lower than the diagonal
        // All entries on the diagonal have to be equal or smaller than
        // the upper bound (Euclidean distance = diagonal path)
        ecNext = i;
        smallerFound = false;

        // column index
        for (std::vector<double>::size_type j = start; j < end; j++)
        {

            double yj = ts2[j];
            double d1, d2, d3;
            d1 = tmp + abs(xi - yj);
            // merge
            d2 = tmpArray[j] + C(xi, ts1[i - 1], yj);
            // split
            d3 = tmpArray[j - 1] + C(yj, xi, ts2[j - 1]);

            // store old entry before overwriting
            tmp = tmpArray[j];
            tmpArray[j] = min(d1, min(d2, d3));

            // PruningExperiments strategy
            double lb = getLowerBound(i, j);
            if ((tmpArray[j] + lb) > upperBound)
            {
                if (!smallerFound)
                    sc = j + 1;
                if (j > ec)
                {
                    std::fill(tmpArray.begin() + j + 1, tmpArray.end(), INF);
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
        std::fill(tmpArray.begin() + 1, tmpArray.begin() + sc, INF);

        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }

    return tmpArray[m];
}


int knn(vector<double> query, const char *sf, int ql)
{
    double bsf = INF;            // best-so-far
    double sequence[ql-1];
    double sval;
    int sclass, bclass;
    FILE *sp = NULL;
    if (NULL == (sp = fopen(sf,"r")))   error(2);
    long long j;
    j = 0;
    while(fscanf(sp,"%lf",&sval) != EOF )
    {
        if(j == 0) sclass = sval;
        else{
            sequence[j] = sval;
        }
        if(j==ql)
        {
            double distance = msmDistPruned(query, vector<double>(sequence, sequence + sizeof sequence / sizeof sequence[0]));
            if(distance < bsf)
            {
                bsf = distance;
                bclass = sclass;
            }
            j=-1;
        }
        j++;
    }

    fclose(sp);
    return bclass;
}

int main(  int argc , char *argv[] )
{
    double qval;
    long long i, nearest;
    int qclass, nclass;
    int tp, qcount;
    double acc;

    if (argc!=4)      error(1);

    double t1,t2;          // timer
    t1 = clock();

    FILE *qp = NULL;
    if (NULL == (qp = fopen(argv[2],"r")))   error(2);

    int ql;                 // length of query
    ql = atoi(argv[3]);
    double query[ql-1];
    i = 0;
    qcount = 0;
    tp = 0;
    while(fscanf(qp,"%lf",&qval) != EOF )
    {
        if(i == 0) qclass = qval;
        else{
            query[i] = qval;
        }
        if(i==ql)
        {
            nclass = knn(vector<double>(query, query + sizeof query / sizeof query[0]), argv[1], ql);
            if(nclass == qclass)   tp++;
            cout << "Query class: "<< qclass << "; 1NN class: "<< nclass << endl;
            cout << "tp: " << tp << endl;
            i=-1;
            qcount++;
        }
        i++;
    }
    fclose(qp);

    t2 = clock();
    acc = tp/qcount;
    cout << "tp: " << tp << endl;
    cout << "qcount: " << qcount << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    return 0;
}
