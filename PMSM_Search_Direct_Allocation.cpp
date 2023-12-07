//
// Created by Chris on 11.11.2023.
//


#include <iostream>
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
#define INF 1e20   // pseudo infinite number for this code



void error(int id)
{
    if(id==1)
    {
        cout << "ERROR: Invalid Number of Arguments!" << endl;
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv" << endl;
    }
    else if ( id == 2 )
        cout << "Error while open file!" << endl;
    exit(1);
}


vector<double> calculateMsmGreedyArray(double *X, double *Y, int m)
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
//double msmDistPruned(const vector<double> &X, const vector<double> &Y)
double msmDistPruned(double *X, double *Y, int m, double bsf)
{

    vector<double> upperBoundArray = calculateMsmGreedyArray(X, Y, m);
    double upperBound = upperBoundArray[0] + 0.0000001;

    auto* ts1 = (double*)malloc(sizeof(double)*(m+2));
    auto* ts2 = (double*)malloc(sizeof(double)*(m+2));
    ts1[0] = INF;
    ts2[0] = INF;

    for (int i=1; i <=m+1; i++) {
        ts1[i] = X[i-1]; // Copy elements from X to ts1
        ts2[i] = Y[i-1]; // Copy elements from X to ts1
    }

    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules

    double *tmpArray;
    int i, j, k;

    tmpArray = (double*)malloc(sizeof(double)*(m+2));
    for(k=0; k<m+2; k++)    tmpArray[k]=INF;
    //vector<double> tmpArray = vector<double>(m + 1, INF);

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

        // column index
        for (j = start; j < end; j++)
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

                double tmplog = tmpArray[j];
                double tmplog2 = upperBoundArray[j];

                upperBound = tmpArray[j] + upperBoundArray[j] + 0.00001;
            }
        }

        // tmpArray = this.fillWithInf(1, sc, tmpArray);
        for(k=1; k<sc; k++)    tmpArray[k]=INF;
        //fill(tmpArray.begin() + 1, tmpArray.begin() + sc, INF);

        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }
    free(ts1);
    free(ts2);

    return tmpArray[m];
}


int main()
{
    FILE *sp;            /// data file pointer
    FILE *qp;            /// query file pointer
    vector<vector<double>> queryfile;
    vector<vector<double>> sequencefile;
    vector<double> qclass;
    vector<double> sclass;
    int nclass, bclass, tp;
    float acc;
    double t1,t2,bsf, distance;
    int m,i,j, query_size, sequence_size;
    double d;

    string dataset, querypath ,sequencepath;
    cout << "Which datset: ";
    cin >> dataset;
    cout << "Datasetlength: ";
    cin >> m;
    cout << "How many queries: ";
    cin >> query_size;
    cout << "How many sequences: ";
    cin >> sequence_size;
    querypath = "data/" + dataset + "/" + dataset + "_TEST.tsv";
    sequencepath = "data/" + dataset + "/" + dataset + "_TRAIN.tsv";

    double** q_file = new double*[query_size+1];
    if( q_file == NULL )
        error(1);
    for(i = 0; i < query_size; i++) {
        q_file[i] = new double[m+1];
    }
    double** s_file = new double*[sequence_size+1];
    if( s_file == NULL )
        error(1);
    for(i = 0; i < sequence_size; i++) {
        s_file[i] = new double[m+1];
    }


    qp = fopen(querypath.c_str(),"r");
    i= 0;
    j=-1;
    while(fscanf(qp,"%lf",&d) != EOF && i < (m+1)*query_size)
    {

        if(i%(m+1)==0){
            j++;
            i=0;
            qclass.push_back(d);
        }
        else {
            q_file[j][i-j-1] = d;
        }
        i++;
    }
    fclose(qp);
    sp = fopen(sequencepath.c_str(),"r");
    i=0;
    j=-1;
    while(fscanf(sp,"%lf",&d) != EOF && i < (m+1)*sequence_size)
    {
        if(i%(m+1)==0){
            j++;
            i=0;
            sclass.push_back(d);
        }
        else {
            s_file[j][i-j-1] = d;
        }
        i++;
    }

    fclose(sp);

    tp = 0;
    t1 = clock();
    for (int i = 0; i < query_size; i++){
        bsf = INF;
        for (int j = 0; j < sequence_size; j++){
            distance = msmDistPruned(q_file[i], s_file[j], m, bsf);
            cout << "distance"<< distance<< endl;
            //knn
            if(distance < bsf)
            {
                bsf = distance;
                bclass = sclass[j];
            }
        }
        if(bclass == qclass[i])   tp++;
    }

    for(int i = 0; i < query_size; i++) {
        delete[] q_file[i];
    }
    delete[] q_file;
    for(int i = 0; i < sequence_size; i++) {
        delete[] s_file[i];
    }
    delete[] s_file;
    t2 = clock();



    acc = tp/query_size;
    cout << "tp: " << tp << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");

    fprintf(rd,"%s, %s, %d, %d, %f, %f secs\n", "PMSM",dataset.c_str(), query_size, m, acc, (t2-t1)/CLOCKS_PER_SEC);
    fclose(rd);
    return 0;
}

