//
// Created by Chris on 11.11.2023.
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define INF 1e20       //Pseudo Infitinte number for this code


using namespace std;

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define C_COST 0.1 // cost for merge and split
#define INF 1e20   // pseudo infinite number for this code



void error(int id)
{
    if(id==1)
    {
        cout << "ERROR: Invalid Number of Arguments!" << endl;
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file bw_step_size cv_split_count" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv  10 10" << endl;
    }
    else if ( id == 2 )
        cout << "Error while open file!" << endl;
    exit(1);
}

bool readFiles(const char *fileToRead,
               vector<vector<double>> &data,
               vector<int> &classValue,
               char separator = '\t') // default separator is tab
{
    ifstream fin(fileToRead);
    string line;

    if (!fin.is_open())
    {
        return false;
    }

    while (getline(fin, line))
    {
        istringstream lineStream(line);
        string value;
        vector<double> auxTS;
        double auxVal;
        bool isclassValue = true;

        while (getline(lineStream, value, separator))
        {
            istringstream(value) >> auxVal;
            if (isclassValue)
            {
                classValue.push_back(auxVal);
                isclassValue = false;
            }
            else
            {
                auxTS.push_back(auxVal);
            }
        }
        data.push_back(auxTS);
    }
    return true;
}


bool readData(const char &queryToRead,
              const char &sequenceToRead,
              vector<vector<double>> &queryfile,
              vector<vector<double>> &sequencefile,
              vector<int> &qclass,
              vector<int> &sclass) {
    if (!readFiles(&sequenceToRead, sequencefile, sclass))
    {
        error(2);
        return 0;
    }
    printf("Sequence data is compose of %d examples with %d observations each\n", sequencefile.size(), sequencefile[0].size());
    if (!readFiles(&queryToRead, queryfile, qclass))
    {
        error(2);
        return 0;
    }
    printf("Query data is compose of %d examples with %d observations each\n\n", queryfile.size(), queryfile[0].size());
    return 1;
}


vector<double> calculateMsmGreedyArray(const vector<double> &X, const vector<double> &Y)
{
    const vector<double>::size_type m = X.size();

    vector<double> greedyArray;
    // greedyArray.reserve(m + 1);
    greedyArray = vector<double>(m + 1, 0);

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

    for (vector<double>::size_type i = 2; i <= m; i++)
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
double msmDistPruned(const vector<double> &X, const vector<double> &Y, const int sakoe_bandwidth)
{

    const vector<double>::size_type m = X.size();

    vector<double> upperBoundArray = calculateMsmGreedyArray(X, Y);
    double upperBound = upperBoundArray[0] + 0.0000001;

    vector<double> ts1 = vector<double>(1, INF);
    vector<double> ts2 = vector<double>(1, INF);

    ts1.reserve(m + 1);
    ts2.reserve(m + 1);

    ts1.insert(ts1.end(), X.begin(), X.end());
    ts2.insert(ts2.end(), Y.begin(), Y.end());


    // Create an array with one extra entry, regarding the whole matrix we initialize
    // MsmDistAStar.Entry [0,0] is set to 0
    // the first row and the first column with inf --> Every entry follows the same computational rules
    vector<double> tmpArray = vector<double>(m + 1, INF);

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
    for (vector<double>::size_type i = 1; i < tmpArray.size(); i++)
    {

        // compute bandwidth regarding the upper bound
        unsigned int local_bandwidth = computeBandwidth(upperBound);
        if (sakoe_bandwidth < local_bandwidth) local_bandwidth = sakoe_bandwidth;

        unsigned int start = max(1,i-sakoe_bandwidth);
        unsigned int end = min(m,i+sakoe_bandwidth);
        //unsigned int start = (local_bandwidth > i) ? sc : max(sc, i - local_bandwidth);

        //unsigned int end = min(i + local_bandwidth + 1, tmpArray.size());

        double xi = ts1[i];
        // the index for the pruned end cannot be lower than the diagonal
        // All entries on the diagonal have to be equal or smaller than
        // the upper bound (Euclidean distance = diagonal path)
        ecNext = i;
        smallerFound = false;

        // column index
        for (vector<double>::size_type j = start; j < end; j++)
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
            if(tmpArray[j] != d1 ) cout << "false" << endl;

            // PruningExperiments strategy
            double lb = getLowerBound(i, j);
            if ((tmpArray[j] + lb) > upperBound)
            {
                if (!smallerFound)
                    sc = j + 1;
                if (j > ec)
                {
                    fill(tmpArray.begin() + j + 1, tmpArray.end(), INF);
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
        fill(tmpArray.begin() + 1, tmpArray.begin() + sc, INF);

        // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
        tmp = INF;
        ec = ecNext;
    }

    return tmpArray[m];
}


int knn(const vector<double> &query, const vector<vector<double>> &sequencefile, const vector<int> &sclass, const int &bandwidth)
{
    double bsf = INF;            // best-so-far
    double distance;
    int bclass;
    for (vector<double>::size_type j = 0; j < sequencefile.size(); j++){
        distance = msmDistPruned(query, sequencefile[j], bandwidth);
        if(distance < bsf)
        {
            bsf = distance;
            bclass = sclass[j];
        }
    }
    return bclass;
}

double crossValidate(vector<vector<double>>& sequenceFile, vector<int>& sclass, int bandwidth, const int cv_split_count) {
    int foldSize = sequenceFile.size() / cv_split_count;
    int tp = 0;
    int pcount = 0;
    double averageAcc =0;

    // Shuffling the indices for random splitting
    vector<int> indices(sequenceFile.size());
    iota(indices.begin(), indices.end(), 0);
    random_device rd;
    mt19937 g(rd());
    shuffle(indices.begin(), indices.end(), g);

    for (int i = 0; i < cv_split_count; ++i) {
        // Creating train and test sets for this fold
        vector<vector<double>> trainSet;
        vector<int> trainClass;
        vector<vector<double>> testSet;
        vector<int> testClass;
        for (int j = 0; j < sequenceFile.size(); ++j) {
            if (j >= i * foldSize && j < (i + 1) * foldSize) {
                testSet.push_back(sequenceFile[indices[j]]);
                testClass.push_back(sclass[indices[j]]);
            } else {
                trainSet.push_back(sequenceFile[indices[j]]);
                trainClass.push_back(sclass[indices[j]]);
            }
        }

        // Testing the classifier
        for (int k = 0; k < testSet.size(); ++k) {
            int predicted = knn(testSet[k], trainSet, trainClass, bandwidth);
            if (predicted == testClass[k]) {
                ++tp;
            }
            ++pcount;
        }
        averageAcc +=  tp * 100 / (double) pcount * (1 / (double) cv_split_count);
        cout << "Acc: " << tp * 100 / (double) pcount << "for iteration: " << i << endl;
        cout << "avAcc: " << averageAcc << endl;
    }

    // Calculating the accuracy
    return averageAcc;
}


int findOptimalBandwidth(vector<vector<double>> sequencefile, vector<int> sclass, const int &bw_step_size, const int &cv_split_count) {
    int bestBandwidth = sequencefile[0].size();
    double bestAccuracy = 0;

    for (int bandwidth = 0; bandwidth <= sequencefile[0].size(); bandwidth += bw_step_size) {
        cout << "bandwidth: " << bandwidth << endl;
        double accuracy = crossValidate(sequencefile, sclass, bandwidth, cv_split_count);
        if (accuracy > bestAccuracy) {
            bestAccuracy = accuracy;
            bestBandwidth = bandwidth;
            cout << "bestBandwidth: " << bestBandwidth << endl;
        }
    }
    return bestBandwidth;
}


int main(  int argc , char *argv[] )
{
    vector<vector<double>> queryfile;
    vector<vector<double>> sequencefile;
    vector<int> qclass;
    vector<int> sclass;
    int nclass, tp;
    float acc;
    double t1,t2;          // timer

    if (argc!=5)      error(1);

    if(!readData(*argv[2],
                 *argv[1],
                 queryfile,
                 sequencefile,
                 qclass,
                 sclass))
        return 0;

    tp = 0;

    int bw_step_size;
    bw_step_size = atoi(argv[3]);
    int cv_split_count;
    cv_split_count = atoi(argv[4]);
    int bandwidth = findOptimalBandwidth(sequencefile, sclass, bw_step_size, cv_split_count);

    t1 = clock();
    for (vector<double>::size_type i = 0; i < queryfile.size(); i++){
        nclass = knn(queryfile[i], sequencefile, sclass, bandwidth);
        if(nclass == qclass[i])   tp++;
        cout << "queryfile: " << i << endl;
    }
    t2 = clock();
    acc = tp/(float)queryfile.size();
    cout << "tp: " << tp << endl;
    cout << "qcount: " << queryfile.size() << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    char *ptr;
    ptr = strtok(argv[1], "/");
    ptr = strtok(NULL, "/");
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");
    fprintf(rd,"%s%d , %s, %d, %d, %f, %f secs\n", "PMSM with learned bandwith: ", bandwidth, ptr, queryfile.size(), queryfile[0].size(), acc, (t2-t1)/CLOCKS_PER_SEC);

    return 0;
}

