//
// Created by Chris on 11.11.2023.
//


#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <cstring>

#define INF 1e20       //Pseudo Infitinte number for this code


using namespace std;

void error(int id)
{
    if(id==1)
    {
        cout << "ERROR: Invalid Number of Arguments!" << endl;
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file bandwidthInPerc" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv  50" << endl;
    }
    else if ( id == 2 )
        cout << "Error while open file!" << endl;
    else if ( id == 2 )
        cout << "Error with Bandwidth" << endl;
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

double MSM_Distance(vector<double> X, vector<double> Y, int bandwidth){
    int i,j;
    int m = X.size();
    int n = Y.size();
    double Cost[m][n];
    double d1, d2, d3;

    int start;
    int end;
    // Initialize
    Cost[0][0] = fabs( X[0] - Y[0] );
    for (i = 1; i < m; i++){
        Cost[i][0] = Cost[i-1][0] + C(X[i], X[i-1], Y[0]);
    }
    for ( j = 1; j < n; j++){
        Cost[0][j] = Cost[0][j-1] + C(Y[j], X[0], Y[j-1]);
    }
    // length of query
    // Main Loop
    for( i = 1; i < m; i++){
        start = max(1,i-bandwidth);
        end = min(n,i+bandwidth);
        if(start>end)
            error(3);
        for (j = start; j < end; j++){
            d1 = Cost[i-1][j-1] + fabs(X[i]-Y[j]);
            d2 = Cost[i-1][j] + C(X[i], X[i-1], Y[j]);
            d3 = Cost[i][j-1] + C(Y[j], X[i], Y[j-1]);
            Cost[i][j] = min( d1, min(d2,d3) );
        }
    }
    // Output
    return Cost[m-1][n-1];
}



int knn(vector<double> query, vector<vector<double>> sequencefile, vector<int> sclass, int bandwidth)
{
    double bsf = INF;            // best-so-far
    double distance;
    int bclass;
    for (vector<double>::size_type j = 0; j < sequencefile.size(); j++){
        distance = MSM_Distance(query, sequencefile[j], bandwidth);
        if(distance < bsf)
        {
            bsf = distance;
            bclass = sclass[j];
        }
    }
    return bclass;
}
/*
double crossValidate(vector<vector<double>> sequencefile, double bandwidth, vector<int> sclass) {
    double averageAccuracy;
    double bsf = INF;            // best-so-far
    double distance;
    int bclass;
    for (vector<double>::size_type j = 0; j < sequencefile.size(); j++){

        distance = MSM_Distance(sequencefile[j], sequencefile[j], bandwidth);
        if(distance < bsf)
        {
            bsf = distance;
            bclass = sclass[j];
        }
    }
    return averageAccuracy;
}
*/
int crossValidate(vector<vector<double>>& sequenceFile, vector<int>& sclass, int bandwidth) {
    int foldSize = sequenceFile.size() / 10;
    int tp = 0;
    int pcount = 0;

    // Shuffling the indices for random splitting
    vector<int> indices(sequenceFile.size());
    iota(indices.begin(), indices.end(), 0);
    random_device rd;
    mt19937 g(rd());
    shuffle(indices.begin(), indices.end(), g);

    for (int i = 0; i < 10; ++i) {
        cout << "cv iteration: " << i << endl;
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
    }

    // Calculating the accuracy
    return tp * 100 / (double) pcount;
}
double findOptimalBandwidth(vector<vector<double>> sequencefile, vector<int> sclass) {
    double bestBandwidth = sequencefile[0].size();
    double bestAccuracy = 0;

    for (int bandwidth = 0; bandwidth <= sequencefile[0].size(); bandwidth += 10) {
        double accuracy = crossValidate(sequencefile, sclass, bandwidth);
        if (accuracy > bestAccuracy) {
            bestAccuracy = accuracy;
            bestBandwidth = bandwidth;
            cout << "bestBandwidth: " << bestBandwidth << endl;
        }
    }
    double accuracy = crossValidate(sequencefile, sclass, sequencefile[0].size());
    if (accuracy > bestAccuracy) {
        bestAccuracy = accuracy;
        bestBandwidth = sequencefile[0].size();
        cout << "bestBandwidth: " << bestBandwidth << endl;
    }
    return bestBandwidth;
}

int main(  int argc , char *argv[] )
{
    vector<vector<double>> queryfile;
    vector<vector<double>> sequencefile;
    vector<int> qclass;
    vector<int> sclass;
    double qval;
    long long i, nearest;
    int nclass;
    int tp, qcount;
    float acc;
    if (argc!=3)      error(1);

    double t1,t2;          // timer
    t1 = clock();


    qcount = 0;
    tp = 0;

    if (!readFiles(argv[1], sequencefile, sclass))
    {
        error(2);
        return 0;
    }
    if (!readFiles(argv[2], queryfile, qclass))
    {

        error(2);
        return 0;
    }


    double bandwidth = findOptimalBandwidth(sequencefile, sclass);
    printf("Read data is compose of %d examples with %d observations each\n\n", queryfile.size(), queryfile[1].size());


    t1 = clock();
    for (vector<double>::size_type i = 0; i < queryfile.size(); i++){
        nclass = knn(queryfile[i], sequencefile, sclass, bandwidth);
        if(nclass == qclass[i])   tp++;
        qcount++;
        cout << "qcount " << qcount << endl;
    }
    t2 = clock();
    acc = (float)tp/(float)qcount;
    cout << "tp: " << tp << endl;
    cout << "qcount: " << qcount << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;    char *ptr;

    ptr = strtok(argv[1], "/");
    ptr = strtok(NULL, "/");
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");
    fprintf(rd,"%s%i, %s, %d, %d, %f, %f secs\n", "MSM learned Sakoe Band: ", bandwidth, ptr, qcount, queryfile[1].size(), acc, (t2-t1)/CLOCKS_PER_SEC);
    return 0;
}