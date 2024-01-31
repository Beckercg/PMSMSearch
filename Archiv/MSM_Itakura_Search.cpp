//
// Created by Chris on 11.11.2023.
//


#include <iostream>
#include <vector>
#include <fstream>

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
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file slope" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv  0.1" << endl;
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

double C(double new_point, double x, double y){
    double c = 0.5; // Change this value if needed.
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

double MSM_Distance(vector<double> &X, vector<double> &Y, double &slope){
    int i,j;
    int m = X.size();
    int n = Y.size();
    double Cost[m][n];
    double d1, d2, d3;

    for ( i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            Cost[i][j] = INF;
        }
    }

    // Initialize
    Cost[0][0] = fabs( X[0] - Y[0] );
    // length of query
    // Main Loop
    for( i = 1; i < m; i++){
        for (j = 1; j < n; j++){

            if (j > (slope * i) || i > (slope * j)) continue;
            d1 = Cost[i-1][j-1] + fabs(X[i]-Y[j]);
            d2 = Cost[i-1][j] + C(X[i], X[i-1], Y[j]);
            d3 = Cost[i][j-1] + C(Y[j], X[i], Y[j-1]);
            Cost[i][j] = min( d1, min(d2,d3) );
        }
    }
    // Output
    return Cost[m-1][n-1];
}


int knn(vector<double> &query, vector<vector<double>> &sequencefile, vector<int> &sclass, double &slope)
{
    double bsf = INF;            // best-so-far
    double distance;
    int bclass;
    for (vector<double>::size_type j = 0; j < sequencefile.size(); j++){
        distance = MSM_Distance(query, sequencefile[j], slope);
        if(distance < bsf)
        {
            bsf = distance;
            bclass = sclass[j];
        }
    }
    return bclass;
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

    if (argc!=4)      error(1);

    tp = 0;
    if(!readData(*argv[2],
                 *argv[1],
                 queryfile,
                 sequencefile,
                 qclass,
                 sclass))
        return 0;
    t1 = clock();
    double slope;
    slope = atoi(argv[3]);
    for (vector<double>::size_type i = 0; i < queryfile.size(); i++){
        nclass = knn(queryfile[i], sequencefile, sclass, slope);
        if(nclass == qclass[i])   tp++;
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
    fprintf(rd,"%s%d, %s, %d, %d, %f, %f secs\n", "MSM with Itakura slope: ", slope, ptr, queryfile.size(), queryfile[1].size(), acc, (t2-t1)/CLOCKS_PER_SEC);
    fclose(rd);
    return 0;
}
