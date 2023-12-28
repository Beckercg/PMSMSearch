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
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file bandwidth(in Percent)" << endl;
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


double ED_Distance(vector<double> X, vector<double> Y){
    int i,j;
    int m = X.size();
    int n = Y.size();
    double distance;
    distance = 0;
    for (i = 0; i < m; i++)
        distance += abs(X[i] - Y[i]);


    return distance;
}


int knn(vector<double> &query, vector<vector<double>> &sequencefile, vector<int> &sclass)
{
    double bsf = INF;            // best-so-far
    double distance;
    int bclass;
    for (vector<double>::size_type j = 0; j < sequencefile.size(); j++){
        distance = ED_Distance(query, sequencefile[j]);
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

    string dataset, querypath ,sequencepath;
    dataset = argv[1];
    querypath = "data/" + dataset + "/" + dataset + "_TEST.tsv";
    sequencepath = "data/" + dataset + "/" + dataset + "_TRAIN.tsv";

    tp = 0;
    if(!readData(*querypath.c_str(),
                 *sequencepath.c_str(),
                 queryfile,
                 sequencefile,
                 qclass,
                 sclass))
        return 0;
    t1 = clock();
    for (vector<double>::size_type i = 0; i < queryfile.size(); i++){
        nclass = knn(queryfile[i], sequencefile, sclass);
        if(nclass == qclass[i])   tp++;
    }
    t2 = clock();
    acc = tp/(float)queryfile.size();
    cout << "tp: " << tp << endl;
    cout << "qcount: " << queryfile.size() << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");
    fprintf(rd,"%s, %s, %d, %d, %f, %f secs\n", "PMSM",dataset.c_str(), queryfile.size(), queryfile[0].size(), sequencefile.size(), sequencefile[0].size(),acc, (t2-t1)/CLOCKS_PER_SEC);
    fclose(rd);
    return 0;
}

