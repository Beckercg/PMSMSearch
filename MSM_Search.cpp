//
// Created by Chris on 11.11.2023.
//


#include <iostream>
#include <fstream>
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
        cout << "Command Usage:   PMSMSearch.exe  data_file  query_file query_length" << endl;
        cout << "For example  :   PMSMSearch.exe  data.tsv   query.tsv  20" << endl;
    }
    else if ( id == 2 )
        cout << "Error while open file!" << endl;
    exit(1);
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

int knn(double *query, const char *sf, int ql)
{
    double bsf = INF;            // best-so-far
    double sequence[ql];
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
            double distance = MSM_Distance(query, ql, sequence, ql);
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
    float acc;

    if (argc!=4)      error(1);

    double t1,t2;          // timer
    t1 = clock();

    FILE *qp = NULL;    //query data
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
            nclass = knn(query, argv[1], ql);
            if(nclass == qclass)   tp++;
            cout << "Query class: "<< qclass << "; 1NN class: "<< nclass << endl;
            i=-1;
            qcount++;
        }
        i++;
    }

    t2 = clock();
    acc = (float)tp/(float)qcount;
    cout << "tp: " << tp << endl;
    cout << "qcount: " << qcount << endl;
    cout << "Accuracy: " << acc << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

    char *ptr;
    ptr = strtok(argv[1], "/");
    ptr = strtok(NULL, "/");
    FILE *rd = NULL;    //result data
    rd = fopen("results.csv", "a");
    fprintf(rd,"%s, %d, %d, %f, %f secs\n", ptr, qcount, ql, acc, (t2-t1)/CLOCKS_PER_SEC);
    fclose(qp);
    return 0;
}

