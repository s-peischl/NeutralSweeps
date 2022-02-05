#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <string.h>
#include "MersenneTwister.h"

using namespace std;

int main()
{
    printf("Hello world \n");
    MTRand mtrand1;
    double r;
    r = mtrand1.randDblExc();
    cout << r << endl;
    int N(20);

    cout << "number of individuals N?" << endl;
    cin >> N;
    cin.ignore();


    int k_rand;
    int t = 0;
    int N_iter;
    int fixation_time;
    int Proba_tau = 0;




    vector<int> Gene(2*N);
    vector<int> New_Gene(2*N);

    string filename;

    cout << "number of iterations?" << endl;
    cin >> N_iter;
    cin.ignore();


    cout << "time to fixation?" << endl;
    cin >> fixation_time;
    cin.ignore();

    vector< vector<int> > dem_hist(fixation_time+1,vector<int>(2*N,0));
    //vector< vector<int> > dem_hist(fixation_time+1,vector<int>(2,0));
    vector<int> avg_dem(fixation_time+1,0);

    string str1 = to_string(N);

    string str2 = to_string(fixation_time);

    filename = "demHistN" + str1 + "tau" + str2 + "Distr.dat";

    cout << filename << endl;

    ofstream files(filename);
    files.close();



    for (int j = 0; j < N_iter; j++)
    {
        t = 0;

        Gene = vector<int>(2*N,0);

        //for (int i = 0; i < 8; i++)
        for (int i = 0; i < 2*N; i++)
        {
            //Gene[i] = 1;
            Gene[i] = i;
        }

        for (int i = 0; i <= fixation_time; i++)
        {
            //dem_hist[i] = vector<int>(2,0);
            dem_hist[i] = vector<int>(2*N,0);
        }


        while (t < fixation_time)             /////////// run until only one allele left on locus 1
        {
            t++;

            for (int i = 0; i < 2*N; i++)    ////////// choose 2N new haplotypes
            {
                k_rand = mtrand1.randInt(2*N-1);
                New_Gene[i] = Gene[k_rand];
                dem_hist[t][New_Gene[i]]++;
            }

            for (int i = 0; i < 2*N; i++)
            {
                Gene[i] = New_Gene[i];
            }
        } /////////// end while loop


        for (int i = 0; i < 2*N; i++)
        {
            //if (dem_hist[t][1] == 2*N)
            if (dem_hist[t][i] == 2*N)
            {
                //if (dem_hist[t-1][1] != 2*N)
                if (dem_hist[t-1][i] != 2*N)
                {
                    ofstream iobitch (filename,ios::app);
                    for (int j = 1; j <= t; j++)
                    {
                        iobitch << dem_hist[j][i] << " ";
                    }
                    iobitch << endl;
                    iobitch.close();
                }
            }
        }
    }

    cout << 1.*Proba_tau/N_iter << endl;





    return 0;

}

//plot average diversity as a function of time and diversity on the linked locus at fixation (absolute and relative to the expectation at the correponding time)


