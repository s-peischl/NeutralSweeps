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

    double rec(0.005);
    int k_rand;
    int t = 0;
    int diversity,old_diversity,diversity2;
    vector<int> div(40*N,0);
    int N_iter(1000000);
    int k_rand2,ni_choose2,coal_time,proba_root0,coal_time_max;
    proba_root0 = 0;
    unsigned long long int mean_coal_time;

    k_rand = 10;
    cout << 4*(k_rand/2)+1-k_rand << endl;
    k_rand = 13;
    cout << 4*(k_rand/2)+1-k_rand << endl;

    int Proba_t = 0;
    int tau;
    vector<int> div2(40*N,0);

    vector< vector<int> > history1(2*N);
    vector< vector<int> > history2(2*N);
    vector< vector<int> > New_history1(2*N);
    vector< vector<int> > New_history2(2*N);

    history1[0].push_back(1);
    history1[0] = history1[1];
    cout << history1[0].size() << endl;

    vector<int> is_there(2*N,0);
    vector<int> is_there2(2*N,0);

    vector<double> theory_time(40*N,0);
    vector<double> theory_proba(40*N,0);
    vector<double> sum_theory_proba(40*N,0);



    string filename;



    cout << "number of iterations?" << endl;
    cin >> N_iter;
    cin.ignore();

    cout << "time to fixation?" << endl;
    cin >> tau;
    cin.ignore();

    vector<unsigned long long int> ancestry(tau+1,0);


    rec = 0;

    string str1 = to_string(N);

    string str = to_string(tau);


/*
    str2.erase(str2.begin()+1);     //////////////// remove the point after the first zero


    string str4;
    str4.append(str2.end()-1,str2.end());

    while (str4 == "0" && str2.length() > 1)        ////////////// remove the zeros after last decimal (always 6 decimals)
    {
        str2.erase(str2.end()-1);
        str4.clear();
        str4.append(str2,str2.length()-1,1);
    }
*/
    filename = "TAdistrN" + str1 + "tau" + str + ".dat";

    cout << filename << endl;

    ofstream diversity_files (filename);
        //diversity_files << 0 << " " << 1 << " " << 0 << " " << 0 << endl;
        diversity_files << 1 << " " << 2*rec*(1-rec) << " " << 0 << " " << 0 << endl;
    diversity_files.close();

    for (int j = 0; j < N_iter; j++)
    {
        t = 0;

        for (int i = 0; i < 2*N; i++)
        {
            history1[i].clear();

        }
        diversity = 2*N;


        while (diversity != 1 && t < tau+1)             /////////// run until only one allele left on locus 1
        {
            t++;

            for (int i = 0; i < 2*N; i++)
            {
                is_there[i] = 0;
            }

            for (int i = 0; i < 2*N; i++)    ////////// choose 2N new haplotypes
            {
                k_rand = mtrand1.randInt(2*N-1);

                New_history1[i] = history1[k_rand];
                New_history1[i].push_back(k_rand);
                is_there[New_history1[i][0]]++;
            }

            diversity = 0;

            for (int i = 0; i < 2*N; i++)
            {
                history1[i] = New_history1[i];

                if (is_there[i] == 2*N)
                {
                    diversity = 1;
                }
            }
        } /////////// end while loop


        if (t == tau)
        {
            coal_time_max = 0;
            for (int i = 0; i < 2*N; i++)
            {
                for (int j = i+1; j < 2*N; j++)
                {

                    coal_time = 1;                //find last common ancestor
                    while (history1[i][history1[i].size()-coal_time] != history1[j][history1[j].size()-coal_time])
                    {
                        coal_time++;
                    }
                    Proba_t++;
                    ancestry[coal_time]++;
                    if (coal_time > coal_time_max)
                    {
                        coal_time_max = coal_time;
                    }
                }
            }

            if (coal_time_max == tau)
            {
                proba_root0++;
            }
        }



        //ancestry[t] += ni_choose2;
    }

/*
    for (int i = 1; i < 40*N; i++)
    {
        theory_proba = vector<double>(i+1,1);
        theory_proba[1] = 1./(2*N-1./i*(2*N-1));

        double sum_theory_proba = theory_proba[1];

        for (int k = 2; k <= i; k++)
        {
            for (int j = 1; j < k; j++)
            {
                theory_proba[k] *= 1-1./(2*N-1.*j/i*(2*N-1));
            }
            theory_proba[k] *= 1./(2*N-1.*k/i*(2*N-1));
            sum_theory_proba += theory_proba[k];
        }

        theory_time[i] = 0;
        for (int k = 1; k <= i; k++)
        {
            theory_time[i] += k*theory_proba[k]/sum_theory_proba;
        }
    }

    for (int i = 1; i < 40*N; i++)
    {
        sum_theory_proba[i] = 0;
        for (int j = 1; j <= i; j++)
        {
            sum_theory_proba[i] += theory_proba[j];
        }

    }

    for (int i = 1; i < 40*N; i++)
    {
        theory_time[i] = 0;
        for (int j = 1; j <= i; j++)
        {
            theory_time[i] += j*theory_proba[j];
        }
        theory_time[i] /= sum_theory_proba[i];
    }
*/

        //ofstream diversity_files ("diversity.dat",ios::app);
        //    diversity_files << t << " " << diversity/2./N << endl;
        //diversity_files.close();
        //choose 2*N new alleles at random
        //calculate number of classes
        //si

        ofstream files (filename,ios::app);
        //int clock = 1;
        //while (survival[clock] > 0)
        for (int clock = 1; clock <= tau; clock++)
        {
            if (Proba_t > 0)
            {
                //files << clock << " " << div[clock]/2./N/survival[clock] << endl;
                //files << clock << " " << div2[clock]/2./N/Proba_t[clock] << " " << (double)div2[clock]/Proba_t[clock]/div[clock]*N_iter << " " << Proba_t[clock] << endl;
                files << clock << " " << Proba_t << " " << 1.*ancestry[clock]/Proba_t << " " << theory_time[clock] << endl;
                //clock++;
                //files << clock << " " << div2[clock] << " " << Proba_t[clock] << endl;
            }

        }
        files.close();

        //for (int i = 0; i < 20; i++)
        //{
        //    cout << survival[i] << endl;
        //}

    cout << proba_root0 << endl;
    return 0;

}

//plot average diversity as a function of time and diversity on the linked locus at fixation (absolute and relative to the expectation at the correponding time)


