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
    double rec(0.005);
    int k_rand;
    int t = 0;
    int diversity,old_diversity,diversity2;
    vector<int> div(40*N,0);
    int N_iter(1000000);
    int k_rand2,ni_choose2,mean_coal_time,coal_time;

    k_rand = 10;
    cout << 4*(k_rand/2)+1-k_rand << endl;
    k_rand = 13;
    cout << 4*(k_rand/2)+1-k_rand << endl;


    vector<unsigned long long int> Proba_t(40*N,0);
    vector<unsigned long long int> div2(40*N,0);
    vector<unsigned long long int> ancestry(40*N,0);

    vector<double> meanSquare(40*N,0);
    vector<double> meanProduct(40*N,0);
    vector<double> mean_coal_timeSquare(40*N,0);

    vector< vector<int> > history1(2*N);
    vector< vector<int> > history2(2*N);
    vector< vector<int> > New_history1(2*N);
    vector< vector<int> > New_history2(2*N);

    history1[0].push_back(1);
    history1[0] = history1[1];
    cout << history1[0].size() << endl;

    vector<int> is_there(2*N,0);
    vector<int> is_there2(2*N,0);



    string filename;

    cout << "number of individuals N?" << endl;
    cin >> N;
    cin.ignore();

    cout << "number of iterations?" << endl;
    cin >> N_iter;
    cin.ignore();

    cout << "recombination rate?" << endl;
    cin >> rec;
    cin.ignore();

    string str1 = to_string(N);

    string str2 = to_string(rec);



    str2.erase(str2.begin()+1);     //////////////// remove the point after the first zero


    string str4;
    str4.append(str2.end()-1,str2.end());

    while (str4 == "0" && str2.length() > 1)        ////////////// remove the zeros after last decimal (always 6 decimals)
    {
        str2.erase(str2.end()-1);
        str4.clear();
        str4.append(str2,str2.length()-1,1);
    }

    filename = "coalN" + str1 + "rec" + str2 + ".dat";

    cout << filename << endl;

    ofstream diversity_files (filename);
        //diversity_files << 0 << " " << 1 << " " << 0 << " " << 0 << endl;
        //diversity_files << 1 << " " << 2*rec*(1-rec) << " " << 0 << " " << 0 << " " << 0 << endl;
    diversity_files.close();

    for (int j = 0; j < N_iter; j++)
    {
        t = 0;

        for (int i = 0; i < 2*N; i++)
        {
            history1[i].clear();
            history2[i].clear();
        }
        diversity = 2*N;
        diversity2 = 2*N;

        while (diversity != 1 && t < 2*N)             /////////// run until only one allele left on locus 1
        {
            t++;

            for (int i = 0; i < 2*N; i++)
            {
                is_there[i] = 0;
            }

            for (int i = 0; i < 2*N; i++)    ////////// choose 2N new haplotypes
            {
                k_rand = mtrand1.randInt(2*N-1);
                r = mtrand1.randDblExc();
                k_rand2 = k_rand;
                if (r < rec)
                {
                    k_rand2 = 4*(k_rand/2)+1-k_rand;
                }

                New_history1[i] = history1[k_rand];
                New_history1[i].push_back(k_rand);
                New_history2[i] = history2[k_rand2];
                New_history2[i].push_back(k_rand2);
                is_there[New_history1[i][0]]++;
            }

            diversity = 0;

            for (int i = 0; i < 2*N; i++)
            {
                history1[i] = New_history1[i];
                history2[i] = New_history2[i];

                if (is_there[i] == 2*N)
                {
                    diversity = 1;
                }
            }
        } /////////// end while loop

        if (diversity == 1)
        {
            for (int i = 0; i < 2*N; i++)
            {
                is_there2[i] = 0;
            }

            ni_choose2 = 0;
            mean_coal_time = 0;

            for (int i = 0; i < 2*N; i++)
            {
                for (int j = i+1; j < 2*N; j++)
                {
                    if (history2[i][0] == history2[j][0])
                    {
                        coal_time = 1;                //find last common ancestor
                        while (history2[i][history2[i].size()-coal_time] != history2[j][history2[j].size()-coal_time])
                        {
                            coal_time++;
                        }
                        ni_choose2++;
                        mean_coal_time += coal_time;
                    }
                }
            }


            Proba_t[t]++;
            div2[t] += N*(2*N-1)-ni_choose2;
            ancestry[t] += mean_coal_time;

            meanSquare[t] += pow(1-1.*ni_choose2/(N*(2*N-1)),2.);
            meanProduct[t] += 1.*mean_coal_time/N/(2*N-1)*(1-1.*ni_choose2/N/(2*N-1));
            mean_coal_timeSquare[t] += 1.*mean_coal_time/N/(2*N-1)*1.*mean_coal_time/N/(2*N-1);
        }
    }


        //ofstream diversity_files ("diversity.dat",ios::app);
        //    diversity_files << t << " " << diversity/2./N << endl;
        //diversity_files.close();
        //choose 2*N new alleles at random
        //calculate number of classes
        //si

        ofstream files (filename,ios::app);
        //int clock = 1;
        //while (survival[clock] > 0)
        for (int clock = 1; clock < 40*N; clock++)
        {
            if (Proba_t[clock] > 0)
            {
                //files << clock << " " << div[clock]/2./N/survival[clock] << endl;
                //files << clock << " " << div2[clock]/2./N/Proba_t[clock] << " " << (double)div2[clock]/Proba_t[clock]/div[clock]*N_iter << " " << Proba_t[clock] << endl;
                files << clock << " " << 1.*div2[clock]/N/Proba_t[clock]/(2*N-1) << " " << Proba_t[clock] << " " << 1.*ancestry[clock]/N/Proba_t[clock]/(2*N-1) << " " << meanSquare[clock]/Proba_t[clock] << " " << mean_coal_timeSquare[clock]/Proba_t[clock] << " " << meanProduct[clock]/Proba_t[clock] << endl;
                //clock++;
                //files << clock << " " << div2[clock] << " " << Proba_t[clock] << endl;
            }

        }
        files.close();

        //for (int i = 0; i < 20; i++)
        //{
        //    cout << survival[i] << endl;
        //}


    return 0;

}

//plot average diversity as a function of time and diversity on the linked locus at fixation (absolute and relative to the expectation at the correponding time)


