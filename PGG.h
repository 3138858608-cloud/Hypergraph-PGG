#pragma once
#ifndef AVPGG_H
#define AVPGG_H

#include"stdc++.h"
#include"HI.h"
#include"rand.h"
#include <algorithm>
using namespace std;
namespace PGG_avpi {

    int c = 1;//cost
    double fc = 0;
    vector<bool> Vectex_status_list;
    vector<bool> temp_Vectex_status_list;
    vector<double> Vectex_payoff_list;
    vector<int> Status_list;
    map<int, double>num_fc_k;
    vector<double>p;//timescale

    void Unconditional_Imitation(double r, HypergraphLib::Hypergraph B)
    {
        for (int u = 0; u < B.N; u++)
        {

            //B.Groups[u]
            //random choose hyperedge
            if (B.Groups[u].size() == 0) continue;
            int v = (int)_rand(0.0, 1.0 * B.Groups[u].size());
            int random_edge = B.Groups[u][v];
            int max_member = u;
            double max_payoff = 0;
            for (auto& v : B.Hyperedges[random_edge])
            {
                if (max_payoff < Vectex_payoff_list[v])
                {
                    max_member = v;
                    max_payoff = Vectex_payoff_list[v];
                }
            }
            if (u == max_member) continue;
            double delta = r * (B.g - 2) + 1;
            if (r > 1)delta = r * B.g - 1;
            double w = (Vectex_payoff_list[max_member] - Vectex_payoff_list[u]) / delta;
            w = min(max(w, 0.0), 1.0);
            if (_rand(0, 1) < w)
                Vectex_status_list[u] = temp_Vectex_status_list[max_member];
        }
    }
    void Betters_Possess_Chance(double r, HypergraphLib::Hypergraph B)
    {
        for (int u = 0; u < B.N; u++)
        {

            //B.Groups[u]
            //random choose hyperedge
            if (B.Groups[u].size() == 0) continue;
            int v = (int)_rand(0.0, 1.0 * B.Groups[u].size());
            int random_edge = B.Groups[u][v];
            //int max_member = u;
            //double max_payoff = 0;
            //for (auto& v : B.Hyperedges[random_edge])
            //{
            //    if (max_payoff < Vectex_payoff_list[v])
            //    {
            //        max_member = v;
            //        max_payoff = Vectex_payoff_list[v];
            //    }
            //}
            //if (u == max_member) continue;
            int member;
            while (1)
            {
                v = (int)_rand(0, B.g);
                member = B.Hyperedges[random_edge][v];
                if (member != u)break;
            }
            double delta = r * (B.g - 2) + 1;
            if (r > 1)delta = r * B.g - 1;
            double w = (Vectex_payoff_list[member] - Vectex_payoff_list[u]) / delta;
            w = min(max(w, 0.0), 1.0);
            if (_rand(0, 1) < w)
                Vectex_status_list[u] = temp_Vectex_status_list[member];
        }
    }
    void Fermi_rule(HypergraphLib::Hypergraph B, double k = 0.1)
    {
        for (int u = 0; u < B.N; u++)
        {

            //B.Groups[u]
            //random choose hyperedge
            /*   1    1
            wij=---X-------------   k=1
                di  e^((fi-fj)/k)
                    1
            wij=-------------  k=0.04
                e^((fi-fj)/k)
            */
            double a = _rand(0, 1);
            if (a > p[u])continue;
            if (B.Groups[u].size() == 0) continue;
            int v = (int)_rand(0.0, 1.0 * B.Groups[u].size());
            int random_edge = B.Groups[u][v];
            int member;
            while (1)
            {
                v = (int)_rand(0, B.g);
                member = B.Hyperedges[random_edge][v];
                if (member != u)break;
            }
            double w = 1.0 / (1.0 + pow(2.71828, ((Vectex_payoff_list[u] - Vectex_payoff_list[member]) / k)));
            //w = w / B.Groups[u].size();
            w = min(max(w, 0.0), 1.0);
            if (_rand(0, 1) < w)
                Vectex_status_list[u] = temp_Vectex_status_list[member];
        }
    }
    void Moran_rule(HypergraphLib::Hypergraph B)
    {
        for (int u = 0; u < B.N; u++)
        {
            map<int, int>asd;
            set<int> neighbors;
            for (auto& v : B.Groups[u])
            {
                for (auto& vv : B.Hyperedges[v])
                {
                    neighbors.insert(vv);
                }
            }

            neighbors.erase(u);
            vector<double>p;
            vector<int>v;
            double sum_p = 0;
            set<int>::iterator it;//使用迭代器 
            for (it = neighbors.begin(); it != neighbors.end(); it++)
            {
                //cout << *it << ' ' << asd[*it] << endl;
                v.push_back(*it);
                p.push_back(Vectex_payoff_list[*it]);
                sum_p += (Vectex_payoff_list[*it]);
            }

            for (int i = 1; i < p.size(); i++)
                p[i] += p[i - 1];
            double r = _rand(0, sum_p);
            for (int i = 0; i < p.size(); i++)
            {

                if (p[i] >= r)
                {
                    Vectex_status_list[u] = temp_Vectex_status_list[v[i]];
                    //cout << 1 << endl;
                    break;
                }
            }
        }
    }
    void show_node_status(int _, HypergraphLib::Hypergraph B)
    {
        if (_ == 3 || _ == 10 || _ == 100 || _ == 1000 || _ == 10000)
        {
            ofstream outfile;
            outfile.open(to_string(_) + '_' + to_string(Hypergraph_assortativity(B)) + ".txt");
            for (int i = 0; i < B.N; i++)
                if (B.Groups[i].size() > 0)
                    outfile << i << ' ' << Vectex_status_list[i] << endl;
            outfile.close();
        }
        //int a = 0, a_c = 0;
        //int b = 0, b_c = 0;
        //int c = 0, c_c = 0;
        //for (int i = 0;i < B.N;i++)
        //{
        //    //if (B.Groups[i].size() > 0)cout << i << ' ' << PGG_avpi::Vectex_status_list[i] << ' ' << B.Groups[i].size() << endl;

        //    int di = B.Groups[i].size();
        //    if (di > 0) {
        //        if (di == 1)
        //        {
        //            a++;
        //            if (PGG_avpi::Vectex_status_list[i])a_c++;
        //        }
        //        else if (di <= 30)
        //        {
        //            b++;
        //            if (PGG_avpi::Vectex_status_list[i])b_c++;
        //        }
        //        else
        //        {
        //            c++;
        //            if (PGG_avpi::Vectex_status_list[i])c_c++;
        //        }
        //    }
        //}
    }

    double run(HypergraphLib::Hypergraph B, double r, int step_num, int init_step_num, int x, int timescale)
    {

        for (int i = 0; i < B.N; i++)
        {
            Vectex_status_list.push_back(0);
            temp_Vectex_status_list.push_back(0);
            Vectex_payoff_list.push_back(0);
            p.push_back(1);
        }
        for (int i = 0; i <= step_num; i++) Status_list.push_back(0);
        for (int i = 0; i < B.N; i++)
            if (_rand(0, 1) > 0.5)
            {
                Vectex_status_list[i] = 1;//1:cooperator,0:defector
                Status_list[1] += Vectex_status_list[i];
            }

        for (int _ = 2; _ <= step_num; _++)
        {

            //update payoff
            for (int i = 0; i < B.N; i++)Vectex_payoff_list[i] = 0;
            for (int l = 0; l < B.L; l++)
            {
                //B.Hyperedges[l]
                double R = r * B.g;
                //num of C
                int cooperator_num = 0;
                for (auto& v : B.Hyperedges[l]) cooperator_num += Vectex_status_list[v];
                //avpi
                double share = R * cooperator_num * c / (B.g * 1.0);
                for (auto& v : B.Hyperedges[l])
                {
                    if (Vectex_status_list[v] == 1)
                        Vectex_payoff_list[v] = Vectex_payoff_list[v] + (share - c) / (1.0 * B.Groups[v].size());
                    else
                        Vectex_payoff_list[v] = Vectex_payoff_list[v] + (share) / (1.0 * B.Groups[v].size());
                }
            }
            //update status
            for (int i = 0; i < B.N; i++)
                temp_Vectex_status_list[i] = Vectex_status_list[i];
            //timescale
            //if (timescale > 0)
            //    for (int i = 0; i < B.N; i++)
            //    {
            //        if (timescale == 1)p[i] = 1.0 / (1.0 + 1 * max(0.0, Vectex_payoff_list[i]));
            //        else if (timescale == 2)p[i] = 1.0 / (1.0 + B.Dergee_Centralites[i] * max(0.0, Vectex_payoff_list[i]));
            //        else if (timescale == 3)p[i] = 1.0 / (1.0 + B.Sub_Centralites[i] * max(0.0, Vectex_payoff_list[i]));
            //        else if (timescale == 4 || timescale == 5 || timescale == 6)p[i] = 1.0 / (1.0 + B.Important_Centralites[i] * max(0.0, Vectex_payoff_list[i]));
            //    }
            if (x == 1) Unconditional_Imitation(r, B);
            if (x == 2) Fermi_rule(B);
            if (x == 3) Moran_rule(B);
            if (x == 4)Betters_Possess_Chance(r, B);

            for (int i = 0; i < B.N; i++)
                if (!B.Groups[i].empty())
                    Status_list[_] += Vectex_status_list[i];
            //cout << _ << ' ' << Status_list[_] << endl;
        }
        for (int i = init_step_num; i <= step_num; i++)
        {
            fc += Status_list[i];
        }
        fc /= (1.0 * ((B.N - B.empty_nodes_num) * (step_num - init_step_num + 1)));

        return fc;
    }
    void clear()
    {
        /*
            int c = 1;//cost
    double fc = 0;
    vector<bool> Vectex_status_list;
    vector<bool> temp_Vectex_status_list;
    vector<double> Vectex_payoff_list;
    vector<int> Status_list;
    map<int, double>num_fc_k;
    vector<double>p;//timescale
        */
        fc = 0.0;
        Vectex_payoff_list.clear();
        temp_Vectex_status_list.clear();
        Vectex_status_list.clear();
        Status_list.clear();
        num_fc_k.clear();
        p.clear();
        vector<double>().swap(Vectex_payoff_list);
        vector<bool>().swap(temp_Vectex_status_list);
        vector<bool>().swap(Vectex_status_list);
        vector<int>().swap(Status_list);
        vector<double>().swap(p);
    }
}
#endif