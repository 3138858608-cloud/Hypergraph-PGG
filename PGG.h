#pragma once
#ifndef AVPGG_H
#define AVPGG_H

#include"stdc++.h"
#include"HI.h"
#include"rand.h"
#include <algorithm>
using namespace std;
#include <iostream>
#include <algorithm>
using namespace std;

void showProgress(size_t current, size_t total)
{
    if (total == 0) return;

    const int width = 50;
    double progress = static_cast<double>(current) / total;
    progress = min(progress, 1.0);

    int filled = static_cast<int>(progress * width);

    cout << "\r[";

    for (int i = 0; i < width; ++i)
        cout << (i < filled ? '#' : '-');

    cout << "] "
        << static_cast<int>(progress * 100) << "%"
        << " (" << current << "/" << total << ")"
        << flush;

    if (current >= total)
        cout << '\n';
}

namespace PGG_avpi {

    int c = 1;//cost
    double fc = 0;
    vector<bool> Vectex_status_list;
    vector<bool> temp_Vectex_status_list;
    vector<double> Vectex_payoff_list;
    vector<int> Status_list;
    map<int, double>num_fc_k;
    vector<double>p;//timescale


    /*
    合作簇
    */
    struct UF {
        vector<int> parent, size;
        UF(int n = 0) { init(n); }
        void init(int n) {
            parent.resize(n);
            size.assign(n, 1);
            for (int i = 0; i < n; ++i) parent[i] = i;
        }
        int find(int x) {
            if (parent[x] != x) parent[x] = find(parent[x]);
            return parent[x];
        }
        void unite(int x, int y) {
            int rx = find(x), ry = find(y);
            if (rx == ry) return;
            if (size[rx] < size[ry]) swap(rx, ry);
            parent[ry] = rx;
            size[rx] += size[ry];
        }
        int get_size(int x) { return size[find(x)]; }
    };
    ofstream fout_summary;
    ofstream fout_cluster;
    // 初始化输出文件
    void InitClusterAnalysis(const string& cin_path, double r, double p, int graph_id, int run_id)
    {
        // 构建子文件夹名：R=xx_p=yy
        char subfolder[128];
        sprintf(subfolder, "R=%.2f_p=%.2f", r*4, p);

        // 完整路径：cin_path/CR/R=xx_p=yy/
        string cr_dir = cin_path + "/ClusterAnalysisResults/" + string(subfolder);
        filesystem::create_directories(cr_dir);

        char filename[512];
        sprintf(filename, "%s/summary_g%02d_r%02d.txt", cr_dir.c_str(), graph_id, run_id);
        fout_summary.open(filename);

        sprintf(filename, "%s/cluster_g%02d_r%02d.txt", cr_dir.c_str(), graph_id, run_id);
        fout_cluster.open(filename);

        // 表头
        fout_summary << "t\tS_C\tmax_ratio\tavg_payoff_cluster\tavg_payoff_boundary\tnum_clusters\n";
        fout_cluster << "t\tsize_sequence\n";
    }
    void AnalyzeCooperativeClusters(int t, const HypergraphLib::Hypergraph& B)
    {
        // ---------- 1. 只用合作者构建并查集 ----------
        UF uf(B.N);

        for (int e = 0; e < B.L; ++e)
        {
            // 收集当前超边中的所有合作者
            vector<int> cooperators;
            for (int v : B.Hyperedges[e])
            {
                if (Vectex_status_list[v] == 1)
                    cooperators.push_back(v);
            }
            // 把同一超边内的合作者两两连接（实际只需链式连接即可）
            for (size_t i = 1; i < cooperators.size(); ++i)
                uf.unite(cooperators[0], cooperators[i]);
        }

        // ---------- 2. 收集所有合作簇的大小 ----------
        vector<int> cluster_sizes;
        vector<char> visited(B.N, 0);
        int max_root = -1;
        int curr_max_size = 0;

        for (int i = 0; i < B.N; ++i)
        {
            if (Vectex_status_list[i] == 1)
            {
                int root = uf.find(i);
                if (!visited[root])
                {
                    visited[root] = 1;
                    int sz = uf.get_size(root);
                    cluster_sizes.push_back(sz);
                    if (sz > curr_max_size)
                    {
                        curr_max_size = sz;
                        max_root = root;
                    }
                }
            }
        }
        sort(cluster_sizes.begin(), cluster_sizes.end(), greater<int>());

        double S_C = (B.N > 0) ? (double)curr_max_size / B.N : 0.0;   // 最大簇占全网络比例
        double C = 0.0;
        for (int i = 0; i < B.N; ++i) if (Vectex_status_list[i] == 1) C += 1.0;
        double max_ratio = (C > 0) ? (double)curr_max_size / C : 0.0;  // 最大簇占所有合作者比例

        // ---------- 3. 标记属于最大合作簇的节点 ----------
        vector<char> in_max_cluster(B.N, 0);
        if (max_root != -1)
        {
            for (int i = 0; i < B.N; ++i)
            {
                if (Vectex_status_list[i] == 1 && uf.find(i) == max_root)
                    in_max_cluster[i] = 1;
            }
        }

        // ---------- 4. 找出与最大合作簇相邻的边界背叛者 ----------
        vector<char> is_boundary_D(B.N, 0);
        for (int e = 0; e < B.L; ++e)
        {
            bool has_max_C = false, has_D = false;
            for (int v : B.Hyperedges[e])
            {
                if (in_max_cluster[v]) has_max_C = true;
                if (Vectex_status_list[v] == 0) has_D = true;
            }
            if (has_max_C && has_D)
            {
                for (int v : B.Hyperedges[e])
                {
                    if (Vectex_status_list[v] == 0)
                        is_boundary_D[v] = 1;
                }
            }
        }

        // ---------- 5. 计算最大合作簇内收益 & 边界背叛者收益 ----------
        double sum_payoff_cluster = 0.0;
        int cnt_cluster = 0;
        double sum_payoff_boundary = 0.0;
        int cnt_boundary = 0;

        for (int i = 0; i < B.N; ++i)
        {
            if (in_max_cluster[i])
            {
                sum_payoff_cluster += Vectex_payoff_list[i];
                cnt_cluster++;
            }
            if (is_boundary_D[i])
            {
                sum_payoff_boundary += Vectex_payoff_list[i];
                cnt_boundary++;
            }
        }

        double avg_payoff_cluster = (cnt_cluster > 0) ? sum_payoff_cluster / cnt_cluster : 0.0;
        double avg_payoff_boundary = (cnt_boundary > 0) ? sum_payoff_boundary / cnt_boundary : 0.0;

        // ---------- 6. 输出 ----------
        fout_summary << t << "\t"
            << S_C << "\t"
            << max_ratio << "\t"
            << avg_payoff_cluster << "\t"
            << avg_payoff_boundary << "\t"
            << cluster_sizes.size() << "\n";

        // 大小分布（只输出 size > 1 的簇，从大到小）
        fout_cluster << t << "\t";
        for (int sz : cluster_sizes)
        {
            if (sz > 1) fout_cluster << sz << " ";
        }
        fout_cluster << "\n";
    }
    void ClusterAnalysisEnd()
    {
        if (fout_summary.is_open()) fout_summary.close();
        if (fout_cluster.is_open()) fout_cluster.close();
    }
    /*
   策略更新规则：
   */
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
    double run(HypergraphLib::Hypergraph B, double r, int step_num, int init_step_num, int x, int timescale,
                bool enable_cluster = false, const string& cin_path = "", int graph_id = 0, int run_id = 0)
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

        if (enable_cluster)
                InitClusterAnalysis(cin_path, r, B.p, graph_id, run_id);
        for (int _ = 2; _ <= step_num; _++)
        {
            if (enable_cluster && _ <= 1000)
            {
                AnalyzeCooperativeClusters(_, B);
            }
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
            if (timescale > 0)
                for (int i = 0; i < B.N; i++)
                {
                    if (timescale == 1)p[i] = 1.0 / (1.0 + 1 * max(0.0, Vectex_payoff_list[i]));
                    else if (timescale == 2)p[i] = 1.0 / (1.0 + B.Dergee_Centralites[i] * max(0.0, Vectex_payoff_list[i]));
                    else if (timescale == 3)p[i] = 1.0 / (1.0 + B.Sub_Centralites[i] * max(0.0, Vectex_payoff_list[i]));
                    else if (timescale == 4)p[i] = 1.0 / (1.0 + B.important_centralites[i] * max(0.0, Vectex_payoff_list[i]));
                }
            if (x == 1) Unconditional_Imitation(r, B);
            if (x == 2) Fermi_rule(B);
            if (x == 3) Moran_rule(B);
            if (x == 4)Betters_Possess_Chance(r, B);

            for (int i = 0; i < B.N; i++)
                if (!B.Groups[i].empty())
                    Status_list[_] += Vectex_status_list[i];
            showProgress(_, step_num);
        }
        for (int i = init_step_num; i <= step_num; i++)
        {
            fc += Status_list[i];
        }
        if (enable_cluster)
            ClusterAnalysisEnd();
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

    // 返回 data 中前10大的元素的原始下标
    vector<int> getTop10Index(const vector<double>& data)
    {
        vector<int> idx(data.size());
        iota(idx.begin(), idx.end(), 0);

        // 先找前100大
        int n100 = min(100, (int)idx.size());

        partial_sort(
            idx.begin(),
            idx.begin() + n100,
            idx.end(),
            [&](int a, int b)
            {
                return data[a] > data[b];
            }
        );

        idx.resize(n100);

        // 前100个本身已经按照从大到小排列
        // 所以前10个就是前10大
        int n10 = min(10, (int)idx.size());

        idx.resize(n10);

        return idx;
    }

}
#endif