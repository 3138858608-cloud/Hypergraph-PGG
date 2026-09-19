#pragma once
#ifndef HI_H
#define HI_H
#include "stdc++.h"
#include "rand.h"
#include<filesystem>
#include <iostream>
#include <fstream>
#include <regex>
#include<string>

using namespace std;
namespace HypergraphLib {
    //Define hypergraph structure
    struct Hypergraph {
        int N;//Number of nodes
        int g;// rank of the hyperedge
        int L;//size of hyperedge
        double mu;
        int empty_nodes_num = 0;
        vector<int> nodes;//nodes set
        vector< vector<int> > Hyperedges; //Hyperedge set
        vector< vector<int> > Groups;//Hyperedge containing the node set
        vector<pair<int, double>>KNN;
        vector<pair<int, double>>RK;
        //vector<double>Dergee_Centralites;
        //vector<double>Sub_Centralites;
        //vector<double>Important_Centralites;
        //超图清空
        void clear()
        {
            nodes.clear();
            Hyperedges.clear();
            Groups.clear();
            empty_nodes_num = 0;
            KNN.clear();
            RK.clear();
        }
        //超图构造
        void URH_make()
        {
            vector<double> c;
            for (int i = 1; i <= N - g + 1; i++) c.push_back(1 - pow((1.0 - i / (1.0 * (N - g + 1))), g));
            for (int i = 0; i < L; i++)
            {
                //step1: Random number
                double r = _rand(0, 1);

                //step2: find min i c_i>= r
                int idx = N - g;
                for (int ii = 0; ii < N - g + 1; ii++)
                {
                    if (c[ii] >= r)
                    {
                        idx = ii;
                        break;
                    }
                }
                //step3: from idx to N  find g-1
                Hyperedges[i] = make_sorted_edge(idx, N, g);
            }
        }
        void SRH_make()
        {
            /*
            Jhun, B., Minjae, J. & Kahng, B. Simplicial SIS model in scale-free uniform  hypergraph. J. Stat. Mech. https://doi.org/10.1088/1742-5468/ab5367 (2019).
            */
            vector<double>w;
            for (int i = 0; i < N; i++)w.push_back(0);
            double sum_w = 0;
            for (int i = 0; i < N; i++)
            {
                w[i] = pow(i + 1.0, -mu);
                if (i > 0)w[i] += w[i - 1];
            }
            for (int i = 0; i < N; i++) w[i] /= w[N - 1];
            for (int i = 0; i < L; i++)
            {
                while (1)
                {
                    set<int> edge;
                    while (edge.size() < g)
                    {
                        double r = _rand(0, 1);
                        for (int ii = 0; ii < N; ii++)
                        {
                            if (w[ii] >= r)
                            {
                                edge.insert(ii);
                                break;
                            }
                        }
                    }
                    vector<int> key(edge.begin(), edge.end());
                    sort(key.begin(), key.end());
                    if (find(Hyperedges.begin(), Hyperedges.end(), key) == Hyperedges.end())
                    {
                        for (auto& v : key)
                            Hyperedges[i].push_back(v);
                        break;
                    }
                }
            }

        }

        //重复节点检查
        bool Check_duplicate_nodes()
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 1; j < g; j++)
                {
                    if (Hyperedges[i][j] == Hyperedges[i][j - 1])
                    {
                        cout << Hyperedges[i][j] << ' ';
                        return true;
                    }
                }
            }
            return false;
        }
        //计算最大连通分量
        /*
        pair<int, vector<int>> largestConnectedComponentWithVertices(Mat& A)
        {

            vector<vector<int>>adj;
            for (int i = 0; i < A.row; i++)
            {
                adj.push_back(vector<int>{});
                for (int j = 0; j < A.col; j++)
                {
                    adj[i].push_back((int)A.val[i][j]);
                }
            }
            vector<bool> visited(N, false);
            int maxSize = 0;
            vector<int> maxComponent;
            auto dfs = [&](int start) -> vector<int> {
                vector<int> component;
                stack<int> st;
                st.push(start);

                while (!st.empty())
                {
                    int u = st.top();
                    st.pop();

                    if (!visited[u])
                    {
                        visited[u] = true;
                        component.push_back(u);

                        for (int v = 0; v < N; v++)
                        {
                            if (adj[u][v] == 1 && !visited[v])
                            {
                                st.push(v);
                            }
                        }
                    }
                }
                return component;
                };
            for (int i = 0; i < N; i++) {
                if (!visited[i]) {
                    vector<int> component = dfs(i);
                    if (component.size() > maxSize) {
                        maxSize = (int)component.size();
                        maxComponent = move(component);
                    }
                }
            }
            return { maxSize, maxComponent };
        }
        */
        //重连函数：用于随机重连超图结构
        void rewiring(int dv, int de)
        {
            if (dv == 1 && de == 1)
            {
                vector<int> node_stub;//000000001111122222
                vector<int> edge_stub;//000111222333

                for (int i = 0; i < N; i++)
                {
                    int ki = (int)Groups[i].size();
                    for (int t = 0; t < ki; t++)
                        node_stub.push_back(i);
                }
                for (int i = 0; i < L; i++)
                {
                    int sj = (int)Hyperedges[i].size();
                    for (int t = 0; t < sj; t++)
                        edge_stub.push_back(i);
                }

                mt19937 rng(random_device{}());
                shuffle(node_stub.begin(), node_stub.end(), rng);
                shuffle(edge_stub.begin(), edge_stub.end(), rng);
                vector<vector<int>> temp_Hyperedges(L);
                int m = (int)node_stub.size();
                for (int i = 0; i < m; i++)
                    temp_Hyperedges[edge_stub[i]].push_back(node_stub[i]);

                Hyperedges.swap(temp_Hyperedges);
                Groups.clear();
                Groups.assign(N, vector<int>());
                for (int i = 0; i < L; i++)
                    for (int v : Hyperedges[i])
                        Groups[v].push_back(i);
            }
        }
        //同配函数：用于构造同配/异配超图结构
        void xs(double p)
        {
            // 用于比较两个节点所属“组大小”的排序规则
            // 注意：this 捕获表明 Groups 是类成员变量
            auto cmp = [this](int a, int b) {
                return Groups[a].size() > Groups[b].size();
                };

            // ================= 重连过程 =================
            // 迭代次数与 p 和系统规模 L 成正比
            // p 越大 → 重连尝试越多
            for (int _ = 0; _ <= 500 * L * abs(p); _++)
            {
                // 概率触发机制：以 |p| 为概率执行一次重连尝试
                if (_rand(0, 1) > abs(p)) continue;

                // 随机选择两个超边
                int e1 = (int)_rand(0, L);
                int e2 = (int)_rand(0, L);

                // 避免选到同一个超边
                if (e1 == e2) continue;

                bool flag = 0;

                // ================= 判断两个超边是否有公共节点 =================
                // 若有交集则跳过（避免破坏结构或重复连接）
                for (int i = 0; i < g; i++)
                    for (int j = 0; j < g; j++)
                        if (Hyperedges[e1][i] == Hyperedges[e2][j])
                            flag = 1;

                if (flag) continue;

                // ================= 合并两个超边的节点集合 =================
                vector<int> v;
                for (int i = 0; i < g; i++)
                {
                    v.push_back(Hyperedges[e1][i]);
                    v.push_back(Hyperedges[e2][i]);
                }

                // 按节点“度/组大小”进行排序（用于同配/异配结构构造）
                sort(v.begin(), v.end(), cmp);

                // =========================================================
                // p > 0 ：同配（assortative mixing）
                // 高度节点倾向连接高度节点
                // =========================================================
                if (p > 0)
                {
                    vector<int> temp_Hyperedge1;
                    vector<int> temp_Hyperedge2;

                    // 前两个最大节点分配给 edge1
                    temp_Hyperedge1.push_back(v[0]);
                    temp_Hyperedge1.push_back(v[1]);

                    // 最小两个节点分配给 edge2
                    temp_Hyperedge2.push_back(v[2 * g - 2]);
                    temp_Hyperedge2.push_back(v[2 * g - 1]);

                    // 中间节点打乱后重新分配
                    vector<int> u;
                    for (int i = 2; i < 2 * g - 2; i++)
                        u.push_back(v[i]);

                    shuffle(u);

                    // 填充剩余节点
                    for (int i = 0; i < g - 2; i++)
                        temp_Hyperedge1.push_back(u[i]);

                    for (int i = g - 2; i < 2 * g - 4; i++)
                        temp_Hyperedge2.push_back(u[i]);

                    // 标准化（排序）用于判重
                    sort(temp_Hyperedge1.begin(), temp_Hyperedge1.end());

                    // 检查是否已存在相同超边
                    auto it1 = std::find(Hyperedges.begin(), Hyperedges.end(), temp_Hyperedge1);
                    if (it1 != Hyperedges.end()) continue;

                    sort(temp_Hyperedge2.begin(), temp_Hyperedge2.end());
                    auto it2 = std::find(Hyperedges.begin(), Hyperedges.end(), temp_Hyperedge2);
                    if (it2 != Hyperedges.end()) continue;

                    // 更新超边结构
                    Hyperedges[e1].clear();
                    Hyperedges[e2].clear();

                    Hyperedges[e1].assign(temp_Hyperedge1.begin(), temp_Hyperedge1.end());
                    Hyperedges[e2].assign(temp_Hyperedge2.begin(), temp_Hyperedge2.end());
                }

                // =========================================================
                // p < 0 ：异配（disassortative mixing）
                // 高度节点倾向连接低度节点
                // =========================================================
                else if (p < 0)
                {
                    vector<int> temp_Hyperedge1;
                    vector<int> temp_Hyperedge2;

                    // 交叉分配高低节点（异配结构）
                    temp_Hyperedge1.push_back(v[0]);
                    temp_Hyperedge2.push_back(v[1]);
                    temp_Hyperedge2.push_back(v[2 * g - 2]);
                    temp_Hyperedge1.push_back(v[2 * g - 1]);

                    // 中间节点打乱
                    vector<int> u;
                    for (int i = 2; i < 2 * g - 2; i++)
                        u.push_back(v[i]);

                    shuffle(u);

                    // 分配剩余节点
                    for (int i = 0; i < g - 2; i++)
                        temp_Hyperedge1.push_back(u[i]);

                    for (int i = g - 2; i < 2 * g - 4; i++)
                        temp_Hyperedge2.push_back(u[i]);

                    // 判重 + 标准化
                    sort(temp_Hyperedge1.begin(), temp_Hyperedge1.end());
                    auto it1 = std::find(Hyperedges.begin(), Hyperedges.end(), temp_Hyperedge1);
                    if (it1 != Hyperedges.end()) continue;

                    sort(temp_Hyperedge2.begin(), temp_Hyperedge2.end());
                    auto it2 = std::find(Hyperedges.begin(), Hyperedges.end(), temp_Hyperedge2);
                    if (it2 != Hyperedges.end()) continue;

                    // 更新结构
                    Hyperedges[e1].clear();
                    Hyperedges[e2].clear();

                    Hyperedges[e1].assign(temp_Hyperedge1.begin(), temp_Hyperedge1.end());
                    Hyperedges[e2].assign(temp_Hyperedge2.begin(), temp_Hyperedge2.end());
                }

                // ================= 更新节点-超边邻接表 =================
                // 清空旧结构
                for (int i = 0; i < N; i++)
                    Groups[i].clear();

                // 重新构建映射：节点 → 所属超边列表
                for (int i = 0; i < L; i++)
                    for (auto& v : Hyperedges[i])
                        Groups[v].push_back(i);
            }
        }

        /***************************************************************************************************************************************/
        //initialization function
        Hypergraph() :N(0), nodes(0), Groups(0), L(0), Hyperedges(0), g(0), mu(0) {}
        Hypergraph(int n, int L0, int Hs, double mu0) 
        {
            N = n;
            g = Hs;
            L = L0;
            mu = mu0;
            /*
                nodes has been initialized as a vector of size N, with element values of 0
                Hyperedges has been initialized as a vector of size L, each element being an empty vector
                x is map selector
                g is the rank of the hyperedge
            */
            while (1)
            {
                for (int i = 0; i < N; i++) nodes.push_back(i), Groups.push_back(vector<int>());
                for (int i = 0; i < L; i++) Hyperedges.push_back(vector<int>());
                cout << "Generating map..." << endl;
                SRH_make();
                for (int i = 0; i < L; i++)
                {
                    for (auto& v : Hyperedges[i])
                        Groups[v].push_back(i);
                }
                for (int i = 0; i < n; i++)
                {
                    if (Groups[i].empty())
                    {
                        empty_nodes_num++;
                    }
                }
                if (Check_duplicate_nodes())
                {
                    cout << "Duplicate nodes found, regenerating map..." << endl;
                    clear();
                }
                else
                    break; // No duplicates found, exit the loop

            }
            cout << "Map generation complete" << endl;
        }
    };
    //超图矩阵计算
    //Mat Incidence_matrix(Hypergraph B)
    //{
    //    /*
    //    关联矩阵是一个大小为∣V∣×∣E∣ 的矩阵，它的第 i 行、第 j 列的元素 1(vi∈ej)/0（other）
    //    */
    //    Mat res(B.N, B.L);
    //    for (int i = 0; i < B.N; i++)
    //        for (auto& v : B.Groups[i])
    //            res.val[i][v] = 1;
    //    return res;
    //}
    //Mat Degree_diagonal_matrix(Hypergraph B)
    //{
    //    /*
    //    超度对角矩阵是一个大小为∣V∣×∣V∣ 的矩阵，它的第 i 行、第 j 列的元素 节点超度(i==j)/0（other）
    //    */
    //    Mat res(B.N, B.N);
    //    for (int i = 0; i < B.N; i++)
    //    {
    //        res.val[i][i] = 1.0 * B.Groups[i].size();
    //    }
    //    return res;
    //}
    //Mat Adjacency_matrix(Hypergraph B)
    //{
    //    /*
    //    邻接矩阵A=I × IT −D
    //    */
    //    Mat res(B.N, B.N);
    //    Mat I = Incidence_matrix(B);
    //    Mat D = Degree_diagonal_matrix(B);
    //    res = M_Multiplication(I, M_T(I));
    //    res = M_Subtraction(res, D);
    //    return res;
    //}
    //计算指标
    double Hypergraph_assortativity(Hypergraph B)
    {
        /*
        Hypergraph assortativity: a dynamical systems perspective
           <k>^2 <kk1>E
        p=---------------- - 1
               <k^2>^2
        g=4
        */
        double avk = 0;
        double avkk = 0;
        for (int i = 0; i < B.N; i++)
        {
            double ki = 1.0 * B.Groups[i].size();
            avk += ki;
            avkk += (ki * ki);
        }
        avk /= (1.0 * B.N);
        avkk /= (1.0 * B.N);
        double kk1 = 0;
        for (int _ = 0; _ < B.L; _++)
        {
            //cout << _ << endl;
            for (int i = 0; i < B.g; i++)
                for (int j = i + 1; j < B.g; j++)
                {
                    double k = 1.0 * B.Groups[B.Hyperedges[_][i]].size();
                    double k1 = 1.0 * B.Groups[B.Hyperedges[_][j]].size();
                    kk1 += (k * k1);
                }
        }
        kk1 /= (B.L * 6.0);
        double ans = avk * avk * kk1;
        ans = ans / (avkk * avkk);
        ans--;
        return ans;
    }
    vector <pair< int, double >> knn(Hypergraph B)
    {
        vector <pair< int, double >> result;
        vector<vector<int>> mkk_matrix;
        for (int i = 0; i <= B.L; i++)
        {
            mkk_matrix.push_back({});
            for (int j = 0; j <= B.L; j++)
                mkk_matrix[i].push_back(0);
        }
        for (int i = 0; i < B.L; i++)
        {
            for (int j = 0; j < B.g; j++)
            {
                for (int jj = j + 1; jj < B.g; jj++)
                {
                    int k = (int)B.Groups[B.Hyperedges[i][j]].size();
                    int k_ = (int)B.Groups[B.Hyperedges[i][jj]].size();
                    //cout << k << ' ' << k_ << endl;
                    mkk_matrix[k][k_]++;
                    mkk_matrix[k_][k]++;
                    //cout << mkk_matrix[k][k_] << endl;
                }
            }
        }
        for (int k = 1; k <= B.L; k++)
        {
            int sum_mkk_ = 0;
            int k_sum_mkk_ = 0;
            for (int k_ = 1; k_ <= B.L; k_++)
                if (mkk_matrix[k][k_] > 0)
                {
                    //cout << k << ' ' << k_ << ' ' << mkk_matrix[k][k_] << ' ' ;
                    sum_mkk_ += mkk_matrix[k][k_];
                    k_sum_mkk_ += k_ * mkk_matrix[k][k_];
                    //cout << sum_mkk_ << ' ' << k_sum_mkk_ << endl;
                }
            if (sum_mkk_ > 0)
                result.push_back(make_pair(k, (1.0 * k_sum_mkk_) / ((1.0 * sum_mkk_))));
        }
        return result;
    }
    vector<pair<int, double>> rk(Hypergraph B)
    {
        vector<pair<int, double>> result;
        vector<double> ri(B.N, 0.0);
        vector<int> nk(B.L + 1, 0);
        for (int i = 0; i < B.N; i++)
        {
            int ki = (int)B.Groups[i].size();
            nk[ki]++;
            if (ki <= 1)
            {
                ri[i] = 0.0;
                continue;
            }

            int overlap_pairs = 0;
            // 枚举节点i所在超边对
            for (int a = 0; a < ki; a++)
            {
                int e1 = B.Groups[i][a];
                for (int b = 0; b < a; b++)
                {
                    int e2 = B.Groups[i][b];

                    bool overlap = false;
                    for (auto& u : B.Hyperedges[e1])
                    {
                        if (overlap)break;
                        for (auto& v : B.Hyperedges[e2])
                        {
                            if (overlap)break;
                            if (u == v && u != i)
                            {
                                overlap = true;
                                //cout << u << ' ' << v << endl;
                                break;
                            }
                        }
                    }

                    if (overlap)
                        overlap_pairs++;
                }
            }

            ri[i] =
                2.0 * overlap_pairs /
                (1.0 * ki * (ki - 1));
        }

        // 计算 r(k)
        for (int k = 1; k <= B.L; k++)
        {
            if (nk[k] == 0)
                continue;
            double rk_value = 0.0;
            for (int i = 0; i < B.N; i++)
            {
                if ((int)B.Groups[i].size() == k)
                    rk_value += ri[i];
            }
            rk_value /= nk[k];

            result.push_back({ k, rk_value });
        }

        return result;
    }
    //vector<double>dergee_centralites(Hypergraph B)
    //{
    //    vector<double>res;
    //    double sum = 0;
    //    for (int i = 0; i < B.N; i++)
    //    {
    //        res.push_back(B.Groups[i].size() * (B.g - 1.0));
    //        sum += res[i];
    //    }
    //    sum /= B.N;
    //    for (int i = 0; i < B.N; i++)
    //    {
    //        res[i] /= sum;
    //        //cout << i<<' '<<res[i] << '\n ';
    //    }
    //    //cout << endl;
    //    return res;
    //}
    //vector<double>sub_centralites(Hypergraph B)
    //{
    //    vector<double> C(B.N, 0.0);
    //    // =========================
    //    // 1. 构造邻接矩阵 A
    //    // =========================
    //    Mat A = Adjacency_matrix(B);
    //    int n = B.N;
    //    // =========================
    //    // 2. 转 Eigen Matrix
    //    // =========================
    //    Eigen::MatrixXd M(n, n);
    //    for (int i = 0; i < n; i++)
    //        for (int j = 0; j < n; j++)
    //            M(i, j) = A.val[i][j];
    //    // =========================
    //    // 3. 特征分解
    //    // =========================
    //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    //    if (solver.info() != Eigen::Success)
    //    {
    //        std::cout << "Eigen decomposition failed!" << std::endl;
    //        return C;
    //    }
    //    Eigen::VectorXd lambda = solver.eigenvalues();
    //    Eigen::MatrixXd U = solver.eigenvectors();
    //    // =========================
    //    // 4. exp(λ)
    //    // =========================
    //    Eigen::VectorXd exp_lambda = lambda.array().exp();
    //    // =========================
    //    // 5. 计算 CSH
    //    // =========================
    //    for (int i = 0; i < n; i++)
    //    {
    //        double sum = 0.0;
    //        for (int j = 0; j < n; j++)
    //        {
    //            double u = U(i, j);
    //            sum += u * u * exp_lambda(j);
    //        }
    //        C[i] = sum;
    //    }
    //    double sum = 0;
    //    for (int i = 0; i < n; i++)
    //        sum += C[i];
    //    sum /= n;
    //    for (int i = 0; i < n; i++)
    //    {
    //        C[i] /= sum;
    //        //cout << i << ' ' << C[i] << '\n ';
    //    }
    //    return C;
    //}
    //vector<double>important_centralites(Hypergraph B, int l)
    //{
    //    vector<double>CI;
    //    vector<int>nodes_with_d;
    //    for (int i = 0; i < B.N; i++)
    //    {
    //        int d = 0;
    //        deque<int>p;
    //        vector<bool>vis;
    //        vector<int>dis;
    //        for (int i = 0; i < B.N; i++)vis.push_back(0), dis.push_back(0);
    //        p.push_back(i);
    //        dis[i] = 0;
    //        vis[i] = 1;
    //        while (!p.empty())
    //        {
    //            int u = p.front();
    //            p.pop_front();
    //            set<int>neb_of_u;
    //            for (int j = 0; j < B.Groups[u].size(); j++)
    //                for (auto v : B.Hyperedges[B.Groups[u][j]])
    //                    if (!vis[v])
    //                        neb_of_u.insert(v);
    //            set<int>::iterator it;
    //            for (it = neb_of_u.begin(); it != neb_of_u.end(); it++)
    //            {
    //                int v = *it;
    //                if (vis[v])continue;
    //                dis[v] = dis[u] + 1;
    //                vis[v] = 1;
    //                p.push_back(v);
    //            }
    //        }
    //        for (int j = 0; j < B.N; j++)
    //            if (dis[j] == l)nodes_with_d.push_back(j);
    //        double sum = 0;
    //        CI.push_back(0);
    //        for (auto j : nodes_with_d)
    //            CI[i] = CI[i] + B.Groups[j].size() - 1;
    //        CI[i] = CI[i] * (B.Groups[i].size() - 1);
    //        nodes_with_d.clear();
    //    }
    //    double sum = 0;
    //    for (int i = 0; i < B.N; i++)
    //        sum += CI[i];
    //    sum /= B.N;
    //    for (int i = 0; i < B.N; i++)
    //        CI[i] /= sum;
    //    return CI;
    //}
    void Hypergraph_file_generation(int n, int G, double MU, string path)
    {
        // ================= 基本参数初始化 =================
        int N = n;              // 节点数
        int g = G;              // 每条超边包含的节点数
        double m = 1 * log(N);  // 控制超边数量的尺度参数
        int L = (int)(N * m) / g;    // 超边数量（由 N、g 决定）
        double mu = MU;         // 超图生成的结构参数

        // ===================================================
        // ===================== 基准图 p = 0 =================
        // ===================================================
        double p0 = 0;
        // 构造初始随机超图（无同配偏好）
        Hypergraph hypergraph0(N, L, g,mu);
        // 计算该超图的同配系数（assortativity）
        double uk = Hypergraph_assortativity(hypergraph0);
        // 输出文件名（p=0 基准图）
        string fileadd = path + "//Hypergraph_N" + to_string(N) +
            "_g" + to_string(g) +
            "_mu" + to_string(mu) +
            "_p" + to_string(p0) + ".txt";

        ofstream outfile;
        outfile.open(fileadd);

        // 写入同配系数
        outfile << uk << endl;
        // 写入超图结构
        for (int i = 0; i < hypergraph0.L; i++)
        {
            outfile << 1 << ' '; // 超边权重
            for (auto& v : hypergraph0.Hyperedges[i])
            {
                outfile << v << ' '; // 输出超边中的节点
            }
            outfile << endl;
        }
        outfile.close();


        // ===================================================
        // ===================== 同配（p > 0）=================
        // ===================================================
        for (double p = 0.01; p <= 0.1; p += 0.01)
        {
            cout << p << endl;
            HypergraphLib::Hypergraph hypergraph;
            hypergraph.clear();
            double p1 = p0;
            while (1)
            {
                HypergraphLib::Hypergraph hypergraph1;
                hypergraph1 = hypergraph0;
                hypergraph1.xs(p);
                p1 = HypergraphLib::Hypergraph_assortativity(hypergraph1);
                if (1 > p1 && p1 > p0)
                {
                    p0 = p1;
                    hypergraph = hypergraph1;
                    break;
                }
            }

            // 输出文件名
            string fileadd = path + "//Hypergraph_N" + to_string(N) +
                "_g" + to_string(g) +
                "_mu" + to_string(mu) +
                "_p" + to_string(p) + ".txt";
            ofstream outfile;
            outfile.open(fileadd);
            // 写入同配系数
            outfile << p1 << '\n';
            // 写入超图结构
            for (int i = 0; i < hypergraph.L; i++)
            {
                //超边权重
                outfile << 1 << ' ';
                for (auto& v : hypergraph.Hyperedges[i])
                {
                    outfile << v << ' ';
                }
                outfile << endl;
            }
            outfile.close();
        }

        // ===================================================
        // ===================== 异配（p < 0）=================
        // ===================================================
        p0 = 0;
        for (double p = -0.01; p >= -0.1; p -= 0.01)
        {
            cout << p << endl;
            HypergraphLib::Hypergraph hypergraph;
            // 计算同配性
            double p1 = p0;
            while (1)
            {
                HypergraphLib::Hypergraph hypergraph1;
                hypergraph1 = hypergraph0;
                hypergraph1.xs(p);
                p1 = HypergraphLib::Hypergraph_assortativity(hypergraph1);
                if (p1 < p0)
                {
                    p0 = p1;
                    hypergraph = hypergraph1;
                    break;
                }
            }
            // 输出文件名
            string fileadd = path + "//Hypergraph_N" + to_string(N) +
                "_g" + to_string(g) +
                "_mu" + to_string(mu) +
                "_p" + to_string(p) + ".txt";
            ofstream outfile;
            outfile.open(fileadd);
            // 写入同配系数
            outfile << p1 << '\n';
            // 写入超图结构
            for (int i = 0; i < hypergraph.L; i++)
            {
                //超边权重
                outfile << 1 << ' ';
                for (auto& v : hypergraph.Hyperedges[i])
                {
                    outfile << v << ' ';
                }
                outfile << endl;
            }
            outfile.close();
        }
    }
    Hypergraph create_graph_from_file(filesystem::directory_entry entry)
    {
        Hypergraph hypergraph;
        hypergraph.clear();
        string fileadd = entry.path().string();
        ifstream in(entry.path());
        string fname = entry.path().filename().string();
        // 解析文件名：Hypergraph_N<NUM>_g<NUM>_mu<NUM>_p<NUM>.txt
        std::regex pattern0("Hypergraph_N(\\d+)_g(\\d+)_mu([0-9.]+)_p([0-9.]+)\\.txt", std::regex::icase);
        std::regex pattern1("Hypergraph_N(\\d+)_g(\\d+)_mu([0-9.]+)_p-([0-9.]+)\\.txt", std::regex::icase);
        std::smatch m;
        int N_file = {};
        int g_file = {};
        double mu_file = {};
        double p_file = {};
        if (std::regex_match(fname, m, pattern0)) {
            N_file = std::stoi(m[1].str());
            g_file = std::stoi(m[2].str());
            mu_file = std::stod(m[3].str());
            p_file = std::stod(m[4].str());
            std::cout << fileadd << "  -> N=" << N_file << " g=" << g_file << " mu=" << mu_file << " p=" << p_file << '\n';
        }
        if (std::regex_match(fname, m, pattern1)) {
            N_file = std::stoi(m[1].str());
            g_file = std::stoi(m[2].str());
            mu_file = std::stod(m[3].str());
            p_file = std::stod(m[4].str());
            std::cout << fileadd << "  -> N=" << N_file << " g=" << g_file << " mu=" << mu_file << " p=-" << p_file << '\n';
        }
        hypergraph.N = N_file;//Number of nodes
        hypergraph.g = g_file;// rank of the hyperedge
        hypergraph.L = (int)(hypergraph.N * 1.0 * log(hypergraph.N)) / hypergraph.g;//size of hyperedge
        hypergraph.mu = mu_file;
        for (int i = 0; i < hypergraph.N; i++)
            hypergraph.nodes.push_back(i), hypergraph.Groups.push_back(vector<int>());
        for (int i = 0; i < hypergraph.L; i++)
            hypergraph.Hyperedges.push_back(vector<int>());
        string line;
        bool flag = 0;
        int l = 0;
        bool duplicate_node_error = false;
        while (getline(in, line)) {
            if (!flag)
            {
                //超边权重
                flag = 1;
                continue;
            }
            //cout << line << std::endl;
            stringstream ss(line);
            int x;
            int x2 = -1;
            bool flag1 = 0;
            while (ss >> x) {
                if (!flag1)
                {

                    flag1 = 1;
                    continue;
                }
                if (x == x2)
                {
                    duplicate_node_error = 1;
                }
                hypergraph.Hyperedges[l].push_back(x);
                x2 = x;
            }
            l++;
        }
        for (int i = 0; i < hypergraph.L; i++)
        {
            for (auto& v : hypergraph.Hyperedges[i])
                hypergraph.Groups[v].push_back(i);
        }
        //empty_nodes_num
        for (int i = 0; i < hypergraph.N; i++)
        {
            if (hypergraph.Groups[i].empty())
            {
                hypergraph.empty_nodes_num++;
            }
        }
        if (duplicate_node_error)
        {
            cout << "Warning: Duplicate nodes found in hyperedges. Please check the input file." << endl;
        }
        //cout << "Hypergraph_assortativity:" << hypergraph.Hypergraph_assortativity() << endl;
        return hypergraph;
    }
    Hypergraph create_real_graph_from_file(filesystem::directory_entry entry)
    {
        Hypergraph hypergraph;
        hypergraph.clear();
        string fileadd = entry.path().string();
        ifstream in(entry.path());
        string fname = entry.path().filename().string();
        string line;
        bool flag = 0;
        int l = 0;
        bool duplicate_node_error = false;
        while (getline(in, line)) {
            if (!flag)
            {
                //超边权重
                flag = 1;
                stringstream ss(line);
                int N;
                int L;
                int g;
                ss >> N >> L >> g;
                hypergraph.N = N;
                hypergraph.L = L;
                hypergraph.g = g;
                cout << N << ' ' << L << ' ' << g << endl;
                for (int i = 0; i < hypergraph.N; i++)
                    hypergraph.nodes.push_back(i), hypergraph.Groups.push_back(vector<int>());
                for (int i = 0; i < hypergraph.L; i++)
                    hypergraph.Hyperedges.push_back(vector<int>());
                continue;
            }
            //cout << line << std::endl;
            stringstream ss(line);
            int x;
            int x2 = -1;
            while (ss >> x) {
                if (x == x2)
                {
                    duplicate_node_error = 1;
                }
                hypergraph.Hyperedges[l].push_back(x);
                x2 = x;
            }
            l++;
        }
        for (int i = 0; i < hypergraph.L; i++)
        {
            for (auto& v : hypergraph.Hyperedges[i])
                hypergraph.Groups[v].push_back(i);
        }
        //empty_nodes_num
        for (int i = 0; i < hypergraph.N; i++)
        {
            if (hypergraph.Groups[i].empty())
            {
                hypergraph.empty_nodes_num++;
            }
        }
        if (duplicate_node_error)
        {
            cout << "Warning: Duplicate nodes found in hyperedges. Please check the input file." << endl;
        }
        //cout << "Hypergraph_assortativity:" << hypergraph.Hypergraph_assortativity() << endl;
        return hypergraph;
    }
}

#endif#pragma once
