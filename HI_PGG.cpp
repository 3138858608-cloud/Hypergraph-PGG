#include <iostream>
#include "stdc++.h"
#include "HI.h"
#include "PGG.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
// ================= 基本参数设置 =================
int N = 10000;          // 网络规模（节点数）
int g = 2;             // 超边大小（后面会在循环中修改）
double m = 1 * log(N);// 用于计算超边数量的系数
int L = (int)(N * m) / g + 1; // 超边数量（由N和g决定）

double mu = 1;         // 参数mu（后面循环修改）
int step_num = 1000000;  // 主演化步数
int init_step_num = 20000; // 初始热身步数
int retime = 100;        // 重复实验次数

int x = 2;  // 网络类型选择：1=URH，2=SRH
int y = 2;  // 更新规则选择：1=UI，2=feimi，3=MOR，4=BPC
int timescale = 0; // 时间尺度选择：0=无时间尺度，1=有时间尺度(n=1), 2 = 有时间尺度(节点度中心性), 3 = 有时间尺度(节点子超图中心性), 4 = 有时间尺度(节点集体影响力)

// 存储不同 r 下的最终合作水平 Fc
vector<double> Fc;
vector<double> r;


void Generate_Hypergraphfiles(int N0,int g0,double mu0,string cin_path0)
{
    //HypergraphLib::Hypergraph_file_generation(1000, 4, 0.25, cin_path);
    HypergraphLib::Hypergraph_file_generation(N0, g0, mu0, cin_path0);
    HypergraphLib::Hypergraph hypergraph;
    hypergraph.clear();
    filesystem::path dir_path(cin_path0);
    for (const auto& entry : filesystem::directory_iterator(dir_path))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".txt")
        {
            cout << "Processing file: " << entry.path() << endl;
            hypergraph = HypergraphLib::create_graph_from_file(entry);
        }
        hypergraph.clear();
    }
}
//void run_by_files()
//{
//    //使用前请通过Generate_Hypergraphfiles()生成超图文件
//    HypergraphLib::Hypergraph hypergraph;
//    hypergraph.clear();
//    filesystem::path dir_path(cin_path);
//    ofstream outfile;
//    outfile.open(cout_path);
//    outfile << cin_path << endl;
//    for (const auto& entry : filesystem::directory_iterator(dir_path))
//    {
//        if (entry.is_regular_file() && entry.path().extension() == ".txt")
//        {
//            cout << "Processing file: " << entry.path() << endl;
//            hypergraph = HypergraphLib::create_graph_from_file(entry);
//            hypergraph.KNN = HypergraphLib::knn(hypergraph);
//            double p = HypergraphLib::Hypergraph_assortativity(hypergraph);
//            string fname = entry.path().filename().string();
//            outfile << fname << endl;
//            outfile << p << endl;
//            // 计算超图指标
//            //hypergraph.KNN = HypergraphLib::knn(hypergraph);
//            //hypergraph.RK = HypergraphLib::rk(hypergraph);
//            for (timescale = 0; timescale <= 0; timescale++)
//            {
//                //if (timescale == 2)hypergraph.Dergee_Centralites = HypergraphLib::dergee_centralites(hypergraph);
//                //if (timescale == 3)hypergraph.Sub_Centralites = HypergraphLib::sub_centralites(hypergraph);
//                //if (timescale == 4)hypergraph.Important_Centralites = HypergraphLib::important_centralites(hypergraph, 1);
//                //if (timescale == 5)hypergraph.Important_Centralites = HypergraphLib::important_centralites(hypergraph, 2);
//                //if (timescale == 6)hypergraph.Important_Centralites = HypergraphLib::important_centralites(hypergraph, 3);
//                r.clear();
//                Fc.clear();
//                // ================= 构造 r 参数扫描序列 =================
//                for (double i = 0.8; i <= 0.7; i += 0.0125)
//                {
//                    r.push_back(i);
//                    Fc.push_back(0);
//                }
//                // ================= 重复实验 =================
//                for (int ii = 0; ii < r.size(); ii++)
//                {
//                    // 遍历所有 r 值，逐个进行仿真
//                    for (int i = 0; i < retime; i++)
//                    {
//                        // 运行演化过程
//                        PGG_avpi::run(hypergraph, r[ii], step_num, init_step_num, y, timescale);
//                        // 输出当前 r 对应的合作水平 fc
//                        cout << r[ii] << ' ' << PGG_avpi::fc << endl;
//                        // 累加结果（用于平均）
//                        Fc[ii] += PGG_avpi::fc;
//                        // 清理本轮仿真状态
//                        PGG_avpi::clear();
//                    }
//                }
//                // ================= 结果归一化 + 输出 =================
//                for (int ii = 0; ii < Fc.size(); ii++)
//                {
//                    // 计算平均合作水平
//                    Fc[ii] /= (1.0 * retime);
//                    cout << r[ii] << ' ' << Fc[ii] << endl;
//                    outfile << r[ii] << ' ' << Fc[ii] << endl;
//                }
//            }
//            // 清理当前超图
//            hypergraph.clear();
//        }
//    }
//    outfile.close();
//}

// ============================================================
// 参数结构体
// ============================================================
struct SimulationParameters
{
    // ==========================================
    // 超图参数
    // ==========================================
    int N = 1000;
    int g = 4;
    double mu = 0.5;
    // 超图输入根目录
    string cin_path;

    // ==========================================
    // 仿真参数
    // ==========================================
	int x = 2;  // 输出内容选择：1=超图指标，2=PGG仿真
    int step_num = 20000;
    int init_step_num = 15000;
    int retime = 1;
    double c = 1.0;
    int update_rule = 1;
    string cout_path;
};

// ============================================================
// 输入超图生成参数
// ============================================================
void InputHypergraphParameters(SimulationParameters& params)
{
    cout << endl;
    cout << "========== 超图生成参数 ==========" << endl;
    cout << "输入节点数 N: ";
    cin >> params.N;
    cout << "输入超边阶数 g: ";
    cin >> params.g;
    cout << "输入度分布参数 mu: ";
    cin >> params.mu;
    cout << "输入超图存储地址: ";
    cin >> params.cin_path;
}

// ============================================================
// 输入PGG仿真参数
// ============================================================
void InputSimulationParameters(SimulationParameters& params)
{
    cout << endl;
    cout << "==============================================" << endl;
    cout << "                 PGG仿真" << endl;
    cout << "==============================================" << endl;

    // --------------------------------------------------------
    // 输入超图地址
    // --------------------------------------------------------
    cout << "输入超图地址: ";
    cin >> params.cin_path;


    // --------------------------------------------------------
    // 选择运行内容
    // --------------------------------------------------------
    int mode;

    cout << endl;
    cout << "选择运行内容：" << endl;
    cout << "1. 超图指标" << endl;
    cout << "2. PGG仿真(使用1生成的超图)" << endl;

    cout << "输入选项 (1/2): ";
    cin >> mode;


    // --------------------------------------------------------
    // 检查输入
    // --------------------------------------------------------
    while (mode != 1 && mode != 2)
    {
        cout << "输入错误，请重新输入 (1/2): ";
        cin >> mode;
    }


    // --------------------------------------------------------
    // 超图指标
    // --------------------------------------------------------
    if (mode == 1)
    {
        params.x = 1;
        cout << endl;
    }

    // --------------------------------------------------------
    // PGG仿真
    // --------------------------------------------------------
    else
    {
        params.x = 2;
        cout << endl;
    }
}
/* run this program using the console pauser or add your own getch, system("pause") or input loop */
int main(int argc, char** argv) {

    SimulationParameters params;
    int mode;
    // ========================================================
    // 主菜单
    // ========================================================
    cout << endl;
    cout << "==============================================" << endl;
    cout << "        Hypergraph Public Goods Game" << endl;
    cout << "==============================================" << endl;

    cout << endl;
    cout << "What do you want?" << endl;
    cout << endl;

    cout << "1. 超图生成" << endl;
    cout << "2. PGG 仿真" << endl;

    cout << endl;
    cout << "输入选项 (1/2): ";
    cin >> mode;
    // ========================================================
    // 模式 1：超图生成
    // ========================================================
    if (mode == 1)
    {
        cout << endl;
        cout << "==============================================" << endl;
        cout << "                 1. 超图生成" << endl;
        cout << "==============================================" << endl;
        cout << "会在cin_path中建立20个文件夹，每个文件夹里有一个超图以及这个超图同配、异配操作后的超图\n";
        // 输入超图参数
        InputHypergraphParameters(params);


        // 参数确认
        cout << endl;
        cout << "========== 参数确认 ==========" << endl;

        cout << "N  = " << params.N << endl;
        cout << "g  = " << params.g << endl;
        cout << "mu = " << params.mu << endl;
        cout << "cin_path = " << params.cin_path << endl;
        // 开始生成
        char choice;
        cout << endl;
        cout << "是否开始生成超图？(y/n): ";
        cin >> choice;
        if (choice != 'y' && choice != 'Y')
        {
            cout << "已取消。" << endl;
            return 0;
        }
        cout << endl;
        cout << "生成超图..." << endl;
        // ====================================================
        // 超图生成
        // ====================================================
        for (int i = 1; i <= 10; i++)
        {
            string folder_path = params.cin_path + "/Hypergraph_" + to_string(i);
            filesystem::create_directories(filesystem::path(folder_path));
            Generate_Hypergraphfiles(params.N, params.g, params.mu, folder_path);
        }
        cout << "超图生成完成。" << endl;
        return 0;
    }
    // ========================================================
    // 模式 2：PGG 仿真
    // ========================================================
    else if (mode == 2)
    {
        cout << endl;
        cout << "==============================================" << endl;
        cout << "                  2. PGG 仿真" << endl;
        cout << "==============================================" << endl;
        // 输入仿真参数
        InputSimulationParameters(params);
        if (params.x == 1)
        {
            // ====================================================
            // 创建结果根目录
            // ====================================================
            filesystem::path result_root =
                filesystem::path(params.cin_path) / "IndicatorResults";
            filesystem::create_directories(result_root);
            params.cout_path = result_root.string();

            // ====================================================
            // 创建10个对应的结果文件夹
            // ====================================================
            for (int i = 1; i <= 10; i++)
            {
                filesystem::path result_folder =
                    result_root / ("Hypergraph_" + to_string(i));

                filesystem::create_directories(result_folder);

                cout << "创建结果目录: "
                    << result_folder << endl;
            }

            cout << endl;
            cout << "超图指标结果目录创建完成。" << endl;

            // ====================================================
            // 开始计算超图指标
            // ====================================================
            cout << endl;
            cout << "超图指标计算。。。" << endl;

            for (int i = 1; i <= 10; i++)
            {
                // ------------------------------------------------
                // 当前输入文件夹
                // ------------------------------------------------
                filesystem::path folder_cin_path =
                    filesystem::path(params.cin_path)
                    / ("Hypergraph_" + to_string(i));
                // ------------------------------------------------
                // 当前输出文件夹
                // ------------------------------------------------
                filesystem::path folder_cout_path =
                    filesystem::path(params.cout_path)
                    / ("Hypergraph_" + to_string(i));
                cout << endl;
                cout << "==============================================" << endl;
                cout << "处理 Hypergraph_" << i << endl;
                cout << "输入目录: " << folder_cin_path << endl;
                cout << "输出目录: " << folder_cout_path << endl;
                cout << "==============================================" << endl;
                // ------------------------------------------------
                // 检查输入目录
                // ------------------------------------------------
                if (!filesystem::exists(folder_cin_path))
                {
                    cout << "输入目录不存在，跳过。" << endl;
                    continue;
                }

                // ------------------------------------------------
                // 遍历当前文件夹中的所有 txt
                // ------------------------------------------------
                for (const auto& entry :
                    filesystem::directory_iterator(folder_cin_path))
                {
                    if (!entry.is_regular_file())
                        continue;

                    if (entry.path().extension() != ".txt")
                        continue;

                    cout << endl;
                    cout << "Processing file: "
                        << entry.path() << endl;

                    // ====================================================
                    // 1. 读取超图
                    // ====================================================
                    HypergraphLib::Hypergraph hypergraph;

                    hypergraph =
                        HypergraphLib::create_graph_from_file(entry);

                    // ====================================================
                    // 2. 计算超图指标
                    // ====================================================
                    hypergraph.KNN = HypergraphLib::knn(hypergraph);
                    hypergraph.RK = HypergraphLib::rk(hypergraph);
                    // ====================================================
                    // 3. 创建输出文件名
                    // ====================================================
                    string input_filename =
                        entry.path().stem().string();

                    filesystem::path output_file =
                        folder_cout_path /
                        (input_filename + "_indicator.txt");

                    // ====================================================
                    // 4. 输出指标
                    // ====================================================
                    ofstream fout(output_file);
                    for (int i = 0; i < hypergraph.KNN.size(); i++)
                    {
                        fout << hypergraph.KNN[i].first << ' ' << hypergraph.KNN[i].second << ' ' << hypergraph.RK[i].second << endl;
                        cout << hypergraph.KNN[i].first << ' ' << hypergraph.KNN[i].second << ' ' << hypergraph.RK[i].second << endl;
                    }
                    fout.close();
                    cout << "指标输出完成: " << output_file << endl;
                    // ====================================================
                    // 5. 清空当前超图
                    // ====================================================
                    hypergraph.clear();
                }
            }

            cout << endl;
            cout << "==============================================" << endl;
            cout << "所有超图指标计算完成。" << endl;
            cout << "==============================================" << endl;
        }
        else
        {
            // ====================================================
            // 创建结果根目录
            // ====================================================
            filesystem::path result_root =
                filesystem::path(params.cin_path) / "SimulationResults";
            filesystem::create_directories(result_root);
            params.cout_path = result_root.string();
            // ====================================================
            // 创建10个对应的结果文件夹
            // ====================================================
            for (int i = 1; i <= 10; i++)
            {
                filesystem::path result_folder =
                    result_root / ("Hypergraph_" + to_string(i));
                filesystem::create_directories(result_folder);
                cout << "创建结果目录: "
                    << result_folder << endl;
            }
            cout << endl;
            cout << "超图指标结果目录创建完成。" << endl;
            // ====================================================
            // 选择仿真结果内容
            // ====================================================
            cout << endl;
            cout << "选择横坐标" << endl;
            cout << "1.超度(不同时刻下,合作频率，收益和度之间的关系，记录1~100，1000，10000)" << endl;
            cout << "2.时间（不同时刻下,合作频率，合作者收益，合作者超度和时间之间的关系，记录1~100000）" << endl;
            cout << "3.协同因子" << endl;
            int choice;
            cin >> choice;
            if (choice == 1)
            {
                
                cout << "输入协同因子(归一化)\n";
                double jay;
                cin >> jay;
                cout << "输入重复次数\n";
                cin >> retime;
                vector<int>T;
                for (int i = 1; i <= 101; i = i + 1)T.push_back(i);
                T.push_back(1001);
                T.push_back(10001);
                cout << "输入Hypergraph范围【l-r】\n";
                int left, right;
                cin >> left >> right;
                for (int i = left; i <= right; i++)
                {
                    // ------------------------------------------------
                    // 当前输入文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cin_path =
                        filesystem::path(params.cin_path)
                        / ("Hypergraph_" + to_string(i));
                    // ------------------------------------------------
                    // 当前输出文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cout_path =
                        filesystem::path(params.cout_path)
                        / ("Hypergraph_" + to_string(i));
                    cout << endl;
                    cout << "==============================================" << endl;
                    cout << "处理 Hypergraph_" << i << endl;
                    cout << "输入目录: " << folder_cin_path << endl;
                    cout << "输出目录: " << folder_cout_path << endl;
                    cout << "==============================================" << endl;
                    // ------------------------------------------------
                    // 检查输入目录
                    // ------------------------------------------------
                    if (!filesystem::exists(folder_cin_path))
                    {
                        cout << "输入目录不存在，跳过。" << endl;
                        continue;
                    }
                    // ------------------------------------------------
                    // 遍历当前文件夹中的所有 txt
                    // ------------------------------------------------
                    for (const auto& entry :
                        filesystem::directory_iterator(folder_cin_path))
                    {
                        if (!entry.is_regular_file())
                            continue;

                        if (entry.path().extension() != ".txt")
                            continue;

                        cout << endl;
                        cout << "Processing file: "
                            << entry.path() << endl;

                        // ====================================================
                        // 1. 读取超图
                        // ====================================================
                        HypergraphLib::Hypergraph hypergraph;

                        hypergraph =
                            HypergraphLib::create_graph_from_file(entry);

                        // ====================================================
                        // 2. 仿真
                        // ====================================================
                        for (auto tt : T)
                        {
                            step_num = tt;
                            init_step_num = tt;
                            r.clear();
                            Fc.clear();
                            {
                                r.push_back(jay);
                                Fc.push_back(0);
                            }
                            map<double, double>fck;
                            map<double, double>pik;
                            map<double, double>numk;
                            for (int i = 0; i < hypergraph.N; i++)
                            {
                                double k = hypergraph.Groups[i].size();
                                fck[k] = 0;
                                pik[k] = 0;
                                numk[k] = 0;
                            }
                            // ================= 重复实验 =================
                            for (int ii = 0; ii < r.size(); ii++)
                            {
                                //t=2 10 100 1000
                                // 遍历所有 r 值，逐个进行仿真
                                for (int i = 0; i < retime; i++)
                                {
                                    // 运行演化过程
                                    PGG_avpi::run(hypergraph, r[ii], step_num, init_step_num, y, timescale);
                                    for (int j = 0; j < hypergraph.N; j++)
                                    {
                                        double k = hypergraph.Groups[j].size();
                                        fck[k] += PGG_avpi::Vectex_status_list[j];
                                        pik[k] += PGG_avpi::Vectex_payoff_list[j];
                                        numk[k] += 1;
                                    }
                                    // 清理本轮仿真状态
                                    PGG_avpi::clear();
                                }
                                for (auto& [k, v] : fck)
                                {
                                    fck[k] /= (1.0 * numk[k]);
                                    pik[k] /= (1.0 * numk[k]);
                                }
                            }
                            // ====================================================
                            // 3. 创建输出文件名
                            // ====================================================
                            string input_filename =
                                entry.path().stem().string();

                            filesystem::path output_file =
                                folder_cout_path /
                                (input_filename + "_Simulation1_K_R" + to_string(4 * jay) + "_T" + to_string(tt - 1) + ".txt");

                            // ====================================================
                            // 4. 输出指标
                            // ====================================================
                            ofstream fout(output_file);
                            for (auto& [k, v] : fck)
                            {
                                fout << k << ' ' << v << ' ' << pik[k] << endl;
                            }
                            fout.close();
                        }
                        // ====================================================
                        // 5. 清空当前超图
                        // ====================================================
                        hypergraph.clear();
                    }
                }
            }
            else if (choice == 2)
            {
                
                cout << "输入协同因子(归一化)\n";
                double jay;
                cin >> jay;
                cout << "输入重复次数\n";
                cin >> retime;
                cout << "输入截至时间\n";
                cin >> step_num;
                init_step_num = step_num;
                cout << "输入Hypergraph范围【l-r】\n";
                int left, right;
                cin >> left >> right;
                for (int i = left; i <= right; i++)
                {
                    // ------------------------------------------------
                    // 当前输入文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cin_path =
                        filesystem::path(params.cin_path)
                        / ("Hypergraph_" + to_string(i));
                    // ------------------------------------------------
                    // 当前输出文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cout_path =
                        filesystem::path(params.cout_path)
                        / ("Hypergraph_" + to_string(i));
                    cout << endl;
                    cout << "==============================================" << endl;
                    cout << "处理 Hypergraph_" << i << endl;
                    cout << "输入目录: " << folder_cin_path << endl;
                    cout << "输出目录: " << folder_cout_path << endl;
                    cout << "==============================================" << endl;
                    // ------------------------------------------------
                    // 检查输入目录
                    // ------------------------------------------------
                    if (!filesystem::exists(folder_cin_path))
                    {
                        cout << "输入目录不存在，跳过。" << endl;
                        continue;
                    }
                    // ------------------------------------------------
                    // 遍历当前文件夹中的所有 txt
                    // ------------------------------------------------
                    for (const auto& entry :
                        filesystem::directory_iterator(folder_cin_path))
                    {
                        if (!entry.is_regular_file())
                            continue;

                        if (entry.path().extension() != ".txt")
                            continue;

                        cout << endl;
                        cout << "Processing file: "
                            << entry.path() << endl;

                        // ====================================================
                        // 1. 读取超图
                        // ====================================================
                        HypergraphLib::Hypergraph hypergraph;

                        hypergraph =
                            HypergraphLib::create_graph_from_file(entry);

                        // ====================================================
                        // 2. 仿真
                        // ====================================================
                        {
                            r.clear();r.push_back(jay);
							vector<double>fct;
           
                            for (int i = 0; i <= step_num; i++)
                            {
								fct.push_back(0);
                            }
                            // ================= 重复实验 =================
                            for (int ii = 0; ii < r.size(); ii++)
                            {
                                // 遍历所有 r 值，逐个进行仿真
                                for (int i = 0; i < retime; i++)
                                {
                                    // 运行演化过程
                                    PGG_avpi::run(hypergraph, r[ii], step_num, init_step_num, y, timescale);
                                    for (int i = 1; i <= step_num; i++)
                                    {
                                        fct[i] += ((PGG_avpi::Status_list[i] * 1.0) / (hypergraph.N * 1.0));
                                    }
                                    
                                    // 清理本轮仿真状态
                                    PGG_avpi::clear();
                                }
                                
                            }
                            // ====================================================
                            // 3. 创建输出文件名
                            // ====================================================
                            string input_filename =
                                entry.path().stem().string();

                            filesystem::path output_file =
                                folder_cout_path /
                                (input_filename + "_Simulation2_T_R" + to_string(4 * jay)+ ".txt");

                            // ====================================================
                            // 4. 输出指标
                            // ====================================================
                            int t = 0;
                            ofstream fout(output_file);
                            for (auto& v : fct)
                            {
                    	        v /= (1.0 * retime);
                    	        fout <<t<<' '<< v << endl;
                                t++;
                            }
                            fout.close();
                            vector<double>().swap(fct);
                        }
                        // ====================================================
                        // 5. 清空当前超图
                        // ====================================================
                        hypergraph.clear();
                    }
                }
            }
            else if (choice == 3)
            {
                cout << "输入协同因子(归一化)\nstart end step\n";
                double start,end,st;
                cin >> start>>end>>st;
                cout << "输入重复次数\n";
                cin >> retime;
                cout << "输入截至时间init_step_num(step_num = init_step_num+5000)\n";
                cin >> init_step_num;
                step_num = init_step_num+5000;
                cout << "输入Hypergraph范围【l-r】\n";
				int left, right;
                cin >> left >> right;
                for (int i = left; i <= right; i++)
                {
                    // ------------------------------------------------
                    // 当前输入文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cin_path =
                        filesystem::path(params.cin_path)
                        / ("Hypergraph_" + to_string(i));
                    // ------------------------------------------------
                    // 当前输出文件夹
                    // ------------------------------------------------
                    filesystem::path folder_cout_path =
                        filesystem::path(params.cout_path)
                        / ("Hypergraph_" + to_string(i));
                    cout << endl;
                    cout << "==============================================" << endl;
                    cout << "处理 Hypergraph_" << i << endl;
                    cout << "输入目录: " << folder_cin_path << endl;
                    cout << "输出目录: " << folder_cout_path << endl;
                    cout << "==============================================" << endl;
                    // ------------------------------------------------
                    // 检查输入目录
                    // ------------------------------------------------
                    if (!filesystem::exists(folder_cin_path))
                    {
                        cout << "输入目录不存在，跳过。" << endl;
                        continue;
                    }
                    // ------------------------------------------------
                    // 遍历当前文件夹中的所有 txt
                    // ------------------------------------------------
                    for (const auto& entry :
                        filesystem::directory_iterator(folder_cin_path))
                    {
                        if (!entry.is_regular_file())
                            continue;

                        if (entry.path().extension() != ".txt")
                            continue;

                        cout << endl;
                        cout << "Processing file: "
                            << entry.path() << endl;

                        // ====================================================
                        // 1. 读取超图
                        // ====================================================
                        HypergraphLib::Hypergraph hypergraph;

                        hypergraph =
                            HypergraphLib::create_graph_from_file(entry);

                        // ====================================================
                        // 2. 仿真
                        // ====================================================
                        {
                            r.clear(); 
                            for(double i=start;i<=end;i+=st)
                            {
                                r.push_back(i);
								Fc.push_back(0);    
                            }
                            // ================= 重复实验 =================
                            for (int ii = 0; ii < r.size(); ii++)
                            {
                                // 遍历所有 r 值，逐个进行仿真
                                for (int i = 0; i < retime; i++)
                                {
                                    // 运行演化过程
                                    PGG_avpi::run(hypergraph, r[ii], step_num, init_step_num, y, timescale);
									Fc[ii] += PGG_avpi::fc;
                                    // 清理本轮仿真状态
                                    PGG_avpi::clear();
                                }

                            }
                            // ====================================================
                            // 3. 创建输出文件名
                            // ====================================================
                            string input_filename =
                                entry.path().stem().string();

                            filesystem::path output_file =
                                folder_cout_path /
                                (input_filename + "_Simulation3.txt");

                            // ====================================================
                            // 4. 输出指标
                            // ====================================================
                            ofstream fout(output_file);
                            for (int i=0;i<r.size();i++)
                            {
                                fout << r[i] << ' ' << Fc[i]/(retime*1.0) << endl;
                            }
                            fout.close();
                        }
                        // ====================================================
                        // 5. 清空当前超图
                        // ====================================================
                        hypergraph.clear();
                    }
                }
            }
        }
    }
    else
    {
        cout << endl;
        cout << "输入错误！请输入 1 或 2。" << endl;

        return 1;
    }
}