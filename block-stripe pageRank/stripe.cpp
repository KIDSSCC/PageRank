#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<iomanip>
using namespace std;
vector<pair<int, int>> block_ranges;
map<int, int>id2Index;
struct NodeData {
    int degree = 0;
    vector<int> tos;
};
// 定义比较函数，用于sort的排序规则
bool compare(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second > b.second;
}
vector<int> get_stripe_sparse_matrix(const string& edges_file, const string& matrix_file, int block_size)
{
    map<int, NodeData> m;
    set<int> nodes;
    ifstream infile(edges_file);
    string line;
    while (getline(infile, line)) {
        stringstream ss(line);
        int from, to;
        ss >> from >> to;
        m[from].degree++;
        m[from].tos.push_back(to);
        nodes.insert(from);
        nodes.insert(to);
    }
    vector<int> Index2id(nodes.begin(), nodes.end());
    int node_index = 0;
    for (int i = 0; i < nodes.size(); i++)
        id2Index[Index2id[i]] = node_index++;

    int nodes_num = Index2id.size();
    int blocks_num = nodes_num / block_size + 1;
    //确定各个块之间的节点范围，保存在block_ranges中
    for (int start = 0; start < nodes_num; start += block_size)
    {
        block_ranges.push_back(make_pair(start, min(start + block_size, nodes_num)));
    }
    for (int i = 0; i < block_ranges.size(); i++)
    {
        //对于每一个块都生成其独有的邻接矩阵文件
        int begin = block_ranges[i].first;
        int end = block_ranges[i].second;
        //需要遍历map，找出所有对当前块中的节点有指向关系的节点。
        map<int, NodeData> Matrix_for_block;
        for (auto line : m)
        {
            int from = line.first;
            int degree = line.second.degree;
            bool scan = false;
            for (auto dest : line.second.tos)
            {
                if (id2Index[dest] >= begin && id2Index[dest] < end)
                {
                    scan = true;
                    Matrix_for_block[from].tos.push_back(dest);
                }
            }
            if (scan)
                Matrix_for_block[from].degree = degree;
        }
        //输出 Matrix_for_block
        stringstream ss;
        ss << i;
        string outfile = matrix_file + ss.str() + ".txt";
        ofstream out(outfile);
        for (auto each : Matrix_for_block)
        {
            out << each.first << " ";
            out << each.second.degree << " ";
            for (int j = 0; j < each.second.tos.size(); j++)
            {
                out << each.second.tos[j] << " ";
            }
            out << "\n";
        }
    }
    return Index2id;

}
void stripe_pageRank(const string& edges_file, const string& sparse_matrix_file, string& r1, string& r2,
    double beta, double epsilon, int block_size)
{
    //得到条纹稀疏矩阵
    vector<int>nodes = get_stripe_sparse_matrix(edges_file, sparse_matrix_file, block_size);
    int nodes_num = nodes.size();
    //写入初始数据
    double init = 1.0 / nodes_num;
    ofstream r_old(r1);
    for (int i = 0; i < nodes_num; i++)
    {
        r_old <<setprecision(10)<< init << "\n";
    }
    r_old.close();

    double delta = 0.0;
    double err = 0.0;
    for (int i = 0; i < 100; i++)
    {
        double r_new_sum = 0.0;
        err = 0.0;

        string middle_file = "Middle.txt";
        ofstream r_middle_file(middle_file, ios::trunc);
        for (int j = 0; j < block_ranges.size(); j++)
        {
            
            vector<double> curr_block(block_ranges[j].second - block_ranges[j].first, 0);
            //在每次迭代中，遍历每一个矩阵块
            stringstream ss;
            ss << j;
            string matrix_file = sparse_matrix_file + ss.str() + ".txt";
            //读取邻接矩阵
            ifstream martix_In(matrix_file);
            ifstream old_In(r1);
            
            string eachline;
            string oldline;

            int r_old_index = 0;
            while (getline(martix_In, eachline))
            {
                stringstream ss(eachline);
                int from, degree;
                ss >> from >> degree;
                int from_index = id2Index[from];

                while (r_old_index != from_index)
                {
                    getline(old_In, oldline);
                    r_old_index++;
                }
                getline(old_In, oldline);
                r_old_index++;

                double old_score;
                stringstream r_ss(oldline);
                r_ss.precision(std::numeric_limits<double>::max_digits10);
                r_ss >> old_score;
                double contribute = beta*old_score / degree;
                //from节点：from,degree,old_score,以及其指向的若干节点
                int dest;
                while (ss >> dest)
                {
                    int dest_index = id2Index[dest];
                    curr_block[dest_index - block_ranges[j].first] += contribute;
                }
            }
            martix_In.close();
            old_In.close();
            //准备进行写入
            for (int k = 0; k < curr_block.size(); k++)
            {
                r_new_sum += curr_block[k];
                r_middle_file<<setprecision(10) << curr_block[k];
                if ((j != block_ranges.size() - 1) || (k != curr_block.size() - 1))
                {
                    r_middle_file << "\n";
                }
            }
        }
        
        r_middle_file.close();
        //到这里计算完了一次迭代中所有的块
        delta = (1 - r_new_sum) / nodes_num;
        ifstream middle_result(middle_file);
        ifstream old_result(r1);
        string middleline;
        string oldline;

        ofstream new_result(r2);
        while (getline(middle_result, middleline))
        {
            getline(old_result, oldline);
            stringstream middle(middleline);
            stringstream old(oldline);
            double middle_score, oldscore;
            middle >> middle_score;
            old >> oldscore;
            err += abs(middle_score + delta - oldscore);
            new_result <<setprecision(10) << middle_score + delta << "\n";
        }
        swap(r1, r2);
        cout << "finish iter: " << i << endl;
        cout << "err is: " << err << endl;
        if (err < epsilon)
        {
            break;
        }
    }
    //此时r1中是最后一次的结果
    vector<double> res;
    ifstream finalresult(r1);
    string r_old_line;
    while (getline(finalresult, r_old_line))
    {
        stringstream r_ss(r_old_line);
        double r_;
        r_ss >> r_;
        res.push_back(r_);
    }
    cout << res.size() << endl;
    vector<pair<int, double>> vec;
    for (size_t i = 0; i < res.size(); ++i) {
        vec.push_back({ nodes[i],res[i] });
    }
    // 使用自定义比较函数对vector中的元素进行降序排序
    sort(vec.begin(), vec.end(), compare);
    ofstream result("result_stripe.txt");
    int num = 0;

    for (vector<pair<int, double>>::iterator it = vec.begin(); it != vec.end(); ++it) {
        result << it->first << "\t" << setprecision(10) << it->second << "\n";
        num++;

        if (num == 100) {
            break;
        }
    }
}
int main() {
    string edges_file = "Data.txt";
    string sparse_matrix_file = "data_sparse";
    string r_old_file = "r.txt";
    string r1 = "r1.txt";
    string r2 = "r2.txt";
    double beta = 0.85;
    double epsilon = 1e-6;
    int block_size = 1000;
    clock_t start = clock();
    stripe_pageRank(edges_file, sparse_matrix_file, r1, r2, beta, epsilon, block_size);
    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "total time: " << elapsed_time << "s." << endl;

}