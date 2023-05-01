#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <string>
#include <iomanip>

using namespace std;

struct NodeData {
    int degree=0;
    vector<int> tos;
};

// 定义比较函数，用于std::sort的排序规则
bool compare(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second > b.second;
}
void pagerank_inmemory(string edges_file, string outfile, double beta, double epsilon)
{
    map<int, NodeData> m;
    set<int> nodes;
    ifstream infile(edges_file);
    string line;
    while (getline(infile, line))
    {
        stringstream ss(line);
        int from_, to;
        ss >> from_ >> to;
        m[from_].degree++;
        m[from_].tos.push_back(to);
        nodes.insert(from_);
        nodes.insert(to);
    }
    //有一个问题需要考虑，是不是有重复的边
    for (auto item : m)
    {
        NodeData data = item.second;

        sort(data.tos.begin(), data.tos.end());
    }
    vector<int> sorted_nodes(nodes.begin(), nodes.end());
    //sorted_nodes存储了所有出度不为0的点。
    map<int, int> id2Serial;
    int idx = 0;
    for (auto id : sorted_nodes)
    {
        id2Serial[id] = idx++;
    }
    int nodes_num = sorted_nodes.size();
    vector<double> r_old(nodes_num, 1.0 / nodes_num);//pageRank值
    int max_iter = 100;
    cout << "iteration start!" << endl;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        double init = 0.0;
        vector<double> r_new(nodes_num, init);
        int r_index = 0;
        map<int, NodeData>::iterator iter_Matrix = m.begin();
        while (iter_Matrix != m.end())
        {
            int from_ = (*iter_Matrix).first;
            int degree = (*iter_Matrix).second.degree;
            int from_idx = id2Serial[from_];
            double r_;
            while (r_index != from_idx)
            {
                r_index++;
            }
            r_ = r_old[r_index] / degree;
            r_index++;
            for (int i = 0; i < degree; ++i)
            {
                int to = (*iter_Matrix).second.tos[i];
                int idx = id2Serial[to];
                r_new[idx] += beta * r_;
            }
            iter_Matrix++;
        }
        
        double r_new_sum = 0;
        for (const auto& r_ : r_new) {
            r_new_sum += r_;
        }
        double delta = (1 - r_new_sum) / nodes_num;
        for (auto& r_ : r_new) {
            r_ += delta;
        }
        double err = 0;
        for (int i = 0; i < r_old.size(); i++)
        {
            err += abs(r_old[i] - r_new[i]);
        }
        for (int i = 0; i < r_old.size(); i++)
        {
            r_old[i] = r_new[i];
        }
        cout << "finish iter " << iter << endl;
        cout << "err " << err << endl;
        if (err < epsilon) {
            cout << "finish at iter " << iter << endl;
            break;
        }
    }
    vector<pair<int, double>>vec;
    for (int i = 0; i < r_old.size(); i++)
    {
        vec.push_back({ sorted_nodes[i],r_old[i] });
    }
    sort(vec.begin(), vec.end(), compare);
    ofstream result(outfile);
    int num = 0;
    for (vector<std::pair<int, double>>::iterator it = vec.begin(); it != vec.end(); ++it) {
        result << it->first << "\t" << setprecision(10) << it->second << "\n";
        num++;
        if (num == 100) {
            break;
        }
    }

}
int main() {
    string edges_file = "Data.txt";
    string outfile = "inmemory.txt";
    double beta = 0.85;
    double epsilon = 1e-6;

    clock_t start = clock();
    pagerank_inmemory(edges_file, outfile, beta, epsilon);
    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "total time: " << elapsed_time << "s." << endl;

    return 0;
}
