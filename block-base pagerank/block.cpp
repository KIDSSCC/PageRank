#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <string>
#include <iomanip>

using namespace std;

struct NodeData {
    int degree;
    vector<int> tos;
};

// 定义比较函数，用于std::sort的排序规则
bool compare(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second > b.second;
}

vector<int> get_sparse_matrix(const string& input_file, const string& sparse_matrix_file) {
    map<int, NodeData> m;
    set<int> nodes;

    ifstream infile(input_file);
    string line;

    while (getline(infile, line)) {
        stringstream ss(line);
        int from_, to;
        ss >> from_ >> to;

        m[from_].degree++;
        m[from_].tos.push_back(to);
    }

    ofstream outfile(sparse_matrix_file);

    for (const auto& item : m) {
        int from_ = item.first;
        NodeData data = item.second;

        sort(data.tos.begin(), data.tos.end());
        outfile << from_ << " " << data.degree << " ";

        for (size_t i = 0; i < data.tos.size(); ++i) {
            outfile << data.tos[i];

            if (i < data.tos.size() - 1) {
                outfile << " ";
            }
        }

        outfile << "\n";
        nodes.insert(from_);

        for (const auto& to : data.tos) {
            nodes.insert(to);
        }
    }
    vector<int> sorted_nodes(nodes.begin(), nodes.end());
    return sorted_nodes;
}

void block_based_pagerank(const string& edges_file, const string& sparse_matrix_file, string& r1_file, string& r2_file, double beta, double epsilon, int bsize) {
    vector<int> nodes = get_sparse_matrix(edges_file, sparse_matrix_file);
    map<int, int> id2idx;
    int idx = 0;
    for (const auto& id : nodes) {
        id2idx[id] = idx++;
    }

    int nodes_num = nodes.size();
    vector<double> r(nodes_num, 1.0 / nodes_num);//pageRank值

    ofstream r1_file_out(r1_file);
    for (const auto& r_ : r) {
        r1_file_out << r_ << "\n";
    }
    r1_file_out.close();

    vector<pair<int, int>> block_ranges;
    for (int start = 0; start < nodes_num; start += bsize) {
        block_ranges.push_back(make_pair(start, min(start + bsize, nodes_num)));
    }

    double delta_new = 0;
    int max_iter = 100;
    cout << "iteration start!" << endl;
    for (int iter = 0; iter < max_iter; ++iter) {
        double delta_old = delta_new;
        delta_new = 0;
        double err = 0;
        double r_new_sum = 0;

        ofstream r2_file_out(r2_file);

        for (const auto& block_range : block_ranges) {
            unordered_map<int, double> r_new;
            for (int i = block_range.first; i < block_range.second; i++) {
                r_new[i] = 0;
            }

            ifstream sparse_matrix(sparse_matrix_file);
            ifstream r1_file_in(r1_file);

            string line;
            string r_line;
            int r_old_file_idx = 0;
            while (getline(sparse_matrix, line)) {
                stringstream ss(line);
                int from, degree;
                ss >> from >> degree;

                int from_idx = id2idx[from];

                double r_;
                while (r_old_file_idx != from_idx)
                {
                    r_old_file_idx++;
                    getline(r1_file_in, r_line);
                }
                getline(r1_file_in, r_line);
                stringstream r_ss(r_line);
                r_ss >> r_;
                r_ += delta_old;
                r_old_file_idx++;

                for (int i = 0; i < degree; ++i) {
                    int to;
                    ss >> to;
                    int idx = id2idx[to];
                    if (block_range.first <= idx && idx < block_range.second) {
                        r_new[idx] += beta * r_ / degree;
                        r_new_sum += beta * r_ / degree;
                    }
                }
            }
            sparse_matrix.close();
            r1_file_in.close();
            // 计算err
            r_old_file_idx = 0;
            ifstream r_old_again(r1_file);
            for (int i = block_range.first; i < block_range.second; i++) {
                double r_;
                while (r_old_file_idx != i)
                {
                    r_old_file_idx++;
                    getline(r_old_again, r_line);
                }
                r_old_file_idx++;
                getline(r_old_again, r_line);
                stringstream r_ss(r_line);
                r_ss >> r_;
                err += abs(r_ - r_new[i]);
            }
            r_old_again.close();
            // 写入r2
            for (int i = block_range.first; i < block_range.second; i++) {
                r2_file_out << r_new[i] << "\n";
            }
        }
        r2_file_out.close();
        delta_new = (1 - r_new_sum) / nodes_num;
        swap(r1_file, r2_file);

        cout << "finish iter " << iter << endl;
        cout << "err " << err << endl;
        if (err < epsilon) {
            cout << "finish at iter " << iter << endl;
            break;
        }
    }

    std::vector<std::pair<int, double>> vec;
    // 计分
    ifstream r_old_again_again(r1_file);
    int r_old_file_idx = 0;
    string r_line;
    for (size_t i = 0; i < r.size(); ++i) {
        double r_;
        while (r_old_file_idx != i)
        {
            r_old_file_idx++;
            getline(r_old_again_again, r_line);
        }
        r_old_file_idx++;
        getline(r_old_again_again, r_line);
        stringstream r_ss(r_line);
        r_ss >> r_;
        vec.push_back({ nodes[i],r_ });
    }
    r_old_again_again.close();

    // 使用自定义比较函数对vector中的元素进行降序排序
    std::sort(vec.begin(), vec.end(), compare);

    ofstream result("result_base.txt");
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
    string sparse_matrix_file = "data_sparse.txt";
    string r_old_file = "r1.txt";
    string r_new_file = "r2.txt";
    int bsize = 1000;
    double beta = 0.85;
    double epsilon = 1e-6;

    clock_t start = clock();
    block_based_pagerank(edges_file, sparse_matrix_file, r_old_file, r_new_file, beta, epsilon, bsize);
    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "total time: " << elapsed_time << "s." << endl;

    return 0;
}