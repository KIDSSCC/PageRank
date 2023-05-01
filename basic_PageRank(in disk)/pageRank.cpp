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

void base_pagerank(const string& edges_file, const string& sparse_matrix_file,
    const string& r_old_file, double beta, double epsilon) {
    vector<int> nodes = get_sparse_matrix(edges_file, sparse_matrix_file);
    map<int, int> id2idx;
    int idx = 0;
    for (const auto& id : nodes) {
        id2idx[id] = idx++;
    }

    int nodes_num = nodes.size();
    vector<double> r(nodes_num, 1.0 / nodes_num);//pageRank值

    ofstream r_old(r_old_file);
    for (const auto& r_ : r) {
        r_old << r_ << "\n";
    }
    r_old.close();

    int max_iter = 100;
    cout << "iteration start!" << endl;
    for (int iter = 0; iter < max_iter; ++iter) {
        ifstream r_file(r_old_file);
        double init = (1 - beta) / nodes_num;
        vector<double> r_new(nodes_num, init);

        ifstream sparse_matrix(sparse_matrix_file);
        string line;
        string r_line;
        int r_index = 0;
        while (getline(sparse_matrix, line)){
            stringstream ss(line);
            int from_, degree;
            ss >> from_ >> degree;

            int from_idx = id2idx[from_];

            double r_;
            while (r_index != from_idx)
            {
                r_index++;
                getline(r_file, r_line);
            }
            getline(r_file, r_line);
            stringstream r_ss(r_line);
            r_ss >> r_;
            r_index++;

        /*    if (iter!=0&&(r_ != r[from_idx]))
            {
                cout << r_ - r[from_idx] << endl;
                cout << r_ << " " << r[from_idx] << endl;
                cout << "精度问题" << endl;

            }*/

            for (int i = 0; i < degree; ++i) {
                int to;
                ss >> to;
                int idx = id2idx[to];
                r_new[idx] += beta * r_ / degree;
            }
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
        ifstream r_old_again(r_old_file);
        for (size_t i = 0; i < r.size(); ++i) {
            double r_;
            r_old_again >> r_;
            r_ = r[i];
            err += abs(r_ - r_new[i]);
        }
        r_old_again.close();

        ofstream r_old_write(r_old_file);
        for (const auto& r_ : r_new) {
            r_old_write << r_ << "\n";
            
        }
        r_old_write.close();

        r = r_new;
        cout << "finish iter " << iter << endl;
        cout << "err " << err << endl;
        if (err < epsilon) {
            cout << "finish at iter " << iter << endl;
            break;
        }
    }

    std::vector<std::pair<int, double>> vec;
    for (size_t i = 0; i < r.size(); ++i) {
        vec.push_back({ nodes[i],r[i] });
    }

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
    string r_old_file = "r.txt";
    double beta = 0.85;
    double epsilon = 1e-6;

    clock_t start = clock();
    base_pagerank(edges_file, sparse_matrix_file, r_old_file, beta, epsilon);
    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "total time: " << elapsed_time << "s." << endl;

    return 0;
}
