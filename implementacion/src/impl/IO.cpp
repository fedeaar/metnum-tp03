#include "../IO.h"


//
// UTILS
//

size_t IO::stolcast(const string &val, const string &msg) {
    size_t res;
    try {
        res = std::stoll(val);
    } catch(std::invalid_argument &ia) {
        throw std::invalid_argument(msg);
    }
    return res;
}


double IO::stodcast(const string &val, const string &msg) {
    double res;
    try {
        res = std::stod(val);
    } catch(std::invalid_argument &ia) {
        throw std::invalid_argument(msg);
    }
    return res;
}


string IO::filename(const string &path) {
    string name;
    // encontrar nombre
    auto delim = path.find_last_of('/');
    if (delim != string::npos) {
        name = path.substr(delim + 1, path.size());
    } else {
        name = path;
    }
    // remover extension
    delim = name.find_last_of('.');
    if (delim != string::npos) {
        name = name.substr(0, delim);
    }
    return name;
}


map<string, string> IO::oparams(int argc,  char** argv) {
    map<string, string> params;
    for (int i = 0; i < argc; ++i) {
        if (argv[i][0] != '-') {
            continue;
        }
        string param = argv[i];
        if (i + 1 < argc && argv[i+1][0] != '-') {
            string val = argv[i+1];
            params[param] = val;
            ++i;
        } else {
            params[param] = "true";
        }
    }
    return params;
}


void IO::skip_lines(ifstream &file, size_t n) {
    size_t i = 0;
    while (i++ < n) {
        file.ignore(unsigned(-1), '\n');
    }
}




//
// FILE HANDLING
//

pair<size_t, size_t> IO::_shape(ifstream &file) {
    auto pos = file.tellg();
    string _line;
    size_t n {}, m {};
    getline(file, _line);
    if (!_line.empty()) {
        ++n;
        for (char c : _line) {
            m += c == ' ';
        }
        ++m;
    }
    while (getline(file, _line)) {
        ++n;
    }
    file.clear();
    file.seekg(pos);
    return {n, m};
}


graph IO::read_grafo(const string &in, size_t start) {
    ifstream file{in, ios::binary};
    if (!file.is_open()) {
        throw std::invalid_argument("no se pudo leer el archivo: " + in + ".");
    }
    skip_lines(file, start);
    // cantidad de nodos
    string _nodos {};
    std::getline(file, _nodos);
    size_t nodos = stodcast(_nodos, "error de formato: linea 1.");
    // cantidad de links
    string _links {};
    std::getline(file, _links);
    size_t links = stodcast(_links, "error de formato: linea 2.");
    // init store
    graph res(nodos, links);
    // in coords
    string _i, _j;
    size_t i, j, k = 2;
    while (file >> _j >> _i) {
        string msg = "error de formato: linea " + std::to_string(k) + ".";
        i = stolcast(_i, msg);
        j = stolcast(_j, msg);
        res.relaciones.emplace_back(coords{i, j});
        ++k;
    }
    file.close();
    return res;
}


Eigen::SparseMatrix<double, Eigen::RowMajor> IO::read_matriz(const string &in, size_t start) {
    ifstream file{in, ios::binary};
    if (!file.is_open()) {
        throw std::invalid_argument("no se pudo leer el archivo: " + in + ".");
    }
    skip_lines(file, start);
    return read_matriz(file);
}
Eigen::SparseMatrix<double, Eigen::RowMajor> IO::read_matriz(ifstream &file) {
    // n, m
    pair<size_t, size_t> shape = _shape(file);
    long long n = shape.first, m = shape.second;
    // init
    std::vector<Eigen::Triplet<double>> t_list;
    size_t i = 0, j = 0, k = 0;
    string _e {};
    double e {};
    while (file >> _e) {
        string msg = "error de formato: linea " + std::to_string(k) + ".";
        ++k;
        e = stodcast(_e, msg);
        t_list.emplace_back(Eigen::Triplet<double>(i, j, e));
        ++j;
        if (j >= m) {
            j = 0;
            ++i;
        }
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> res{n, m};
    res.setFromTriplets(t_list.begin(), t_list.end());
    file.close();
    return res;
}


Eigen::VectorXd IO::read_vector(const string &in, size_t start) {
    ifstream file{in, ios::binary};
    if (!file.is_open()) {
        throw std::invalid_argument("no se pudo leer el archivo: " + in + ".");
    }
    skip_lines(file, start);
    return read_vector(file);
}
Eigen::VectorXd IO::read_vector(ifstream &file) {
    // n
    pair<size_t, size_t> shape = _shape(file);
    long long n = shape.first;
    // init
    Eigen::VectorXd res(n);
    long long i = 0;
    string _e {};
    double e {};
    while (file >> _e) {
        string msg = "error de formato: linea " + std::to_string(i) + ".";
        e = stodcast(_e, msg);
        res.coeffRef(i) = e;
        ++i;
    }
    file.close();
    return res;
}


void IO::write_matriz(const string &out, const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int precision) {
    ofstream file{out};
    write_matriz(file, mat, precision);
}
void IO::write_matriz(ofstream &file, const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int precision) {
    file << std::setprecision(precision) << std::fixed << mat << endl;
    file.close();
}


void IO::write_vector(const string &out, const Eigen::VectorXd &vec, int precision) {
    ofstream file{out};
    write_vector(file, vec, precision);
}
void IO::write_vector(ofstream &file, const Eigen::VectorXd &vec, int precision) {
    file << std::setprecision(precision) << std::fixed << vec << endl;
    file.close();
}


void IO::write_time(const string &out,  const vector<pair<string, chrono::microseconds>> &time_measures) {
    ofstream file {out};
    file << "medida,microseconds\n";
    for (auto &v: time_measures) {
        file << v.first << "," << v.second.count() << '\n';
    }
    file.close();
}
