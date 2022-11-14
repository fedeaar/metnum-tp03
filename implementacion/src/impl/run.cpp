#include "../run.h"


//
// PARAMS
//

params get(int argc,  char** argv) {

    string error;
    params res;

    // path
    res.string_params["in_path"] = argv[1];
    string file_name = IO::filename(argv[1]);
    res.string_params["file_name"] = file_name;

    // p value
    error = "error: el valor p no se pudo interpretar como doble.";
    res.double_params["p_val"] = IO::stodcast(argv[2], error);


    // params opcionales
    map <string, string> op = IO::oparams(argc, argv);

    // formato
    res.string_params["formato"] = op.count("-f") ? op.at("-f") : "matriz";

    // metodo
    res.string_params["metodo"] = op.count("-m") ? op.at("-m") : "EG";

    // tol
    error = "error: el valor de tolerancia no se pudo interpretar como doble.";
    res.double_params["tol"] = op.count("-tol") ? IO::stodcast(op.at("-tol"), error) : 1e-4;

    // iter
    error = "error: el valor de iteraciones no se pudo interpretar como entero.";
    res.size_t_params["iter"] = op.count("-iter") ? IO::stodcast(op.at("-iter"), error) : 1e5;

    // out dir
    string out_path = op.count("-o") ? op.at("-o") : "./";
    if (out_path[out_path.length() - 1] != '/') {
        out_path += '/';
    }
    res.string_params["out_path"] = out_path;

    // out name
    res.string_params["out_name"] = op.count("-as") ? op.at("-as") : file_name;

    // precision
    error = "error: el valor de precision no se pudo interpretar como entero.";
    res.size_t_params["precision"] = op.count("-p") ? IO::stodcast(op.at("-p"),  error) : 15;

    // verbose
    res.bool_params["verbose"] = op.count("-v");

    // time
    res.bool_params["time"] = op.count("-time");

    return res;
}




//
// RUN
//

void run(const params &program) {

    bool verbose = program.bool_params.at("verbose");

    // read p_val
    double p_val = program.double_params.at("p_val");

    // read tol
    double tol = program.double_params.at("tol");

    // read iter
    size_t iter = program.size_t_params.at("iter");

    // make matriz pagerank
    string in_path =  program.string_params.at("in_path");
    if (verbose) std::cout << "leyendo el archivo (como grafo): " + in_path << "\n";
    pagerank::IO::in_file data = pagerank::IO::read_in(in_path, p_val);

    // metodo
    string metodo = program.string_params.at("metodo");
    pagerank::metodo pr_metodo = metodo == "GS" ? pagerank::GS : metodo == "J" ? pagerank::J : pagerank::EG;

    // calcular
    if (verbose) std::cout << "resolviendo pagerank (" << metodo << ")...\n";
    auto inicio_make = chrono::high_resolution_clock::now();
    Eigen::SparseMatrix<double> A = pagerank::make(data);
    auto fin_make = chrono::high_resolution_clock::now();
    auto time_make = chrono::duration_cast<chrono::microseconds>(fin_make - inicio_make);
    if (verbose) std::cout << "tiempo de ejecucion etapa make (microsegundos): " << time_make.count() << '\n';

    auto inicio_solve = chrono::high_resolution_clock::now();
    Eigen::VectorXd res = pagerank::solve(A, pr_metodo, tol, iter);
    auto fin_solve = chrono::high_resolution_clock::now();
    auto time_solve = chrono::duration_cast<chrono::microseconds>(fin_solve - inicio_solve);
    if (verbose) std::cout << "tiempo de ejecucion etapa solve (microsegundos): " << time_solve.count() << '\n';

    // write
    string out = program.string_params.at("out_path") + program.string_params.at("out_name");
    size_t precision = program.size_t_params.at("precision");
    if (verbose) std::cout << "guardando resultados en: " << out << ".out (si el path existe)\n";
    pagerank::IO::write_out(out + ".out", {res, p_val}, (int) precision);

    // write time
    if (program.bool_params.at("time")) {
        string time_out = out + ".time";
        if (verbose) std::cout << "guardando tiempo de ejecucion en: " + time_out << " (si el path existe)" << endl;
        IO::write_time(time_out, {{"make", time_make}, {"solve", time_solve}});
    }
}
