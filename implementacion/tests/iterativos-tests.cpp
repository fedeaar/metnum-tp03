#include "../dependencies/gtest-1.8.1/gtest.h"
#include <Eigen/Sparse>
#include "../src/IO.h"
#include "../src/pagerank.h"

#if METODO == 0
pagerank::metodo M = pagerank::EG;
#elif METODO == 1
pagerank::metodo M = pagerank::GS;
#else
pagerank::metodo M = pagerank::J;
#endif


//
// INIT
//

class IterativosTests : public testing::Test {
protected:
    string basedir;
    double epsilon;
    size_t iter;
    void SetUp() override {
        basedir = "../../catedra/tests-pagerank/";
        epsilon = 1e-4;
        iter = 1e5;
    }

    void pagerank_test(const string &in, const string &out);
};


void IterativosTests::pagerank_test(const string &in, const string &out) {

    grafo g = IO::read_grafo(basedir + in);
    pagerank::IO::out_file expected = pagerank::IO::read_out(basedir + out);
    Eigen::SparseMatrix<double> A = pagerank::make({g, expected.p_val});
    Eigen::VectorXd res = pagerank::solve(A, M);
    EXPECT_LE((res - expected.solucion).norm(), epsilon);
}




//
// TESTS
//

TEST_F(IterativosTests, test_aleatorio) {
    pagerank_test("test_aleatorio.txt", "test_aleatorio_desordenado.txt.out");
}


TEST_F(IterativosTests, test_aleatorio_desordenado) {
    pagerank_test("test_aleatorio_desordenado.txt", "test_aleatorio_desordenado.txt.out");
}


TEST_F(IterativosTests, test_trivial) {
    pagerank_test("test_trivial.txt", "test_trivial.txt.out");
}


TEST_F(IterativosTests, test_completo) {
    pagerank_test("test_completo.txt", "test_completo.txt.out");
}


TEST_F(IterativosTests, test_sin_links) {
    pagerank_test("test_sin_links.txt", "test_sin_links.txt.out");
}


TEST_F(IterativosTests, test_15_segundos) {
    pagerank_test("test_15_segundos.txt", "test_15_segundos.txt.out");
}


TEST_F(IterativosTests, test_30_segundos) {
    pagerank_test("test_30_segundos.txt", "test_30_segundos.txt.out");
}
