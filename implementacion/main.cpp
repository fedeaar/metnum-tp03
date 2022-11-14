#include "src/run.h"


//
// MAIN
//

int main(int argc,  char** argv) {

    if (argc < 3 || argc > 17) {
        cout << "error: cantidad invalida de parametros.\n" <<
             "expected: [source] [p_value]\n"  <<
             "          -m      (metodo)            ['EG' | 'GS' | 'J'],  default = 'EG'\n" <<
             "          -tol    (tolerancia GS, J)  [double],             default = 1e-4\n" <<
             "          -iter   (max iter. GS, J)   [entero],             default = 1e5\n"  <<
             "          -o      (out dir)           [string],             default = './'\n" <<
             "          -as     (save as)           [string],             default = (nombre del source)\n" <<
             "          -p      (precision)         [uint(0, 15)],        default = 15\n" <<
             "          -time   (medir tiempos)                           default = false\n" <<
             "          -v      (verbose)                                 default = false" << endl;
        return -1;
    }

    // RUN
    params program = get(argc, argv);
    run(program);

    return 0;
}
