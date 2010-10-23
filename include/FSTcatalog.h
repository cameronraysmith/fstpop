/**/
    #ifndef FSTCATALOG_H
    #define FSTCATALOG_H
    #include <vector>
    #include <cmath>

    //------------OpenFST-------------//
    #include <fst/fstlib.h>
    #include <fst/extensions/far/farlib.h>

    //--------------GSL----------------//
    #include <gsl/gsl_math.h>
    #include <gsl/gsl_eigen.h>
    #include <gsl/gsl_vector.h>
    #include <gsl/gsl_matrix.h>
    #include <gsl/gsl_blas.h>

    #include "MMatrix.h"


    using namespace fst;
    typedef list< pair<StdVectorFst, int> > FSTlist;
    typedef vector< MMatrix > MatrixGroup;
    extern ofstream fstpoplog;

    class FSTcatalog
    {
        public:
            FSTlist popID;
            vector<double> popdist; // container for population type distribution...a
            MatrixGroup intxnNet;
            int N; // population size

            FSTcatalog(vector<StdVectorFst> V); // constructor
            virtual ~FSTcatalog();

            void update(StdVectorFst result, int d, int T1type, int T2type);
            double ncomplexity();
            double scomplexity();
            void printpop(string outpdf);
        protected:
        private:
    };

    #endif // FSTCATALOG_H
/**/
