/**/
    #ifndef FSTCATALOG_H
    #define FSTCATALOG_H
    #include <vector>
    #include <fst/fstlib.h>
    #include "MMatrix.h"

    using namespace fst;
    typedef list< pair<StdVectorFst, int> > FSTlist;
    typedef vector< MMatrix > MatrixGroup;

    class FSTcatalog
    {
        public:
            FSTlist popID;
            vector<double> popdist; // container for population type distribution...a
            MatrixGroup intxnNet;

            FSTcatalog(vector<StdVectorFst> V);
            virtual ~FSTcatalog();

            void update(StdVectorFst result, int d, int T1type, int T2type);
        protected:
        private:
    };

    #endif // FSTCATALOG_H
/**/
