/**/
    #ifndef MMATRIX_H
    #define MMATRIX_H

   // #include <vector>
    #include <bitset>
    #include <tnt/tnt_array2d.h>

    using namespace TNT;
//    using namespace std;

    class MMatrix : public Array2D<bool>
    {
        public:
            MMatrix(int x, int y);
            virtual ~MMatrix();
            void removerow(int i);
            void removecolumn(int i);
            void addrow();
            void addcolumn();
           // void rownorm();
        protected:
        private:
    };

    #endif // MMATRIX_H
/**/
