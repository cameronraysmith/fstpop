/**/
    #ifndef MMATRIX_H
    #define MMATRIX_H

    #include <tnt/tnt_array2d.h>

    using namespace TNT;

    class MMatrix : public Array2D<double>
    {
        public:
            MMatrix(int x, int y);
            virtual ~MMatrix();
            inline Array2D<double> & operator=(const Array2D<double> &A);
            void removerow(int i);
            void removecolumn(int i);
            void addrow();
            void addcolumn();
        protected:
        private:
    };

    #endif // MMATRIX_H
/**/
