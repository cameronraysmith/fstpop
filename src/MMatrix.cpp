    #include "MMatrix.h"

    MMatrix::MMatrix(int x, int y) : Array2D<double> (x,y) // use constructor from Array2D
    {
        MMatrix &M = *this;
        for (int i=0; i<x; i++)
        {
            for (int j=0; j<y; j++)
            {
                M[i][j]=0;
            }
        }
    }

    MMatrix::~MMatrix()
    {
        //dtor
    }
/*
    MMatrix & MMatrix::operator=(const MMatrix &A)
    {
        ref(*(dynamic_cast<const Array2D<double>* > (&A)));
        return *this;
    }

    MMatrix & MMatrix::operator=(const Array2D<double> &A)
    {
        ref(A);
        return *this;
    }


    template <class T>
    Array2D<T> & Array2D<T>::ref(const Array2D<T> &A)
    {
        if (this != &A)
        {
            v_ = A.v_;
            data_ = A.data_;
            m_ = A.m_;
            n_ = A.n_;

        }
        return *this;
    }
*/
    void MMatrix::removerow(int i)
    {
        MMatrix &M = *this;

        int numrows = M.dim1()-1;
        MMatrix out(numrows, M.dim2());

        for(int j=0; j<i; j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                out[j][k] = M[j][k];
            }
        }

        for(int j=i; j<numrows; j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                out[j][k] = M[j+1][k];
            }
        }
        M = out;
    }

    void MMatrix::removecolumn(int i)
    {
        MMatrix &M = *this;

        int numcols = M.dim2()-1;
        MMatrix out(M.dim1(), numcols);

        for(int j=0; j<i; j++)
        {
            for (int k=0; k<M.dim1(); k++)
            {
                out[k][j] = M[k][j];
            }
        }

        for(int j=i; j<numcols; j++)
        {
            for (int k=0; k<M.dim1(); k++)
            {
                out[k][j] = M[k][j+1];
            }
        }
        M = out;
    }

    void MMatrix::addrow()
    {
        MMatrix &M = *this;

        MMatrix out(M.dim1()+1, M.dim2());

        for(int j=0; j<M.dim1(); j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                out[j][k] = M[j][k];
            }
        }

        M = out;
    }

    void MMatrix::addcolumn()
    {
        MMatrix &M = *this;

        MMatrix out(M.dim1(), M.dim2()+1);

        for(int j=0; j<M.dim1(); j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                out[j][k] = M[j][k];
            }
        }

        M = out;
    }

    void MMatrix::rownorm()
    {
        MMatrix &M = *this;

        MMatrix out(M.dim1(), M.dim2());

        std::vector<double> rowsum;

        for(int j=0; j<M.dim1(); j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                rowsum[j]+= M[j][k];
            }
        }

        for(int j=0; j<M.dim1(); j++)
        {
            for (int k=0; k<M.dim2(); k++)
            {
                out[j][k] = M[j][k]/rowsum[j];
            }
        }

        M = out;
    }
