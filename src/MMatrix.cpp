/**/
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

    Array2D<double> & MMatrix::operator=(const Array2D<double> &A)
    {
        return ref(A);
    }

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
        M = out.copy();
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
        M = out.copy();
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

        M = out.copy();
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

        M = out.copy();
    }
/**/
