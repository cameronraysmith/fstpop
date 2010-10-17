/**/
    #include "FSTcatalog.h"

    FSTcatalog::FSTcatalog(vector<StdVectorFst> V)
    {
        N=V.size(); // store population size

        FSTlist::iterator it;

        popID.push_back( make_pair(V[0],1) );
        for (unsigned int i=1; i<V.size(); i++)
        {
            bool IDswitch = 1;
            for (it=popID.begin(); it!=popID.end() ; it++)
            {
                StdVectorFst ProT = (*it).first; //ProT = prototype from list

                if (RandEquivalent(V[i],ProT,100,0))
                    {
                        (*it).second = (*it).second + 1;
                        IDswitch = 0;
                        break;
                    }
            }
            if (IDswitch)
            {
                popID.push_back( make_pair(V[i],1) );
            }
        }
        cout << "popID size: " << (int) popID.size() << "\n\n";

        //compute popdist
        for (it=popID.begin(); it!=popID.end() ; it++)
        {
            popdist.push_back((*it).second);
        }
        ostream_iterator< double > output( cout, " " );
        copy(popdist.begin(),popdist.end(), output); cout << "\n\n\n";

        //initialize the group of interaction matrices
       // MMatrix ZeroArray(popID.size(),popID.size());


        for (unsigned int i=0; i<popID.size(); i++)
        {
        MMatrix ZeroArray(popID.size(),popID.size());
        intxnNet.push_back(ZeroArray);
        }
        cout << "number of matrices in intxnNet: " << (int) intxnNet.size() << endl;
    }

    FSTcatalog::~FSTcatalog()
    {
        //dtor
    }

    void FSTcatalog::update(StdVectorFst result, int d, int T1type, int T2type)
    {
        FSTcatalog &C = *this; // corresponds to FSTcatalog popInfo in main.cpp
        FSTlist::iterator it;


        //--------update C.popID with the composition result --------//
        bool IDswitch = 1;
        int counter = 0;

        StdVectorFst ProT;

        for (it=C.popID.begin(); it!=C.popID.end() ; it++)
        {
            ProT = (*it).first; //ProT = prototype from list

            if (RandEquivalent(result,ProT,100,0))
                {
                    (*it).second = (*it).second + 1;
                    C.intxnNet[counter][T1type][T2type] = 1;

                    IDswitch = 0;
                    break;
                }
            counter++;
        }



        if (IDswitch)
        {
            C.popID.push_back( make_pair(result,1) );
            MMatrix ZeroArray(popID.size(),popID.size());
            C.intxnNet.push_back(ZeroArray);
            for (unsigned int i=0; i<C.intxnNet.size(); i++)
            {
                C.intxnNet[i].addrow();
                C.intxnNet[i].addcolumn();
            }
            int lastind = C.intxnNet.size()-1;
            C.intxnNet[lastind][T1type][T2type]= 1;
        }

        //-------update C.popID based upon the individual randomly selected for removal
        it=C.popID.begin();
        advance(it,d);
        if((*it).second==1)
        {
            C.popID.erase(it);

            C.intxnNet.erase(C.intxnNet.begin()+d);
            for (unsigned int i=0; i<C.intxnNet.size(); i++)
            {
                C.intxnNet[i].removerow(d);
                C.intxnNet[i].removecolumn(d);
            }
        }else {
            (*it).second = (*it).second - 1;
        }

        //-------update C.popdist----------------
        int popsum=0;
        C.popdist.clear();
        for (it=C.popID.begin(); it!=C.popID.end() ; it++)
        {
            C.popdist.push_back((*it).second);
            popsum += (*it).second;
        }
        cout << "population size is: " << popsum << endl;
    }

    double FSTcatalog::ncomplexity()
    {
        FSTcatalog &C = *this;
        double CmuG = 0;
        vector<double> vjk;
        double Vi;
        int count = 0;

        for (unsigned int i=0; i<C.intxnNet.size(); i++)
        {
            vjk.clear();
            Vi = 0;

            for (int j=0; j<C.intxnNet[0].dim1(); j++)
            {
                for (int k=0; k<C.intxnNet[0].dim2(); k++)
                {
                    if (C.intxnNet[i][j][k])
                    {
                        vjk.push_back(C.popdist[j]/C.N*C.popdist[k]/C.N);
                        Vi += C.popdist[j]/C.N*C.popdist[k]/C.N;
                        count++;
                        cout << "true indices from intxNet: " << i << ",  " << j << ",  " << k << endl;
                    }
                }
            }
            for(unsigned int l=0; l<vjk.size(); l++)
            {
                CmuG += -(vjk[l]/Vi*log2(vjk[l]/Vi));
            }
        }
        cout << "number of 1s in intxnNet: " << count << endl;
        cout << "interaction network complexity (CmuG): " << CmuG << endl;
        return CmuG;
    }

    double FSTcatalog::scomplexity()
    {
        FSTcatalog &C = *this;
        FSTlist::iterator it;
        double avgCmu = 0;

        for (it=C.popID.begin(); it!=C.popID.end() ; it++)
        {
            StdVectorFst F = (*it).first;
            ArcSort(&F, StdILabelCompare());
            int totst = F.NumStates();
            //cout << "totst is: " << totst << endl;
            gsl_matrix * T = gsl_matrix_calloc(totst, totst);

            if (totst < 1)
            {
                cout << "attempted to calculate complexity for machine with " << totst << " states" << endl;
            }

            for (int i=0; i<totst ; i++)
            {
                for (int l=1; l<=2 ; l++)
                {
                    Matcher<StdVectorFst> matcher(F, MATCH_INPUT);
                    matcher.SetState(i); // i marks the current state
                    if ((matcher.Find(l)))
                    {
                        for (; !matcher.Done(); matcher.Next())
                        {
                        const StdArc &arc = matcher.Value();
                        int ds = arc.nextstate; // store the destination state
                        gsl_matrix_set (T,i,ds,1); // need to row normalize this gsl_matrix before using to calculate eigenvectors for stat comp
                        }
                    }
                }
            }

            for (unsigned int i=0; i<(*T).size1; i++)
            {
                gsl_vector_view row = gsl_matrix_row(T, i);
                double rowsum = gsl_blas_dasum(&row.vector);
                if(rowsum!=0)
                {
                    gsl_vector_scale(&row.vector, 1/rowsum);
                }
            }

            // transpose transition matrix to compute left eigenvectors
            gsl_matrix_transpose(T);

            gsl_vector_complex *eval = gsl_vector_complex_alloc (totst);
            gsl_matrix_complex *evec = gsl_matrix_complex_alloc (totst, totst);

            gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (totst);

            gsl_eigen_nonsymmv (T, eval, evec, w);


            gsl_eigen_nonsymmv_free (w);

            gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

            vector<double> pvec;

            //   cout << "number of evals: " << eval.size << endl;
            //cout << "number of evals: " << (*eval).size << endl;

            for(unsigned int i=0; i<(*eval).size; i++)
            {
                gsl_complex tempeval = gsl_vector_complex_get(eval, i);
                //cout << "Re(tempeval) = " << tempeval.dat[0] << endl;
                //cout << "Im(tempeval) = " << tempeval.dat[1] << endl;

                if ((tempeval.dat[0])==1 && !tempeval.dat[1]) //selects for eigenvalue 1 + 0i // && !tempeval.dat[1]
                {
                   // printf("eval = %g + %gi\n",GSL_REAL(tempeval),GSL_IMAG(tempeval));
                    gsl_vector_complex * compe1vec = gsl_vector_complex_alloc ((*evec).size1);
                    gsl_vector * e1vec = gsl_vector_alloc ((*evec).size1);
                    gsl_matrix_complex_get_col(compe1vec, evec, i);
                    gsl_vector_view e1vecview = gsl_vector_complex_real(compe1vec);
                    gsl_vector_memcpy(e1vec, &e1vecview.vector);

                   // for(int i=0; i<(*e1vec).size; i++) printf("e1vec of %d = %f \n",i,gsl_vector_get(e1vec, i));

                    double vecsum = gsl_blas_dasum(e1vec);
                    if(vecsum!=0)
                    {
                        gsl_vector_scale(e1vec, 1/vecsum);
                    }

                    for(unsigned int j=0; j<(*e1vec).size; j++)
                    {
                        pvec.push_back(fabs(gsl_vector_get(e1vec, j)));
                        //cout << "prob: " << pvec[j] << endl;
                    }

                    gsl_vector_complex_free(compe1vec);
                    gsl_vector_free(e1vec);
                    break;
                }
            }

            double Cmu = 0;
            for(unsigned int l=0; l<pvec.size(); l++)
            {
                if(pvec[l]!=0.0)
                {
                    Cmu += -(pvec[l]*log2(pvec[l]));
                }
      //          printf("pvec of %d = %f\n", l, pvec[l]);
            }

            avgCmu += ((*it).second*Cmu)/(C.N);
      //      cout << "individual Cmu = : " << Cmu << endl;
      //      cout << "avg individual Cmu = : " << avgCmu << endl << endl;


            gsl_vector_complex_free(eval);
            gsl_matrix_complex_free(evec);
            gsl_matrix_free(T);
        }
        cout << "avg individual complexity (Cmu): " << avgCmu << endl << endl;
        return avgCmu;
    }

/**/
