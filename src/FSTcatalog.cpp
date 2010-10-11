/**/
    #include "FSTcatalog.h"

    FSTcatalog::FSTcatalog(vector<StdVectorFst> V)
    {
        FSTlist::iterator it;

        popID.push_back( make_pair(V[0],1) );
        for (int i=1; i<V.size(); i++)
        {
            bool IDswitch = 1;
            for (it=popID.begin(); it!=popID.end() ; it++)
            {
                StdVectorFst ProT = (*it).first; //ProT = prototype from list

                if (RandEquivalent(V[i],ProT,10,0))
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
        MMatrix ZeroArray(popID.size(),popID.size());

        for (unsigned int i=0; i<popID.size(); i++)
        {
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
        FSTcatalog &C = *this;
        FSTlist::iterator it;

        //--------update C.popID with the composition result --------//
        bool IDswitch = 1;
        int counter = 0;
        for (it=C.popID.begin(); it!=C.popID.end() ; it++)
        {
            StdVectorFst ProT = (*it).first; //ProT = prototype from list

            if (RandEquivalent(result,ProT,10,0))
                {
                    (*it).second = (*it).second + 1;

                    C.intxnNet[counter][T1type][T2type]= 1;

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
            for (int i=0; i<C.intxnNet.size(); i++)
            {
                C.intxnNet[i].addrow();
                C.intxnNet[i].addcolumn();
            }
            C.intxnNet[C.intxnNet.size()][T1type][T2type]= 1;
        }

        //-------update C.popID based upon the individual randomly selected for removal
        it=C.popID.begin();
        advance(it,d);
        if((*it).second==1)
        {
            C.popID.erase(it);

            C.intxnNet.erase(C.intxnNet.begin()+d);
            for (int i=0; i<C.intxnNet.size(); i++)
            {
                C.intxnNet[i].removerow(d);
                C.intxnNet[i].removecolumn(d);
            }
        }else {
            (*it).second = (*it).second - 1;
        }

    }

/**/
