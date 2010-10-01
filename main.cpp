//------------standard--------------//
/*-----
/usr/include/c++/4.4
/usr/include/c++/4.4/tr1
/usr/include
----*/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <vector>
#include <bitset>
#include <algorithm>
#include <numeric> // for accumulate
#include <iterator>

//-------------openFST--------------//
/*--------
http://openfst.org/
statically linked shared libraries for openfst
/usr/local/lib/libfst.so
/usr/local/lib/libfstfar.so
/usr/local/lib/libfstcompact.so

linker search directories for dynamically linked libraries
/usr/local/lib
/usr/local/lib/fst
----------*/
#include <fst/fstlib.h>
#include <fst/compact-fst.h>
#include <fst/extensions/far/farlib.h>

//--------------TNT----------------//
/*--------
http://math.nist.gov/tnt/index.html
pure template library
headers at /usr/include/tnt
--------- */
#include <tnt/tnt_array2d.h>

//--------------GSL----------------//
/*-------------
http://www.gnu.org/software/gsl/manual/html_node/
statically linked shared libraries for gsl
/usr/lib/libgsl.so
/usr/lib/libgslcblas.so

headers at /usr/include/gsl
--------*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define N 50 // set population size
#define GEN 100 // set number of generations

using namespace fst;
using namespace TNT;
using namespace std;

//-------initialize global variables for gsl RNG-----------//
const gsl_rng_type * T;
gsl_rng * r;


//------generate n integers between 0 and range---//
vector<int> gslrandgen(int n, int range)

{
    int i;
    vector<int> uout (n);
    //uout = new int [n];
    for (i = 0; i < n; i++)
     {
       uout[i] = gsl_rng_uniform_int (r, range+1);
       //printf ("%d, ", uout[i]);
     }

    return uout;
    //gsl_rng_free (r);
}

//-------generate a random 2-state FST------------//
StdVectorFst randfst()
{
    //bool combo[]={0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,1,1,0,1,0,0,0,1,0,0,1,1,0,1,0};
    //vector<bool> arcbit; arcbit.assign(combo, combo+36);
    bitset<36> arcbit (string("000000010010010001010110100010011010"));
    vector<int> rout = gslrandgen(2, 8);
    vector<int> sout = gslrandgen(4, 1);

    //bitset<64> arcbit (string("0000000100100011010001010110011110001001101010111100110111101111"));
    //vector<int> rout = gslrandgen(4, 15);

    StdVectorFst Ra; // A vector FST is a general mutable FST
    // Adds state 0 to the initially empty FST and make it the start state.
    Ra.AddState();   // 1st state will be state 0 (returned by AddState)
    Ra.AddState();   // 2nd state

    //-------set start and final states randomly to 0 or 1 ------//
    vector<int> rfin = gslrandgen(2, 1);
    if (rfin[0]) {Ra.SetStart(0);}   else {Ra.SetStart(1);}
    if (rfin[1]) {Ra.SetFinal(1,0);} else {Ra.SetFinal(0,0);} // SetFinal (StateID, weight)

    // Adds arcs exiting state 0 and returning to state 0
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[0]*4+0]) { Ra.AddArc(0, StdArc(1, 1, 0, sout[0])); } // 1st arg is src state ID
    if (arcbit[rout[0]*4+1]) { Ra.AddArc(0, StdArc(1, 2, 0, sout[0])); }
    if (arcbit[rout[0]*4+2]) { Ra.AddArc(0, StdArc(2, 1, 0, sout[1])); }
    if (arcbit[rout[0]*4+3]) { Ra.AddArc(0, StdArc(2, 2, 0, sout[1])); }


    // Adds arcs exiting state 1 and returning to state 1
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[1]*4+0]) { Ra.AddArc(1, StdArc(1, 1, 0, sout[2])); } // 1st arg is src state ID
    if (arcbit[rout[1]*4+1]) { Ra.AddArc(1, StdArc(1, 2, 0, sout[2])); }
    if (arcbit[rout[1]*4+2]) { Ra.AddArc(1, StdArc(2, 1, 0, sout[3])); }
    if (arcbit[rout[1]*4+3]) { Ra.AddArc(1, StdArc(2, 2, 0, sout[3])); }

/*    // Adds arcs exiting state 0 and entering state 1
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[2]*4+0]) { Ra.AddArc(0, StdArc(1, 1, 0, 1)); } // 1st arg is src state ID
    if (arcbit[rout[2]*4+1]) { Ra.AddArc(0, StdArc(1, 2, 0, 1)); }
    if (arcbit[rout[2]*4+2]) { Ra.AddArc(0, StdArc(2, 1, 0, 1)); }
    if (arcbit[rout[2]*4+3]) { Ra.AddArc(0, StdArc(2, 2, 0, 1)); }

    // Adds arcs exiting state 1 and returning to state 1
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[3]*4+0]) { Ra.AddArc(1, StdArc(1, 1, 0, 0)); } // 1st arg is src state ID
    if (arcbit[rout[3]*4+1]) { Ra.AddArc(1, StdArc(1, 2, 0, 0)); }
    if (arcbit[rout[3]*4+2]) { Ra.AddArc(1, StdArc(2, 1, 0, 0)); }
    if (arcbit[rout[3]*4+3]) { Ra.AddArc(1, StdArc(2, 2, 0, 0)); }
*/

    bool stateconnect = 0;
    int i = Ra.Start();
    for (ArcIterator<StdVectorFst> aiter(Ra, i); !aiter.Done(); aiter.Next())
    {
      const StdArc &arc = aiter.Value();
      int ds = arc.nextstate; // store the destination state
      if (i!=ds){stateconnect=1;}
    }

    if (Verify(Ra) && stateconnect) {return Ra;} else {return randfst();}
}

/*template<class Arc>
int GetStates(const Fst<Arc> &fst) {

    // Count states
    int ns = 0;
    for (StateIterator< Fst<Arc> > siter(fst);
        !siter.Done();
        siter.Next())
        ++ns;
    return ns;
}*/

/*double scomplexity(StdVectorFst C)
{
    ArcSort(&C, StdILabelCompare());
    int totst = C.NumStates();
    Array2D<double> T(totst, totst);

    for (int i=0; i<=totst ; i++)
    {
        for (int l=1; l<=2 ; l++)
        {
            Matcher matcher(C, MATCH_INPUT);
            matcher.SetState(i); // i marks the current state
            if (matcher.Find(l))
            for (; !matcher.Done(); matcher.Next())
            {
            const StdArc &arc = matcher.Value();
            int ds = arc.nextstate; // store the destination state
            T[i][ds] = 1; // need to row normalize this Array2D before using to calculate eigenvectors for stat comp
            }
        }
    }
}
*/

//double ncomplexity(list<pair<StdVectorFst, int> > U, vector<int> PI) // Unique, Population Index
// intending to pass (popID, nVT)


int main()
{
// uses MT19937 generator of Makoto Matsumoto and Takuji Nishimura RNG by defualt
// Mersenne prime period of 2^19937 - 1 (about 10^6000) and is equi-distributed in 623 dimensions.
// http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_ran_discrete_t * g;
//---------do not call gslrandgen or randfst before this point-----//


    int fstnum = 3;
    string fstnames[fstnum];

    fstnames[0]="T1";
    fstnames[1]="T2";
    fstnames[2]="T3";
//------------------Generate [population] of FSTs--------------------//
    vector<StdVectorFst> VT (N); // Vector of Transducers (VT) is a container for the population
    //vector<int> nVT (N);
    //vector<int> randfstind(0);
    vector<int> ran_popdist (3);
    //Array2D<bool> popID(N,N);
    //map<int, int> popID;
    //map<int, int>::iterator it;
    vector<int> popdist; // container for population type distribution
    vector<double> popfreq; // container for population type frequency = popdist/N
    double psize = N; // double version of N for computing popfreq

    typedef list< pair<StdVectorFst, int> > FSTcatalog;
    FSTcatalog popID;
    FSTcatalog::iterator it;

    typedef vector< Array2D<double> > MatrixGroup;
    MatrixGroup intxnNet;

    for (int n=0; n<N; n++)
    {
        VT[n] = randfst();
    }

    //determine the number of unique individuals in the population

    popID.push_back( make_pair(VT[0],1) );
    for (int i=1; i<N; i++)
    {
        bool IDswitch = 1;
        for (it=popID.begin(); it!=popID.end() ; it++)
        {
            StdVectorFst ProT = (*it).first; //ProT = prototype from list

            if (RandEquivalent(VT[i],ProT,10,0))
                {
                    (*it).second = (*it).second + 1;
                    IDswitch = 0;
                    break;
                }
        }
        if (IDswitch)
        {
            popID.push_back( make_pair(VT[i],1) );
        }
    }
    cout << "popID size: " << (int) popID.size() << "\n\n";


    for (it=popID.begin(); it!=popID.end() ; it++)
    {
        popdist.push_back((*it).second);
        popfreq.push_back((*it).second/psize);
    }
    ostream_iterator< int > output( cout, " " );
    ostream_iterator< double > outd( cout, " " );
    copy(popdist.begin(),popdist.end(), output); cout << "\n\n\n";
    copy(popfreq.begin(),popfreq.end(), outd);cout << "\n\n\n";

    Array2D< double > ZeroArray(popID.size(),popID.size(), 0.0);

    for (int i=0; i<popID.size(); i++)
    {
        intxnNet.push_back(ZeroArray);
    }
    cout << "number of matrices in intxnNet: " << (int) intxnNet.size() << endl;

//-------------------Select two FSTs at random for composition------------//
    StdVectorFst compres; // Container for composition result
    StdVectorFst result; // Container for minimized composition result.

    StdVectorFst T1; int T1type; double T1freq;
    StdVectorFst T2; int T2type; double T2freq;
    int docounter = 0;
    int dolimit = 100;

    do
    {
        // generate 3 random integers from 0 to N-1
        // 1: first machine in composition
        // 2: second machine in compositon
        // 3: machine scheduled for replacement
        g = gsl_ran_discrete_preproc (popID.size(), &popfreq[0]); // can pass &vector[0] to function expecting an array

//??????
        for (int i=0; i<3; i++) {ran_popdist[i]= gsl_ran_discrete (r, g);} // will this go from 0 to popID.size()??

        it=popID.begin();
        advance(it,ran_popdist[0]);
            T1 = (*it).first;
            T1type = ran_popdist[0];
            T1freq = (*it).second/psize;

        it=popID.begin();
        advance(it,ran_popdist[1]);
            T2 = (*it).first;
            T2type = ran_popdist[1];
            T2freq = (*it).second/psize;

        T1.Write("onestate/T1.fst");
        T2.Write("onestate/T2.fst");

        // The FSTs must be sorted along the dimensions they will be joined.
        // In fact, only one needs to be so sorted.
        // This could have instead been done for "model.fst" when it was created.
        ArcSort(&T1, StdOLabelCompare());
        ArcSort(&T2, StdILabelCompare());

        // Create the composed FST.
        Compose(T1, T2, &result);
        //Minimize(&compres, &result);

        //cout << "number of states in result:" << GetStates(result) << endl;
        //cout << "result start state is: " << result.Start() << endl;

        result.Write("onestate/T3.fst");

        //cout << "composition result is:" << Verify(result) << endl;

        //prevent infinite loop
        if (docounter > dolimit)
        {
            cout << "# of composition failures exceeded limit: " << dolimit << endl;
            break;
        }
        docounter++;

    } while ((Verify(result) == 0) | (result.Start() == -1));

    //int d = randfstind[2]; // index of machine scheduled for replacement
    int d = ran_popdist[2]; // index of machine type scheduled for removal

    bool IDswitch = 1;
    int counter = 0;
    for (it=popID.begin(); it!=popID.end() ; it++)
    {
        StdVectorFst ProT = (*it).first; //ProT = prototype from list

        if (RandEquivalent(result,ProT,10,0))
            {
                (*it).second = (*it).second + 1;
                IDswitch = 0;
                //need to extract new population frequency vector from popID and update T1freq/T2freq
                popfreq[counter] = (*it).second/psize;
                if (counter==T1type) {T1freq = (T1freq*psize+1)/N;}
                if (counter==T2type) {T2freq = (T2freq*psize+1)/N;}

                intxnNet[counter][T1type][T2type]= T1freq*T2freq;
                //need to iterate through row "counter" and column "counter" of each matrix and update frequencies

                break;
            }
        counter++;
    }
    if (IDswitch)
    {
        popID.push_back( make_pair(result,1) );
        popfreq.push_back(1/psize);
        //need to add a row and a column to each matrix
        //need to add a row and a column to ZeroArray
        intxnNet.push_back(ZeroArray);
        intxnNet[intxnNet.size()][T1type][T2type]= T1freq*T2freq;

    }


    it=popID.begin();
    advance(it,d);
        if((*it).second==1)
        {
            popID.erase(it);
       //     intxnNet.erase(intxnNet.begin()+d-1); //segmentation fault
       //     popfreq.erase(popfreq.begin()+d-1);
            //need to erase row d and column d from each matrix in intxnNet

        }else {
            (*it).second = (*it).second - 1;
            //extract new population frequency vector from popID
       //     popfreq[d] = (*it).second/psize;
            //need to iterate through row d and column d of each matrix and update frequencies
        }


 /*   int counter = 0;
    for (it=popID.begin(); it!=popID.end() ; it++)
    {
        if (d==counter)
        {
            if((*it).second==1)
            {
                popID.erase(it);
            }else {
                (*it).second = (*it).second - 1;
            }
        }
    }
*/


    cout << "number of composition failures: " << (docounter-1) << endl;
    cout << "fst eliminated is #: " << d << endl;
    cout << "T1 start state is " << T1.Start() << endl;
    cout << "T2 start state is " << T2.Start() << endl;
    cout << "result start state is " << result.Start() << endl;


    //cout << "size of randfstind: " << (int) randfstind.size();
//---------compute and store interaction network complexity-------//



//---------print interaction network---------------//
    int cmdtst = system(NULL);
    if ( cmdtst != 0 )
    {
        for ( int i=0; i<fstnum; i++ )
        {

            string sh1="fstdraw onestate/.fst onestate/.dot";
            string sh2="dot -Tps onestate/.dot > onestate/.ps";
            string sh3="evince onestate/.ps &";
            string fstname = fstnames[i];
            sh1.insert(17,fstname); sh1.insert(33,fstname);
            sh2.insert(18,fstname); sh2.insert(36,fstname);
            sh3.insert(16,fstname);
            system(sh1.c_str()); system(sh2.c_str()); system(sh3.c_str());

        }

    } else { cout << "No command interpreter available \n"; };

    gsl_rng_free (r);
    gsl_ran_discrete_free (g);
    return 0;

}
