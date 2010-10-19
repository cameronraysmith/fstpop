//------------standard--------------//
/*-----
/usr/include/c++/4.4
/usr/include/c++/4.4/tr1
/usr/include
----*/
#include <iostream>
#include <fstream>
#include <ostream>
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
//#include <tnt/tnt_array2d.h>

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
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>

//------------classes------------//
#include "FSTcatalog.h"
//#include "MMatrix.h"

#define POP 100 // set population size
#define GEN 10000 // set number of generations

using namespace fst;
//using namespace TNT;
using namespace std;


//-------initialize global variables -----------//
    //---for GSL global RNG
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

    // Adds arcs exiting state 0
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[0]*4+0]) { Ra.AddArc(0, StdArc(1, 1, 0, sout[0])); } // 1st arg is src state ID
    if (arcbit[rout[0]*4+1]) { Ra.AddArc(0, StdArc(1, 2, 0, sout[0])); }
    if (arcbit[rout[0]*4+2]) { Ra.AddArc(0, StdArc(2, 1, 0, sout[1])); }
    if (arcbit[rout[0]*4+3]) { Ra.AddArc(0, StdArc(2, 2, 0, sout[1])); }


    // Adds arcs exiting state 1
    // Arc constructor args: ilabel, olabel, weight, dest state ID.
    if (arcbit[rout[1]*4+0]) { Ra.AddArc(1, StdArc(1, 1, 0, sout[2])); } // 1st arg is src state ID
    if (arcbit[rout[1]*4+1]) { Ra.AddArc(1, StdArc(1, 2, 0, sout[2])); }
    if (arcbit[rout[1]*4+2]) { Ra.AddArc(1, StdArc(2, 1, 0, sout[3])); }
    if (arcbit[rout[1]*4+3]) { Ra.AddArc(1, StdArc(2, 2, 0, sout[3])); }

    bool stateconnect = 0;
    bool tofinal = 0;
    for (int i = 0; i<=1; i++)
    {
        for (ArcIterator<StdVectorFst> aiter(Ra, i); !aiter.Done(); aiter.Next())
            {
              const StdArc &arc = aiter.Value();
              int ds = arc.nextstate; // store the destination state
              if (ds==rfin[1]) {tofinal = 1;}
              if (i!=ds){stateconnect=1;}
            }
    }
    if ((Verify(Ra) && stateconnect && tofinal)) { return Ra;} else {randfst(); return Ra;}
}

// this function prints progress bar
// pre: percentage complete
// post: progress bar
void progressbar(int percent)
{
    //static stringstream bars;
    static int z = 0, w = 0, q = 100;
    string slash[4];
    slash[0] = "\\";
    slash[1] = "-";
    slash[2] = "/";
    slash[3] = "|";
    //bars << "|";
    cout << "\r"; // carriage return back to beginning of line
    cout << "[" ;
    for (int i=0; i<w; i++) {cout << "|";} //<< bars.str()
    cout << slash[z];
    for (int i=0; i<q; i++) {cout << " ";}
    cout << "]" << percent << " %"; // print the bars and percentage
    cout.flush();
    z++; w++; q--; // increment to make the slash appear to rotate
    if(z == 4)
    z = 0; // reset slash animation
}

int main()
{

// uses MT19937 generator of Makoto Matsumoto and Takuji Nishimura RNG by defualt
// Mersenne prime period of 2^19937 - 1 (about 10^6000) and is equi-distributed in 623 dimensions.
// http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    gsl_rng_default_seed = 0;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_ran_discrete_t * g;
//---------do not call gslrandgen or randfst before this point-----//

    // set probability of generating a random FST (phi_in)
    // and probability of generating an FST via composition (1-phi_in)
    double phi_in = 0.5;

//---initiate log files
    ofstream fstpoplog;
    fstpoplog.open("fstpop.log");

    stringstream foutname;
    foutname << "complexity_" << phi_in << "_" << gsl_rng_default_seed << ".dat";
    string s = foutname.str();
    ofstream ofile;
    ofile.open(s.c_str());


//------------------Generate [population] of FSTs--------------------//
    vector<StdVectorFst> VT (POP); // Vector of Transducers (VT) is a container for the population

    vector<int> ran_popdist (3);

    for (int n=0; n<POP; n++)
    {
        VT[n] = randfst();
    }

    //determine the number of unique individuals in the population...in constructor of FSTcatalog

    FSTcatalog popInfo(VT);

//-------------------Select two FSTs at random for composition------------//
  //  StdVectorFst compres; // Container for composition result
    StdVectorFst result; // Container for minimized composition result.

    StdVectorFst T1; int T1type; //double T1freq;
    StdVectorFst T2; int T2type; //double T2freq;
    int dolimit = 100;
    int docounter = 0;
    vector<double> CmuG;
    vector<double> avgCmu;
    FSTlist::iterator it;

    for(int x=0; x<GEN; x++)
    {
        if(!(gsl_ran_bernoulli(r, phi_in)))
        {
            docounter = 0;
            do
            {
                // generate 3 random integers from 0 to N-1
                // 1: first machine in composition
                // 2: second machine in compositon
                // 3: machine scheduled for replacement
                g = gsl_ran_discrete_preproc (popInfo.popID.size(), &popInfo.popdist[0]); // can pass &vector[0] to function expecting an array
                for (int i=0; i<3; i++) {ran_popdist[i]= gsl_ran_discrete (r, g);} // will this go from 0 to popID.size()??

                    it=popInfo.popID.begin();
                    advance(it,ran_popdist[0]);
                        T1 = (*it).first;
                        T1type = ran_popdist[0];
                        //T1freq = (*it).second/psize;

                    it=popInfo.popID.begin();
                    advance(it,ran_popdist[1]);
                        T2 = (*it).first;
                        T2type = ran_popdist[1];
                        //T2freq = (*it).second/psize;

                    fstpoplog << "T1 type is: " << T1type << ", T2 type is: " << T2type << endl;
                    //T1.Write("onestate/T1.fst");
                    //T2.Write("onestate/T2.fst");

                    // The FSTs must be sorted along the dimensions they will be joined.
                    // In fact, only one needs to be so sorted.
                    // This could have instead been done for "model.fst" when it was created.
                    ArcSort(&T1, StdOLabelCompare());
                    ArcSort(&T2, StdILabelCompare());

                    // Create the composed FST.
                    Compose(T1, T2, &result);

                    fstpoplog << "number of states in result:" << result.NumStates() << endl;
                    //fstpoplog << "result start state is: " << result.Start() << endl;

                    //result.Write("onestate/T3.fst");



                //------prevent infinite loop------

                if (docounter > dolimit)
                {
                    fstpoplog << "# of composition failures exceeded limit: " << dolimit << endl;
                    break;
                }
                docounter++;

            } while ((Verify(result) == 0) | (result.Start() == -1));
        } else {
            g = gsl_ran_discrete_preproc (popInfo.popID.size(), &popInfo.popdist[0]); // can pass &vector[0] to function expecting an array
            for (int i=0; i<3; i++) {ran_popdist[i]= gsl_ran_discrete (r, g);} // will this go from 0 to popID.size()??

            result = randfst();
            T1type = -1;
            T2type = -1;
            fstpoplog << "result is a random 2-state machine" << endl;
        }

        int d = ran_popdist[2]; // index of machine type scheduled for removal
        fstpoplog << "type index scheduled for removal: " << d << endl;
        popInfo.update(result, d, T1type, T2type);

        //CmuG.push_back(popInfo.ncomplexity());
        //avgCmu.push_back(popInfo.scomplexity());
        ofile << popInfo.ncomplexity() << "\t" << popInfo.scomplexity() << "\n";


        double gennum = GEN;
        int genmod = int(gennum/100);

        if (x % genmod == 0)
        {
            int percentcomplete = int(x/gennum*100);
            progressbar(percentcomplete);
        }

    }

    gsl_rng_free (r);
    gsl_ran_discrete_free (g);
    fstpoplog.close();
    ofile.close();
    return 0;
}
