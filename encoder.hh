#ifndef ENCODER_HH
#define ENCODER_HH

#include <vector>
#include <string>
#include "gen_types.hh"
#include "coder.hh"

using std::vector;
using std::string;

class Encoder: public Coder {

   public:

     Encoder( 
        Image< uint8_t >* Ydecoded,
        vector< Image< uint8_t >* >* cwarpedOK,
        Image< int >* Crestr);
     ~Encoder();

     void computeMinMaxW();
     void encode(string file);
     void find_predictors();
     void read_predictors(string filename);
     void writeMetadata(string predictor_fn, string mask_fn, string minmax_fn);
     vector< vector< int > > mergeRegions();

   private:

     static double fastOLS(vector< vector< uint8_t > >& A, vector< uint8_t >& d, int Ms, vector< int >* thetav, vector< int >* sparsev);
     static void fastOLS_looper(vector<int> inds, vector< vector < vector< uint8_t > > >& A_all, vector< vector< uint8_t > >& d_all, int Ms, vector< vector< int >* >* thetav, vector< vector< int >* >* sparsev);

};
#endif
