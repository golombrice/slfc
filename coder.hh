#ifndef CODER_HH
#define CODER_HH

#include <vector>
#include <string>
#include "gen_types.hh"
#include "image.hh"

using std::vector;
using std::string;

class Coder {

   public:

     Coder( 
        Image< uint8_t >* Ydecoded,
        vector< Image< uint8_t >* >* cwarpedOK,
        Image< int >* Crestr);
     ~Coder();

     void SetErrorImagep(Image<int>* im);
     Image<int>* GetErrorImagep();
     void SetErrorImagec(Image<int>* im);
     Image<int>* GetErrorImagec();

   protected:


     Image< uint8_t >* Ydecoded_;
     vector< Image< uint8_t >* >* cwarpedOK_;
     Image< int >* Crestr_;
     vector< vector < vector < int >* >* >* SPARSEW_;
     vector< vector < vector < int >* >* >* THETAW_;
     vector< vector< vector< double >* >* >* MINMAXW_;
     Image< int >* ErrorImagep_;
     Image< int >* ErrorImagec_;

     int N_regions_;
     int NC_;
     int NR_;

     double GeneralPredictionp(const int icomp, const int iir, const int iic, vector< int >& Sq);
     double WarpedPredictionp(const int icomp, const int iir, const int iic, vector< int >& Sq);

     int computeContext(const int icomp, const int iir, const int iic);
     bool fixBoundary(vector< int >& Sq, int iir, int iic, int iR);

     void compute_Nregions();

};
#endif
