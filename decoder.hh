#ifndef DECODER_HH
#define DECODER_HH

#include <vector>
#include <string>
#include "coder.hh"
#include "gen_types.hh"

using std::vector;
using std::string;

class Decoder: public Coder {

   public:

     Decoder( 
        Image< uint8_t >* Ydecoded,
        vector< Image< uint8_t >* >* cwarpedOK,
        Image< int >* Crestr);
     ~Decoder();

     void decode(string file);
     void readMetadata(string predictor_fn, string mask_fn, string minmax_fn);

   private:

};
#endif
