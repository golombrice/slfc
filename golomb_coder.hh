#ifndef GOLOMB_CODER_HH
#define GOLOMB_CODER_HH 

#include <vector>
#include <string>
#include "bit_output.hh"
#include "bit_input.hh"
#include "gen_types.hh"

using namespace std;

class GolombCoder {

   public:
      GolombCoder(std::string filename, bool decode);
      ~GolombCoder();

      void encode_symbols(vector< int >& symbols, int N_bits);
      void decode_symbols(vector< int >& symbols, int N_bits);
   private:

      bool decode_;
      BitOutput* bitoutputter;
      BitInput* bitinputter;
      int bit_counter;


};

#endif
