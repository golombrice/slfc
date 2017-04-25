#ifndef ADAPTIVE_BINARY_CODER_HH
#define ADAPTIVE_BINARY_CODER_HH 

#include <vector>
#include <string>
#include "bit_output.hh"
#include "bit_input.hh"
#include "gen_types.hh"

using namespace std;

class AdaptiveBinaryCoder {

   public:
      AdaptiveBinaryCoder(std::string filename, int k, int N, bool decode);
      ~AdaptiveBinaryCoder();

      void encode_symbol(int symbol);
      int decode_symbol();
      void initialize_binary(int k, int N);

   private:

      bool decode_;
      vector< int > freq;
      vector< int > cum_freq;
      BitOutput* bitoutputter;
      BitInput* bitinputter;
      long high;
      long low;
      long value;
      int bits_to_follow;

      void bit_plus_follow(int bit);
      void updateFreqz(int symbol);

};

#endif
