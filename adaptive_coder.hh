#ifndef ADAPTIVE_CODER_HH
#define ADAPTIVE_CODER_HH 

#include <vector>
#include <string>
#include "bit_output.hh"
#include "bit_input.hh"
#include "gen_types.hh"

using std::vector;
using std::string;

class AdaptiveCoder {

   public:
      AdaptiveCoder(std::string filename, vector< vector < vector < double >* >* >* MINMAXW, bool decode);
      ~AdaptiveCoder();

      void encode_symbol(int symbol, int iW, int context);
      int decode_symbol(int iW,int icon);
   private:

      bool decode_;
      vector< vector< vector<uint32_t> > > freq;
      vector< vector< vector<uint32_t> > > cum_freq;
      BitOutput* bitoutputter;
      BitInput* bitinputter;
      vector< vector< vector< double >* >* >* MINMAXW_;
      uint64_t high;
      uint64_t low;
      uint64_t value;
      int bits_to_follow;

      void bit_plus_follow(int bit);
      void updateFreqz(int symbol, int iW, int context);
      void initialize_with_ones();

};

#endif
