#ifndef BIT_OUTPUT_HH 
#define BIT_OUTPUT_HH

#include <vector>
#include <string>
#include <cstdio>

using namespace std;

class BitOutput {

   public:
      BitOutput(std::string filename);
      ~BitOutput();
      void output_bit(int bit);

   private:

      int buffer;		/* Bits buffered for output                 */
      int bits_to_go;		/* Number of bits free in buffer            */
      FILE* outputFile;
      int n_bits_;
};

#endif
