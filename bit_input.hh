#ifndef BIT_INPUT_HH_
#define BIT_INPUT_HH_

#include <vector>
#include <string>
#include <cstdio>

using namespace std;

class BitInput {

   public:
      BitInput(std::string filename);
      ~BitInput();
      int input_bit();

   private:

      int buffer;		/* Bits buffered for output                 */
      int bits_to_go;		/* Number of bits free in buffer            */
      FILE* outputFile;
};

#endif
