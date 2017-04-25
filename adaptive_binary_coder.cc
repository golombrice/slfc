#include "adaptive_binary_coder.hh"
#include <iostream>

AdaptiveBinaryCoder::AdaptiveBinaryCoder(std::string filename, int k1, int N1, bool decode): high(Top_value), low(0), bits_to_follow(0), decode_(decode), value(0)
{   
  if(decode) {
   bitinputter = new BitInput(filename);
  }
  else {
   bitoutputter = new BitOutput(filename);
  }

  initialize_binary(k1,N1);

   if( decode )
    for (int i = 1; i<=Code_value_bits; i++) {	/* code value.              */
        value = 2*value+bitinputter->input_bit();
    }
}

void AdaptiveBinaryCoder::initialize_binary(int k,int N) {
   if( freq.size() == 0 ) {
     freq.push_back(0);
     freq.push_back(k);
     freq.push_back(N-k);
     cum_freq.push_back(N);
     cum_freq.push_back(k);
     cum_freq.push_back(0);
   } else {
     freq[0]=(0);
     freq[1]=(k);
     freq[2]=(N-k);
     cum_freq[0]=(N);
     cum_freq[1]=(k);
     cum_freq[2]=(0);
   }
}

int AdaptiveBinaryCoder::decode_symbol()
{
    int symbol = 0;			/* Symbol decoded                           */
    long range = (long)(high-low)+1;
    int cum = 					/* Find cum freq for value. */
      (((long)(value-low)+1)*cum_freq[0]-1)/range;
    for (symbol = 1; cum_freq[symbol]>cum; symbol++) ; /* Then find symbol. */
    high = low +				/* Narrow the code region   */
      (range*cum_freq[symbol-1])/cum_freq[0]-1;	/* to that allotted to this */
    low = low + 				/* symbol.                  */
      (range*cum_freq[symbol])/cum_freq[0];
    while(true) {					/* Loop to get rid of bits. */
        if (high<Half) {
            /* nothing */			/* Expand low half.         */
        } 
        else if (low>=Half) {			/* Expand high half.        */
            value -= Half;
            low -= Half;			/* Subtract offset to top.  */
            high -= Half;
        }
        else if (low>=First_qtr			/* Expand middle half.      */
              && high<Third_qtr) {
            value -= First_qtr;
            low -= First_qtr;			/* Subtract offset to middle*/
            high -= First_qtr;
        }
        else break;				/* Otherwise exit loop.     */
        low = 2*low;
        high = 2*high+1;			/* Scale up code range.     */
        value = 2*value+bitinputter->input_bit();		/* Move in next input bit.  */
    }

  updateFreqz(symbol);

  return symbol;
}


void AdaptiveBinaryCoder::encode_symbol(int symbol)
{
   // ENCODE
   /* Size of the current code region          */
   long range = (long)(high-low)+1;
   high = low +				/* Narrow the code region   */
      (range*cum_freq[symbol-1])/cum_freq[0]-1;	/* to that allotted to this */
   low = low + 				/* symbol.                  */
      (range*cum_freq[symbol])/cum_freq[0];
   if( cum_freq[symbol-1] == cum_freq[symbol])
      cout << "ermor" << endl;
   while(true) {					/* Loop to output bits.     */
      if (high<Half) {
         bit_plus_follow(0);			/* Output 0 if in low half. */
      } 
      else if (low>=Half) {			/* Output 1 if in high half.*/
         bit_plus_follow(1);
         low -= Half;
         high -= Half;			/* Subtract offset to top.  */
      }
      else if (low>=First_qtr			/* Output an opposite bit   */
            && high<Third_qtr) {		/* later if in middle half. */
         bits_to_follow += 1;
         low -= First_qtr;			/* Subtract offset to middle*/
         high -= First_qtr;
      }
      else break;				/* Otherwise exit loop.     */
      low = 2*low;
      high = 2*high+1;			/* Scale up code range.     */
   }


  updateFreqz(symbol);

}

void AdaptiveBinaryCoder::updateFreqz(int symbol) {
   freq[symbol] -= 1;				/* Increment the frequency  */
   int running_symbol = symbol;
   while (running_symbol>0) {				
      running_symbol -= 1;					
      cum_freq[running_symbol] -= 1;				/* Increment the frequency  */
   }
}

void AdaptiveBinaryCoder::bit_plus_follow(int bit)
{   
   bool estim = false;
   if (estim){
      //nrBiti++;
   }
   else {
      bitoutputter->output_bit(bit);				/* Output the bit.          */
   }
   while (bits_to_follow>0) {
      if (estim) {
         //nrBiti++;
      }
      else {
         bitoutputter->output_bit(!bit);			/* Output bits_to_follow    */
      }
      bits_to_follow -= 1;			/* opposite bits. Set       */
   }						/* bits_to_follow to zero.  */
}

AdaptiveBinaryCoder::~AdaptiveBinaryCoder() {
    if( decode_ ) {
      delete bitinputter;
    } else {
      bits_to_follow += 1;			/* Output two bits that     */
      if (low<First_qtr) 
        bit_plus_follow(0);	/* select the quarter that  */
      else 
        bit_plus_follow(1);			/* the current code range   */
      delete bitoutputter;
    }
}
