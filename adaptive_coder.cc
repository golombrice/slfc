#include "adaptive_coder.hh"
#include <iostream>

AdaptiveCoder::AdaptiveCoder(std::string filename, vector< vector< vector < double >* >* >* MINMAXW, bool decode): high(Top_value), low(0), bits_to_follow(0), decode_(decode), value(0)
{   
  if(decode) {
   bitinputter = new BitInput(filename);
  }
  else {
   bitoutputter = new BitOutput(filename);
  }

   MINMAXW_ = MINMAXW;

   initialize_with_ones();

   if( decode )
    for (int i = 1; i<=Code_value_bits; i++) {	/* code value.              */
        value = 2*value+bitinputter->input_bit();
    }
}

void AdaptiveCoder::initialize_with_ones() {
   int N_symbols = 513;
   int N_regions = MINMAXW_->at(0)->size();
   int N_con = MINMAXW_->size();

   for( int i = 0; i < N_regions; ++i ) {
      vector< vector < uint32_t > > one_region;
      for( int j = 0; j < N_con; ++j ) {

         vector< uint32_t > one_context(N_symbols);
         one_region.push_back( one_context );
      }
      freq.push_back( one_region);
      cum_freq.push_back( one_region);
   }

   int vmax;
   for(int k=0;k<N_regions;k++)
   {
      for(int i=0; i < N_con; i++)
      {
         vmax = (int) (MINMAXW_->at(i)->at(k)->at(1)-MINMAXW_->at(i)->at(k)->at(0))+1;
         for (int j=0; j<=vmax; j++)
         {	
            freq[k][i][j] = 1;
            cum_freq[k][i][j] = vmax-j;  // cum_freq2[i][j] = My_no_of_symbols-j;
         }
         freq[k][i][0] = 0;
      }
   }
}

int AdaptiveCoder::decode_symbol(int iW,int icon)
{
    int symbol = 0;			/* Symbol decoded                           */
    uint64_t range = (uint64_t)(high-low)+1;
    int cum = 					/* Find cum freq for value. */
      (((uint64_t)(value-low)+1)*cum_freq[iW][icon][0]-1)/range;
    for (symbol = 1; cum_freq[iW][icon][symbol]>cum; symbol++) ; /* Then find symbol. */
    high = low +				/* Narrow the code region   */
      (range*cum_freq[iW][icon][symbol-1])/cum_freq[iW][icon][0]-1;	/* to that allotted to this */
    low = low + 				/* symbol.                  */
      (range*cum_freq[iW][icon][symbol])/cum_freq[iW][icon][0];
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

  updateFreqz(symbol,iW,icon);

  return symbol;
}


void AdaptiveCoder::encode_symbol(int symbol,int iW,int icon)
{
   // ENCODE
   /* Size of the current code region          */
   uint64_t range = (uint64_t)(high-low)+1;
   high = low +				/* Narrow the code region   */
      (range*cum_freq[iW][icon][symbol-1])/cum_freq[iW][icon][0]-1;	/* to that allotted to this */
   low = low + 				/* symbol.                  */
      (range*cum_freq[iW][icon][symbol])/cum_freq[iW][icon][0];
   if( cum_freq[iW][icon][symbol-1] == cum_freq[iW][icon][symbol])
      cout << "ermor" << endl;

   //cout << symbol << " " << high << " " << low << " " << range << endl;
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

  updateFreqz(symbol,iW,icon);

}

void AdaptiveCoder::updateFreqz(int symbol, int iW, int icon) {
   int vmax = 0;
   vmax = (int) (MINMAXW_->at(icon)->at(iW)->at(1)-MINMAXW_->at(icon)->at(iW)->at(0))+1;
   if ( (cum_freq[iW][icon][0]>=2000))
   //if ( ((icon==0) & (cum_freq[iW][icon][0]>=1000)) | ((icon>0) & (cum_freq[iW][icon][0]>= 2000)) )  // 16383  100
   {		                            
      int cum = 0;				         
      cum = 0;
      for (int i = vmax; i>=0; i--) 
      {	                           
         freq[iW][icon][i] = (freq[iW][icon][i]+1)/2;	
         cum_freq[iW][icon][i] = cum; 			
         cum += freq[iW][icon][i];
      }
   }

   freq[iW][icon][symbol] += 1;				/* Increment the frequency  */
   int running_symbol = symbol;
   while (running_symbol>0) {				
      running_symbol -= 1;					
      cum_freq[iW][icon][running_symbol] += 1;				/* Increment the frequency  */
   }
}

void AdaptiveCoder::bit_plus_follow(int bit)
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

AdaptiveCoder::~AdaptiveCoder() {
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
