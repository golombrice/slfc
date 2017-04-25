#include <iostream>
#include "bit_output.hh"

BitOutput::BitOutput(std::string filename): buffer(0), bits_to_go(8), n_bits_(0)
{   
   outputFile = fopen(filename.c_str(), "wb");
}


/* OUTPUT A BIT. */
void BitOutput::output_bit(int bit)
{   buffer >>= 1; if (bit) buffer |= 0x80;	/* Put bit in top of buffer.*/
    bits_to_go -= 1;
    if (bits_to_go==0) {			/* Output buffer if it is   */
        //putc(buffer,stdout);			/* now full.                */
        fwrite(&buffer,1,1,outputFile);
        n_bits_ += 8;

        bits_to_go = 8;
    }
}


/* FLUSH OUT THE LAST BITS. */
BitOutput::~BitOutput()
{   
	buffer = (buffer>>bits_to_go);
	fwrite(&buffer,1,1,outputFile);
   fclose(outputFile);
}
