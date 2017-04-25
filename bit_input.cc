#include <iostream>
#include "bit_input.hh"

BitInput::BitInput(std::string filename): buffer(0), bits_to_go(0)
{   
   outputFile = fopen(filename.c_str(), "rb");
}


BitInput::~BitInput()
{   
}

int BitInput::input_bit()
{   int t;
    bits_to_go -= 1;
    if (bits_to_go<0) {				/* Read the next byte if no */
        if(fread(&buffer,1,1,outputFile) != 1) {
          //printf("Can not read!\n");
        }
		//buffer = getc(stdin);			/* bits are left in the     */
        bits_to_go = 7;				/* buffer. Return anything  */
    }			   			/* after end-of-file.       */
    t = buffer&1;
    buffer >>= 1;				/* Return the next bit from */
    return t;					/* the bottom of the byte.  */
}
