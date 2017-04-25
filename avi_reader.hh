#ifndef AVI_READER_HH
#define AVI_READER_HH 

#include <string>
#include <cstdio>
#include "image.hh"
using std::string;

class AVIReader{
   public:

      // Constructs an object reading an "Uncompressed AVI" (see e.g. Matlab's VideoWriter).
      // Takes the filename and size of one frame as input.
      AVIReader(string filename, int NR, int NC): NR_(NR), NC_(NC), is_avi_file_(false) {

        if(filename.at(filename.size()-1) == 'i' && filename.at(filename.size()-2) == 'v' && filename.at(filename.size()-3) == 'i') {
          is_avi_file_ = true;
        }

        if(is_avi_file_) {
           avi_file_ = fopen(filename.c_str(), "rb");
           rewind(avi_file_);
           char temp[4];
           long bytes = 0;
           while(fread(temp, sizeof(char), 4, avi_file_) == 4) {
              bytes += 4;
              if(temp[0]=='m' && temp[1]=='o' && temp[2] == 'v' && temp[3] == 'i') {
                 break;
              }
           }
           if(temp[0]=='m' && temp[1]=='o' && temp[2] == 'v' && temp[3] == 'i') {
              start_byte_ = bytes;
           }
           else {
              cout << "AVIReader: data not found!" << endl;
           }
        }
        else {
           avi_file_ = fopen(filename.c_str(), "rb");
        }
      }

      // This class deletes the images it creates.
      ~AVIReader() {
         fclose(avi_file_);
         for( int n = 0; n < images_.size(); ++n ) {
            delete images_[n];
         }
      }

      // Returns a pointer to n-th (indexing from zero) view (frame) in the sequence.
      Image< uint8_t>* get_view(int n) {

        if( is_avi_file_ ) {
         // start_byte_ denotes the starting of frames 
         // and one frame takes the amount of bytes given by the second term. 
         long offset = start_byte_+n*(NR_*(3*NC_+1)+8);
         fseek(avi_file_, offset, SEEK_SET);

         // Allocate new image and store the pointer.
         Image<uint8_t>* A_ij = new Image<uint8_t>(NR_,NC_,3);
         images_.push_back(A_ij);

         // Read 7 bytes of metadata before pixels.
         uint8_t one_element = 0;
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         fread(&one_element,1,1,avi_file_);
         // The order goes from the last row to the first, first column to last. 
         // Every row has an extra byte in the beginning.
         for( int iir = NR_-1; iir >= 0; --iir) {
            fread(&one_element,1,1,avi_file_);
            for( int iic = 0; iic < NC_; ++iic ) {
               fread(A_ij->getp(iir,iic,2), 1, 1, avi_file_);
               //cout << A_ij->get(iir,iic,2) << " ";
               fread(A_ij->getp(iir,iic,1), 1, 1, avi_file_);
               fread(A_ij->getp(iir,iic,0), 1, 1, avi_file_);
            }
         }
         fread(&one_element,1,1,avi_file_);

         return A_ij;
        }
        else {
         long offset = n*3*NR_*NC_;
         fseek(avi_file_, offset, SEEK_SET);

         // Allocate new image and store the pointer.
         Image<uint8_t>* A_ij = new Image<uint8_t>(NR_,NC_,3);
         images_.push_back(A_ij);
         for( int iir = 0; iir < NR_; ++iir) {
            for( int iic = 0; iic < NC_; ++iic ) {
               fread(A_ij->getp(iir,iic,0), 1, 1, avi_file_);
               fread(A_ij->getp(iir,iic,1), 1, 1, avi_file_);
               fread(A_ij->getp(iir,iic,2), 1, 1, avi_file_);
            }
         }

         return A_ij;
        }
      }
      
   private:
      vector< Image< uint8_t >* > images_;
      long start_byte_;
      bool is_avi_file_;
      FILE* avi_file_;
      int NR_;
      int NC_;
};

#endif
