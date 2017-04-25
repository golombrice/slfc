#ifndef IMAGE_HH
#define IMAGE_HH

#include <iostream>
#include <vector>
#include <stdint.h>
#include "gen_types.hh"

using namespace std;

template < class type > class Image {
  public:
    //template< class type >
    inline  Image(int height,int width,int nc): height_(height), width_(width), nc_(nc) {
      if( nc_ == 1 ) {
        r_ = new vector< type >(width*height);
      }
      if( nc_ == 3 ) {
        r_ = new vector< type >(width*height);
        g_ = new vector< type >(width*height);
        b_ = new vector< type >(width*height);
      }
    }
    //template< class  type >
    inline void set(int x,int y, int c, type val) {
      switch(c){
        case 0:
          r_->at(sub2ind(x,y)) = val;
          break; 
        case 1:
          g_->at(sub2ind(x,y)) = val;
          break; 
        case 2:
          b_->at(sub2ind(x,y)) = val;
          break; 
      }
    }


    //template< class  type >
    inline type get(int x,int y, int c) {
      switch(c){
        case 0:
          return r_->at(sub2ind(x,y));
        case 1:
          return g_->at(sub2ind(x,y));
        case 2:
          return b_->at(sub2ind(x,y));
      }
      return 0;
    }

    inline type* getp(int x,int y, int c) {
      switch(c){
        case 0:
          return &r_->operator[](sub2ind(x,y));
        case 1:
          return &g_->operator[](sub2ind(x,y));
        case 2:
          return &b_->operator[](sub2ind(x,y));
      }
      return 0;
    }

    //template< class  type >
    inline MatSize size()
    {
      MatSize size = {height_, width_};
      return size;
    }

    // Image(int height,int width,int nc);
    // void set(int x,int y, int c, type val);
    // type get(int x,int y, int c);
    // void set(int ind, int val);
    // MatSize size(); 
    inline void writeImage(string filename, string operation="wb") {
       FILE* decode_file = fopen(filename.c_str(), operation.c_str());

       if( nc_ == 1 ) {
         for(int i=0;i<height_;i++) {
           for(int l = 0; l < width_; ++l ) {
             type temp = r_->at(sub2ind(i,l));
             fwrite(&temp, sizeof(type), 1, decode_file);
           }
         }
       }
       if( nc_ == 3 ) {
          for(int i=0;i<height_;i++) {
             for(int l = 0; l < width_; ++l ) {
                type temp = r_->at(sub2ind(i,l));
                fwrite(&temp, sizeof(type), 1, decode_file);
                temp = g_->at(sub2ind(i,l));
                fwrite(&temp, sizeof(type), 1, decode_file);
                temp = b_->at(sub2ind(i,l));
                fwrite(&temp, sizeof(type), 1, decode_file);
             }
          }
       }
       fclose(decode_file);
    }

    inline void readImage(string filename) {
       int temp = 0;
       FILE* THIS_IMAGE = fopen(filename.c_str(), "rb");
       for(int i=0;i<height_;i++) {
          for( int k = 0; k < width_; ++k ) {
             fread(&temp, sizeof(type), 1, THIS_IMAGE);
             r_->at(sub2ind(i,k)) = temp;
          }
       }
       if( nc_ == 3 ) {
          for(int i=0;i<height_;i++) {
             for( int k = 0; k < width_; ++k ) {
                fread(&temp, sizeof(type), 1, THIS_IMAGE);
                g_->at(sub2ind(i,k)) = temp;
             }
          }
          for(int i=0;i<height_;i++) {
             for( int k = 0; k < width_; ++k ) {
                fread(&temp, sizeof(type), 1, THIS_IMAGE);
                b_->at(sub2ind(i,k)) = temp;
             }
          }
       }
       fclose(THIS_IMAGE);
    }

  private:
    int width_;
    int height_;
    int nc_;
    vector< type >* r_;
    vector< type >* g_;
    vector< type >* b_;

    //template< class  type >
    inline int sub2ind(int x, int y) {
      return x + y*height_;
    }

    //template< class  type >
    inline void ind2sub(int ind, int& x, int& y) {
    }



    //int sub2ind(int x, int y);
    //void ind2sub(int ind, int& x, int& y);
};

#endif


