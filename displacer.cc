#include "displacer.hh"
#include "bit_output.hh"
#include "bit_input.hh"
#include "gen_types.hh"

Displacer::Displacer(Image<uint8_t>* im1, Image<uint8_t>* im2, Image<int>* BqF): 
  im1_(im1), im2_(im2), BqF_(BqF), BqFSV_(0), best_disp_h_(0), best_disp_v_(0) {
    best_disp_h_ = new vector< int >();
    best_disp_v_ = new vector< int >();
}

Displacer::Displacer(Image<int>* BqF, string displacements_file): 
  im1_(0), im2_(0), BqF_(BqF), BqFSV_(0), best_disp_h_(0), best_disp_v_(0) {
    best_disp_h_ = new vector< int >();
    best_disp_v_ = new vector< int >();
    read_displacements(displacements_file);
}

Displacer::~Displacer() {
  delete best_disp_h_;
  delete best_disp_v_;
}

Image< int >* Displacer::get_BqFSV() {
  return BqFSV_;
}

vector<int>* Displacer::get_displacements_h() {
  return best_disp_h_;
}

void Displacer::reorder_displacements(vector<int>& indices) {
  best_disp_h_new_ = new vector<int>();
  best_disp_v_new_ = new vector<int>();
  for( int i = 0; i < indices.size(); ++i ) {
    best_disp_h_new_->push_back(best_disp_h_->at(indices(i)));
    best_disp_v_new_->push_back(best_disp_h_->at(indices(i)));
  }
  delete best_disp_h_;
  delete best_disp_v_;
  best_disp_h_ = best_disp_h_new_;
  best_disp_v_ = best_disp_v_new_;
}

void Displacer::write_displacements(string displacements_file) {
  //FILE* fid = fopen(displacements_file.c_str(), "wb");
  //for(int i = 0; i < best_disp_h_->size(); ++i ) {
  //  fwrite(&best_disp_h_->at(i), sizeof(int), 1, fid);
  //  fwrite(&best_disp_v_->at(i), sizeof(int), 1, fid);
  //}
  BitOutput bo(displacements_file);
  for( int bit = 0; bit < 10; ++bit) {
    bo.output_bit(best_disp_h_->size() & (1 << bit) );
  }

  for(int i = 0; i < best_disp_h_->size(); ++i ) {
    for( int bit = 0; bit < 5; ++bit) {
      bo.output_bit((best_disp_h_->at(i)+8) & (1 << bit) );
    }
    for( int bit = 0; bit < 5; ++bit) {
      bo.output_bit((best_disp_v_->at(i)+8) & (1 << bit) );
    }
  }
  //fclose(fid);

}

void Displacer::read_displacements(string displacements_file) {
  BitInput bo(displacements_file);
  int N = 0;
  for( int bit = 0; bit < 10; ++bit) {
    N += (bo.input_bit() << bit);
  }
  for( int i = 0; i < N; ++i ) {
    best_disp_h_->push_back(0);
    best_disp_v_->push_back(0);
    for( int bit = 0; bit < 5; ++bit) {
      best_disp_h_->at(i) += (bo.input_bit() << bit);
    }
    best_disp_h_->at(i) -= 8;
    for( int bit = 0; bit < 5; ++bit) {
      best_disp_v_->at(i) += (bo.input_bit() << bit);
    }
    best_disp_v_->at(i) -= 8;
  }
  
 // FILE* fid = fopen(displacements_file.c_str(), "rb");
 // int temp = 0;
 // while(fread(&temp, sizeof(int), 1, fid)==1) {
 //   best_disp_h_->push_back(temp);
 //   fread(&temp, sizeof(int), 1, fid);
 //   best_disp_v_->push_back(temp);
 // }
 // fclose(fid);
}

void Displacer::search_displacements() {

  int NR = im1_->size().height;
  int NC = im1_->size().width;

  // Find number of regions
  int N_regions = 0;
  for(int iir = 0; iir < NR; iir++){  
    for(int iic = 0; iic < NC; iic++){
      int iR = BqF_->get(iir,iic,0);
      if( iR > N_regions ) {
        N_regions = iR;
      }
    }
  }

  // Displacement grid
  vector< int > displace_opts;
  int min_displacement = -8;
  for( int i = min_displacement; i <= 8; ++i ) {
    displace_opts.push_back(i);
  }
  int N_opts = displace_opts.size();
  vector< vector< vector < double > > > MSE_displacements;

  // Initialize for each region an array with displacement grid
  for(int i = 0; i < N_regions; ++i ) {
    vector< vector <double> > MSE_vh;
    for( int j = 0; j < displace_opts.size(); ++j ) {
      vector< double > MSE_one_region(displace_opts.size());
      MSE_vh.push_back(MSE_one_region);
    }
    MSE_displacements.push_back(MSE_vh);
  }

  // Compute MSE
  vector< double > pix_count(N_regions);
  for(int iir = 0; iir < NR; iir++){  
    for(int iic = 0; iic < NC; iic++){
      int iR = BqF_->get(iir,iic,0);
      if( iR < 1 ) {
        continue;
      }
      bool pixok = iir > 8 && iic > 8 && iir < NR-9 && iic < NC-9;
      for (int icomp=0;icomp<3;icomp++){
        double pix1 = (double)im1_->get(iir,iic,icomp);
        //bool pixok = true;
        //for(int i = 0; i < N_opts; ++i ) {
        //  int disp_h = displace_opts.at(i);
        //  for(int j = 0; j < N_opts; ++j ) {
        //    int disp_v = displace_opts[j];
        //    int ih = iic+disp_h; 
        //    int iv = iir+disp_v;
        //    if( ih < 1 || iv < 1 || ih > NC-1 || iv > NR-1 ) { 
        //      pixok = false;
        //    }
        //  }
        //}
        if(pixok) {
          pix_count.at(iR-1) += 1;
          for(int i = 0; i < N_opts; ++i ) {
            int disp_h = displace_opts[i];
            for(int j = 0; j < N_opts; ++j ) {
              int disp_v = displace_opts[j];
              int ih = iic+disp_h; 
              int iv = iir+disp_v;
              double pix2 = (double)im2_->get(iv,ih,icomp);
              //double err = (pix2-pix1)*(pix2-pix1);
              //if(verbose) cout << iir << " " << iic << " " << iv << " " << ih << " " << pix1 << " " << pix2 << " " << err << endl;
              MSE_displacements[iR-1][i][j] += (pix2-pix1)*(pix2-pix1);
              //MSE_displacements_v.at(iR-1).at(j) += err;
            }
          }
        }
      }
    }
  }

  //vector< int > best_disp_h_(N_regions);
  //vector< int > best_disp_v_(N_regions);


  //FILE* outfile = fopen(filename.c_str(), "wb");
  for( int ir = 0; ir < N_regions; ++ir ) {
    double min_h = 1e20;
    int min_ih = 0; int min_iv = 0;
    for( int i = 0; i < displace_opts.size(); ++i ) {
      for( int j = 0; j < displace_opts.size(); ++j ) {
      if(verbose) cout << MSE_displacements.at(ir).at(i).at(j) << endl;
        if( MSE_displacements.at(ir).at(i).at(j)/pix_count.at(ir) <= min_h ) {
          min_ih = i; 
          min_iv = j;
          min_h = MSE_displacements.at(ir).at(i).at(j)/pix_count.at(ir);
        }
      }
    }
    best_disp_h_->push_back(displace_opts.at(min_ih));
    best_disp_v_->push_back(displace_opts.at(min_iv));
    int region = ir+1;
    //fwrite(&region,1,sizeof(int),outfile);
   // fwrite(&best_disp_v_.at(ir),1,sizeof(int),outfile);
   // fwrite(&best_disp_h_.at(ir),1,sizeof(int),outfile);
  }
  //fclose(outfile);

}

void Displacer::propagateRegions() {

  int NR = BqF_->size().height;
  int NC = BqF_->size().width;

  int N_regions = 0;
  for(int iir = 0; iir < NR; iir++){  
    for(int iic = 0; iic < NC; iic++){
      int iR = BqF_->get(iir,iic,0);
      if( iR > N_regions ) {
        N_regions = iR;
      }
    }
  }

  // initialize
  BqFSV_ = new Image< int >(NR,NC,1);
  for( int iir = 0; iir < NR; ++iir ) {
    for( int iic = 0; iic < NC; ++iic ) {
        BqFSV_->set(iir, iic, 0, -1);
    } 
  }
  // displace
  //for(int iR = N_regions-1; iR >= 0; --iR ) {
  //  int disp_h = displacements_h->at(iR);
  //  int disp_v = displacements_v->at(iR);
  //  for( int iir = 1; iir < NR - 1; ++iir ) {
  //    for( int iic = 1; iic < NC - 1; ++iic ) {
  //      if( iR == BqF_->get(iir,iic,0) ) {
  //        int iirp = iir+disp_v; int iicp = iic+disp_h;
  //        if( iirp >= 1 && iirp <= NR-1 && iicp >= 1 && iicp <= NC-1) {
  //          BqFSV_->set(iirp, iicp, 0, iR);
  //          N_unknown_pixels--;
  //        }
  //      }
  //    } 
  //  }
  //}

  //cout << "displace" << endl;
  for( int iir = 0; iir < NR; ++iir ) {
    for( int iic = 0; iic < NC; ++iic ) {
      int iR = BqF_->get(iir,iic,0);
      if( iR < 1 ) {
        continue;
      }
      int disp_h = best_disp_h_->at(iR-1);
      int disp_v = best_disp_v_->at(iR-1);
      int iirp = iir+disp_v; int iicp = iic+disp_h;
      if( iirp >= 0 && iirp <= NR-1 && iicp >= 0 && iicp <= NC-1) {
        int this_region = BqFSV_->get(iirp,iicp,0);
        if( this_region == -1 || this_region > iR ) {
          BqFSV_->set(iirp, iicp, 0, iR);
        }
      }
    }
  } 
  int N_unknown_pixels = (NR)*(NC);
  for( int iir = 0; iir < NR; ++iir ) {
    for( int iic = 0; iic < NC; ++iic ) {
      if( BqFSV_->get(iir,iic,0) != -1 ) {
        N_unknown_pixels--;
      }
    }
  }
  //cout << "median fill" << endl;
  // median fill
  for( int iir = 0; iir < NR; ++iir ) {
    for( int iic = 0; iic < NC; ++iic ) {
      if( BqFSV_->get(iir,iic,0) == -1 ) {
        int N_ok = 0;
        int max_depth = -1;
        int N_known_depth = 0;
        for( int i = 0; i < 9; ++i ) {
          int inr = iir+dhv[i][0];
          int inc = iic+dhv[i][1];
          if( inr >= 0 && inr <= NR-1 && inc >= 0 && inc <= NC-1) {
            N_ok++;
            int this_depth = BqFSV_->get(inr,inc,0);
            if( this_depth > max_depth ) {
              max_depth = this_depth;
            }
            if( this_depth != -1 ) {
              N_known_depth++;
            } 
          }
        }
        if( N_ok > 3 && N_known_depth > 3) {
          BqFSV_->set(iir, iic, 0, max_depth);
          N_unknown_pixels--;
        }
      }
    } 
  }
  //cout << "fill holes" << endl;
  // fill holes 
  //for(int iR = N_regions-1; iR >= 0; --iR ) {
  //  for( int iir = 0; iir < NR; ++iir ) {
  //    for( int iic = 0; iic < NC; ++iic ) {
  //      if( N_unknown_pixels == 0) {
  //        break;
  //      }
  //      if( iR == BqF_->get(iir,iic,0) ) {
  //        for( int iirp = iir-8; iirp <= iir+8; ++iirp ) {
  //          for( int iicp = iic-8; iicp <= iic+8; ++iicp ) {
  //            if( iirp >= 0 && iirp <= NR-1 && iicp >= 0 && iicp <= NC-1) {
  //              if( BqFSV_->get(iirp,iicp,0) == -1) {
  //                BqFSV_->set(iirp, iicp, 0, iR);
  //                N_unknown_pixels--;
  //              }
  //            }
  //          }
  //        }
  //      }
  //    } 
  //  }
  //}
  for( int iir = 0; iir < NR; ++iir ) {
    for( int iic = 0; iic < NC; ++iic ) {
      if( N_unknown_pixels == 0) {
        break;
      }
      if( BqFSV_->get(iir,iic,0) == -1) {
        int iR = -1;
        for( int iirp = iir-8; iirp <= iir+8; ++iirp ) {
          for( int iicp = iic-8; iicp <= iic+8; ++iicp ) {
            if( iirp >= 0 && iirp <= NR-1 && iicp >= 0 && iicp <= NC-1) {
              int iR1 = BqF_->get(iirp,iicp,0);
              if( iR1 > iR) {
                iR = iR1;
              }

            }
          }
        }
        BqFSV_->set(iir, iic, 0, iR);
        N_unknown_pixels--;
      }
    }
  }
}
