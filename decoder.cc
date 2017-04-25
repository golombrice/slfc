#include "decoder.hh"
#include "adaptive_coder.hh"
#include "golomb_coder.hh"
#include "adaptive_binary_coder.hh"

Decoder::Decoder( 
        Image< uint8_t >* Ydecoded,
        vector< Image< uint8_t >* >* cwarpedOK,
        Image< int >* Crestr):
    Coder(Ydecoded, cwarpedOK, Crestr) {
}

Decoder::~Decoder() {
  // nothing?
}

void Decoder::decode(string outputfilename) {

  //Image< uint8_t > Ydecodedvrt(NR_,NC_,3); for (int icomp=0;icomp<3;icomp++){
  //for(int iic = 0; iic < NC_; iic++){
  //    for(int iir = 0; iir < NR_; iir++){
  //      //Ydecodedvrt.set(iir,iic,icomp,Ydecoded_->get(iir,iic,icomp));
  //      Ydecoded_->set(iir,iic,icomp,uint8_t(0));
  //    }
  //  }
  //}
  //for (int icomp=0;icomp<3;icomp++){
  //  int iir = 0;
  //  for(int iic = 0; iic < NC_; iic++){
  //    Ydecoded_->set(iir,iic,icomp,uint8_t(246));
  //  }
  //  iir = NR_-1;
  //  for(int iic = 0; iic < NC_; iic++){
  //    Ydecoded_->set(iir,iic,icomp,uint8_t(246));
  //  }
  //  int iic = 0;
  //  for(int iir = 0; iir < NR_; iir++){
  //    Ydecoded_->set(iir,iic,icomp,uint8_t(246));
  //  }
  //  iic = NC_-1;
  //  for(int iir = 0; iir < NR_; iir++){
  //    Ydecoded_->set(iir,iic,icomp,uint8_t(246));
  //  }
  //}

  AdaptiveCoder coder(outputfilename, MINMAXW_, true);

  vector< int > Sq(9);

  for (int icomp=0;icomp<3;icomp++){
    for(int iir = 0; iir < NR_; iir++){  
      for(int iic = 0; iic < NC_; iic++){

  //for (int icomp=0;icomp<1;icomp++){
  //  for(int iir = 0; iir < 1; iir++){  
  //    for(int iic = 0; iic < 1; iic++){
        bool general_prediction = true;

        for(int i=0;i<9;i++) {
          Sq[i] = i;
        }

        int iR = Crestr_->get(iir,iic,0);
        if(verbose) printf("%s\t%i\t%i\t%i\n","iR, iir, iic", iR, iir, iic);
        int iW = iR;

        general_prediction = fixBoundary(Sq, iir, iic, iR);
        if( iR <= 0 ) {
          general_prediction = true;
        }

        double ypred = 0;

        if(general_prediction){
          ypred = GeneralPredictionp(icomp, iir, iic, Sq);
          if(verbose) printf("%s\n","general prediction!");
        }
        else {
          ypred = WarpedPredictionp(icomp, iir, iic, Sq);
          if(verbose) printf("%s %f\n","warped prediction! y=", ypred);
        }

        int icon = computeContext(icomp, iir, iic);
        //bool general_prediction = true;

        //for(int i=0;i<9;i++) {
        //  Sq[i] = i;
        //}

        //int iR = Crestr_->get(iir,iic,0);
        //if (verbose) printf("%s\t%i\t%i\t%i\n","iR, iir, iic", iR, iir, iic);
        //int iW = iR;

        //general_prediction = fixBoundary(Sq, iir, iic, iR);
        //if( iR <= 0 ) {
        //  general_prediction = true;
        //}

        //double ypred = 0;

        //if(general_prediction){
        //  ypred = GeneralPredictionp(icomp, iir, iic, Sq);
        //  //printf("%s\n","general prediction!");
        //}
        //else {
        //  ypred = WarpedPredictionp(icomp, iir, iic, Sq);
        //  //printf("%s %f\n","warped prediction! y=", ypred);
        //}

        //int icon = computeContext(icomp, iir, iic);

        double symb_to_enc = 0; //ypred - Ydecoded_->get(iir,iic,icomp);
        if( !general_prediction ) {
          int symb_to_enc1 = coder.decode_symbol(iW-1,icon);
          symb_to_enc1 -= 1;
          symb_to_enc = (double)symb_to_enc1 + MINMAXW_->at(icon)->at(iW-1)->at(0);
        }
        else {
          int symb_to_enc1 = coder.decode_symbol(N_regions_,icon);
          symb_to_enc1 -= 1;
          symb_to_enc = (double)symb_to_enc1 + MINMAXW_->at(icon)->at(N_regions_)->at(0);
        }
        ErrorImagec_->set(iir,iic,icomp,symb_to_enc);
        uint8_t temp_Ydecoded = (uint8_t)(ypred - symb_to_enc);
        //cout << iir << " " << iic << " " << icomp <<  "pred: " << ypred << "corr: " << (int) symb_to_enc << "tot: " << (int)temp_Ydecoded << endl;
        Ydecoded_->set(iir,iic,icomp, temp_Ydecoded);
      }
    }
  }
}

void Decoder::readMetadata(string predictor_fn, string mask_fn, string minmax_fn) {
  // INITIALIZE
  int M = 25;

 // for(int j=0;j<3;j++) {
 //   //SPARSEW_.push_back( new vector < vector < int >* >);
 //   for(int i=0;i<N_regions_+1;i++) {
 //     //SPARSEW_->at(j)->push_back( new vector < int > );
 //     for(int k = 0; k < M; ++k ) {
 //       SPARSEW_->at(j)->at(i)->push_back(0);
 //     }
 //   }
 // }
 // for(int j=0;j<3;j++) {
 //   //THETAW_.push_back( new vector < vector < int >* >);
 //   for(int i=0;i<N_regions_+1;i++) {
 //     //THETAW_->at(j)->push_back( new vector < int > );
 //     for(int k = 0; k < M; ++k ) {
 //       THETAW_->at(j)->at(i)->push_back( 0 );
 //     }
 //   }
 // }

 // for(int j=0;j<=N_regions_;j++) {
 //   //MINMAXW_.push_back( new vector < vector < double >* >);
 //   for(int i=0;i<N_CON;i++) {
 //     //MINMAXW_->at(j)->push_back( new vector < double > );
 //     for(int k = 0; k < 2; ++k ) {
 //       if( k == 0 )
 //         MINMAXW_->at(i)->at(j)->push_back( 255 );
 //       if( k == 1 )
 //         MINMAXW_->at(i)->at(j)->push_back( -255 );
 //     }
 //   }
 // }

  // READ

  GolombCoder* gcoder = new GolombCoder(predictor_fn, 1);
  cout << predictor_fn << endl;
  for(int i = 0; i < 3; ++i ) {
    for( int j = 0; j < N_regions_+1; ++j ) {
      vector< int > symbols;
      gcoder->decode_symbols(symbols,5);
      //cout << symbols.size() << endl;
      for( int k = 0; k < symbols.size(); ++k ) {
        //cout << symbols[k] << " " << symbols[k] - (int)( THETAW_.at(i)->at(j)->at(k)*pow(2,10)) << endl;
        //THETAW_.at(i)->at(j)->at(k) = ((double)symbols[k])/((double)pow(2,10));
        THETAW_->at(i)->at(j)->push_back(symbols[k]);
        //cout << "THETA: " << THETAW_.at(i)->at(j)->at(k) << endl;
      }
    }
  }
  delete gcoder;

  cout << mask_fn << endl;
  AdaptiveBinaryCoder* acoder = new AdaptiveBinaryCoder(mask_fn, 1, 25, 1);
  int N_regressors = 0;
  for(int i = 0; i < 3; ++i ) {
    for( int j = 0; j < N_regions_+1; ++j ) {
      if( j==N_regions_ ) {
        N_regressors = 1+9*cwarpedOK_->size();
      }
      else { 
        //N_regressors = 1+4+9*cwarpedOK_->size();
        if( i == 0 ) {
        N_regressors = 1+4+9*cwarpedOK_->size();
        }
        if( i == 1 ) {
        N_regressors = 1+4+9+9*cwarpedOK_->size();
        }
        if( i == 2 ) {
        N_regressors = 1+4+2*9+9*cwarpedOK_->size();
        }
      }
      acoder->initialize_binary(THETAW_->at(i)->at(j)->size(), N_regressors); 
      int n_nnz = 0;
      for( int k = 0; k < N_regressors; ++k ) {
        int symbol_this = acoder->decode_symbol();
        if(symbol_this == 1+1) {
          SPARSEW_->at(i)->at(j)->push_back(k+1);
          n_nnz++;
        }
      }
    }
  }
  delete acoder;

  //gcoder = new GolombCoder(minmax_fn, 1);
  //for(int i = 0; i < N_CON; ++i ) {
  //  for( int j = 0; j < N_regions_+1; ++j ) {
  //    vector< int > symbols;
  //    gcoder->decode_symbols(symbols,2);
  //    for( int k = 0; k < 2; ++k ) {
  //      //cout << symbols.at(k) << endl;
  //      MINMAXW_->at(j)->at(i)->at(k) = (double)symbols.at(k);
  //      //cout << (int) MINMAXW_->at(j)->at(i)->at(k) << endl;
  //    }
  //  }
  //}
  //delete gcoder;

  cout << minmax_fn << endl;
  gcoder = new GolombCoder(minmax_fn, 1);
  vector< int > symbols;
  gcoder->decode_symbols(symbols,15);
  int running_index = 0;
  for(int i = 0; i < N_CON; ++i ) {
    for( int j = 0; j < N_regions_+1; ++j ) {
      for( int k = 0; k < 2; ++k ) {
        //cout << symbols.at(k) << endl;
        //MINMAXW_->at(j)->at(i)->at(k) = (double)symbols.at(k);
        MINMAXW_->at(i)->at(j)->at(k) = ((double)symbols.at(running_index));
        running_index++;
        //cout << (int) MINMAXW_->at(i)->at(j)->at(k) << endl;
      }
    }
  }
  delete gcoder;
}
