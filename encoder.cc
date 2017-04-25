#include <cmath>
#include <algorithm>
#include <thread>
#include "encoder.hh"
#include "adaptive_coder.hh"
#include "adaptive_binary_coder.hh"
#include "golomb_coder.hh"
//#include "lars_main.h"

Encoder::Encoder( 
    Image< uint8_t >* Ydecoded,
    vector< Image< uint8_t >* >* cwarpedOK,
    Image< int >* Crestr):
  Coder(Ydecoded, cwarpedOK, Crestr) {
  }

Encoder::~Encoder() {
  // nothing?
}

void Encoder::encode(string outputfilename) {

  AdaptiveCoder coder(outputfilename, MINMAXW_, false);
  int N_regions = SPARSEW_->at(0)->size()-1; 
  //int no_inW = N_regions;
  int NR = Ydecoded_->size().height;
  int NC = Ydecoded_->size().width;

  vector< int > Sq(9);

  for (int icomp=0;icomp<3;icomp++){
    for(int iir = 0; iir < NR; iir++){  
      for(int iic = 0; iic < NC; iic++){

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

        double symb_to_enc = ypred - Ydecoded_->get(iir,iic,icomp);
        //cout << iir << " " << iic << " " << ypred << " " << symb_to_enc << endl;
        ErrorImagec_->set(iir,iic,icomp,symb_to_enc);
        if( !general_prediction ) {
          if(verbose) cout << MINMAXW_->at(icon)->at(iW-1)->at(0) << endl;
          int symb_to_enc1 = (int) (symb_to_enc - MINMAXW_->at(icon)->at(iW-1)->at(0));
          symb_to_enc1 = symb_to_enc1+1;
          //cout << iir << " " << iic << endl;
          coder.encode_symbol(symb_to_enc1,iW-1,icon);
        }
        else {
          int symb_to_enc1 = (int) (symb_to_enc - MINMAXW_->at(icon)->at(N_regions)->at(0));
          symb_to_enc1 = symb_to_enc1+1;
          //cout << iir << " " << iic << endl;
          coder.encode_symbol(symb_to_enc1,N_regions,icon);
        }

      }
    }
  }
}



void Encoder::computeMinMaxW() {
  int NR = Ydecoded_->size().height;
  int NC = Ydecoded_->size().width;
  int N_regions = SPARSEW_->at(0)->size()-1; 

  vector< int > Sq(9);

  for (int icomp=0;icomp<3;icomp++){
    for(int iir = 0; iir < NR; iir++){  
      for(int iic = 0; iic < NC; iic++){

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
          if(verbose) cout << ypred << "Yd:" << (int)Ydecoded_->get(iir,iic,icomp) << endl;
        }
        else {
          ypred = WarpedPredictionp(icomp, iir, iic, Sq);
          if(verbose) cout << ypred << "Yd:" << (int)Ydecoded_->get(iir,iic,icomp) << endl;
        }

        int icon = computeContext(icomp, iir, iic);

        if(!general_prediction) {
          double symb_to_enc = ypred - Ydecoded_->get(iir,iic,icomp);
          if( symb_to_enc > MINMAXW_->at(icon)->at(iW-1)->at(1) ) {
            MINMAXW_->at(icon)->at(iW-1)->at(1) = symb_to_enc; 
          }
          if( symb_to_enc < MINMAXW_->at(icon)->at(iW-1)->at(0) ) {
            MINMAXW_->at(icon)->at(iW-1)->at(0) = symb_to_enc; 
          }
        } else {
          double symb_to_enc = ypred - Ydecoded_->get(iir,iic,icomp);
          if( symb_to_enc > MINMAXW_->at(icon)->at(N_regions)->at(1) ) {
            MINMAXW_->at(icon)->at(N_regions)->at(1) = symb_to_enc; 
          }
          if( symb_to_enc < MINMAXW_->at(icon)->at(N_regions)->at(0) ) {
            MINMAXW_->at(icon)->at(N_regions)->at(0) = symb_to_enc; 
          }
        }
      }
    }
  }

  for( int i = 0; i <= N_regions; ++i ) {
    for( int j = 0; j < N_CON; ++j ) {
      if( floor(MINMAXW_->at(j)->at(i)->at(0)+.5) == 500 ) {
        MINMAXW_->at(j)->at(i)->at(0) = 0;
      }
      if( floor(MINMAXW_->at(j)->at(i)->at(1)+.5) == -500 ) {
        MINMAXW_->at(j)->at(i)->at(1) = 0;
      }
    }
  }
}

void Encoder::read_predictors(string filename) {
  int Ms = 25;
  int NR = Ydecoded_->size().height;
  int NC = Ydecoded_->size().width;
  int N_regions = SPARSEW_->at(0)->size();

  cout << filename << endl;

  FILE* fid = fopen(filename.c_str(), "rb");

  for (int icomp=0;icomp<3;icomp++){

      for(int iR = 1; iR < N_regions; ++iR) {

        int nnz = 0;
        fread(&nnz, sizeof(int), 1, fid);

        for( int i = 0; i < nnz; ++i) {
          int temp = 0;
          fread(&temp, sizeof(int), 1, fid);
          THETAW_->at(icomp)->at(iR-1)->push_back(temp);
        }

        for( int i = 0; i < nnz; ++i) {
          int temp = 0;
          fread(&temp, sizeof(int), 1, fid);
          SPARSEW_->at(icomp)->at(iR-1)->push_back(temp);
        }

        //fastOLS(A_all.at(iR-1),d_all.at(iR-1),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
      }

      int nnz = 0;
      fread(&nnz, sizeof(int), 1, fid);

      for( int i = 0; i < nnz; ++i) {
        int temp = 0;
        fread(&temp, sizeof(int), 1, fid);
        THETAW_->at(icomp)->at(N_regions-1)->push_back(temp);
      }

      for( int i = 0; i < nnz; ++i) {
        int temp = 0;
        fread(&temp, sizeof(int), 1, fid);
        SPARSEW_->at(icomp)->at(N_regions-1)->push_back(temp);
      }
  }
  fclose(fid);
}

void Encoder::find_predictors() {
  int Ms = 25;
  int NR = Ydecoded_->size().height;
  int NC = Ydecoded_->size().width;
  int N_regions = SPARSEW_->at(0)->size();

  //  vector< vector <  uint8_t > > A_gen;
  //  vector< uint8_t > d_gen;
  // vector< uint8_t > v;
  // v.push_back(3);
  // v.push_back(2);
  // v.push_back(1);
  // A_gen.push_back(v);
  // v.clear();
  // v.push_back(1);
  // v.push_back(21);
  // v.push_back(5);
  // A_gen.push_back(v);
  // v.clear();
  // v.push_back(3);
  // v.push_back(0);
  // v.push_back(12);
  // A_gen.push_back(v);
  // v.clear();
  // v.push_back(5);
  // v.push_back(51);
  // v.push_back(3);
  // A_gen.push_back(v);
  // v.clear();
  // v.push_back(1);
  // v.push_back(1);
  // v.push_back(1);
  // A_gen.push_back(v);

  // d_gen.push_back(2);
  // d_gen.push_back(1);
  // d_gen.push_back(41);
  // d_gen.push_back(6);
  // d_gen.push_back(19);

  // fastOLS(A_gen,d_gen,25,THETAW_->at(0)->at(0), SPARSEW_->at(0)->at(0));
  // return;

  for (int icomp=0;icomp<3;icomp++){
    vector< vector< vector <  uint8_t > > > A_all(N_regions);
    vector< vector< uint8_t > > d_all(N_regions);
    for(int iR = -1; iR < N_regions; ++iR ) {

      for(int iir = 0; iir < NR; iir++){  
        for(int iic = 0; iic < NC; iic++){
          int iR_c = Crestr_->get(iir,iic,0);
          if( iR_c != iR ) {
            continue;
          }

          bool general_prediction = true;
          vector< int > Sq(9);
          for(int i=0;i<9;i++) {
            Sq[i] = i;
          }

          general_prediction = fixBoundary(Sq, iir, iic, iR);
          if( iR <= 0 ) {
            general_prediction = true;
          }

          if( !general_prediction ){
            d_all.at(iR-1).push_back(Ydecoded_->get(iir,iic,icomp));
            vector < uint8_t > v;
            v.push_back(1);
            for( int nk = 0; nk < 4; ++nk ) {
              v.push_back(Ydecoded_->get(iir+dhv[Sq[dcaus[nk]]][0] , iic+dhv[Sq[dcaus[nk]]][1],icomp));
            }
            if( icomp > 0 ) {
              for( int nk = 0; nk < 9; ++nk ) {
                v.push_back(Ydecoded_->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp-1));
              }
            }
            if( icomp > 1 ) {
              for( int nk = 0; nk < 9; ++nk ) {
                v.push_back(Ydecoded_->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp-2));
              }
            }
            for( int ni = 0; ni < cwarpedOK_->size(); ++ni ) {
              for( int nk = 0; nk < 9; ++nk ) {
                v.push_back(cwarpedOK_->at(ni)->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp));
              }
            }
            A_all.at(iR-1).push_back( v );
          }
        }
      }
    }

    bool parallel = false;

    string method = "OOMP";

    if( !parallel ) {
      //cout << "not parallel" << endl;
      if( method == "OOMP") {
        for(int iR = 1; iR < N_regions; ++iR) {
          fastOLS(A_all.at(iR-1),d_all.at(iR-1),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
        }
      }
      if( method == "lars" ) {
        for(int iR = 1; iR < N_regions; ++iR) {
         // main_lars_wrapper(A_all.at(iR-1),d_all.at(iR-1),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
        }
      }
    } 
    else {
      int iR = 1;
      int N_threads = 4;
      vector< int > inds1;
      vector< int > inds2;
      vector< int > inds3;
      vector< int > inds4;
      while(iR <= 1*(N_regions/N_threads)) {
        inds1.push_back(iR);
        ++iR;
      }
      while(iR <= 2*(N_regions/N_threads)) {
        inds2.push_back(iR);
        ++iR;
      }
      while(iR <= 3*(N_regions/N_threads)) {
        inds3.push_back(iR);
        ++iR;
      }
      while(iR <= 4*(N_regions/N_threads)) {
        inds4.push_back(iR);
        ++iR;
      }

      std::thread t1(&Encoder::fastOLS_looper, inds1, std::ref(A_all),std::ref(d_all),Ms,THETAW_->at(icomp), SPARSEW_->at(icomp));
      std::thread t2(&Encoder::fastOLS_looper, inds2, std::ref(A_all),std::ref(d_all),Ms,THETAW_->at(icomp), SPARSEW_->at(icomp));
      std::thread t3(&Encoder::fastOLS_looper, inds3, std::ref(A_all),std::ref(d_all),Ms,THETAW_->at(icomp), SPARSEW_->at(icomp));
      std::thread t4(&Encoder::fastOLS_looper, inds4, std::ref(A_all),std::ref(d_all),Ms,THETAW_->at(icomp), SPARSEW_->at(icomp));
      //  std::thread t2(&Encoder::fastOLS, std::ref(A_all.at(iR-1)),std::ref(d_all.at(iR-1)),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
      //  std::thread t3(&Encoder::fastOLS, std::ref(A_all.at(iR-1)),std::ref(d_all.at(iR-1)),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
      //  std::thread t4(&Encoder::fastOLS, std::ref(A_all.at(iR-1)),std::ref(d_all.at(iR-1)),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
      t1.join();
      t2.join();
      t3.join();
      t4.join();
      for( iR = N_threads*(N_regions/N_threads)+1; iR < N_regions; ++iR) {
        fastOLS(A_all.at(iR-1),d_all.at(iR-1),Ms,THETAW_->at(icomp)->at(iR-1), SPARSEW_->at(icomp)->at(iR-1));
      }
    }

    vector< vector <  uint8_t > > A_gen;
    vector< uint8_t > d_gen;

    for(int iir = 0; iir < NR; iir++){  
      for(int iic = 0; iic < NC; iic++){
        int iR = Crestr_->get(iir,iic,0);

        bool general_prediction = true;
        vector< int > Sq(9);
        for(int i=0;i<9;i++) {
          Sq[i] = i;
        }

        general_prediction = fixBoundary(Sq, iir, iic, iR);
        if( iR <= 0 ) {
          general_prediction = true;
        }

        if(general_prediction){
          try {
          d_gen.push_back(Ydecoded_->get(iir,iic,icomp));
          vector < uint8_t > v;
          v.push_back(1);
          for( int ni = 0; ni < cwarpedOK_->size(); ++ni ) {
            for( int nk = 0; nk < 9; ++nk ) {
              v.push_back(cwarpedOK_->at(ni)->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp));
            }
          }
          A_gen.push_back( v );
          } //catch(...)
           catch (const std::exception& ex) {
             cout << ex.what();
              cout << "ermor: " << iir << " " << iic << endl;
              cout << cwarpedOK_->at(0)->size().width << endl;
              for( int i = 0; i < 9; ++i ) {
                cout << Sq[i] << endl;
              }
              return;
          }
        }
      }
    }
    //cout << THETAW_->at(icomp)->size() << " " << N_regions-1 << endl;
    //cout << "general prediction for " <<  A_gen.size() << endl;
      if( method == "OOMP") {
        fastOLS(A_gen,d_gen,Ms,THETAW_->at(icomp)->at(N_regions-1), SPARSEW_->at(icomp)->at(N_regions-1));
      }
      if( method == "lars" ) {
        //main_lars_wrapper(A_gen,d_gen,Ms,THETAW_->at(icomp)->at(N_regions-1), SPARSEW_->at(icomp)->at(N_regions-1));
      }
  }

}

void Encoder::fastOLS_looper(vector<int> inds, vector< vector < vector< uint8_t > > >& A_all, vector< vector< uint8_t > >& d_all, int Ms, vector< vector< int >* >* thetav, vector< vector< int >* >* sparsev) {
  for(int i = 0; i < inds.size(); ++i) {
    int iR = inds.at(i);
    fastOLS( (A_all.at(iR-1)),(d_all.at(iR-1)),Ms,thetav->at(iR-1), sparsev->at(iR-1));
  }
}

double Encoder::fastOLS(vector< vector< uint8_t > >& A, vector< uint8_t >& d, int Ms, vector< int >* thetav, vector< int >* sparsev) {
  int N = A.size();
  if( N == 0 ) {
    return 0;
  }
  int M = A[0].size();
  if( M == 0 ) {
    return 0;
  }
  if( Ms > M ) {
    Ms = M;
  }
  if( Ms > N ) {
    Ms = N;
  }
  bool all_d_zero = true;
  for( int i = 0; i < d.size(); ++i ) {
    if( d[i] != 0) {
      all_d_zero=false;
      break;
    }
  }
  if( all_d_zero ) {
    return 0;
  }

  //cout << N << " " << M << " " << d.size() << Ms << endl;

  // B = [AA Yd]'*[AA Yd];
  vector< vector< double > > B;
  B.reserve(M+1);
  vector< vector< double > > Bo;
  Bo.reserve(M+1);
  for( int i = 0; i < M; ++i ) {
    vector< double > this_row;
    this_row.reserve(M);
    for( int j = 0; j < M; ++j ) {
      double this_val = 0;
      for( int n = 0; n < N; ++n ) {
        this_val += (double)(A[n][i]*A[n][j]);
      }
      this_row.push_back(this_val);
    }

    double this_val = 0;
    for( int n = 0; n < N; ++n ) {
      this_val += (double)(A[n][i]*d[n]);
    }
    this_row.push_back(this_val);

    B.push_back(this_row);
    Bo.push_back(this_row);
  }

  // last row
  vector< double > this_row;
  this_row.reserve(M);
  for( int j = 0; j < M; ++j ) {
    double this_val = 0;
    for( int n = 0; n < N; ++n ) {
      this_val += (double)(d[n]*A[n][j]);
    }
    this_row.push_back(this_val);
  }
  double this_val = 0;
  for( int n = 0; n < N; ++n ) {
    this_val += (double)(d[n]*d[n]);
  }
  this_row.push_back(this_val);

  B.push_back(this_row);
  Bo.push_back(this_row);

//
//    for(int i = 0; i < A.size(); ++i ) {
//      for(int j = 0; j < A[0].size(); ++j ) {
//        cout << (int)A[i][j] << " ";
//      }
//      cout << endl;
//    }

  // main loop
  vector< vector < double > > C;
  C.reserve(N);
  for( int i = 0; i < N; ++i ){
    vector< double > row(M+1);
    for( int j = 0; j < M; ++j ) {
      if( i == j ) {
        row.at(j)=(1);
      } else {
        row.at(j)=(0);
      }
    }
    C.push_back(row);
  }
    //for(int i = 0; i < C.size(); ++i ) {
    //  for(int j = 0; j < C[0].size(); ++j ) {
    //    cout << C[i][j] << " ";
    //  }
    //  cout << endl;
    //}

  vector< int > Regr;
  for( int i = 0; i < M; ++i ) {
    Regr.push_back(i);
  }
  double crit = B[M][M];

  double Nd = (double)N;
  const double MDL_CONST_TERM = (Nd/2.0)*log2(2.0*3.14159265359*2.71828182846/Nd);

  double MDL_min = 1e20;
  int p_min = 0;


  for( int p = 0; p < Ms; ++p ) {
    //cout << p << endl;
    double sigerr_max = -1e20;
    int max_ind = -1;
    for( int j = p; j < M; ++j ) {
      double C1 = B[j][M]/B[j][j];
      double sigerr = C1*C1*B[j][j]/B[M][M];
      if( sigerr > sigerr_max ) {
        sigerr_max = sigerr;
        max_ind = j;
      }
    }
    if(max_ind==-1){
      break;
    }
    crit -= sigerr_max*B[M][M];
    int temp = Regr[p];
    Regr[p] = Regr[max_ind];
    Regr[max_ind] = temp;

    for( int i = p; i < M+1; ++i ) {
      double tempd = B[i][p];
      B[i][p] = B[i][max_ind];
      B[i][max_ind] = tempd; 
    }
    for( int i = p; i < M+1; ++i ) {
      double tempd = B[p][i];
      B[p][i] = B[max_ind][i];
      B[max_ind][i] = tempd; 
    }

    for( int i = 0; i <= p-1; ++i ) {
      double tempd = C[i][p];
      C[i][p] = C[i][max_ind];
      C[i][max_ind] = tempd; 
    }

    for( int i = p + 1; i < M+1; ++i ) {
      C[p][i]= B[p][i]/B[p][p];
    }
    for(int j = p+1; j < M; ++j) {
      for(int k = j; k < (M+1); ++k ) {
        B[j][k] = B[j][k]-C[p][j]*C[p][k]*B[p][p];
        B[k][j] = B[j][k];
      }
    }

    // backsubstitution
    vector<double> theta(p+1);
    for(int k = p; k >=0; --k) {
      double sp = 0;
      for( int ki = k+1; ki <= p; ++ki ) {
          sp += C[k][ki]*theta[ki];
      }
      theta[k] = C[k][M]-sp;
    }
    // quantization
    double sumthetas = 0;
    //cout << "theta: ";
    for( int k = 0; k <= p; ++k ) {
      theta[k] = floor(theta[k]*1024+.5)/1024;
      //cout << theta[k] << " " << std::floor(fabs(theta[k])) << " ";
      sumthetas+= std::floor(fabs(theta[k]))+1;
    } 
    //cout << endl;

    double rr = Bo[M][M];
    for( int k = 0; k <= p; ++k ) {
      rr -= 2*Bo[Regr[k]][M]*theta[k];
    }
    double tempv[p+1];
    for( int i = 0; i <= p; ++i ) {
      tempv[i] = 0;
      for( int j = 0; j <= p; ++j ) {
        tempv[i] += theta[j]*Bo[Regr[i]][Regr[j]];
      }
    }
    for( int i = 0; i <= p; ++i ) {
      rr += theta[i]*tempv[i];
    }

    if( rr < 1 ) {
      rr = 1;
    }
    double MDL = Nd/2*log2(rr) + MDL_CONST_TERM + (p+1)+sumthetas + 10*(p+1) + MDL_MASK_BIT_COST[p][M-1];
    if(MDL < MDL_min) {
      MDL_min = MDL;
      p_min = p;
    }
    //cout << p << " " << MDL << " " << (rr) << " " << MDL_CONST_TERM << " " << sumthetas << " " << MDL_MASK_BIT_COST[p][M-1] << endl;
    if( rr < 1.000001 || rr > 1e10 ) {
      //cout << "breaking prematurely" << endl;
      break;
    }
    //MDL(p,:) = [N/2*log2(e_sigma(p)) + MDL_const_term, p+sum(floor(abs(theta_this))+1) + N_bits*p + log2(exp(gammaln(M+1)-(gammaln(M-p+1)+gammaln(p+1))))]; % lookup table

  }

  int p = p_min;
  double theta[p+1];
  for(int k = p; k >=0; --k) {
    double sp = 0;
    for( int ki = k+1; ki <= p; ++ki ) {
      sp += C[k][ki]*theta[ki];
    }
    theta[k] = C[k][M]-sp;
  }
  if( thetav != 0 ) {
    for( int k = 0; k <= p; ++k ) {
      thetav->push_back(floor(theta[k]*1024.0+.5));
      sparsev->push_back(Regr[k]+1);
    } 
  }

  //cout << "MIN order: " << p_min+1 << endl;
  //cout << "MDL: " << MDL_min << endl;

  return MDL_min;

}

void Encoder::writeMetadata(string predictor_fn, string mask_fn, string minmax_fn) {

  //for(int i = 0; i < 3; ++i ) {
  //  for( int j = 0; j < no_inW+1; ++j ) {
  //    vector<int> index_vec;
  //    for (int f = 0; f < no_SPARSEW.at(i).at(j); ++f) { 
  //      index_vec.push_back(f); 
  //    }
  //    sort(index_vec.begin(), index_vec.end(),[&](int a, int b) { return SPARSEW[i]->at(j)->at(a) < SPARSEW[i]->at(j)->at(b); });

  //    vector< int > symbols;
  //    for( int k = 0; k < no_SPARSEW.at(i).at(j); ++k ) {
  //      symbols.push_back( (int)(THETAW.at(i)->at(j)->at(index_vec[k])*pow(2,10)) );
  //    }
  //  }
  //}

  int N_regions = SPARSEW_->at(0)->size()-1; 

  GolombCoder* gcoder = new GolombCoder(predictor_fn, 0);
  int N_regressors = 0;
  AdaptiveBinaryCoder* acoder = new AdaptiveBinaryCoder(mask_fn, 1, 25, 0);
  for(int i = 0; i < 3; ++i ) {
    for( int j = 0; j < N_regions+1; ++j ) {
      if( j==N_regions ) {
        N_regressors = 1+9*cwarpedOK_->size();
      }
      else { 
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
      acoder->initialize_binary(SPARSEW_->at(i)->at(j)->size(), N_regressors); 
      vector< int > mask(N_regressors);
      for( int k = 0; k < SPARSEW_->at(i)->at(j)->size(); ++k ) {
        mask.at(SPARSEW_->at(i)->at(j)->at(k)-1) = 1;
      }
      vector< int > symbols;
      for( int k = 0; k < N_regressors; ++k ) {
        int symbol_this = mask.at(k)+1;
        acoder->encode_symbol(symbol_this);
        if(symbol_this == 2){
          int f = 0;
          for( f = 0; f < SPARSEW_->at(i)->at(j)->size(); ++f ) {
            if(SPARSEW_->at(i)->at(j)->at(f)-1 == k) {
              break;
            }
          }
          symbols.push_back( THETAW_->at(i)->at(j)->at(f));
          //cout << j << THETAW_->at(i)->at(j)->at(f) << endl;

        }
      }
      gcoder->encode_symbols(symbols,5);
    }
  }
  delete gcoder;
  delete acoder;

  //gcoder = new GolombCoder(minmax_fn, 0);
  //for(int i = 0; i < N_CON; ++i ) {
  //  for( int j = 0; j < no_inW+1; ++j ) {
  //    vector< int > symbols;
  //    for( int k = 0; k < 2; ++k ) {
  //      symbols.push_back( (int)(MINMAXW[j]->at(i)->at(k)) );
  //    }
  //    gcoder->encode_symbols(symbols,2);
  //  }
  //}
  //delete gcoder;
  gcoder = new GolombCoder(minmax_fn, 0);
  vector< int > symbols;
  for(int i = 0; i < N_CON; ++i ) {
    for( int j = 0; j < N_regions+1; ++j ) {
      for( int k = 0; k < 2; ++k ) {
        //cout << i << " " << j << " " << k << " " << (int)MINMAXW_->(i)->at(j)->at(k) << endl;
        symbols.push_back( (int)floor(MINMAXW_->at(i)->at(j)->at(k)+.5) );
      }
    }
  }
  gcoder->encode_symbols(symbols,15);
  delete gcoder;
}

vector< vector< int > > Encoder::mergeRegions() {

  cout << "heps" << endl;
  // regions-color-matrix
  vector< vector< vector < vector < uint8_t > > > > all_A(N_regions_, vector<vector<vector< uint8_t > > >(3)); // = new vector< vector < vector< uint8_t > > >*;
  // regions-color-vector
  vector< vector< vector < uint8_t > > > all_d(N_regions_, vector<vector< uint8_t > >(3));
  for(int icomp = 0; icomp < 3; ++icomp) {
    for( int iir = 0; iir < NR_; ++iir) {
      for(int iic = 0; iic < NC_; ++iic) {

        int iR = Crestr_->get(iir,iic,0);
        //cout << iir << " " << iR << " " << N_regions_ <<  endl;
        bool general_prediction = true;
        vector< int > Sq(9);
        for(int i=0;i<9;i++) {
          Sq[i] = i;
        }

        general_prediction = fixBoundary(Sq, iir, iic, iR);
        if( iR <= 0 ) {
          general_prediction = true;
        }

        if( general_prediction ){
          continue;
        }

        iR--;

        all_d.at(iR).at(icomp).push_back(Ydecoded_->get(iir,iic,icomp));
        vector < uint8_t > v;
        v.push_back(1);
        for( int nk = 0; nk < 4; ++nk ) {
          v.push_back(Ydecoded_->get(iir+dhv[Sq[dcaus[nk]]][0] , iic+dhv[Sq[dcaus[nk]]][1],icomp));
        }
        if( icomp > 0 ) {
          for( int nk = 0; nk < 9; ++nk ) {
            v.push_back(Ydecoded_->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp-1));
          }
        }
        if( icomp > 1 ) {
          for( int nk = 0; nk < 9; ++nk ) {
            v.push_back(Ydecoded_->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp-2));
          }
        }
        for( int ni = 0; ni < cwarpedOK_->size(); ++ni ) {
          for( int nk = 0; nk < 9; ++nk ) {
            v.push_back(cwarpedOK_->at(ni)->get(iir+dhv[Sq[nk]][0], iic+dhv[Sq[nk]][1],icomp));
          }
        }
        all_A.at(iR).at(icomp).push_back( v );
      }
    }
  }

  vector<double> crit_matrix(N_regions_*N_regions_, -1e5);

  for( int iir = 0; iir < NR_-1; ++iir) {
    for(int iic = 0; iic < NC_-1; ++iic) {
      //cout << iir << " " << iic << endl;
      int iR1 = Crestr_->get(iir,iic,0)-1;
      
      int iRn1 = Crestr_->get(iir,iic+1,0)-1;
      int iRn2 = Crestr_->get(iir+1,iic,0)-1;

      if( iR1 < 0 ) {
        continue;
      }

      //cout << iR1 << " " << iRn1 << " " << iRn2 << endl;

      if( iRn1 > 0 && iRn1 != iR1 && fabs(crit_matrix.at(iR1+iRn1*N_regions_) - (-1e5)) < 1e-6) {
      //cout << "computing criterion for right neighbor" << endl;
        double crit_val_1 = 0;
        double crit_val_2 = 0;
        double crit_val_c = 0;
        for(int icomp = 0; icomp < 3; ++icomp ) {
          double MDL1 = fastOLS(all_A.at(iR1).at(icomp), all_d.at(iR1).at(icomp), 25, 0, 0);
          //cout << all_A.at(iR1).at(icomp).size() << " " << MDL1 << endl;
          double MDL2 = fastOLS(all_A.at(iRn1).at(icomp), all_d.at(iRn1).at(icomp), 25, 0, 0);
          //cout << MDL2 << endl;
          vector< vector< uint8_t > > A_temp = all_A.at(iR1).at(icomp);
          A_temp.insert(A_temp.end(), all_A.at(iRn1).at(icomp).begin(), all_A.at(iRn1).at(icomp).end());
          vector< uint8_t > d_temp = all_d.at(iR1).at(icomp);
          d_temp.insert(d_temp.end(), all_d.at(iRn1).at(icomp).begin(), all_d.at(iRn1).at(icomp).end());
          double MDLc = fastOLS(A_temp, d_temp, 25, 0, 0);
          //cout << MDLc << endl;
          //crit_val+=(MDL1+MDL2)-MDLc;
          crit_val_1+=MDL1;
          crit_val_2+=MDL2;
          crit_val_c+=MDLc;
        }
        crit_matrix.at(iR1+iRn1*N_regions_) = crit_val_1+crit_val_2-crit_val_c;
        crit_matrix.at(iRn1+iR1*N_regions_) = crit_val_1+crit_val_2-crit_val_c;
      }
      if( iRn2 > 0 && iRn2 != iR1 && fabs(crit_matrix.at(iR1+iRn2*N_regions_) - (-1e5)) < 1e-6 ) {
      //cout << "computing criterion for bottom neighbor" << endl;
        //double crit_val = 0;
        double crit_val_1 = 0;
        double crit_val_2 = 0;
        double crit_val_c = 0;
        for(int icomp = 0; icomp < 3; ++icomp ) {
          double MDL1 = fastOLS(all_A.at(iR1).at(icomp), all_d.at(iR1).at(icomp), 25, 0, 0);
          double MDL2 = fastOLS(all_A.at(iRn2).at(icomp), all_d.at(iRn2).at(icomp), 25, 0, 0);
          vector< vector< uint8_t > > A_temp = all_A.at(iR1).at(icomp);
          A_temp.insert(A_temp.end(), all_A.at(iRn2).at(icomp).begin(), all_A.at(iRn2).at(icomp).end());
          vector< uint8_t > d_temp = all_d.at(iR1).at(icomp);
          d_temp.insert(d_temp.end(), all_d.at(iRn2).at(icomp).begin(), all_d.at(iRn2).at(icomp).end());
          double MDLc = fastOLS(A_temp, d_temp, 25, 0, 0);
          //crit_val+=(MDL1+MDL2)-MDLc;
          crit_val_1+=MDL1;
          crit_val_2+=MDL2;
          crit_val_c+=MDLc;
        }
        crit_matrix.at(iR1+iRn2*N_regions_) = crit_val_1+crit_val_2-crit_val_c;
        crit_matrix.at(iRn2+iR1*N_regions_) = crit_val_1+crit_val_2-crit_val_c;
        //crit_matrix.at(iR1+iRn2*N_regions_) = crit_val;
        //crit_matrix.at(iRn2+iR1*N_regions_) = crit_val;
      }
    }
  }

  cout << "starting merging" << endl;

  vector< vector <int> > merged_regions(N_regions_, vector<int>() );
  for(int i = 0; i < N_regions_; ++i ) {
    merged_regions.at(i).push_back(i);
  }
  double max_val = 0;
  do {
    max_val = *std::max_element(crit_matrix.begin(), crit_matrix.end());
    int max_ind = std::distance(crit_matrix.begin(), std::max_element(crit_matrix.begin(), crit_matrix.end()));
    int iR1 = max_ind/N_regions_;
    int iR2 = max_ind-iR1*N_regions_;
    cout << "merging " << iR1 << " and " << iR2 << " with max val " << max_val << endl;
    merged_regions.at(iR1).insert(merged_regions.at(iR1).end(), merged_regions.at(iR2).begin(), merged_regions.at(iR2).end());
    merged_regions.at(iR2).clear();
    vector< int > all_neighbors;
    for(int iRn = 0; iRn < N_regions_; ++iRn) {
      if( iRn == iR1 || iRn == iR2 ) {
        continue;
      }
      if(fabs(crit_matrix.at(iR1+iRn*N_regions_) - (-1e5)) > 1e-6 || fabs(crit_matrix.at(iR2+iRn*N_regions_) - (-1e5)) > 1e-6) {
        all_neighbors.push_back(iRn);
      }
    }
    int iRnew = iR1;
    for(int iR = 0; iR < N_regions_; ++iR) {
      crit_matrix.at(iR+iR1*N_regions_) = -1e5;
      crit_matrix.at(iR1+iR*N_regions_) = -1e5;
      crit_matrix.at(iR+iR2*N_regions_) = -1e5;
      crit_matrix.at(iR2+iR*N_regions_) = -1e5;
    }
    for(int icomp = 0; icomp < 3; ++icomp) {
      all_A.at(iR1).at(icomp).insert(all_A.at(iR1).at(icomp).end(), all_A.at(iR2).at(icomp).begin(), all_A.at(iR2).at(icomp).end());
      all_d.at(iR1).at(icomp).insert(all_d.at(iR1).at(icomp).end(), all_d.at(iR2).at(icomp).begin(), all_d.at(iR2).at(icomp).end());
    }
    for(int icomp = 0; icomp < 3; ++icomp) {
      all_A.at(iR2).at(icomp).clear();
      all_d.at(iR2).at(icomp).clear();
    }
    for(int iRi = 0; iRi < all_neighbors.size(); ++iRi) {
      int iRn1 = all_neighbors.at(iRi);
      double crit_val_1 = 0;
      double crit_val_2 = 0;
      double crit_val_c = 0;
      for(int icomp = 0; icomp < 3; ++icomp ) {
        double MDL1 = fastOLS(all_A.at(iR1).at(icomp), all_d.at(iR1).at(icomp), 25, 0, 0);
        double MDL2 = fastOLS(all_A.at(iRn1).at(icomp), all_d.at(iRn1).at(icomp), 25, 0, 0);
        vector< vector< uint8_t > > A_temp = all_A.at(iR1).at(icomp);
        A_temp.insert(A_temp.end(), all_A.at(iRn1).at(icomp).begin(), all_A.at(iRn1).at(icomp).end());
        vector< uint8_t > d_temp = all_d.at(iR1).at(icomp);
        d_temp.insert(d_temp.end(), all_d.at(iRn1).at(icomp).begin(), all_d.at(iRn1).at(icomp).end());
        double MDLc = fastOLS(A_temp, d_temp, 25, 0, 0);
        //crit_val+=(MDL1+MDL2)-MDLc;
        crit_val_1+=MDL1;
        crit_val_2+=MDL2;
        crit_val_c+=MDLc;
      }
      double crit_val = crit_val_1+crit_val_2-crit_val_c;
      crit_matrix.at(iR1+iRn1*N_regions_) = crit_val;
      crit_matrix.at(iRn1+iR1*N_regions_) = crit_val;
    }

  } while(max_val > 0);


  return merged_regions;

}
