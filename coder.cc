#include "coder.hh"
#include <cmath>

Coder::Coder(Image< uint8_t >* Ydecoded, vector< Image< uint8_t >* >* cwarpedOK, Image< int >* Crestr):
  Ydecoded_(Ydecoded), cwarpedOK_(cwarpedOK), Crestr_(Crestr) 
{
   compute_Nregions();
   //cout << "N_regs: " << N_regions_ << endl;
   SPARSEW_ = new vector< vector< vector < int >* >* >();
   for(int icomp = 0; icomp < 3; ++icomp) {
      vector< vector < int >* >* one_comp = new vector< vector < int >* >;
      for(int n = 0; n < N_regions_+1; ++n ) {
         vector< int >* one_region = new vector< int >;
         one_comp->push_back(one_region);
      }
      SPARSEW_->push_back(one_comp);
   }
   THETAW_ = new vector< vector< vector < int >* >* >();
   for(int icomp = 0; icomp < 3; ++icomp) {
      vector< vector < int >* >* one_comp = new vector< vector < int >* >;
      for(int n = 0; n < N_regions_+1; ++n ) {
         vector< int >* one_region = new vector< int >;
         one_comp->push_back(one_region);
      }
      THETAW_->push_back(one_comp);
   }
   MINMAXW_ = new vector< vector < vector < double >* >* >() ;
   for(int icomp = 0; icomp < N_CON; ++icomp) {
      vector< vector < double >* >* one_comp = new vector< vector < double >* >;
      for(int n = 0; n < N_regions_+1; ++n ) {
         vector< double >* one_region = new vector< double >;
         one_region->push_back(500);
         one_region->push_back(-500);
         one_comp->push_back(one_region);
      }
      MINMAXW_->push_back(one_comp);
   }

   NR_ = Crestr_->size().height;
   NC_ = Crestr_->size().width;
}

void Coder::SetErrorImagep(Image<int>* im) {
  ErrorImagep_ = im;
}

Image<int>* Coder::GetErrorImagep() {
  return ErrorImagep_;
}

void Coder::SetErrorImagec(Image<int>* im) {
  ErrorImagec_ = im;
}

Image<int>* Coder::GetErrorImagec() {
  return ErrorImagec_;
}

Coder::~Coder() {
  // deletions could happen here
   for( int i = 0; i < SPARSEW_->size(); ++i ) {
      for( int j = 0; j < SPARSEW_->at(i)->size(); ++j ) {
         delete SPARSEW_->at(i)->at(j);
      }
   }
   for( int i = 0; i < THETAW_->size(); ++i ) {
      for( int j = 0; j < THETAW_->at(i)->size(); ++j ) {
         delete THETAW_->at(i)->at(j);
      }
   }
   for( int i = 0; i < MINMAXW_->size(); ++i ) {
      for( int j = 0; j < MINMAXW_->at(i)->size(); ++j ) {
         delete MINMAXW_->at(i)->at(j);
      }
   }
}

void Coder::compute_Nregions() {
  N_regions_ = 0;
  int NR = Crestr_->size().height;
  int NC = Crestr_->size().width;

  for(int iir = 0; iir < NR; iir++){  
    for(int iic = 0; iic < NC; iic++){
      int iR = Crestr_->get(iir,iic,0);
      //cout << iR << " ";
      if( iR > N_regions_ ) {
        N_regions_ = iR;
      }
    }
  }
}

int Coder::computeContext(const int icomp, const int ir, const int ic)
{
  double grV = 0;
  double grH = 0;
  double gr45 = 0;
  double gr135 = 0;

  int NR = Crestr_->size().height;
  int NC = Crestr_->size().width;

  if( ir == 0 || ic == 0 || ir == NR-1 ||  ic == NC-1 ) {
    return 0;
  }

//  int this_val = 0;
//  this_val = log2(abs(ErrorImagep_->get(ir,ic,icomp))+1);
//  return this_val;

return 0;

	if (icomp == 0)
	{
		grV = abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp) - *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp) ); //5-4
		grH = ( abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp) - *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) )
			+ abs( *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) - *Ydecoded_->getp(ir+dhv[7][0],ic+dhv[7][1],icomp) ) )/2; //5-2,2-8
		gr45 = abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp) - *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) ); //4-2
	}
	else
	{
		grV =  ( abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp) - *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp) )  //5-4
			+ abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp-1) - *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp-1) ) //5-4
			+ abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp-1) - *Ydecoded_->getp(ir+dhv[5][0],ic+dhv[5][1],icomp-1) ) //4-6
			+ abs( *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp-1)- *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) ) //2-1
			+ abs( *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) - *Ydecoded_->getp(ir+dhv[2][0],ic+dhv[2][1],icomp-1) ) //1-3
			+ abs( *Ydecoded_->getp(ir+dhv[7][0],ic+dhv[7][1],icomp-1) - *Ydecoded_->getp(ir+dhv[6][0],ic+dhv[6][1],icomp-1) ) //8-7
			+ abs( *Ydecoded_->getp(ir+dhv[6][0],ic+dhv[6][1],icomp-1) - *Ydecoded_->getp(ir+dhv[8][0],ic+dhv[8][1],icomp-1) ) )/7; //7-9

		grH = ( abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp) - *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) ) //5-2
			+ abs( *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) - *Ydecoded_->getp(ir+dhv[7][0],ic+dhv[7][1],icomp) )	//2-8
			+ abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp-1) - *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp-1) ) //5-2
			+ abs( *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp-1) - *Ydecoded_->getp(ir+dhv[7][0],ic+dhv[7][1],icomp-1) ) //2-8
			+ abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp-1) - *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) ) //4-1
			+ abs( *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) - *Ydecoded_->getp(ir+dhv[6][0],ic+dhv[6][1],icomp-1) ) //1-7
			+ abs( *Ydecoded_->getp(ir+dhv[5][0],ic+dhv[5][1],icomp-1) - *Ydecoded_->getp(ir+dhv[2][0],ic+dhv[2][1],icomp-1) ) //6-3
			+ abs( *Ydecoded_->getp(ir+dhv[2][0],ic+dhv[2][1],icomp-1) - *Ydecoded_->getp(ir+dhv[8][0],ic+dhv[8][1],icomp-1) ) )/8; //3-9

		gr45 =  ( abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp) - *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp) ) //4-2
			+ abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp-1)- *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp-1) ) //4-2
			+ abs( *Ydecoded_->getp(ir+dhv[5][0],ic+dhv[5][1],icomp-1) - *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) ) //6-1
			+ abs( *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) - *Ydecoded_->getp(ir+dhv[7][0],ic+dhv[7][1],icomp-1) ) //1-8
			+ abs( *Ydecoded_->getp(ir+dhv[2][0],ic+dhv[2][1],icomp-1) - *Ydecoded_->getp(ir+dhv[6][0],ic+dhv[6][1],icomp-1) ) )/5; //3-7

		gr135 =  ( abs( *Ydecoded_->getp(ir+dhv[3][0],ic+dhv[3][1],icomp-1) - *Ydecoded_->getp(ir+dhv[2][0],ic+dhv[2][1],icomp-1) ) //4-3
			+ abs( *Ydecoded_->getp(ir+dhv[4][0],ic+dhv[4][1],icomp-1) - *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) ) //5-1
			+ abs( *Ydecoded_->getp(ir+dhv[0][0],ic+dhv[0][1],icomp-1) - *Ydecoded_->getp(ir+dhv[8][0],ic+dhv[8][1],icomp-1) ) //1-9
			+ abs( *Ydecoded_->getp(ir+dhv[1][0],ic+dhv[1][1],icomp-1) - *Ydecoded_->getp(ir+dhv[6][0],ic+dhv[6][1],icomp-1) ) )/4; //2-7
	}

	int icon = 0;

	if( grV>=2 )
		icon = icon+1;
	icon = icon*2;
	if( grH>=2 )
		icon = icon+1;
	icon = icon*2;
	if( gr45>=2 )
		icon = icon+1;
	icon = icon*2;
	if( gr135>=2 )
		icon = icon+1;

  return icon;
}

double Coder::GeneralPredictionp(const int icomp, const int iir, const int iic, vector< int >& Sq)
{
  int N_regions = SPARSEW_->at(icomp)->size()-1;
  vector< int >* sparsityW = SPARSEW_->at(icomp)->at(N_regions);
  vector< int >* thetaW = THETAW_->at(icomp)->at(N_regions);
  int no_sparse = sparsityW->size();
 // cout << no_sparse << endl;
 // for(int i = 0; i < no_sparse; ++i ) {
 //   cout << sparsityW->at(i) << " " << thetaW->at(i) << endl;
 // }

    double maxv = 0;
    double temp = 0;
    double regressor = 0;
    double ypred = 0;

    const int M_max = 100;
    int ik = 0;

    vector< uint8_t* > XXXY(M_max);

    uint8_t constant = 1;

    XXXY[0] = &constant;

    ik = 1;

    int n_im = 0;
    for(int n_im=0;n_im<cwarpedOK_->size();++n_im) {
        for(int i=ik;i<ik+9;i++) {
            //cout << (int)cwarpedOK_->at(n_im)->get(iir+dhv[Sq[i-ik]][0],iic+dhv[Sq[i-ik]][1], icomp) << endl;
            XXXY[i] = cwarpedOK_->at(n_im)->getp(iir+dhv[Sq[i-ik]][0],iic+dhv[Sq[i-ik]][1], icomp);
            //cout << "ic" << icomp << endl;
            //cout << XXXY[i] << endl;
            //cout << "Sq " << iir+dhv[Sq[i-ik]][0] << endl;
        }
        ik=ik+9;
    }

    int ypredi = 0;
    for(int i=0;i<no_sparse;i++) {
      //cout << thetaW->at(i) << " " << (int)*XXXY[sparsityW->at(i)-1] << " ";
        ypredi += *XXXY[sparsityW->at(i)-1]*thetaW->at(i);
    }
    //cout << endl;
    //ypred = ((double)ypredi);
    ypred = ((double)ypredi)/1024.0;

    ypred = floor(ypred+0.5);

    if(0>ypred)
        ypred = 0;
    if(ypred > 255)
        ypred = 255;

    return ypred;
}

double Coder::WarpedPredictionp(const int icomp, const int iir, const int iic, vector< int >& Sq) 
{
  int iW = Crestr_->get(iir,iic,0);
  vector< int >* sparsityW = SPARSEW_->at(icomp)->at(iW-1);
  vector< int >* thetaW = THETAW_->at(icomp)->at(iW-1);
  int no_sparse = sparsityW->size();

  int ypredi = 0;
  int ik = 0;

  for(int ii = 0; ii < sparsityW->size(); ++ii ) {
    int i = sparsityW->operator[](ii)-1;
    if( i == 0) {
      ypredi += thetaW->operator[](ii);
      continue;
    }
    if( i < 5 ) {
      ypredi += *Ydecoded_->getp(iir+dhv[Sq[dcaus[i-1]]][0] , iic+dhv[Sq[dcaus[i-1]]][1], icomp)*thetaW->operator[](ii);
      continue;
    }
    if( icomp == 0 ) {
      //if( i >= 5 ) {
      int n_im = (i-5)/9;
      ik = n_im*9 + 5;
      int iSq = Sq[i-ik];
      ypredi += *cwarpedOK_->operator[](n_im)->getp(iir+dhv[iSq][0],iic+dhv[iSq][1], icomp)*thetaW->operator[](ii);
      //}
      continue;
    }
    if( icomp > 0 && i < 14 ) {
      ypredi += *Ydecoded_->getp(iir+dhv[Sq[i-5]][0] , iic+dhv[Sq[i-5]][1], icomp-1)*thetaW->operator[](ii);
      continue;
    }
    if( icomp == 1 ) {
      int n_im = (i-14)/9;
      ik = n_im*9 + 14;
      int iSq = Sq[i-ik];
      ypredi += *cwarpedOK_->operator[](n_im)->getp(iir+dhv[iSq][0],iic+dhv[iSq][1], icomp)*thetaW->operator[](ii);
      continue;
    }
    if( icomp > 1 && i < 23 ) {
      ypredi += *Ydecoded_->getp(iir+dhv[Sq[i-14]][0] , iic+dhv[Sq[i-14]][1], icomp-2)*thetaW->operator[](ii);
      continue;
    }
    int n_im = (i-23)/9;
    ik = n_im*9 + 23;
    int iSq = Sq[i-ik];
    ypredi += *cwarpedOK_->operator[](n_im)->getp(iir+dhv[iSq][0],iic+dhv[iSq][1], icomp)*thetaW->operator[](ii);
    continue;
  }
  double ypred= ((double)ypredi)/1024.0;
  ypred = floor(ypred+0.5);

  if(0>ypred)
    ypred = 0;
  if(ypred > 255)
    ypred = 255;

  return ypred;
}
//double Coder::WarpedPredictionp(const int icomp, const int iir, const int iic, vector< int >& Sq) 
//{
//  int iW = Crestr_->get(iir,iic,0);
//  vector< int >* sparsityW = SPARSEW_->at(icomp)->at(iW-1);
//  vector< int >* thetaW = THETAW_->at(icomp)->at(iW-1);
//  int no_sparse = sparsityW->size();
//
//    double maxv = 0;
//    double temp = 0;
//    double regressor = 0;
//    double ypred = 0;
//
//    const int M_max = 100;
//    int ik = 0;
//
//    static vector< uint8_t* > XXXY(M_max);
//
//    uint8_t constant = 1;
//    XXXY[0] = &constant;
//
//    for(int i=1;i<5;i++)
//        XXXY[i] = Ydecoded_->getp(iir+dhv[Sq[dcaus[i-1]]][0] , iic+dhv[Sq[dcaus[i-1]]][1], icomp);
//    ik = 5;
//
//    //cout << cwarpedOK_.size() << endl;
//
//    int n_im = 0;
//    for(int n_im=0;n_im<cwarpedOK_->size();++n_im) {
//        for(int i=ik;i<ik+9;i++) {
//            XXXY[i] = cwarpedOK_->at(n_im)->getp(iir+dhv[Sq[i-ik]][0],iic+dhv[Sq[i-ik]][1], icomp);
//            //cout << "ic" << icomp << endl;
//            //cout << (int)cwarpedOK_[n_im]->get(iir+dhv[Sq[i-ik]][0],iic+dhv[Sq[i-ik]][1], icomp) << endl;
//            //cout << XXXY[i] << endl;
//            //cout << "Sq " << iir+dhv[Sq[i-ik]][0] << endl;
//        }
//        ik=ik+9;
//    }
//
//    //cout << "theta: ";
//    int ypredi = 0;
//    for(int i=0;i<no_sparse;i++) {
//       // cout << xxxy[sparsityw->at(i)-1] << " " << thetaw->at(i) << endl;
//       // cout << sparsityW->at(i) << " " << thetaW->at(i) << " " << (int)*XXXY[sparsityW->at(i)-1] << " ";
//        ypredi += *XXXY[sparsityW->at(i)-1]*thetaW->at(i);
//    }
//    //cout << endl;
//   // ypred= ((double)ypredi);
//    ypred= ((double)ypredi)/1024.0;
//
//    //cout << no_sparse << endl;
//
//    ypred = floor(ypred+0.5);
//
//    if(0>ypred)
//        ypred = 0;
//    if(ypred > 255)
//        ypred = 255;
//
//    return ypred;
//}

bool Coder::fixBoundary(vector< int >& Sq, int iir, int iic, int iR)
{
    bool general_prediction = true;
    int first_one = 1;
    int with_problems = 0;
    int ivok = 999;

   // if( iir == 0 || iic == 0 || iir == Crestr_->size().height-1 || iic == Crestr_->size().width-1 ) {
   //   return true;
   // }

    for(int i=1;i<9;i++){
        int iirc = iir+dhv[i][0];
        int iicc = iic+dhv[i][1];
        if(iirc >= 0 && iicc >=0 && iirc < NR_ && iicc < NC_ && Crestr_->get(iirc,iicc,0) == iR) {
          if(first_one==1 && (i==1 || i==3 || i==4 || i==7)){
            ivok = i;
            first_one = 0;
          }
        }
        else
        {
            with_problems = 1;
        }
    }
    if(ivok==999){
      general_prediction = true;
      for(int i=0;i<9;i++){
        int iirc = iir+dhv[i][0];
        int iicc = iic+dhv[i][1];
        if(iirc >= 0 && iicc >=0 && iirc < NR_ && iicc < NC_ ) {
          Sq[i] = i;
        }
        else
        {
          Sq[i] = 0;
        }
      }
    }
    else
    {
        general_prediction = false;

        for(int i=0;i<9;i++) {
            Sq[i] = i;
        }

        if(with_problems==1){
          for(int i=1;i<9;i++){
            int iirc = iir+dhv[i][0];
            int iicc = iic+dhv[i][1];
            if(iirc >= 0 && iicc >=0 && iirc < NR_ && iicc < NC_ && Crestr_->get(iirc,iicc,0) == iR) {
            //if(Crestr_->get(iir+dhv[i][0],iic+dhv[i][1],0) == iR){
              Sq[i] = i;
            }
            else
            {
              Sq[i] = ivok;
            }
          }
        }
    }

    return general_prediction;
}

