#include <cstdlib> 
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdint.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
//#include <chrono>
#include "image.hh"
#include "golomb_coder.hh"
#include "gen_types.hh"
#include "decoder.hh"
#include "encoder.hh"
#include "displacer.hh"
#include "input_parser.hh"
#include "avi_reader.hh"
#include "openjpeg.h"
//#include <openjpeg-2.2/openjpeg.h>
#include "cerv/cerv.h"
#include "segm/msImageProcessor.h"
//#define MINIZ_HEADER_FILE_ONLY
#include "miniz.c"

using std::vector;
using std::string;
using std::cout;
using std::endl;

typedef std::pair<double,int> mypair;
struct compare_pairs {
  inline bool operator()( const mypair& l, const mypair& r)
  { return l.first < r.first; }
};

void compress_with_jp2(Image<uint8_t>* A_ij, string filename, int precision=8);
void compress_with_jp2_BqF(Image<int>* A_ij, string filename, int precision=16);
Image<uint8_t>* decompress_with_jp2(string filename);
Image<int>* decompress_with_jp2_BqF(string filename);
void cerv_compress(Image< int >* BqF, string fn, string fn_labels);
void cerv_decompress(Image< int >* BqF, string fn, string fn_labels);
Image<int>* segment(Image< uint8_t >* A, const int h_s, const int h_r);
Image<int>* order_segmentation_by_depth(Image< int >* BqF, vector< Image< uint8_t >* >& views);
void merge_center(Image< int >* BqF, vector< vector< int > >& merged_regions);
bool create_zip(vector<string>& files_to_archive, string zip_filename);
bool zip_extract(string zip_filename);

string append_string(string base, int n) {
  ostringstream ss;
  ss << base << n;
  return ss.str();
}

int main(int argc, char** argv) 
{
  int NR = 434;
  int NC = 625;


  InputParser parser(argc, argv);
  bool decode = parser.cmdOptionExists("-d");

  string dir_name = "";
  if(decode) {
    string input_file = parser.getCmdOption("-i");
    dir_name = "tmp_" + input_file;
  }
  else {
    string input_file = parser.getCmdOption("-i");
    ostringstream ss;
    ss << input_file << ".slfc";
    string compressed_file = ss.str();
    if( parser.cmdOptionExists("-o") ) {
      compressed_file = parser.getCmdOption("-o");
    }
    dir_name = "tmp_" + compressed_file;
  }

  string displacement_fn = dir_name + "/displacements_";
  string pred_fn = dir_name + "/pred_";
  string residuals_fn = dir_name + "/residuals_";
  string minmax_fn = dir_name + "/minmax_";
  string mask_fn = dir_name + "/mask_";
  string center_fn = dir_name + "/center.jp2";
  string partition_fn = dir_name + "/partition";
  string partition_labels_fn = dir_name + "/partition_labels";


  system(string("mkdir \"" + dir_name + "\"").c_str());
  
  if(decode) {
    string input_file = parser.getCmdOption("-i");
    string decoded_fn = input_file + "_decoded.yuv";
    zip_extract(input_file);
    vector< Image< uint8_t >* > views;
    Image<uint8_t>* A_icjc = decompress_with_jp2(center_fn);
    
    //Image<uint8_t>* A_icjc = new Image<uint8_t>(NR,NC,3);
    //A_icjc->readImage("center.txt");
    views.push_back(A_icjc);
    A_icjc->writeImage(decoded_fn);
    cout << (int)A_icjc->get(0,0,0) << endl;
    Image< int >* BqF = new Image< int >(NR,NC,1);
    cerv_decompress(BqF, partition_fn, partition_labels_fn);

    Image< int >* ErrorImagep = new Image<int>(NR,NC,3);
    Image< int >* ErrorImagec = new Image<int>(NR,NC,3);
    for(int n = 1; n < 225; ++n ) {
      cout << n << endl;
      string displacement_fn_n = append_string(displacement_fn, n);

      // 2. search displacements
      Displacer* displacer = new Displacer(BqF, displacement_fn_n);
      displacer->propagateRegions();
      Image<int>* BqFSV = displacer->get_BqFSV();
      //BqFSV->writeImage("BqFSV_dec.txt");

      // Read regressors
      vector< Image< uint8_t >* >* A_regressors = new vector< Image< uint8_t >* >();
      for(int i = 0; i < 5; ++i)  {
        int iii = A_REGR_INDS[n-1][i];
        if( iii == 0 ) {
          continue;
        }
        Image< uint8_t >* A_r = views.at(iii-1);
        A_regressors->push_back(A_r);
      }
      Image< uint8_t >* A_ij = new Image< uint8_t >(NR,NC,3);
      Decoder decoder(A_ij, A_regressors, BqFSV);
      decoder.SetErrorImagep(ErrorImagep);
      decoder.SetErrorImagec(ErrorImagec);

      string pred_fn_n = append_string(pred_fn, n);
      string mask_fn_n = append_string(mask_fn, n);
      string minmax_fn_n = append_string(minmax_fn, n);
      string residuals_fn_n = append_string(residuals_fn, n);
      string decoded_fn_n = append_string(decoded_fn, n);

      decoder.readMetadata(pred_fn_n, mask_fn_n, minmax_fn_n); 
      decoder.decode(residuals_fn_n);
      delete ErrorImagep;
      ErrorImagep = decoder.GetErrorImagec();
      ErrorImagec = new Image<int>(NR,NC,3);
      A_ij->writeImage(decoded_fn, "ab");
      views.push_back(A_ij);
    }
    system(string("rm -rf \"" + dir_name + "\"").c_str());
  }
  
  else {
    // ENCODE
    string input_file = parser.getCmdOption("-i");
    ostringstream ss;
    ss << input_file << ".slfc";
    string compressed_file = ss.str();
    if( parser.cmdOptionExists("-o") ) {
      compressed_file = parser.getCmdOption("-o");
    }
    vector<string> files_to_be_archived;
    files_to_be_archived.push_back(center_fn);
    files_to_be_archived.push_back(partition_fn);
    files_to_be_archived.push_back(partition_labels_fn);

    AVIReader avi_reader(input_file, NR, NC);
    vector< Image< uint8_t >* > views;
    cout << "Reading views from disk to memory." << endl; 
    for(int n = 0; n < 225; ++n) {
      views.push_back(avi_reader.get_view(n));
    }

    // get center view and encode it with JPEG 2000.
    Image< uint8_t >* A_icjc = views.at(0);
    compress_with_jp2(A_icjc, string(center_fn));

    Image<int>* BqF = 0;
    //Image<int>* BqFs = 0;
    if( parser.cmdOptionExists("-p") ) {
      string BqF_filename = parser.getCmdOption("-p");
      BqF = new Image< int >(NR,NC,1);
      BqF->readImage(BqF_filename);
    }
    else {
      // segment center view with mean shift
      int h_s = 8; 
      int h_r = 8;
      if( parser.cmdOptionExists("-mhs") ) {
        h_s = atoi( parser.getCmdOption("-mhs").c_str() );
      }
      if( parser.cmdOptionExists("-mhr") ) {
        h_r = atoi( parser.getCmdOption("-mhr").c_str() );
      }
      cout << "Segmenting center view with mean shift parameters " << h_s << " " << h_r << endl;
      BqF = segment(A_icjc, h_s, h_r);

    }

    vector< Displacer* > displacers;
    displacers.reserve(225);
    for(int n = 1; n < 225; ++n ) {
      cout << "Searching displacements for view " << n << endl;
      // 1. read image A_ij
      Image< uint8_t >* A_ij= views.at(n);
      // 2. search displacements
      Displacer* displacer = new Displacer(A_icjc, A_ij, BqF);
      displacers.at(n) = displacer;
      displacer->search_displacements();
    }

    if( !parser.cmdOptionExists("-p") ) {
      cout << "Reordering the labels of the segmentation according to depth" << endl;
      
      Image< int >* BqFn = order_segmentation_by_depth(BqF, views, displacers);
      delete BqF;
      BqF = BqFn;
    }

    bool optimize_regions = parser.cmdOptionExists("-opt");
    if( optimize_regions ) {
      vector< Image< uint8_t >* >* A_regressors = new vector< Image< uint8_t >* >();
      int n = 7;
      Image< uint8_t >* A_ij= views.at(n);
      for(int i = 0; i < 5; ++i)  {
        int iii = A_REGR_INDS[n-1][i];
        if( iii == 0 ) {
          continue;
        }
        Image< uint8_t >* A_r= avi_reader.get_view(iii-1);
        A_regressors->push_back(A_r);
      }
      Displacer* displacer = new Displacer(A_icjc, A_ij, BqF);
      displacer->search_displacements();
      displacer->propagateRegions();
      Image<int>* BqFSV = displacer->get_BqFSV();
      Encoder encoder(A_ij, A_regressors, BqFSV);
      cout << "Optimizing regions." << endl;
      vector< vector< int > > merged_regions = encoder.mergeRegions();
      merge_center(BqF, merged_regions);
    }

    if( parser.cmdOptionExists("-s") ) {
      BqF->writeImage(input_file + "_BqFSV.txt");
    }

    // Read partition and compress it with cerv.
    cout << "Compress center partition with CERV" << endl;
    cerv_compress(BqF, partition_fn, partition_labels_fn);

    Image< int >* ErrorImagep = new Image<int>(NR,NC,3);
    Image< int >* ErrorImagec = new Image<int>(NR,NC,3);

    for(int n = 1; n < 225; ++n ) {
      cout << "Encoding view " << n << endl;
      // 1. read image A_ij
      Image< uint8_t >* A_ij= views.at(n);
      string displacement_fn_n = append_string(displacement_fn, n);
      displacer = displacers.at(n);
      displacer->write_displacements(displacement_fn_n);
      // 3. get BqFSV
      displacer->propagateRegions();
      Image<int>* BqFSV = displacer->get_BqFSV();
      if( parser.cmdOptionExists("-s") ) {
        BqFSV->writeImage(input_file + "_BqFSV.txt", "ab");
        delete displacer; 
        continue;
      }
      // 4. read cwarpedOK
      vector< Image< uint8_t >* >* A_regressors = new vector< Image< uint8_t >* >();
      for(int i = 0; i < 5; ++i)  {
        int iii = A_REGR_INDS[n-1][i];
        if( iii == 0 ) {
          continue;
        }
        Image< uint8_t >* A_r = views.at(iii-1);
        A_regressors->push_back(A_r);
      }
      // 5. search predictors
      Encoder encoder(A_ij, A_regressors, BqFSV);
      encoder.SetErrorImagep(ErrorImagep);
      encoder.SetErrorImagec(ErrorImagec);
      //cout << "find predictors" << endl;
      if( parser.cmdOptionExists("-c")) {
        string predfile = parser.getCmdOption("-c");
        ostringstream ssp;
        ssp << predfile << n+1;
        predfile = ssp.str();
        encoder.read_predictors(predfile);
      }
      else {
        encoder.find_predictors();
      }
      //cout << "compute minmax" << endl;
      encoder.computeMinMaxW();
      // 6. write residuals and metadata
      string pred_fn_n = append_string(pred_fn, n);
      string mask_fn_n = append_string(mask_fn, n);
      string minmax_fn_n = append_string(minmax_fn, n);
      string residuals_fn_n = append_string(residuals_fn, n);
      //cout << "write residuals" << endl;
      encoder.encode(residuals_fn_n);
      delete ErrorImagep;
      ErrorImagep = encoder.GetErrorImagec();
      ErrorImagec = new Image<int>(NR,NC,3);
      //cout << "write header" << endl;
      encoder.writeMetadata(pred_fn_n, mask_fn_n, minmax_fn_n); 
      files_to_be_archived.push_back(displacement_fn_n);
      files_to_be_archived.push_back(pred_fn_n);
      files_to_be_archived.push_back(minmax_fn_n);
      files_to_be_archived.push_back(mask_fn_n);
      files_to_be_archived.push_back(residuals_fn_n);
      delete displacer;
    }
    if( !parser.cmdOptionExists("-s") ) {
      create_zip(files_to_be_archived, compressed_file);
    }
    system(string("rm -rf \"" + dir_name + "\"").c_str());
  }
 return EXIT_SUCCESS; 
}

void merge_center(Image< int >* BqF, vector< vector< int > >& merged_regions) {
  for(int iir = 0; iir < BqF->size().height; ++iir) {
    for(int iic = 0; iic < BqF->size().width; ++iic) {
      int iR = BqF->get(iir,iic,0)-1;
      cout << iir << " " << iic << " " << iR << endl;
      if( iR < 0 || merged_regions.at(iR).size() > 0 ) {
        continue;
      }

      for( int i = 0; i < merged_regions.size(); ++i) {
       if( std::find(merged_regions.at(i).begin(), merged_regions.at(i).end(), iR) != merged_regions.at(i).end()) {
         BqF->set(iir,iic,0, i+1);
         break;
       }
      }
    }
  }
}

Image< int >* order_segmentation_by_depth(Image< int >* BqF, vector< Image< uint8_t >* >& views, vector< Displacer* >& displacers ) {
    // find displacements for a certain column of views

    //cout << "searching displacements" << endl;
    int NR = BqF->size().height;
    int NC = BqF->size().width;

    Image< uint8_t >* A_icjc = views.at(0);
    vector< Displacer* > displacers;
    vector< vector< double >* > mean_disp_for_all;
    vector< int > x_opts;
    for(int jv = -5; jv <= 5; ++jv) {
      x_opts.push_back(jv);
    }
    for(int jvi = 0; jvi < x_opts.size(); ++jvi ) {
      int jv = x_opts.at(jvi)+8;

      vector<double>* mean_displacement = new vector< double >;
      for(int iv = 3; iv <= 13; ++iv ) {
        //cout << iv << " " << jv << endl;
        int n = 0;
        if( iv == 8 && jv == 8 ) {
          n = 0;
        } else {
          n = get_n(iv,jv);
        }
        Image< uint8_t >* A_ij= views.at(n);
        //Displacer* displacer = new Displacer(A_icjc, A_ij, BqF);
        //displacer->search_displacements();
        //displacers.push_back(displacer);
        Displacer* displacer = displacers.at(n);
        vector<int>* disp_h = displacer->get_displacements_h();
        if(mean_displacement->size() == 0){
          for( int ii = 0; ii < disp_h->size(); ++ii ) {
            mean_displacement->push_back(disp_h->at(ii));
          }
        }
        else {
          for( int ii = 0; ii < disp_h->size(); ++ii ) {
            mean_displacement->at(ii) += disp_h->at(ii);
          }
        }
      }
      for( int ii = 0; ii < mean_displacement->size(); ++ii ) {
        mean_displacement->at(ii) /= 11.0;
      }
      mean_disp_for_all.push_back(mean_displacement);
    }
    // fit line
    vector< double > ks(mean_disp_for_all.at(0)->size());
    for(int j = 0; j < mean_disp_for_all.size(); ++j) {
      for(int i = 0; i < mean_disp_for_all.at(j)->size(); ++i) {
        ks.at(i) += x_opts.at(j)*mean_disp_for_all.at(j)->at(i);
      }
    }
    double k_c = 0;
    for(int j = 0; j < x_opts.size(); ++j ) {
      k_c += x_opts[j]*x_opts[j];
    }
    vector< mypair > ks_p;
    for(int i = 0; i < ks.size(); ++i ) {
      ks_p.push_back( mypair(1/k_c*ks.at(i)*x_opts.at(x_opts.size()-1), i+1) );
    }
    std::sort(ks_p.begin(), ks_p.end(), compare_pairs());
   // for(int i = 0; i < ks_p.size(); ++i ) {
   //   cout << ks_p[i].first << " " << ks_p[i].second << endl;
   // }
   // for(int i = 0; i < mean_disp_for_all.at(0)->size(); ++i ) {
   // }

      // sort the displacements
    //  vector< mypair > mean_displacement_p;
    //  for(int i = 0; i < mean_displacement.size(); ++i ) {
    //    mean_displacement_p.push_back( mypair(mean_displacement.at(i), i) );
    //  }
    //  std::sort(mean_displacement_p.begin(), mean_displacement_p.end(), compare_pairs());
    //  for(int i = 0; i < mean_displacement_p.size(); ++i ) {
    //    cout << mean_displacement_p[i].first << " " << mean_displacement_p[i].second << endl;
    //  }

    Image< int >* BqFn = new Image< int >(NR,NC,1);
    for(int ir = 0; ir < NR; ++ir) {
      for(int ic = 0; ic < NC; ++ic) {
        int iR = BqF->get(ir,ic,0);
        int iRn = -1;
        for(int i = 0; i < ks_p.size(); ++i ) {
          if( ks_p.at(i).second == iR ) {
            iRn = i;
            break;
          }
        }
        BqFn->set(ir,ic,0,iRn);
      }
    }
    //BqFn->writeImage("labels2.txt");

   // for(int i = 0; i < displacers.size(); ++i) {
   //   if(displacers.at(i) != 0) {
   //     delete displacers.at(i);
   //   }
   // }

    vector<int> sort_indices;
    for(int i = 0; i < ks_p.size(); ++i ) {
      sort_indices.push_back(ks_p.at(i).second-1);
    }
    for(int i = 0; i < 225; ++i ) {
      displacers.at(i)->reorder_displacements(sort_indices);
    }

    return BqFn;
}

// JPEG 2000 stuff
/**
sample error debug callback expecting no client object
*/
static void error_callback(const char *msg, void *client_data) {
	(void)client_data;
	fprintf(stdout, "[ERROR] %s", msg);
}
/**
sample warning debug callback expecting no client object
*/
static void warning_callback(const char *msg, void *client_data) {
	(void)client_data;
	fprintf(stdout, "[WARNING] %s", msg);
}
/**
sample debug callback expecting no client object
*/
static void info_callback(const char *msg, void *client_data) {
	(void)client_data;
	//fprintf(stdout, "[INFO] %s", msg);
}

void compress_with_jp2(Image<uint8_t>* A_ij, string filename, int precision) {
   int NR = A_ij->size().height;
   int NC = A_ij->size().width;
  // JP2 stuff
  opj_cparameters_t l_param;
  opj_set_default_encoder_parameters(&l_param);
  l_param.tcp_numlayers = 1;
  l_param.cp_fixed_quality = 1;
  //l_param.tcp_distoratio[0] = 20;
  /* tile definitions parameters */
  /* position of the tile grid aligned with the image */
  l_param.cp_tx0 = 0;
  l_param.cp_ty0 = 0;
  /* tile size, we are using tile based encoding */
  l_param.tile_size_on = OPJ_TRUE;
  l_param.cp_tdx = (NC);
  l_param.cp_tdy = (NR);

  /** number of resolutions */
  //l_param.numresolution = 6;

  /** progression order to use*/
  /** OPJ_LRCP, OPJ_RLCP, OPJ_RPCL, PCRL, CPRL */
  l_param.prog_order = OPJ_LRCP;



  opj_image_cmptparm_t l_params [3];
  opj_image_cmptparm_t * l_current_param_ptr;
  l_current_param_ptr = l_params;
  for (int i=0;i<3;++i) {
     l_current_param_ptr->dx = 1;
     l_current_param_ptr->dy = 1;

     l_current_param_ptr->h = (OPJ_UINT32)(NR);
     l_current_param_ptr->w = (OPJ_UINT32)(NC);

        l_current_param_ptr->sgnd = 0;
     l_current_param_ptr->prec = (OPJ_UINT32)precision;

     l_current_param_ptr->x0 = 0;
     l_current_param_ptr->y0 = 0;

     ++l_current_param_ptr;
  }

  opj_codec_t * l_codec = opj_create_compress(OPJ_CODEC_JP2);

  opj_image_t * l_image = opj_image_tile_create(3,l_params,OPJ_CLRSPC_SRGB);
  l_image->x0 = 0;
  l_image->y0 = 0;
  l_image->x1 = (OPJ_UINT32)(NC);
  l_image->y1 = (OPJ_UINT32)(NR);
  l_image->color_space = OPJ_CLRSPC_SRGB;

  opj_setup_encoder(l_codec,&l_param,l_image);

  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename.c_str(),OPJ_FALSE);

  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_start_compress(l_codec,l_image,l_stream);

  int l_data_size;
  OPJ_BYTE* l_data;
  if(precision==8) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    for( int icomp = 0; icomp < 3; ++icomp) {
      for(int iir = 0; iir < A_ij->size().height; ++iir) {
        for(int iic = 0; iic < A_ij->size().width ; ++iic) {
          l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
          ++i;
        }
      }
    }
  }
  if(precision==16) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)2; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    for(int iir = 0; iir < A_ij->size().height; ++iir) {
      for(int iic = 0; iic < A_ij->size().width ; ++iic) {
        uint16_t val = (uint16_t)A_ij->get(iir,iic,0);
        l_data[i] = (OPJ_BYTE)(val >> 8);
        l_data[i+1] = (OPJ_BYTE)(val & 0x00ff);
        i += 2;
      }
    }
  }
  int i=0;
  opj_write_tile(l_codec,i,l_data,l_data_size,l_stream);

  opj_end_compress(l_codec,l_stream);
  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(l_image);
  free(l_data);
}

void compress_with_jp2_BqF(Image<int>* A_ij, string filename, int precision) {
   int NR = A_ij->size().height;
   int NC = A_ij->size().width;
  // JP2 stuff
  opj_cparameters_t l_param;
  opj_set_default_encoder_parameters(&l_param);
  l_param.tcp_numlayers = 1;
  l_param.cp_fixed_quality = 1;
  //l_param.tcp_distoratio[0] = 20;
  /* tile definitions parameters */
  /* position of the tile grid aligned with the image */
  l_param.cp_tx0 = 0;
  l_param.cp_ty0 = 0;
  /* tile size, we are using tile based encoding */
  l_param.tile_size_on = OPJ_TRUE;
  l_param.cp_tdx = (NC);
  l_param.cp_tdy = (NR);

  /** number of resolutions */
  //l_param.numresolution = 6;

  /** progression order to use*/
  /** OPJ_LRCP, OPJ_RLCP, OPJ_RPCL, PCRL, CPRL */
  l_param.prog_order = OPJ_LRCP;



  opj_image_cmptparm_t l_params [1];
  opj_image_cmptparm_t * l_current_param_ptr;
  l_current_param_ptr = l_params;
     l_current_param_ptr->dx = 1;
     l_current_param_ptr->dy = 1;

     l_current_param_ptr->h = (OPJ_UINT32)(NR);
     l_current_param_ptr->w = (OPJ_UINT32)(NC);

     l_current_param_ptr->sgnd = 0;
     l_current_param_ptr->prec = (OPJ_UINT32)precision;

     l_current_param_ptr->x0 = 0;
     l_current_param_ptr->y0 = 0;

  opj_codec_t * l_codec = opj_create_compress(OPJ_CODEC_JP2);
  cout << "create compress" << endl;

  opj_image_t * l_image = opj_image_tile_create(1,l_params,OPJ_CLRSPC_GRAY);
  l_image->x0 = 0;
  l_image->y0 = 0;
  l_image->x1 = (OPJ_UINT32)(NC);
  l_image->y1 = (OPJ_UINT32)(NR);
  l_image->color_space = OPJ_CLRSPC_GRAY;

  opj_setup_encoder(l_codec,&l_param,l_image);
  cout << "setup encoder" << l_codec << endl;

  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename.c_str(),OPJ_FALSE);
  cout << "create stream " << l_stream <<  endl;

  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  cout << l_codec << l_image << l_stream << endl;
  opj_start_compress(l_codec,l_image,l_stream);
  cout << "start compress" << endl;

  int l_data_size;
  OPJ_BYTE* l_data;
  if(precision==8) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)3; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    for( int icomp = 0; icomp < 3; ++icomp) {
      for(int iir = 0; iir < A_ij->size().height; ++iir) {
        for(int iic = 0; iic < A_ij->size().width ; ++iic) {
          l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
          ++i;
        }
      }
    }
  }
  if(precision==16) {
    l_data_size = (OPJ_UINT32)(NC)* (OPJ_UINT32)(NR)* (OPJ_UINT32)2; 
    l_data = (OPJ_BYTE*) malloc(l_data_size * sizeof(OPJ_BYTE));
    int i = 0;
    for(int iir = 0; iir < A_ij->size().height; ++iir) {
      for(int iic = 0; iic < A_ij->size().width ; ++iic) {
        uint16_t val = (uint16_t)A_ij->get(iir,iic,0);
        l_data[i] = (OPJ_BYTE)(val >> 8);
        l_data[i+1] = (OPJ_BYTE)(val & 0x00ff);
        i += 2;
      }
    }
  }
  int i=0;
  opj_write_tile(l_codec,i,l_data,l_data_size,l_stream);

  opj_end_compress(l_codec,l_stream);
  opj_stream_destroy(l_stream);
  opj_destroy_codec(l_codec);
  opj_image_destroy(l_image);
  free(l_data);
}

Image<int>* decompress_with_jp2_BqF(string filename) {
   opj_image_t * l_image;
  opj_dparameters_t l_param;
  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename.c_str(),OPJ_TRUE);
  /* Set the default decoding parameters */
  opj_set_default_decoder_parameters(&l_param);

  /* */
  //l_param.decod_format = JP2_CFMT; 

  /** you may here add custom decoding parameters */
  /* do not use layer decoding limitations */
  l_param.cp_layer = 0;

  /* do not use resolutions reductions */
  l_param.cp_reduce = 0;
  opj_codec_t* l_codec = opj_create_decompress(OPJ_CODEC_JP2);
  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_setup_decoder(l_codec, &l_param);
  opj_read_header(l_stream, l_codec, &l_image);
  OPJ_UINT32 da_x0=0;
  OPJ_UINT32 da_y0=0;
  OPJ_UINT32 da_x1=l_image->x1;
  OPJ_UINT32 da_y1=l_image->y1;
  opj_set_decode_area(l_codec, l_image, da_x0, da_y0,da_x1, da_y1);
  OPJ_UINT32 l_tile_index, l_data_size, l_nb_comps;
  OPJ_INT32 l_current_tile_x0, l_current_tile_y0, l_current_tile_x1, l_current_tile_y1, l_go_on;
  opj_read_tile_header( l_codec,
        l_stream,
        &l_tile_index,
        &l_data_size,
        &l_current_tile_x0,
        &l_current_tile_y0,
        &l_current_tile_x1,
        &l_current_tile_y1,
        &l_nb_comps,
        &l_go_on);

  OPJ_BYTE* l_data = (OPJ_BYTE *) malloc(l_data_size);
  opj_decode_tile_data(l_codec,l_tile_index,l_data,l_data_size,l_stream);
        /** now should inspect image to know the reduction factor and then how to behave with data */
  int i = 0;
  cout << (int)l_data[0] << endl;
  //cout << (int)l_data[0] << " ";
  Image<int>* A = new Image<int>(l_image->y1,l_image->x1,1);
  for(int iir = 0; iir < A->size().height; ++iir) {
    for(int iic = 0; iic < A->size().width ; ++iic) {
      //l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
      uint8_t val1 = l_data[i];
      uint8_t val2 = l_data[i+1];
      A->set(iir,iic,0,(val1<<8)+val2);
      i+=2;
    }
  }

  free(l_data);

  return A;

}

Image<uint8_t>* decompress_with_jp2(string filename) {
   opj_image_t * l_image;
  opj_dparameters_t l_param;
  opj_stream_t* l_stream = opj_stream_create_default_file_stream(filename.c_str(),OPJ_TRUE);
  /* Set the default decoding parameters */
  opj_set_default_decoder_parameters(&l_param);

  /* */
  //l_param.decod_format = JP2_CFMT; 

  /** you may here add custom decoding parameters */
  /* do not use layer decoding limitations */
  l_param.cp_layer = 0;

  /* do not use resolutions reductions */
  l_param.cp_reduce = 0;
  opj_codec_t* l_codec = opj_create_decompress(OPJ_CODEC_JP2);
  opj_set_info_handler(l_codec, info_callback,00);
  opj_set_warning_handler(l_codec, warning_callback,00);
  opj_set_error_handler(l_codec, error_callback,00);
  opj_setup_decoder(l_codec, &l_param);
  opj_read_header(l_stream, l_codec, &l_image);
  OPJ_UINT32 da_x0=0;
  OPJ_UINT32 da_y0=0;
  OPJ_UINT32 da_x1=l_image->x1;
  OPJ_UINT32 da_y1=l_image->y1;
  opj_set_decode_area(l_codec, l_image, da_x0, da_y0,da_x1, da_y1);
  OPJ_UINT32 l_tile_index, l_data_size, l_nb_comps;
  OPJ_INT32 l_current_tile_x0, l_current_tile_y0, l_current_tile_x1, l_current_tile_y1, l_go_on;
  opj_read_tile_header( l_codec,
        l_stream,
        &l_tile_index,
        &l_data_size,
        &l_current_tile_x0,
        &l_current_tile_y0,
        &l_current_tile_x1,
        &l_current_tile_y1,
        &l_nb_comps,
        &l_go_on);

  OPJ_BYTE* l_data = (OPJ_BYTE *) malloc(l_data_size);
  opj_decode_tile_data(l_codec,l_tile_index,l_data,l_data_size,l_stream);
        /** now should inspect image to know the reduction factor and then how to behave with data */
  int i = 0;
  cout << (int)l_data[0] << endl;
  //cout << (int)l_data[0] << " ";
  Image<uint8_t>* A = new Image<uint8_t>(l_image->y1,l_image->x1,3);
  for( int icomp = 0; icomp < 3; ++icomp) {
     for(int iir = 0; iir < A->size().height; ++iir) {
     //for(int iir = A->size().height-1; iir > 0; --iir) {
        for(int iic = 0; iic < A->size().width ; ++iic) {
           //l_data[i] = (OPJ_BYTE)A_ij->get(iir,iic,icomp); /*rand();*/
           A->set(iir,iic,icomp,l_data[i]);
           ++i;
        }
     }
  }

  free(l_data);

  return A;

}

// callbacks for mean shift
void bgLogVar(const char* msg, va_list list) {
} 
bool stop_flag = false;
int percentDone = 0;
// </ meanshift >

Image<int>* segment(Image< uint8_t >* A, const int h_s, const int h_r) {
  int NR = A->size().height;
  int NC = A->size().width;
  uint8_t* img = (uint8_t*)malloc((NR*NC*3)*sizeof(uint8_t));
  int ir = 0;
  for(int i = 0; i < NR; ++i ) {
    for(int j = 0; j < NC; ++j ) {
      for(int icomp = 0; icomp < 3; ++icomp) {
        img[ir] = A->get(i,j,icomp);
        ++ir;
      }
    }
  }
  msImageProcessor ms;
  ms.DefineImage(img,COLOR,NR,NC);
  ms.Segment(h_s,h_r,256,NO_SPEEDUP);
  
  int* labelshere = new int[NR*NC];
  float* modeshere = 0;
  int* modePointCountshere = 0;
  int return_regions = ms.GetRegions(&labelshere, &modeshere, &modePointCountshere);
  //cout << "regions " << return_regions << endl;
  Image<int>* BqFi = new Image<int>(NR,NC,1);
  ir = 0;
  for(int i = 0; i < NR; ++i ) {
    for(int j = 0; j < NC; ++j ) {
      BqFi->set(i,j,0,labelshere[ir]+1);
      ++ir;
    }
  }
  //BqFi->writeImage("labels.txt");

  delete[] labelshere;
  delete[] modeshere;
  delete[] modePointCountshere;
  free(img);

  return BqFi;
}

void cerv_compress(Image< int >* BqF, string partition_fn, string partition_labels_fn) {
  int NR = BqF->size().height;
  int NC = BqF->size().width;
    vector<int> ind2label(1000,-50);
    int** img = (int**)malloc((NR)*sizeof(int*));
    for(int i = 0; i < BqF->size().height; ++i ) {
      img[i] = (int*)malloc((NC)*sizeof(int));
      for(int j = 0; j < BqF->size().width; ++j ) {
        img[i][j] = BqF->get(i,j,0)+2;
      }
    }
    char fn[500];
    strcpy(fn, partition_fn.c_str());
    cerv_encode(img, NR, NC, fn);

    int** img2 = (int**)malloc((NR)*sizeof(int*));
    for(int i = 0; i < BqF->size().height; ++i ) {
      img2[i] = (int*)malloc((NC)*sizeof(int));
      for(int j = 0; j < BqF->size().width; ++j ) {
        img2[i][j] = 0;
      }
    }

    cerv_decode(img2, NR, NC, fn);

    int L_max = -5;
    for(int i = 0; i < BqF->size().height; ++i ) {
      for(int j = 0; j < BqF->size().width; ++j ) {
        ind2label.at(img2[i][j]) = BqF->get(i,j,0); 
        if( img2[i][j] > L_max ) {
          L_max = img2[i][j];
        }
      }
    }
    FILE* partition_label_file = fopen(partition_labels_fn.c_str(), "wb");
    for( int i = 0; i <= L_max; ++i ) {
      short temp = ind2label.at(i);
      fwrite(&temp, sizeof(short), 1, partition_label_file);
    }
    fclose(partition_label_file);

    for(int i = 0; i < BqF->size().height; ++i ) {
      free(img2[i]); 
    }
    free(img2);


}

void cerv_decompress(Image< int >* BqFd, string partition_fn, string partition_labels_fn) {
  int NR = BqFd->size().height;
  int NC = BqFd->size().width;

  int** img2 = (int**)malloc((NR)*sizeof(int*));
  for(int i = 0; i < BqFd->size().height; ++i ) {
    img2[i] = (int*)malloc((NC)*sizeof(int));
    for(int j = 0; j < BqFd->size().width; ++j ) {
      img2[i][j] = 0;
    }
  }
  char fn[500];
  strcpy(fn, partition_fn.c_str());
  cerv_decode(img2, NR, NC, fn);

  vector< int > ind2label;
  FILE* partition_label_file = fopen(partition_labels_fn.c_str(), "rb");
  short temp = 0;
  while(fread(&temp, sizeof(short), 1, partition_label_file)==1) {
    ind2label.push_back(temp);
  }
  fclose(partition_label_file);

  for(int i = 0; i < BqFd->size().height; ++i ) {
    for(int j = 0; j < BqFd->size().width; ++j ) {
      int this_label = ind2label.at(img2[i][j]);
      BqFd->set(i,j,0,this_label);
    }
  }
  //BqFd->writeImage("BqF_decoded.txt");

    for(int i = 0; i < BqFd->size().height; ++i ) {
      free(img2[i]); 
    }
    free(img2);
}

void readBqF(const int NR, const int NC, Image< uint8_t >& BqF, string filename) {
  int temp = 0;
  FILE* THIS_IMAGE = fopen(filename.c_str(), "rb");
  for(int i=0;i<NR;i++) {
    for( int k = 0; k < NC; ++k ) {
      fread(&temp, sizeof(int), 1, THIS_IMAGE);
      BqF.set(i,k,0,temp);
    }
  }
  fclose(THIS_IMAGE);
}


void readImage(const int NR, const int NC, Image< uint8_t >& Ydecoded, string imagename, const int image_index_i, const int image_index_j) {
  double temp = 0;
  ostringstream stringStream;
  stringStream << imagename << "_" << image_index_i << "_" << image_index_j;
  string imagefile = stringStream.str();
  stringStream.str("");
  FILE* THIS_IMAGE = fopen(imagefile.c_str(), "rb");

  for(int j=0;j<3;j++) {
    for(int i=0;i<NR;i++) {
      for( int k = 0; k < NC; ++k ) {
        fread(&temp, sizeof(double), 1, THIS_IMAGE);
        Ydecoded.set(i,k,j,uint8_t(temp));
      }
    }
  }
  fclose(THIS_IMAGE);
}


bool create_zip(vector<string>& files_to_archive, string zip_filename) {
  //const char* pZip_filename = "kakka.zip";
  mz_zip_archive zip;
  memset(&zip, 0, sizeof(zip));
  if (!mz_zip_writer_init_file(&zip, zip_filename.c_str(), 0))
  {
    return false;
  }

  for(int i = 0; i < files_to_archive.size(); ++i) {
    mz_bool success = mz_zip_writer_add_file(&zip, files_to_archive.at(i).c_str(), files_to_archive.at(i).c_str(), 0, 0, MZ_BEST_SPEED);
    if (!success)
    {
      mz_zip_writer_end(&zip);
      remove(zip_filename.c_str());
      return false;
    }
  }

  if (!mz_zip_writer_finalize_archive(&zip))
  {
    mz_zip_writer_end(&zip);
    remove(zip_filename.c_str());
    return false;
  }

  mz_zip_writer_end(&zip);

  printf("Created zip file %s\n", zip_filename.c_str()); 
  return true;
}

bool zip_extract(string zip_filename)
{
  mz_zip_archive zip;
  memset(&zip, 0, sizeof(zip));
  if (!mz_zip_reader_init_file(&zip, zip_filename.c_str(), 0))
  {
    return false;
  }

  for (int i = 0; i < mz_zip_reader_get_num_files(&zip); i++)
  {
    char pFilename[400];
    mz_uint bytes_written =  mz_zip_reader_get_filename(&zip, i, pFilename, 400);
    mz_zip_reader_extract_to_file(&zip, i, pFilename, 0);
    printf("Extracted file \"%s\"\n", pFilename);
  }

  mz_zip_reader_end(&zip);

  return true;
}
