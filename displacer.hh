#ifndef DISPLACER_HH
#define DISPLACER_HH

#include <vector>
#include <string>
#include "image.hh"
using std::vector;
using std::string;

class Displacer {

  public:
    Displacer(Image<uint8_t>* im1, Image<uint8_t>* im2, Image<int>* BqF);
    Displacer(Image<int>* BqF, string displacements_file);
    ~Displacer();
    void propagateRegions();
    vector<int>* get_displacements_h();
    void search_displacements();

    void write_displacements(string displacements_file); 
    void read_displacements(string displacements_file); 
    
    Image< int >* get_BqFSV();

  private:
    Image<uint8_t>* im1_;
    Image<uint8_t>* im2_;
    Image<int>* BqF_;
    Image<int>* BqFSV_;
    vector< int >* best_disp_h_;
    vector< int >* best_disp_v_;

    
};

#endif
