#ifndef READ_LINKER_OUTPUT_H
#define READ_LINKER_OUTPUT_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "variant_site.h"
#include "build_hap_matrix.h"
#include "sv_junction.h"

//////////////// definitions //////////////////
#define min_sv_reads 10

/////////////// functions /////////////////////
void read_het_coverage( std::string coverageFile, variant_graph& vgraph );
void read_hap_solution( std::string hapFile, coord_dictionary& pdict );
std::string split_string_first( std::string s, std::string delimiter, int choice);
bool is_file_exist( std::string filename );
void read_sv_file( std::string svFile, std::vector<sv_entry>& sv_list );

#endif  // READ_LINKER_OUTPUT

