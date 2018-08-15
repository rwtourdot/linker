#ifndef RUN_BX_BIN_H
#define RUN_BX_BIN_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "bin_assembly.h"
#include "bin_reference.h"
#include "read_bam.h"
#include "write_linker_output.h"

//////////////// definitions //////////////////
typedef std::unordered_map<std::string,int> contig_dict;
typedef std::unordered_map<std::string,contig_node> contig_bx;
//std::unordered_map<std::string,variant_node> variant_graph;
//std::unordered_map<std::string,read_tree> read_graph;
//std::vector<vcf_entry> vcf_vector;
//std::map<std::string,int> map_str_int;

/////////////// functions /////////////////////
static void parse_bx_bin_options( int argc, char** argv );
void run_bin_bxtag( int argc, char** argv );

#endif  // RUN_BX_BIN_H



