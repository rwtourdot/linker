#ifndef WRITE_LINKER_OUTPUT_H
#define WRITE_LINKER_OUTPUT_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "variant_site.h"
#include "bin_assembly.h"
#include "bin_reference.h"

/////////////// functions /////////////////////
void write_bin_matrix( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize );
void write_het_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice );
void write_hap_solution( std::unordered_map<std::string,variant_node>& var_dict, std::string hapsolutionFile, coord_dictionary& pdict, std::string chr_choice );
void write_link_network( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice );
void write_phased_sv( std::vector<sv_entry>& sv_list, std::string outputFile );
void write_bin_matrix_genome( full_map& chr_map, std::string outputFile, int binsize );
void write_cn_phased( cn_map& chromosome_map, std::string outputFile, std::vector<int>& good_bins, std::vector<int>& merged_bins );
void write_het_coverage_filter( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile );
void write_het_bx_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice );
void write_bin_bx_cov( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize );
void write_hic_links( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice );

#endif  // WRITE_LINKER_OUTPUT_H
