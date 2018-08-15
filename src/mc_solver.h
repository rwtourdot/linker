#ifndef MC_SOLVER_H
#define MC_SOLVER_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "map_matrix.h"
#include "sub_matrix.h"

///////////// solver_cutoffs //////
#define solver_loops 30                 // 10  // {val} the number of spin flip block flip loops
#define pos_diff_cutoff 100000          // {val} if the maximum delta genome distance - band width

/////////////// functions /////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict,map_matrix<double> diff_matrix,map_matrix<int> num_matrix_second ); ///// fix this
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void energy_sum_min_max( coord_dictionary& pdict );
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima );
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype );
static double dot_product( vector<int> a, vector<double> b );
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b );
void call_blocks(coord_dictionary& pdict, std::string technology);

/////////////////
void solver_recursive(std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second );
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 );
static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );

#endif  // MC_SOLVER_H
