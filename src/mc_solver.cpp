#include "mc_solver.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
	int t = clock();
        static const std::size_t length = pdict.num_paired;
        map_matrix<int> num_matrix_third(length);
        cout << "========" << endl;
        cout << "monte carlo solver: " << endl;
        length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
        vector<bool> flip_maxima(length); 
	for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
        for (int i = 0; i < solver_loops; i++) {
                cout << "******************* solver loop     " << i << endl;
                t = clock();
                single_spin_flip_map( pdict, diff_matrix, num_matrix_third );  // num_matrix_second
                if (i > 0) { get_flip_positions( pdict, flip_maxima ); }
                t = clock() - t;
                cout << "spin flip --------- time: " << t << endl;
                t = clock();
                cout << "-- finished flip spins     " << endl;
                block_flip_map( pdict, diff_matrix, num_matrix_third, flip_maxima );  // num_matrix_second
                t = clock() - t;
                cout << "block flip -------- time: " << t << endl;
                t = clock();
                cout << "-- finished block flip     " << endl;
                energy_sum_min_max( pdict );
        }
        set_switch_energy( pdict, flip_maxima );   /// create a new haplotype block class and find the segments
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
        int t = clock();
	static const std::size_t length = pdict.num_paired;
	map_matrix<int> num_matrix_third(length);
        cout << "========" << endl;
        cout << "monte carlo solver: " << endl;
	length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
        //vector<bool> flip_maxima(pdict.num_paired);
        //for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
        for (int i = 0; i < solver_loops; i++) {
                cout << "******************* solver loop     " << i << endl;
                t = clock();
		single_spin_flip_recursive( pdict, diff_matrix, num_matrix_third );
                t = clock() - t;
                cout << "spin flip --------- time: " << t << endl;
                t = clock();
		block_flip_recursive( pdict, diff_matrix, num_matrix_third );
                t = clock() - t;
                cout << "block flip -------- time: " << t << endl;
                energy_sum_min_max( pdict );
        }
        //set_switch_energy(pdict,flip_maxima);   /// create a new haplotype block class and find the segments
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;      //cout << i << " " << m << " " << diff_pos << " " << nmatrix(i,m) << endl;
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        if (diff_pos < pos_diff_cutoff) { nmatrix2.add_to(i,m,nmatrix(i,m)); }
                }
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void call_blocks(coord_dictionary& pdict,std::string technology) {
	int switch_cutoff = -200;
	int block_num = 0;
	if ( technology == "tenx" ) { switch_cutoff = -200; }
	else if ( technology == "pacbio" ) { switch_cutoff = -50; }
        for (int i = 0; i < pdict.num_paired; i++) {
		if (pdict.switchE[i] > switch_cutoff) { block_num += 1; }
		pdict.block[i] = block_num;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                if (!flip_maxima[i]) { pdict.switchE[i] = -10000.0; }
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> neighbors{i-1,i+1};  //cout << "  " << i << "  ";
                double energy = pdict.deltaE[i];  bool check = true;
		if ( energy < -200 ) {
                for (int j=0; j < neighbors.size(); j++) {   //cout << neighbors[j] << "  ";
                        if ( neighbors[j] >= 0 && neighbors[j] < pdict.num_paired ) {
				if ( pdict.deltaE[neighbors[j]] > energy ) { check = false; } 
			}
                } }
                flip_maxima[i] = check;
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void energy_sum_min_max( coord_dictionary& pdict ) {
        double sum = 0.0; double max = 0.0; double min = 0.0;
        for (int i = 0; i < pdict.num_paired; i++) {
                sum += pdict.deltaE[i];
                if (pdict.deltaE[i] > max) { max = pdict.deltaE[i]; }
                if (pdict.deltaE[i] < min) { min = pdict.deltaE[i]; }
        }
        cout << "-- energy sum: " << sum << " minimum: " << min << " maximum: " << max << endl;
}

// added negative signs
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int init_i = pdict.haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        initial_energy += -1.0*init_i*pdict.haplotype[m]*diff_matrix(i,m);  //final_energy += final_i*pdict.haplotype[m]*diff_matrix(i,m);
                        //cout << i << " " << m << " " << pdict.haplotype[i] << " " << diff_pos << " " << nmatrix(i,m) << " " << diff_matrix(i,m) << endl;
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                //cout << "init e " << initial_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
                }
        }
        double energyA = initial_energy;
        double energyB = initial_energy;
        vector<int> loop_haplotype = pdict.haplotype;
        for (int i = 0; i < pdict.num_paired; i++) {
                double save_energyB = energyB;
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                int hap_spin  = loop_haplotype[i];
                int flip_spin = -1*loop_haplotype[i];
                double spin_i_energy = 0.0;
                double spin_f_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        spin_i_energy += -1.0*hap_spin*switch_haplotype[m]*diff_matrix(i,m);  /// may need to add negative sign
                        spin_f_energy += -1.0*flip_spin*switch_haplotype[m]*diff_matrix(i,m);
                }
                double switch_energy = spin_f_energy - spin_i_energy;
                energyA = energyA + switch_energy;
                double diff_energy = energyA - energyB;
                if ( diff_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++;
                }
                else {
                        energyB = save_energyB;
                        energyA = energyA;
                }
                pdict.switchE[i] = (-1.0)*diff_energy;
        }
        pdict.haplotype = loop_haplotype;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0; clock_t t;
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int l = 0; bool cross_i = true; int spin_pos;  double init_energy,final_energy,diff_energy;
                for(auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        if( m>i && cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; l++; }
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[m] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << m << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[m] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[m] << "\t double_pos " << diff_pos2 << endl;      //////////////////////////
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(m); l++;}
                }
                if(cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                submatrix_opt tempsub;  tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                //cout << "flip  position: " << i << " sparse_map_size " << sparse_map.size() << " spin_pos " << spin_pos << "  total " << l << endl;
                final_sub_haplotype[spin_pos] = final_sub_haplotype[spin_pos]*(-1);
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
        cout << "-- number of spin flips: " << nsf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima ) {
        int nbf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                if (flip_maxima[i]) {                   //if (pdict.switchE[i] > -switch_cut) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int flip_pos,l=0;   double init_energy,final_energy,diff_energy;
                for ( int j = pdict.flip_low_bound[i]; j <= pdict.flip_up_bound[i]; j++ ) {
                        if (j == i) { flip_pos = l; }
                        int diff_pos = abs(pdict.sorted_paired_positions[j] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[j] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << j << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[j] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[j] << "\t double_pos " << diff_pos2 << endl;
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(j); l++; }
                }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                //cout << "block position: " << i << " sparse_map_size " << sparse_map.size() << " flip_pos " << flip_pos << "  total " << l << endl;
                submatrix_opt tempsub; tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                for (int j = 0; j < flip_pos; j++) { final_sub_haplotype[j] = final_sub_haplotype[j]*(-1); }
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { for (int l = 0; l < i; l++) { pdict.haplotype[l] = pdict.haplotype[l]*(-1); } nbf++; }
                pdict.switchE[i] = (-1.0)*diff_energy;
                }    //}
        }
        cout << "-- number of block flips: " << nbf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype ) {
        vector<double> temp_vector = dot_product_matrix(haplotype,tempsub.smat);
        double energy_return = (-0.5)*dot_product(haplotype,temp_vector);
        return energy_return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double dot_product( vector<int> a, vector<double> b ) {
        double dot = 0.0; for (int i = 0; i < a.size(); i++) { dot += (double)(a[i])*(b[i]); };
        return dot;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b ) {
        vector<double> dot_vector;
        for (int i = 0; i < a.size(); i++) {
                double i_value = 0.0; for (int j = 0; j < b.size(); j++) { i_value += (double)a[j]*b[i][j]; }
                dot_vector.push_back(i_value);
        }
        return dot_vector;
};


