#include "write_linker_output.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(coverageFile);
        for (auto& it : var_dict) {
		if ( it.second.var ) {
                ofile << it.first << "\t" << it.second.pos << "\t";
                ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
                for (auto& it2 : it.second.base_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
		ofile << it.second.total_bases << "\t";
                ofile << endl;
		}
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_bx_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice ) {
	// print each bx tag
        ofstream ofile; ofile.open(coverageFile);
        for (auto& it : var_dict) {
                if ( it.second.var ) {
                char rbase = it.second.ref_base[0]; char vbase = it.second.var_base[0];
                std::set<std::string> ref_bx_set = it.second.base_dict_set[rbase];
                std::set<std::string> var_bx_set = it.second.base_dict_set[vbase];
                ofile << it.first << "\t" << it.second.pos << "\t";
                ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
                for (auto& it2 : it.second.unique_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
                ofile << it.second.unique_total << "\t";
		ofile << it.second.ref_base << ":";
		for ( std::set<std::string>::iterator it4=ref_bx_set.begin(); it4!=ref_bx_set.end(); ++it4 ) { 
			ofile << *it4;
			if (std::distance(ref_bx_set.begin(), it4) < (ref_bx_set.size()-1) ) { ofile << ","; }
		}
		ofile << "\t";
		ofile << it.second.var_base << ":";
                for ( std::set<std::string>::iterator it5=var_bx_set.begin(); it5!=var_bx_set.end(); ++it5 ) { 
			ofile << *it5;
			if (std::distance(var_bx_set.begin(), it5) < (var_bx_set.size()-1) ) { ofile << ","; }
		}
                ofile << endl;
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_coverage_filter( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile ) {
        ofstream ofile; ofile.open(coverageFile);
        for (auto& it : var_dict) { 
		//cout << it.first << "  " << it.second.filter << "  " << it.second.var << endl;
                if ( it.second.var ) {
		if ( it.second.filter ) {
        		ofile << it.first << "\t" << it.second.pos << "\t";
        		ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
        		for (auto& it2 : it.second.base_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
        		ofile << it.second.total_bases << "\t";
        		ofile << endl;
		}
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hic_links( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                int pos1 = it.second.pos;
                std::string base1 = it.first;
                if (it.second.unique_total > 0) {
                for (auto it2 : it.second.connections) {  //int pos2 = 100;  //var_dict[it2.first].pos;
                        int pos2 = var_dict[it2.first].pos;
			int nconnections = it2.second;
                        std::string base2 = it2.first;
			std::vector<std::string> v1 = var_dict[it2.first].connected_reads_long_form;
			std::vector<std::string> v2 = var_dict[it.first].connected_reads_long_form;
                        ofile << pos1 << "\t" << base1 << "\t" << pos2 << "\t" << base2 << "\t" << nconnections << "\t";
			bool first = true; ofile << "[";
			for (int j = 0; j < v1.size(); j++) {
				if ( std::find(v2.begin(), v2.end(), v1[j]) != v2.end() ) {
					if (!first) { ofile << ','; } 
					first = false; 
					ofile << v1[j]; 
				}
			}
			ofile << "]" << endl;
                }
                }
        }
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_link_network( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice ) {   // cout << " writing network file " << endl;
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                ofile << chr_choice << "_" << it.first << "  " << it.second.pos << "  ["; //<< endl;
                bool first = true;
                for (auto it2 : it.second.connections) { if (!first) { ofile << ','; } first = false; ofile << chr_choice << "_" << it2.first; }
                ofile << "] [";
                bool first_num = true;
                for (auto it2 : it.second.connections) { if (!first_num) { ofile << ','; } first_num = false; ofile << it2.second; }
                ofile << "] {"; std::string base_output;
                for (auto it2 : it.second.base_dict) {
                        if(!base_output.empty()) { base_output += ","; }
                        string ssbase(1,it2.first);
                        base_output += ssbase + ":" + std::to_string(it2.second);
                }
                ofile << base_output << "} paired " << std::boolalpha << it.second.paired << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hap_solution( std::unordered_map<std::string,variant_node>& var_dict, std::string hapsolutionFile, coord_dictionary& pdict, std::string chr_choice ) {
        ofstream ofile; ofile.open(hapsolutionFile);
        std::unordered_map<int,int> opposite; opposite[0] = 1; opposite[1] = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int p = pdict.sorted_paired_positions[i];
                int m = pdict.ref_index[p];
                ofile << i << "\t" << p << "\t" << chr_choice << "_" << pdict.paired_dict[p][opposite[m]] << "\t" << chr_choice << "_" << pdict.paired_dict[p][m] << "\t" << pdict.haplotype[i] << "\t" << pdict.base_number[p][opposite[m]] << "\t" << pdict.base_number[p][m] << "\t" << pdict.deltaE[i] << "\t" << pdict.switchE[i] << "\t" << pdict.block[i] << "\t" << pdict.span_bound[i] << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_matrix( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : bx_map) {
                for (auto& it2 : it.second.connections) { 
			ofile << contig_name << "\t" << it.first*binsize << "\t" << "-" << "\t" << contig_name << "\t" << it2.first*binsize << "\t" << "+" << "\t" << it2.second << endl; 
		}
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_bx_cov( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : bx_map) {
        	ofile << contig_name << "\t" << it.first*binsize << "\t" << it.second.num_bx << "\t" << it.second.num_unique_bx << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_matrix_genome( full_map& chr_map, std::string outputFile, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : chr_map) {
                for (auto& it2 : it.second) {
                	for (auto& it3 : it2.second.chr_connections) {
				for (auto& it4 : it3.second) {
                        		ofile << it.first << "\t" << it2.first*binsize << "\t" << "-" << "\t" << it3.first << "\t" << it4.first*binsize << "\t" << "+" << "\t" << it4.second << endl;
				}
                	}
        	}
	}
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_phased_sv( std::vector<sv_entry>& sv_list, std::string outputFile ) {
	ofstream ofile; ofile.open(outputFile);
        for (int i = 0; i < sv_list.size(); i++) {
                ofile << sv_list[i].chr_one << ":" << sv_list[i].pos_one << "  " << sv_list[i].chr_two << ":" << sv_list[i].pos_two << endl;
                for (auto& it : sv_list[i].hap_end1) { ofile << "0 " << sv_list[i].chr_one << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end1[it.first] << endl; }
                for (auto& it : sv_list[i].hap_end2) { ofile << "1 " << sv_list[i].chr_two << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end2[it.first] << endl; }
                ofile << "######" << endl;
        }
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_cn_phased( cn_map& chromosome_map, std::string outputFile, std::vector<int>& good_bins, std::vector<int>& merged_bins ) {
	ofstream ofile; ofile.open(outputFile);
	cout << "write output " << endl;
	int k = 0;
	int p = 0;
	for (auto ind : good_bins) {
		for (int j = 0; j < chromosome_map[ind].num_hets; j++) {
			int hap = chromosome_map[ind].switch_spin*chromosome_map[ind].het_hap[j];
			int orig_hap = chromosome_map[ind].het_hap[j];
			int hapA_cov = 0; 
			int hapB_cov = 0;
			int orig_hapA_cov = 0; 
			int orig_hapB_cov = 0; 
			if (hap == 1) {  hapA_cov = chromosome_map[ind].ref_cov[j]; hapB_cov = chromosome_map[ind].var_cov[j]; }
			if (hap == -1) { hapB_cov = chromosome_map[ind].ref_cov[j]; hapA_cov = chromosome_map[ind].var_cov[j]; }
                        if (orig_hap == 1) {  orig_hapA_cov = chromosome_map[ind].ref_cov[j]; orig_hapB_cov = chromosome_map[ind].var_cov[j]; }
                        if (orig_hap == -1) { orig_hapB_cov = chromosome_map[ind].ref_cov[j]; orig_hapA_cov = chromosome_map[ind].var_cov[j]; }
			ofile << k << "\t" << ind << "\t" << chromosome_map[ind].het_pos[j] << "\t" << chromosome_map[ind].ref_base[j] << "\t" << chromosome_map[ind].var_base[j] << "\t" << hap << "\t" << hapA_cov <<  "\t" << hapB_cov << "\t" << orig_hap << "\t" << orig_hapA_cov <<  "\t" << orig_hapB_cov << "\t" << chromosome_map[ind].het_block[j] << "\t" << chromosome_map[ind].switch_spin << "\t" << merged_bins[p] << endl;
			k += 1;
		}
		p+=1;
	}
	cout << "number of bins: " << p << endl;
	cout << "number of hets: " << k << endl;
        ofile.close();
}














                        //cout << pos1 << "\t" << base1 << "\t" << pos2 << "\t" << base2 << "\t" << nconnections << endl;
                        //cout << pos1 << "\t" << pos2 << "\t" << " reads1" << endl;
                        //for (int j = 0; j < v1.size(); j++) { cout << v1[j] << " " << v1[j].length() << endl;}
                        //for ( std::vector<std::string>::iterator it5 = v1.begin(); it5 != v1.end(); ++it5 ) { cout << *it5 << " " << *it5.length() << endl; }
                        //cout << pos1 << "\t" << pos2 << "\t" << " reads2" << endl;
                        //for (int k = 0; k < v2.size(); k++) { cout << v2[k] << " " << v2[k].length() << endl;}
                        //for ( std::vector<std::string>::iterator it6 = v2.begin(); it6 != v2.end(); ++it6 ) { cout << *it6 << " " << *it6.length() << endl; }







                        //bool first = true;
                        //cout << "===" << endl;
                        //for (auto rn : read_names) { if (!first) { ofile << ','; } first = false; ofile << rn; }  //cout << rn << ",";





