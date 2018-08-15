#include "read_linker_output.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool is_file_exist(std::string filename ) {
	std::ifstream infile(filename);
	return infile.good();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string split_string_first( std::string s, std::string delimiter, int choice ) {
	std::string return_string;
	if ( choice == 0 ) { return_string = s.substr(0, s.find(delimiter)); }
	if ( choice == 1 ) { return_string = s.substr(s.find(delimiter)+1, s.length()); }
	return return_string;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_het_coverage( std::string coverageFile, variant_graph& vgraph ) { 
	if (!is_file_exist(coverageFile)) { cout << " file not found " << coverageFile << endl; exit(1); }
        std::ifstream inputFile(coverageFile); std::string line;
	while (getline(inputFile,line)) {
		std::istringstream ss(line);
		int pos,tot_cov;
		std::string het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase;	
		ss >> het_name >> pos >> ref_var >> ibase >> dbase >> gbase >> cbase >> abase >> tbase >> tot_cov;
		//cout << coverageFile << "\t" << pos << " " << abase << " " << tbase << " " << cbase << " " << gbase << endl;
		std::string ref = split_string_first(ref_var,":",0);
		std::string var = split_string_first(ref_var,":",1);
		int anum = std::stoi(split_string_first(abase,"|",1));
		int tnum = std::stoi(split_string_first(tbase,"|",1));
		int cnum = std::stoi(split_string_first(cbase,"|",1));
		int gnum = std::stoi(split_string_first(gbase,"|",1));
		int dnum = std::stoi(split_string_first(dbase,"|",1));
		int inum = std::stoi(split_string_first(ibase,"|",1));
		std::string variant_id = std::to_string(pos)+"_"+ref+"_"+var;
		std::string reference_id = std::to_string(pos)+"_"+ref+"_"+ref;
		//cout << pos << " " << anum << " " << tnum << " " << cnum << " " << gnum << endl;
		//cout << variant_id << " " << reference_id << " " << het_name << "  " << pos << "  " << ref_var << "  " << ref << "  " << var << "  " << tot_cov << "  " << anum << endl;
                variant_node v_node1,v_node2;
                vgraph[variant_id] = v_node1;  
		vgraph[reference_id] = v_node2;
                vgraph[variant_id].reset_values(pos,true,var,ref,anum,tnum,cnum,gnum,dnum,inum);
                vgraph[reference_id].reset_values(pos,false,var,ref,anum,tnum,cnum,gnum,dnum,inum);
		vgraph[variant_id].paired = true;
		vgraph[reference_id].paired = true;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_hap_solution( std::string hapFile, coord_dictionary& pdict ) { 
        if (!is_file_exist(hapFile)) { cout << " file not found " << hapFile << endl; exit(1); }
	std::unordered_map<int,int> index_map;
	for (int i = 0; i < pdict.num_paired; i++) { index_map[pdict.sorted_paired_positions[i]] = i; }
        std::ifstream inputFile(hapFile); std::string line;
        while (getline(inputFile,line)) {
                std::istringstream ss(line);
                int index,pos,hap,hapA_cov,hapB_cov,block,width;
                std::string ref_name,het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase,spinE,switchE;
                ss >> index >> pos >> ref_name >> het_name >> hap >> hapA_cov >> hapB_cov >> spinE >> switchE >> block >> width;
                std::string chrom = split_string_first(ref_name,"_",0);
		int pind = index_map[pos];
		pdict.haplotype[pind] = hap;
		pdict.block[pind] = block;
		pdict.deltaE[pind] = std::stod(spinE);
		pdict.switchE[pind] = std::stod(switchE);
		pdict.reload_bool[pind] = true;   //cout << pdict.reload_bool[pind] << endl;
                cout << index << "  " << pind << "  " << pos << "  " << ref_name << "  " << chrom << "  " << hap << "  " << spinE << "  " << switchE << " "  << block << endl;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_sv_file( std::string svFile, std::vector<sv_entry>& sv_list ) { 
        if (!is_file_exist(svFile)) { cout << " file not found " << svFile << endl; exit(1); }
        std::ifstream inputFile(svFile); std::string line;
	int i=0;
        while (getline(inputFile,line)) {
		sv_entry sv_temp;
                std::istringstream ss(line);
                int raindex,pos1,str1,pos2,str2,totalcount;
                std::string chr1,chr2;
                ss >> raindex >> chr1 >> pos1 >> str1 >> chr2 >> pos2 >> str2 >> totalcount;
		if ( chr1.length() > 1 ) {
			if ( totalcount > min_sv_reads ) { 
				sv_temp.set_locations(chr1,pos1,chr2,pos2,str1,str2,totalcount);
				sv_list.push_back(sv_temp);
                		//cout << chr1 << "  " << pos1 << "  " << str1 << "  " << chr2 << "  " << pos2 << "  " << str2 << "  " << totalcount << endl;
			}
		}
		//if ( i > 40 ) { break; }
		i++;
        }
};

