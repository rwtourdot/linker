#include "run_recover.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr20";
        static std::string input_scaffold_file = "./output/hap_full_scaffold_oct7_K562_chr20.dat";
        static std::string input_graph_file = "./output/graph_variant_oct3_K562_bam1_tenx_chr20.dat";
        static std::string id_string = "default";
        static std::string technology = "tenx";
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:i:g:c:n:";
static const struct option longopts[] = {
        { "help",        no_argument, NULL, 'h' },
        { "scaf-file",   no_argument, NULL, 'i' },
        { "graph-file",  no_argument, NULL, 'g' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "tech",  	 no_argument, NULL, 'e' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *RECOVER_USAGE_MESSAGE =
"Usage: linker recover [OPTION] -i /output/hap_full_scaffold_*chr4.dat -g /output/graph_variant_*chr4.dat -c chr4 \n\n"
"\n"
"********* this command is set up for hg38 *********"
"\n"
"  Options\n"
"  -i,      input scaffold path \n"
"  -g,      input graph file \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -e,      technology ( tenx, hic, nanopore ) \n"
"  -n,      id string for output files \n"
"\n";


///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_recover_options( int argc, char** argv ) {
        bool die = false; //bool vcf_load = false; //bool cov_load = false;  //if(argc <= 2) { die = true; }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << RECOVER_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_scaffold_file; break;
                case 'g': arg >> opt::input_graph_file; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'e': arg >> opt::technology; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
        //if (opt::input_bam_file.length() == 0) { die = true; }
        //if (opt::input_vcf_file.length() == 0) { die = true; }
        //if (!opt::vcf_load && !opt::cov_load) { die = true; }
        if (die) {
          std::cerr << "\n" << RECOVER_USAGE_MESSAGE;
          exit(1);
        }
        //parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running recover ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== input scaf file  === " << opt::input_scaffold_file << endl;
        cout << "== input graph file  === " << opt::input_graph_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "== technology  === " << opt::technology << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

void recover_links( std::string chr_choice, coord_dictionary& pdict, variant_graph& vgraph, std::string technology, coord_dictionary& rdict, std::string scaffoldsolutionFile ) {
	std::map<int,std::string> ref_position_map;
	std::map<int,std::string> alt_position_map;
        for (auto& it : vgraph) {
		if (it.second.var) {
			ref_position_map[it.second.pos] = it.second.ref_base;
			alt_position_map[it.second.pos] = it.second.var_base;
		}
        }
	std::unordered_map<int,std::string> ab_hap = {{1,"A"},{-1,"B"}};
	int max_delta_p;
	if ( technology == "tenx" ) { max_delta_p = 100000; }
	else if ( technology == "hic" ) { max_delta_p = 10000000; }
	else if ( technology == "nanopore" ) { max_delta_p = 1000000; }
	cout << "maximum delta pos: " << max_delta_p << endl;
        //#############################################
	ofstream ofile; ofile.open(scaffoldsolutionFile);
        //#############################################
	for (auto& it : ref_position_map) {
		int pos1 = it.first;
                int ref_numA,ref_numB,alt_numA,alt_numB;
                ref_numA = 0; ref_numB = 0; 
		alt_numA = 0; alt_numB = 0;
		std::string ref_base = ref_position_map[pos1];
		std::string alt_base = alt_position_map[pos1];
                std::string ref_hash = std::to_string(pos1) + "_" + ref_base + "_" + ref_base; // pdict.ref_handle[i];
                std::string alt_hash = std::to_string(pos1) + "_" + ref_base + "_" + alt_base;
		//std::string ref_hash = ref_position_map[pos1]; // pdict.ref_handle[i];
		//std::string alt_hash = alt_position_map[pos1];
		cout << ref_hash << "\t" << alt_hash << "\t" << pdict.double_positions.size() << endl;
		int hap1 = 0;
		ptrdiff_t i = get_index_var( pdict.double_positions, pos1 );
		bool in_scaffold = false;
		if (i < pdict.double_positions.size() ) { hap1 = pdict.haplotype[i]; in_scaffold = true; }
                for (auto it2 : vgraph[ref_hash].connections) { /////// reference base
                        ptrdiff_t j = get_index_var( pdict.double_positions, vgraph[it2.first].pos );
                        //cout << "ref: " << ref_hash << "\t" << pos1 << "\t" << hap1 << "\t" << it2.first << "\t" << vgraph[it2.first].pos << "\t" << j << endl;
                        if (j < pdict.double_positions.size() && i != j) {
                                std::string temp_string = split_string_first(it2.first,"_",1);
                                std::string ref_base = split_string_first(temp_string,"_",0);
                                std::string cnx_base = split_string_first(temp_string,"_",1);
                                int is_variant = -1; if (ref_base != cnx_base) { is_variant = 1; }
                                int pos2 = pdict.double_positions[j];
                                int hap2 = pdict.haplotype[j];
                                int call_allele = hap2*is_variant;
				int diffp = std::abs( pos1 - pos2 );
				if ( diffp < max_delta_p ) {
                                	if ( ab_hap[call_allele] == "A" ) { ref_numA += 1; }
                               		if ( ab_hap[call_allele] == "B" ) { ref_numB += 1; }
				}  
				//cout << pos2 << "\t" << hap2 << "\t" << is_variant << "\t" << call_allele << "\t" << ab_hap[call_allele] << endl;
                        }
                }
                for (auto it2 : vgraph[alt_hash].connections) { /////// alternate base
                        ptrdiff_t j = get_index_var( pdict.double_positions, vgraph[it2.first].pos );
                        //cout << "alt: " << alt_hash << "\t" << pos1 << "\t" << hap1 << "\t" << it2.first << "\t" << vgraph[it2.first].pos << "\t" << j << endl;
                        if (j < pdict.double_positions.size() && i != j) {
                                std::string temp_string = split_string_first(it2.first,"_",1);
                                std::string ref_base = split_string_first(temp_string,"_",0);
                                std::string cnx_base = split_string_first(temp_string,"_",1);
                                int is_variant = -1; if (ref_base != cnx_base) { is_variant = 1; }
                                int pos2 = pdict.double_positions[j];
                                int hap2 = pdict.haplotype[j];
                                int call_allele = hap2*is_variant;
				int diffp = std::abs( pos1 - pos2 );
				if ( diffp < max_delta_p ) {
                                	if ( ab_hap[call_allele] == "A" ) { alt_numA += 1; }
                                	if ( ab_hap[call_allele] == "B" ) { alt_numB += 1; }
				}  
				//cout << pos2 << "\t" << hap2 << "\t" << is_variant << "\t" << call_allele << "\t" << ab_hap[call_allele] << endl;
                        }
                }
                ofile << pos1 << "\t" << hap1 << "\t" << ref_numA << "\t" << ref_numB << "\t" << alt_numA << "\t" << alt_numB << endl;
		cin.get();
        }
	ofile.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_recover( int argc, char** argv ) {
        parse_recover_options(argc, argv);
        std::string scaffoldsolutionFile = "./output/hap_recovered_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        cout << "- output recovered file: " << scaffoldsolutionFile << endl;
        cout << endl;
        // ############################### load hap solution and hic links ###########################
        coord_dictionary pdict;
	cout << "reading scaffold file: " << opt::input_scaffold_file << endl;
        read_scaffold( opt::input_scaffold_file, pdict );
        //#############################################
        read_graph rgraph;
        variant_graph vgraph;
	cout << "reading graph file: " << opt::input_graph_file << endl;
        read_variant_graph_file( opt::input_graph_file, opt::chr_choice, vgraph, rgraph );
	cout << "linking hashes" << endl;
        link_hashes( vgraph, rgraph );
	prune_graph( vgraph );
	for ( auto& it : vgraph ) { it.second.unique_hash(); }
        coord_dictionary rdict;
	bool paired;
	initialize_pdict( vgraph, rdict, paired );
	cout << "making recovered dict" << endl;
	recover_links( opt::chr_choice, pdict, vgraph, opt::technology, rdict, scaffoldsolutionFile );
        //#############################################
        return;
};



        //for (int l = 0; l < rdict.num_paired; l++) {
        //      int position = rdict.double_positions[l];
        //      cout << l << "\t" << rdict.num_paired << "\t"<< position << endl;
        //        ref_position_map[position] = rdict.ref_handle[l];
        //        alt_position_map[position] = rdict.alt_handle[l];
        //}


        //for (auto& it : vgraph) {
        //      ref_position_map[it.second.pos] = it.second.ref_base;
        //      alt_position_map[it.second.pos] = it.second.var_base;
        //}

        //coord_dictionary pdict;


        //
        ///////////// write haplotype output
        //cout << " writing scaffold file: " << scaffoldsolutionFile << endl;
        //write_scaffold( opt::chr_choice, scaffoldsolutionFile, bl_dict, pdict, centromere_pos );



        //for (int i = 0; i < pdict.num_paired; i++) {
                //int pos1 = pdict.double_positions[i];
                //int hap1 = pdict.haplotype[i];
                //int ref_hap = hap1;
                //int alt_hap = hap1*-1;
                //std::string ref_hash = pdict.ref_handle[i];
                //std::string alt_hash = pdict.alt_handle[i];


