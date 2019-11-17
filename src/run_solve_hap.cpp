#include "run_solve_hap.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr20";
	static std::string input_graph_file = "./output/graph_variant_jun10_BL1954_tenx_chr20.dat";
        static std::string id_string = "default";
	static int start_bound = 0;
	static int end_bound = 300000000;
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:i:c:n:";
static const struct option longopts[] = {
	{ "help",        no_argument, NULL, 'h' },
        { "graph-file",    no_argument, NULL, 'i' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *PHASER_USAGE_MESSAGE =
"Usage: linker solve [OPTION] -i /output/graph_variant_*.dat -c chr20 -n trial \n\n"
"\n"
"  Options\n"
"  -i,      input graph_variant file \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_solve_hap_options( int argc, char** argv ) {
	bool die = false; //bool vcf_load = false; //bool cov_load = false;
	//if(argc <= 2) { die = true; }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << PHASER_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
		case 'h': die = true; break;
                case 'i': arg >> opt::input_graph_file; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
	if (die) {
	  std::cerr << "\n" << PHASER_USAGE_MESSAGE;
	  exit(1);
	}
	//parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running link phaser ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== input graph  === " << opt::input_graph_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_solve_hap( int argc, char** argv ) {
        parse_solve_hap_options(argc, argv);
        std::string output_network_file = "./output/variant_network_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        std::string hapsolutionFile = "./output/hap_solution_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        cout << "- output net file: " << output_network_file << endl;
        cout << "- output hap file: " << hapsolutionFile << endl;
        cout << endl;
        coord_dictionary pdict;
        read_graph rgraph;
        variant_graph vgraph;
	bool paired;
        //#################### start of code ##########################
	cout << "loading graph " << endl;
	read_variant_graph_file(opt::input_graph_file,opt::chr_choice,vgraph,rgraph);	
	cout << "linking graph " << endl;
	//#############################################################
        link_hashes(vgraph,rgraph);
        prune_graph(vgraph);
        initialize_pdict(vgraph,pdict,paired);
	//#############################################################
        write_link_network(vgraph,output_network_file,opt::chr_choice);
        static const std::size_t length = pdict.num_paired;
	///////////// create global datastructures for haplotype solver
        map_matrix<int> num_matrix_second(length);
	map_matrix<double> diff_matrix(length);
        initialize_solver(vgraph,pdict,diff_matrix,num_matrix_second);
	solver_recursive(vgraph,pdict,diff_matrix,num_matrix_second);
	//solver(vgraph,pdict,diff_matrix,num_matrix_second);
	///////////// write haplotype output
	cout << " solver finished " << endl;
	int switch_cutoff = -700;
	call_blocks(pdict,switch_cutoff);
	cout << " called haplotype blocks " << endl;
        write_hap_solution(vgraph,hapsolutionFile,pdict,opt::chr_choice);
        return;
};


