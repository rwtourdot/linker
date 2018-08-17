#include "run_bx_bin.h"

///////////////////////////////////////////////
namespace opt {
        static std::string chr_choice = "chr20";
        static std::string technology = "tenx";
        static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string id_string = "default";
	static int binsize = 10000;
};

static const char* shortopts = "ho:i:c:e:s:n:";
static const struct option longopts[] = {
        { "bam-file",    no_argument, NULL, 'i' },
        { "technology",  no_argument, NULL, 'e' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "id_string",   no_argument, NULL, 'n' },
	{ "binsize",     no_argument, NULL, 'b' }
};

///////////////////////////////////////////////
static const char *COVERAGE_USAGE_MESSAGE =
"Usage: linker bx_bin [OPTION] -i /path/to/bamfile.bam \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -e,      long read tech (tenx,pacbio,illumina) \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
"  -b,      binsize - raw base number (default 10000 - 10kb) \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_bx_bin_options( int argc, char** argv ) {
        bool die = false;
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'e': arg >> opt::technology; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'n': arg >> opt::id_string; break;
		case 'b': arg >> opt::binsize; break;
                }
        }
        if (die) {
          std::cerr << "\n" << COVERAGE_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running linker het_coverage ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== technology === " << opt::technology << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "== binsize    === " << opt::binsize << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void run_bin_bxtag(int argc, char** argv) {
        parse_bx_bin_options(argc, argv);
        std::string output_file = "./output/bx_unique_binsize_" + std::to_string(opt::binsize) + "_" + opt::id_string + "_" + opt::chr_choice + ".dat";
	cout << " - output bx file: " << output_file << endl;
        //#################### some global structures ################
        //std::map<std::string,int> chr_str_map;   //coord_dictionary pdict;
	contig_dict cdict;
        contig_bx cbx;
        bx_map gen_map;
        full_map chr_map;
	contig_name_length(opt::input_bam_file,cdict);
        //contig_name_map(opt::input_bam_file,chr_str_map);
        //int chromosome = chr_str_map[opt::chr_choice];
        //#################### start of code ##########################
        //if (opt::technology == "pacbio" || opt::technology == "tenx" || opt::technology == "illumina") {
        //        connect_up_variants_bam_pileup(vvec,opt::input_bam_file,chromosome,vgraph,rgraph,opt::technology);
        //        if(opt::multiple_bams) { connect_up_variants_bam_pileup(vvec,opt::second_input_bam_file,chromosome,vgraph,rgraph,opt::technology); }
        //}
        if (opt::technology == "tenx") {
		create_contig_dict_nolinks(opt::input_bam_file,cdict,cbx,opt::chr_choice,opt::binsize,output_file,gen_map,chr_map);
        }
	else { cout << "error: not a valid technology choice [tenx] " << endl; return; }
        //else { write_het_coverage(vgraph,coverageFile,opt::chr_choice); }
	write_bin_bx_cov( gen_map, output_file, opt::chr_choice, opt::binsize );
        return;
};



        //vcf_vector vvec;
        //read_graph rgraph;
        //variant_graph vgraph;



		//connect_up_variants_bam_pileup(vvec,opt::input_bam_file,chromosome,vgraph,rgraph,opt::technology);

        //        calc_coverage_unique_hash(vgraph);
        //        write_het_bx_coverage(vgraph,coverageFile,opt::chr_choice);



//        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
//        static bool multiple_bams = false;
//        static std::string second_input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";



//        { "vcf-file",    no_argument, NULL, 'v' },
//        { "bam-file2",   no_argument, NULL, 's' },
//        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
