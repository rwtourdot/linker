//####################################
#include "linker.h"

static const char *LINKER_USAGE_MESSAGE =
"############### running linker ############### \n"
"Program: Linker \n"
"Version: 2.1 \n"
"Usage: linker <command> [options]\n\n"
"Commands:\n"
"	phase		phase germline haplotype from linked read sample \n"
"	cn_phase	further phase genome due to LOH or aneuplody in genome \n"
"	sv_phase	phase sv call on long reads \n"
"	hic_phase	phase split reads in hic \n"
"	coverage	extract coverage of each het site from bam file \n"
"	matrix		bin genome and create linked read heatmap \n"
" 	filter		filter normal variants by copy fraction \n"
"	bx_bin          get unique bx coverage from 10X \n"
"\n";


int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cerr << LINKER_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << LINKER_USAGE_MESSAGE;
      return 0;
    } else if (command == "phase") { 
	run_link_phaser(argc-1, argv+1);
    } else if (command == "cn_phase") { 
	run_cn_phaser(argc-1, argv+1);
    } else if (command == "coverage") { 
	run_het_coverage(argc-1, argv+1);
    } else if (command == "matrix") {  	
	run_bin_matrix(argc-1, argv+1);
    } else if (command == "sv_phase") { 
	run_sv_phaser(argc-1, argv+1); 	
    } else if (command == "filter") {
        run_variant_filtering(argc-1, argv+1);
    } else if (command == "bx_bin") {
        run_bin_bxtag(argc-1, argv+1);
    } else if (command == "hic_phase") {
        run_hic_phaser(argc-1, argv+1);
    }
    else {
      std::cerr << LINKER_USAGE_MESSAGE;
      return 0;
    }
  }

  std::cerr << "done with linker" << std::endl;
  return 0;
}
