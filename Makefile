PROG=linker
CC=c++ -o3 -std=c++11                           # linux c++
#CC=clang++ -std=c++11 -stdlib=libc++ -lz       # macos c++
#CC=c++ -g -Wall -Wextra -std=c++11   #g++

#LDIR=/usr/lib/x86_64-linux-gnu/
LDIR=./packages/bamtools/lib/
IDIR2=./htslib/
LDIR2=./packages/htslib/lib/
IDIR=./packages/bamtools/include/bamtools/
IDIR3=./packages/bamtools/lib/

CFLAGS=-I$(IDIR) -I$(IDIR2)
LFLAGS=-L$(LDIR) -L$(LDIR2) -Wl,-rpath,$(IDIR2),-rpath,$(IDIR3) #-L$(LDIR3)
LBAMTOOLS=-lbamtools
LHTS=-lhts
#LCURL=-lcurl
LIBRARIES=$(LHTS) $(LBAMTOOLS) #$(LHTS) -fPIC
#CCPLUS=$(CC) $(LFLAGS) $(CFLAGS)

SRCDIR:=src
#BINDIR:=bin
#TARGET:=./$(BINDIR)/$(PROG)
TARGET:=./$(PROG)
SOURCES=$(shell find ./$(SRCDIR)/*.cpp)
OBJECTS=$(addsuffix .o,$(basename $(SOURCES)) )

##########################
all: $(TARGET)

##########################
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LFLAGS) $(CFLAGS) $(LIBRARIES)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(LFLAGS) $(CFLAGS) -c -o $@ $<

print-% : ; @echo $* = $($*)

##########################
clean:
	-@rm $(OBJECTS)
	-@rm -v $(TARGET)

distclean:
	-@rm $(OBJECTS)
	-@rm -v $(TARGET)
	-@rm ./output/variant_network_*
	-@rm ./output/hap_solution_*
	-@rm ./output/het_coverage_*
	-@rm ./output/bxmatrix_*
	-@rm ./output/filtered_coverage_*
	-@rm ./output/sv_phased_*
	-@rm ./output/cn_phased_*
	-@rm ./output/hic_links_*
	-@rm ./output/phased_sv_sites.dat

##########################
batch:
	make clean
	make all

























##########################
#$(PROG).o: $(PROG).cpp $(PROG).h
#	$(CCPLUS) -c $(PROG).cpp $(LIBRARIES) -o $(PROG).o

##########################
#run_bin_matrix.o: run_bin_matrix.cpp run_bin_matrix.h bin_reference.h
#	$(CCPLUS) -c run_bin_matrix.cpp $(LIBRARIES) -o run_bin_matrix.o

##########################
#run_cn_phaser.o: run_cn_phaser.cpp run_cn_phaser.h
#	$(CCPLUS) -c run_cn_phaser.cpp $(LIBRARIES) -o run_cn_phaser.o

##########################
#run_het_coverage.o: run_het_coverage.cpp run_het_coverage.h read_tree.h variant_site.h
#	$(CCPLUS) -c run_het_coverage.cpp $(LIBRARIES) -o run_het_coverage.o

##########################
#run_link_phaser.o: run_link_phaser.cpp run_link_phaser.h map_matrix.h sub_matrix.h coord_dict.h variant_site.h
#	$(CCPLUS) -c run_link_phaser.cpp $(LIBRARIES) -o run_link_phaser.o

##########################
#write_linker_output.o: write_linker_output.cpp write_linker_output.h coord_dict.h variant_site.h
#	$(CCPLUS) -c write_linker_output.cpp $(LIBRARIES) -o write_linker_output.o

##########################
#read_vcf.o: read_vcf.cpp read_vcf.h variant_site.h coord_dict.h variant_site.h
#	$(CCPLUS) -c read_vcf.cpp $(LIBRARIES) -o read_vcf.o

##########################
#read_bam.o: read_bam.cpp read_bam.h variant_site.h read_tree.h bin_reference.h
#	$(CCPLUS) -c read_bam.cpp $(LIBRARIES) -o read_bam.o

##########################
#mc_solver.o: mc_solver.cpp mc_solver.h coord_dict.h map_matrix.h sub_matrix.h
#	$(CCPLUS) -c mc_solver.cpp $(LIBRARIES) -o mc_solver.o

##########################
#build_hap_matrix.o: build_hap_matrix.cpp build_hap_matrix.h coord_dict.h map_matrix.h variant_site.h read_tree.h
#	$(CCPLUS) -c build_hap_matrix.cpp $(LIBRARIES) -o build_hap_matrix.o

##########################
#bin_assembly.o: bin_assembly.cpp bin_assembly.h bin_reference.h read_tree.h
#	$(CCPLUS) -c bin_assembly.cpp $(LIBRARIES) -o bin_assembly.o






#$(PROG): $(PROG).o bin_reference.o coord_dict.o map_matrix.o read_tree.o sub_matrix.o variant_site.o read_vcf.o read_bam.o bin_assembly.o write_linker_output.o mc_solver.o build_hap_matrix.o run_bin_matrix.o run_cn_phaser.o run_het_coverage.o run_link_phaser.o
#	$(CCPLUS) $(PROG).o bin_reference.o coord_dict.o map_matrix.o read_tree.o sub_matrix.o variant_site.o read_vcf.o read_bam.o bin_assembly.o write_linker_output.o mc_solver.o build_hap_matrix.o run_bin_matrix.o run_cn_phaser.o run_het_coverage.o run_link_phaser.o -o $(PROG) $(LIBRARIES)

#%.o: %.cpp $(DEPS)
#	$(CCPLUS) -c $*.cpp -o $*.o $(LIBRARIES)


#	$(CC) $(LFLAGS) $(CFLAGS) $(PROG).cpp -o $(PROG) -lbamtools -lhts


	#%.o: %.c
	#	$(CC) $(LFLAGS) $(CFLAGS) -o $@ -c $<

	#all: linker.cpp $(DEPS)     #utils.cpp  hash_link.cpp
	#	$(CC) $(LFLAGS) $(CFLAGS) $(PROG).cpp -o $(PROG) -lbamtools -lhts


#DEPS= read_tree.h map_matrix.h bin_reference.h coord_dict.h sub_matrix.h variant_site.h

# Set the LD_LIBRARY_PATH [*]
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. #/usr/lib/x86_64-linux-gnu/
#export LD_LIBRARY_PATH=.

# gaaaah -lbamtools needs to be after the .cpp file(s)

#@rm -fv ./*.out
#$(CC) $(LFLAGS) $(CFLAGS) $(PROG).cpp hash_link.cpp -o $(PROG) -lbamtools -lhts
