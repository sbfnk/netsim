#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lboost_program_options-gcc41 -lpng -lpngwriter -lfreetype

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = network.cc
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = network
MANIP       = pngmanip
#--------------------------------------------------------------------

all:  bin/$(EXEC) bin/$(MANIP)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

bin/$(MANIP): out/pngmanip.o
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $< -o $@

out/%.o: src/%.cc include/%.hh
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) src/*.cc

htmldoc:
	doxygen

pdfdoc: htmldoc
	$(MAKE) -C doc/latex pdf

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f bin/* out/*
	@echo ... done

cleaner:
	@echo '... cleaning also *~'
	@rm -f bin/* out/* *~ .#*
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

out/network.o: include/lattice_generator.hh include/graph_structure.hh \
	include/erdos_renyi_generator2.hh include/graph_io.hh \
	include/albert_barabasi_generator.hh include/cluster_coeffs.hh \
	include/assortativity.hh include/graph_statistics.hh \
	include/path_length.hh include/breadth_first_search.hpp \
	include/degree_overlap.hh include/tree_generator.hh \
        include/community_generator.hh include/community_structure.hh
