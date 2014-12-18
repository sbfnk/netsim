#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude -I/home/seb/.local/include

LINKER      = g++
LDFLAGS     = -L/home/seb/.local/lib -lboost_program_options -lboost_iostreams \
              -lpngwriter -lpng -lfreetype

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
	$(LINKER) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo

bin/$(MANIP): out/pngmanip.o
	$(LINKER) $(OPTFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

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
        include/community_generator.hh include/nearest_infected.hh \
        include/community_structure.hh
