#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lpng -lpngwriter -lfreetype -lboost_program_options 

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

#SOURCES     = simulate.cc Bin.cc InfoSIRS.cc Model.cc
SOURCES     = simulate.cc Bin.cc InfoSIRS.cc ProtectiveSIRS.cc \
	VaccinationSIRS.cc Model.cc
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = simulate
AVG         = average_runs
#--------------------------------------------------------------------

all:  bin/$(EXEC) bin/$(AVG)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

bin/$(AVG): out/average_runs.o
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $< -o $@

out/%.o : src/%.cc include/%.hh
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

out/InfoSIRS.o: include/InfoSIRS.hh include/Model.hh

out/ProtectiveSIRS.o: include/ProtectiveSIRS.hh include/Model.hh

out/Vaccination.o: include/Vaccination.hh include/Model.hh

out/simulate.o: include/Vertex.hh include/Simulator.hh \
	include/GillespieSimulator.hh \
	include/ChrisSimulator.hh include/Tree.hh include/Bin.hh \
	include/InfoSIRS.hh include/Model.hh include/lattice_generator.hh \
	include/graph_structure.hh include/erdos_renyi_generator2.hh \
	include/visualise_graph.hh include/albert_barabasi_generator.hh \
	include/cluster_coeffs.hh include/assortativity.hh \
	include/ProtectiveSIRS.hh include/VaccinationSIRS.hh \
	include/graph_statistics.hh include/path_length.hh \
	include/breadth_first_search.hpp include/degree_overlap.hh
