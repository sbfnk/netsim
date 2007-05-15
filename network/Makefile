#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lboost_program_options -lgsl -lgslcblas

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

#SOURCES     = simulate.cc Bin.cc InfoSIRS.cc Model.cc
SOURCES     = simulate.cc Bin.cc InfoSIRS.cc ProtectiveSIRS.cc Model.cc
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

out/simulate.o: include/Vertex.hh include/Simulator.hh \
	include/GillespieSimulator.hh \
	include/ChrisSimulator.hh include/Tree.hh include/Bin.hh \
	include/InfoSIRS.hh include/Model.hh include/lattice_generator.hh \
	include/graph_structure.hh include/erdos_renyi_generator2.hh \
	include/visualize_graph.hh include/albert_barabasi_generator.hh \
	include/cluster_coeffs.hh
