#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lboost_program_options -lgsl -lgslcblas

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = simulate.cc Bin.cc InfoSIRS.cc Model.cc
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
	$(MAKEDEPEND) $(DEPENDFLAGS) *.cc

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f $(OBJECTS) $(EXEC)
	@echo ... done

cleaner:
	@echo '... cleaning also *~'
	@rm -f $(OBJECTS) *~
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

InfoSIRS.o: InfoSIRS.hh Model.hh

main.o: Vertex.hh Simulator.hh GillespieSimulator.hh ChrisSimulator.hh \
        Tree.hh Bin.hh InfoSIRS.hh Model.hh lattice_generator.hh \
        graph_structure.hh erdos_renyi_generator2.hh visualize_graph.hh
