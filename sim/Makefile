#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -g

LINKER      = g++
LDFLAGS     = -lboost_program_options

#------ Intel icpc  --------------------------------------------------

#CC         = icpc
#CPPFLAGS     = -g -Wall
#OPTFLAGS   = -O3 -ipo -fnsplit -ip -ftz -funroll-loops
#OPTFLAGS   = -static
#OPTFLAGS   = -fast / -mp / 
#OPTFLAGS   = 

#LINKER      = icpc
#LDFLAGS   = -lm -lscs -lmkl 
#LDFLAGS   = -i-dynamic/static
#LDFLAGS     = 

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = main.cc Bin.cc Model.cc
OBJECTS     = $(SOURCES:.cc=.o)

EXEC        = simulate
AVG         = average_runs
#--------------------------------------------------------------------

all:  $(EXEC) $(AVG)

$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

$(AVG): average_runs.o
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) average_runs.o -o $@

%.o : %.cc
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

main.o: Vertex.hh GillespieSimulator.hh Tree.hh Bin.hh Model.hh \
        lattice_generator.hh graph_structure.hh erdos_renyi_generator2.hh \
        visualize_graph.hh

