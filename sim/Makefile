#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -I/home/seb/.local/include/ -I.. -Iinclude -O3 -g 

LINKER      = g++
LDFLAGS     = -L/home/seb/.local/lib -lpngwriter -lpng -lfreetype \
              -lboost_system -lboost_program_options -lboost_filesystem \
              -lboost_iostreams  -g

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = simulate.cc Bin.cc Vertex.cc Model.cc StatRecorder.cc
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = simulate
AVG         = average_runs
#--------------------------------------------------------------------

all:  bin/$(EXEC) bin/$(AVG)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo

bin/$(AVG): out/average_runs.o
	$(LINKER) $(OPTFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

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

out/simulate.o: include/Vertex.hh include/Simulator.hh \
	include/EpiSimulator.hh include/GillespieSimulator.hh \
	include/GroupFormSimulator.hh include/Model.hh \
	include/GroupFormModel.hh include/EpiRewireSimulator.hh \
	include/ChrisSimulator.hh include/Tree.hh include/Bin.hh \
	include/InfoSIRS.hh include/DimInfoSIRS.hh include/EpiModel.hh \
	include/ProtectiveSIRS.hh include/VaccinationSIRS.hh \
	include/SingleSIRS.hh include/SIRInfo.hh \
        include/sim_statistics.hh include/RewireSimulator.hh \
        include/GroupFormModel.hh include/NetEvoSimulator.hh \
	include/EvoModel.hh \
	../network/include/graph_structure.hh \
	../network/include/graph_io.hh \
	../network/include/nearest_infected.hh \
	../network/include/Edge.hh \
	../network/include/community_structure.hh \
	../network/include/graph_structure.hh

