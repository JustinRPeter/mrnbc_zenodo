
EXECUTABLE=mrnbc

FC=gfortran
#FCFLAGS=-Wall -c -Wno-tabs -Wextra
FCFLAGS=-c -Wno-tabs
#FCFLAGS=-Wtabs -Wextra
SOURCES=MRNBC_main.f MRNBC_QMM_rev_2.f MRNBC_stat.f
CONFIGS=data.inc dimension.inc temp.inc
OBJECTS=$(SOURCES:.f=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(OBJECTS) -o $@

.f.o: 
	$(FC) $(FCFLAGS) $< 



