f90compiler = gfortran
optimize = -O2
# debug = -g
LD = gfortran
LDFLAGS = $(optimize)
LIB_HDF4 = -L$(HOME)/Documents/local/hdf4/lib -lm -lmfhdf -ldf -lz -lsz -ljpeg

######## The targets ########

equiSources: definitionsModule.o utilities.o dustModule.o stellarPopulationModule.o \
	stellarBetaTable.o uniformTable.o uvbBetaTable.o coll_rates.o colh2diss.o calc_rates.o \
	rotateIndicesModule.o transportRoutinesModule.o equiSources.o
	$(LD) $(LDFLAGS) $^ $(LIB_HDF4) -o $@
bin2hdf4: definitionsModule.o rotateIndicesModule.o bin2hdf4.o
	$(LD) $(LDFLAGS) $^ $(LIB_HDF4) -o $@
hdf42bin: hdf42bin.o
	$(LD) $(LDFLAGS) $^ $(LIB_HDF4) -o $@
readCellArray: definitionsModule.o rotateIndicesModule.o readCellArray.o
	$(LD) $(LDFLAGS) $^ $(LIB_HDF4) -o $@
convertFormats: convertFormats.o
	$(LD) $(LDFLAGS) $^ $(LIB_HDF4) -o $@

######## Rules to create object files ########

.SUFFIXES: .o .f90 .f77 .f
.f90.o:
	$(f90compiler) $(optimize) $(debug) -c $< -o $@
.f.o:
	$(f90compiler) $(optimize) $(debug) -c $< -o $@

######## Special targets ########

clean:
	@/bin/rm -rf *.o *.mod *~ core.* equiSources bin2hdf4 hdf42bin readCellArray convertFormats
