COMPFLAGS     =

CXXFLAGS      = $(COMPFLAGS) $(FlameManCXXOpts) -I$(FlameManLib) -D$(FlameManMach) -D$(FlameManCVODEVERSION) -I$(FlameManCVODEDIR)/include

CXX           = $(FlameManCXX)


HDRS	      = BroadeningFactors.h \
		FlameMaster.h \
		TCountDiffFlamePhys.h \
		TCountDiffPhysEigen.h \
		TCountDiffFlameMix.h \
		TCountDiffFlameCont.h \
		TCountDiffFlameSim.h \
		TFlameSheet.h \
		T0DIsoChor.h \
		TTransFlamelet.h \
		TTransFlameSolver.h \
		TCountPremFlameSim.h \
		TCountPremFlamePhys.h \
		TUnstrPremFlamePhys.h \
		ScanManStructs.h

#		TTrans1DIsoChor.h \
#		TTrans1DIsoChorSolver.h \
#		TTransCoalParticle.h \
#		TTransCoalParticleSolver.h \

LD	      = $(CXX)

LDFLAGS	      = $(COMPFLAGS) $(FlameManLDCCOpts) $(FlameManLDFortOpts)

LIBS          = $(FlameManLib)/NewtonLib.a \
	$(FlameManLib)/FMLib.a \
	$(FlameManLib)/libAM.a \
	$(FlameManLib)/alligator.a \
	$(FlameManCVODEDIR)/lib/libsundials_cvode.a \
	$(FlameManCVODEDIR)/lib/libsundials_nvecserial.a \
	$(FlameManLIBF77)

#	$(FlameManLib)/ddasslLib.a \


MAKEFILE      = Makefile

OBJS	      = BroadeningFactors.o \
		FlameMasterMain.o \
		FlameUtilities.o \
		TFlame.o \
		TInputData.o \
		TProperties.o \
		TReaction.o \
		SteadyStates.o \
		TSpecies.o \
		TSoot.o \
		TPAH.o \
		TCountDiffFlameMix.o \
		TCountDiffFlameCont.o \
		TCountDiffFlameSim.o \
		TCountDiffFlamePhys.o \
		TCountDiffPhysEigen.o \
		TFlameSheet.o \
		T0DIsoChor.o \
		TTransFlamelet.o \
		TTransFlameSolver.o \
		TCountPremFlameSim.o \
		TCountPremFlamePhys.o \
		TUnstrPremFlamePhys.o \
		lex.yy.o

#		TTrans1DIsoChor.o \
#		TTrans1DIsoChorSolver.o \
#		TTransCoalParticle.o \
#		TTransCoalParticleSolver.o \

PROGRAM       = FlameMaster

SRCS	      = BroadeningFactors.C \
		FlameMasterMain.C \
		FlameUtilities.C \
		TFlame.C \
		TInputData.C \
		TProperties.C \
		TReaction.C \
		SteadyStates.C \
		TSpecies.C \
		TSoot.C \
		TPAH.C \
		TCountDiffFlameMix.C \
		TCountDiffFlameCont.C \
		TCountDiffFlameSim.C \
		TCountDiffFlamePhys.C \
		TCountDiffPhysEigen.C \
		TFlameSheet.C \
		T0DIsoChor.C \
		TTransFlamelet.C \
		TTransFlameSolver.C \
		TCountPremFlameSim.C \
		TCountPremFlamePhys.C \
		TUnstrPremFlamePhys.C \
		lex.yy.C

#		TTrans1DIsoChor.C \
#		TTrans1DIsoChorSolver.C \
#		TTransCoalParticle.C \
#		TTransCoalParticleSolver.C \

$(PROGRAM):     $(OBJS) $(LIBS) $(MAKEFILE)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(OBJS) $(LIBS) $(LDFLAGS) -o $(PROGRAM)
		cp $(PROGRAM) $(FlameManLib)
		@echo "done"

touch:;		@touch $(SRCS) $(HDRS);rm lex.yy.C
toucho:;	@touch lex.yy.C
		@sleep 1
		@touch $(OBJS)

lex.yy.C: ReadAddData.flex $(MAKEFILE) FlameMaster.h
			$(FLEXPATH)flex -8 -t -i -s ReadAddData.flex > lex.yy.C
BroadeningFactors.o: BroadeningFactors.h
FlameMasterMain.o: $(HDRS)
#FlameMasterMain.o: FlameMaster.h
TSoot.o: FlameMaster.h ScanManStructs.h
TPAH.o: FlameMaster.h ScanManStructs.h
FlameUtilities.o: FlameMaster.h ScanManStructs.h
TCountDiffFlamePhys.o: FlameMaster.h TCountDiffFlamePhys.h ScanManStructs.h
TCountDiffPhysEigen.o: FlameMaster.h TCountDiffPhysEigen.h ScanManStructs.h
TCountDiffFlameMix.o: FlameMaster.h TCountDiffFlameMix.h ScanManStructs.h
TCountDiffFlameCont.o: FlameMaster.h TCountDiffFlameCont.h TCountDiffFlameSim.h ScanManStructs.h
TCountDiffFlameSim.o: FlameMaster.h TCountDiffFlameSim.h ScanManStructs.h
TFlameSheet.o: FlameMaster.h TFlameSheet.h ScanManStructs.h 
T0DIsoChor.o: FlameMaster.h T0DIsoChor.h ScanManStructs.h 
TTransFlamelet.o: FlameMaster.h TTransFlamelet.h TTransFlameSolver.h ScanManStructs.h 
TTransFlameSolver.o: FlameMaster.h TTransFlameSolver.h ScanManStructs.h 
TCountPremFlameSim.o: FlameMaster.h TCountPremFlameSim.h ScanManStructs.h 
TCountPremFlamePhys.o: FlameMaster.h TCountPremFlamePhys.h ScanManStructs.h 
TUnstrPremFlamePhys.o: FlameMaster.h TUnstrPremFlamePhys.h ScanManStructs.h 
TFlame.o: FlameMaster.h ScanManStructs.h 
TInputData.o: FlameMaster.h ScanManStructs.h 
TProperties.o: FlameMaster.h ScanManStructs.h
TReaction.o: FlameMaster.h BroadeningFactors.h ScanManStructs.h
SteadyStates.o: FlameMaster.h
TSpecies.o: FlameMaster.h ScanManStructs.h 
#TTrans1DIsoChor.o: FlameMaster.h TTrans1DIsoChor.h TTrans1DIsoChorSolver.h ScanManStructs.h 
#TTrans1DIsoChorSolver.o: FlameMaster.h TTrans1DIsoChorSolver.h ScanManStructs.h 
#TTransCoalParticle.o: FlameMaster.h TTransCoalParticle.h TTransCoalParticleSolver.h ScanManStructs.h 
#TTransCoalParticleSolver.o: FlameMaster.h TTransCoalParticleSolver.h ScanManStructs.h 
lex.yy.o: FlameMaster.h ScanManStructs.h
$(OBJS) : $(MAKEFILE)

.SUFFIXES: $(.SUFFIXES) .cp .C

.cp.C:
	mv $< $*.C
.C.o:
	$(CXX) $(CXXFLAGS) -c $<
.cp.o:
	mv $< $*.C;$(CXX) $(CXXFLAGS) -c $*.C

