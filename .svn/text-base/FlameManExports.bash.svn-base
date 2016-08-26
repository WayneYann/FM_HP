### change the location of the FlameMaster package here
export FlameManPath=~/FlameMaster

### do not change the following couple of lines
export FlameManSource=$FlameManPath/FlameManSource
export FlameManLib=$FlameManPath/FlameManLibs
export myData=$FlameManPath/FlameManData
alias FlameMaster=$FlameManLib/FlameMaster
alias LT='$FlameManLib/ListTool -M'
export PATH=$PATH:$FlameManLib

### You have to choose between two possible options for FLEXPATH depending
### on weather flex is installed on your machine and in the current path or not.
### You can find out by typing 'which flex' in the command line. If you get a valid
### path as a result, this means flex is installed. Then choose the first option.
### If flex cannot be found, choose the second one. Then, flex will be installed
### during the FlameMaster installation. You will do this following the procedure
### outlined in the Readme file.
###
### Option 1: flex already exists on your machine (don't insert the path in the FLEXPATH definition)
export FLEXPATH=""
### Option 2: flex will be installed to $FlameManLib
#export FLEXPATH=$FlameManLib/bin/


### choose one of the following sets

# following for Linux - GNU compiler
export FlameManMach=LINUXGXX
export FlameManCC=gcc
export FlameManCCOpts="-O3"
export FlameManCXX=g++
export FlameManCXXOpts="-O3"
export FlameManYACC="yacc"
export FlameManLDCCOpts=""
export FlameManRANLIB='ranlib $(LIBRARY)'
export FlameManCVODEDIR="/opt/cvode/"
export FlameManCVODEVERSION="SUNDIALS24"

# following for Mac OSX - INTEL compiler
#export FlameManMach=DIGITAL
#export FlameManCC=icc
#export FlameManCCOpts="-O3"
#export FlameManCXX=icc
#export FlameManCXXOpts="-D MACOSX"
#export FlameManYACC="yacc"
#export FlameManLDCCOpts="-lm -lstdc++"
#export FlameManCVODEDIR="/opt/sundials/intel"
#export FlameManCVODEVERSION="SUNDIALS24"

# following for Mac OSX - GNU compiler
#export FlameManMach=DIGITAL
#export FlameManCC=cc
#export FlameManCCOpts="-O3"
#export FlameManCXX=CC
#export FlameManCXXOpts="-O3"
#export FlameManYACC="yacc"
#export FlameManLDCCOpts="-lm"
#export FlameManRANLIB='ranlib $(LIBRARY)'
#export FlameManCVODEDIR="/opt/sundials/intel"
#export FlameManCVODEVERSION="SUNDIALS24"
