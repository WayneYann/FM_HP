### change the location of the FlameMaster package here
setenv FlameManPath ~/FlameMaster

### do not change the following couple of lines
setenv FlameManSource $FlameManPath/FlameManSource
setenv FlameManLib $FlameManPath/FlameManLibs
setenv myData $FlameManPath/FlameManData
alias FlameMaster $FlameManLib/FlameMaster
alias LT '$FlameManLib/ListTool -M'
set path=($path $FlameManLib)

### You have to choose between two possible options for FLEXPATH depending
### on weather flex is installed on your machine and in the current path or not.
### You can find out by typing 'which flex' in the command line. If you get a valid
### path as a result, this means flex is installed. Then choose the first option.
### If flex cannot be found, choose the second one. Then, flex will be installed
### during the FlameMaster installation. You will do this following the procedure
### outlined in the Readme file.
###
### Option 1: flex already exists on your machine (don't insert the path in the FLEXPATH definition)
setenv FLEXPATH ""
### Option 2: flex will be installed to $FlameManLib
#setenv FLEXPATH $FlameManLib/bin/

# CVODE or the entire sundials library has to be installed. The path has to be provided below


### choose one of the following sets

# following for Linux - GNU compiler
setenv FlameManMach LINUXGXX
setenv FlameManCC gcc
setenv FlameManCCOpts "-ansi -O3"
setenv FlameManCXX g++
setenv FlameManCXXOpts "-x c++ -O3"
setenv FlameManYACC yacc
setenv FlameManLDCCOpts ""
setenv FlameManRANLIB 'ranlib $(LIBRARY)'
setenv FlameManCVODEDIR "/opt/sundials/intel"
setenv FlameManCVODEVERSION "SUNDIALS24"

# following for Mac OSX - INTEL compiler
#setenv FlameManMach DIGITAL
#setenv FlameManCC icc
#setenv FlameManCCOpts "-O3"
#setenv FlameManCXX icc
#setenv FlameManCXXOpts "-D MACOSX"
#setenv FlameManYACC yacc
#setenv FlameManLDCCOpts "-lm -lstdc++"
#setenv FlameManCVODEDIR "/opt/sundials/intel"
#setenv FlameManCVODEVERSION "SUNDIALS24"

# following for Mac OSX - GNU compiler
#setenv FlameManMach DIGITAL
#setenv FlameManCC cc
#setenv FlameManCCOpts "-O3"
#setenv FlameManCXX CC
#setenv FlameManCXXOpts "-O3"
#setenv FlameManYACC yacc
#setenv FlameManLDCCOpts "-lm"
#setenv FlameManRANLIB 'ranlib $(LIBRARY)'
#setenv FlameManCVODEDIR "/opt/sundials/intel"
#setenv FlameManCVODEVERSION "SUNDIALS24"
