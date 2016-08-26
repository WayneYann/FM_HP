#!/bin/csh 
 
echo "Install FlameMan to "$FlameManLib
echo ""

echo "### Installation instructions are given in the Readme file."
echo "### "
echo "### This installation script always ensures a full build of all"
echo "### libraries and programs. This can be helpful in case you make"
echo "### mistakes or have to experiment with the settings in the file"
echo "### ~/.FlameManExports"
echo "### If you are sure that the installation of a particular part was "
echo "### correct, full build can be prevented by commenting the respective"
echo "### 'rm *.o' or 'make touch' statements."

#install bison only if yacc does not work
#cd bison-1.25
#echo ""
#echo "start installing bison"
#mkdir Installation
#configure --prefix=$FlameManLib
#make
#chmod 755 install-sh
#make install
#cd ..

cd $FlameManSource/..
cd FlameManLibsSource

#cd DasslLib
#echo ""
#echo "start installing DasslLib"
#rm *.o
#make
#cd ..

cd Alligator
echo ""
echo "start installing Alligator"
rm *.o
make
cd ..
cd ArrayMan
echo ""
echo "start installing ArrayMan"
rm *.o
make
cd ..
cd FMLibs
echo ""
echo "start installing FMLibs"
make touch
make
cd ..
cd NewtonLib
echo ""
echo "start installing NewtonLib"
make touch
make
cd ../..
cd FlameManSource
cp SteadyStateDefault.cp SteadyStates.C
echo ""
echo "start installing FlameMan"
make touch
make
cd ..
cd FlameManTools
cd CreateBinFile
echo  ""
echo "start installing CreateBinFile"
make touch
make
chmod 755 ScanIt.script
./ScanIt.script >&! Scan.out
cd ..
echo ""
echo "start installing ListTool package"
cd ListTool
cd WSS
echo ""
echo "start installing WSS"
rm *.o
make 
cd ..
cd List
echo ""
echo "start installing List"
rm *.o
make
cd ..
cd FLReader
echo ""
echo "start installing FLReader"
rm *.o
make
cd ..
cd ListTool
echo ""
echo "start installing ListTool"
rm *.o
make
cd ../../..
cd ScanManSource
echo ""
echo "start installing ScanMan"
make touch
make
cd ..
cd ScanManRun
chmod 755 ScanMan.script
./ScanMan.script H2.allstar3
./ScanMan.script CH4.72
./ScanMan.script C3H8.82
./ScanMan.script CH4.Igni73
cd HeptaneShort
$FlameManLib/CreateBinFile -i nHeptane.allstarnew_oks.thermo -m nHeptane.allstarnew_oks.trans -o nHeptane.allstarnew_oks.thermo.bin
$FlameManLib/ScanMan -i nHeptane.allstarnew_oks.mech -t nHeptane.allstarnew_oks.thermo.bin -S
mv nHeptane.allstarnew_oks.pre $myData
cd ../..

echo ""
echo done
