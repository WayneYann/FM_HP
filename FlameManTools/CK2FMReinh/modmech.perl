#!/usr/bin/perl -w 
#modmech13.perl Translated from previous C++ program
#modmech14.perl The allowed limit for the common temperature in the
#               Chemkin thermodata file is changed from (500-1500) to
#               (50-5000) 
#modmech15.perl The additional specifier in front of a species name is
#               converted to uppercase e.g.   N-C3H7  (before n-C3H7)
#modmech16.perl -T can also be used to translate species in the FlameMaster 
#               mechanism (not thermodata though)
#modmech17.perl now the temperature can be specified at which the Steady State
#               reactions are evaluated ($tempss) 
#Version 18.0   the help output was updated

$tempss = 1000;
$maxspecies = 1200;
$gotocorrect = 1.0;
$maxuppertempforfit = 2100.0;
$factornumeric = 1000.0;
$program_name = "modmech.pl";
$program_version = "version 18.0 (2003/10/06)";
$maxnewspecieslength = 25;

$fmechnameout = "newmechfile";  # name of the new mechanism file
$fintmechnameout = "intmechfile";
$fsteadystatesnameout = "new.ss";  # name of the new .ss file
$fthermnameout = "newthermofile";   # name of the new thermdat file
$fspecnameout = "speciestranslated";
$increaseafactor=0.5;
# name of the new thermodata file
#int i,j;

sub usage {
    $exit_value = $_[0];
    #print "exitvalue: $exit_value\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "modmech.perl - $program_version\n";
    }
	printf STDERR"usage:\n$program_name [-h] [-T translationtablefile] [-i species_to_include] [-x species_to_exclude] [-t thermodatafile] [-r transportdatafile] [-o mechanismoutputfile] ";
    printf STDERR"mechanism_file_name\n\n\t-h: help. \n\t\tgenerates a help message and exits\n\n\t";
    printf STDERR"-i species_to_include\n\t\t";
    printf STDERR "filename where species are listed or species string\n\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "\t\tReactions where any of these species appear will be used,\n\t\t";
	printf STDERR "unless one of the species of the reaction is also to be excluded.\n\n";
    } 
    printf STDERR "\t-x species_to_exclude\n\t\t";
    printf STDERR "filename where species are listed or species string\n\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "\t\tThe species files contain species to be in/excluded, one species per line.\n\t\t";
	printf STDERR "Alternatively a single species can be given in the command line.\n\t\t";
	printf STDERR "In general species names are case-sensitive.\n\t\t";
	printf STDERR "The wildcard '*' can be used at the beginning or end of an expression\n\t\t";
	printf STDERR "In case of simultaneous translation of names (thermodata file given),\n\t\t";
	printf STDERR "the species to be in/excluded are new names.\n\n"; 
    } 
    printf STDERR "\t-t thermodatafile\n\t\t";
    printf STDERR "filename where species original names are listed.\n\t\t";
    printf STDERR "They will be translated into the origname-chemicalname\n\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "\t\tThe thermodata file can be in CHEMKIN format.\n\t\t";
	printf STDERR "or a file of species containing origname and newname per line.\n\n"; 
    } 
    printf STDERR "\t-r transportdatafile\n\t\t";
    printf STDERR "filename where species with their transportdata are listed.\n\t\t";
    printf STDERR "They contain original name, epsilon/k and sigma\n\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "\t\tThe transportdatafile has to be in CHEMKIN format.\n\t\t";
	printf STDERR "Each line contains speciesname, a number, epsilon/k in kelvin and\n\t\t";
	printf STDERR "sigma in angstrom all space or tab seperated.\n\n"; 
    } 
    printf STDERR "\t-o mechanismoutputfile\n\t\t";
    printf STDERR "filename where the new mechanism in FlameMaster format is written to.\n\t\t";
    printf STDERR "If not given 'newmechfile' is used.\n\n";
    
    printf STDERR "\tmechanism_file_name\n\t\t";
    printf STDERR "the name of the mechanism file to process\n";
    if ($exit_value==0) {   #it was called with the -h option
	printf STDERR "\t\tThis file has currently the following restrictions:\n\t\t";
	printf STDERR "A reaction block starts at a line where the 1. character is a number\n\t\t";
	printf STDERR "The reaction itself must be within the 1. line\n\t\t";
	printf STDERR "A reaction block must be written in subsequent lines and\n\t\t";
	printf STDERR "ends after the line where a '}' appears.\n\t\t";
	printf STDERR "A reaction block, that is commented out for the purpose of excluding species, \n\t\t";
	printf STDERR "starts with '#' followed by the number of the reaction\n\t\t";
	printf STDERR "All other single line comments must start with '#' followed by any\n\t\t";
	printf STDERR "character that is not a number, e.g. use '##' for permanent comments\n\t\t";
	printf STDERR "Comments over a range are allowed, if '/*' and '*/' appear on single lines\n\n";
    } 
    printf STDERR "output:\n\t";
    printf STDERR "A new mechanism in FlameMaster format is generated.\n\t";
    printf STDERR "with option -t also the file 'newthermofile' is generated.\n\n";
    printf STDERR "examples:\n\t";
    printf STDERR "$program_name -i includespecs.txt -x \"C4\\S*\" N-C7H16.mech\n\t";
    printf STDERR "$program_name -x \"\\S*C7\\S*\" -i \"\\S*\" -o N-C7H16.mech N-C7H16.mec\n\t";
    printf STDERR "\texcludes species that contain the the string \"C7\" \n\t";
    printf STDERR "\tsurrounded by any character\n\t";
    printf STDERR "$program_name -x \"\\S*H13\" -i \"\\S*\" -T manualtranslationtable.txt -t nc7h16.thermo -r nc7h16.trans N-C7H16.mech\n\t";
    printf STDERR "$program_name -x H2O2 -o CH4new.mech CH4.mech\n\n";
    exit;
}



sub readtransportdata { 
#(char *speciesorigname, double epsioverk[], double sigma[])
#    my ($chtransfile) = @_;
    my $endoffile = 0;
    my $i=0;
    my $line1="";
    my $epsi;
    my $sig;
    my $origname="";
    my $currentline="";
    my @token;

    if (!open FILETRANSIN, $chtransfile) { 
	print "\ncannot open transport data file: $chtransfile\n";
	die;
    }
    print "\nreading from transportdata file $chtransfile";
    my  $linenr=0;
    while (defined ($currentline = <FILETRANSIN>)) {
	chomp($currentline);
	my $validline=1;
	if ($currentline =~ /^!/ ) {$validline = 0;}
	else {
	    @token = split /\s+/, $currentline; 
	    #print "\ntoken: $token[0] $token[1] $token[2] $token[3] $token[4] ";
	    if ($#token < 0) {$validline=0;} 
	    elsif ($#token<4) { 
		print "\nline ",$linenr+1," has only ",$#token+1," token instead of 4"; 
		$validline=0;
	    }
#here I should check if token2 and token3 are indeed numerical values
	}    
	if ($validline) {
	    $origname=$token[0];
	    $epsi=$token[2];
	    $sig=$token[3];
	    $speciesintransportdata[$i] = $origname;
	    $epsioverk[$i]=$epsi;
	    $sigma[$i]=$sig;
	    printf "\nreading transport data: '$origname' epsilon: $epsi sigma: $sig";
	    $i++;
	}
	$linenr++;
    } # for all lines in FILETRANSIN
	$speciesintransportdata[$i] = "end";
    close(FILETRANSIN);
}


sub readspeciesinfofile { #this can be either a chemkinthermodatafile or a species-translationtable
#(FILE *FILETHERMIN,FILE *FILETRANSLATEIN,char *speciesorigname,char *speciesnewname,char *speciestoshow,int $maxnewspecieslength)   
    my $endoffile = 0;
    my $speciesfoundintranslation;
    my $itransl=0;
    my $i=0;
    my ($j,$k,$l);
    my $line1="";
    my $helperline="";
    my $line2="";
    my $line3="";
    my $line4="";
    my @Ahi;
    my @alo;
    my @Afit; #6,10
    my @bfit;
    my @fitx; #[16]
    my @fitcp; #[16]
    my @fitcpwrong; #[16];
    my ($lowerlimit,$upperlimit,$chcommtemp,$commtemp,$testT,$cplo,$cphi,$cperr,$H0RTlo,$H0RThi,$H0RTerr,$S0Rlo,$S0Rhi,$S0Rerr);
    my ($cpfit,$cporig,$H0RTfit,$H0RTorig,$H0RTorigwrong,$S0Rfit,$S0Rorig,$S0Rorigwrong);
    my $fitupper;
    my ($molmass,$epsioverkfinal,$sigmafinal);
    #my (@epsioverk,@sigma); #[$maxspecies];
    #my @speciesintransportdata; #[$maxspecies][$specieslength];
    my $origname;
    my $newname;
    my $chemicalname;
    my $newchemicalname;
    my $orignameuppercase;
    my $chemicalnameuppercase;
    my $helper;
    my @digitcharacter=(0,0);
    my $currentline;
    my $parameter;
    my @token; #[100][$tokenlength];
    my @numberatoms; #[4];
    my $currentatomnumberstring;
    my @atom;
    
    if ($chtransfile_exist) {
	#print "\nReading transportdata";
	readtransportdata;
    } else { #no transportfile
	print "\nNo transportdatafile was specified, transport parameter will be -1";
    }
    if (!open FILETHERMOUT, ">$fthermnameout") { 
	print "\ncannot open thermodata file for writing\n";
	exit;
    }
    
    if ($translatefile_exist) {
	if (!open FILETRANSLATEIN, $translatefile) { 
	    print "\ncannot open translation file: $translatefile\n";
	    die;
	}
	print "\nReading translationfile";
	my $linenr=0;
	while (defined ($currentline = <FILETRANSLATEIN>)) {
	    chomp($currentline);
	    my $validline=1;
	    if ($currentline =~ /^\#/ ) {$validline = 0;}
	    else {
		@token = split /\s+|=/, $currentline; 
		#print "\ntoken: $token[0] $token[1] $token[2] $token[3] $token[4] ";
		if ($#token < 0) {$validline=0;} 
		elsif ($#token<1) { 
		    print "\nline ",$linenr+1," has only ",$#token+1," token instead of 2"; 
		    $validline=0;
		} elsif ($#token>1) {
		    if ($token[2] !~ /^\#/) {
			print "\nline ",$linenr+1," has ",$#token+1," token instead of 2, use '#' before the 3. token";
			$validline = 0;
		    }
		}
	    }
	    if ($validline) {
		$speciesorigname[$itransl]=$token[0];
		$speciesnewname[$itransl]=$token[1];
		$itransl++;
	    }
	    $linenr++;
	} # all lines in the translation file
    } # only if there was a translation file
    
if ($chthermfile_exist) {        
    print "\nReading thermodata file";
    my $linenr=0;
    my $thermodatanr=0;
    if (!defined ($currentline = <FILETHERMIN>)) { 
	print "\nerror reading thermodata file\n";
	die;
    }
    if($currentline !~ /thermo/i) {
	print "\nthermodata file not in chemkin format, first line should start with 'thermo'";
    } else { #a valid chemkin thermodata file
	print "\nreading thermofile";
	$currentline = <FILETHERMIN>;
	$linenr++;
	my $endoffile=0;
	while (!$endoffile) {
	    my $endofempty=0;
	    while (!$endofempty) {
		$endofempty=1;
		if (!defined ($currentline = <FILETHERMIN>)) { 
		    $endoffile=1; 
		}
		else {
		    chomp ($currentline);
		    if ($currentline =~ /^!/) { $endofempty=0; }
		    if ($currentline =~ /^\s+$/) { $endofempty=0; }
		    $linenr++;
		}
	    }
	    $line1 = $currentline;
	    if ($line1 =~ /end/i) { $endoffile=1; } 
	    else {
		if (length($line1)<80) {
		    $endcharacter = " ";
		} else {
		    $endcharacter = substr($line1,79,1);
		}
		if ($endcharacter ne "1") {
		    print "\nerror: '1' expected in thermodata file at position 80 of line $linenr";
		    print "\n$line1";
		    $endoffile=1; 
		}
		if (!defined ($line2 = <FILETHERMIN>)) { $endoffile=1; }
		$linenr++;
		if (!defined ($line3 = <FILETHERMIN>)) { $endoffile=1; }
		$linenr++;
		if (!defined ($line4 = <FILETHERMIN>)) { $endoffile=1; }
		$linenr++;  
		
		if (!$endoffile) {
		    $origname = substr ($line1,0,16);
		    @dummy = split (/\s+/, $origname);
		    $origname = $dummy[0]; #this should be much easier
		    #print "\norigname '$origname'";
		    
		    $commtemp = 1000.0;
		    $lowerlimit = substr($line1,45,10);
		    $upperlimit = substr($line1,55,10);
		    if ($upperlimit>$maxuppertempforfit) { $upperlimit=$maxuppertempforfit; }
		    $chcommtemp = substr($line1,65,8); #attention: this is according to chemkin specs. Pitz thermo was 1 char longer
		    $Ahi[1]=substr($line2,0,15);
		    $Ahi[2]=substr($line2,15,15);
		    $Ahi[3]=substr($line2,30,15);
		    $Ahi[4]=substr($line2,45,15);
		    $Ahi[5]=substr($line2,60,15);
		    $Ahi[6]=substr($line3,0,15);
		    $Ahi[7]=substr($line3,15,15);
		    $alo[1]=substr($line3,30,15);
		    $alo[2]=substr($line3,45,15);
		    $alo[3]=substr($line3,60,15);
		    $alo[4]=substr($line4,0,15);
		    $alo[5]=substr($line4,15,15);
		    $alo[6]=substr($line4,30,15);
		    $alo[7]=substr($line4,45,15);
		    if($origname eq $speciestoshow) {
			print "\n$speciestoshow HT";
			for ($j=1;$j<=7;$j++) { print "\n$Ahi[$j]"; }
			print "\n$speciestoshow LT";
			for ($j=1;$j<=7;$j++) { print "\n$alo[$j]"; }
			die;
		    }
		    # The following statement makes sure that the common temp
                    # makes sense. Otherwise 1000K is assumed.
		    if ($chcommtemp eq "        ") { $chcommtemp = 1000.0;}
		    if ($chcommtemp < 50.0) { 
			print "\n$origname:\t $commtemp K is used as the common temp. instead of $chcommtemp K";
		    } elsif ($chcommtemp > 5000.0) { 
			print "\n$origname:\t $commtemp K is used as the common temp. instead of $chcommtemp K"; 
		    } else { # the same as chemkin common temp
			$commtemp = $chcommtemp;
		    }  
		    #$testT = $commtemp; #could be checked for the comm temp.
		    $testT = 1000.0;     #but we assume it is ok and also look
                                         #at 1000K. If it agrees there then we
                                         #do not refit but just ignore the 
                                         #common temp and use 1000K as common 
                                         #temp.
		    $cplo=calccpR(\@alo,\$testT);
		    $cphi=calccpR(\@Ahi,\$testT);
		    if ($cplo>0.00001) { 
			$cperr=($cphi-$cplo)/$cplo;
		    } else {
			print "\n no cp data available";
			$cperr=1;
		    }
		    $H0RTlo=calcH0RT(\@alo,\$testT);
		    $H0RThi=calcH0RT(\@Ahi,\$testT);
		    $H0RTerr=($H0RThi-$H0RTlo);
		    
		    $S0Rlo=calcS0R(\@alo,\$testT);
		    $S0Rhi=calcS0R(\@Ahi,\$testT);
		    $S0Rerr=($S0Rhi-$S0Rlo);
		    
		    if ((abs($cperr) > 0.002 || abs($H0RTerr) > 0.02 || abs($S0Rerr) > 0.02)) {
			printf "\n%s: Dev at %.0fK: cp(rel.): %.3f, H0/RT: %.3f, S0/R: %.3f, comm. T.: %.0fK ",$origname,$testT,$cperr,$H0RTerr,$S0Rerr,$commtemp; 
			if ($commtemp>1000.001) {
			    $fitupper=1;
			    print "\n...refit upper thd down to 1000K";
			} elsif ($commtemp<999.999) {
			    $fitupper=0;
			    print "\n...refit lower thd up to 1000K";
			} else { #commtemp = 1000
			    $fitupper=-1;
			    print "\nthermodata for $origname are not continuous at the cutofftemp. of 1000K";
			    printf STDERR"\nthermodata for $origname are not continuous at the cutofftemp. of 1000K\npress enter to continue";
			    read (STDIN,$_,1); 
			}
			if ($fitupper!=-1) {
			    #print "\nfitting....";
			    #generate points with old thermofunctions to be fit later
			    my $t=1000.0;
			    my $numberfitpoints = 15;
			    my $numberpar = 2;
			    for ($j=1;$j<=$numberfitpoints;$j++) {
				$fitx[$j]=$t/$factornumeric;
				if ($t<$commtemp) {
				    $fitcp[$j]=calccpR(\@alo,\$t);
				} else { #t>=commtemp
				    $fitcp[$j]=calccpR(\@Ahi,\$t);
				}
				# This part is just temporary. It is used to get thermodata some-
				# where in between the original nonfitted (with wrong common T)
				# and the new fitted thermodata
				if ($t<1000.0) {
				    $fitcpwrong[$j]=calccpR(\@alo,\$t);
				    #fitcpwrong[$j]=alo[1]+alo[2]*t+alo[3]*pow(t,2)+alo[4]*pow(t,3)+alo[5]*pow(t,4);
				} else { # t>=1000
				    $fitcpwrong[$j]=calccpR(\@Ahi,\$t);
				    #fitcpwrong[$j]=Ahi[1]+Ahi[2]*t+Ahi[3]*pow(t,2)+Ahi[4]*pow(t,3)+Ahi[5]*pow(t,4);
				}
				$fitcp[$j]=$gotocorrect*$fitcp[$j]+(1.0-$gotocorrect)*$fitcpwrong[$j];
				# until here temporary
				
				#print "\nfitcp[%i] @T=%f: %f\t",$j,$t,$fitcp[$j]);
				if ($fitupper) {
				    $t+=($upperlimit-1000.0)/($numberfitpoints-1);
				} else { #fitlower
				    $t-=(1000.0-$lowerlimit)/($numberfitpoints-1);
				}  
			    }
			    
			    #get the 1. derivative at 1000 to match it later
			    $t=1000.0;
			    my $derivHi=$Ahi[2]+2*$Ahi[3]*$t+3*$Ahi[4]*$t**2+4*$Ahi[5]*$t**3;
			    my $derivlo=$alo[2]+2*$alo[3]*$t+3*$alo[4]*$t**2+4*$alo[5]*$t**3;
			    my $fitcpprime;
			    if ($fitupper) {
				$fitcpprime=$gotocorrect*$derivlo+(1.0-$gotocorrect)*$derivHi;
			    } else { #fitlower
				$fitcpprime=$gotocorrect*$derivHi+(1.0-$gotocorrect)*$derivlo;
			    }
			    $fitcpprime*=$factornumeric;
			    
			    #	   printf("\n");
			    #generate the fit equation, these come from minimizing least square
			    #of deviations. They look like A*x=b.
			    #because we match the first point and slope of first point
			    #there will be only 3 parameters left for the least square-fit
			    my @factor;
			    for ($j=2;$j<=$numberfitpoints;$j++) {
				$factor[1]=1*$fitx[1]**2-2*$fitx[1]*$fitx[$j]+$fitx[$j]**2;
				$factor[2]=2*$fitx[1]**3-3*$fitx[1]**2*$fitx[$j]+$fitx[$j]**3;
				$factor[3]=3*$fitx[1]**4-4*$fitx[1]**3*$fitx[$j]+$fitx[$j]**4;
				$factor[$numberpar+1]=$fitcp[$j]-$fitcp[1]+($fitx[1]-$fitx[$j])*$fitcpprime;
				
				for ($l=1;$l<=$numberpar;$l++) {
				    for ($k=1;$k<=$numberpar;$k++) {
					if ($j==2) { $Afit[$k][$l]=0.0; }
					$Afit[$k][$l]+=$factor[$k]*$factor[$l];
				    }
				    if ($j==2) { $bfit[$l]=0.0; }
				    $bfit[$l]+=$factor[$numberpar+1]*$factor[$l];
				}
			    }
			    # A is the matrix built from the Gaussian sums. They are
			    # now more complicated than just [x],[x2],[xx]... because 
			    # the value and slope of the first point comes in
			    # x is the set of polynomial parameters we want to obtain 
			    # b comes also from Gaussian sums [y],[yx],[yx2]... 
			    #now get the answer
			    (*parfit)=solvelineqs(\@Afit,\@bfit,\$numberpar);
			    $parfit[5]=$parfit[3];
			    $parfit[5]=0.0; 
			    $parfit[4]=$parfit[2];
#  	   $parfit[4]=0.0;  
			    $parfit[3]=$parfit[1];
			    $parfit[2]=$fitcpprime-2*$parfit[3]*$fitx[1]-3*$parfit[4]*$fitx[1]**2-4*$parfit[5]*$fitx[1]**3;
			    $parfit[2]=$fitcpprime-2*$parfit[3]*$fitx[1]-3*$parfit[4]*$fitx[1]**2; 
#  	   $parfit[2]=$fitcpprime-2*$parfit[3]*$fitx[1]; 
			    $parfit[1]=$fitcp[1]-$fitx[1]*$fitcpprime+1*$parfit[3]*$fitx[1]**2+2*$parfit[4]*$fitx[1]**3+3*$parfit[5]*$fitx[1]**4;
			    $parfit[1]=$fitcp[1]-$fitx[1]*$fitcpprime+1*$parfit[3]*$fitx[1]**2+2*$parfit[4]*$fitx[1]**3;
#  	   $parfit[1]=$fitcp[1]-$fitx[1]*$fitcpprime+1*$parfit[3]*$fitx[1]**2;
			    $parfit[2]=$parfit[2]/$factornumeric**1;
			    $parfit[3]=$parfit[3]/$factornumeric**2;
			    $parfit[4]=$parfit[4]/$factornumeric**3;
			    $parfit[5]=$parfit[5]/$factornumeric**4;
			    
			    #get H0 with equating orig and fit at 1000K
			    $testT=1000.0;
			    if (!$fitupper) {
				$H0RTorig=calcH0RT(\@Ahi,\$testT);
				$parfit[6]=$Ahi[6];
				$S0Rorig=calcS0R(\@Ahi,\$testT);
				$parfit[7]=$Ahi[7];
				$H0RTorigwrong=calcH0RT(\@alo,\$testT);
				$S0Rorigwrong=calcS0R(\@alo,\$testT);
			    } else { #fitupper
				$H0RTorig=calcH0RT(\@alo,\$testT);
				$parfit[6]=$alo[6];
				$S0Rorig=calcS0R(\@alo,\$testT);
				$parfit[7]=$alo[7];
				$H0RTorigwrong=calcH0RT(\@Ahi,\$testT);
				$S0Rorigwrong=calcS0R(\@Ahi,\$testT);
			    }
			    $H0RTfit=calcH0RT(\@parfit,\$testT);
			    $parfit[6]+=(($gotocorrect*$H0RTorig+(1-$gotocorrect)*$H0RTorigwrong)-$H0RTfit)*$testT; # correction of parameter 6
			    $S0Rfit=calcS0R(\@parfit,\$testT);
			    $parfit[7]+=(($gotocorrect*$S0Rorig+(1-$gotocorrect)*$S0Rorigwrong)-$S0Rfit); # correction of parameter 7
			    
			    $testT=1000.0;
			    my $numbercheckpoints = 10;
			    for($j=1;$j<=$numbercheckpoints;$j++) {
				# first the orig chemkin
				if ($testT<$commtemp) {
				    $cporig=calccpR(\@alo,\$testT);
				} else { #testT>=commtemp
				    $cporig=calccpR(\@Ahi,\$testT);
				}
				#print "\ncporig @ T=%f: %f\t",testT,cporig);
				# then my fit
				$cpfit=calccpR(\@parfit,\$testT);
				#print "cpfit @ T=%f: %f",testT,cpfit);
				$cperr=($cporig-$cpfit)/$cporig;
				
				if ($testT<$commtemp) {
				    $H0RTorig=calcH0RT(\@alo,\$testT);
				} else { # testT>=commtemp
				    $H0RTorig=calcH0RT(\@Ahi,\$testT);
				}
				$H0RTfit=calcH0RT(\@parfit,\$testT);
				$H0RTerr=($H0RTorig-$H0RTfit);
				
				if ($testT<$commtemp) {
				    $S0Rorig=calcS0R(\@alo,\$testT);
				} else { #testT>=commtemp
				    $S0Rorig=calcS0R(\@Ahi,\$testT);
				}
				$S0Rfit=calcS0R(\@parfit,\$testT);
				$S0Rerr=($S0Rorig-$S0Rfit);
				
				if (abs($cperr) > 0.005  || abs($H0RTerr) > 0.05 || abs($S0Rerr) > 0.05) {
				    printf "\n%s:\tfit at %.1f K: cp(rel.): %.3f, H0/RT: %.3f, S0/R: %.3f, given comm. temp.: %.1f",$origname,$testT,$cperr,$H0RTerr,$S0Rerr,$commtemp; 
				}
				if ($fitupper) {
				    $testT+=($upperlimit-1000.0)/($numbercheckpoints-1);
				} else { #fitlower
				    $testT-=(1000.0-$lowerlimit)/($numbercheckpoints-1);
				}  
			    }
			    if ($fitupper) {
				for ($j=1;$j<=7;$j++) { $Ahi[$j]=$parfit[$j]; }
			    } else { #fitlower
				for ($j=1;$j<=7;$j++) { $alo[$j]=$parfit[$j]; }
			    }
			}   
		    } #if a refit was necessary
		    
		    #now we build the correct chemical name from the chemkin thermodata atom information
		    $chemicalname ="";
		    my $totalnumberatoms=0;
		    for ($j=0;$j<4;$j++) {
			$atom[$j]=substr($line1,24+$j*5,2);
			$atom[$j]=uc($atom[$j]);
			@atomjunk = split(/\s+/, $atom[$j]);
			$atom[$j] = $atomjunk[0];  #check if this could be easier
			#print "\natom: $atom[$j]";
			$numberatoms[$j]=substr($line1,26+$j*5,3);
			#print "\nnumberatomsstring: $numberatomsstring[$j]";
			#@dummy = split (/\s+/,$numberatomsstring[$j]);
			#$numberatoms[$j]=$dummy[0];
			if ($numberatoms[$j] eq "   ") { $numberatoms[$j]=0 };
			#     print "\natom: '%s' number: '%i'",@atom,$numberatoms[$j];
			$numberatoms[$j] =~ s/\s+//g;
			if ($numberatoms[$j] > 0) {
			    $chemicalname=$chemicalname.$atom[$j];
			    if ($numberatoms[$j] > 1) { $chemicalname=$chemicalname.$numberatoms[$j]; }
			}
			$totalnumberatoms+=$numberatoms[$j];
		    }
		    #print "\n chemical name: '$chemicalname'";
		    # now we will build the new name suitable for flamemaster. Before that we
		    # check if it is already in the translationtable b/c then we take that name
		    # and don't worry any further.
		    $speciesfoundintranslation = 0;
		    for(my $itable=0;$itable<$itransl && (!$speciesfoundintranslation);$itable++) {
			if ($speciesorigname[$itable] eq $origname) {
			    $newname=$speciesnewname[$itable];
			    $speciesfoundintranslation = 1;
			    printf "\nalready found in translation: '%s' => `%s`",$origname,$newname;
			}
		    }
		    if (!$speciesfoundintranslation) {     
			$newchemicalname="";
			my $statchemsearch = 's';
			my $actualchar;
			my $currentatom;
			$currentatom="";
			my $atomintableuppercase;
			my $currentatomuppercase;
			my ($currentatomnumber,$currentatomnumberfactor);
#       my  *pstart;
			printf "\n'%s'    =>    \t",$origname;
			$p=length($origname)-1;
			for(;($statchemsearch ne 'y') && ($statchemsearch ne 'n');$p--) {
			    #       print "\nnewchemicalname: '%s'",$newchemicalname;
			    $actualchar = substr($origname,$p,1);
			    #print "\nactualchar: '%c'",$actualchar;
			    if ($statchemsearch eq 's') {
				$currentatomnumberstring="";
				if($actualchar=~ /\d/) { 
				    $statchemsearch='d'; 
				    $currentatomnumberfactor=1;
				    $currentatomnumber=0;
				} elsif ($actualchar!~/\d/) {
				    $statchemsearch='a';
				    $currentatomnumber=1;
				} else { # not found
				    $statchemsearch='n';
				}
			    }
			    if ($statchemsearch eq 'd') {
				if($actualchar=~/\d/) { 
				    $currentatomnumber+=$currentatomnumberfactor*$actualchar;
				    $currentatomnumberfactor*=10;
				    $helper=$currentatomnumberstring;
				    $currentatomnumberstring=$actualchar;
				    $currentatomnumberstring=$currentatomnumberstring.$helper;
				} elsif($actualchar!~/\d/) {
				    $statchemsearch='a';
				    $currentatom="";
				} else { # not found
				    $statchemsearch='n';
				}
			    }
			    if ($statchemsearch eq 'a') {
				$helper=$currentatom;
				$currentatom=$actualchar;
				$currentatom=$currentatom.$helper;
				#	   strncat(currentatom,&actualchar,1);
				my $withinatomtable = 1;
				for($j=0;$withinatomtable;$j++) {
				    $currentatomuppercase=uc($currentatom);
				    $atomintableuppercase=uc($atom[$j]);
				    
				    if ($atomintableuppercase eq $currentatomuppercase) {
					$withinatomtable=0;
					$numberatoms[$j]-=$currentatomnumber;
					$currentatom="";
					if($numberatoms[$j]<0) {
					    $statchemsearch='n';
					} else { #atoms there
					    $totalnumberatoms-=$currentatomnumber;
					    $helper=$atom[$j];
					    $helper=$helper.$currentatomnumberstring.$newchemicalname;
					    $newchemicalname=$helper;
					    if($totalnumberatoms<=0) {
						$statchemsearch='y';
					    } else { #still searching
						$statchemsearch='s';
					    }
					}
				    }
				    if ($j>=3) { $withinatomtable=0; }
				}
			    }
			    if (($statchemsearch ne 'y') && (($p==0) || (length($currentatom)>=2))) { $statchemsearch='n'; } 
			    #	 print "\nnewchemicalname: '%s'",$newchemicalname);
			    #	 print "\nstatchem: '%c'",$statchemsearch);
			}
			if ($statchemsearch eq 'y') {
			    #print "\nnewchemicalname: $newchemicalname";
			    my $oldchars=$p+1;
			    $newname=substr ($origname,0,$oldchars);
			    #print "\nnewname: $newname";
			    if ($oldchars && ( substr($origname,$oldchars-1) ne '-')) { $newname=$newname."-"; }
			    # the following part was translated from c++, with Perl now it could be made much easier
			    for (my $iorigname=0;$iorigname<$oldchars-1;$iorigname++) {
				if (substr($newname,$iorigname,1)eq'-') { substr($newname,$iorigname,1)='-'; } 
				if (substr($newname,$iorigname,1)eq',') { substr($newname,$iorigname,1)='J'; }
				if (substr($newname,$iorigname,1)eq'(') { substr($newname,$iorigname,1)='G'; }
				if (substr($newname,$iorigname,1)eq')') { substr($newname,$iorigname,1)='G'; }
				if (substr($newname,$iorigname,1)eq'#') { substr($newname,$iorigname,1)='U'; }
				if (substr($newname,$iorigname,1)eq'*') { substr($newname,$iorigname,1)='D'; }
				if (substr($newname,$iorigname,1)eq'.') { substr($newname,$iorigname,1)='X'; }
			    }
			    $newname=uc($newname).$newchemicalname;
			} else { # chemicalname not found in origname
			    #print "\nchemicalname: '%s'",$chemicalname);
			    $newname=$origname;
			    my $origlength=length($newname);
			    for (my $iorigname=0;$iorigname<$origlength;$iorigname++) {
				if (substr($newname,$iorigname,1)eq'-') { substr($newname,$iorigname,1)='-'; }
				if (substr($newname,$iorigname,1)eq',') { substr($newname,$iorigname,1)='J'; }
				if (substr($newname,$iorigname,1)eq'(') { substr($newname,$iorigname,1)='G'; }
				if (substr($newname,$iorigname,1)eq')') { substr($newname,$iorigname,1)='G'; }
				if (substr($newname,$iorigname,1)eq'#') { substr($newname,$iorigname,1)='U'; }
				if (substr($newname,$iorigname,1)eq'*') { substr($newname,$iorigname,1)='D'; }
				if (substr($newname,$iorigname,1)eq'.') { substr($newname,$iorigname,1)='X'; }
			    }
			    $newname=uc($newname)."-".$chemicalname;
			}
			if (length($newname)>$maxnewspecieslength) {
			    $digitcharacter[1]+=1;
			    if($digitcharacter[1]>65+25) {
				$digitcharacter[1]=65+0;
				$digitcharacter[0]+=1;
				if($digitcharacter[0]>65+25) {
				    $digitcharacter[1]=65+0;
				    $digitcharacter[0]=65+0;
				}
			    }
			    $newname=uc(chr($digitcharacter[0]).chr($digitcharacter[1])."-".$chemicalname);
			}
			print "'$newname'";
			$speciesorigname[$itransl]=$origname;
			$speciesnewname[$itransl]=$newname;
			$itransl++;
		    } # if species not already found in translation table
		    
		    $molmass=-1;
		    $epsioverkfinal=-1;
		    $sigmafinal=-1;
		    if ($chtransfile_exist) {
			$stillsearching=1;
			for (my $lineintransdata=0;$stillsearching;$lineintransdata++) {
			    my $testname;
			    #printf "compare '%s' with '%s'\n",$origname,$speciesintransportdata[$lineintransdata];
			    if ($speciesintransportdata[$lineintransdata] eq $origname) {
				$stillsearching=0;
				$epsioverkfinal=$epsioverk[$lineintransdata];
				$sigmafinal=$sigma[$lineintransdata];
				#	   print "\nfound transportdata for %i: '%s'='%s', epsi: %f, sigma: %f",$lineintransdata,$origname,$speciesintransportdata[$lineintransdata],$epsioverkfinal,$sigmafinal);
			    }
			    $helperline=lc($speciesintransportdata[$lineintransdata]);
			    #trim(strtolower(helperline));
			    if (($helperline eq "end") || ($lineintransdata > $maxspecies)) {
				$stillsearching=0;
				printf STDERR "\n$origname not in transportdata file, -1 is used"; 
			    }
			}
		    }
		    printf FILETHERMOUT "%s % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e \n",$newname,$Ahi[1],$Ahi[2],$Ahi[3],$Ahi[4],$Ahi[5],$Ahi[6],$Ahi[7],$alo[1],$alo[2],$alo[3],$alo[4],$alo[5],$alo[6],$alo[7],$molmass,$epsioverkfinal,$sigmafinal;
		} # if 1 was found at the end of the chemkin thermodata line
	    } # for all lines in FILETHERMIN
	} # if chemkin format 
    } # if there is a thermodatafile
    close(FILETHERMOUT);
}
    $speciesorigname[$itransl]="end";
    if (!open FILESPECTRANSLOUT, ">$fspecnameout") { 
	print "\ncannot open output species translation file for writing\n";
	die;
    }
    
    for (my $j=0;$j<$itransl;$j++) { printf FILESPECTRANSLOUT "%s\t%s\n",$speciesorigname[$j],$speciesnewname[$j]; }
    close(FILESPECTRANSLOUT);
} #end sub readspeciesinfofile

sub calccpR {
# double par[], float T)
    local (*par,*T) = @_;;
    my $cp=$par[1]+$par[2]*$T+$par[3]*$T**2+$par[4]*$T**3+$par[5]*$T**4;
    return $cp;
}

sub calcH0RT {
#( double par[], float T)
    local (*par,*T) = @_;;
    my $H=$par[1]+$par[2]/2*$T+$par[3]/3*$T**2+$par[4]/4*$T**3+$par[5]/5*$T**4+$par[6]/$T;
    return $H;
}

sub calcS0R {
#( double par[], float T)
    local (*par,*T) = @_;;
    my $S=$par[1]*log($T)+$par[2]*$T+$par[3]/2*$T**2+$par[4]/3*$T**3+$par[5]/4*$T**4+$par[7];
    return $S;
}



sub solvelineqs {
#(double a[][10],double binp[],double x[],int n) {
  local (*a, *binp, *n) =@_;	
  my $msize = 10;
  #	float a[msize][msize];
#  my @inva[msize][msize];
#  my @l[msize][msize];
#  my @u[msize][msize];
  my @y;
  my @b;
  my ($i,$j,$k,$irow,$jcol);
  my $sum;
  my @x;
  # Initializing the arrays
  for($i=1;$i<=$n;$i++)
    {
      for($j=1;$j<=$n;$j++)
	{
	  $l[$i][$j]=0;
	  $u[$i][$j]=0;
	  $y[$i]=0;
	  $b[$i]=$binp[$i];
	  $x[$i]=0;
	}
    }
  # first column of l = first column of a
  for($i=1;$i<=$n;$i++) { $l[$i][1] = $a[$i][1]; }
  # first row of u can be obtained from a and l
  for($j=1;$j<=$n;$j++) { $u[1][$j] = $a[1][$j]/$a[1][1]; }
  for($j=2;$j<=$n; $j++) 
    {
      for($i=$j; $i<=$n; $i++) {	
	$sum = 0.0;	
	for ($k=1;$k<=$j-1;$k++) { $sum = $sum + $l[$i][$k]*$u[$k][$j]; }	
	$l[$i][$j] = $a[$i][$j]-$sum;
      }
      $u[$j][$j]= 1;
      $irow = $j;
      for($jcol=$irow+1;$jcol<=$n;$jcol++)
	{
	  $sum=0.0;
	  for($k=1;$k<=$irow-1;$k++) { $sum = $sum + $l[$irow][$k]*$u[$k][$jcol]; }
	  $u[$irow][$jcol] = ($a[$irow][$jcol] - $sum)/$l[$j][$j];
	}
    }
#   printf STDOUT "The array read is A = \n"; 
#   printarray($a,$n); 

#   printf STDOUT "The array L is \n"; 
#   printarray($l,$n); 
#   printf STDOUT "The array U is  \n"; 
#   printarray($u,$n);	 
  for($i=1;$i<=$n;$i++) 
    {
      for($j=1;$j<=$n;$j++) 
	{
	    if ( $i==$j) { $b[$j]= 1; }
	    else { $b[$j] = 0; }
	}
      (*y)=forsub(\@l,\@b,\@y,\$n);
      (*x)=backsub(\@u,\@y,\@x,\$n);
      for($j=1;$j<=$n;$j++) { $inva[$j][$i] = $x[$j]; }
    }
  #printf STDOUT "The inverse array is \n";
  #printarray(inva,n);
  # compute x = inva * binp
  for ($i=1;$i<=$n;$i++) {
    $x[$i]=0.0;
    for ($j=1;$j<=$n;$j++) {
      $x[$i]+=$inva[$j][$i]*$binp[$j];
      #printf("inva[%i][%i]: %g\n",j,i,inva[$j][$i]);
    }
    #printf("a[%i]=%g\n",i,x[$i]);
  }
  return (\@x);
}

sub printarray {
#( double a[][10], int n)
#{
    local (*a, *n) =@_;	
    my ($i,$j);
    for ( $i =1; $i<=$n; $i++)
    {
	for($j = 1; $j<=$n; $j++) { 
	   printf "%7.3f\t",$a[$i][$j];
	}
	print "\n";
    }
}


sub forsub {
#( double l[][10], double b[], double y[], int n)
  local (*l, *b, *y, *n) =@_;	
  my $i=0; my $j=0;
  my $sum=0;
  for($j=1;$j<=$n;$j++)
  {
      my $sum = 0;
      for($i=1;$i<=($j-1);$i++) {  $sum = $sum + $l[$j][$i]*$y[$i]; }
      $y[$j] = ($b[$j]-$sum)/$l[$j][$j];
  }
  return (\@y); 
}

sub backsub {
# double u[][10], double y[], double x[], int n) # does backsubstitution
  local (*u, *y, *x, *n) =@_;	
  my $i=0; my $j=0;		# counter variables
  my $sum = 0;
  
  for($j=$n;$j>=1;$j--)
  {
      my $sum = 0;
      for($i=$j+1;$i<=$n;$i++) { $sum = $sum + $x[$i]*$u[$j][$i]; }
      $x[$j] = ($y[$j] - $sum)/$u[$j][$j];
  }
  return (\@x);
}


use Env (HOST);

$specincl_on = 0; # flag when parsing the commandline
$specexcl_on = 0;
$specsteadystates_on = 0;
$chthermfile_on = 0;
$translatefile_on = 0;
$chtransfile_on = 0;
$showaspecies_on = 0;
$speciesmaxchar_on = 0;
$steadystatesfileout_on = 0;
$mechoutfile_on = 0;
$increasea_on = 0;

$specincl_exist = 0;  # if there is a species or a file to include species
$specexcl_exist = 0;
$specsteadystates_exist = 0;
$chthermfile_exist = 0;
$translatefile_exist = 0;
$speciestoshow_exist = 0;
$increasea_exist = 0;
#$refitifnot1000K = 1;

$currentline="";
$mechfile="";



$steadystatesnotfound=0;
$countreactions=0;
# The following section reads the commandline arguments, mainly several
# filenames. It remembers then if there was a CHEMKIN thermodatafile 
# specified, a CHEMKIN transportdatafile and so on. The main argument is 
# the FlameMaster mechanismfile
#print "number args: '$#ARGV'\n";
for ( my $ix=0; $ix <= $#ARGV; ++$ix) { #reads all arguments from the commandline string by string
    my $pchar = $ARGV[$ix];
    $sw1=substr($pchar,0,1);
      #print "sw1 '$sw1'\n";
  SWITCH1: {
      if($sw1 eq '-') {
	  $sw2=substr($pchar,1,1);
	SWITCH2: {
	    #print "sw2 '$sw2'\n";
	    if($sw2 eq 'h') { usage(0); }
	    if($sw2 eq 'i') {
		$specincl_on = 1;
		last SWITCH2; 
	    }
	    if($sw2 eq 'x') {
		$specexcl_on = 1;
		last SWITCH2; 
	    }
	    if($sw2 eq 'X') {
		$specsteadystates_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 't') {
		$chthermfile_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'T') {
		$translatefile_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'r') {
		$chtransfile_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 's') {
		$showaspecies_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'S') {
		$steadystatesfileout_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'c') {
		$speciesmaxchar_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'o') {
		$mechoutfile_on = 1;
		last SWITCH2;
	    }
	    if($sw2 eq 'a') {
		$increasea_on = 1;
		last SWITCH2;
	    }
	    if ($sw2 !~ 'h|i|x|X|t|T|r|s|S|c|o|a') {
		printf STDERR "$program_name: error: unrecognized option: - $pchar\n\n";
		usage(-1);
	    }
	    last SWITCH2;
	}
      }
      if ($sw1 ne '-') {
	  if ($specincl_on) {
	      $specincl_on = 0;
	      $specincl_exist = 1; # so that we know later if there are species to include
	      $specincl = $pchar;
	  } 
	  elsif ($specexcl_on) {
	      $specexcl_on = 0;
	      $specexcl_exist = 1;
	      $specexcl = $pchar;
	  }
	  elsif ($specsteadystates_on) {
	      $specsteadystates_on = 0;
	      $specsteadystates_exist = 1;
	      $specsteadystates = $pchar;
	  }
	  elsif ($chthermfile_on) {
	      $chthermfile_on = 0;
	      $chthermfile_exist = 1;
	      $chthermfile = $pchar;
	  } 
	  elsif ($translatefile_on) {
	      $translatefile_on = 0;
	      $translatefile_exist = 1;
	      $translatefile = $pchar;
	  } 
	  elsif ($chtransfile_on) {
	      $chtransfile_on = 0;
	      $chtransfile_exist = 1;
	      $chtransfile = $pchar;
	  }
	  elsif ($showaspecies_on) {
	      $showaspecies_on = 0;
	      $speciestoshow_exist = 1;
	      $speciestoshowinp = $pchar;
	  }
	  elsif ($speciesmaxchar_on) {
	      $speciesmaxchar_on = 0;
	      $maxnewspecieslength = $pchar;
	  }
	  elsif ($steadystatesfileout_on) {
	      $steadystatesfileout_on = 0;
	      $fsteadystatesnameout=$pchar;
	  }
	  elsif ($increasea_on) {
	      $increasea_on = 0;
	      $increasea_exist=1;
	      $increasereactionlabel=$pchar;
	  }
	  elsif ($mechoutfile_on) {
	      $mechoutfile_on = 0;
	      $fmechnameout=$pchar;
	  }
	  else { $mechfile = $pchar; }
	  last SWITCH1;
      }
  }
}

if ( $mechfile eq "" ) {
    printf STDERR "$program_name: error: no file specified for processing.\n\n";
    usage(-3);
}

if (!open FILEMECHIN, $mechfile) { 
    printf "\ncannot open mechanism file '%s' \n", $mechfile;
    exit;
}

if ($chthermfile_exist) {
    if (!open FILETHERMIN, $chthermfile) { 
	print "\ncannot open thermodata file: '$chthermfile'\n";
	exit;
    }
}

$speciestoshow="";
if ($speciestoshow_exist) {
    if ($chthermfile_exist) {
	$speciestoshow=$speciestoshowinp;
	printf "\nspecies to show: '%s'",$speciestoshow;
    } else { # no chemkin thermodata file was specified
	print "\nNo thermo data file was specified to read the species from";
	die;
    }
}

# The following section reads the species to be included and then the 
# species to be excluded. Each time it is checked if a species was given
# directly at the commandline or in a file. In the latter case the species
# from the file are read in to an array.
$numberspecincl=1; # if species given in the comandline it will stay 1
if ($specincl_exist) {
    if (!open FILESPECINCL, $specincl) { 
	#printf STDERR "\nmanually included species: '$specincl'"; # exist, 
	#use the species from commandline
	$specinclnames[0]=$specincl;
    }
    else {
	my $endoffile = (!defined ($currentline = <FILESPECINCL>)); 
        for($numberspecincl=0;!$endoffile;) {
	    #read sp from species include file
	    @p = split (/\s+/,$currentline); # takes the first word
	    if (!defined ($p[0])) {$p[0]="";}
	    if (($p[0] ne "") && ($p[0] !~ /^\#/)) { # only if line is not empty
		$specinclnames[$numberspecincl]=$p[0];
		#print "\n",$specinclnames[$numberspecincl];
		$numberspecincl++;
	    }
	    $endoffile = (!defined ($currentline = <FILESPECINCL>)); 
	}
	printf STDERR "\nnumber of included species: $numberspecincl"; 
    }
    for($i=0;$i<$numberspecincl;$i++) {
	#$specinclnames[$i] =~ s/\*/\.\*/g;  #for perl '*' is replaced with '.*'
	print "\nincluded species #",($i+1),": '",$specinclnames[$i],"'";
	}
}
my $numberspecexcl=1;
if ($specexcl_exist) { # same as block above
    if (!open FILESPECEXCL, $specexcl) { 
	#printf STDERR "\nmanually excluded species: '$specexcl'";   
	$specexclnames[0]=$specexcl;
    }
    else {
	my $endoffile = (!defined ($currentline = <FILESPECEXCL>)); 
	for($numberspecexcl=0;!$endoffile;) {
	    @p= split (/\s+/,$currentline); #maybe easier    
	    if (!defined ($p[0])) {$p[0]="";}
	    if (($p[0] ne "") && ($p[0] !~ /^\#/)) { # only if line is not empty
		$specexclnames[$numberspecexcl]=$p[0];
		$numberspecexcl++;
	    }
	    $endoffile = (!defined ($currentline = <FILESPECEXCL>)); 
	}
	printf STDERR "\nnumber of excluded species: $numberspecexcl"; 
    }
    for($i=0;$i<$numberspecexcl;$i++) {
	#$specexclnames[$i] =~ s/\*/\.\*/g;
	print "\nexcluded species #",($i+1),": '",$specexclnames[$i],"'";
	}
}
my $numberspecsteadystates=1;
if ($specsteadystates_exist) { # same as block above
    if (!open FILESTEADYSTATESPECIES, $specsteadystates) { 
	#printf STDERR "\nmanually chosen steady states species: '$specsteatystates'";   
	$specsteadystatesnames[0]=$specsteadystates;
    }
    else {
	my $endoffile = (!defined ($currentline = <FILESTEADYSTATESPECIES>)); 
	for($numberspecsteadystates=0;!$endoffile;) {
	    @p= split (/\s+/,$currentline); #maybe easier    
	    if (!defined ($p[0])) {$p[0]="";}
	    if ($p[0] ne "" && $p[0] !~ /^\#/) { # only if line is not empty
		$specsteadystatesnames[$numberspecsteadystates]=$p[0];
		$numberspecsteadystates++;
	    }
	    $endoffile = (!defined ($currentline = <FILESTEADYSTATESPECIES>)); 
	}
	printf STDERR "\nnumber of steady state species: $numberspecsteadystates"; 
    }
#    for(my $i=0;$i<$numberspecsteadystates;$i++) {
#	$p=" ".$specsteadystatesnames[$i]." ";
#	$p=~s/ \*//;
#	$p=~s/\* //;
#	$specsteadystatesnames[$i]=$p;
#	printf STDOUT "\nsteady state species \#",$i+1,": '",$specsteadystatesnames[$i],"'";
#    }
}
if (!open FILEMECHOUT, ">$fintmechnameout") { 
    print "\ncannot open mechanism output file\n";
    die;
}
print "\nchthermfile: $chthermfile_exist, translatefile: $translatefile_exist";  
if ($chthermfile_exist || $translatefile_exist ) { readspeciesinfofile; }
#(FILETHERMIN,FILETRANSLATEIN,*speciesorigname,*speciesnewname,speciestoshow,$maxnewspecieslength); 
# If a CHEMKIN thermodata file was specified, then the species there are read
# and matched with species in the mechanism, then converted to a new name on
# basis of the here constructed chemical name.
# Also in addition a converted FlameMaster thermodatafile is written, 
# using, if specified, a CHEMKIN transport data file, for adding the 
# 2 transport coefficients. 

# Now we go through the mechanism, line by line. First we check which type
# of line we are in and then the line is converted.
#bool endoffile = false;
my $statusofparsing = 0;  
# 0: looking for new reaction line 
# 1: within a follow up line of an orig. commented reaction with the mission to uncomment it
# 2: within a follow up line of an active reaction line with the mission to comment it out
# 3: within a follow up line of a reaction and leaving it as it is;
# 4: within the steady states species section of the mechanism (firstline)
# 5: within the steady states species section of the mechanism
for (my $linenr=0; defined ($currentline = <FILEMECHIN>); $linenr++) {
    #print "\ncurrentline: '$currentline'";
    $firstlineofcommentreaction=0;
    $firstlineofreaction=0; 
    if ($statusofparsing == 0) {  # it is looking for a new reaction
	$isthirdbodyline=0;
	if ($currentline =~ /let M/i ) {$isthirdbodyline = 1;}
	$issteadystates=0;
	if ($currentline =~ /let steady states be/i ) {$issteadystates = 1;}
	if ($issteadystates) { $statusofparsing=4; }
	$firstlineofcommentreaction = ($currentline=~/^\#\d+\w?:/);
	# this means it is a valid reaction just commented out for the purpose of 
	# excluding species. Note that there can be other comments or reactions 
	# that are meant to be commented out permanentely.
	$firstlineofreaction = ($currentline=~/^\d+\w?:/);
	# this means it is a valid reaction 
	$usethisreaction=0;
	if ($firstlineofreaction) { 
	    @dummy=split(/:/,$currentline); 
	    $reactionlabel[$countreactions]=$dummy[0]; #could be easier
	    $usethisreaction=1;
	    #print "\n(3)ok: $currentline";
	}
	if ($firstlineofcommentreaction) { 
	    @dummy=split(/\#|:/,$currentline); #this could fail if "#2 4f:" or "#24 f :" or similar
	    $reactionlabel[$countreactions]=$dummy[1];
	    $usethisreaction=-1;
	    #print "\n(4)comment: $currentline";
	    #print "\nreactionlabel '",$reactionlabel[$countreactions],"'";
	}
	#print "\nreactionlabel nr $countreactions: ",$reactionlabel[$countreactions];
	$paragraph=$currentline;
#	$reactiontext="";
#	if (($firstlineofreaction || $firstlineofcommentreaction || $isthirdbodyline) && ($chthermfile_exist || $translatefile_exist)) {
#	}
	if ($firstlineofreaction || $firstlineofcommentreaction) {$statusofparsing=1 }
	##print "\n$currentline";
    }
    
    elsif ($statusofparsing == 4) { # start line for steady states
	printf FILEMECHOUT "let steady states be \n";
	for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	    printf FILEMECHOUT "$specsteadystatesnames[$iss] \n";
	}
	if ($currentline=~/\./) {
	    $paragraph=".\n";
	    $statusofparsing=0;
	} else { # not the end yet
	    $paragraph="";
	    $statusofparsing=5;
	}
    } 
    elsif ($statusofparsing == 5) { # follow up lines of steady states
	if ($currentline=~/\./) {
	    $paragraph=".\n";
	    $statusofparsing=0;
	} else { # not the end yet
	    $paragraph="";
	}
    } 
    else { # follow up lines will stay the same, but we still check where the reaction ends
	if ($firstlineofcommentreaction) { printf STDERR "\nerror: I am looking for the end of a reaction ('}'),\n but line ",$linenr+1," looks like the beginning of a new commented reaction"; }
	if ($firstlineofreaction) { printf STDERR "\nerror: I am looking for the end of a reaction ('}'),\n but line ",$linenr+1," looks like a beginning of a reaction"; }
	$paragraph=$paragraph.$currentline;
    }
    #printf STDOUT "currentline: $currentline";
    #this is done no matter which line
    
    if ((($statusofparsing<4) && ($currentline =~ /\}/)) || $isthirdbodyline) {
	# In this case, we check first if there are species to be translated.
	# This is done by going through all species in the array, that were
	# read from the thermodata file. Each species then is searched for in
	# the current line, and replaced, until no more is found in the line
	if ($chthermfile_exist || $translatefile_exist) {
	    for (my $lineinspectable=0;$speciesorigname[$lineinspectable] ne "end";$lineinspectable++) {
		#printf("specorig: '%s'",speciesorigname[lineinspectable]);
		$pold=0;
		$newline=""; # a new line will be built while searching the line 
		#the following could be completely rewritten to use Perl's capabilities
		$oldname = $speciesorigname[$lineinspectable];
		$oldname =~ s/\*/\\\*/g;
		$oldname =~ s/\./\\\./g;
		$oldname =~ s/\(/\\\(/g;
		$oldname =~ s/\)/\\\)/g;
		$oldname =~ s/\#/\\\#/g;
		#printf "\noldname: $oldname"; 
		for(;$paragraph=~/($oldname)/g;) {
		    $posfound=length($`);
		    $lengthmatch = length($&);
		    $before = substr($paragraph,$posfound-1,1);
		    $before2 = substr($paragraph,$posfound-2,1);
		    $after = substr($paragraph,$posfound+$lengthmatch,1);
		    $after2 = substr($paragraph,$posfound+$lengthmatch+1,1);
		    my $speciesfrontok=0;
		    if(($after =~ (/ |\t|\[|]|{|\+/)) || ($after2 eq "->")) {
			if(($before2 =~ (/ |\t|\[|\]|:|\+|>/)) && ($before=~/\d/)) {
			    printf STDERR "\nassuming $before $speciesorigname[$lineinspectable] in line:\n%s",$paragraph;
			    $speciesfrontok=2;
			}
			if($before =~ / |\t|\[|\]|:|\+|>/) { 
			    $speciesfrontok=1; 
			}
			if($speciesfrontok) {
			    # here a species was found and it was surrounded by possible 
			    # delimiters or a single digit in front
			    $newline=$newline.substr($paragraph,$pold,$posfound-$pold);
			    if($speciesfrontok==2) { $newline=$newline." "; }
			    $newline=$newline.$speciesnewname[$lineinspectable];
			} else { # the portion until here stays the same
			    $newline=$newline.substr($paragraph,$pold,$posfound-$pold+$lengthmatch);
			}
		    } else { # the portion until here stays the same
			$newline=$newline.substr($paragraph,$pold,$posfound-$pold+$lengthmatch);
		    }
		    $pold=pos $paragraph;
		}
		$newline=$newline.substr($paragraph,$pold);
		$paragraph=$newline;
	    }
	}
	
	if (!$isthirdbodyline) {
	    @dummy = split(/\{/,$paragraph);  #}
	    $reactionpart=$dummy[0];
	    #print "\n(1)reactionpart: ",$reactionpart;
	    
	    if ($usethisreaction == -1 && $specincl_exist) { #test if species to be included only if reaction was commented out
		for (my $counter=0;($counter<$numberspecincl) && ($usethisreaction!=2);$counter++) {
		    $exitsearchinline=0;
		    for(;$paragraph=~/($specinclnames[$counter])/g && $exitsearchinline==0;) {
			#printf "\nnumberspecincl: $numberspecincl counter: $counter specincl: '%s'",$specinclnames[$counter];
			$posfound=length($`);
			$lengthmatch = length($&);
			$match=substr($paragraph,$posfound,$lengthmatch);
			$before = substr($paragraph,$posfound-1,1);
			$before2 = substr($paragraph,$posfound-2,1);
			$after = substr($paragraph,$posfound+$lengthmatch,1);
			$after2 = substr($paragraph,$posfound+$lengthmatch+1,1);
			my $speciesfrontok=0;
			if(($after =~ (/ |\t|\[|]|{|\+/)) || ($after2 eq "->")) {
			    if(($before2 =~ (/ |\t|\[|\]|:|\+|>/)) && ($before=~/\d/)) {
				printf STDERR "\nassuming $before $match in line:\n%s",$paragraph;
				$speciesfrontok=2;
			    }
			    if($before =~ / |\t|\[|\]|:|\+|>/) { 
				$speciesfrontok=1; 
			    }
			    if($speciesfrontok>0) {
				# here a species was found and it was surrounded by possible 
				# delimiters or a single digit in front
				$usethisreaction=2;
				$foundspecies=$match;
				$exitsearchinline=1;
			    }
			}
		    }
		}
	    }
	    if ($usethisreaction != -1 && $specexcl_exist) { #test if species to be excluded only if reaction was not commented out
		for (my $counter=0;($counter<$numberspecexcl) && ($usethisreaction!=-2);$counter++) {
		    $exitsearchinline=0;
		    for(;$paragraph=~/($specexclnames[$counter])/g && $exitsearchinline==0;) {
			$posfound=length($`);
			$lengthmatch = length($&);
			$match=substr($paragraph,$posfound,$lengthmatch);
			$before = substr($paragraph,$posfound-1,1);
			$before2 = substr($paragraph,$posfound-2,1);
			$after = substr($paragraph,$posfound+$lengthmatch,1);
			$after2 = substr($paragraph,$posfound+$lengthmatch+1,1);
			my $speciesfrontok=0;
			if(($after =~ (/\s|\[|]|{|\+/)) || ($after.$after2 eq "->")) {
			    if(($before2 =~ (/\s|\[|\]|:|\+|>/)) && ($before=~/\d/)) {
				printf STDERR "\nassuming $before $match in line:\n%s",$paragraph;
				$speciesfrontok=2;
			    }
			    if($before =~ /\s|\[|\]|:|\+|>/) { 
				$speciesfrontok=1; 
			    }
			    if($speciesfrontok>0) {
				# here a species was found and it was surrounded by possible 
				# delimiters or a single digit in front
				$usethisreaction=-2;
				$foundspecies=$match;
				$exitsearchinline=1;
			    }
			}
		    }
		}
	    }
	    
	    if (($paragraph =~ /^\#\d+\w?:/) && ($usethisreaction>0)) {
		$paragraph =~ s/^\#//;
		$paragraph =~ s/\n\#/\n/g;
		printf STDOUT "\nACTIVATED REACTION: '$foundspecies' \t$paragraph";
	    }
	    if (($paragraph =~ /^\d+\w?:/) && ($usethisreaction<0)) {
		$paragraph = "#".$paragraph;
		$paragraph =~ s/\n/\n\#/g;
		$paragraph =~ s/\n\#$/\n/g;
		printf STDOUT "\n#DISABLED REACTION: '$foundspecies' \t$paragraph";
	    }
	    
	    if ($usethisreaction>0) {
		#print "\n(2)paragraph $paragraph";
		if ($specsteadystates_exist) {
		    # Here we check the steady state species, if any of them appears in the reaction
		    # This is done by going through all given steady state species in the array.
		    # Each species then is searched for in the current line. 
		    # If it was found...
		    my $linebeforespecies;
		    $linebeforespecies=$currentline;
		    #printf("check steady states in reaction: %i\n",countreactions);
		    
		    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
			$steadystatespecfound[$iss]=0;
			
			for(;$paragraph=~/($specsteadystatesnames[$iss])/g;) {
			    $posfound=length($`);
			    $lengthmatch = length($&);
			    $before = substr($paragraph,$posfound-1,1);
			    $before2 = substr($paragraph,$posfound-2,1);
			    $after = substr($paragraph,$posfound+$lengthmatch,1);
			    $after2 = substr($paragraph,$posfound+$lengthmatch+1,1);
			    my $speciesfrontok=0;
			    if(($after =~ (/ |\t|\[|]|{|\+/)) || ($after2 eq "->")) {
				if(($before2 =~ (/ |\t|\[|\]|:|\+|>/)) && ($before=~/\d/)) {
				    printf STDERR "\nassuming $before $specsteadystatesnames[$iss] in line:\n%s",$paragraph;
				    $speciesfrontok=2;
				}
				if($before =~ / |\t|\[|\]|:|\+|>/) { 
				    $speciesfrontok=1; 
				}
				if($speciesfrontok>0) {
				    # here a species was found and it was surrounded by possible 
				    # delimiters or a single digit in front
				    $linebeforespecies=substr($reactionpart,0,$posfound);
				    if($linebeforespecies=~">|=") {
					$steadystatespecfound[$iss]=1;
					# it is a product
				    } else {   # it is a reactant
					$steadystatespecfound[$iss]=-1;
					print "\nreactant '",$specsteadystatesnames[$iss],"' in $countreactions";
					#now we know which of the steady state spec. were found as reactants in the reaction
				    }
				}
			    }
			}
		    }
		    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
			$numberreactionforss[$iss]=0;
		    }
		    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
			if($steadystatespecfound[$iss]!=0) {

			    $reactionforss[$iss][$numberreactionforss[$iss]]=$countreactions;
			    #print "\n$reactionforss[$iss][$numberreactionforss[$iss]]",$reactionforss[$iss][$numberreactionforss[$iss]];
			    $stoichcoeff[$iss][$numberreactionforss[$iss]]=$steadystatespecfound[$iss];
			    $numberreactionforss[$iss]++; 
			}
		    }
		}
		
		if ($increasea_exist || $specsteadystates_exist) {		    
		    @dummy = split(/\{|\}/,$paragraph);  #}
		    $reactionpart=$dummy[0];
		    $parameterpart=$dummy[1];
		    $partafterparameter=$dummy[2];
		    #print "\nparameterpart $countreactions: '$parameterpart'";
		    $correcteda=0;
		    my $freqfactor=-1e7;
		    my $tempcoeff=-1e7;
		    my $actenergy=-1e7;
		    #$reactiontext=$reactiontext.$currentline;  #it adds lines as long a reaction lasts
		    @tokens=split(/\s+/,$parameterpart);
		    #print "\ntokens: '$tokens[0]'$tokens[1]'$tokens[2]'$tokens[3]'$tokens[4]'$tokens[5]'$tokens[6]'$tokens[7]'$tokens[8]'";
		    #print "$tokens[9]'$tokens[10]'$tokens[11]'$tokens[12]'$tokens[13]'$tokens[14]'$tokens[15]'$tokens[16]'$tokens[17]'";
		    for (my $tokeninline=0;$tokeninline<=($#tokens-1);$tokeninline++) {
			#printf("tokens[tokeninline]: '%s'\n",tokens[tokeninline]);
			#    if (!defined ($tokens[$tokeninline+1])) { $tokens[$tokeninline+1] = ""; } 
			if ($tokens[$tokeninline+1] eq "=") {
			    if(uc($tokens[$tokeninline]) eq "A") { 
				if($increasea_exist && ($increasereactionlabel eq $reactionlabel[$countreactions])) {
				    print "\nreaction: $paragraph";
				    print "\nchanged A from ",$tokens[$tokeninline+2];
				    $tokens[$tokeninline+2]=sprintf "%.6g",$tokens[$tokeninline+2]*$increaseafactor;    
				    print " to ",$tokens[$tokeninline+2];
				    $correcteda=1;
				}
				$freqfactor=$tokens[$tokeninline+2]; 
			    }
			    if(uc($tokens[$tokeninline]) eq "N") { $tempcoeff=$tokens[$tokeninline+2]; }
			    if(uc($tokens[$tokeninline]) eq "E") { $actenergy=$tokens[$tokeninline+2]; }
			    if(uc($tokens[$tokeninline]) eq "AI") { 
				if($increasea_exist && ($increasereactionlabel eq $reactionlabel[$countreactions])) {
				    print "\nchanged AI from ",$tokens[$tokeninline+2];
				    $tokens[$tokeninline+2]=sprintf "%.6g", $tokens[$tokeninline+2]*$increaseafactor;    
				    print "to ",$tokens[$tokeninline+2];
				}
				$freqfactor=$tokens[$tokeninline+2]; 
			    }
			    if(uc($tokens[$tokeninline]) eq "NI") { $tempcoeff=$tokens[$tokeninline+2]; }
			    if(uc($tokens[$tokeninline]) eq "EI") { $actenergy=$tokens[$tokeninline+2]; }
			}
		    }
		    if ($freqfactor<-1e6) { print "\nfreqfactor not found in reaction \n$parameterpart"; }
		    if ($tempcoeff<-1e6) { print "\ntempcoeff not found in reaction \n$parameterpart"; }
		    if ($actenergy<-1e6) { print "\nactenergy not found in reaction \n$parameterpart"; }
		    #evaluate reaction rate for $tempss
		    #print "\nfreqfactor: ",$freqfactor;
		    #print "\ntempcoeff: ",$tempcoeff;
		    #print "\nactenergy: ",$actenergy;
		    if ($correcteda==1) {
			$paragraph=$reactionpart." { ";
			for (my $tokeninline=0;$tokeninline<=$#tokens;$tokeninline++) {
			    $paragraph=$paragraph.$tokens[$tokeninline]." ";
			}
			$paragraph=$paragraph."}".$partafterparameter;
		    }
		    $ratecoeff[$countreactions]=$freqfactor*$tempss**$tempcoeff*exp(-$actenergy/0.0083143/$tempss);
#printf("reaction: %i %s ratecoeff: %f\n",countreactions-1,reactionlabel[countreactions-1],ratecoeff[countreactions-1]);
		    $countreactions++;  #finally we know that we have one valid reaction more
		}
	    }
	    
	}
	
	$statusofparsing=0; 
    }
    if (($statusofparsing==1) && ($currentline !~ /\}/))  {

    } else {
	printf FILEMECHOUT "$paragraph";
    }
    #printf("%s",currentline);
} # for all lines in mechfile
close(FILEMECHIN);


if ($specsteadystates_exist) { # now put the steady states in file
    if (!open FILESTEADYSTATESOUT, ">$fsteadystatesnameout") { 
	print "\ncannot open steady states output file $fsteadystatesnameout\n";
	die;
    }
    
    print "\nSteady States:";
    # sort the reactions that belong to one species by (reactionrate*stoichcoeff) so that the most
    # negative is left. That is the one with the fastest consumption.
    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	if ($numberreactionforss[$iss]) {
	    my $switchbufferdouble;
	    for(my $sortrun=0;$sortrun<$numberreactionforss[$iss];$sortrun++) {
		for($ireac=0;$ireac<$numberreactionforss[$iss]-1;$ireac++) {
		    if ($stoichcoeff[$iss][$ireac]*$ratecoeff[$reactionforss[$iss][$ireac]]>$stoichcoeff[$iss][$ireac+1]*$ratecoeff[$reactionforss[$iss][$ireac+1]]) {
			#switch them
			$switchbufferint=$stoichcoeff[$iss][$ireac];
			$stoichcoeff[$iss][$ireac]=$stoichcoeff[$iss][$ireac+1];
			$stoichcoeff[$iss][$ireac+1]=$switchbufferint;
			$switchbufferint=$reactionforss[$iss][$ireac];
			$reactionforss[$iss][$ireac]=$reactionforss[$iss][$ireac+1];
			$reactionforss[$iss][$ireac+1]=$switchbufferint;
		    }
		}
	    }
	    #fprintf(FILESTEADYSTATESOUT,"%s\tw%s\n",specsteadystatesnames[iss],reactionlabel[reactionforss[iss][0]]);
	    print "\n'$specsteadystatesnames[$iss]'";
	    printf "chosen: %s %g",$reactionlabel[$reactionforss[$iss][0]],$ratecoeff[$reactionforss[$iss][$0]];
	    for ($ireac=0;$ireac<$numberreactionforss[$iss];$ireac++) {
		print "\n",$reactionlabel[$reactionforss[$iss][$ireac]];
		#printf("%g \n",ratecoeff[reactionforss[iss][ireac]]*stoichcoeff[iss][ireac]);
	    }
	    my $newiss=$iss-$steadystatesnotfound;
	    for($ireac=0;$ireac<$numberreactionforss[$iss];$ireac++) {
		$reactionforss[$newiss][$ireac]=$reactionforss[$iss][$ireac];
		$stoichcoeff[$newiss][$ireac]=$stoichcoeff[$iss][$ireac];
	    }
	    $specsteadystatesnames[$newiss]=$specsteadystatesnames[$iss];
	    $numberreactionforss[$newiss]=$numberreactionforss[$iss];
	} else { #no reactions found for that species 
	    printf "\n%s was not a reactant in the mechanism",$specsteadystatesnames[$iss];
	    $steadystatesnotfound++;
	}
    }
    $numberspecsteadystates -= $steadystatesnotfound;
    # now comes the big, big routine for sorting the steady state species according to their
    # dependency of each other. If they are not coupled this should be a fail-safe method
    # If they are coupled there might be a chance to uncouple them.
    # First we determine for each species how often its reaction is also used by other ss specs
    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	$numberdepending[$iss]=0;
	for ($isstocheck=$iss+1;$isstocheck<$numberspecsteadystates;$isstocheck++) {
	    for($ireac=0;$ireac<$numberreactionforss[$isstocheck];$ireac++) {
		if($reactionforss[$isstocheck][$ireac]==$reactionforss[$iss][0]) { $numberdepending[$iss]++; }
	    }
	}
	printf "\ninit. depend. for species %s: %i",$specsteadystatesnames[$iss],$numberdepending[$iss];
    }
    
    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	$numberserving[$iss]=0;
	for ($isstocheck=0;$isstocheck<$iss;$isstocheck++) {
	    for($ireac=0;$ireac<$numberreactionforss[$iss];$ireac++) {
		#printf("%s,%s\n",reactionlabel[reactionforss[isstocheck][0]],reactionlabel[reactionforss[iss][ireac]]);
		if($reactionforss[$isstocheck][0]==$reactionforss[$iss][$ireac]) { $numberserving[$iss]++; }
	    }
	}
	printf "\ninit. serving. for species %s: %i",$specsteadystatesnames[$iss],$numberserving[$iss];
    }
    
    # now we perform a bubble sort on the steady state species. The species are compared by 
    # dependencies. Two species are switched if the dependency of the second species is
    # lower than the one of the first species. A relative dependency of the two species to
    # each other must be subtracted before applying the criteria.
    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	#printf("\n");
	for ($isstocheck=0;$isstocheck<$numberspecsteadystates-1;$isstocheck++) {
	    #printf("spec: %s dep: %i\n",specsteadystatesnames[isstocheck+1],numberdepending[isstocheck+1]);
	    # first we check the dependency of each other
	    $dep1on2=0;
	    for($ireac=0;$ireac<$numberreactionforss[$isstocheck+1];$ireac++) {
		#printf("%s %s\n",reactionlabel[reactionforss[isstocheck][0]],reactionlabel[reactionforss[isstocheck+1][ireac]]);
		if($reactionforss[$isstocheck][0]==$reactionforss[$isstocheck+1][$ireac]) { $dep1on2=1; }
	    }
	    $dep2on1=0;
	    for($ireac=0;$ireac<$numberreactionforss[$isstocheck];$ireac++) {
		#printf("%s %s\n",reactionlabel[reactionforss[isstocheck+1][0]],reactionlabel[reactionforss[isstocheck][ireac]]);
		if($reactionforss[$isstocheck+1][0]==$reactionforss[$isstocheck][$ireac]) { $dep2on1=1; }
	    }
	    #printf("dep1on2, dep2on1  %i  %i\n",dep1on2,dep2on1);
	    my $doswitch=-1;
	    if(!$dep1on2 && $dep2on1) {
		if ($numberdepending[$isstocheck]-$numberserving[$isstocheck]>$numberdepending[$isstocheck+1]-$numberserving[$isstocheck+1]) { $doswitch=1; } else { $doswitch=0; }
	    }
	    if(!$dep1on2 && $dep2on1) {
		# leave it for now (maybe I find a better criteria later)
		$doswitch=0;
	    }
	    if($dep1on2 && !$dep2on1) {
		{ $doswitch=1; }
	    }
	    if($dep1on2 && $dep2on1) {
		# now starts the real mess. lets just leave it for now
		$doswitch=0;
		printf "\nSteady State Species %s and %s depend on each other",$specsteadystatesnames[$isstocheck],$specsteadystatesnames[$isstocheck+1];
	    }
	    if ($doswitch) {
		for($ireac=0;$ireac<max($numberreactionforss[$isstocheck],$numberreactionforss[$isstocheck+1]);$ireac++) {
		    $switchbufferint=$reactionforss[$isstocheck+1][$ireac];
		    $reactionforss[$isstocheck+1][$ireac]=$reactionforss[$isstocheck][$ireac];
		    $reactionforss[$isstocheck][$ireac]=$switchbufferint;
		    $switchbufferint=$stoichcoeff[$isstocheck+1][$ireac];
		    $stoichcoeff[$isstocheck+1][$ireac]=$stoichcoeff[$isstocheck][$ireac];
		    $stoichcoeff[$isstocheck][$ireac]=$switchbufferint;
		}
		$switchbufferchar=$specsteadystatesnames[$isstocheck+1];
		$specsteadystatesnames[$isstocheck+1]=$specsteadystatesnames[$isstocheck];
		$specsteadystatesnames[$isstocheck]=$switchbufferchar;
		$switchbufferint=$numberreactionforss[$isstocheck+1];
		$numberreactionforss[$isstocheck+1]=$numberreactionforss[$isstocheck];
		$numberreactionforss[$isstocheck]=$switchbufferint;
		$switchbufferint=$numberdepending[$isstocheck+1];
		$numberdepending[$isstocheck+1]=$numberdepending[$isstocheck];
		$numberdepending[$isstocheck]=$switchbufferint;
		$switchbufferint=$numberserving[$isstocheck+1];
		$numberserving[$isstocheck+1]=$numberserving[$isstocheck];
		$numberserving[$isstocheck]=$switchbufferint;
		# we have to update the dependencies here
		if($dep1on2) {
		    if($numberdepending[$isstocheck+1]>0) { $numberdepending[$isstocheck+1]--; }
		    if($numberserving[$isstocheck]>0) { $numberserving[$isstocheck]--; }
		}
		if($dep2on1) {
		    $numberdepending[$isstocheck]++;
		    $numberserving[$isstocheck+1]++;
		}
	    }
	}
	printf "\ndependencies for species %s: %i",$specsteadystatesnames[$iss],$numberdepending[$iss];
    }
    
    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
	printf FILESTEADYSTATESOUT "%s\tw%s\n",$specsteadystatesnames[$iss],$reactionlabel[$reactionforss[$iss][0]];
	printf "\n%s",$specsteadystatesnames[$iss];
    }
}
#now read and write the mechanism again, this time with the right order of
#steady state species
close(FILEMECHOUT);

if (!open FILEMECHIN, $fintmechnameout) { 
    printf "\ncannot open mechanism file %s \n", $mechfile;
    die;
}
if (!open FILEMECHOUT, ">$fmechnameout") { 
    print "\ncannot open output file '$fmechnameout'\n";
    die;
}

for (my $linenr=0; defined ($currentline = <FILEMECHIN>); $linenr++) {
    if($specsteadystates_exist) {
	my $issteadystates=0; 
	if ($currentline=~"^let steady states be") { $issteadystates=1; }
	if ($issteadystates) { # start line for steady states
	    printf FILEMECHOUT "let steady states be \n";
	    for ($iss=0;$iss<$numberspecsteadystates;$iss++) {
		$currentline = <FILEMECHIN>; #basically dummy reading
		printf FILEMECHOUT "%s \n",$specsteadystatesnames[$iss];
	    }
	    $currentline = <FILEMECHIN>;
	    for (my $issnotfound=0;$issnotfound<$steadystatesnotfound;$issnotfound++) { $currentline = <FILEMECHIN>; }
	    if ($currentline=~/\./) {
		$currentline=".\n";
	    } else { # no end found
		print "\n. not found at end of steady state species\n";
		die;
	    }
	}
    } 
    printf FILEMECHOUT $currentline;
    #printf("%s",currentline);
} # for all lines in mechfile
close(FILEMECHIN);
close(FILEMECHOUT);

# main

print "\n";
print STDERR "\n";















