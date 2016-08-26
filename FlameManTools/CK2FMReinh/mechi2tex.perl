#!/usr/bin/perl -w 
(-e  $ARGV[0]) or print " $ARGV[0] not found \n
usage will see\n ";
$fnameroot = $ARGV[0];
$fnameroot =~  s/(\S+\.)(\w+)/$1/g;
print "fnameroot  $fnameroot  \n";
#    use File::Copy;
#    copy("$ARGV[0]","tmp.$$");


open Infile0, "$ARGV[0]" or die ;
$ilines = 0;
while (defined ($line[$ilines] = <Infile0>)) {chomp($line[$ilines++]) };
close(Infile0);
#$nlines = $ilines-1;

$ilines = 0;
while  ($line[$ilines++] !~ /ELEMENTS/) { };
$ilines+=2;
$i = 0;
while  ($line[$ilines++] =~ /\d+\./) { 
    @token = split  /\s+/, $line[$ilines-1];
    $element[$i] = $token[2];
    print "$element[$i]\n";
    $nelement = $i;
    $i++;
};
    print "$element[3]\n";

while  ($line[$ilines++] !~ /REACTIONS CONSIDERED/) { };
$ilines--;
#print  "$ilines : \t $line[$ilines]\n";
# found beginning of reaction data

$ireac = 1;
while  ($line[$ilines++] !~ /=/) { };
$ilines--;
#print  "$ilines : \t $line[$ilines]\n";
while ($line[$ilines] =~ /=/) { # all reactions must have =, some will have =>
    @token0 = split /\. /, $line[$ilines],2;
    #$No[$ireac] =  $token0[0];
    @token = split /\s+/, $token0[1];
    $reac[$ireac] = $token[0];
    $forwardonly[$ireac] = 0;
    if (($reac[$ireac] =~ /=>/) && ($reac[$ireac] !~ /</)) { 
	$forwardonly[$ireac] = 1;
    }  #check for forward only
    
    @sides =  split /=/, $token[0];
    $sides[$ireac][0] = $sides[0];
    $sides[$ireac][1] = $sides[1];

#    print "left:  $sides[$ireac][0]\t right:  $sides[$ireac][1]\n";
    $kA[$ireac] =  $token[1];	# 
    $kb[$ireac] =  $token[2];	# 
    $kE[$ireac] =  $token[3];	# 

    $low[$ireac] = 0;
    $centering[$ireac] = 0;
    $enhanced[$ireac] = 0;

    while ($line[$ilines+1] !~ /=/ and $line[$ilines+1] !~ /NOTE/) { #check simultaneously for centering,3rd body,and low pressure

	if ($line[$ilines+1] =~ /Low/) {
	    #$ilines++;
	    @token = split /\s+/, $line[$ilines+1];
	    $low[$ireac] = 1;
	    $kAlow[$ireac] =  $token[4];
	    $kblow[$ireac] =  $token[5];
	    $kElow[$ireac] =  $token[6];
	} 
	
        #------------- check for reverse reaction
	
	if ($line[$ilines+1] =~ /Reverse Arrhenius coe/) {
	    #$ilines++;
	    $reverse[$ireac] = 1;
	    @token = split /\s+/, $line[$ilines+1];
	    $kArev[$ireac] =  $token[4];
            #	print " token[4] $token[4]\n";
	    $kbrev[$ireac] =  $token[5];
	    $kErev[$ireac] =  $token[6];
            #	print " token[5] $token[5]\n";
	}

        #------------- i'm sure i will work on the centering 
	
	if ($line[$ilines+1] =~ /T&H/) {
	    $TandH[$ireac] = 1;
            #	$ilines++;
	    @token = split /\s+/, $line[$ilines+1];
	    $centering[$ireac] = $#token-2;
            #	print " TH[$ireac] : $centering[$ireac]\n";
	    for ($iTH = 0; $iTH < $centering[$ireac] ;$iTH++ ) {$centering[$ireac][$iTH] =  $token[$iTH+3]};
	} 
	if ($line[$ilines+1] =~ /TROE centering:/) {
	    $TROE[$ireac] = 1;
            #	$ilines++;
	    @token = split /\s+/, $line[$ilines+1];
	    $centering[$ireac] = $#token-2;
            #	print " TH[$ireac] : $centering[$ireac]\n";
	    for ($iTH = 0; $iTH < $centering[$ireac] ;$iTH++ ) {$centering[$ireac][$iTH] =  $token[$iTH+3]};
	} 
	if ($line[$ilines+1] =~ /SRI centering:/) {
	    $SRI[$ireac] = 1;
            #	$ilines++;
	    @token = split /\s+/, $line[$ilines+1];
	    $centering[$ireac] = $#token-2;
            #	print " TH[$ireac] : $centering[$ireac]\n";
	    for ($iTH = 0; $iTH < $centering[$ireac] ;$iTH++ ) {$centering[$ireac][$iTH] =  $token[$iTH+3]};
	} 
        #------------- i'm sure i will work on the centering 
	
	if ($line[$ilines+1] =~ /Reverse Arrh/) {
        #	$ilines++;
	    @token = split /\s+/, $line[$ilines+1];
	    $reverse[$ireac] = $#token-3;
        #	print " reverse[$ireac] : $reverse[$ireac]\n";
	    for ($i= 0; $i < $reverse[$ireac] ;$i++ ) {$reverse[$ireac][$i] =  $token[$i+3]};
	} 
	$incrilines=1;  
	$ienhanced = 0;
	while ($line[$ilines+1] =~ /Enhanced/) {
	    $incrilines=0;  
	    $ilines++;
	    @token = split /\s+/, $line[$ilines];
	    $enhanced[$ireac] += 1;	# 
	    $enhSpec[$ireac][$ienhanced] =  $token[1];
	    $enhVal[$ireac][$ienhanced] =  $token[4]; # 
	    $ienhanced++;		 
	} ;
	while ($line[$ilines+1] =~ /duplicate/) {
	    $incrilines=0;  
	    $ilines++;
	}
	$ilines+=$incrilines;
    } #end of while next reaction not found 
    
    print "$ireac \t$reac[$ireac]\t";
    if($low[$ireac]) {print "k_{\\inf}"};
    print "\t$kA[$ireac]\t$kb[$ireac]\t$kE[$ireac]\n";
    if ($low[$ireac]) {print "\t\t\t\tk_0\t$kAlow[$ireac]\t$kblow[$ireac]\t$kElow[$ireac]\n"};	 
    for ($i = 0; $i < $centering[$ireac]; $i++) {
	printf "\ta_$i = %2.4f", $centering[$ireac][$i];
	$troefactors[$ireac][$i]=$centering[$ireac][$i];
    }				# 
    if ( $centering[$ireac]) {print "\n"};
    for ($ienhanced =0; $ienhanced < $enhanced[$ireac];$ienhanced++){
	if ($ienhanced ==0) {print "\tthird-body efficiencies "};
	printf " / %s %2.1f",$enhSpec[$ireac][$ienhanced],$enhVal[$ireac][$ienhanced];
    }
    if( $enhanced[$ireac]){print " /\n"};
    $ilines++;
    $ireac++;			# 
}				# 

# -----------LaTeX output-----------------------------------------------
if ($#ARGV >0) {
    $latexname = $fnameroot ;
    $latexname .= "tex";
    print "printing LaTeX table into file $latexname\n";

    $nreac = $ireac-1;
    open Outfile, ">$latexname" or die ;  

    print Outfile  "\\documentclass[11pt]{article}\n";
    print Outfile  "\\pagestyle{empty}\n";
    print Outfile  "\\textheight 9.00in\n";
    print Outfile  "\\textwidth 6.5in\n";
    print Outfile  "\\voffset -3.0cm\n";
    print Outfile  "\\oddsidemargin 0.00in\n";
    print Outfile  "\n";

    print Outfile  "\%\\hline\n";
    print Outfile  "\%\\end{tabular}\n";
    print Outfile  "\%\\end{table}\n";
    print Outfile  "\%\\begin{table}\n";
    print Outfile  "\%\\begin{tabular}{|r l c l l l|} \\hline \n";
    print Outfile  "\%No& Reaction &   & \$A_n\$& \$b_n\$ & \$E_n\$\\\\ \\hline \n";
    print Outfile  "\n";

    print Outfile  "\\begin{document}\n";
    print Outfile  "\\begin{table}[hbt]\n";
    print Outfile  "\\begin{tabular}{|r l c l l l|} \\hline \n";
    print Outfile  "No& Reaction &   & \$A_n\$& \$b_n\$ & \$E_n\$\\\\ \\hline \n";
    
    for ($ireac=1;$ireac <= $nreac; $ireac++)  { 
#    print Outfile "$ireac &\t$reac[$ireac]&\t";
	$_ = $reac[$ireac];
#	$reac[$ireac] =~ s/([A-Z])(\d+)/$1\$_{$2}\$/g;
	s/([A-Z])(\d+)/$1\$_{$2}\$/g;
#	$skipnext=0;
#	if($sides[$ireac][0]==$sides[$ireac+1][1] && $sides[$ireac][1]==$sides[$ireac+1][0]) {
#	    $skipnext=1;
#	    if ($kA[$ireac]<1.0)
        s/=/\$\\rightleftharpoons\$ /;
	$texreac[$ireac] = $_;
	print Outfile "$ireac &\t$texreac[$ireac]&\t";
        
	if($low[$ireac]) {print Outfile "\$k_{\\infty}\$"};
	print Outfile "&\t$kA[$ireac]&\t$kb[$ireac]&\t$kE[$ireac]\\\\ \n";
	if ($low[$ireac]) {
	    printf Outfile "\t&\t&\t\$k_0\$&\t%.2E&%5.2f&\t%7.1f\\\\ \n",$kAlow[$ireac],$kblow[$ireac],$kElow[$ireac]
	    };			 

	if ( $centering[$ireac]) {	 
	    print Outfile "\t&\\multicolumn{5}{l|}{";
	    for ($i = 0; $i < $centering[$ireac]; $i++) {
		printf Outfile "\t\$a_$i\$ = %2.4G", $centering[$ireac][$i];  
	    }				
	    print Outfile "}\\\\ \n";
	} 
	if( $enhanced[$ireac]){
	    print Outfile "\t&\\multicolumn{5}{l|}{";
	    for ($ienhanced =0; $ienhanced < $enhanced[$ireac];$ienhanced++){
		if ($ienhanced ==0) {print Outfile "third-body efficiencies "};
		$_ = $enhSpec[$ireac][$ienhanced];
		s/([A-Z])(\d+)/$1\$_{$2}\$/g;
		$texenhSpec[$ireac][$ienhanced] = $_;
		printf  Outfile " / %s %2.1f",$texenhSpec[$ireac][$ienhanced],$enhVal[$ireac][$ienhanced];
	    }
	    print Outfile " /}\\\\ \n";
	};		       
    }
    print Outfile "\\hline\n";
    print Outfile  "\\end{tabular} \n";
    print Outfile  "\\end{table} \n";
    print Outfile  "\\end{document}\n";
};




# -----------ScanMan output-----------------------------------------------
if ($#ARGV >1) {
    $ScanManname = $fnameroot ;
    $ScanManname .= "mech";
    print "printing ScanManMechanism into file $ScanManname\n";

    $nreac = $ireac-1;
    open Outfile, ">$ScanManname" or die ;  

    print Outfile  "let allowed atoms be ";
    for ($i=0;$i<$nelement;$i++) {
	print Outfile  " $element[$i],";
    }
    print Outfile  " $element[$nelement].\n";
    my $flag = 0;
    foreach (@element){if ($_ =~ /ar/i){$flag = 1}}
    if ($flag){
	print Outfile  "let additional species be N2, AR.\n";
    }else{
	print Outfile  "let additional species be N2.\n";
    }
    print Outfile  "let temperature exponent be n_k.\n";
    print Outfile  "let order of reaction be n.\n";
    print Outfile  "let units for A be [ cm^(3(n-1)) / ( s * mole^(n-1) * K^n_k ) ].\n";
    print Outfile  "let units for E be [ kJ / mole ].\n";
    print Outfile  "\n";

    # init vars
    $enhancedix = 0;
    for ($ireac=1;$ireac <= $nreac; $ireac++)  { 
#    print Outfile "$ireac &\t$reac[$ireac]&\t";
	$reac[$ireac] =~ s/\(\+[mM]\)/\+M/g;	# Here (+m) is changed to +M 
	$reac[$ireac] =~ s/\+m/\+M/g;	#  Here +m is changed to +M
	$reac[$ireac] =~ s/\+/ \+ /g;	# 
	$reac[$ireac] =~ s/<?=+>?/ -> /g;	# 	
	@sides = split /->/, $reac[$ireac];     #
	if ( $enhanced[$ireac]){ # 
	    $reacenhix[$ireac] = ++$enhancedix;
	    $reac[$ireac] =~ s/M/M$reacenhix[$ireac]/g;
	};
	if ($forwardonly[$ireac]) { print Outfile "${ireac}: \t$reac[$ireac]\t";}
	else { print Outfile "${ireac}f: \t$reac[$ireac]\t";}
        $kE[$ireac]*=0.004184;   # in kJ/mole
        if ($kElow[$ireac]) {$kElow[$ireac]*=0.004184}; # in kJ/mole
	if($low[$ireac]) {	# 
	    printf Outfile "\{ ai = %.3E ni = %5.3f Ei = %7.3f \n",$kA[$ireac],$kb[$ireac],$kE[$ireac];
	    printf Outfile "\t\t\t a = %.3E n = %5.3f E = %7.3f \n",$kAlow[$ireac],$kblow[$ireac],$kElow[$ireac];
	}			# 
	else {			# 
	    printf Outfile "\{ a = %.3E n = %5.3f E = %7.3f ",$kA[$ireac],$kb[$ireac],$kE[$ireac];
	};			# 

# centering: 

    	if ($TROE[$ireac]) {
	    printf Outfile "\n\tfca = %.3E fcta = %.3E fcb = %.3E fctb = %.3E",
	    1-$troefactors[$ireac][0],$troefactors[$ireac][1],$troefactors[$ireac][0],$troefactors[$ireac][2];
	    if ($troefactors[$ireac][3]) {      # (regular) 4 coef TROE form
		printf Outfile " fcc = 1 fctc = %.3E ",$troefactors[$ireac][3];
	    }
	    printf Outfile "\n";
	}
#    	if ($TROE[$ireac]) {
#	    printf Outfile "\tfc = (1-%s) * exp(-T/%s) + %s * exp(-T/%s)",
#	    $troefactors[$ireac][0],$troefactors[$ireac][1],$troefactors[$ireac][0],$troefactors[$ireac][2];
#	    if ($troefactors[$ireac][3]) {      # (regular) 4 coef TROE form
#		printf Outfile " + exp(-%s/T)",$troefactors[$ireac][3];
#	    }
#	    printf Outfile ";\n";
#	}
	elsif ($SRI[$ireac]) {
	     printf Outfile "\tfc = 1.0;\n\#warning: SRI centering changed\n";
	 }
	     
	elsif  ($TandH[$ireac]) { # 
	    printf Outfile "\tfc = %s + (%s *T);}\n",$centering[$ireac][0],$centering[$ireac][1];
	}
	elsif  ($centering[$ireac] == 1) {
	    printf Outfile "\tfc = %s ;\n",$centering[$ireac][0];
	}
	elsif ($low[$ireac]){
	    printf Outfile "\tfcc = 1.0 fctc = 0.0\n";  # default lindemann form 
	}
	print Outfile "} \n";


	if ($reverse[$ireac]) {
	  $backwardreaction = "${ireac}b: \t$sides[1] -> $sides[0]\t";
	if ( $enhanced[$ireac]){ # 
	    $backwardreaction =~ s/M/M$reacenhix[$ireac]/g;
	};
	    print Outfile "$backwardreaction";
	    $kErev[$ireac]*=0.004184;   # in kJ/mole
    	if ($TROE[$ireac]) {
	    printf Outfile "\{ ai = %.3E ni = %5.3f Ei = %7.3f \}\n",$kArev[$ireac],$kbrev[$ireac],$kErev[$ireac];
	  } else {
	    printf Outfile "\{ a = %.3E n = %5.3f E = %7.3f \}\n",$kArev[$ireac],$kbrev[$ireac],$kErev[$ireac];
	  }
	}
    }				# 
 

   # print third body efficiencies
    print Outfile "let M = 1.0[OTHER].\n"; 
    for ($ireac=1;$ireac <= $nreac; $ireac++)  { 
	if( $enhanced[$ireac]){	# 
	    for ($ienhanced =0; $ienhanced < $enhanced[$ireac];$ienhanced++){
		if ($ienhanced ==0) {print Outfile "let M$reacenhix[$ireac] = "};
		printf  Outfile "%2.2f [%s] + ",$enhVal[$ireac][$ienhanced],$enhSpec[$ireac][$ienhanced];
	    }			# 
	    print Outfile "1.0 [OTHER]. \n";
	}
    }
#    for ($j=1;$j <= 100; $j++)  {
#      for ($k=0;$k <= 3; $k++)  {
#      print "centering no. $j $k = $troefactors[$j][$k] \n"; 
#}
#} 

};




