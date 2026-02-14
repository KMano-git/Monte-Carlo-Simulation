#!/usr/bin/perl

$home = $ENV{"HOME"};
print "home = $home\n";

$drfm = "$home/sol2d.v2#0/src/pstNEW";
print "from Dir : $drfm\n";

chomp( $drto = `pwd` );
print "  to Dir : $drto\n";

chdir $drfm;
@temp = `ls *.f`;
@tbfl = ();
foreach $fnam (@temp){
    chomp($fnam);
    push( @tbfl, $fnam );
}

chdir $drto;
@tbex = ();
foreach $fnam (@tbfl){
  if( -f $fnam ){
     print "Find $fnam\n"; 
     push( @tbex, $fnam );
  }
  else{
    print "copy  $fnam\n";
###    system(" cp $drfm/$fnam ."); 
  }
} 

print "\n";
print "\n";
$dffl = "\@diff";
foreach $fnam (@tbex){
    system("sdiff -l -w 180  $fnam  $drfm/$fnam > $dffl");
  }
print "diff files with the same name  see $dffl\n";

