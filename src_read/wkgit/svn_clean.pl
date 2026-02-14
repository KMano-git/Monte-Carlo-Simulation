#!/usr/bin/perl

## delete dir. of .svn
## delete Makefile
## shift  PerlM/Myst

## current dir. :  sonicA

$prj = "$ENV{'PRJ_sol2d'}";

$sft = "PerlM/Myst"
chdir $prj;




##@del = ();
##push( @del, ".svn" );
##push( @del, "Makefile" );
##push( @del, "Makeopt" );
##push( @del, "*.o" );
##push( @del, "Myst" );


$flst = "\@flst";

chdir $prj;
system("find . -name $del > $flst");
##system("cat $flst");


open(FLST,"$flst");
$ii = 0;
while(<FLST>){
  $line = $_;
  chomp($line);
  $tab[$ii] = "rm -rf $line";
  if( $ii < 6 || $ii%10 == 0 ){
      printf "   %3d:  %s\n", $ii, $tab[$ii];
  }
  $ii = $ii + 1;
}
close FLST;

