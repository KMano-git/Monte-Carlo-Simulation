#!/usr/bin/perl

##  data ==> make subroutine to print data
##  &mksub("dt_git,"git_mesg.f");
 
sub mksub{
    local ($dat, $sub) = @_;


#Header KH 20121009
    open(FSUB,"> $sub");
    open(FHED,"$sub.head");
    while(<FHED>){
	$line = $_;
	chomp($line);
	next if( length($line) <= 0 ) ;
	printf FSUB "$line\n";
    }
    close  FHED;

##Data read write
    open(FDAT,"$dat");
    while(<FDAT>){
	$line = $_;
	chomp($line);
	next if( length($line) <= 0 ) ;
	printf FSUB "write(n6,61)  \"$line\"\n";
    }
    close  FDAT;

##Closer
    printf FSUB "write(n6,61)  \"   \"\n";
    printf FSUB "   \n";
    printf FSUB "return\n";
    printf FSUB "end\n";
    close  FSUB;
}
1;
