#!/usr/bin/perl

##  current directory   prj/src/wksvn

push( @INC, "./wkgit" );
require "mksub.pl";

$lod = $ARGV[0];

chdir "./wkgit";
$wkd = `pwd`;
chomp($wkd);
$prj = "../../";   ## <== ~/PRJ/sonicA
$prj = `cd $prj && pwd`;
chomp($prj);

##print "CHK  pwd = $wkd\n";
##print "CHK  prj = $prj\n";
##$lod = "/home/g9/a003859/PRJ/sonicA/exe/slimCS10/LOD/ld_name";

$dat = "$wkd/dt_git";
$sub = "git_mesg.f";
$tim = `date`;

system("rm -rf $dat");
system("rm -rf *.o");

system("echo '   ' > $dat");
system("echo '=== Load module ===' >>  $dat");
system("echo 'prj : $prj'  >> $dat");
system("echo 'lod : $lod'  >> $dat");
system("echo 'dat : $tim'  >> $dat");
system("echo '   '  >> $dat");

system("echo '=== git branch  ===' >> $dat");
system("cd $prj && git branch  >> $dat");
system("echo '   '  >> $dat");
system("echo '=== git log  ===' >> $dat");
system("cd $prj && git log -1  >> $dat");
system("echo '   '  >> $dat");
system("echo '=== git status  ===' >> $dat");
system("cd $prj && git status -s  >> $dat");
system("echo '   '  >> $dat");
#system("echo '=== git diff  ===' >> $dat");
#system("cd $prj && git diff  >> $dat");
#system("echo '   '  >> $dat");
system("cat $dat");

&mksub($dat, $sub);

