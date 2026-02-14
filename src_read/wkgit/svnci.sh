#!/bin/csh

set langb = `printenv LANG`
set topprj = $PRJ_sol2d   # ~/PRJ/sonicA
set log = "lgmsg"
set logpth = "$topprj/$log"
set logtmp = "@lgmsg"
set inf = "@info"
set drnm = $topprj
set pwd = `pwd`
set langb = `printenv LANG`
setenv LANG C

## make lgmsg
cd $drnm
svn log -r HEAD >! @1
set new = `grep ^r @1 | cut -f1 -d' ' | sed -n 's/r//p'`
##echo "last revisin = [$new]"

svn log -r $new\:1 >! @1
grep ^20 @1 >! $log

## revision number
cd $drnm
set url = `svn info | grep ^URL | sed 's/.* //'`
set Myrev = `svn info      | grep ^Revision | cut -f2 -d: | sed 's/ *//'`
set Cmrev = `svn info $url | grep ^Revision | cut -f2 -d: | sed 's/ *//'`
echo "Revision : Work = $Myrev  Repository = $Cmrev"

##svn st -u
##goto LB100

echo "  "
echo "topprj = $topprj"
echo "wk-dir = $pwd"

echo "  "
##echo -n ">> commit for dir (top:<rtn>/w:wk/src/q) ==> "
##
##set ans = $<
##if( "$ans" == "" || $ans == "t" ) then
##  set drnm = $topprj
##else if( "$ans" == "w" ) then  ## <===
##  set drnm = $pwd
##else if( "$ans" == "q" ) then
##  echo "Command was cancelled."
##  exit
##else
##  set drnm = "$topprj/$ans"
##  if(! -d $drnm ) then
##    echo "No found dir : $drnm"
##    exit
##  endif
##endif

set drnm = $pwd

LB100:
cd   $drnm
echo "commit for dir : "`pwd`
if( $drnm != $topprj ) then
  svn st -u
endif
##echo -n ">> enter execution ? (<rtn>/q:quit) ==> "
##set ans = $<
##if( "$ans" == "q" ) then
##  echo "Command was cancelled."
##  exit
##endif

## revision
cd $topprj
svn info > $inf

set no = $Cmrev
rm -f @info
@ no++

set nomx = 4
set mj = `echo $no | wc -c`
@ mj--
@ i = $mj
set rev = "r$no"
while( $i < $nomx )
 set rev = "$rev "
 @ i++
end
set rev = "$rev"': '
set dat = `date +"20%y %m/%d"`


cp  $log  $logtmp
echo "$dat  $rev" > $log
cat $logtmp >> $log
rm -f $logtmp

mv lgmsg_org XXlgmsg
edt $log   # emacs -nw $log
mv XXlgmsg lgmsg_org

set comd = "svn commit -m "
echo -n ">> svn commit O.K. (<rtn>:y/n) ==> "
set ans = $<
if( "$ans" == "n" ) then
  echo "command was cancelled."
  exit
endif

set mesg = `head -1 $log`
set mesg = `echo $mesg | sed 's/\n//'`
set mj = `echo "$mesg" | wc -c`
@ mj--
set mesg = `echo $mesg | cut -c1-$mj`

echo "COMD  = [""$comd "'"'"$mesg"'"'"]"
cd $drnm
$comd  "$mesg"

## revision number
cd $topprj
set url = `svn info | grep ^URL | sed 's/.* //'`
set Myrev = `svn info      | grep ^Revision | cut -f2 -d: | sed 's/ *//'`
set Cmrev = `svn info $url | grep ^Revision | cut -f2 -d: | sed 's/ *//'`
echo "URL : $url"
echo "Revision : Work = $Myrev  Repository = $Cmrev"

echo "  "
echo "Do not forget [svn update] command."
echo "  "
