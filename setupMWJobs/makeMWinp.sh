#!/bin/bash

cat inp.list | while read temp pres spec edown ntrial tfin einit moli
   do
     echo  "$name"
# Make Temperature and subsequent pressure folders, make all corresponding .dat files and copy
# DensData folder into folders
    mkdir "$temp"K
    cd "$temp"K
    mkdir $pres
    cd "$pres"
    rm "$pres".dat
    cat ../../a.dat >> "$pres".dat
    echo "$temp $temp">> "$pres".dat
    echo "1 !no. pressures" >> "$pres".dat
    echo  "$pres" >> "$pres".dat
    cat ../../b.dat >> "$pres".dat
    echo "1 4.46 390 1 $edown 0.0 0.0 0.0 0.0 0.0 0.0 0.0 !exp down model" >> "$pres".dat
    echo "'LJ'"  >> "$pres".dat
    cat ../../c.dat >> "$pres".dat
    echo "$ntrial  'TIME' $tfin 'THERMAL'       $moli       3       $einit !24.66  !Ntrials, Tspec, Tread, KEYTEMP, Molinit, IR(negl" >> "$pres".dat
    echo " " >> "$pres".dat
    cp -r ../../DensData ./
#
# Make and copy qsub scripts into all folders
    rm "MW$temp.$pres.$spec"
    cat ../../sh1 >>   "MW$temp.$pres.$spec"
    echo "cd `pwd`" >> "MW$temp.$pres.$spec"
    echo "cp $pres.dat multiwell.dat"  >> "MW$temp.$pres.$spec"
    echo "cp /home/edames/multiwell-2014/bin/multiwell ./" >> "MW$temp.$pres.$spec"
    echo "./multiwell" >> "MW$temp.$pres.$spec"
    echo "cp multiwell.out $pres.out"  >> "MW$temp.$pres.$spec"
    echo "cp multiwell.rate $pres.rate" >> "MW$temp.$pres.$spec"
    echo "rm multiwell.dat" >> "MW$temp.$pres.$spec"
    echo "exit 0" >> "MW$temp.$pres.$spec"   
    qsub "MW$temp.$pres.$spec"    
      
#go back up    
    cd ../../       
   done
    echo "make sure a,b,c.dat and sh1 are correct!" 
exit 0
