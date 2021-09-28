#!/bin/bash

# Rules of Usage
# Accept at most 5 input arguments
# Argument $1 must be rows [mandatory]
# Argument $2 must be cols [mandatory]
# If data file name is to be input, it must be placed at $3
# Argument $3-5 (optional) can either be "3d", "r2" or "r3" 
# Hardcoded strings: "3d" - plot a 3D map; "r2" - replot every 2 secs; "r3" - replot every 3 sec.
# Press Ctrl-C to stop gnuplot if you use "r2"/"r3".

# Examples of Usage
# ./plot.sh 200 200			Plot jacobi.dat as a colormap of size 200 x 200
# ./plot.sh 200 200 r2			Plot jacobi.dat as a colormap of size 200 x 200, replot every 2 seconds 
# ./plot.sh 200 200 r3			Plot jacobi.dat as a colormap of size 200 x 200, replot every 3 seconds
# ./plot.sh 200 200 jacobi_seq.dat 3d	Plot jacobi_seq.dat as a 3D colormap of size 200 x 200
# ./plot.sh 200 200 3d r3		Plot jacobi.dat as a 3D colormap of size 200 x 200, replot every 3 seconds

rows=$1
cols=$2

file="jacobi.dat"
if [ "$3" != "" -a "$3" != "3d" -a "$3" != "r2" -a "$3" != "r3" ] ; then
   file=$3
fi

refresh="pause -1 \"Hit return to exit\""
if   [ "$3" = "r2" -o "$4" = "r2" -o "$5" = "r2" ] ; then
   refresh="pause 2; replot; reread;"
elif [ "$3" = "r3" -o "$4" = "r3" -o "$5" = "r3" ] ; then
   refresh="pause 3; replot; reread;"
fi

style="set pm3d map"
if [ "$3" = "3d" -o "$4" = "3d" -o "$5" = "3d" ] ; then
   style="set pm3d at b"
   view3d="yes"
fi

# Generate the gnuplot script
echo "$style"                               > jacobi.p
if [ "$rows" = "$cols" ] ; then
   echo "set size square"                  >> jacobi.p
fi
echo "set xrange [0:$rows]"                >> jacobi.p
echo "set yrange [$cols:0]"                >> jacobi.p
if [ "$view3d" = "yes" ] ; then
   echo "set view 70, 35"                  >> jacobi.p
fi
echo "set palette rgbformulae 22,13,-31"   >> jacobi.p
echo "splot \"$file\""		           >> jacobi.p
echo "$refresh"                            >> jacobi.p

gnuplot jacobi.p
