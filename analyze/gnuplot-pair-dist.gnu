df = "pair-dist.dat"

set term postscript eps enhanced color font "Times-Roman,25"


set pm3d map
set palette rgb 21,22,23
#set size square 1,1
set size square 0.8,1.0

set lmargin at screen 0.1
set rmargin at screen 0.65
set bmargin at screen 0.17
set tmargin at screen 0.98

set xlabel "x"
set ylabel "y"

set output "pair-dist.eps"
splot "pair-dist.dat" notitle

set xlabel "k_x / 2{/Symbol p}"
set ylabel "k_y / 2{/Symbol p}"

set xrange [-6:6]
set yrange [-6:6]

set output "pair-dist-fft.eps"
splot "pair-dist-fft.dat" notitle
