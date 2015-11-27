# Gnuplot
#
# plot 1D data for several stages
#
#

# set xrange [-0.5:0.5]
# set yrange [-1:1]

plot "DATA/st000000.txt" using 1:2 w lp pt 7 ps 1.5, \
     "DATA/st000092.txt" using 1:2 w lp pt 7 ps 1.5

pause -1 "press [Enter] key to quit " # prompt