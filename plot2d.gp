# Gnuplot
#
# plot 2D data
#
#

##### 3D lines
# set xrange [-0.5:0.5]
# set yrange [-0.5:0.5]
# set dgrid3d 64,64
# set hidden3d
# splot "DATA/st000000.txt" using 1:2:3 with lines
# splot "DATA/st000183.txt" using 1:2:3 with lines

##### 2D image
set pm3d map
set size ratio -1
splot "DATA/st000183.txt" using 1:2:3 with image

pause -1 "press [Enter] key to quit " # prompt