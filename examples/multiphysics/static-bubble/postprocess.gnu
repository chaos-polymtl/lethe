# Set the output to a png file
set terminal pngcairo  enhanced font "arial,14" fontscale 1.0 size 600, 400


# The file we'll write to
set output ARG1."/L2Error.png"

# Set the label of the axis
set xlabel "Time [s]"
set ylabel "||u||_{L2} [-]"

# Set logarithmic
set logscale y

# Set the position of the legend
set key bottom left

# Enable the grid
set grid

# Plot the CD and CL evolutions
currentPath=ARG1."/L2Error.dat"
plot currentPath every ::1 using 1:2 with lines notitle

set terminal x11 enhanced font "arial,10" size 600, 400
set output
replot

pause -1
