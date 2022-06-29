# Set the output to a png file
set terminal pngcairo

# The file we'll write to
set output ARG1."/CL-CD.png"

# Set the label of the axis
set xlabel "Time [s]"
set ylabel "CL and CD [-]"

# Set the position of the legend
set key bottom left

# Enable the grid
set grid

# Plot the CD and CL evolutions
currentPath=ARG1."/force.00.dat"
plot [0:200] [-1.0:1.5] currentPath using 1:(2*$2) with lines title "CD", currentPath using 1:(2*$3) with lines title "CL"

set terminal x11
set output
replot

pause -1
