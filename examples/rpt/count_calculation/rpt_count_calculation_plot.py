"""
Name   : rpt_count_calculation_plot.py
Author : Amishga Alphonius
Date   : 16-06-2022
Desc   : Generates plots for the data resulting from the count calculation example
"""


from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

header = []
rows = []
position_x = []
position_y = []
position_z = []
counts = []

# Verify if an argument was added when running the script
assert(len(sys.argv)>1), "The .csv containing the data to be read wasn't specified as an argument when running the script."
	

# Opening and reading .csv file, and converting positions from 'm' to 'cm'
with open(sys.argv[1],'r') as data_file:
	csvReader = csv.reader(data_file, delimiter=',')
	header = next(csvReader)
	for row in csvReader:
		position_x.append(float(row[0])*100)
		position_y.append(float(row[1])*100)
		position_z.append(float(row[2])*100)
		counts.append(float(row[4]))
data_file.close()

# Generate plot

# Scenario 1
if (position_x[0]!=0 and position_y[0]!=position_x[0]):
	plt.plot(position_x, counts, linestyle="none", color="#CB5D09", marker='o', markersize=3, markerfacecolor="#CB5D09")
	plt.xlabel("x (cm)")
	plt.ylabel("Photon count")
	plt.grid(True)
	plt.savefig("./plots/results_horizontalx.png", dpi=300)
	plt.show()	
	
# Scenario 2
elif (position_x[0]==0 and  position_y[0]!=0):
	plt.plot(position_y, counts, linestyle="none", color="#CB5D09", marker='o', markersize=3, markerfacecolor="#CB5D09")
	plt.xlabel("y (cm)")
	plt.ylabel("Photon count")
	plt.grid(True)
	plt.savefig("./plots/results_horizontaly.png", dpi=300)
	plt.show()
	
# Scenario 3
elif (position_x[0]==0 and  position_y==position_x):
	plt.plot(position_z, counts, linestyle="none", color="#CB5D09", marker='o', markersize=1, markerfacecolor="#CB5D09")
	plt.xlabel("z (cm)")
	plt.ylabel("Photon count")
	plt.grid(True)
	plt.savefig("./plots/results_vertical.png", dpi=300)
	plt.show()
	
# Scenario 4
elif (position_x[0]!=0 and position_y==position_x):	
	fig = plt.figure()	
	ax = plt.axes(projection="3d")
	color_map = plt.get_cmap("turbo")
	p = ax.scatter3D(position_x, position_y, position_z, c=counts, cmap=color_map);
	ax.set_xlabel("x (cm)")
	ax.set_ylabel("y (cm)")
	ax.set_zlabel("z (cm)");
	ax.axes.set_xlim3d(left=-10, right=10) 
	ax.axes.set_ylim3d(bottom=-10, top=10) 
	ax.axes.set_zlim3d(bottom=0, top=30)
	cbar = fig.colorbar(p, location="left")
	cbar.set_label("Photon count", rotation=90)
	plt.savefig("./plots/results_diagonal.png", dpi=300)
	plt.show()

# For any other set of particle positions
else :
	fig = plt.figure()	
	ax = plt.axes(projection="3d")
	color_map = plt.get_cmap("turbo")
	p = ax.scatter3D(position_x, position_y, position_z, c=counts, cmap=color_map);
	ax.set_xlabel("x (cm)")
	ax.set_ylabel("y (cm)")
	ax.set_zlabel("z (cm)");
	cbar = fig.colorbar(p, location="left", pad=0.1)
	cbar.set_label("Photon count", rotation=90)
	plt.show()

	
