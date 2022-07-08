"""
Name   : rpt_parameter_tuning_plot.py
Author : Amishga Alphonius
Date   : 23-06-2022
Desc   : Generates plots with the data resulting from the parameter tuning example
"""

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys


header = []
rows = []
position_x = []
exp_counts = []
cal_counts = []

# Verify if an argument was added when running the script
assert(len(sys.argv)>1), "The .csv containing the data to be read wasn't specified as an argument when running the script."
	

# Opening and reading .csv file, and converting positions from 'm' to 'cm'
with open("counts_exp.csv",'r') as data_file:
	csvReader = csv.reader(data_file, delimiter=',')
	header = next(csvReader)
	for row in csvReader:
		exp_counts.append(float(row[0]))
data_file.close()

with open(sys.argv[1],'r') as data_file:
	csvReader = csv.reader(data_file, delimiter=',')
	header = next(csvReader)
	for row in csvReader:
		position_x.append(float(row[0])*100)
		cal_counts.append(float(row[4]))
data_file.close()

cal_counts = np.array(cal_counts)
exp_counts = np.array(exp_counts)


# Generate plots
# Experimental and calculated counts comparison
plt.plot(position_x, cal_counts, linestyle="dashed", color="#fc710d", linewidth=2, label="tuned parameters")
plt.plot(position_x, exp_counts, linestyle="none", color="#000000", marker='o', ms=5, mfc="None", label="experimental")
plt.xlabel("Photon count")
plt.ylabel("x (cm)")
plt.legend()
plt.grid(True)
plt.show()

# Linear fit graph
plt.plot(exp_counts, cal_counts, linestyle="None", color="#000000", marker='o', ms=5, mfc="None")
cal_counts, exp_counts = cal_counts.reshape(-1,1), exp_counts.reshape(-1,1)
plt.plot(cal_counts, LinearRegression().fit(cal_counts, exp_counts).predict(cal_counts), linestyle="dashed", color="#fc710d", linewidth=2)
plt.annotate("R^2 = {:.4f}".format(r2_score(cal_counts, exp_counts)), (490, 650), color="#fc710d")
plt.xlabel("Experimental photon count")
plt.ylabel("Calculated photon count")
plt.grid(True)
plt.show()
