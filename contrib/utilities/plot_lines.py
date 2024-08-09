import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sys

import pandas as pd
from datetime import datetime

colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
#Plot font and colors
font = {'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

# read text file into pandas DataFrame and create 
# header with names.
df = pd.read_csv(sys.argv[1], sep=" ", 
                 names=["Date", "Code","Test"])

# transform date to datetime
date_list=[]
for j in df["Date"]: 
  date_list.append((datetime.strptime(j,'%Y-%m-%d')))

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot_date(date_list,df["Code"]/1000,ms=4,color=colors[0],label="Code")
plt.plot_date(date_list,df["Test"]/1000,ms=4,color=colors[1],label="Tests")
plt.xlabel("Date")
plt.ylabel("Number of lines $\\cdot 10^3$")

plt.legend()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))
plt.tight_layout()
plt.savefig("lines.png",dpi=300)
plt.show()