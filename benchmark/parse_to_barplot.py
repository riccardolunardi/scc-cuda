import os
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import traceback
import csv

FOLDER_PATH = 'F:/network-benchmark/final/'
graph_name = "soclivejournal1_not_u"

with open(f'./benchmark/result/{graph_name}.csv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    # Initialize empty lists for each column
    versions_runned = []
    avg = []
    std = []
    
    # Iterate over the rows of the CSV file
    for row in reader:
        # Append the values in each column to the corresponding list
        if row[0] == "Versione":
            continue

        versions_runned.append(row[0])
        avg.append(float(row[1]))
        std.append(float(row[2]))

print(versions_runned, avg, std, sep="\n")

plt.figure()
fig, ax = plt.subplots()

plt.bar(versions_runned, avg, yerr=std, error_kw=dict(lw=5, capsize=5, capthick=3))

# Add labels and title
ax.set_xlabel('Versione', labelpad=-3)
ax.set_ylabel('Time (ms)')

name = graph_name
name = "SOC-LiveJournal1"
if "not_u" in graph_name:
    name = graph_name.capitalize().replace("_not_u", " (Nodi U opposti)")

ax.set_title(name, fontdict={'fontsize': 20, 'fontweight': 'heavy'}, pad=15)
ax.set_xticks(np.arange(len(versions_runned)))
ax.set_xticklabels([v.replace(" ", "\n") for v in versions_runned], fontsize=8)

# Show the plot
plt.savefig(f'./benchmark/result/{graph_name}.png')