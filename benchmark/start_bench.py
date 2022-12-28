import os
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np

def print_stdout(output):
    print(output.stdout.decode('utf-8').replace('\r', ''))

print("Compiling CUDA code...")
compiler_output = subprocess.run(["nvcc", "-Xcompiler", "/openmp", "-DDEBUG_FINAL=1", "-DOMP_MIN_NODES=100000", "-DWARMUP=5", ".\cuda\scc_runner.cu", "-o", "./build/scc.exe"], capture_output=True)
print_stdout(compiler_output)

FOLDER_PATH = 'F:/network-benchmark/final/'

files = os.listdir(FOLDER_PATH)
for file in files:
    if file.endswith('.txt') and file.startswith('twitter'): # or file.startswith('amazon')
        print(f"Benchmarking {file}...")
        
        output = subprocess.run(['./build/scc.exe', FOLDER_PATH + file, '5', '0'], capture_output=True)
        print_stdout(output)
        output_string = output.stdout.decode('utf-8')

        graph_names = re.findall(r"\/(\w+)\.txt", output_string)
        
        avg = []
        std = []
        for time_data in re.findall(r"\d+\.\d+,\d+\.\d+", output_string):
            avg_run, std_run = time_data.split(",")
            avg.append(float(avg_run))
            std.append(float(std_run))

        versions_runned = re.findall(r"Versione \d", output_string)
        
        print(versions_runned)
        print(avg)
        print(std)

        x = np.arange(len(versions_runned))  # the label locations
        width = 0.35  # the width of the bars
    
        plt.bar(versions_runned, avg, label='Values 1')
        plt.bar(versions_runned, std, bottom=avg, label='Values 2')
        fig, ax = plt.subplots()
        rects1 = ax.bar(x - width/2, avg, width, label='Values 1', yerr=std, error_kw=dict(lw=5, capsize=5, capthick=3))

        # Add labels and title
        ax.set_xlabel('Versions')
        ax.set_ylabel('Values')
        ax.set_title('Bar Plot Example')
        ax.set_xticks(x)
        ax.set_xticklabels(versions_runned)

        # Add the legend
        ax.legend()

        # Show the plot
        plt.show()

        data_run = zip(versions_runned, avg, std)

        for a,b,c in data_run:
            print(f"{a}, {b}, {c}")
        
        break