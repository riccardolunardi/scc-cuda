import os
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np

def print_stdout(output):
    print(output.stdout.decode('utf-8').replace('\r', ''))

print("Compiling CUDA code...")
compiler_output = subprocess.run(["nvcc", "-Xcompiler", "/openmp", "-DDEBUG_FINAL=0", "-DOMP_MIN_NODES=100000", "-DWARMUP=5", ".\cuda\scc_runner.cu", "-o", "./build/scc.exe"], capture_output=True)
print_stdout(compiler_output)

FOLDER_PATH = 'C:/Users/rikyn/Desktop/final/'

files = os.listdir(FOLDER_PATH)

for file in files:
    if file.endswith('.txt'):
        print(f"Benchmarking {file}...")
        
        try:
            output = subprocess.run(['./build/scc.exe', FOLDER_PATH + file, '25', '1'], capture_output=True)
            print_stdout(output)
            output_string = output.stdout.decode('utf-8')

            graph_names = re.findall(r"\/(\w+)\.txt", output_string)
            
            avg = []
            std = []
            for time_data in re.findall(r"\d+\.\d+,\d+\.\d+", output_string):
                avg_run, std_run = time_data.split(",")
                avg.append(float(avg_run))
                std.append(float(std_run))

            versions_runned = re.findall(r"Versione \d -(.+)-", output_string)
            
            print(versions_runned)
            print(avg)
            print(std)

            with open(f'./benchmark/result/{graph_names[0]}.csv', 'w') as f:
                f.write(f"Versione\tAvg\tStd\n")
                for i in range(len(versions_runned)):
                    f.write(f"{versions_runned[i]}\t{avg[i]}\t{std[i]}\n")

            for i in range(0,2):
                plt.figure()
                fig, ax = plt.subplots()

                plt.bar(versions_runned[i:], avg[i:], yerr=std[i:], error_kw=dict(lw=5, capsize=5, capthick=3))

                # Add labels and title
                ax.set_xlabel('Versione', labelpad=-3)
                ax.set_ylabel('Time (ms)')
                ax.set_title(graph_names[0].capitalize(), fontdict={'fontsize': 20, 'fontweight': 'heavy'}, pad=15)
                ax.set_xticks(np.arange(len(versions_runned[i:])))
                ax.set_xticklabels([v.replace(" ", "\n") for v in versions_runned[i:]], fontsize=8)

                # Show the plot
                plt.savefig(f'./benchmark/result/{graph_names[0]}_{i}.png')
        except Exception as e:
            print("Error: ", e, " - skipping file...")