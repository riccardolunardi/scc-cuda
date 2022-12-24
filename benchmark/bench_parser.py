SEPARATOR = ";"
SIMBOLO_INIZIO_DECIMALI = ","

headers = ["Nodes","Avg(v0)","Std(v0)","Avg(v1)","Std(v1)","Avg(v2)","Std(v2)","Avg(v3)","Std(v3)","Avg(v4)","Std(v4)","Avg(v5)","Std(v5)","Avg(v6)","Std(v6)","Avg(v7)","Std(v7)"]
headers = SEPARATOR.join(headers)
csv_res = []
csv_line = []

with open("benchmark/benchmarks_partial_results.txt") as f:
    for line in f.readlines():
        if line.startswith("Number"):
            if csv_line:
                csv_res.append(csv_line)

            csv_line = []
            csv_line.append(line.split(" ")[-1].rstrip())
        elif line.startswith("Versione"):
            bench = line.split("-")[-1].split(",")
            bench = [format(float(b.rstrip()), ".5f").replace(".", f"{SIMBOLO_INIZIO_DECIMALI}") for b in bench]
            csv_line.extend(bench)

with open("benchmarks_results.csv", "w") as f:
    f.write(headers+"\n" + "\n".join([";".join(line) for line in csv_res]))

print("Done!")