help:
	@make -qpRr | egrep -v "^Makefile:" | egrep -e '^[^\.\#].*:$$' | sed -e 's~:~~g' | sort

generate:
	@ python3 ./gen/generate.py

cpp-compile:
	@ g++ -std=c++11 standalone.cpp -o ./build/standalone.exe

cpp-run:
	@ ./build/standalone.exe F:/network-benchmark/final/twitter.txt

cuda-compile:
	@ nvcc -Xcompiler /openmp -DDEBUG_FINAL=1 -DOMP_MIN_NODES=100000 -DWARMUP=0 ./cuda/scc_runner.cu -o ./build/scc.exe

cuda-run:
	@ ./build/scc.exe F:/network-benchmark/final/twitter.txt 1