help:
	@make -qpRr | egrep -v "^Makefile:" | egrep -e '^[^\.\#].*:$$' | sed -e 's~:~~g' | sort

compile:
	@ g++ -o main.out main.cpp

compile-win:
	@ g++ -o ./build/main_cpp.exe main.cpp

run:
ifdef n
	@ ./main.out "./samples/sample_graph$(n)"
else
	@ echo "Specifica il numero di grafo su cui vuoi testare lo script"
endif

run-win:
ifdef n
	@ ./main.exe "./samples/sample_graph$(n)"
else
	@ echo "Specifica il numero di grafo su cui vuoi testare lo script"
endif


cuda-compile:
ifdef version
	@ nvcc ./cuda/$(version).cu -o ./build/$(version).exe
else
	@ nvcc ./cuda/naive.cu -o ./build/naive.exe
endif

cuda-run:
ifdef version 
ifdef test
	@ ./build/$(version).exe ./samples/$(test)
else
	@ ./build/$(version).exe ./samples/killer
endif
else
	@ nvcc ./cuda/naive.cu -o ./samples/killer
endif

cuda-compile-omp:
	@ nvcc -Xcompiler /openmp ./cuda/sccv5_openmp.cu -o ./build/sccv5_openmp.exe

cuda-compile-all:
	@ nvcc ./cuda/sccv1_naive.cu -o ./build/sccv1_naive.exe
	@ nvcc ./cuda/sccv2_status.cu -o ./build/sccv2_status.exe
	@ nvcc ./cuda/sccv3_streams.cu -o ./build/sccv3_streams.exe
	@ nvcc ./cuda/sccv4_pinned.cu -o ./build/sccv4_pinned.exe
	@ nvcc -Xcompiler /openmp ./cuda/sccv5_openmp.cu -o ./build/sccv5_openmp.exe

cuda-test-all:
	@ .\build\sccv1_naive.exe .\samples\mid_tests\sample_test_scc_fewu
	@ .\build\sccv1_naive.exe .\samples\final_tests\sample_test_scc_moreu

	@ .\build\sccv2_status.exe .\samples\mid_tests\sample_test_scc_fewu
	@ .\build\sccv2_status.exe .\samples\final_tests\sample_test_scc_moreu

	@ .\build\sccv3_streams.exe .\samples\mid_tests\sample_test_scc_fewu
	@ .\build\sccv3_streams.exe .\samples\final_tests\sample_test_scc_moreu

	@ .\build\sccv4_pinned.exe .\samples\mid_tests\sample_test_scc_fewu
	@ .\build\sccv4_pinned.exe .\samples\final_tests\sample_test_scc_moreu
	
	@ .\build\sccv5_openmp.exe .\samples\mid_tests\sample_test_scc_fewu
	@ .\build\sccv5_openmp.exe .\samples\final_tests\sample_test_scc_moreu

profile-cpp:
	@ g++ -std=c++11 -pg -no-pie .\main.cpp -o main_prof.exe
	@ .\main_prof.exe .\samples\final_tests\sample_test_scc_moreu
	@ gprof .\main_prof.exe > profiled_cpp.txt

generate:
	@ python3 ./gen/generate.py

cuda-compile-test:
	@ nvcc -Xcompiler /openmp .\cuda\scc_runner.cu -o ./build/scc.exe
	@ ./build/scc.exe .\samples\final_tests\sample_test_gscc_fewu 2 0
	@ ./build/scc.exe .\samples\mid_tests\sample_test_scc_fewu 2 0
	@ ./build/scc.exe .\samples\final_tests\sample_test_scc_fewu 2 0

cuda-test-only:
	@ ./build/scc.exe .\samples\mid_tests\sample_test_scc_fewu 10 0