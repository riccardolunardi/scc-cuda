nvcc -Xcompiler /openmp -DDEBUG_FINAL=0 -DOMP_MIN_NODES=100000 .\cuda\scc_runner.cu -o "./build/scc.exe"
echo > ./benchmarks_partial_results.txt
for /F %%f in ('dir /B /O:N "%cd%\samples\iter_tests\sample_test_iter*"') do %cd%\build\scc.exe %cd%\samples\iter_tests\%%f 25 1 >> ./benchmarks_partial_results.txt
