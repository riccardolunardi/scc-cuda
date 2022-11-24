help:
	@make -qpRr | egrep -v "^Makefile:" | egrep -e '^[^\.\#].*:$$' | sed -e 's~:~~g' | sort
compile:
	@ g++ -o main.out main.cpp
run:
	@ ./main.out ./samples/sample_graph
run1:
	@ ./main.out ./samples/sample_graph1
run2:
	@ ./main.out ./samples/sample_graph2
run3:
	@ ./main.out ./samples/sample_graph3
run4:
	@ ./main.out ./samples/sample_graph4
run4T:
	@ ./main.out ./samples/sample_graph4T
run5:
	@ ./main.out ./samples/sample_graph5
run5B:
	@ ./main.out ./samples/sample_graph5B
run6:
	@ ./main.out ./samples/sample_graph6
run6B:
	@ ./main.out ./samples/sample_graph6B
run7A:
	@ ./main.out ./samples/sample_graph7A
run7B:
	@ ./main.out ./samples/sample_graph7B
run7C:
	@ ./main.out ./samples/sample_graph7C
run7D:
	@ ./main.out ./samples/sample_graph7D
run_fewu:
	@ ./main.out ./samples/mid_tests/sample_test_gscc_fewu
