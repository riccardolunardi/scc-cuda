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
run7B_compare:
	@ ./main_avanti.out ./samples/sample_graph7B > gen/avanti_mini.txt
	@ ./main_indietro.out ./samples/sample_graph7B > gen/indietro_mini.txt
run7C:
	@ ./main.out ./samples/sample_graph7C
run7D:
	@ ./main.out ./samples/sample_graph7D
run_gscc_fewu_compare:
	@ ./main_avanti.out ./samples/mid_tests/sample_test_gscc_fewu > gen/avanti_mini.txt
	@ ./main_indietro.out ./samples/mid_tests/sample_test_gscc_fewu > gen/indietro_mini.txt
run_gscc_moreu_compare:
	@ ./main_avanti.out ./samples/mid_tests/sample_test_gscc_moreu > gen/avanti_mini.txt
	@ ./main_indietro.out ./samples/mid_tests/sample_test_gscc_moreu > gen/indietro_mini.txt
run_pscc_fewu_compare:
	@ ./main_avanti.out ./samples/mid_tests/sample_test_scc_fewu > gen/avanti_mini.txt
	@ ./main_indietro.out ./samples/mid_tests/sample_test_scc_fewu > gen/indietro_mini.txt
run_pscc_moreu_compare:
	@ ./main_avanti.out ./samples/mid_tests/sample_test_scc_moreu > gen/avanti_mini.txt
	@ ./main_indietro.out ./samples/mid_tests/sample_test_scc_moreu > gen/indietro_mini.txt
run_killer_compare:
	@ ./main_avanti.out ./gen/killer > gen/avanti_mini.txt
	@ ./main_indietro.out ./gen/killer > gen/indietro_mini.txt
run_mini_compare:
	@ ./main_avanti.out ./gen/mini > gen/avanti_mini.txt
	@ ./main_indietro.out ./gen/mini > gen/indietro_mini.txt
compile_is_checked:
	@ g++ -o is_checked.out is_checked.cpp
run_is_checked:
	@ ./is_checked.out