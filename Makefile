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