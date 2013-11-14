a.out: main.o SparseMat.o Utilities.o
	g++ -g -L/usr/lib64 main.o SparseMat.o Utilities.o -o a.out -lpapi

main.o: main_comp_QR.cpp SparseMat.h Utilities.h
	g++ -g -c -I/usr/include main_comp_QR.cpp -o main.o

SparseMat.o: SparseMat.cpp SparseMat.h
	g++ -g -c SparseMat.cpp -o SparseMat.o

Utilities.o: Utilities.cpp Utilities.h
	g++ -g -c Utilities.cpp -o Utilities.o

clean:
	rm -rf main.o SparseMat.o Utilities.o Utilities.h.gch SparseMat.h.gch a.out
