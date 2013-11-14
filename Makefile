a.out: main.o SparseMat.o Utilities.o GMRES_sol.h
	g++ -g main.o SparseMat.o Utilities.o GMRES_sol.h -o a.out

main.o: main.cpp SparseMat.h Utilities.h
	g++ -g -c main.cpp -o main.o

SparseMat.o: SparseMat.cpp SparseMat.h
	g++ -g -c SparseMat.cpp -o SparseMat.o

Utilities.o: Utilities.cpp Utilities.h GMRES_sol.h SparseMat.h
	g++ -g -c Utilities.cpp -o Utilities.o

clean:
	rm -rf main.o SparseMat.o Utilities.o Utilities.h.gch SparseMat.h.gch a.out
