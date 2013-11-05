a.out: main.o SparseMat.o SparseVec.o
	g++ -g main.o SparseMat.o SparseVec.o -o a.out

main.o: main.cpp SparseVec.h SparseMat.h
	g++ -g -c main.cpp -o main.o

SparseMat.o: SparseMat.cpp SparseMat.h SparseVec.h
	g++ -g -c SparseMat.cpp -o SparseMat.o

SparseVec.o: SparseVec.cpp SparseVec.h
	g++ -g -c SparseVec.cpp -o SparseVec.o

clean:
	rm -rf main.o SparseMat.o SparseVec.o SparseVec.h.gch SparseMat.h.gch a.out
