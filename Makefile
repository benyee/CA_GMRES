a.out: main.o SparseMat.o Utilities.o GMRES_sol.h
	g++ -g main.o SparseMat.o Utilities.o GMRES_sol.h -o a.out

<<<<<<< HEAD
main.o: main_comp_QR.cpp SparseMat.h Utilities.h
	g++ -g -c -I/usr/include main_comp_QR.cpp -o main.o
=======
main.o: main.cpp SparseMat.h Utilities.h
	g++ -g -c main.cpp -o main.o
>>>>>>> 6302d8f4f0b3bcc4b3e3abf2d2b17f59df03902f

SparseMat.o: SparseMat.cpp SparseMat.h
	g++ -g -c SparseMat.cpp -o SparseMat.o

Utilities.o: Utilities.cpp Utilities.h GMRES_sol.h SparseMat.h
	g++ -g -c Utilities.cpp -o Utilities.o

clean:
	rm -rf main.o SparseMat.o Utilities.o Utilities.h.gch SparseMat.h.gch a.out
