CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : VariationalMC.o random.o
	$(CC) random.o VariationalMC.o -o VariationalMC.exe
MD_MC.o : VariationalMC.cpp
	$(CC) -c VariationalMC.cpp -o VariationalMC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VariationalMC.exe
