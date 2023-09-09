CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o classes.o
	$(CC) random.o classes.o main.o -o main.exe
main.o : main.cpp globals.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
classes.o : classes.cpp classes.h
	$(CC) -c classes.cpp -o classes.o $(CFLAGS)
clean :
	rm *.o main.exe
