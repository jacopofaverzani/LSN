CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o
	$(CC) random.o funzioni.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.hpp random.h
	$(CC) -c funzioni.cpp funzioni.hpp
clean :
	rm *.o *.gch *:Zone.Identifier
