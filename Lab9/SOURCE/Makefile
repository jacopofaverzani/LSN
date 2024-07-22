CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: main.exe
main.exe : genetic.o cities.o main.o random.o
	$(CC) genetic.o cities.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp genetic.hpp cities.hpp path.hpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
genetic.o : genetic.cpp genetic.hpp random.hpp cities.hpp path.hpp
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
cities.o : cities.cpp cities.hpp
	$(CC) -c cities.cpp -o cities.o $(CFLAGS)
random.o : random.cpp random.hpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 
remove:
	rm ../OUTPUT/*.*
