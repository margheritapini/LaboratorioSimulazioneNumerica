CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: main.exe
main.exe : Genetic_Algorithm.o main.o city.o random.o
	$(CC) Genetic_Algorithm.o city.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp Genetic_Algorithm.h city.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
Genetic_Algorithm.o : Genetic_Algorithm.cpp Genetic_Algorithm.h random.h city.h
	$(CC) -c Genetic_Algorithm.cpp -o Genetic_Algorithm.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
