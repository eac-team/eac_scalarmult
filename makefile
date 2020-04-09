CC = gcc
FLAGS = -Wall
GMP = -lgmp

all: eac_measure

eac_measure: measure.c zaddu.o 
	$(CC) $(FLAGS) $^ -o $@ $(GMP) -O3

test_zaddu: test_zaddu.c zaddu.o
	$(CC) $(FLAGS) $^ -o $@ $(GMP)

zaddu.o: zaddu.c
	$(CC) -c $<



measure: eac_measure
	./eac_measure

test: test_zaddu
	./test_zaddu
