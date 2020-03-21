CC = gcc
FLAGS = -Wall
GMP = -lgmp


eac_measure: measure.c zaddu.o
	$(CC) $(FLAGS) $^ -o $@ $(GMP) -O3

zaddu.o: zaddu.c
	$(CC) -c $<

measure: eac_measure
	./eac_measure

