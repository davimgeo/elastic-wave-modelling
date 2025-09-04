CFLAGS = -std=c99 -O3 -march=native -Wall -fopenmp -mavx2 -mfma

SRCS = src/main.c src/bin.c debug/debug.c src/fd.c 

all:
	gcc $(CFLAGS) $(SRCS) -o run.out -lm

run: all
	./run.out
	
	rm run.out

clean:
	rm run.out
