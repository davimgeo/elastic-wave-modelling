CFLAGS = -std=c99 -O3 -g -march=native -Wall -fopenmp -mavx2 -mfma


SRCS = src/main.c src/bin.c src/wavelet.c src/fd.c 

all:
	gcc $(CFLAGS) $(SRCS) -o run.out -lm

run: all
	./run.out
	
	rm run.out

clean:
	rm run.out

clean_snap:
	rm data/output/snapshots/*.bin
