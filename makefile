CFLAGS = -std=c99 -O2 -Wall -Werror -fopenmp

SRCS = src/run.c src/bin.c debug/debug.c src/fd_test.c src/wavelet.c

all:
	gcc $(CFLAGS) $(SRCS) -o run.out -lm

run: all
	./run.out
	
	rm run.out

clean:
	rm run.out
