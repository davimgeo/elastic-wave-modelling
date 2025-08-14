CFLAGS = -std=c99 -Wall -Werror

SRCS = run.c src/bin.c debug/debug.c src/fd.c src/wavelet.c

all:
	gcc $(CFLAGS) $(SRCS) -o run.out -lm

run: all
	./run.out
	
	rm run.out

clean:
	rm run.out