CFLAGS = -std=c99 -Wall -Werror

SRCS = main.c src/read_bin.c debug/debug.c src/FD.c

all:
	gcc $(CFLAGS) $(SRCS) -o main.out -lm

run: all
	./main.out
	
	rm main.out

clean:
	rm main.out