CFLAGS = -std=c99 -O0 -g -march=native -Wall -fopenmp -mavx2 -mfma

SRCS = $(wildcard src/*.c)

TARGET = run.out

all:
	gcc $(CFLAGS) $(SRCS) -o $(TARGET) -lm

run: all
	./$(TARGET)
	
	rm $(TARGET)

clean:
	rm $(TARGET)

clean_snap:
	rm data/output/snapshots/*.bin

