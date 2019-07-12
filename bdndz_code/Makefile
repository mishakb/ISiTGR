CC = gcc
default: dndz

clean:
	rm -f *.o *.x

dndz: dndz_func.c dndz_func_scale.c dndz_bin.c smith2.c linalg.c read.c min1d.c \
    utils_func.c utils_func_scale.c utils_bin.c
	$(CC) smith2.c -c -o smith2.o -w -O3
	$(CC) dndz_func.c smith2.o -o dndz_func.x -lm -w -O3 -DSAVE_TABLES
	$(CC) dndz_func_scale.c smith2.o -o dndz_func_scale.x -lm -w -O3
	$(CC) dndz_bin.c smith2.o -o dndz_bin.x -lm -w -O3

