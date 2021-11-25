CC=gcc
PFLAGS=-lpthread

qsort: qsortseq.c
	$(CC) -O2 -o qsort.out qsortseq.c $(PFLAGS)
gaus: gaussian.c
	$(CC) -O2 -o gaus.out gaussian.c $(PFLAGS)
all: qsort gaus


clean:
	rm -rf *.out
run: qsort gaus
	./qsort.out
	./gaus.out
