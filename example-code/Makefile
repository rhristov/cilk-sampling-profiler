CC=gcc
CFLAGS= -g
CILK=g++
CILKFLAGS= -g -w -O3 -fcilkplus -lcilkrts -lprofiler -Wl,-rpath,/usr/local/lib
LDFLAGS= -L/usr/local/lib
#AR=ar

all: main

main : main.cpp tfk_utils.h Makefile
	$(CILK) $(CILKFLAGS) $@.cpp $(LDFLAGS) -o $@

clean :
	rm -f main *~

run :
	./run.sh $(file) $(type) $(nworkers)
