CC=gcc
CFLAGS= -g
CILK=g++
CILKFLAGS= -g -w -O3 -L/usr/local/lib -lcilkrts -fcilkplus -lprofiler -ldl -Wl,-rpath,/usr/local/lib
#LDFLAGS= -L$(CURDIR)
#AR=ar

all: main

main : main.cpp tfk_utils.h Makefile
	$(CILK) $(CILKFLAGS) $@.cpp $(LDFLAGS) -o $@

clean :
	rm -f main *~

run :
	./run.sh $(file) $(type) $(nworkers)
