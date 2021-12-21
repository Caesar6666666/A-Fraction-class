.PHONY:test

task:fraction.o
	ar -src libfraction.a fraction.o
	g++ -shared -fpic -o libfraction.so fraction.o

fraction.o:fraction.cpp
	g++ -c fraction.cpp

test:
	g++ test.cpp -o test -I. -L. -lfraction