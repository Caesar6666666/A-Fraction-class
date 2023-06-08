.PHONY:test

sharelib:fraction_share.o
	g++ -shared -fpic -o libfraction.so fraction_share.o

staticlib:fraction_static.o
	ar -src libfraction.a fraction_static.o

fraction_static.o:fraction.cpp
	g++ -c fraction.cpp -o fraction_static.o

fraction_share.o:fraction.cpp
	g++ -fpic -c fraction.cpp -o fraction_share.o

test:
	g++ -std=c++20 -O2 test.cpp -o test -lprofiler