.PHONY:test

test:
	g++ -std=c++23 -march=native -D__disable_gcc_builtin__ -D__disable_check_overflow__ -O3 test.cpp -o test