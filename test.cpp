#include "fraction.h"
#include <format>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

#define frac(x) (x).get_sign() == 1 ? ' ' : '-', (x).get_numer(), (x).get_denom()

int main() {
    default_random_engine e;
    e.seed(chrono::system_clock::now().time_since_epoch().count());
    auto start = chrono::high_resolution_clock::now();
    for (int i = 1; i <= 10000000; ++i) {
        auto num1 = e();
        auto num2 = e();
        auto den1 = e();
        auto den2 = e();
        auto f1 = caesar::Fraction<unsigned long long>(num1, den1, -1);
        auto f2 = caesar::Fraction<unsigned long long>(num2, den2, 1);
        auto &&f3 = f1 + f2;
        if (fabs(f3.ConvToFloat() - (((double)num1 / den1 * -1) + ((double)num2 / den2)))
            > fabs(0.00001 * f3.ConvToFloat()))
            cout << format("Error: {}{}/{} + {}{}/{} = {}{}/{}", frac(f1), frac(f2), frac(f3)) << endl;
    }
    auto end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
    return 0;
}
// -732564845/821930932 - 173400416/394589155 = -1/4 -0.25 -0.451828