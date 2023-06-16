#include "fraction.h"
#include <cstdlib>
#include <ctime>
#include <format>
#include <functional>
#include <iostream>
#include <random>
using namespace std;

#define frac(x) (x).get_sign() == 1 ? ' ' : '-', (x).get_numer(), (x).get_denom()

int main() {
    // mersenne_twister_engine<unsigned long long, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29,
    //                         0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL,
    //                         43, 6364136223846793005ULL>
    //     e1;
    default_random_engine e;
    e.seed(time(NULL));
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
    // cout << caesar::Fraction<uint>(1,2) << endl;
    // cout << caesar::Fraction<uint64_t>(732564845, 821930932, -1)
    //             + caesar::Fraction<uint64_t>(173400416, 394589155,1)
    //      << endl;
    // cout << -(double)732564845 / 821930932 + (double)173400416 / 394589155 << endl;
    // cout << caesar::gcd(146538977639188263ull,324325031926242460ull) << endl;
    return 0;
}
// -732564845/821930932 - 173400416/394589155 = -1/4 -0.25 -0.451828