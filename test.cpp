#include "fraction.h"
#include <cstdlib>
#include <functional>
#include <gperftools/profiler.h>
#include <iostream>
#include <random>
#include <sys/time.h>
using namespace std;
template <typename T>
std::istream &operator>>(std::istream &is, caesar::Fraction<T> &u) {
    T _a, _b;
    char ch;
    is >> _a >> ch >> _b;
    u = caesar::Fraction<T>(_a, _b);
    return is;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const caesar::Fraction<T> &u) {
    if (u.get_sign() > 0)
        os << u.get_numer() << '/' << u.get_denom();
    else
        os << '-' << u.get_numer() << '/' << u.get_denom();
    return os;
}

double operator-(const timeval &t1, const timeval &t2) {
    return t1.tv_sec - t2.tv_sec + (double)(t1.tv_usec - t2.tv_usec) / 1000000;
}

int main() {
    // mersenne_twister_engine<unsigned long long, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29,
    //                         0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43,
    //                         6364136223846793005ULL>
    //     e1;
    default_random_engine e;
    e.seed(time(NULL));
    for (int i = 1; i <= 10000000; ++i) {
        auto num1 = e();
        auto num2 = e();
        auto den1 = e();
        auto den2 = e();
        int8_t sign1 = (i / 2) % 2 == 0 ? 1 : -1;
        int8_t sign2 = i % 2 == 0 ? 1 : -1;
        auto f1 = caesar::Fraction<unsigned long long>(num1, den1, sign1);
        auto f2 = caesar::Fraction<unsigned long long>(num2, den2, sign2);
        auto &&f3 = f1 / f2;
        if (fabs(f3.ConvToFloat() - (((double)num1 / den1 * sign1) / ((double)num2 / den2 * sign2)))
            > fabs(0.00001 * f3.ConvToFloat()))
            cout << f1 << " / " << f2 << " = " << f3 << " " << f3.ConvToFloat() << " "
                 << ((((double)num1 / den1 * sign1) / ((double)num2 / den2 * sign2))) << endl;
    }
    // cout << caesar::Fraction<uint>(1,2) << endl;
    // cout << caesar::Fraction<uint64_t>(732564845, 821930932, -1)
    //             + caesar::Fraction<uint64_t>(173400416, 394589155,1)
    //      << endl;
    // cout << -(double)732564845 / 821930932 + (double)173400416 / 394589155 << endl;
    //cout << caesar::gcd(146538977639188263ull,324325031926242460ull) << endl;
    return 0;
}
// -732564845/821930932 - 173400416/394589155 = -1/4 -0.25 -0.451828