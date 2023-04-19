#include <cstdlib>
#include <functional>
#include <iostream>
#include "fraction.h"
#include <sys/time.h>
#include <random>
#include <gperftools/profiler.h>
using namespace std;
template<typename T>
std::istream& operator >> (std::istream& is,caesar::Fraction<T>& u) {
    T _a,_b;
    char ch;
    is >> _a >> ch >> _b;
    u = caesar::Fraction<T>(_a,_b);
    return is;
}

template<typename T>
std::ostream& operator << (std::ostream& os,const caesar::Fraction<T>& u) {
    if(u.mark() > 0) os << u.get_numer() << '/' << u.get_denom();
    else os << '-' << u.get_numer() << '/' << u.get_denom();
    return os;
}

double operator - (const timeval& t1, const timeval& t2) {
    return t1.tv_sec - t2.tv_sec + (double)(t1.tv_usec - t2.tv_usec) / 1000000;
}

int main() {
    default_random_engine e;
    caesar::Fraction<unsigned long long>f1(5,2,1);
    caesar::Fraction<unsigned long long>f2(7,3,1);
    for(int i = 1;i <= 100000000;++ i) {
        auto&& f3 = f1 + f2;
        f1 = f2;
        f2 = f3;
        //if(fabs(f3.ConvToFloat() - (double)(-a * b1 + a1 * b) / (b * b1)) > 1e-6) cout << f1 << endl << f2 << endl<< "f3 计算值:" << f3.ConvToFloat() << endl << "理论值:" << (double)(a * b1 + a1 * b) / (b * b1) << endl << "error" << endl;
    }
    cout << f2;
    return 0;
}