#include <iostream>
#include "fraction.h"

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
    if(u.mark() > 0) os << u.num() << '/' << u.deno();
    else os << '-' << u.num() << '/' << u.deno();
    return os;
}

int main() {
    int op;
    while(std::cin >> op) {
        caesar::Fraction<long long> f1;
        caesar::Fraction<int> f2;
        std::cin >> f1 >> f2;
        std::cout << f1 + f2 << std::endl;
        std::cout << f1 - f2 << std::endl;
        std::cout << f1 * f2 << std::endl;
        std::cout << f1 / f2 << std::endl;
        std::cout << (f1 < f2) << std::endl;
        std::cout << (f1 > f2) << std::endl;
        std::cout << (f1 == f2) << std::endl;
    }
}