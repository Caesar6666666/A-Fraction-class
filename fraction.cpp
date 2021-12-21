#include "fraction.h"

Fraction::Fraction() {
    op = 1;
    a = 0;
    b = 1;
}

Fraction::Fraction(long _a = 0, long _b = 1) {
    op = _a >= 0 ? (_b > 0 ? 1 : -1) : (_b > 0 ? -1 : 1);
    a = _a > 0 ? _a : -_a;
    b = _b > 0 ? _b : -_b;
    if(b == 0) {
        exit(-1);
    }
    Simplify();
}

Fraction::Fraction(ull _a, ull _b, int _op):a(_a),b(_b),op(_op) {
    if(_b == 0) {
        exit(-1);
    }
    Simplify();
}

Fraction::Fraction(const Fraction& u) {
    a = u.a;
    b = u.b;
    op = u.op;
}

std::istream& operator >> (std::istream& os,Fraction& u) {
    long _a,_b;
    char ch;
    os >> _a >> ch >> _b;
    u = Fraction(_a,_b);
    return os; 
}

void Fraction::Simplify() {
    ull g = gcd(a, b);
    a = a / g;
    b = b / g;
}

double Fraction::transform_to_float() {
    return (double)a / b * op;
}

ull Fraction::gcd(ull _x,ull _y) {
    return _y == 0 ? _x : gcd(_y, _x % _y);
}

std::ostream& operator << (std::ostream& os,const Fraction& u) {
    if(u.op > 0) os << u.a << '/' << u.b;
    else os << '-' << u.a << '/' << u.b;
    return os;
}

Fraction Fraction::operator + (const Fraction& u) {
    ull b1 = b * u.b;
    long a1 = a * u.b * op + u.a * b * u.op;
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 > 0 ? 1 : -1);
}

Fraction Fraction::operator - (const Fraction& u) {
    ull b1 = b * u.b;
    long a1 = a * u.b * op - u.a * b * u.op;
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 > 0 ? 1 : -1);
}

Fraction Fraction::operator * (const Fraction& u) {
    ull&& m1 = gcd(a,u.b);
    ull&& m2 = gcd(b,u.a);
    ull&& a1 = a / m1 * u.a / m2;
    ull&& b1 = b / m2 * u.b / m1;
    return Fraction(a1,b1,op * u.op);
}

Fraction Fraction::operator / (const Fraction& u) {
    if(u.a == 0) exit(-1);
    ull&& m1 = gcd(a,u.a);
    ull&& m2 = gcd(b,u.b);
    ull&& a1 = a / m1 * u.b / m2;
    ull&& b1 = b / m2 * u.a / m1;
    return Fraction(a1,b1,op * u.op);
}
