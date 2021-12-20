#include "fraction.h"

Fraction::Fraction(long _a = 0, long _b = 1) {
    op = _a > 0 ? (_b > 0 ? 1 : -1) : (_b > 0 ? -1 : 1);
    if(_b == 0) op = 0;
    a = _a > 0 ? _a : -_a;
    b = _b > 0 ? _b : -_b;
    if(b == 0) {
        std::cout << "ERROR" << std::endl;
        exit(0);
    }
    Simplify();
}

Fraction::Fraction(ull _a,ull _b,int _op):a(_a),b(_b),op(_op){
    Simplify();
}

std::istream& operator >> (std::istream& os,Fraction& u) {
    ull _a,_b;
    char ch;
    os >> _a >> ch >> _b;
    u.a = _a;
    u.b = _b;
    u.Simplify();
    return os; 
}

void Fraction::Simplify() {
    ull g = gcd(a, b);
    a = a / g;
    b = b / g;
}

double Fraction::change_to_float() {
    if(op == 0) {
        std::cout << "ERROR!";
        exit(0);
    }
    return (double)a / b * op;
}

ull Fraction::gcd(ull _x,ull _y) {
    return _y == 0? _x : gcd(_y, _x % _y);
}

std::ostream& operator << (std::ostream& os,const Fraction& u) {
    if(u.op > 0) os << u.a << '/' << u.b;
    else if(u.op < 0) os << '-' << u.a << '/' << u.b;
    else if(u.op == 0) {
        os << "ERROR!";
        exit(0);
    }
    return os;
}

Fraction Fraction::operator + (const Fraction& u) {
    if(op == 0 || u.op == 0) {
        std::cout << "ERROR" << std::endl;
        exit(0);
    }
    ull b1 = b * u.b;
    long a1 = a * u.b * op + u.a * b * u.op;
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 > 0 ? 1 : -1);
}

Fraction Fraction::operator - (const Fraction& u) {
    if(op == 0 || u.op == 0) {
        std::cout << "ERROR" << std::endl;
        exit(0);
    }
    ull b1 = b * u.b;
    long a1 = a * u.b * op - u.a * b * u.op;
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 > 0 ? 1 : -1);
}

Fraction Fraction::operator * (const Fraction& u) {
    if(op == 0 || u.op == 0) {
        std::cout << "ERROR" << std::endl;
        exit(0);
    }
    ull&& a1 = a * u.a;
    ull&& b1 = b * u.b;
    return Fraction(a1,b1,op * u.op);
}

Fraction Fraction::operator / (const Fraction& u) {
    if(op == 0 || u.op == 0) {
        std::cout << "ERROR" << std::endl;
        exit(0);
    }
    ull&& a1 = a * u.b;
    ull&& b1 = b * u.a;
    return Fraction(a1,b1,op * u.op);
}
