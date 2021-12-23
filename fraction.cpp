#include "fraction.h"

caesar::Fraction::Fraction() {
    this->op = 1;
    this->a = 0;
    this->b = 1;
}

caesar::Fraction::Fraction(long long _a, long long _b) {
    this->op = _a >= 0 ? (_b > 0 ? 1 : -1) : (_b > 0 ? -1 : 1);
    this->a = _a > 0 ? _a : -_a;
    this->b = _b > 0 ? _b : -_b;
    if(this->b == 0) throw(-1);
    Simplify();
}

caesar::Fraction::Fraction(ull _a, ull _b, int _op):a(_a),b(_b),op(_op) {
    if(_b == 0) throw(-1);
    Simplify();
}

caesar::Fraction::Fraction(const Fraction& u) {
    this->a = u.a;
    this->b = u.b;
    this->op = u.op;
}

void caesar::Fraction::Simplify() {
    ull g = gcd(this->a, this->b);
    this->a = this->a / g;
    this->b = this->b / g;
}

ull caesar::Fraction::num() const {
    return this->a;
}

ull caesar::Fraction::deno() const {
    return this->b;
}

int caesar::Fraction::mark() const {
    return this->op;
}

double caesar::Fraction::transform_to_float() const {
    return (double) this->a / this->b * op;
}

template<typename T>
T caesar::gcd(T _x,T _y) {
    return _y == 0 ? _x : gcd(_y, _x % _y);
}

caesar::Fraction caesar::Fraction::operator + (const Fraction& u) const {
    ull b1 = this->b * u.num();
    long a1 = this->a * u.deno() * this->op + u.num() * this->b * u.mark();
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 >= 0 ? 1 : -1);
}

caesar::Fraction caesar::Fraction::operator - (const Fraction& u) const {
    ull b1 = this->b * u.deno();
    long a1 = this->a * u.deno() * this->op - u.num() * this->b * u.mark();
    return Fraction(a1 > 0 ? a1 : -a1, b1, a1 >= 0 ? 1 : -1);
}

caesar::Fraction caesar::Fraction::operator * (const Fraction& u) const {
    ull&& m1 = gcd(this->a,u.deno());
    ull&& m2 = gcd(this->b,u.num());
    ull&& a1 = this->a / m1 * u.num() / m2;
    ull&& b1 = this->b / m2 * u.deno() / m1;
    return Fraction(a1,b1,this->op * u.mark());
}

caesar::Fraction caesar::Fraction::operator / (const Fraction& u) const {
    if(u.num() == 0) throw(-1);
    ull&& m1 = gcd(this->a,u.num());
    ull&& m2 = gcd(this->b,u.deno());
    ull&& a1 = this->a / m1 * u.deno() / m2;
    ull&& b1 = this->b / m2 * u.num() / m1;
    return Fraction(a1,b1,this->op * u.mark());
}

bool caesar::Fraction::operator < (const Fraction& u) const {
    if(this->op == u.mark() && this->a == u.num() && this->b == u.deno()) return false;
    return this->transform_to_float() - u.transform_to_float() < 0;
}

bool caesar::Fraction::operator > (const Fraction& u) const {
   if(this->op == u.mark() && this->a == u.num() && this->b == u.deno()) return false;
    return this->transform_to_float() - u.transform_to_float() > 0;
}

bool caesar::Fraction::operator == (const Fraction& u) const {
    if(this->op == u.mark() && this->a == u.num() && this->b == u.deno()) return true;
    return false;
}