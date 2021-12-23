#ifndef __FRACTION_H__
#define __FRACTION_H__

typedef unsigned long long ull;
namespace caesar {

class Fraction {
private:
    int op;
    unsigned long long a,b;

public:
    Fraction();

    Fraction(long long a,long long b);

    Fraction(const Fraction& u);

    Fraction(ull _a,ull _b,int _op);

    double transform_to_float() const;

    ull num() const;

    ull deno() const;

    int mark() const;

    void Simplify();

    Fraction operator * (const Fraction& u) const;

    Fraction operator / (const Fraction& u) const;

    Fraction operator + (const Fraction& u) const;

    Fraction operator - (const Fraction& u) const;

    bool operator < (const Fraction& u) const;

    bool operator > (const Fraction& u) const;

    bool operator == (const Fraction& u) const;
};

template<typename T>
T gcd(T _x,T _y);

}
#endif