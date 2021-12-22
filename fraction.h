#ifndef __FRACTION_H__
#define __FRACTION_H__
#include <iostream>

typedef unsigned long long ull;

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
 
    friend std::istream& operator >> (std::istream& os,Fraction& u); 

    friend std::ostream& operator << (std::ostream& os,const Fraction& u);

    void Simplify();

    ull gcd(ull _x,ull _y);

    Fraction operator * (const Fraction& u);

    Fraction operator / (const Fraction& u);

    Fraction operator + (const Fraction& u);

    Fraction operator - (const Fraction& u);

    bool operator < (const Fraction& u);

    bool operator > (const Fraction& u);
    
    bool operator == (const Fraction& u);
};
#endif