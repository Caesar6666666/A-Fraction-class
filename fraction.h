#ifndef __FRACTION_H__
#define __FRACTION_H__
#include <iostream>

typedef unsigned long long ull;

class Fraction {
private:
    int op;
    unsigned long long a,b;

public:
    Fraction(long a,long b);

    Fraction(ull _a,ull _b,int _op);

    double change_to_float();
 
    friend std::istream& operator >> (std::istream& os,Fraction& u); 

    void Simplify();

    ull gcd(ull _x,ull _y);

    Fraction operator * (const Fraction& u);

    Fraction operator / (const Fraction& u);

    Fraction operator + (const Fraction& u);

    Fraction operator - (const Fraction& u);

    friend std::ostream& operator << (std::ostream& os,const Fraction& u);
};
#endif