#ifndef __FRACTION_H__
#define __FRACTION_H__

namespace caesar {

    template<typename T>
    class Fraction {
    private:
        int op;
        T a,b;

    public:
        Fraction();

        Fraction(T _a,T _b);   // a/b

        template<typename T1>
        Fraction(const Fraction<T1>& u);

        Fraction(T _a,T _b,int _op);    //op * a/b

        double ConvToFloat() const;

        T num() const;

        T deno() const;

        int mark() const;

        void Simplify();
    };

    template<typename T1,typename T2>
    T1 gcd(const T1& _x,const T2& _y) {
        return _y == 0 ? _x : gcd(_y, _x % _y);
    }

    template<typename T>
    Fraction<T>::Fraction() {
        this->op = 1;
        this->a = 0;
        this->b = 1;
    }

    template<typename T>
    Fraction<T>::Fraction(T _a, T _b) {
        this->op = _a >= 0 ? (_b > 0 ? 1 : -1) : (_b > 0 ? -1 : 1);
        this->a = _a > 0 ? _a : -_a;
        this->b = _b > 0 ? _b : -_b;
        if(this->b == 0) throw(-1);
        Simplify();
    }

    template<typename T>
    Fraction<T>::Fraction(T _a, T _b, int _op) {
        if(_b == 0) throw(-1);
        this->a = _a > 0 ? _a : -_a;
        this->b = _b > 0 ? _b : -_b;
        this->op = _op;
        Simplify();
    }

    template<typename T>
    template<typename T1>
    Fraction<T>::Fraction(const Fraction<T1>& u) {
        this->a = u.num();
        this->b = u.deno();
        this->op = u.mark();
    }

    template<typename T>
    void Fraction<T>::Simplify() {
        T g = gcd<T>(this->a, this->b);
        this->a = this->a / g;
        this->b = this->b / g;
    }

    template<typename T>
    T Fraction<T>::num() const {
        return this->a;
    }

    template<typename T>
    T Fraction<T>::deno() const {
        return this->b;
    }

    template<typename T>
    int Fraction<T>::mark() const {
        return this->op;
    }

    template<typename T>
    double Fraction<T>::ConvToFloat() const {
        return (double) this->a / this->b * op;
    }

    template<typename T1,typename T2>
    Fraction<typename std::common_type_t<T1, T2>> operator + (const Fraction<T1>& u,const Fraction<T2>& v) {
        T1&& b1 = u.deno() * v.deno();
        T1&& a1 = u.num() * v.deno() * u.mark() + v.num() * u.deno() * v.mark();
        return {a1 > 0 ? a1 : -a1, b1, a1 >= 0 ? 1 : -1};
    }

    template<typename T1,typename T2>
    Fraction<typename std::common_type_t<T1, T2>> operator - (const Fraction<T1>& u,const Fraction<T2>& v) {
        T1&& b1 = u.deno() * v.deno();
        T1&& a1 = u.num() * v.deno() * u.mark() - v.num() * u.deno() * v.mark();
        return {a1 > 0 ? a1 : -a1, b1, a1 >= 0 ? 1 : -1};
    }

    template<typename T1,typename T2>
    Fraction<typename std::common_type_t<T1, T2>> operator * (const Fraction<T1>& u,const Fraction<T2>& v) {
        T1&& m1 = gcd<T1>(u.num(),v.deno());
        T1&& m2 = gcd<T1>(u.deno(),v.num());
        T1&& a1 = u.num() / m1 * v.num() / m2;
        T1&& b1 = u.deno() / m2 * v.deno() / m1;
        return {a1,b1,u.mark() * v.mark()};
    }

    template<typename T1,typename T2>
    Fraction<typename std::common_type_t<T1, T2>> operator / (const Fraction<T1>& u,const Fraction<T2>& v) {
        if(v.num() == 0) throw(-1);
        T1&& m1 = gcd(u.num(),v.num());
        T1&& m2 = gcd(u.deno(),v.deno());
        T1&& a1 = u.num() / m1 * v.deno() / m2;
        T1&& b1 = u.deno() / m2 * v.num() / m1;
        return {a1,b1,u.mark() * v.mark()};
    }

    template<typename T1,typename T2>
    bool operator < (const Fraction<T1>& u,const Fraction<T2>& v) {
        return u.ConvToFloat() - v.ConvToFloat() < 0;
    }

    template<typename T1,typename T2>
    bool operator > (const Fraction<T1>& u,const Fraction<T2>& v) {
        return u.ConvToFloat() - v.ConvToFloat() > 0;
    }

    template<typename T1,typename T2>
    bool operator == (const Fraction<T1>& u,const Fraction<T2>& v) {
        if(u.mark() == v.mark() && u.num() == v.num() && u.deno() == v.deno()) return true;
        return false;
    }

}
#endif