#ifndef __FRACTION_H__
#define __FRACTION_H__
#include <gperftools/profiler.h>
#include <limits>
#include <type_traits>
#pragma once
namespace caesar {
using namespace std;
template <typename T>
class Fraction {
private:
    int op;
    T numer, denom;

public:
    Fraction();

    Fraction(const T &_numer, const T &_denom); // a/b

    template <typename T1>
    Fraction(const Fraction<T1> &u);

    Fraction(const T &_numer, const T &_denom, const int &_op); // op * a/b

    [[nodiscard]] auto ConvToFloat() const -> double;

    auto get_numer() const -> T;

    auto get_denom() const -> T;

    [[nodiscard]] auto mark() const -> int;

    void Simplify();
};

template <typename T1, typename T2>
auto gcd(const T1 &_x, const T2 &_y) -> typename std::common_type_t<T1, T2> {
    std::common_type_t<T1, T2> x = _x, y = _y;
    if (!x || !y) return x > y ? x : y;
    for (int t; t = x % y; x = y, y = t)
        ;
    return y;
}

template <typename T1, typename T2>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> bool {
    return _y > numeric_limits<typename std::common_type_t<T1, T2>>::max() - _x
           || _y < numeric_limits<typename std::common_type_t<T1, T2>>::min() - _x;
}

template <typename T1, typename T2>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> bool {
    return _y > numeric_limits<typename std::common_type_t<T1, T2>>::max() / (double)_x
           || _y < numeric_limits<typename std::common_type_t<T1, T2>>::min() / (double)_x;
}

template <typename T>
Fraction<T>::Fraction() {
    this->op = 1;
    this->numer = 0;
    this->denom = 1;
}

template <typename T>
Fraction<T>::Fraction(const T &_numer, const T &_denom) {
    this->op = _numer >= 0 ? (_denom > 0 ? 1 : -1) : (_denom > 0 ? -1 : 1);
    this->numer = _numer > 0 ? _numer : -_numer;
    this->denom = _denom > 0 ? _denom : -_denom;
    if (this->denom == 0) throw(-1);
    Simplify();
}

template <typename T>
Fraction<T>::Fraction(const T &_numer, const T &_denom, const int &_op) {
    if (_denom == 0) throw(-1);
    this->numer = _numer > 0 ? _numer : -_numer;
    this->denom = _denom > 0 ? _denom : -_denom;
    this->op = _op;
    Simplify();
}

template <typename T>
template <typename T1>
Fraction<T>::Fraction(const Fraction<T1> &u) {
    this->numer = u.get_numer();
    this->denom = u.get_denom();
    this->op = u.mark();
}

template <typename T>
void Fraction<T>::Simplify() {
    if (this->numer == 0)
        this->denom = 1;
    else {
        T g = gcd(this->numer, this->denom);
        this->numer = this->numer / g;
        this->denom = this->denom / g;
    }
}

template <typename T>
auto Fraction<T>::get_numer() const -> T {
    return this->numer;
}

template <typename T>
auto Fraction<T>::get_denom() const -> T {
    return this->denom;
}

template <typename T>
auto Fraction<T>::mark() const -> int {
    return this->op;
}

template <typename T>
inline auto Fraction<T>::ConvToFloat() const -> double {
    return (double)this->numer / this->denom * op;
}

template <typename T1, typename T2>
auto operator+(const Fraction<T1> &u, const Fraction<T2> &v) {
    ProfilerStart("cpp_demo_perf.prof");
    if (u.mark() * v.mark() == 1) {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)
                || check_add_overflow(u.get_numer() * v1_denom, v.get_numer() * u1_denom))
                throw(-1);
            auto t = u.get_numer() * v1_denom + v.get_numer() * u1_denom;
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, u.mark()};
        } else {
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())
                || check_add_overflow(u.get_numer() * v.get_denom(), v.get_numer() * u.get_denom()))
                throw(-1);
            auto t = u.get_numer() * v.get_denom() + v.get_numer() * u.get_denom();
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, u.mark()};
        }
    } else {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            int op = u.mark();
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom))
                throw(-1);
            auto t = u.get_numer() * v1_denom - v.get_numer() * u1_denom;
            if (t < 0) t = -t, op = -op;
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, op};
        } else {
            int op = u.mark();
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom()))
                throw(-1);
            auto t = u.get_numer() * v.get_denom() - v.get_numer() * u.get_denom();
            if (t < 0) t = -t, op = -op;
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, op};
        }
    }
    ProfilerStop();
}

template <typename T1, typename T2>
auto operator-(const Fraction<T1> &u, const Fraction<T2> &v)
    -> Fraction<typename std::common_type_t<T1, T2>> {
    if (u.mark() * v.mark == 1) {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            int op = u.mark();
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom))
                throw(-1);
            auto t = u.get_numer() * v1_denom - v.get_numer() * u1_denom;
            if (t < 0) t = -t, op = -op;
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, op};
        } else {
            int op = u.mark();
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom()))
                throw(-1);
            auto t = u.get_numer() * v.get_denom() - v.get_numer() * u.get_denom();
            if (t < 0) t = -t, op = -op;
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, op};
        }
    } else {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)
                || check_add_overflow(u.get_numer() * v1_denom, v.get_numer() * u1_denom))
                throw(-1);
            auto t = u.get_numer() * v1_denom + v.get_numer() * u1_denom;
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, u.mark()};
        } else {
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())
                || check_add_overflow(u.get_numer() * v.get_denom(), v.get_numer() * u.get_denom()))
                throw(-1);
            auto t = u.get_numer() * v.get_denom() + v.get_numer() * u.get_denom();
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, u.mark()};
        }
    }
}

template <typename T1, typename T2>
auto operator*(const Fraction<T1> &u, const Fraction<T2> &v)
    -> Fraction<typename std::common_type_t<T1, T2>> {
    auto d1 = gcd(u.get_numer(), v.get_denom());
    auto d2 = gcd(u.get_denom(), v.get_numer());
    auto u_numer = u.get_numer() / d1;
    auto u_denom = u.get_denom() / d2;
    auto v_numer = v.get_numer() / d2;
    auto v_denom = v.get_denom() / d1;
    if (check_mul_overflow(u_numer, v_numer) || check_mul_overflow(u_denom, v_denom)) throw(-1);
    return Fraction<typename std::common_type_t<T1, T2>>{u_numer * v_numer, u_denom * v_denom,
                                                         u.mark() * v.mark()};
}

template <typename T1, typename T2>
auto operator/(const Fraction<T1> &u, const Fraction<T2> &v)
    -> Fraction<typename std::common_type_t<T1, T2>> {
    if (v.get_numer() == 0) throw(-1);
    auto d1 = gcd(u.get_numer(), v.get_numer());
    auto d2 = gcd(u.get_denom(), v.get_denom());
    auto u_numer = u.get_numer() / d1;
    auto u_denom = u.get_denom() / d2;
    auto v_numer = v.get_numer() / d1;
    auto v_denom = v.get_denom() / d2;
    if (check_mul_overflow(u_numer, v_denom) || check_mul_overflow(u_denom, v_numer)) throw(-1);
    return Fraction<typename std::common_type_t<T1, T2>>{u_numer * v_denom, u_denom * v_numer,
                                                         u.mark() * v.mark()};
}

template <typename T1, typename T2>
auto operator<(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.ConvToFloat() - v.ConvToFloat() < 0;
}

template <typename T1, typename T2>
auto operator>(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.ConvToFloat() - v.ConvToFloat() > 0;
}

template <typename T1, typename T2>
auto operator==(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    if (u.mark() == v.mark() && u.get_numer() == v.get_numer() && u.get_denom() == v.get_denom()) return true;
    return false;
}

} // namespace caesar
#endif