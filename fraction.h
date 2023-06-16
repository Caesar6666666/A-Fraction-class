#ifndef __FRACTION_H__
#define __FRACTION_H__
#pragma once

// #include <gperftools/profiler.h>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <type_traits>

#define MathUtils_SignBit(x) (((signed char *)&x)[sizeof(x) - 1] >> 7 | 1)
namespace caesar {
using namespace std;
template <typename T>
concept UnsignedType = std::is_unsigned_v<T>;

template <typename T>
concept UnsignedInt = UnsignedType<T> && std::is_same_v<T, unsigned int>;

template <typename T>
concept UnsignedLong = UnsignedType<T> && std::is_same_v<T, unsigned long>;

template <typename T>
concept UnsignedLongLong = UnsignedType<T> && std::is_same_v<T, unsigned long long>;
template <UnsignedType T>
class Fraction {
private:
    int8_t sign;
    T numer, denom;
    double val;

public:
    Fraction();

    template <typename T1>
    Fraction(const T1 &_numer, const T1 &_denom); // a/b

    template <UnsignedType T1>
        requires(sizeof(T1) <= sizeof(T))
    Fraction(const Fraction<T1> &u);

    Fraction(const T &_numer, const T &_denom, const int8_t &_sign); // op * a/b

    auto ConvToFloat() const -> double;

    auto get_numer() const -> T;

    auto get_denom() const -> T;

    auto get_sign() const -> int8_t;

    void Simplify();
};


template <UnsignedType T1, UnsignedType T2>
auto gcd(T1 a, T2 b) -> std::common_type_t<T1, T2> {
    if (!a || !b) return a > b ? a : b;
    for (std::common_type_t<T1, T2> t; t = a % b; a = b, b = t)
        ;
    return b;
}

#ifdef __GNUC__
// https://www.luogu.com.cn/blog/maysoul/solution-p5435
template <UnsignedType T1, UnsignedType T2>
    requires UnsignedInt<std::common_type_t<T1, T2>>
auto gcd(T1 a, T2 b) -> std::common_type_t<T1, T2> {
    int az = __builtin_ctz(a), bz = __builtin_ctz(b);
    int z = min(az, bz);
    int dif;
    b >>= bz;
    while (a) {
        a >>= az;
        dif = b - a;
        az = __builtin_ctz(dif);
        if (a < b) b = a;
        if (dif < 0)
            a = -dif;
        else
            a = dif;
    }
    return b << z;
}

template <UnsignedType T1, UnsignedType T2>
    requires UnsignedLong<std::common_type_t<T1, T2>>
auto gcd(T1 a, T2 b) -> std::common_type_t<T1, T2> {
    int az = __builtin_ctzl(a), bz = __builtin_ctzl(b);
    int z = min(az, bz);
    long dif;
    b >>= bz;
    while (a) {
        a >>= az;
        dif = b - a;
        az = __builtin_ctzl(dif);
        if (a < b) b = a;
        if (dif < 0)
            a = -dif;
        else
            a = dif;
    }
    return b << z;
}

template <UnsignedType T1, UnsignedType T2>
    requires UnsignedLongLong<std::common_type_t<T1, T2>>
auto gcd(T1 a, T2 b) -> std::common_type_t<T1, T2> {
    int az = __builtin_ctzll(a), bz = __builtin_ctzll(b);
    int z = min(az, bz);
    int64_t dif;
    b >>= bz;
    while (a) {
        a >>= az;
        dif = b - a;
        az = __builtin_ctzll(dif);
        if (a < b) b = a;
        if (dif < 0)
            a = -dif;
        else
            a = dif;
    }
    return b << z;
}
#endif

template <typename T1, typename T2>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    if(_x > std::numeric_limits<std::common_type_t<T1, T2>>::max() - _y) {
        throw std::overflow_error("addition overflow");
    } else {
        return _x + _y;
    }
}

template <typename T1, typename T2>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    if(_x > std::numeric_limits<std::common_type_t<T1, T2>>::max() / _y) {
        throw std::overflow_error("multiplication overflow");
    } else {
        return _x * _y;
    }
}


#ifdef __GNUC__
template <typename T1, typename T2>
    requires UnsignedInt<std::common_type_t<T1, T2>>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned int ans;
    __builtin_uadd_overflow(_x, _y, &ans);
    return ans;
}

template <typename T1, typename T2>
    requires UnsignedLong<std::common_type_t<T1, T2>>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned long ans;
    __builtin_uaddl_overflow(_x, _y, &ans);
    return ans;
}

template <typename T1, typename T2>
    requires UnsignedLongLong<std::common_type_t<T1, T2>>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned long long ans;
    __builtin_uaddll_overflow(_x, _y, &ans);
    return ans;
}

template <typename T1, typename T2>
    requires UnsignedInt<std::common_type_t<T1, T2>>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned int ans;
    __builtin_umul_overflow(_x, _y, &ans);
    return ans;
}

template <typename T1, typename T2>
    requires UnsignedLong<std::common_type_t<T1, T2>>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned long ans;
    __builtin_umull_overflow(_x, _y, &ans);
    return ans;
}

template <typename T1, typename T2>
    requires UnsignedLongLong<std::common_type_t<T1, T2>>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    unsigned long long ans;
    __builtin_umulll_overflow(_x, _y, &ans);
    return ans;
}
#endif

template <UnsignedType T>
Fraction<T>::Fraction() {
    this->sign = 1;
    this->numer = 0;
    this->denom = 1;
    this->val = 0;
}

template <UnsignedType T>
template <typename T1>
Fraction<T>::Fraction(const T1 &_numer, const T1 &_denom) {
    if (_denom == 0) {
        cout << "divided by 0" << endl;
        throw("input error");
    }
    this->sign = MathUtils_SignBit(numer) * MathUtils_SignBit(denom);
    this->numer = _numer * MathUtils_SignBit(numer);
    this->denom = _denom * MathUtils_SignBit(denom);
    Simplify();
}

template <UnsignedType T>
Fraction<T>::Fraction(const T &_numer, const T &_denom, const int8_t &_sign) {
    if ((_sign != 1 && _sign != -1) || _numer < 0 || _denom <= 0) {
        cout << format("input error: op = {}, numer = {}, denom = {}", _sign, _numer, _denom) << endl;
        throw("input error");
    }
    this->sign = _sign;
    this->numer = _numer;
    this->denom = _denom;
    Simplify();
}

template <UnsignedType T>
template <UnsignedType T1>
    requires(sizeof(T1) <= sizeof(T))
Fraction<T>::Fraction(const Fraction<T1> &u) {
    this->numer = u.get_numer();
    this->denom = u.get_denom();
    this->sign = u.get_sign();
    this->val = u.ConvToFloat();
}

template <UnsignedType T>
void Fraction<T>::Simplify() {
    if (this->numer == 0) {
        this->denom = 1;
        this->sign = 1;
        this->val = 0;
    } else {
        auto g = gcd(this->numer, this->denom);
        this->numer = this->numer / g;
        this->denom = this->denom / g;
        this->val = (double)this->numer / this->denom * this->sign;
    }
}

template <UnsignedType T>
auto Fraction<T>::get_numer() const -> T {
    return this->numer;
}

template <UnsignedType T>
auto Fraction<T>::get_denom() const -> T {
    return this->denom;
}

template <UnsignedType T>
auto Fraction<T>::get_sign() const -> int8_t {
    return this->sign;
}

template <UnsignedType T>
inline auto Fraction<T>::ConvToFloat() const -> double {
    return this->val;
}

template <UnsignedType T1, UnsignedType T2>
auto operator+(const Fraction<T1> &u, const Fraction<T2> &v) {
    // ProfilerStart("cpp_demo_perf.prof");
    if (u.get_sign() == v.get_sign()) {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            auto t = check_add_overflow(check_mul_overflow(u.get_numer(), (v.get_denom() / d1)),
                                        check_mul_overflow(v.get_numer(), (u.get_denom() / d1)));
            auto d2 = gcd(t, d1);
            return Fraction<std::common_type_t<T1, T2>>(
                t / d2, check_mul_overflow((u.get_denom() / d1), (v.get_denom() / d2)), u.get_sign());
        } else {
            auto t = check_add_overflow(check_mul_overflow(u.get_numer(), v.get_denom()),
                                        check_mul_overflow(v.get_numer(), u.get_denom()));
            return Fraction<std::common_type_t<T1, T2>>(t, check_mul_overflow(u.get_denom(), v.get_denom()),
                                                        u.get_sign());
        }
    } else {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            auto u1_numer = check_mul_overflow(u.get_numer(), (v.get_denom() / d1));
            auto v1_numer = check_mul_overflow(v.get_numer(), (u.get_denom() / d1));
            if (u1_numer > v1_numer) {
                auto t = u1_numer - v1_numer;
                auto d2 = gcd(t, d1);
                return Fraction<std::common_type_t<T1, T2>>(
                    t / d2, check_mul_overflow((u.get_denom() / d1), (v.get_denom() / d2)), u.get_sign());
            } else {
                auto t = v1_numer - u1_numer;
                auto d2 = gcd(t, d1);
                return Fraction<std::common_type_t<T1, T2>>(
                    t / d2, check_mul_overflow((u.get_denom() / d1), (v.get_denom() / d2)), -u.get_sign());
            }
        } else {
            auto u1_numer = check_mul_overflow(u.get_numer(), v.get_denom());
            auto v1_numer = check_mul_overflow(v.get_numer(), u.get_denom());
            if (u1_numer > v1_numer) {
                auto t = u1_numer - v1_numer;
                return Fraction<std::common_type_t<T1, T2>>(
                    t, check_mul_overflow(u.get_denom(), v.get_denom()), u.get_sign());
            } else {
                auto t = v1_numer - u1_numer;
                return Fraction<std::common_type_t<T1, T2>>(
                    t, check_mul_overflow(u.get_denom(), v.get_denom()), -u.get_sign());
            }
        }
    }
    // ProfilerStop();
}

template <UnsignedType T>
auto operator-(const Fraction<T> &u) -> Fraction<T> {
    return Fraction<T>(u.get_numer(), u.get_denom(), -u.get_sign());
}

template <UnsignedType T1, UnsignedType T2>
auto operator-(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    return u + (-v);
}

template <UnsignedType T1, UnsignedType T2>
auto operator*(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    auto d1 = gcd(u.get_numer(), v.get_denom());
    auto d2 = gcd(u.get_denom(), v.get_numer());
    auto u_numer = u.get_numer() / d1;
    auto u_denom = u.get_denom() / d2;
    auto v_numer = v.get_numer() / d2;
    auto v_denom = v.get_denom() / d1;
    return Fraction<typename std::common_type_t<T1, T2>>{check_mul_overflow(u_numer, v_numer),
                                                         check_mul_overflow(u_denom, v_denom),
                                                         u.get_sign() * v.get_sign()};
}

template <UnsignedType T1, UnsignedType T2>
auto operator/(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    if (v.get_numer() == 0) {
        cout << format("overflow at operator '/': u = {}{}/{}, v = {}{}/{}", u.get_sign() == -1 ? '-' : ' ',
                       u.get_numer(), u.get_denom(), v.get_sign() == -1 ? '-' : ' ', v.get_numer(),
                       v.get_denom())
             << endl;
    }
    auto d1 = gcd(u.get_numer(), v.get_numer());
    auto d2 = gcd(u.get_denom(), v.get_denom());
    auto u_numer = u.get_numer() / d1;
    auto u_denom = u.get_denom() / d2;
    auto v_numer = v.get_numer() / d1;
    auto v_denom = v.get_denom() / d2;
    return Fraction<typename std::common_type_t<T1, T2>>{check_mul_overflow(u_numer, v_denom),
                                                         check_mul_overflow(u_denom, v_numer),
                                                         u.get_sign() * v.get_sign()};
}

template <UnsignedType T1, UnsignedType T2>
auto operator<(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.ConvToFloat() - v.ConvToFloat() < 0;
}

template <UnsignedType T1, UnsignedType T2>
auto operator>(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.ConvToFloat() - v.ConvToFloat() > 0;
}

template <UnsignedType T1, UnsignedType T2>
auto operator==(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    if (u.get_sign() == v.get_sign() && u.get_numer() == v.get_numer() && u.get_denom() == v.get_denom())
        return true;
    return false;
}

} // namespace caesar
#endif
