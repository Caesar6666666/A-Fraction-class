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
class Fraction {
private:
    int8_t sign;
    T numer, denom;
    double val;

public:
    Fraction();

    Fraction(const T &_numer, const T &_denom); // a/b

    template <typename T1>
    Fraction(const Fraction<T1> &u);

    Fraction(const T &_numer, const T &_denom, const int8_t &_sign); // op * a/b

    [[nodiscard]] auto ConvToFloat() const -> double;

    auto get_numer() const -> T;

    auto get_denom() const -> T;

    [[nodiscard]] auto get_sign() const -> int8_t;

    void Simplify();
};

// template <typename T1, typename T2>
// auto gcd(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
//     std::common_type_t<T1, T2> x = _x, y = _y;
//     if (!x || !y) return x > y ? x : y;
//     for (uint64_t t; t = x % y; x = y, y = t)
//         ;
//     return y;
// }

template <typename T>
class gcd {
public:
    gcd(T x, T y) : x_(x), y_(y) {
    }

    T operator()() {
        T a = x_, b = y_;
        if (!a || !b) return a > b ? a : b;
        for (T t; t = a % b; a = b, b = t)
            ;
        return b;
    }

private:
    T x_;
    T y_;
};

template <>
class gcd<uint32_t> {
public:
    gcd(uint32_t a, uint32_t b) : a_(a), b_(b) {
    }

    uint32_t operator()() {
        uint32_t _x = a_, _y = b_;
        int az = __builtin_ctz(a_), bz = __builtin_ctz(b_);
        int z = min(az, bz);
        int dif;
        _y >>= bz;
        while (_x) {
            _x >>= az;
            dif = _y - _x;
            az = __builtin_ctz(dif);
            if (_x < _y) _y = _x;
            if (dif < 0)
                _x = -dif;
            else
                _x = dif;
        }
        return _y << z;
    }

private:
    uint32_t a_;
    uint32_t b_;
};


template <>
class gcd<uint64_t> {
public:
    gcd(uint64_t a, uint64_t b) : a_(a), b_(b) {
    }

    uint64_t operator()() {
        uint64_t _x = a_, _y = b_;
        int az = __builtin_ctz(a_), bz = __builtin_ctz(b_);
        int z = min(az, bz);
        int dif;
        _y >>= bz;
        while (_x) {
            _x >>= az;
            dif = _y - _x;
            az = __builtin_ctz(dif);
            if (_x < _y) _y = _x;
            if (dif < 0)
                _x = -dif;
            else
                _x = dif;
        }
        return _y << z;
    }

private:
    uint64_t a_;
    uint64_t b_;
};

auto gcd1(uint64_t _x, uint64_t _y) {
    auto az = __builtin_ctzll(_x), bz = __builtin_ctzll(_y);
    int z = min(az, bz);
    int dif;
    _y >>= bz;
    while (_x) {
        _x >>= az;
        dif = _y - _x;
        az = __builtin_ctzll(dif);
        if (_x < _y) _y = _x;
        if (dif < 0)
            _x = -dif;
        else
            _x = dif;
    }
    return _y << z;
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
    this->sign = 1;
    this->numer = 0;
    this->denom = 1;
    this->val = 0;
}

template <typename T>
Fraction<T>::Fraction(const T &_numer, const T &_denom) {
    if (_denom == 0) {
        cout << "divided by 0" << endl;
        throw("input error");
    }
    this->sign = MathUtils_SignBit(numer) * MathUtils_SignBit(denom);
    this->numer = _numer * MathUtils_SignBit(numer);
    this->denom = _denom * MathUtils_SignBit(denom);
    Simplify();
}

template <typename T>
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

template <typename T>
template <typename T1>
Fraction<T>::Fraction(const Fraction<T1> &u) {
    this->numer = u.get_numer();
    this->denom = u.get_denom();
    this->sign = u.get_sign();
    this->val = u.ConvToFloat();
}

template <typename T>
void Fraction<T>::Simplify() {
    if (this->numer == 0) {
        this->denom = 1;
        this->sign = 1;
        this->val = 0;
    } else {
        auto g = gcd(this->numer, this->denom)();
        this->numer = this->numer / g;
        this->denom = this->denom / g;
        this->val = (double)this->numer / this->denom * this->sign;
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
auto Fraction<T>::get_sign() const -> int8_t {
    return this->sign;
}

template <typename T>
inline auto Fraction<T>::ConvToFloat() const -> double {
    return this->val;
}

template <typename T1, typename T2>
auto operator+(const Fraction<T1> &u, const Fraction<T2> &v) {
    // ProfilerStart("cpp_demo_perf.prof");
    if (u.get_sign() * v.get_sign() == 1) {
        auto d1 = gcd(u.get_denom(), v.get_denom())();
        if (d1 > 1) {
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)
                || check_add_overflow(u.get_numer() * v1_denom, v.get_numer() * u1_denom)) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            auto t = u.get_numer() * v1_denom + v.get_numer() * u1_denom;
            auto d2 = gcd(t, d1)();
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2,
                                                                 u.get_sign()};
        } else {
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())
                || check_add_overflow(u.get_numer() * v.get_denom(), v.get_numer() * u.get_denom())) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            auto t = u.get_numer() * v.get_denom() + v.get_numer() * u.get_denom();
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom)();
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, u.get_sign()};
        }
    } else {
        auto d1 = gcd(u.get_denom(), v.get_denom())();
        if (d1 > 1) {
            int8_t sign = u.get_sign();
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            std::common_type_t<T1, T2> t;
            if (u.get_sign() * u.ConvToFloat() > v.get_sign() * v.ConvToFloat()) {
                t = u.get_numer() * v1_denom - v.get_numer() * u1_denom;
            } else {
                t = v.get_numer() * u1_denom - u.get_numer() * v1_denom;
                sign = -sign;
            }
            auto d2 = gcd(t, d1)();
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, sign};
        } else {
            int8_t sign = u.get_sign();
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            std::common_type_t<T1, T2> t;
            if (u.get_sign() * u.ConvToFloat() > v.get_sign() * v.ConvToFloat()) {
                t = u.get_numer() * v.get_denom() - v.get_numer() * u.get_denom();
            } else {
                t = v.get_numer() * u.get_denom() - u.get_numer() * v.get_denom();
                sign = -sign;
            }
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom)();
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, sign};
        }
    }
    // ProfilerStop();
}

template <typename T1, typename T2>
auto operator-(const Fraction<T1> &u, const Fraction<T2> &v)
    -> Fraction<typename std::common_type_t<T1, T2>> {
    if (u.get_sign() * v.get_sign() == 1) {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            int8_t sign = u.get_sign();
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            std::common_type_t<T1, T2> t;
            if (u.get_sign() * u.ConvToFloat() > v.get_sign() * v.ConvToFloat()) {
                t = u.get_numer() * v1_denom - v.get_numer() * u1_denom;
            } else {
                t = v.get_numer() * u1_denom - u.get_numer() * v1_denom;
                sign = -sign;
            }
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2, sign};
        } else {
            int8_t sign = u.get_sign();
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            std::common_type_t<T1, T2> t;
            if (u.get_sign() * u.ConvToFloat() > v.get_sign() * v.ConvToFloat()) {
                t = u.get_numer() * v.get_denom() - v.get_numer() * u.get_denom();
            } else {
                t = v.get_numer() * u.get_denom() - u.get_numer() * v.get_denom();
                sign = -sign;
            }
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, sign};
        }
    } else {
        auto d1 = gcd(u.get_denom(), v.get_denom());
        if (d1 > 1) {
            auto v1_denom = v.get_denom() / d1;
            auto u1_denom = u.get_denom() / d1;
            if (check_mul_overflow(u.get_numer(), v1_denom) || check_mul_overflow(v.get_numer(), u1_denom)
                || check_add_overflow(u.get_numer() * v1_denom, v.get_numer() * u1_denom)) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            auto t = u.get_numer() * v1_denom + v.get_numer() * u1_denom;
            auto d2 = gcd(t, d1);
            t /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, u1_denom * v.get_denom() / d2,
                                                                 u.get_sign()};
        } else {
            if (check_mul_overflow(u.get_numer(), v.get_denom())
                || check_mul_overflow(v.get_numer(), u.get_denom())
                || check_add_overflow(u.get_numer() * v.get_denom(), v.get_numer() * u.get_denom())) {
                cout << format("overflow at operator '+': u = {}{}/{}, v = {}{}/{}",
                               u.get_sign() == -1 ? '-' : ' ', u.get_numer(), u.get_denom(),
                               v.get_sign() == -1 ? '-' : ' ', v.get_numer(), v.get_denom())
                     << endl;
            }
            auto t = u.get_numer() * v.get_denom() + v.get_numer() * u.get_denom();
            auto denom = u.get_denom() * v.get_denom();
            auto d2 = gcd(t, denom);
            t /= d2;
            denom /= d2;
            return Fraction<typename std::common_type_t<T1, T2>>{t, denom, u.get_sign()};
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
    if (check_mul_overflow(u_numer, v_numer) || check_mul_overflow(u_denom, v_denom)) {
        cout << format("overflow at operator '*': u = {}{}/{}, v = {}{}/{}", u.get_sign() == -1 ? '-' : ' ',
                       u.get_numer(), u.get_denom(), v.get_sign() == -1 ? '-' : ' ', v.get_numer(),
                       v.get_denom())
             << endl;
    }
    return Fraction<typename std::common_type_t<T1, T2>>{u_numer * v_numer, u_denom * v_denom,
                                                         u.get_sign() * v.get_sign()};
}

template <typename T1, typename T2>
auto operator/(const Fraction<T1> &u, const Fraction<T2> &v)
    -> Fraction<typename std::common_type_t<T1, T2>> {
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
    if (check_mul_overflow(u_numer, v_denom) || check_mul_overflow(u_denom, v_numer)) {
        cout << format("overflow at operator '/': u = {}{}/{}, v = {}{}/{}", u.get_sign() == -1 ? '-' : ' ',
                       u.get_numer(), u.get_denom(), v.get_sign() == -1 ? '-' : ' ', v.get_numer(),
                       v.get_denom())
             << endl;
    }
    return Fraction<typename std::common_type_t<T1, T2>>{u_numer * v_denom, u_denom * v_numer,
                                                         u.get_sign() * v.get_sign()};
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
    if (u.get_sign() == v.get_sign() && u.get_numer() == v.get_numer() && u.get_denom() == v.get_denom())
        return true;
    return false;
}

} // namespace caesar
#endif
