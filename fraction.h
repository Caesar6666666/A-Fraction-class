#ifndef __FRACTION_H__
#define __FRACTION_H__
#pragma once

#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <type_traits>
#include <numeric>
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

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator+(const Fraction<T1> &u, const Fraction<T2> &v);

    template <UnsignedType T1>
    friend auto operator-(const Fraction<T1> &u) -> Fraction<T>;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator-(const Fraction<T1> &u, const Fraction<T2> &v)
        -> Fraction<std::common_type_t<T1, T2>>;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator*(const Fraction<T1> &u, const Fraction<T2> &v)
        -> Fraction<std::common_type_t<T1, T2>>;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator/(const Fraction<T1> &u, const Fraction<T2> &v)
        -> Fraction<std::common_type_t<T1, T2>>;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator==(const Fraction<T1> &u, const Fraction<T2> &v) -> bool;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator<(const Fraction<T1> &u, const Fraction<T2> &v) -> bool;

    template <UnsignedType T1, UnsignedType T2>
    friend auto operator>(const Fraction<T1> &u, const Fraction<T2> &v) -> bool;
};

#ifndef __disable_check_overflow__
template <typename T1, typename T2>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    if (_x > std::numeric_limits<std::common_type_t<T1, T2>>::max() - _y) {
        throw std::overflow_error("addition overflow");
    } else {
        return _x + _y;
    }
}

template <typename T1, typename T2>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    if (_x > std::numeric_limits<std::common_type_t<T1, T2>>::max() / _y) {
        throw std::overflow_error("multiplication overflow");
    } else {
        return _x * _y;
    }
}

#ifndef __disable_gcc_builtin__
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
#endif

#else
template <typename T1, typename T2>
inline auto check_add_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    return _x + _y;
}

template <typename T1, typename T2>
inline auto check_mul_overflow(const T1 &_x, const T2 &_y) -> std::common_type_t<T1, T2> {
    return _x * _y;
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
    this->sign = MathUtils_SignBit(_numer) * MathUtils_SignBit(_denom);
    this->numer = _numer > 0 ? _numer : ~(T)_numer + 1;
    this->denom = _denom > 0 ? _denom : ~(T)_denom + 1;
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
    if (u.sign == v.sign) {
        auto d1 = gcd(u.denom, v.denom);
        if (d1 > 1) {
            auto t = check_add_overflow(check_mul_overflow(u.numer, (v.denom / d1)),
                                        check_mul_overflow(v.numer, (u.denom / d1)));
            auto d2 = gcd(t, d1);
            return Fraction<std::common_type_t<T1, T2>>(
                t / d2, check_mul_overflow((u.denom / d1), (v.denom / d2)), u.sign);
        } else {
            auto t = check_add_overflow(check_mul_overflow(u.numer, v.denom),
                                        check_mul_overflow(v.numer, u.denom));
            return Fraction<std::common_type_t<T1, T2>>(t, check_mul_overflow(u.denom, v.denom), u.sign);
        }
    } else {
        auto d1 = gcd(u.denom, v.denom);
        if (d1 > 1) {
            auto u1_numer = check_mul_overflow(u.numer, (v.denom / d1));
            auto v1_numer = check_mul_overflow(v.numer, (u.denom / d1));
            if (u1_numer > v1_numer) {
                auto t = u1_numer - v1_numer;
                auto d2 = gcd(t, d1);
                return Fraction<std::common_type_t<T1, T2>>(
                    t / d2, check_mul_overflow((u.denom / d1), (v.denom / d2)), u.sign);
            } else {
                auto t = v1_numer - u1_numer;
                auto d2 = gcd(t, d1);
                return Fraction<std::common_type_t<T1, T2>>(
                    t / d2, check_mul_overflow((u.denom / d1), (v.denom / d2)), -u.sign);
            }
        } else {
            auto u1_numer = check_mul_overflow(u.numer, v.denom);
            auto v1_numer = check_mul_overflow(v.numer, u.denom);
            if (u1_numer > v1_numer) {
                auto t = u1_numer - v1_numer;
                return Fraction<std::common_type_t<T1, T2>>(t, check_mul_overflow(u.denom, v.denom), u.sign);
            } else {
                auto t = v1_numer - u1_numer;
                return Fraction<std::common_type_t<T1, T2>>(t, check_mul_overflow(u.denom, v.denom), -u.sign);
            }
        }
    }
}

template <UnsignedType T>
auto operator-(const Fraction<T> &u) -> Fraction<T> {
    return Fraction<T>(u.numer, u.denom, -u.sign);
}

template <UnsignedType T1, UnsignedType T2>
auto operator-(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    return u + (-v);
}

template <UnsignedType T1, UnsignedType T2>
auto operator*(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    auto d1 = gcd(u.numer, v.denom);
    auto d2 = gcd(u.denom, v.numer);
    auto u_numer = u.numer / d1;
    auto u_denom = u.denom / d2;
    auto v_numer = v.numer / d2;
    auto v_denom = v.denom / d1;
    return Fraction<typename std::common_type_t<T1, T2>>{
        check_mul_overflow(u_numer, v_numer), check_mul_overflow(u_denom, v_denom), u.sign * v.sign};
}

template <UnsignedType T1, UnsignedType T2>
auto operator/(const Fraction<T1> &u, const Fraction<T2> &v) -> Fraction<std::common_type_t<T1, T2>> {
    if (v.numer == 0) {
        cout << format("overflow at operator '/': u = {}{}/{}, v = {}{}/{}", u.sign == -1 ? '-' : ' ',
                       u.numer, u.denom, v.sign == -1 ? '-' : ' ', v.numer, v.denom)
             << endl;
    }
    auto d1 = gcd(u.numer, v.numer);
    auto d2 = gcd(u.denom, v.denom);
    auto u_numer = u.numer / d1;
    auto u_denom = u.denom / d2;
    auto v_numer = v.numer / d1;
    auto v_denom = v.denom / d2;
    return Fraction<typename std::common_type_t<T1, T2>>{
        check_mul_overflow(u_numer, v_denom), check_mul_overflow(u_denom, v_numer), u.sign * v.sign};
}

template <UnsignedType T1, UnsignedType T2>
auto operator<(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.val - v.val < 0;
}

template <UnsignedType T1, UnsignedType T2>
auto operator>(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    return u.val - v.val > 0;
}

template <UnsignedType T1, UnsignedType T2>
auto operator==(const Fraction<T1> &u, const Fraction<T2> &v) -> bool {
    if (u.sign == v.sign && u.numer == v.numer && u.denom == v.denom) return true;
    return false;
}

} // namespace caesar
#endif
