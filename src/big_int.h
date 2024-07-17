#pragma once

#include <cmath>
#include <bit>
#include <ranges>
#include <random>
#include <compare>
#include <concepts>
#include <numeric>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <numbers>
#include <optional>

#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
#include <intrin.h>
#pragma intrinsic(_umul128)
#pragma intrinsic(_udiv128)
#define BIG_INT_HAS_UMUL128_UDIV128
#endif

#ifdef __SIZEOF_INT128__
#if (defined(__clang__) && !defined(_WIN32)) || \
    (defined(__GNUC__) && !defined(__clang__) && !defined(__CUDACC__))
#define BIG_INT_HAVE_INTRINSIC_INT128
#endif
#endif

#ifdef _MSC_VER
#define BIG_INT_FORCE_INLINE __forceinline
#elif defined(__clang__) || defined(__GNUC__)
#define BIG_INT_FORCE_INLINE inline __attribute__((always_inline))
#else
#define BIG_INT_FORCE_INLINE inline
#endif

namespace spaceless {

namespace detail {
static constexpr BIG_INT_FORCE_INLINE uint64_t umul128(uint64_t lhs, uint64_t rhs, uint64_t *high_product);
static constexpr BIG_INT_FORCE_INLINE uint64_t udiv128(uint64_t high_dividend, uint64_t low_dividend, uint64_t divisor, uint64_t *remainder);
static constexpr BIG_INT_FORCE_INLINE uint64_t high_32(uint64_t x) noexcept { return x >> 32; }
static constexpr BIG_INT_FORCE_INLINE uint64_t low_32(uint64_t x) noexcept { return x & ~uint32_t(0); }
static constexpr BIG_INT_FORCE_INLINE uint64_t high_n(uint64_t x, int n) noexcept { return n ? x >> (64 - n) : 0; }
static constexpr BIG_INT_FORCE_INLINE uint64_t low_n(uint64_t x, int n) noexcept { return n ? x << (64 - n) >> (64 - n) : 0; }
static constexpr BIG_INT_FORCE_INLINE void add(uint64_t *lhs, uint64_t rhs, uint64_t *carry) noexcept { *lhs += rhs; *carry += (*lhs < rhs); }
static constexpr BIG_INT_FORCE_INLINE void minus(uint64_t *lhs, uint64_t rhs, uint64_t *carry) noexcept { *carry += (*lhs < rhs); *lhs -= rhs; }
static constexpr BIG_INT_FORCE_INLINE std::strong_ordering flip(const std::strong_ordering &order) noexcept {
	if (order == std::strong_ordering::greater) return std::strong_ordering::less;
	if (order == std::strong_ordering::less) return std::strong_ordering::greater;
	return order; // std::strong_ordering::equal
}
}

class big_int final {
public:
	// from/to string
	enum class base : uint32_t {
		auto_detect = 0,
		binary = 2,
		octal = 8,
		hex = 16,
		decimal = 10,
	};
	enum class to_string_option : uint32_t {
		none = 0,
		showbase = 1 << 0,
		uppercase = 1 << 1,
		showpos = 1 << 2,
	};
	static constexpr std::optional<big_int> from_string(std::string_view sv, base *b_output = nullptr);
	static constexpr std::optional<big_int> from_string(std::string_view sv, base b) { return from_string_impl(sv.data(), sv.length(), &b); }
	constexpr std::string to_string(base b = base::decimal, to_string_option option = to_string_option::showbase) const { return to_string_impl(uint32_t(b), option); }
	// constructors
	constexpr big_int() : data_({ 0 }) {}
	constexpr big_int(const big_int &) = default;
	constexpr big_int(big_int &&) noexcept = default;
	template<std::signed_integral T>
	constexpr big_int(T num) : big_int(int64_t(num)) {}
	template<std::unsigned_integral T>
	constexpr big_int(T num) : big_int(uint64_t(num)) {}
	constexpr big_int(uint64_t num) : data_({ num }) {}
	constexpr big_int(int64_t num) : big_int(num >= 0 ? uint64_t(num) : (uint64_t(~num) + 1)) { sign_ = num < 0; }
	constexpr big_int(std::vector<uint64_t> d, bool sign = false) noexcept : sign_(sign), data_(std::move(d)) { trim_leading_zeros(); if (is_zero(data_)) sign_ = false; }
	template<size_t E>
	constexpr big_int(std::span<const uint64_t, E> d, bool sign = false) : big_int(std::vector(d.begin(), d.end()), sign) {}
	explicit constexpr big_int(std::string_view sv) : big_int(from_string(sv).value()) {}
	constexpr ~big_int() noexcept = default;
	// assignments
	constexpr big_int &operator = (const big_int &rhs) = default;
	constexpr big_int &operator = (big_int &&) noexcept = default;
	// unary operators, also see sign-related functions
	constexpr const big_int &operator + () const & { return *this; }
	constexpr big_int operator + () && noexcept { return std::move(*this); }
	constexpr big_int operator - () const & { return negation(); }
	constexpr big_int operator - () && noexcept { return std::move(*this).negation(); }
	// binary operators
	friend constexpr big_int operator + (const big_int &lhs, const big_int &rhs) { return lhs.copy_for_add_minus(rhs) + rhs; }
	friend constexpr big_int operator + (const big_int &lhs, big_int &&rhs) { return std::move(rhs) + lhs; }
	friend constexpr big_int operator + (big_int &&lhs, const big_int &rhs) { return std::move(lhs += rhs); }
	friend constexpr big_int operator + (big_int &&lhs, big_int &&rhs) { return std::move(lhs += std::move(rhs)); }
	friend constexpr big_int operator - (const big_int &lhs, const big_int &rhs) { return lhs.copy_for_add_minus(rhs) - rhs; }
	friend constexpr big_int operator - (const big_int &lhs, big_int &&rhs) { return std::move(rhs).negation() + lhs; }
	friend constexpr big_int operator - (big_int &&lhs, const big_int &rhs) { return std::move(lhs -= rhs); }
	friend constexpr big_int operator - (big_int &&lhs, big_int &&rhs) { return std::move(lhs -= std::move(rhs)); }
	friend constexpr big_int operator * (const big_int &lhs, const big_int &rhs) { return { multiply_impl(lhs.data_, rhs.data_), lhs.sign_ != rhs.sign_ }; }
	friend constexpr big_int operator * (const big_int &lhs, big_int &&rhs) { return std::move(rhs *= lhs); }
	friend constexpr big_int operator * (big_int &&lhs, const big_int &rhs) { return std::move(lhs *= rhs); }
	friend constexpr big_int operator * (big_int &&lhs, big_int &&rhs) { return std::move(lhs *= std::move(rhs)); }
	friend constexpr big_int operator / (const big_int &lhs, const big_int &rhs) { return lhs.copy_with_reserve(lhs.data_.size() + 1) / rhs; }
	friend constexpr big_int operator / (big_int &&lhs, const big_int &rhs) { return std::move(lhs).divide_with_reminder(rhs); }
	friend constexpr big_int operator % (const big_int &lhs, const big_int &rhs) { return lhs.copy_with_reserve(lhs.data_.size() + 1) % rhs; }
	friend constexpr big_int operator % (big_int &&lhs, const big_int &rhs) { big_int ret; std::move(lhs).divide_with_reminder(rhs, &ret); return ret; }
	// bitwise binary operators
	constexpr big_int operator ~ () const & { return ~big_int(*this); }
	constexpr big_int operator ~ () && { bitwise_not(); return std::move(*this); }
#ifdef __RESHARPER__
	friend constexpr big_int operator & (big_int, big_int) noexcept { return {}; } constexpr big_int &operator &= (const big_int &) { return *this; }
	friend constexpr big_int operator | (big_int, big_int) noexcept { return {}; } constexpr big_int &operator |= (const big_int &) { return *this; }
	friend constexpr big_int operator ^ (big_int, big_int) noexcept { return {}; } constexpr big_int &operator ^= (const big_int &) { return *this; }
#else
#define BIG_INT_BITWISE_OPERATORS(op) \
	friend constexpr big_int operator op (const big_int &lhs, const big_int &rhs) { return std::move(lhs.data_.size() >= rhs.data_.size() ? big_int(lhs) op= rhs : big_int(rhs) op= lhs); } \
	friend constexpr big_int operator op (const big_int &lhs, big_int &&rhs) { return std::move(rhs op= lhs); } \
	friend constexpr big_int operator op (big_int &&lhs, const big_int &rhs) { return std::move(lhs op= rhs); } \
	friend constexpr big_int operator op (big_int &&lhs, big_int &&rhs) noexcept { return std::move(lhs.data_.size() >= rhs.data_.size() ? lhs op= rhs : rhs op= lhs); } \
	constexpr big_int &operator op= (const big_int &rhs) { bitwise_operation_impl(rhs, [](uint64_t &x, uint64_t y) noexcept { x op= y; }); return *this; }
	BIG_INT_BITWISE_OPERATORS(&)
	BIG_INT_BITWISE_OPERATORS(|)
	BIG_INT_BITWISE_OPERATORS(^)
#undef BIG_INT_BITWISE_OPERATORS
#endif
	constexpr big_int operator << (uint64_t rhs) const & { return copy_with_reserve(data_.size() + bits_to_digits(rhs)) << rhs; }
	constexpr big_int operator << (uint64_t rhs) && { return std::move(*this <<= rhs); }
	constexpr big_int operator >> (uint64_t rhs) const & { return big_int(*this) >> rhs; }
	constexpr big_int operator >> (uint64_t rhs) && noexcept { return std::move(*this >>= rhs); }
	// assignment operators
	constexpr big_int &operator += (const big_int &rhs) { add_minus_impl(rhs, false); return *this; }
	constexpr big_int &operator += (big_int &&rhs) { return data_.capacity() >= rhs.data_.capacity() ? *this += rhs : *this = std::move(rhs += *this); }
	constexpr big_int &operator -= (const big_int &rhs) { add_minus_impl(rhs, true); return *this; }
	constexpr big_int &operator -= (big_int &&rhs) { return data_.capacity() >= rhs.data_.capacity() ? *this -= rhs : *this = std::move(rhs -= *this).negation(); }
	constexpr big_int &operator *= (const big_int &rhs);
	constexpr big_int &operator *= (big_int &&rhs) { return data_.capacity() >= rhs.data_.capacity() ? *this *= rhs : *this = std::move(rhs *= *this); }
	constexpr big_int &operator /= (const big_int &rhs) { return *this = std::move(*this) / rhs; }
	constexpr big_int &operator %= (const big_int &rhs) { return *this = std::move(*this) % rhs; }
	constexpr big_int &operator <<= (uint64_t shift) { left_shift_impl(&data_, shift); return *this; }
	constexpr big_int &operator >>= (uint64_t shift) noexcept { right_shift_impl(&data_, shift); return *this; }
	// increment and decrement operators
	constexpr big_int &operator ++ () { return *this += 1; }
	constexpr big_int &operator -- () { return *this -= 1; }
	[[nodiscard]] constexpr big_int operator ++ (int) { auto ret = *this; ++*this; return ret; }
	[[nodiscard]] constexpr big_int operator -- (int) { auto ret = *this; --*this; return ret; }
	// relational operators
	friend constexpr std::strong_ordering operator <=> (const big_int &lhs, const big_int &rhs) noexcept {
		if (lhs.sign_ != rhs.sign_) return lhs.sign_ ? std::strong_ordering::less : std::strong_ordering::greater;
		auto ret = compare_data(lhs.data_, rhs.data_);
		return lhs.sign_ ? detail::flip(ret) : ret;
	}
	template<std::integral T>
	friend constexpr std::strong_ordering operator <=> (const big_int &lhs, T t) noexcept {
		if (t < 0 != lhs.sign_ || lhs.data_.size() > 1) return lhs.sign_ ? std::strong_ordering::less : std::strong_ordering::greater;
		return T(lhs) <=> t;
	}
	template<std::integral T> friend constexpr std::strong_ordering operator <=> (T t, const big_int &rhs) noexcept { return detail::flip(rhs <=> t); }
	template<std::integral T> friend constexpr bool operator == (const big_int &lhs, T t) noexcept { return (lhs <=> t) == 0; }
	template<std::integral T> friend constexpr bool operator == (T t, const big_int &rhs) noexcept { return (rhs <=> t) == 0; }
	constexpr bool operator == (const big_int &rhs) const noexcept = default;
	// bitwise functions
	constexpr uint64_t count_trailing_zeros() const;
	constexpr uint64_t popcount() const noexcept { return std::ranges::fold_left(data_, 0ull, [](auto sum, auto d) noexcept { return sum + std::popcount(d); }); }
	constexpr uint64_t highest_bit() const noexcept { return data_.size() * 64 - std::countl_zero(data_.back()); }
	constexpr void bitwise_not() noexcept { for (auto &x : data_) x = ~x; }
	// sign-related functions
	constexpr BIG_INT_FORCE_INLINE bool sign() const noexcept { return sign_; }
	constexpr BIG_INT_FORCE_INLINE void set_sign(bool s) noexcept { if (*this) sign_ = s; }
	constexpr BIG_INT_FORCE_INLINE void absolutize() noexcept { sign_ = false; }
	[[nodiscard]] constexpr big_int absolute() const & { return big_int(*this).absolute(); }
	[[nodiscard]] constexpr big_int absolute() && noexcept { absolutize(); return std::move(*this); }
	constexpr BIG_INT_FORCE_INLINE void negate() noexcept { set_sign(!sign_); }
	[[nodiscard]] constexpr big_int negation() const & { return big_int(*this).negation(); }
	[[nodiscard]] constexpr big_int negation() && noexcept { negate(); return std::move(*this); }
	// arithmetic functions
	[[nodiscard]] constexpr big_int power(uint64_t n) const;
	static constexpr big_int factorial(uint64_t n);
	static constexpr big_int gcd(big_int x, big_int y);
	static constexpr big_int lcm(const big_int &x, const big_int &y) { return x / gcd(x, y) * y; }
	constexpr big_int divide_with_reminder(const big_int &divisor, big_int *reminder = nullptr) const & { return big_int(*this).divide_with_reminder(divisor, reminder); }
	constexpr big_int divide_with_reminder(const big_int &divisor, big_int *reminder = nullptr) &&;
	// conversion operators
	explicit constexpr BIG_INT_FORCE_INLINE operator bool() const noexcept { return !is_zero(data_); }
	template<std::integral T> explicit constexpr operator T () const noexcept { return T(int64_t(*this)); }
	explicit constexpr operator uint64_t () const noexcept { return data_.front(); }
	explicit constexpr operator int64_t () const noexcept { auto ret = int64_t(uint64_t(*this)); return sign_ ? -ret : ret; }
	template<std::floating_point T>
	explicit constexpr operator T () const noexcept {
		const T exp2_64 = T(std::exp2(64));
		T ret = 0; for (auto x : std::ranges::views::reverse(data_)) ret = ret * exp2_64 + T(x); return sign_ ? -ret : ret;
	}
	explicit constexpr operator std::string() const { return to_string(); }
	template<typename T> constexpr T to() const noexcept(T(big_int{})) { return T(*this); }
	// other helper functions
	constexpr void reserve(uint64_t bits) { data_.reserve(bits_to_digits(bits)); }
	constexpr BIG_INT_FORCE_INLINE const std::vector<uint64_t> &data() const noexcept { return data_; }

	template<typename G>
	static big_int random(uint64_t bits, G &&rng) { return random_impl(bits, std::forward<decltype(rng)>(rng)); }
	static big_int random(uint64_t bits) { thread_local std::mt19937_64 rng(std::random_device{}()); return random(bits, rng); }

	friend auto &operator >> (std::istream &is, big_int &x) { std::string s; is >> s; if (auto v = from_string(s); v) x = std::move(v.value()); return is; }
	friend auto &operator << (std::ostream &os, const big_int &x) { auto [b, o] = get_ios_base_and_options(os); return os << x.to_string(b, o); }
private:
	static constexpr uint64_t B = std::numeric_limits<uint64_t>::max();
	static constexpr BIG_INT_FORCE_INLINE bool is_zero(const std::vector<uint64_t> &data) noexcept { return data.size() == 1 && !data.front(); }
	static constexpr BIG_INT_FORCE_INLINE uint64_t bits_to_digits(uint64_t bits) noexcept { return (bits + 63) / 64; }
	constexpr big_int copy_with_reserve(uint64_t size) const {
		std::vector<uint64_t> data(size);
		data = data_;
		return { std::move(data), sign_ };
	}
	constexpr big_int copy_for_add_minus(const big_int &rhs) const { return copy_with_reserve(std::max(data_.size(), rhs.data_.size()) + 1); }
	static std::pair<base, to_string_option> get_ios_base_and_options(const std::ios_base &ios) noexcept {
		base base = base::decimal;
		auto flags = ios.flags();
		if (flags & std::ios_base::oct) base = base::octal;
		else if (flags & std::ios_base::hex) base = base::hex;
		to_string_option options = to_string_option::none;
		if (flags & std::ios_base::uppercase) options = to_string_option(int(options) | int(to_string_option::uppercase));
		if (flags & std::ios_base::showbase) options = to_string_option(int(options) | int(to_string_option::showbase));
		if (flags & std::ios_base::showpos) options = to_string_option(int(options) | int(to_string_option::showpos));
		return std::make_pair(base, options);
	}
	template<typename Range1, typename Range2>
	static constexpr std::strong_ordering compare_data(const Range1 &lhs, const Range2 &rhs) noexcept {
		if (lhs.size() != rhs.size()) return lhs.size() <=> rhs.size();
		return std::lexicographical_compare_three_way(std::ranges::rbegin(lhs), std::ranges::rend(lhs), std::ranges::rbegin(rhs), std::ranges::rend(rhs));
	}
	template<typename Range>
	static constexpr auto high_n_digits(const Range &range, uint64_t n) { return std::ranges::subrange(range.begin() + range.size() - n, range.end()); }

	constexpr void trim_leading_zeros() noexcept { trim_leading_zeros(&data_); }
	static constexpr void trim_leading_zeros(std::vector<uint64_t> *data) noexcept { while (data->size() > 1 && !data->back()) data->pop_back(); }

	constexpr void add_minus_impl(const big_int &rhs, bool minus) { add_minus_impl(&data_, &sign_, rhs.data_, rhs.sign_ ^ minus); }

	static constexpr void add_minus_impl(std::vector<uint64_t> *lhs, bool *lhs_sign, const std::vector<uint64_t> &rhs, bool rhs_sign) {
		uint64_t n = lhs->size(), m = rhs.size();
		auto l = [&] { return std::span(lhs->data(), n); };
		auto r = [&] { return std::span(rhs.data(), m); };
		if (*lhs_sign == rhs_sign) {
			lhs->resize(std::max(n, m) + 1);
			if (n >= m) add_minus_impl<false>(l(), r(), *lhs);
			else add_minus_impl<false>(r(), l(), *lhs);
		}
		else if (auto order = compare_data(*lhs, rhs); order == std::strong_ordering::equal)
			*lhs = { 0 }, *lhs_sign = false;
		else if (order == std::strong_ordering::greater)
			add_minus_impl<true>(l(), r(), *lhs);
		else
			lhs->resize(m), add_minus_impl<true>(r(), l(), *lhs), *lhs_sign = rhs_sign;
		trim_leading_zeros(lhs);
	}
	template<bool minus>
	static constexpr void add_minus_impl(std::span<const uint64_t> l, std::span<const uint64_t> r, std::span<uint64_t> o) noexcept {
#define BIG_INT_OVERLOADED(func) [](auto &&...args) noexcept { func(std::forward<decltype(args)>(args)...); }
		if constexpr (!minus) add_minus_impl(l, r, o, BIG_INT_OVERLOADED(detail::add));
		else add_minus_impl(l, r, o, BIG_INT_OVERLOADED(detail::minus));
#undef BIG_INT_OVERLOADED
	}
	template<typename Func>
	static constexpr void add_minus_impl(std::span<const uint64_t> l, std::span<const uint64_t> r, std::span<uint64_t> o, Func &&func) noexcept {
		uint64_t carry = 0;
		for (uint64_t i = 0; i < r.size(); i++) {
			uint64_t overflow = 0;
			uint64_t x = l[i];
			func(&x, r[i], &overflow);
			func(&x, carry, &overflow);
			o[i] = x;
			carry = overflow;
		}
		auto j = r.size();
		while (carry && j < o.size()) {
			uint64_t overflow = 0, x = j < l.size() ? l[j] : 0;
			func(&x, carry, &overflow);
			carry = overflow;
			o[j++] = x;
		}
		if (auto min_size = std::min(l.size(), o.size()); j < min_size) std::ranges::copy_n(l.begin() + j, min_size - j, o.begin() + j);
	}

	static constexpr std::vector<uint64_t> multiply_impl(const std::vector<uint64_t> &lhs, const std::vector<uint64_t> &rhs) {
		std::vector<uint64_t> ret;
		multiply_impl(&ret, lhs, rhs);
		return ret;
	}

	static constexpr void multiply_1_impl(std::vector<uint64_t> *output, const std::vector<uint64_t> &lhs, uint64_t rhs) {
		uint64_t n = lhs.size();
		output->resize(n + 1);
		output->front() = 0;
		const uint64_t *l = lhs.data();
		for (uint64_t i = 0; i < n; i++) {
			uint64_t high = 0;
			uint64_t low = detail::umul128(l[i], rhs, &high);
			detail::add(&low, (*output)[i], &high);
			(*output)[i] = low;
			(*output)[i + 1] = high;
		}
		trim_leading_zeros(output);
	}

	static constexpr void multiply_naive_impl(std::vector<uint64_t> *output, const std::vector<uint64_t> &lhs, const std::vector<uint64_t> &rhs) {
		uint64_t n = lhs.size(), m = rhs.size();
		output->resize(n + m);
		std::ranges::fill(*output, 0);
		const uint64_t *l = lhs.data();
		const uint64_t *r = rhs.data();
		for (uint64_t i = 0; i < n; i++) {
			uint64_t carry = 0;
			uint64_t *o = output->data() + i;
			for (uint64_t j = 0; j < m; j++) {
				uint64_t high = 0;
				uint64_t low = detail::umul128(l[i], r[j], &high);
				detail::add(&low, carry, &high);
				carry = high;
				detail::add(&low, *o, &carry);
				*o++ = low;
			}
			*o = carry;
		}
		trim_leading_zeros(output);
	}

	static constexpr void multiply_impl(std::vector<uint64_t> *output, const std::vector<uint64_t> &lhs, const std::vector<uint64_t> &rhs) {
		if (lhs.size() == 1) multiply_1_impl(output, rhs, lhs.front());
		else if (rhs.size() == 1) multiply_1_impl(output, lhs, rhs.front());
		else multiply_naive_impl(output, lhs, rhs);
	}

	static constexpr std::vector<uint64_t> divide_1_impl(std::vector<uint64_t> dividend, uint64_t divisor, uint64_t *remainder) {
		if (dividend.size() == 1 && dividend.front() < divisor) {
			*remainder = dividend.front();
			return { 0 };
		}
		uint64_t n = dividend.size();
		*remainder = 0;
		for (uint64_t i = n - 1; i < n; --i)
			dividend[i] = detail::udiv128(*remainder, dividend[i], divisor, remainder);
		trim_leading_zeros(&dividend);
		return dividend;
	}

	static constexpr std::vector<uint64_t> divide_knuth_impl(std::vector<uint64_t> dividend, std::vector<uint64_t> divisor, std::vector<uint64_t> *output_remainder) {
		if (auto o = compare_data(dividend, divisor); o < 0) {
			if (output_remainder) *output_remainder = std::move(dividend);
			return{ 0 };
		}
		else if (o == 0) {
			if (output_remainder) *output_remainder = { 0 };
			return{ 1 };
		}
		uint64_t m = dividend.size() - divisor.size(), n = divisor.size();
		int shift = std::countl_zero(divisor.back());
		left_shift_impl(&dividend, shift);
		dividend.resize(n + m + 1);
		left_shift_impl(&divisor, shift);
		std::vector<uint64_t> ret(m + 1), tmp;
		tmp.reserve(n + 1);
		for (auto j = m; j <= m; --j) {
			uint64_t r = 0, r_carry = 0;
			uint64_t q = detail::udiv128(dividend[j + n], dividend[j + n - 1], divisor.back(), &r);
			auto test = [&] {
				if (q == B) return true;
				uint64_t high, low = detail::umul128(q, divisor[n - 2], &high);
				return std::tie(high, low) > std::tie(r, dividend[j + n - 2]);
			};
			while (!r_carry && test()) --q, detail::add(&r, divisor.back(), &r_carry);
			multiply_1_impl(&tmp, divisor, q);
			auto before = dividend[j + n];
			auto dividend_span = std::span(dividend.data() + j, n + 1);
			add_minus_impl<true>(dividend_span, tmp, dividend_span);
			ret[j] = q;
			if (before < dividend[j + n]) --ret[j], add_minus_impl<false>(dividend_span, divisor, dividend_span);
		}
		if (output_remainder) {
			*output_remainder = std::move(dividend);
			output_remainder->resize(n);
			right_shift_impl(output_remainder, shift);
		}
		trim_leading_zeros(&ret);
		return ret;
	}

	static constexpr std::vector<uint64_t> divide_impl(std::vector<uint64_t> dividend, const std::vector<uint64_t> &divisor_data, std::vector<uint64_t> *output_remainder) {
		if (divisor_data.size() == 1) {
			uint64_t remainder = 0;
			auto ret = divide_1_impl(std::move(dividend), divisor_data.front(), &remainder);
			if (output_remainder) *output_remainder = { remainder };
			return ret;
		}
		return divide_knuth_impl(std::move(dividend), divisor_data, output_remainder);
	}

	static constexpr void left_shift_impl(std::vector<uint64_t> *output, uint64_t shift) {
		if (is_zero(*output)) return;
		auto m = shift % 64;
		bool shifted_digit = std::countl_zero(output->back()) < m;
		if (auto n = shift / 64) {
			output->resize(output->size() + n + shifted_digit);
			std::ranges::shift_right(*output, n);
			std::ranges::fill_n(output->begin(), n, 0);
		}
		else if (shifted_digit) output->emplace_back();
		if (m) {
			for (uint64_t i = output->size() - 1; i != ~0ull; --i) {
				if (i != output->size() - 1) (*output)[i + 1] |= detail::high_n((*output)[i], m);
				(*output)[i] <<= m;
			}
		}
		trim_leading_zeros(output);
	}
	static constexpr void right_shift_impl(std::vector<uint64_t> *output, uint64_t shift) noexcept {
		if (is_zero(*output)) return;
		if (auto n = shift / 64) {
			if (n >= output->size()) { *output = { 0 }; return; }
			std::ranges::shift_left(*output, n);
			output->resize(output->size() - n);
		}
		if (auto m = shift % 64)
			for (uint64_t i = 0; i < output->size(); ++i) {
				if (i) (*output)[i - 1] |= (*output)[i] << (64 - m);
				(*output)[i] >>= m;
			}
		trim_leading_zeros(output);
	}

	template<typename Op>
	constexpr void bitwise_operation_impl(const big_int &rhs, Op &&op) {
		uint64_t n = rhs.data_.size(), m = data_.size();
		if (n > m) data_.resize(m = n);
		for (uint64_t i = 0; i < n; ++i) op(data_[i], rhs.data_[i]);
		for (uint64_t i = n; i < m; ++i) op(data_[i], uint64_t(0));
	}

	// string conversion literals
	static constexpr uint64_t sc_d = 10000000000000000000ull;
	static constexpr int sc_bunch_digits = 19; // log10(sc_d)
	static constexpr std::optional<big_int> from_string_impl(const char *s, uint64_t length, base *base_in_out) {
		uint32_t b = base_in_out ? uint32_t(*base_in_out) : 0;
		big_int ret;
		bool sign = false;
		if (!length) return ret;
		auto skip = [&](int cnt) { s += cnt; length -= cnt; };
		if (*s == '+') skip(1);
		else if (*s == '-') skip(1), sign = true;
		if (!b && length > 2) {
			auto starts_with = [&](const char *prefix) { return std::strncmp(s, prefix, 2) == 0; };
			if (starts_with("0b") || starts_with("0B")) b = 2, skip(2);
			else if (starts_with("0o") || starts_with("0O")) b = 8, skip(2);
			else if (starts_with("0x") || starts_with("0X")) b = 16, skip(2);
			else b = 10;
		}
		else if (!b) b = 10;
		if (base_in_out) *base_in_out = base(b);
		auto char_to_num = [](char c) -> uint64_t {
			if ('0' <= c && c <= '9') return c - '0';
			if ('a' <= c && c <= 'z') return c - 'a' + 10;
			if ('A' <= c && c <= 'Z') return c - 'A' + 10;
			return B;
		};
		auto str_to_num = [&char_to_num, &b](const char *begin, const char *end) -> uint64_t {
			uint64_t ret = 0;
			while (begin != end) {
				auto num = char_to_num(*(begin++));
				if (num >= b) return B;
				ret *= 10, ret += num;
			}
			return ret;
		};
		if (std::has_single_bit(b)) {
			int bits = std::bit_width(b) - 1;
			ret.data_.resize(bits_to_digits(length * bits));
			uint64_t idx = 0, shift = 0;
			for (const char *cur = s + length - 1; cur >= s; --cur) {
				auto num = char_to_num(*cur);
				if (num >= b) return std::nullopt;
				ret.data_[idx] |= num << shift;
				if ((shift += bits) >= 64) {
					shift -= 64, ++idx;
					if (shift) ret.data_[idx] |= num >> (bits - shift);
				}
			}
			ret.trim_leading_zeros();
		}
		else {
			const auto d(sc_d);
			for (const char *cur = s + length % sc_bunch_digits; cur <= s + length; cur += sc_bunch_digits) {
				auto num = str_to_num(cur - s < sc_bunch_digits ? s : cur - sc_bunch_digits, cur);
				if (num >= sc_d) return std::nullopt;
				ret *= d, ret += num;
			}
		}
		ret.sign_ = sign;
		return ret;
	}

	static constexpr int FastBaseTenDigits(uint64_t x) {
		constexpr auto log2_to_log10 = []() constexpr {
			std::array<int, 64> ret = { 1 };
			for (auto x = 2ull, i = 1ull; i < 64; ++i, x *= 2)
				for (auto t = x - 1; t; t /= 10)
					++ret[i];
			return ret;
		}();
		constexpr auto exp_10 = []() constexpr {
			std::array<uint64_t, 20> ret = { 1 };
			for (int i = 1; i < 20; i++) ret[i] = ret[i - 1] * 10;
			return ret;
		}();
		int log2 = 63 - std::countl_zero(x);
		int log10 = log2_to_log10[log2];
		return log10 + (x >= exp_10[log10]);
	}

	constexpr std::string to_string_impl(uint32_t b, to_string_option option) const {
		auto has_option = [&](to_string_option flag) { return (int(option) & int(flag)) == int(flag); };
		auto num_to_char = [&](auto num) -> char {
			if (0 <= num && num <= 9) return num + '0';
			if (10 <= num && num < b) return num - 10 + (has_option(to_string_option::uppercase) ? 'A' : 'a');
			return 0;
		};
		if (std::has_single_bit(b)) {
			std::string ret = sign_ ? "-" : (has_option(to_string_option::showpos) ? "+" : "");
			if (has_option(to_string_option::showbase)) {
				if (b == 2) ret += "0b";
				else if (b == 8) ret += "0o";
				else ret += "0x";
			}
			int bits = std::bit_width(b) - 1;
			uint64_t idx = data_.size() - 1;
			int shift = (data_.size() * 64 - 1) / bits * bits % 64;
			auto valid = [&]() { return idx < data_.size(); }; // notice that idx is unsigned
			auto next = [&]() { shift -= bits; if (shift < 0) shift += 64, --idx; };
			auto cur = [&]() {
				int p = 64 - shift;
				if (p < bits && idx != data_.size() - 1)
					return detail::high_n(data_[idx], p) | (detail::low_n(data_[idx + 1], bits - p) << p);
				return detail::low_n(detail::high_n(data_[idx], p), bits);
			};
			while (valid() && !cur()) next();
			ret.reserve(ret.size() + shift / bits + 64 / bits * idx + 1);
			if (!valid()) return ret += '0';
			while (valid()) { ret += num_to_char(cur()); next(); }
			return ret;
		}
		std::string ret;
		const big_int d(sc_d);
		big_int r = *this;
		r.sign_ = false;
		while (r) {
			big_int cur;
			r.data_ = divide_impl(std::move(r.data_), d.data_, &cur.data_);
			int digits = FastBaseTenDigits(int64_t(cur));
			std::string num = std::string(sc_bunch_digits - digits, '0') + std::to_string(uint64_t(cur));
			std::ranges::reverse(num);
			ret += num;
		}
		while (!ret.empty() && ret.back() == '0') ret.pop_back();
		if (ret.empty()) return has_option(to_string_option::showpos) ? "+0" : "0";
		if (sign_) ret += '-';
		else if (has_option(to_string_option::showpos)) ret += '+';
		std::ranges::reverse(ret);
		return ret;
	}
	template<typename G>
	static big_int random_impl(uint64_t bits, G &&rng) {
		using dist_t = std::uniform_int_distribution<uint64_t>;
		dist_t dist;
		big_int ret;
		ret.data_.resize(bits_to_digits(bits));
		for (uint64_t i = 0; i + 1 < ret.data_.size(); ++i)
			ret.data_[i] = dist(rng);
		ret.data_.back() = dist_t(0ull, (1ull << (bits % 64)) - 1)(rng);
		ret.trim_leading_zeros();
		return ret;
	}

	bool sign_ = false;
	std::vector<uint64_t> data_;
};

constexpr std::optional<big_int> big_int::from_string(std::string_view sv, base* b_output) {
	base b = base::auto_detect;
	auto ret = from_string_impl(sv.data(), sv.size(), &b);
	if (b_output) *b_output = b;
	return ret;
}

constexpr big_int &big_int::operator *= (const big_int &rhs) {
	if (data_.capacity() >= data_.size() + rhs.data_.size() && *this && rhs) {
		auto lhs = data_;
		if (std::addressof(rhs) == this) multiply_impl(&data_, lhs, lhs);
		else multiply_impl(&data_, lhs, rhs.data_);
		sign_ ^= rhs.sign_;
		return *this;
	}
	return *this = *this * rhs;
}

constexpr uint64_t big_int::count_trailing_zeros() const {
	uint64_t ret = 0ull;
	for (auto d : data_) {
		ret += std::countr_zero(d);
		if (d) break;
	}
	return ret;
}

constexpr big_int big_int::power(uint64_t n) const {
	if (n == 0 || *this == 1) return 1; if (n == 1 || *this == 0) return *this;
	big_int y(1);
	y.reserve(highest_bit() * n);
	std::optional<big_int> cur;
	while (true) {
		if (n % 2) y *= cur ? *cur : *this;
		if (n /= 2)
			if (cur) *cur *= *cur;
			else cur = *this * *this;
		else break;
	}
	return y;
}

constexpr big_int big_int::factorial(uint64_t n) {
	big_int ret_int = 1;
	if (!std::is_constant_evaluated()) {
		auto esti = std::llround(double(n + 1) * (std::log2(n) - std::numbers::log2e));
		if (esti < 63) esti = 63;
		ret_int.reserve(esti + 64);
	}
	for (auto i = n; i > 1; --i) ret_int *= i;
	return ret_int;
}

constexpr big_int big_int::gcd(big_int x, big_int y) {
	if (x.sign_ || y.sign_) { return gcd(std::move(x).absolute(), std::move(y).absolute()); }
	uint64_t pow = 0;
	while (true) {
		if (!x) return std::move(y <<= pow); if (!y) return std::move(x <<= pow); if (x == 1 || y == 1) return big_int(1) << pow;
		if (auto xzs = x.count_trailing_zeros(), yzs = y.count_trailing_zeros(); xzs || yzs) { pow += std::min(xzs, yzs); x >>= xzs; y >>= yzs; continue; }
		if (x < y) std::swap(x, y);
		x -= y;
	}
}

constexpr big_int big_int::divide_with_reminder(const big_int &divisor, big_int *reminder) && {
	if (std::addressof(divisor) == this) [[unlikely]] { if (reminder) *reminder = 0; return 1; }
	big_int ret = { divide_impl(std::move(data_), divisor.data_, reminder ? &reminder->data_ : nullptr), sign_ != divisor.sign_ };
	if (reminder) { reminder->sign_ = false; if (*reminder) reminder->sign_ = sign_; }
	return ret;
}

namespace big_int_literals {
	constexpr big_int operator ""_bi (unsigned long long int x) { return big_int(x); }
	constexpr big_int operator ""_bi (const char *s, std::size_t len) { return *big_int::from_string({ s, len }); }
}

namespace detail {

static constexpr BIG_INT_FORCE_INLINE uint64_t umul128_naive_impl(uint64_t lhs, uint64_t rhs, uint64_t *high_product) {
	uint64_t h1 = high_32(lhs), l1 = low_32(lhs);
	uint64_t h2 = high_32(rhs), l2 = low_32(rhs);
	uint64_t ret = l1 * l2;
	*high_product = h1 * h2;
	// deal with (l1 * h2 + l2 * h1) << 32
	uint64_t m_carry = 0;
	uint64_t m = l1 * h2;
	add(&m, l2 * h1, &m_carry);
	uint64_t hm = high_32(m), lm = low_32(m);
	uint64_t carry = 0;
	add(&ret, lm << 32, &carry);
	*high_product += hm;
	*high_product += carry;
	*high_product += m_carry << 32;
	return ret;
}

static constexpr BIG_INT_FORCE_INLINE uint64_t udiv128_naive_impl(uint64_t high_dividend, uint64_t low_dividend, uint64_t divisor, uint64_t *remainder) {
	if (high_dividend >= divisor) { return *remainder = std::numeric_limits<uint64_t>::max(); }
	int shift = std::countl_zero(divisor);
	divisor <<= shift;
	high_dividend <<= shift;
	if (shift) high_dividend |= high_n(low_dividend, shift);
	low_dividend <<= shift;
	auto low_dividend_high = high_32(low_dividend), low_dividend_low = low_32(low_dividend);
	auto divisor_high = high_32(divisor), divisor_low = low_32(divisor);
	auto get_q = [&](auto h, auto num) {
		auto q = h / divisor_high, r = h % divisor_high;
		if (auto c1 = q * divisor_low, c2 = r << 32 | num; c1 > c2) q -= (c1 - c2 > divisor) ? 2 : 1;
		return q;
	};
	auto get_ram = [&](auto h, auto num, auto q) { return (h << 32 | num) - q * divisor; };
	auto q1 = get_q(high_dividend, low_dividend_high);
	auto rem = get_ram(high_dividend, low_dividend_high, q1);
	auto q0 = get_q(rem, low_dividend_low);
	*remainder = get_ram(rem, low_dividend_low, q0) >> shift;
	return q1 << 32 | q0;
}

static constexpr BIG_INT_FORCE_INLINE uint64_t umul128(uint64_t lhs, uint64_t rhs, uint64_t *high_product) {
	if (std::is_constant_evaluated()) return umul128_naive_impl(lhs, rhs, high_product);
#ifdef BIG_INT_HAS_UMUL128_UDIV128
	return _umul128(lhs, rhs, high_product);
#elif defined(BIG_INT_HAVE_INTRINSIC_INT128)
	using T = unsigned __int128;
	T res = T(lhs) * T(rhs);
	*high_product = uint64_t(ret >> 64);
	return uint64_t(res & ~uint64_t(0));
#else
	return naive_impl();
#endif
}

static constexpr BIG_INT_FORCE_INLINE uint64_t udiv128(uint64_t high_dividend, uint64_t low_dividend, uint64_t divisor, uint64_t *remainder) {
	if (std::is_constant_evaluated()) return udiv128_naive_impl(high_dividend, low_dividend, divisor, remainder);
#ifdef BIG_INT_HAS_UMUL128_UDIV128
	return _udiv128(high_dividend, low_dividend, divisor, remainder);
#elif defined(BIG_INT_HAVE_INTRINSIC_INT128)
	using T = unsigned __int128;
	T dividend = (T(high_dividend) << 64) | low_dividend;
	*remainder = dividend % divisor;
	return uint64_t(dividend / divisor);
#else
	return naive_impl();
#endif
}

}

}

#ifdef BIG_INT_HAVE_INTRINSIC_INT128
#undef BIG_INT_HAVE_INTRINSIC_INT128
#endif

#ifdef BIG_INT_HAS_UMUL128_UDIV128
#undef BIG_INT_HAS_UMUL128_UDIV128
#endif
