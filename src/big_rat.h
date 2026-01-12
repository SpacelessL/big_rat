#pragma once

#include "big_int.h"

namespace spaceless {

class big_rat final {
public:
	using base = big_int::base;
	using to_string_option = big_int::to_string_option;

	static constexpr std::optional<big_rat> from_string(const std::string_view &sv, base b = base::auto_detect) { return from_string_impl(sv, b); }
	constexpr std::string to_string(base b = base::decimal, to_string_option option = to_string_option::showbase) const { return numerator_.to_string(b, option) + '/' + denominator_.to_string(b, option); }

	struct repeating_decimal {
		base b = base::auto_detect;
		std::string integer_part, non_repeating_fraction_part, repeating_fraction_part;

		constexpr std::string to_string() const;
		static constexpr repeating_decimal from_string(const std::string_view &sv, base b = base::auto_detect);
	};
	static constexpr std::optional<big_rat> from_repeating_decimal(const repeating_decimal &rd) { return from_repeating_decimal_impl(rd); }
	// constructors
	constexpr big_rat() : big_rat(0) {}
	constexpr big_rat(const big_rat &) = default;
	constexpr big_rat(big_rat &&) noexcept = default;
	explicit constexpr big_rat(big_int num) : big_rat(std::move(num), 1, {}) {}
	constexpr big_rat(big_int num, big_int denom) : big_rat(std::move(num), std::move(denom), {}) { simplify(); }
	template<std::integral T>
	constexpr big_rat(T num) : big_rat(big_int(num)) {}
	template<std::floating_point T>
	requires(std::numeric_limits<T>::is_iec559 && (std::is_same_v<T, float> || std::is_same_v<T, double>))
	constexpr big_rat(T val) : big_rat(std::isnan(val) || std::isinf(val) ? big_rat(0) : from_floating_point(val)) {}
	constexpr ~big_rat() noexcept = default;
	// assignments
	constexpr big_rat &operator = (const big_rat &) = default;
	constexpr big_rat &operator = (big_rat &&) noexcept = default;
	// unary operators, also see sign-related functions
	constexpr const big_rat &operator + () const & { return *this; }
	constexpr big_rat operator - () const & { return negation(); }
	constexpr big_rat operator + () && noexcept { return std::move(*this); }
	constexpr big_rat operator - () && noexcept { return std::move(*this).negation(); }
	// binary operators
	friend constexpr big_rat operator + (const big_rat &lhs, const big_rat &rhs) { return { lhs.numerator_ * rhs.denominator_ + rhs.numerator_ * lhs.denominator_, lhs.denominator_ * rhs.denominator_ }; }
	friend constexpr big_rat operator - (const big_rat &lhs, const big_rat &rhs) { return { lhs.numerator_ * rhs.denominator_ - rhs.numerator_ * lhs.denominator_, lhs.denominator_ * rhs.denominator_ }; }
	friend constexpr big_rat operator * (const big_rat &lhs, const big_rat &rhs) { return { lhs.numerator_ * rhs.numerator_, lhs.denominator_ * rhs.denominator_ }; }
	friend constexpr big_rat operator / (const big_rat &lhs, const big_rat &rhs) { return { lhs.numerator_ * rhs.denominator_, lhs.denominator_ * rhs.numerator_ }; }
	friend constexpr big_rat operator + (const big_rat &lhs, const big_int &rhs) { return { lhs.numerator_ + rhs * lhs.denominator_, lhs.denominator_ }; }
	friend constexpr big_rat operator - (const big_rat &lhs, const big_int &rhs) { return { lhs.numerator_ - rhs * lhs.denominator_, lhs.denominator_ }; }
	friend constexpr big_rat operator * (const big_rat &lhs, const big_int &rhs) { return { lhs.numerator_ * rhs, lhs.denominator_ }; }
	friend constexpr big_rat operator / (const big_rat &lhs, const big_int &rhs) { return { lhs.numerator_, lhs.denominator_ * rhs }; }
	friend constexpr big_rat operator + (big_rat &&lhs, const big_int &rhs) { auto num = lhs.numerator_ + rhs * lhs.denominator_; return { std::move(num), std::move(lhs.denominator_) }; }
	friend constexpr big_rat operator - (big_rat &&lhs, const big_int &rhs) { auto num = lhs.numerator_ - rhs * lhs.denominator_; return { std::move(num), std::move(lhs.denominator_) }; }
	friend constexpr big_rat operator * (big_rat &&lhs, const big_int &rhs) { auto num = lhs.numerator_ * rhs; return { std::move(num), std::move(lhs.denominator_) }; }
	friend constexpr big_rat operator / (big_rat &&lhs, const big_int &rhs) { auto denom = lhs.denominator_ * rhs; return { std::move(lhs.numerator_), std::move(denom) }; }
	template<std::integral T> friend constexpr big_rat operator + (const big_rat &lhs, T rhs) { return lhs + big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator - (const big_rat &lhs, T rhs) { return lhs - big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator * (const big_rat &lhs, T rhs) { return lhs * big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator / (const big_rat &lhs, T rhs) { return lhs / big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator + (big_rat &&lhs, T rhs) { return std::move(lhs) + big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator - (big_rat &&lhs, T rhs) { return std::move(lhs) - big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator * (big_rat &&lhs, T rhs) { return std::move(lhs) * big_int(rhs); }
	template<std::integral T> friend constexpr big_rat operator / (big_rat &&lhs, T rhs) { return std::move(lhs) / big_int(rhs); }
	friend constexpr big_rat operator + (const big_int &lhs, const big_rat &rhs) { return rhs + lhs; }
	friend constexpr big_rat operator - (const big_int &lhs, const big_rat &rhs) { return -(rhs - lhs); }
	friend constexpr big_rat operator * (const big_int &lhs, const big_rat &rhs) { return rhs * lhs; }
	friend constexpr big_rat operator / (const big_int &lhs, const big_rat &rhs) { return (rhs / lhs).inverse(); }
	friend constexpr big_rat operator + (const big_int &lhs, big_rat &&rhs) { return std::move(rhs) + lhs; }
	friend constexpr big_rat operator - (const big_int &lhs, big_rat &&rhs) { return -(std::move(rhs) - lhs); }
	friend constexpr big_rat operator * (const big_int &lhs, big_rat &&rhs) { return std::move(rhs) * lhs; }
	friend constexpr big_rat operator / (const big_int &lhs, big_rat &&rhs) { return (std::move(rhs) / lhs).inverse(); }
	template<std::integral T> friend constexpr big_rat operator + (T lhs, const big_rat &rhs) { return rhs + lhs; }
	template<std::integral T> friend constexpr big_rat operator - (T lhs, const big_rat &rhs) { return -(rhs - lhs); }
	template<std::integral T> friend constexpr big_rat operator * (T lhs, const big_rat &rhs) { return rhs * lhs; }
	template<std::integral T> friend constexpr big_rat operator / (T lhs, const big_rat &rhs) { return (rhs / lhs).inverse(); }
	template<std::integral T> friend constexpr big_rat operator + (T lhs, big_rat &&rhs) { return std::move(rhs) + lhs; }
	template<std::integral T> friend constexpr big_rat operator - (T lhs, big_rat &&rhs) { return -(std::move(rhs) - lhs); }
	template<std::integral T> friend constexpr big_rat operator * (T lhs, big_rat &&rhs) { return std::move(rhs) * lhs; }
	template<std::integral T> friend constexpr big_rat operator / (T lhs, big_rat &&rhs) { return (std::move(rhs) / lhs).inverse(); }
	// assignment operators
	constexpr big_rat &operator += (const big_rat &rhs) { return *this = *this + rhs; }
	constexpr big_rat &operator -= (const big_rat &rhs) { return *this = *this - rhs; }
	constexpr big_rat &operator *= (const big_rat &rhs) { return *this = *this * rhs; }
	constexpr big_rat &operator /= (const big_rat &rhs) { return *this = *this / rhs; }
	constexpr big_rat &operator += (const big_int &rhs) { return *this = std::move(*this) + rhs; }
	constexpr big_rat &operator -= (const big_int &rhs) { return *this = std::move(*this) - rhs; }
	constexpr big_rat &operator *= (const big_int &rhs) { return *this = std::move(*this) * rhs; }
	constexpr big_rat &operator /= (const big_int &rhs) { return *this = std::move(*this) / rhs; }
	template<std::integral T> constexpr big_rat &operator += (T rhs) { return *this += big_int(rhs); }
	template<std::integral T> constexpr big_rat &operator -= (T rhs) { return *this -= big_int(rhs); }
	template<std::integral T> constexpr big_rat &operator *= (T rhs) { return *this *= big_int(rhs); }
	template<std::integral T> constexpr big_rat &operator /= (T rhs) { return *this /= big_int(rhs); }
	// increment and decrement operators
	constexpr big_rat &operator ++ () { numerator_ += denominator_; simplify(); return *this; }
	constexpr big_rat &operator -- () { numerator_ -= denominator_; simplify(); return *this; }
	[[nodiscard]] constexpr big_rat operator ++ (int) { auto ret = *this; ++*this; return ret; }
	[[nodiscard]] constexpr big_rat operator -- (int) { auto ret = *this; --*this; return ret; }
	// relational operators
	constexpr bool operator == (const big_rat &) const noexcept = default;
	constexpr std::strong_ordering operator <=> (const big_rat &rhs) const noexcept { return compare_impl(*this, rhs); }
	// sign-related functions
	constexpr bool sign() const noexcept { return numerator_.sign(); }
	constexpr void set_sign(bool s) noexcept { numerator_.set_sign(s); }
	constexpr void absolutize() noexcept { numerator_.absolutize(); }
	[[nodiscard]] constexpr big_rat absolute() const & { return big_rat(*this).absolute(); }
	[[nodiscard]] constexpr big_rat absolute() && noexcept { absolutize(); return std::move(*this); }
	constexpr void negate() noexcept { numerator_.negate(); }
	[[nodiscard]] constexpr big_rat negation() const & { return big_rat(*this).negation(); }
	[[nodiscard]] constexpr big_rat negation() && noexcept { negate(); return std::move(*this); }
	constexpr void invert() noexcept { std::swap(numerator_, denominator_); }
	[[nodiscard]] constexpr big_rat inverse() const & { return big_rat(*this).inverse(); }
	[[nodiscard]] constexpr big_rat inverse() && noexcept { invert(); return std::move(*this); }
	// arithmetic functions
	[[nodiscard]] constexpr big_int floor() const;
	[[nodiscard]] constexpr big_int ceil() const;
	[[nodiscard]] constexpr big_int round() const;
	[[nodiscard]] constexpr big_int trunc() const { return numerator_ / denominator_; }
	[[nodiscard]] constexpr big_rat power(int64_t n) const;
	[[nodiscard]] constexpr big_rat limit_denominator(const big_int &max_denominator = 1000000) const;
	// conversion operators
	explicit constexpr operator bool() const noexcept { return bool(numerator_); }
	template<std::floating_point T>
	requires(std::numeric_limits<T>::is_iec559)
	explicit constexpr operator T() const { return to_floating_point<T>(); }
	explicit constexpr operator big_int() const { return numerator_ / denominator_; }
	template<typename T>
	constexpr T to() const { return T(*this); }
	// other helper functions
	constexpr const big_int &numerator() const noexcept { return numerator_; }
	constexpr const big_int &denominator() const noexcept { return denominator_; }

	friend std::istream &operator >> (std::istream &is, big_rat &x) { std::string s; is >> s; auto v = from_string(s); if (v) x = std::move(v.value()); return is; }
	friend std::ostream &operator << (std::ostream &os, const big_rat &x) { os << x.numerator() << '/' << x.denominator(); return os; }
private:
	struct DoNotSimplify {};
	constexpr big_rat(big_int num, big_int denom, DoNotSimplify) : numerator_(std::move(num)), denominator_(std::move(denom)) {}

	constexpr void simplify() {
		if (denominator_.sign()) numerator_.negate(), denominator_.negate();
		auto gcd = big_int::gcd(numerator_, denominator_); numerator_ /= gcd; denominator_ /= gcd;
	}

	static constexpr std::optional<big_rat> from_string_impl(const std::string_view &sv, base b) {
		auto get = [b](std::string_view sv) { return from_repeating_decimal(repeating_decimal::from_string(sv, b)); };
		if (auto p = sv.find('/'); p != std::string::npos) {
			if (p == 0 || p == sv.length() - 1) return std::nullopt;
			auto num = get(sv.substr(0, p)), denom = get(sv.substr(p + 1));
			if (num && denom) return num.value() / denom.value();
			return std::nullopt;
		}
		return get(sv);
	}

	static constexpr std::optional<big_rat> from_repeating_decimal_impl(const repeating_decimal &rd) {
		base b = rd.b;
		auto get_big_int = [&](const std::string &string) -> std::optional<big_int> {
			if (string.empty()) return big_int(0);
			return b == base::auto_detect ? big_int::from_string(string, &b) : big_int::from_string(string, b);
		};
		auto integer = get_big_int(rd.integer_part); if (!integer) return std::nullopt;
		auto nr = get_big_int(rd.non_repeating_fraction_part); if (!nr) return std::nullopt;
		auto r = get_big_int(rd.repeating_fraction_part); if (!r) return std::nullopt;
		if (b == base::auto_detect) b = base::decimal;
		auto bb = big_int(uint32_t(b));
		auto nr_denom = bb.power(rd.non_repeating_fraction_part.length());
		nr->set_sign(integer->sign());
		return big_rat(*integer) + big_rat(*nr, nr_denom)
			+ (*r ? big_rat(*r, nr_denom * (bb.power(rd.repeating_fraction_part.length()) - 1)) : big_rat(0));
	}

	static constexpr std::strong_ordering compare_impl(const big_rat &lhs, const big_rat &rhs) noexcept {
		if (lhs.sign() != rhs.sign()) return lhs.sign() ? std::strong_ordering::less : std::strong_ordering::greater;
		auto lmax = lhs.numerator_.highest_bit() + rhs.denominator_.highest_bit(), rmax = lhs.denominator_.highest_bit() + rhs.numerator_.highest_bit();
		if (lmax - 1 > rmax) return lhs.sign() ? std::strong_ordering::less : std::strong_ordering::greater;
		if (rmax - 1 > lmax) return lhs.sign() ? std::strong_ordering::greater : std::strong_ordering::less;
		return lhs.numerator_ * rhs.denominator_ <=> lhs.denominator_ * rhs.numerator_;
	}

	static big_rat from_floating_point(float x) { return from_floating_point<8, 23, uint32_t>(&x); }
	static big_rat from_floating_point(double x) { return from_floating_point<11, 52, uint64_t>(&x); }
	template<unsigned int exponent, unsigned int significand, typename T>
	static big_rat from_floating_point(void *data) {
		static_assert((exponent + significand + 1) % 32 == 0);
		T x = *static_cast<const T *>(data);
		bool sign = x >> (exponent + significand);
		int64_t exp = (x << 1) >> (significand + 1);
		exp -= (1 << (exponent - 1)) - 1;
		bool denormalized = exp == -127;
		if (denormalized) exp = -126;
		T frac = x & ((1ull << significand) - 1);
		if (!denormalized) frac |= 1ull << significand;
		auto ret = big_rat(big_int(frac), big_int(1) << (significand - exp));
		return sign ? -std::move(ret) : ret;
	}
	template<typename T>
	T to_floating_point() const { return T(numerator_) / T(denominator_); }

	big_int numerator_, denominator_;
};

namespace big_rat_literals {
	constexpr big_rat operator ""_br (unsigned long long int x) { return { x }; }
	constexpr big_rat operator ""_br (long double x) { return { double(x) }; }
	constexpr big_rat operator ""_br (const char *s, std::size_t len) { return big_rat::from_string({ s, len }).value(); }
}

constexpr std::string big_rat::repeating_decimal::to_string() const {
	std::string ret = integer_part;
	if (ret.empty()) ret = "0";
	if (non_repeating_fraction_part.empty() && repeating_fraction_part.empty()) return ret;
	ret += '.'; ret += non_repeating_fraction_part;
	if (!repeating_fraction_part.empty()) ret += '(' + repeating_fraction_part + ')';
	return ret;
}

constexpr big_rat::repeating_decimal big_rat::repeating_decimal::from_string(const std::string_view &sv, base b) {
	repeating_decimal ret;
	ret.b = b;
	if (auto dot_p = sv.find('.'); dot_p != std::string::npos) {
		ret.integer_part = sv.substr(0, dot_p);
		if (dot_p == sv.size() - 1) return ret;
		if (auto p = sv.find('(', dot_p); p != std::string::npos) {
			ret.non_repeating_fraction_part = sv.substr(dot_p + 1, p - dot_p - 1);
			ret.repeating_fraction_part = sv.substr(p + 1, sv.size() - p - 2);
		}
		else
			ret.non_repeating_fraction_part = sv.substr(dot_p + 1);
	}
	else ret.integer_part = sv;
	return ret;
}

constexpr big_rat big_rat::power(int64_t n) const {
	if (!n) return 1;
	big_rat ret(numerator_.power(std::abs(n)), denominator_.power(std::abs(n)), {});
	if (n < 0) ret.invert();
	return ret;
}

constexpr big_rat big_rat::limit_denominator(const big_int &max_denominator) const {
	if (denominator_ <= max_denominator)
		return *this;
	big_int p0 = 0, q0 = 1, p1 = 1, q1 = 0, n = numerator_.absolute(), d = denominator_;
	while (true) {
		big_int b, a = n.divide_with_reminder(d, &b);
		auto q2 = q0 + a * q1;
		if (q2 > max_denominator)
			break;
		p0 += a * p1;
		std::swap(p0, p1);
		q0 = std::move(q1);
		q1 = std::move(q2);
		n = std::move(d);
		d = std::move(b);
	}

	auto k = (max_denominator - q0) / q1;
	auto bound1 = big_rat(std::move(p0) + k * p1, std::move(q0) + k * q1);
	auto bound2 = big_rat(std::move(p1), std::move(q1));
	bound1.set_sign(sign());
	bound2.set_sign(sign());
	return (bound2 - *this).absolute() <= (bound1 - *this).absolute() ? bound2 : bound1;
}

constexpr big_int big_rat::floor() const {
	big_int rem, quot = numerator_.divide_with_reminder(denominator_, &rem);
	if (sign() && rem) --quot;
	return quot;
}

constexpr big_int big_rat::ceil() const {
	big_int rem, quot = numerator_.divide_with_reminder(denominator_, &rem);
	if (!sign() && rem) ++quot;
	return quot;
}

constexpr big_int big_rat::round() const {
	big_int rem, quot = numerator_.divide_with_reminder(denominator_, &rem);
	rem.absolutize();
	if ((rem << 1) >= denominator_) { if (sign()) --quot; else ++quot; }
	return quot;
}

}
