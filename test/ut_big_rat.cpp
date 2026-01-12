#include "big_rat.h"
#include "external/ut.hpp"

using namespace spaceless;
using namespace big_int_literals;
using namespace big_rat_literals;

static void unit_test() {
	using namespace boost::ut;
	using namespace boost::ut::bdd;

	"big_rat"_test = [] {
		feature("construction") = [] {
			expect(eq(big_rat(3, 4).numerator(), 3_bi));
			expect(eq(big_rat(3, 4).denominator(), 4_bi));
			expect(eq(big_rat(6, 8).numerator(), 3_bi));
			expect(eq(big_rat(6, 8).denominator(), 4_bi));
			expect(eq(big_rat(-6, 8).numerator(), -3_bi));
			expect(eq(big_rat(6, -8).numerator(), -3_bi));
			expect(eq(big_rat(42).numerator(), 42_bi));
			expect(eq(big_rat(42).denominator(), 1_bi));
			expect(eq(big_rat(0.5).numerator(), 1_bi));
			expect(eq(big_rat(0.5).denominator(), 2_bi));
		};
		feature("string") = [] {
			expect(eq(big_rat::from_string("3/4")->numerator(), 3_bi));
			expect(eq(big_rat::from_string("3/4")->denominator(), 4_bi));
			expect(eq(big_rat(3, 4).to_string(), std::string("3/4")));
			auto seventh = big_rat::from_repeating_decimal(big_rat::repeating_decimal::from_string("0.(142857)"));
			expect(eq(seventh->numerator(), 1_bi));
			expect(eq(seventh->denominator(), 7_bi));
		};
		feature("arithmetic") = [] {
			auto a = big_rat(1, 2), b = big_rat(1, 3);
			expect(eq(a + b, big_rat(5, 6)));
			expect(eq(a - b, big_rat(1, 6)));
			expect(eq(a * b, big_rat(1, 6)));
			expect(eq(a / b, big_rat(3, 2)));
			expect(eq(a + 1, big_rat(3, 2)));
			expect(eq(a * 2, big_rat(1, 1)));
			expect(eq(1 + a, big_rat(3, 2)));
			expect(eq(2 * a, big_rat(1, 1)));
			auto c = a; c += b; expect(eq(c, big_rat(5, 6)));
			c = a; c -= b; expect(eq(c, big_rat(1, 6)));
			c = a; c *= b; expect(eq(c, big_rat(1, 6)));
			c = a; c /= b; expect(eq(c, big_rat(3, 2)));
			c = a; ++c; expect(eq(c, big_rat(3, 2)));
			c = a; --c; expect(eq(c, big_rat(-1, 2)));
		};
		feature("comparison") = [] {
			expect(big_rat(1, 2) == big_rat(2, 4));
			expect(big_rat(1, 2) < big_rat(2, 3));
			expect(big_rat(1, 2) > big_rat(1, 3));
			expect(big_rat(-1, 2) < big_rat(1, 2));
		};
		feature("sign") = [] {
			auto a = big_rat(3, 4);
			expect(!a.sign());
			expect((-a).sign());
			expect(eq(a.negation(), big_rat(-3, 4)));
			expect(eq((-a).absolute(), a));
			auto b = a; b.negate(); expect(eq(b, big_rat(-3, 4)));
			b = -a; b.absolutize(); expect(eq(b, a));
		};
		feature("invert") = [] {
			auto a = big_rat(3, 4);
			expect(eq(a.inverse(), big_rat(4, 3)));
			auto b = a; b.invert(); expect(eq(b, big_rat(4, 3)));
		};
		feature("floor_ceil_round_trunc") = [] {
			expect(eq(big_rat(5, 2).floor(), 2_bi));
			expect(eq(big_rat(5, 2).ceil(), 3_bi));
			expect(eq(big_rat(5, 2).round(), 3_bi));
			expect(eq(big_rat(5, 2).trunc(), 2_bi));
			expect(eq(big_rat(-5, 2).floor(), -3_bi));
			expect(eq(big_rat(-5, 2).ceil(), -2_bi));
			expect(eq(big_rat(-5, 2).round(), -3_bi));
			expect(eq(big_rat(-5, 2).trunc(), -2_bi));
			expect(eq(big_rat(4, 2).floor(), 2_bi));
			expect(eq(big_rat(4, 2).ceil(), 2_bi));
			expect(eq(big_rat(4, 2).round(), 2_bi));
			expect(eq(big_rat(1, 3).floor(), 0_bi));
			expect(eq(big_rat(1, 3).ceil(), 1_bi));
			expect(eq(big_rat(1, 3).round(), 0_bi));
		};
		feature("power") = [] {
			expect(eq(big_rat(2, 3).power(3), big_rat(8, 27)));
			expect(eq(big_rat(2, 3).power(-2), big_rat(9, 4)));
			expect(eq(big_rat(2, 3).power(0), big_rat(1)));
		};
		feature("limit_denominator") = [] {
			auto pi = big_rat::from_string("314159265358979/100000000000000").value();
			expect(eq(pi.limit_denominator(100), big_rat(311, 99)));
			expect(eq(pi.limit_denominator(1000), big_rat(355, 113)));
		};
		feature("conversion") = [] {
			expect(bool(big_rat(1, 2)));
			expect(!bool(big_rat(0)));
			expect(eq(big_rat(3, 2).to<double>(), 1.5));
			expect(eq(big_int(big_rat(7, 3)), 2_bi));
		};
	};
}

int main(int argc, char **argv) {
	unit_test();
	return 0;
}
