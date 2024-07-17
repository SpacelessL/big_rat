#include "big_int.h"
#include "external/ut.hpp"

#include <fstream>

using namespace spaceless;
using namespace big_int_literals;

template<typename T>
auto &operator << (std::ostream &os, const std::optional<T> &x) { if (!x) return os << "[nullopt]"; return os << *x; }

static void unit_test() {
	using namespace boost::ut;
	using namespace boost::ut::bdd;

	"big_int"_test = [] {
		feature("basic") = [] {
			{
				auto a = "0x167fffffffffffffffffffffffffffffa60000000000000000000000000000001549400000000000000000000000000229b5fffffffffffffffffffffffffff97fbf4000000000000000000000000008b352fffffffffffffffffffffffffffba01"_bi;
				auto b = "0x7e8ffffffffffffffffffffffffffffc0b80000000000000000000000000000be31ffffffffffffffffffffffffffff017800000000000000000000000000007ff9"_bi;
				auto c = "0x167fffffffffffffffffffffffffffffa5ffffffffffffffffffffffffffffff96b940000000000000000000000000061e35ffffffffffffffffffffffffffed9c9f40000000000000000000000000189bd2fffffffffffffffffffffffffff3a08"_bi;
				expect(eq(a - b, c));
			}
			auto a = "18658316586846546587789749"_bi;
			auto b = a.power(13);
			auto c = "33215588430304983195541539375682464839985654461394098514621059173680386919467827771862952765123987087013273061130144991703574128633748404479096130602468054549085269859006532267535358492450584752280012118799821038295237387310743703453485208812716091333318782521196132832473357275721628178731752174526818361800007516545562521391749"_bi;
			expect(eq(b, c));
			auto d = "96416446134622"_bi;
			expect(lt(d, a));
			expect(eq(b / a, a.power(12)));
			expect(eq(b % a, 0_i));
			expect(eq((b + d) / a, a.power(12)));
			expect(eq((b + d) % a, d));
			expect(eq((2_bi).power(1677215), (1_bi << 1677215)));

			expect(eq(a + a, a * 2));
			expect(eq(b + b, b * 2));
			expect(eq(c + c, c * 2));
			expect(eq(d + d, d * 2));

			expect(eq(a + a, a << 1));
			expect(eq(b + b, b << 1));
			expect(eq(c + c, c << 1));
			expect(eq(d + d, d << 1));

			expect(eq(a - a, 0));
			expect(eq(b - b, 0));
			expect(eq(c - c, 0));
			expect(eq(d - d, 0));

			expect(eq(a / a, 1));
			expect(eq(b / b, 1));
			expect(eq(c / c, 1));
			expect(eq(d / d, 1));

			expect(eq(a % a, 0));
			expect(eq(b % b, 0));
			expect(eq(c % c, 0));
			expect(eq(d % d, 0));

			auto test_self_op = [](big_int x) {
				big_int y = x;
				expect(eq(y + y, x + x));
				expect(eq(y + y, 2 * x));
				expect(eq(y + y, x << 1));
				expect(eq(y - y, 0));
				expect(eq(y * y, x * x));
				expect(eq(y * y, x.power(2)));
				expect(eq(y / y, 1));
				expect(eq(y % y, 0));
				y = x; y += y; expect(eq(y, x + x));
				y = x; y -= y; expect(eq(y, 0));
				y = x; y *= y; expect(eq(y, x * x));
				y = x; y /= y; expect(eq(y, 1));
				y = x; y %= y; expect(eq(y, 0));
			};
			test_self_op(a);
			test_self_op(b);
			test_self_op(c);
			test_self_op(d);
		};
		feature("gcd") = [] {
			auto a = "0x2cffffffffffffffffffffffffffffff4c000000000000000000000000000000b5"_bi;
			auto b = "0x1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"_bi;
			auto c = "0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"_bi;
			auto d = c - b;
			auto e = "0x7ffffffffffffffffffffe0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"_bi;
			expect(eq(d, e));
			expect(eq(big_int::gcd(a * b, a * c), a));
			expect(eq(big_int::gcd(b * a, b * c), b));
			expect(eq(big_int::gcd(c * a, c * b), c));
			expect(eq(big_int::gcd(a * b * c, a * c), a * c));
			expect(eq(big_int::gcd(a * b * c, a * b), a * b));
			expect(eq(big_int::gcd(a * b * c, b * c), b * c));
			expect(eq(big_int::gcd(a * a, a * b), a));
			expect(eq(big_int::gcd(a * a, a * c), a));
			expect(eq(big_int::gcd(b * b, b * a), b));
			expect(eq(big_int::gcd(b * b, b * c), b));
			expect(eq(big_int::gcd(c * c, c * a), c));
			expect(eq(big_int::gcd(c * c, c * b), c));
		};
		feature("add_minus") = [] {
			auto a = "0b111111111111111111"_bi;
			auto b = "0b1111111111111111111"_bi;
			auto c = "0b1000000000000000000000000000000000000000000000000000000000000000000"_bi;
			auto d = "0b111111111111111111111111111111111111111111111111111111111111111111"_bi;
			expect(eq(a - a, 0_i));
			expect(eq(b - b, 0_i));
			expect(eq(c - 1, d));
			expect(eq(d + 1, c));
			expect(eq(d + 8, c + 7));
			expect(eq(-c + 1, -d));
			expect(eq(-c + d, -1));
			expect(eq(d - c, -1));
			expect(eq(-d + c, 1));
			expect(eq(c - d, 1));
			expect(eq(c + d, (c << 1) - 1));
			a += a;
			expect(eq(a, b - 1));
			a -= a;
			expect(eq(a, 0));
		};
		feature("string") = [] {
			auto test_string = [](std::string dec, std::string bin, std::string oct, std::string hex) {
				expect(big_int::from_string(dec).has_value());
				expect(eq(big_int::from_string(dec), big_int::from_string(bin)));
				expect(eq(big_int::from_string(dec), big_int::from_string(oct)));
				expect(eq(big_int::from_string(dec), big_int::from_string(hex)));
				expect(eq(big_int::from_string(dec)->to_string(big_int::base::decimal), dec));
				expect(eq(big_int::from_string(bin)->to_string(big_int::base::binary), bin));
				expect(eq(big_int::from_string(oct)->to_string(big_int::base::octal), oct));
				expect(eq(big_int::from_string(hex)->to_string(big_int::base::hex), hex));
			};
			test_string("18658316586846546587789749"
				, "0b111101101111000011010101111001110001011110101110000101010011101101111010100110110101"
				, "0o7557032571613656052355724665"
				, "0xf6f0d5e717ae153b7a9b5");
			test_string("33215588430304983195541539375682464839985654461394098514621059173680386919467827771862952765123987087013273061130144991703574128633748404479096130602468054549085269859006532267535358492450584752280012118799821038295237387310743703453485208812716091333318782521196132832473357275721628178731752174526818361800007516545562521391749"
				, "0b101000000100001010110011111011111001010111110110111011010011111101001101011010101000000100011000010100101111000101000011111111100011001001001111001010010101000010011010001101100110000011111101111001110001101110011100101111100011111101101110010011110110000011101100100111111010011101011011101110000111101001001001000000110000111010101110110101000010101111101111000010000101001100110111110100111000101010110010000000111100010111100111000010110010110101011110010111111111111101011000000001110110001011010100010100011010100011010001110001010011111110100011101101001110001100100011100111111111111101111111010001011111101101001101101111000000010111110000011011101010010011001000001110011001100011000010100001100101110000010010100011101001001011000111101110100000011000011000111011101011001011100011100111111110111111000110111011110010001110001001100011101001111111110000011001010111101001011001111111011011001000101010001011000111111111000111111101101101101100100000001011001001110100011101101011011100110001001000010011011000010111001010001110001111111010010000010111101101001000000101001010000101"
				, "0o5004126373712766732375153250043024570503774311171225023215460375716156345743755623660354477235335607511100607256650257570205146764705262007427470262653627777530016613242432432161237643551614434777757721375515570027603352231016314302414560224351130756403030735313434776770673621611435177603127513177331052130777077555544013116435533461102330271216177220275510051205"
				, "0xa042b3ef95f6ed3f4d6a811852f143fe324f29509a3660fde71b9cbe3f6e4f60ec9fa75bb87a49030eaed42bef085337d38ab203c5e70b2d5e5fff580762d451a8d1c53fa3b4e3239fff7f45fb4dbc05f06ea4c83998c2865c128e92c7ba0618eeb2e39fefc6ef23898e9ff0657a59fdb22a2c7fc7f6db202c9d1dadcc484d85ca38fe905ed205285");
		};
		feature("shift") = [] {
			for (int k = 1; k <= 128; k++) {
				big_int a = "18658316586846546587789749"_bi * 6416896889068409808_bi, b = a;
				auto c = 1_bi << k;
				for (int i = 0; i < 1000; i++) {
					a <<= k;
					b *= c;
					expect(eq(a, b));
				}
				for (int i = 0; i <= 2000; i++) {
					a >>= k;
					b /= c;
					expect(eq(a, b));
				}
				expect(a == 0_i);
				expect(b == 0_i);
			}
		};
		feature("umul128") = [] {
			std::mt19937_64 mt(std::mt19937_64::default_seed);
			using dist_t = std::uniform_int_distribution<uint64_t>;
			dist_t dist;
			for (int i = 0; i < 10000000; i++) {
				uint64_t a = dist(mt), b = dist(mt);
				uint64_t c, d;
				uint64_t e = spaceless::detail::umul128(a, b, &c);
				uint64_t f = spaceless::detail::umul128_naive_impl(a, b, &d);
				expect(eq(c, d));
				expect(eq(e, f));
			}
		};
		feature("udiv128") = [] {
			std::mt19937_64 mt(std::mt19937_64::default_seed);
			using dist_t = std::uniform_int_distribution<uint64_t>;
			dist_t dist;
			for (int i = 0; i < 10000000; i++) {
				uint64_t a = dist(mt), b = dist(mt), c = 0;
				while (a >= c) c = dist(mt);
				uint64_t d, e;
				uint64_t f = spaceless::detail::udiv128(a, b, c, &d);
				uint64_t g = spaceless::detail::udiv128_naive_impl(a, b, c, &e);
				expect(eq(d, e));
				expect(eq(f, g));
			}
		};
		feature("compare_with_python") = [] {
			std::ifstream fin("test_cases.txt");
			std::string op;
			while (std::getline(fin, op)) {
				std::string at, bt, ct;
				std::getline(fin, at);
				std::getline(fin, bt);
				std::getline(fin, ct);
				big_int a = big_int::from_string(at).value();
				big_int b = big_int::from_string(bt).value();
				big_int c = big_int::from_string(ct).value();
				auto calc = [&op](auto &&a, auto &&b) {
					if (op == "ADD")
						return std::forward<decltype(a)>(a) + std::forward<decltype(b)>(b);
					if (op == "SUB")
						return std::forward<decltype(a)>(a) - std::forward<decltype(b)>(b);
					if (op == "MUL")
						return std::forward<decltype(a)>(a) * std::forward<decltype(b)>(b);
					if (op == "DIV")
						return std::forward<decltype(a)>(a) / std::forward<decltype(b)>(b);
					if (op == "MOD")
						return std::forward<decltype(a)>(a) % std::forward<decltype(b)>(b);
					return big_int(0);
				};
				auto calc_eq = [&op](auto &&a, auto &&b) {
					if (op == "ADD")
						return std::forward<decltype(a)>(a) += std::forward<decltype(b)>(b);
					if (op == "SUB")
						return std::forward<decltype(a)>(a) -= std::forward<decltype(b)>(b);
					if (op == "MUL")
						return std::forward<decltype(a)>(a) *= std::forward<decltype(b)>(b);
					if (op == "DIV")
						return std::forward<decltype(a)>(a) /= std::forward<decltype(b)>(b);
					if (op == "MOD")
						return std::forward<decltype(a)>(a) %= std::forward<decltype(b)>(b);
					return big_int(0);
				};
				expect(eq(c, calc(a, b)));
				expect(eq(c, calc(big_int(a), b)));
				expect(eq(c, calc(a, big_int(b))));
				expect(eq(c, calc(big_int(a), big_int(b))));
				auto test_calc_eq = [&](auto &&a, auto &&b) {
					auto a_bk = a, b_bk = b;
					expect(eq(c, calc_eq(std::forward<decltype(a)>(a), std::forward<decltype(b)>(b))));
					a = a_bk;
					b = b_bk;
				};
				test_calc_eq(a, b);
				test_calc_eq(big_int(a), b);
				test_calc_eq(a, big_int(b));
				test_calc_eq(big_int(a), big_int(b));
			}
		};
	};
}

int main(int argc, char **argv) {
	unit_test();
	return 0;
}
