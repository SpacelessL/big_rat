#include "big_rat.h"
#include "external/ut.hpp"
#include <numbers>

using namespace boost::ut;
using namespace bdd;
using namespace spaceless;
using namespace big_int_literals;
using namespace big_rat_literals;

static void unit_test() {
	big_rat e = "0.(010309278350515463917525773195876288659793814432989690721649484536082474226804123711340206185567)"_br;
	//expect(eq(e.numerator(), 1));
	//expect(eq(e.denominator(), 97));
	e = -"3.1415926535897932"_br;
	std::cout << e.to_string() << std::endl;
	auto e1 = e.limit_denominator(1000);
	auto e2 = e.limit_denominator(100);
	auto e3 = e.limit_denominator(10);
	std::cout << e1.to_string() << std::endl;
	std::cout << e2.to_string() << std::endl;
	std::cout << e3.to_string() << std::endl;
	std::cout << "-1.1"_br.to_string() << std::endl;
	std::cout << "-1.1"_br.limit_denominator().to_string() << std::endl;
	std::cout << big_rat(-std::cos(std::numbers::pi / 3)).to_string() << std::endl;
	std::cout << big_rat(-std::cos(std::numbers::pi / 3)).limit_denominator().to_string() << std::endl;
}

int main(int argc, char **argv) {
	unit_test();
	return 0;
}
