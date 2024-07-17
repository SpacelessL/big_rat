#include "big_rat.h"
#include <numbers>
#define ANKERL_NANOBENCH_IMPLEMENT
#include "external/nanobench.h"

using namespace spaceless;
using namespace big_int_literals;
using namespace big_rat_literals;

static big_rat calculate_pi() {
	big_rat ret = 0;
	/*for (uint64_t k = 0; k < 3; ++k) {
		big_rat x = big_int::factorial(6 * k) *(13591409 + 545140134 * k)
			/ (big_int::factorial(3 * k) *big_int::factorial(k).power(3) *big_rat((640320_bi).power(6 * k + 3)).sqrt());
		if (k % 2) ret -= x;
		else ret += x;
	}
	ret *= 12;*/
	return std::move(ret).inverse();
}

static void benchmark() {
	using namespace ankerl::nanobench;
	Bench bench;
	bench.minEpochTime(std::chrono::milliseconds(100));
	bench.run("calculate_pi", []() {
		doNotOptimizeAway(calculate_pi());
	});
}

int main(int argc, char **argv) {
	big_rat pi = calculate_pi();
	std::cout << pi << std::endl;
	std::cout << double(pi) << std::endl;
	std::cout << big_int::factorial(100) << std::endl;
	//std::cout << big_rat(2).sqrt() << std::endl;
	benchmark();
	return 0;
}
