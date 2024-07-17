#include "big_int.h"
#define ANKERL_NANOBENCH_IMPLEMENT
#include <numbers>
#include <variant>
#include <cstdlib>

#include "external/nanobench.h"

using namespace spaceless;
using namespace big_int_literals;

static void benchmark() {
	using namespace ankerl::nanobench;
	Bench bench;
	bench.minEpochTime(std::chrono::milliseconds(100));

	std::vector bench_sizes{ 30, 100, 300, 1000, 3000, 10000 };
	auto bench_sizes_for_add_minus = bench_sizes | std::ranges::views::transform([](int x) { return x * x; });

	for (int i : bench_sizes_for_add_minus)
		bench.complexityN(i).run("add", [&]() {
			auto a = big_int::random(i), b = big_int::random(i);
			doNotOptimizeAway(a + b);
		});
	for (int i : bench_sizes_for_add_minus)
		bench.complexityN(i).run("minus", [&]() {
			auto a = big_int::random(i), b = big_int::random(i);
			doNotOptimizeAway(a - b);
		});

	for (int i : bench_sizes)
		bench.complexityN(i).run("multiply", [&]() {
			auto a = big_int::random(i), b = big_int::random(i);
			doNotOptimizeAway(a * b);
		});
	for (int i : bench_sizes)
		bench.complexityN(i).run("divide", [&]() {
			auto a = big_int::random(i), b = big_int(0);
			while (!b) b = big_int::random(i / 2);
			doNotOptimizeAway(a / b);
		});
	for (int i : bench_sizes)
		bench.complexityN(i).run("mod", [&]() {
			auto a = big_int::random(i), b = big_int(0);
			while (!b) b = big_int::random(i / 2);
			doNotOptimizeAway(a % b);
		});
}

int main(int argc, char **argv) {
	return 0;
}
