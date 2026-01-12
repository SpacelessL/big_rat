# Big Rat<span style="color:lightgrey">~~ional~~</span>

A **WIP** C++23 single/double-header(s) library aimed at handling arbitrary-sized integer and rational numbers.
Though it currently only uses $O(n^{2})$ algorithms for multiplication and division, its performance can be comparable to [GMP](https://gmplib.org/) for $numbers \leq 2^{10000}$. This library lacks many unit tests, especially for the big rational part. It's still in its extremely early stages, so don't use it in any (especially production) environments.

## Features
1. **Readable and Fast**: The aim is to maintain code readability while ensuring efficiency. This library is meant to be a more accessible alternative to others like mini-GMP, balancing speed with simplicity.
2. **Modern C++ Features**: Move semantics, constexpr, attributes, user-defined literals. This library is also a playground for these modern C++ features.
3. **Versatile Conversion**: Supports conversion from/to binary, octal, hexadecimal, and decimal. Plus, `big_rat` can handle repeating decimals (e.g., $0.\overline{142857} \to \frac{1}{7}$) and convert from IEEE 754 floating points without any precision loss.

## Usage
```cpp
#include "big_int.h"
using namespace spaceless;

// Construction
big_int a = 12345678901234567890_bi;
big_int b = big_int::from_string("FF", big_int::base::hex).value();
big_int c("123456789012345678901234567890");

// Arithmetic
big_int result = a * b + c;
big_int fact = big_int::factorial(100);
big_int exp = a.power(10);

// Number theory
big_int g = big_int::gcd(a, b);
big_int l = big_int::lcm(a, b);

// Bitwise
size_t bits = a.highest_bit();
size_t ones = a.popcount();
a.bitwise_not();                  // in-place

// Sign (in-place)
a.negate();                       // flip sign
a.absolutize();                   // make positive

// String conversion
std::string hex = a.to_string(big_int::base::hex);
std::string bin = a.to_string(big_int::base::binary);
```

```cpp
#include "big_rat.h"
using namespace spaceless;

// Construction
big_rat a = 3.14_br;              // lossless from literal
big_rat b(3, 4);                  // 3/4
big_rat c = big_rat(0.1);         // exact IEEE 754 representation

// Arithmetic
big_rat result = a * b / c;
big_rat pow = b.power(-2);        // 16/9

// In-place operations
b.invert();                       // 3/4 -> 4/3
b.negate();                       // 4/3 -> -4/3

// Components
big_int num = result.numerator();
big_int den = result.denominator();
big_int truncated = result.trunc();

// Approximation
big_rat pi_approx = big_rat(3.14159265358979).limit_denominator(1000); // 355/113

// Repeating decimals: 0.142857142857... = 1/7
auto rd = big_rat::repeating_decimal::from_string("0.(142857)");
big_rat seventh = big_rat::from_repeating_decimal(rd).value();
```

## Todo
1.  Expand unit testing and include code coverage.
2.  Transition to my own bench and ut framework in [sus](https://github.com/SpacelessL/sus), when they are implemented.
3.  Add benchmark tables.
4.  Implement floor, ceil, and round functions in `big_rat`.
5.  Enhance `big_rat` to floating-point conversion.
6.  Refactor duplicated binary operators in `big_rat` for cleaner code.
7.  Implement multiplication algorithms (Karatsuba, Toom–Cook, FFT), and improve division as well.
8.  Integrate exception support.
9.  Include more comments for clarity.
10. Multi-platform support.

# License

This project is licensed under the MIT License - see the LICENSE file for details. 
