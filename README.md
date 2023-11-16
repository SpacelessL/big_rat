# Big Rat<span style="color:lightgrey">~~ional~~</span>

A **WIP** C++23 single/double-header(s) library aimed at handling arbitrary-sized integer and rational numbers. 
Though it currently only uses $O(n^{2})$ algorithms for multiplication and division, its performance can be comparable to [GMP](https://gmplib.org/) for $numbers \leq 2^{10000}$. This library lacks many unit tests, especially for the big rational part. It's still in its extremely early stages, so don't use it in any (especially production) environments.

## Usage
To get started, simply include `big_int.h` or `big_rat.h` in your project.

## Features
1. **Readable and Fast**: The aim is to maintain code readability while ensuring efficiency. This library is meant to be a more accessible alternative to others like mini-GMP, balancing speed with simplicity.
2. **Modern C++ Features**: move semantics, constexpr, attributes, user-defined literals. This library is also a playground for these modern C++ features.
3. **Versatile Conversion**: Supports conversion from/to binary, octal, hexadecimal, and decimal. Plus, `big_rat` can handle repeating decimals (e.g., $0.\overline{142857} \to \frac{1}{7}$) and convert from IEEE 754 floating points without any precision loss.

## Todo
1.  Expand unit testing and include code coverage.
2.  Transition to Catch2 from ut++ and microbench.
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
