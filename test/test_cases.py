import random
import sys
import tqdm

def random_bigint(bits):
    return random.randint(-2 ** bits, 2 ** bits)

def to_str(x):
    return random.choice([str, bin, oct, hex])(x)

def output(f, x):
    f.write(f"{to_str(x)}\n")
    
def output_result(f, op, a, b, c):
    f.write(f"{op}\n")
    output(f, a)
    output(f, b)
    output(f, c)

def gen(n = 10000):
    with open("test_cases.txt", "w") as f:
        for _ in tqdm.tqdm(range(n)):
            bits = random.choice([10, 35, 65, 257, 1025, 10000, 100000])
            a = random_bigint(bits)
            b = random_bigint(bits)
            output_result(f, "ADD", a, b, a + b)
            output_result(f, "SUB", a, b, a - b)
            output_result(f, "MUL", a, b, a * b)
            c = random_bigint(bits // 2)
            if c:
                d = abs(a) // abs(c)
                r = abs(a) % abs(c)
                d = -d if (a < 0) != (c < 0) else d
                r = -r if a < 0 else r
                output_result(f, "DIV", a, c, d)
                output_result(f, "MOD", a, c, r)

if __name__ == "__main__":
    sys.set_int_max_str_digits(100000)
    gen()
