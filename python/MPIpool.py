# MPIpool.py

from mpi4py.futures import MPIPoolExecutor
import math
import textwrap

def calc_primes(range_tuple):
    """Calculate all the prime numbers in the given range."""
    low, high = range_tuple
    print(low)
    if low <= 2 < high:
        primes = [2]
    else:
        primes = []

    start = max(3,low)   # Don't start below 3
    if start % 2 == 0:   # Make sure start is odd, i.e. skip evens
        start += 1

    for num in range(start, high, 2):  # increment by 2's, i.e. skip evens
        if all(num % i != 0 for i in range(3, int(math.sqrt(num)) + 1, 2)):
            primes.append(num)

    return primes


def determine_subranges(fullrange, num_subranges):
    """
    Break fullrange up into smaller sets of ranges that cover all
    the same numbers.
    """
    subranges = []
    inc = fullrange[1] // num_subranges
    for i in range(fullrange[0], fullrange[1], inc):
        subranges.append( (i, min(i+inc, fullrange[1])) )
    return( subranges )


if __name__ == '__main__':
    # fullrange = (0, 10000000)
    fullrange = (0, 2000)
    num_subranges = 1000
    subranges = determine_subranges(fullrange, num_subranges)

    executor = MPIPoolExecutor()
    prime_sets = executor.map(calc_primes, subranges)
    executor.shutdown()

    # flatten the list of lists
    primes = [p for plist in prime_sets for p in plist]
    print(textwrap.fill(str(primes),80))