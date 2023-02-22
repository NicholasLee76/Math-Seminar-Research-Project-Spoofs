import sympy as sp
import numpy as np
import math
import timeit as tt


def searchRatio(m, tot, cmin):
    # Searches for values of c
    sols = []
    for k in range(math.ceil(m / tot),  # Bounds on k
                   math.floor((cmin / (cmin - 1)) * (m / tot) - 1 / ((cmin - 1) * tot)) + 1):
        if not (k * tot - 1) % (k * tot - m):  # Checks whether the second expression divides the first
            sols.append([m, tot, k, (k * tot - 1) // (k * tot - m)])
    return sols


def pfilter(primes, m):
    # Filters out primes which violate the px+1 restriction (as well as the prime factors of m)
    return [p for p in primes if math.gcd(m, p - 1) == 1 and m % p]


def ngreater(x, y):
    return y > 0 and x > y


def timedTest(n, factors=[], omega=4, cmin=3):
    tinit = tt.default_timer()

    sols = []
    masks = 0
    hc = (2 * (cmin - 1) / cmin) ** 2

    row = [[1, 1, hc, 0, 0]]  # Format: m, phi(m), h rest value***, omega(m), len

    for p in factors:
        row[0][0] *= p
        row[0][1] *= p - 1
        row[0][2] *= ((p - 1) / p) ** 2
        row[0][3] += 1

    primes = pfilter(list(sp.primerange(3, sp.prime(n))), row[0][0])
    size = len(primes)

    trow = tt.default_timer()

    tqconds = []
    tbigifs = []
    tsearchorappends = []
    tsearches = []
    tsearches1 = []
    tappends = []
    tappends1 = []

    while row:
        q = 0
        x = row.pop()

        tqconds.append(tt.default_timer())

        if x[4] < size:
            q = primes[x[4]]

        tbigifs.append(tt.default_timer())

        # Max prime restriction, max omega restriction, H restriction
        if not q or x[3] == omega or ngreater(q, 2 * ((omega - x[3]) / (x[2] - 1) + 1)):
            tsearchorappends.append(tt.default_timer())
            if x[0] > 1:  # Filters out m = 1
                tsearches1.append(tt.default_timer())

                c = searchRatio(x[0], x[1], cmin)
                if c:
                    sols.append(c)
                    print(c)
                tsearches.append(tt.default_timer())
                masks += 1
        else:
            tsearchorappends.append(tt.default_timer())

            tappends1.append(tt.default_timer())
            row.append([x[0], x[1], x[2], x[3], x[4] + 1])

            if math.gcd(x[0], q - 1) == 1:
                row.append([x[0] * q, x[1] * (q - 1), x[2] * ((q - 1) / q) ** 2, x[3] + 1, x[4] + 1])
            tappends.append(tt.default_timer())

    tfin = tt.default_timer()

    return sols, masks, trow + tinit, tfin - trow, np.subtract(tbigifs, tqconds), np.subtract(tsearchorappends,
                                                                                              tbigifs), np.subtract(
        tsearches, tsearches1), np.subtract(tappends, tappends1)


b = timedTest(100, factors=[3, 5, 17], omega=6)

print(f'Total Algorithm run time: {b[2] + b[3]} seconds')
print(f'    Initializing: {b[2]} seconds')
print(f'    Calculation: {b[3]} seconds for {b[1]} masks (Average of {b[3] / b[1]} seconds per mask)')
print(f'        Qcond evaluation: {np.sum(b[4])} (Average of {np.mean(b[4])} seconds per mask)')
print(f'        Big if statement evaluation: {np.sum(b[5])} (Average of {np.mean(b[5])} seconds per mask)')
print(f'        C searching: {np.sum(b[6])} (Average of {np.mean(b[6])} seconds per mask)')
print(f'        Appending values: {np.sum(b[7])} (Average of {np.mean(b[7])} seconds per mask)')
print(b[0])
