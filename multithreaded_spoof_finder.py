import sympy as sp
import numpy as np
import math
import timeit as tt
import threading
import queue

# start queues
input_queue = queue.Queue()
c_queue = queue.Queue()
output_queue = queue.Queue()


def generate_queue(qu):
    while True:
        data = qu.get()
        row = data[0]
        omega = data[1]
        cmin = data[2]
        primes = data[3]
        q = 0
        if row:
            x = row.pop()
        else:
            print('bruh')
            check = False
            return
        if x[4] < len(primes):
            q = primes[x[4]]
        if not q or x[3] == omega or ngreater(q, 2 * ((omega - x[3]) / (x[2] - 1) + 1)):
            if x[0] > 1:
                c_queue.put([x[0], x[1], cmin])
        else:
            row.append([x[0], x[1], x[2], x[3], x[4] + 1])
            if math.gcd(x[0], q - 1) == 1:
                row.append([x[0] * q, x[1] * (q - 1), x[2] * ((q - 1) / q) ** 2, x[3] + 1, x[4] + 1])
            input_queue.put([row, omega, cmin, primes])

        qu.task_done()


def search_thread(q):
    while True:
        data = q.get()
        sols = []
        m = data[0]
        tot = data[1]
        cmin = data[2]
        for k in range(math.ceil(m / tot),  # Bounds on k
                       math.floor((cmin / (cmin - 1)) * (m / tot) - 1 / ((cmin - 1) * tot)) + 1):
            if not (k * tot - 1) % (k * tot - m):  # Checks whether the second expression divides the first
                sols.append([m, tot, k, (k * tot - 1) // (k * tot - m)])
        output_queue.put(sols)
        q.task_done()


def output_thread(q):
    while True:
        c = output_queue.get()
        if c:
            print(c)
        q.task_done()


def pfilter(primes, m):
    # Filters out primes which violate the px+1 restriction (as well as the prime factors of m)
    return [p for p in primes if math.gcd(m, p - 1) == 1 and m % p]


def ngreater(x, y):
    return y > 0 and x > y


def main(n, factors, omega=4, cmin=3):
    # init step
    t_init = tt.default_timer()

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

    t1 = threading.Thread(target=generate_queue, args=(input_queue,), daemon=True).start()
    t2 = threading.Thread(target=search_thread, args=(c_queue,), daemon=True).start()
    t3 = threading.Thread(target=output_thread, args=(output_queue,), daemon=True).start()

    global check
    check = True
    # generate and search steps
    while True:
        input_queue.put([row, omega, cmin, primes])  # goes into gen thread
        if not check:
            break

    t1.join()
    t2.join()
    t3.join()
    t_final = tt.default_timer()
    print(f"total runtime = {t_final}")
    print(f"init runtime = {t_init}")
    print(f"algorithm runtime = {t_final-t_init}")


if __name__ == '__main__':
    main(1000, factors=[3, 5, 17], omega=6)
