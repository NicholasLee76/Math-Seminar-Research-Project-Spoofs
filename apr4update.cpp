#include <iostream>
#include <unordered_set>
#include <vector>
#include <stack>
#include <ctime>
#include <fstream>
using namespace std;
using namespace std::chrono;

int gcd(int a, int b) {
    // Euclidean Algorithm
    int h = max(a, b);
    int l = min(a, b);
    while (l > 0) {
        int r = h % l;
        h = l;
        l = r;
    }
    return h;
}

bool modRestriction(long long m, int q) {
    // px+1 restriction
    return (gcd(m, q - 1) == 1);
}

void generate_filtered_primes(int n, vector<int> &primes, vector<int> &factors) {
    // Sieve of Eratosthenes
    vector<bool> isPrime(n + 1);
    for (int i = 0; i <= n; i++) isPrime[i] = true;

    for (int p = 2; p * p <= n; p++) {
        // if isPrime[p] is not changed, then it is a prime
        if (isPrime[p]) {
            // update multiples of p
            for (int i = 2 * p; i <= n; i += p) isPrime[i] = false;
        }
    }

    for (int i = 0; i <= n; i++) {
        for (int f : factors) {
            // pfilter
            if (isPrime[i] && (i == f || gcd(i, f - 1) != 1 || gcd(i - 1, f) != 1)) isPrime[i] = false;
        }

        // if p is an odd prime, add it to vector
        if (isPrime[i] && i > 2) {
            primes.push_back(i);
        }
    }
}

struct Node {
    long long m; // double check everything not gonna be an overflow issue (>2 billion)
    long long phi;
    float h_c;
    int w; // omega
    int min_index;
    float prime_max; // h_rest

    Node(long long m_in, long long p, float h, int w_in, int mini, int max_w) {
        m = m_in;
        phi = p;
        h_c = h;
        w = w_in;
        min_index = mini;

        prime_max = 2.0 * ((((float)(max_w - w)) / (h_c - 1.0)) + 1.0);
    }

    // When you call generate_child, one of two things happens:
    // (1) It loops through primes starting at min_index until it finds one that passes the modRestriction
    // (2) It fails the omega restriction, the hrestriction, or it runs out of primes and thus is terminal,
    // returning a nullptr instead of a new Node

    Node* generate_child(int numPrimes, vector<int> &primes, int max_w) {
        if (w == max_w) return nullptr;

        int q;
        do {
            if (min_index == numPrimes) return nullptr;

            q = primes[min_index];
            min_index++;

            if (prime_max > 0 && (float)q > prime_max)  return nullptr;
        } while (!modRestriction(m, q)); // keeps iterating until q passes px+1 restriction

        return new Node(m*q, phi*(q-1), h_c*((float)(q-1)/(float)q)*((float)(q-1)/(float)q), w+1, min_index, max_w);
    }
};

struct Solution {
    long long m, phi, k, c;

    Solution(long long m_in, long long p, long long k_in, long long c_in) {
        m = m_in;
        phi = p;
        k = k_in;
        c = c_in;
    }

    void print(ofstream &myfile) {
        myfile << "n = " << c << " * " << m << " for k = " << k << endl;
    }
};

Node* generate_root(vector<int> &factors, int max_w, float c_frac) {
    // initialize root node
    long long m_root = 1;
    long long phi_root = 1;
    float h_sqrt = 2.0 / c_frac;
    float h_root = h_sqrt * h_sqrt;
    int w_root = 0;

    for (int p : factors) {
        m_root *= p;
        phi_root *= (p-1);
        h_root *= ((float)(p-1) / (float)p) * ((float)(p-1) / (float)p);
        w_root++;
    }

    return new Node(m_root, phi_root, h_root, w_root, 0, max_w);
}

void search_c(long long m, long long phi, float c_frac, vector<Solution*> &sols, ofstream &myfile) {
    long long lbound = ceil(((float)m / (float)phi));
    long long ubound = floor(c_frac * ((float)m / (float)phi));
    for (long long k = lbound; k <= ubound; k++) {
        if ((k * phi - 1) % (k * phi - m) == 0) {
            long long c = (k * phi - 1)/(k * phi - m);
            Solution* sol = new Solution(m, phi, k, c);
            sols.push_back(sol);
            sol->print(myfile);
        }
    }
}

void dfs(Node* root, int numPrimes, vector<int> &primes, int max_w, float c_frac, vector<Solution*> &sols) {
    long long count = 0;
    long long maxStackSize = 0;
    stack<Node*> s;
    s.push(root);

    // write solutions.txt to file
    ofstream myfile;
    myfile.open("solutions.txt");

    while (!s.empty()) {
        Node* curr = s.top();
        Node* child = curr->generate_child(numPrimes, primes, max_w);

        if (s.size() > maxStackSize) maxStackSize = s.size();

        if (child == nullptr) {
            count++;
            if (curr->m > 1) {
                search_c(curr->m, curr->phi, c_frac, sols, myfile);
            }
            if (count % 10000000 == 0) {
                cout << "Checked " << count / 1000000 << " million nodes (mss of " << maxStackSize << ")..." << endl;
            }
            s.pop();
            delete(curr);
        }

        else {
            s.push(child);
        }
    }
    cout << "Total count: " << count << " nodes" << endl;
    cout << "Max stack size: " << maxStackSize << " dormant nodes" << endl;
    myfile.close();
}


int main() {
    // auto start = high_resolution_clock::now(); // start time
    // initialize input parameters
    long long maxPrime;
    int max_w, min_c;
    vector<int> factors; // prime factors to include in final m value
    int min_index_init = 0;
    cin >> maxPrime;
    cin >> max_w;
    cin >> min_c;
    int num_fac;
    cin >> num_fac;

    for (int i = 0; i < num_fac; i++) {
        int fac;
        cin >> fac;
        factors.push_back(fac);
    }

    time_t start = time(NULL);
    float c_frac = (float)min_c / (float)(min_c - 1);

    Node* root = generate_root(factors, max_w, c_frac);

    // generate list of primes and do p-filter
    vector<int> primes;
    generate_filtered_primes(maxPrime, primes, factors);
    int numPrimes = primes.size();

    cout << "---------------" << endl;

    cout << "Number of primes testing: " << numPrimes << endl;

    vector<Solution*> sols;
    dfs(root, numPrimes, primes, max_w, c_frac, sols);

    // auto stop = high_resolution_clock::now();
    time_t end = time(NULL);
    // auto duration = duration_cast<seconds>(stop - start);
    cout << "The program terminated after " << end - start << " seconds" << endl;
    return 0;
}
