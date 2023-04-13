#include <iostream>
#include <unordered_set>
#include <vector>
#include <stack>
#include <ctime>
#include <fstream>
using namespace std;
using namespace std::chrono;

//int gcd(int a, int b) {
//    // Euclidean Algorithm
//    int h = max(a, b);
//    int l = min(a, b);
//    while (l > 0) {
//        int r = h % l;
//        h = l;
//        l = r;
//    }
//    return h;
//}
//
//bool is_power(long long n, int b) {
//    while (n >= b && n % b == 0) {
//        n /= b;
//    }
//    return (n == 1);
//}
//
//bool modRestriction(long long m, int p) {
//    // px+1 restriction
//    return (gcd(m, p-1) == 1);
//}
//
//void generate_primes(int n, vector<int> &primes) {
//    // Sieve of Eratosthenes
//    vector<bool> isPrime(n+1);
//    for (int i=0; i<=n; i++) {
//        isPrime[i] = true;
//    }
//    for (int p=2; p*p<=n; p++) {
//        // If isPrime[p] is not changed, then it is a prime
//        if (isPrime[p]) {
//            // Update all multiples of p
//            for (int i=p*2; i<=n; i+=p) {
//                isPrime[i] = false;
//            }
//        }
//    }
//    for (int i=0; i<=n; i++) {
//        if (isPrime[i] && i > 2) {
//            primes.push_back(i);
//        }
//    }
//}
//
//void pfilter(vector<int> &newPrimes, vector<int> &primes, vector<int> &factors) {
//    for (int p: primes) {
//        bool valid = true;
//        for (int f: factors) {
//            if (p == f || gcd(p, f-1) != 1 || gcd(p-1, f) != 1) {
//                valid = false;
//            }
//        }
//        if (valid) {
//            newPrimes.push_back(p);
//        }
//    }
//}
//
//struct Node {
//    long long m; // double check everything not gonna be an overflow issue (>2 billion)
//    long long phi;
//    float h_c;
//    int w; // omega
//    vector<int> possiblePrimes;
//
//    Node(long long m_in, long long p, float h, int w_in, vector<int> &pps) {
//        m = m_in;
//        phi = p;
//        h_c = h;
//        w = w_in;
//        possiblePrimes = pps; // possibly sus
//    }
//
////    Node* generate_left_child() {
////        return new Node(m, phi, h_c, w, index+1);
////    }
//
////    Node* generate_right_child(int q) {
////        // check px+1 restriction
////        if (modRestriction(m, q)) {
////            return new Node(m*q, phi*(q-1), h_c*((float)(q-1)/(float)q)*((float)(q-1)/(float)q), w+1, index+1);
////        }
////        else {
////            return nullptr;
////        }
////    }
//
//    void generate_child_primes(vector<int> &primes, vector<int> &childPrimes, int q, int max_w) {
//        float h_rest = 2.0*((((float)(max_w-w))/(h_c-1.0))+1.0);
//        for (int p: primes) {
//            if (p > q && (p < h_rest || h_rest < 0) && modRestriction(m,p)) {
//                childPrimes.push_back(p);
//            }
//        }
//    }
//
//    Node* generate_child(vector<int> &primes, int max_w) {
//        // last one in vector
//        int q = possiblePrimes[possiblePrimes.size()-1];
//        possiblePrimes.pop_back();
//        vector<int> childPrimes;
//        generate_child_primes(primes, childPrimes, q, max_w);
//        return new Node(m*q, phi*(q-1), h_c*((float)(q-1)/(float)q)*((float)(q-1)/(float)q), w+1, childPrimes);
//    }
//};
//
//struct Solution {
//    long long m;
//    long long phi;
//    long long k;
//    long long c;
//
//    Solution(long long m_in, long long p, long long k_in, long long c_in) {
//        m = m_in;
//        phi = p;
//        k = k_in;
//        c = c_in;
//    }
//
//    void print() {
//        cout << "n = " << c << " * " << m << " for k = " << k << endl;
//    }
//};
//
//bool ngreater(int x, float y) {
//    // the kyle mystery function
//    return (y>0 && (float)x > y);
//}
//
//void search_c(long long m, long long phi, int min_c, vector<Solution*> &sols) {
//    long long lbound = ceil(((float)m/(float)phi));
//    long long ubound = floor(((float)min_c/(float)(min_c-1))*((float)m/(float)phi));
//    for (long long k=lbound; k<=ubound; k++) {
//        if ((k*phi - 1) % (k*phi - m) == 0) {
//            long long c = (k*phi - 1)/(k*phi - m);
//            Solution* sol = new Solution(m, phi, k, c);
//            sols.push_back(sol);
//            sol->print();
//        }
//    }
//}
//
//void dfs(Node* root, int numPrimes, vector<int> &primes, int max_w, int min_c, vector<Solution*> &sols) {
//    // this primes is the newPrimes list in main()
//    long long count = 0;
//    stack<Node*> s;
//    s.push(root);
//    while (!s.empty()) {
//        Node* curr = s.top();
//        // SMALL IF
//        if (curr->w == max_w || curr->possiblePrimes.empty()) {
//            count++; // sus?
//            if (curr->m > 1) {
//                search_c(curr->m, curr->phi, min_c, sols);
////                if (is_power(count, 10)) {
////                    cout << "   current count: 10^" << log10(count) << endl;
////                }
//                if (count % 100000 == 0) {
//                    cout << "    Checking " << count << " node..." << endl;
//                }
//            }
//            s.pop();
//        }
//        else {
//            Node* child = curr->generate_child(primes, max_w);
//            s.push(child);
//        }
//    }
//    cout << "count = " << count << endl;
//}
//
//int main() {
//    auto start = high_resolution_clock::now(); // start time
//    // initialize input parameters
//    long long numPrimesInit;
//    int max_w, min_c;
//    vector<int> factors; // prime factors guaranteed to have in final m value
//
//    cin >> numPrimesInit;
//    cin >> max_w, cin >> min_c;
//    int num_fac;
//    cin >> num_fac;
//    for (int i=0; i<num_fac; i++) {
//        int fac;
//        cin >> fac;
//        factors.push_back(fac);
//    }
//
//    // initialize root node
//    long long m_root = 1;
//    long long phi_root = 1;
//    float sqrt = 2.0*((float)(min_c-1)/(float)min_c);
//    float h_root = sqrt*sqrt;
//    int w_root = 0;
//
//    for (int p: factors) {
//        m_root *= p;
//        phi_root *= (p-1);
//        sqrt = ((float)(p-1)/(float)p);
//        h_root *= sqrt*sqrt;
//        w_root++;
//    }
//
//    // generate list of primes and do p-filter
//    vector<int> primes;
//    generate_primes(numPrimesInit, primes);
//    vector<int> newPrimes; // apply p-filter
//    pfilter(newPrimes, primes, factors);
//    int numPrimes = newPrimes.size();
//    Node* root = new Node(m_root, phi_root, h_root, w_root, newPrimes);
//    // could add hrest
//
//    vector<Solution*> sols;
//    dfs(root, numPrimes, newPrimes, max_w, min_c, sols);
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//    cout << "The program terminated after " << duration.count()/1000000000.0 << " seconds" << endl;
//    return 0;
//}

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

    void print(ofstream &myfile, time_t start) { // do need to pass time_t by ref?
        time_t curr = time(NULL);
        myfile << "n = " << c << " * " << m << " for k = " << k << ", in " << curr-start << " seconds" << endl;
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

void search_c(long long m, long long phi, float c_frac, vector<Solution*> &sols, ofstream &myfile, time_t start) {
    long long lbound = ceil(((float)m / (float)phi));
    long long ubound = floor(c_frac * ((float)m / (float)phi));
    for (long long k = lbound; k <= ubound; k++) {
        if ((k * phi - 1) % (k * phi - m) == 0) {
            long long c = (k * phi - 1)/(k * phi - m);
            Solution* sol = new Solution(m, phi, k, c);
            sols.push_back(sol);
            sol->print(myfile, start);
        }
    }
}

void dfs(Node* root, int numPrimes, vector<int> &primes, int max_w, float c_frac, vector<Solution*> &sols, time_t start) {
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
                search_c(curr->m, curr->phi, c_frac, sols, myfile, start);
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
    dfs(root, numPrimes, primes, max_w, c_frac, sols, start);

    // auto stop = high_resolution_clock::now();
    time_t end = time(NULL);
    // auto duration = duration_cast<seconds>(stop - start);
    cout << "The program terminated after " << end - start << " seconds" << endl;
    return 0;
}
