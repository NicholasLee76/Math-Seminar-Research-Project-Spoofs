#include <iostream>
#include <math.h>
using namespace std;

int sigma(int m) {
    int s = m+1;
    for (int i=2; i<m; i++) {
        if (m % i == 0) {
            s += i;
        }
    }
    return s;
}

bool check_quasi(int m, int s, int c) {
    return s*(c*c+c+1) == 2*m*(c+1)+1 && c > 1;
}

bool check_almost(int m, int s, int c) {
    return s*(c*c+c+1) == 2*m*(c+1)-1 && c > 1;
}

int main() {
    int limit = 10000000;
    for (int m=1; m<limit; m++) {
        if (m % 200000 == 0) {
            cout << "Checking " << m << "... " << endl;
        }
        int s = sigma(m);
        int c_top = ceil(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));
        int c_bot = floor(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));

        // quasi
        if (check_quasi(m, s, c_top)) {
            cout << "c = " << c_top << ", m = " << m << ", quasi";
        }

        if (check_quasi(m, s, c_bot)) {
            cout << "c = " << c_bot << ", m = " << m << ", quasi";
        }

        // almost
        if (check_almost(m, s, c_top)) {
            cout << "c = " << c_top << ", m = " << m << ", almost";
        }

        if (check_almost(m, s, c_bot)) {
            cout << "c = " << c_bot << ", m = " << m << ", almost";
        }
    }

    return 0;
}
