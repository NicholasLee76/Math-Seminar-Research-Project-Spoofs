#include <iostream>
#include <math.h>
using namespace std;

int sigma(int m) {
    if (m == 1) {
        return 1;
    }
    int s = m+1;
    for (int i=2; i<m; i++) {
        if (m % i == 0) {
            s += i;
        }
    }
    return s;
}

bool check_quasi(int m, int s, int c, bool check1) {
    bool check = s*(c*c+c+1) == 2*m*c*c+1;
    if (check1) return check;
    return check && c > 1;
}

bool check_almost(int m, int s, int c, bool check1) {
    bool check = s*(c*c+c+1) == 2*m*c*c-1;
    if (check1) return check;
    return check && c > 1;
}

int main() {
    int limit = 10000000;
    bool check1 = true;
    for (int m=1; m<limit; m++) {
        // n = mc^2
        if (m % 200000 == 0) {
            cout << "Checking " << m << "... " << endl;
        }
        int s = sigma(m);
        int c_top = ceil(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));
        int c_bot = floor(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));

        // quasi
        if (check_quasi(m, s, c_top, check1)) {
            cout << "c = " << c_top << ", m = " << m << ", quasi" << endl;
        }

        if (c_top != c_bot && check_quasi(m, s, c_bot, check1)) {
            cout << "c = " << c_bot << ", m = " << m << ", quasi" << endl;
        }

        // almost
        if (check_almost(m, s, c_top, check1)) {
            cout << "c = " << c_top << ", m = " << m << ", almost" << endl;
        }

        if (c_top != c_bot && check_almost(m, s, c_bot, check1)) {
            cout << "c = " << c_bot << ", m = " << m << ", almost" << endl;
        }
    }

    return 0;
}
