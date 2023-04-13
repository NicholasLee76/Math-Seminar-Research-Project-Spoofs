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

int main() {
    int limit = 10000000;
    for (int m=1; m<limit; m++) {
        if (m % 200000 == 0) {
            cout << "Checking " << m << "th value" << endl;
        }
        int s = sigma(m);
        int c_top = ceil(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));
        int c_bot = floor(((2*m-s) + sqrt((s-2*m)*(s-2*m) - 4*s*(s-2*m+1)))/(2*s));

        // quasi
        if (s*(c_top*c_top+c_top+1) == 2*m*(c_top+1)+1 && c_top > 1) {
            cout << "c = " << c_top << ", m = " << m << ", quasi";
        }

        if (s*(c_bot*c_bot+c_bot+1) == 2*m*(c_bot+1)+1 && c_bot > 1) {
            cout << "c = " << c_bot << ", m = " << m << ", quasi";
        }

        // almost
        if (s*(c_top*c_top+c_top+1) == 2*m*(c_top+1)-1 && c_top > 1) {
            cout << "c = " << c_top << ", m = " << m << ", almost";
        }

        if (s*(c_bot*c_bot+c_bot+1) == 2*m*(c_bot+1)-1 && c_bot > 1) {
            cout << "c = " << c_top << ", m = " << m << ", almost";
        }
    }

    return 0;
}
