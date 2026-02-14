#include <bits/stdc++.h>
using namespace std;

int main() {
    int in;
    cin >> in;
    int a = in / 1000;
    int b = (in / 100) % 10;
    int c = (in / 10) % 10;
    int d = in % 10;
    
    // 演算子を試す
    bool flag = false;
    char op1_char, op2_char, op3_char;
    for (int i = 0; i < 2; i++) { // op1
        for (int j = 0; j < 2; j++) { // op2
            for (int k = 0; k < 2; k++) { // op3
                int num = a;
                if (i == 0) {
                    num += b;
                } else {
                    num -= b;
                }
                if (j == 0) {
                    num += c;
                } else {
                    num -= c;
                }
                if (k == 0) {
                    num += d;
                } else {
                    num -= d;
                }
                if (num == 7) {
                    flag = true;
                }
                if (flag) {
                    if (i == 0) {
                        op1_char = '+';
                    } else {
                        op1_char = '-';
                    }
                    if (j == 0) {
                        op2_char = '+';
                    } else {
                        op2_char = '-';
                    }
                    if (k == 0) {
                        op3_char = '+';
                    } else {
                        op3_char = '-';
                    }
                    cout << a << op1_char << b << op2_char << c << op3_char << d << "=7" << endl;
                    break;
                }
            }
            if (flag) {
                break;
            }
        }
        if (flag) {
            break;
        }
    }
    
}