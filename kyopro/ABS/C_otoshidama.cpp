#include <bits/stdc++.h>
using namespace std;

int main() {
    int N, Y;
    cin >> N >> Y;

    // n枚ぴったりじゃないとだめだったわ
    int x = 0, y = 0, z = 0;
    for (int i = 0; i <= N; i++) {  // 10000の枚数
        for (int j = 0; j <= N - i; j++) {  // 5000の枚数
            x = i;
            y = j;
            z = N - i - j;
            if (10000*x + 5000*y + 1000*z == Y) {
                cout << x << " " << y << " " << z << endl;
                return 0;
            }
        }
    }
    cout << "-1 -1 -1" << endl;
}