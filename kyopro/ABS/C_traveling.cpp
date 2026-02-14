#include <bits/stdc++.h>
using namespace std;

bool canMove(int t1, int t2, int x1, int x2, int y1, int y2) {
    int dt = t2 - t1;
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    // 2秒余ればループできる
    if (dt < dx + dy) return false;
    if ((dt - (dx + dy)) % 2) return false;
    return true;
}

int main() {
    int n;
    cin >> n;
    vector<int> t(n), x(n), y(n);
    for (int i = 0; i < n; i++) {
        cin >> t.at(i) >> x.at(i) >> y.at(i);
    }

    // スタート
    int time = 0, x0 = 0, y0 = 0;
    // 1ステップごとに移動できるか判定できる関数を用意
    if (!canMove(time, t.at(0), x0, x.at(0), y0, y.at(0))) {
        cout << "No" << endl;
        return 0;
    } else if (n == 1) {
        cout << "Yes" << endl;
        return 0;
    } else {
        time = t.at(0);
    }
    for (int i = 1; i < n; i++) {
        if (canMove(time, t.at(i), x.at(i-1), x.at(i), y.at(i-1), y.at(i))) {
            time = t.at(i);
        }else {
            cout << "No" << endl;
            return 0;
        }
    }
    cout << "Yes" << endl;
}