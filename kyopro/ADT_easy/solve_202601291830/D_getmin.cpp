#include <bits/stdc++.h>
using namespace std;

int main() {
    int Q;
    cin >> Q;

    vector<int> balls(0);
    for (int i = 0; i < Q; i++) {
        int type;
        cin >> type;
        if (type == 1) {
            int x;
            cin >> x;
            balls.push_back(x);
            sort(balls.begin(), balls.end());
        }
        else if (type == 2) {
            cout << balls.at(0) << endl;
            balls.erase(balls.begin());
        }
    }
}

// getlineを使わずに、逐次入力を受け取る方式に変更