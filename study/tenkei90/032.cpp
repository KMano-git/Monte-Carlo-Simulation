#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;
    vector<vector<int>> A(N, vector<int>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cin >> A.at(i).at(j);
        }
    }
    int M;
    cin >> M;
    vector<vector<bool>> bad(N, vector<bool>(N, false));
    for (int i = 0; i < M; i++) {
        int u, v;
        cin >> u >> v;
        u--;
        v--;
        bad.at(u).at(v) = true;
        bad.at(v).at(u) = true;
    }

    vector<int> p(N);
    iota(p.begin(), p.end(), 0);

    int ans = 1e9;
    do {
        bool ok = true;
        int current_time = 0;
        for (int i = 0; i < N; i++) {
            current_time += A.at(p.at(i)).at(i);
            if (i < N - 1 && bad.at(p.at(i)).at(p.at(i + 1))) {
                ok = false;
                break;
            }
        }
        if (ok) {
            ans = min(ans, current_time);
        }
    } while (next_permutation(p.begin(), p.end()));

    if (ans == 1e9) {
        cout << -1 << endl;
    } else {
        cout << ans << endl;
    }

    return 0;
}