#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    vector<int> P(n);
    for (int i = 0; i < n; i++) {
        cin >> P.at(i);
    }
    vector<int> Q(n);
    for (int i = 0; i < n; i++) {
        cin >> Q.at(i);
    }

    // いわゆる逆引きをする(indexのズレに注意)
    vector<int> rev_P(n);
    for (int i = 0; i < n; i++) {
        rev_P.at(P.at(i) - 1) = i;
    }
    vector<int> rev_Q(n);
    for (int i = 0; i < n; i++) {
        rev_Q.at(Q.at(i) - 1) = i;
    }
    // これで、逆関数を用いる場合は全部インデックス表記で
    // 元の関数はインデックス表記でないので注意

    for (int i = 0; i < n; i++) {
        cout << Q.at(P.at(rev_Q.at(i)) - 1) << " ";
    }
    cout << endl;

}