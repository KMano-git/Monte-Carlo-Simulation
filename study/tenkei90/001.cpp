#include <bits/stdc++.h>
using namespace std;

// 最小ピースの長さが x 以上になるように K 回切れるか？
bool canCut(int x, vector<int> &A, int N, int L, int K) {
    int cuts = 0;      // 切った回数
    int last = 0;      // 最後に切った位置
    
    for (int i = 0; i < N; i++) {
        // 今の位置で切ると、前のピースの長さは A[i] - last
        if (A.at(i) - last >= x) {
            cuts++;
            last = A.at(i);
            if (cuts == K) break; // K回切ったらループを抜ける
            // K+1回切って条件を満たさないことがある
        }
    }
    // 最後のピース（last から L まで）も x 以上必要
    if (L - last < x) return false;
    
    return cuts >= K;
}

int main() {
    int N, L, K;
    cin >> N >> L >> K;
    vector<int> A(N);
    for (int i = 0; i < N; i++) {
        cin >> A.at(i);
    }

    // 二分探索 切れ目の長さを0~Lで二分探索
    int left = 0, right = L;
    while (right - left > 1) {
        int mid = (left + right) / 2;
        if (canCut(mid, A, N, L, K)) {
            left = mid;   // midでOK → もっと大きくできるかも
        } else {
            right = mid;  // midでNG → もっと小さくする
        }
    }
    
    cout << left << endl;
}