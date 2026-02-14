#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

// B x B の行列を表す型
using Matrix = vector<vector<ll>>;

// 行列の積を計算
Matrix multiply(const Matrix &A, const Matrix &B, int size) {
    Matrix C(size, vector<ll>(size, 0));
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            if (A[i][k] == 0) continue;  // 枝刈り
            for (int j = 0; j < size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % MOD;
            }
        }
    }
    return C;
}

// 行列のN乗を計算（ダブリング）
Matrix matrix_power(Matrix A, ll N, int size) {
    // 単位行列で初期化
    Matrix result(size, vector<ll>(size, 0));
    for (int i = 0; i < size; i++) {
        result[i][i] = 1;
    }
    
    while (N > 0) {
        if (N & 1) {
            result = multiply(result, A, size);
        }
        A = multiply(A, A, size);
        N >>= 1;
    }
    return result;
}

int main() {
    ll N, B, K;
    cin >> N >> B >> K;
    vector<int> c(K);
    for (int i = 0; i < K; i++) {
        cin >> c[i];
    }
    
    // 遷移行列を作成
    // A[i][j] = 「余りがjの状態から数字を1つ追加して余りがiになる」方法の数
    // 余りがjのとき、数字dを追加すると余りは (j*10 + d) % B になる
    Matrix A(B, vector<ll>(B, 0));
    for (int j = 0; j < B; j++) {
        for (int k = 0; k < K; k++) {
            int next = (j * 10 + c[k]) % B;
            A[next][j]++;
        }
    }
    
    // 行列をN乗
    Matrix result = matrix_power(A, N, B);
    
    // 初期状態: 余り0から始まる（0桁の数は余り0）
    // 答え: N桁後に余り0になる場合の数 = result[0][0]
    cout << result[0][0] << endl;
    
    return 0;
}
