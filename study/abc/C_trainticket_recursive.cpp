#include <bits/stdc++.h>
using namespace std;

// 再帰関数で演算子の組み合わせを試す
// pos: 現在の位置（0〜3）, sum: 現在の合計, ops: 演算子の履歴, digits: 各桁の数字
bool solve(int pos, int sum, string ops, vector<int> &digits) {
    // ベースケース：4桁すべて処理した
    if (pos == 4) {
        if (sum == 7) {
            // 答えを出力
            cout << digits[0] << ops[0] << digits[1] << ops[1] << digits[2] << ops[2] << digits[3] << "=7" << endl;
            return true;
        }
        return false;
    }
    
    // 最初の桁は足すだけ
    if (pos == 0) {
        return solve(pos + 1, digits[0], ops, digits);
    }
    
    // '+' を試す
    if (solve(pos + 1, sum + digits[pos], ops + '+', digits)) {
        return true;
    }
    
    // '-' を試す
    if (solve(pos + 1, sum - digits[pos], ops + '-', digits)) {
        return true;
    }
    
    return false;
}

int main() {
    int in;
    cin >> in;
    
    vector<int> digits(4);
    digits[0] = in / 1000;
    digits[1] = (in / 100) % 10;
    digits[2] = (in / 10) % 10;
    digits[3] = in % 10;
    
    string ops = "";
    solve(0, 0, ops, digits);
}
