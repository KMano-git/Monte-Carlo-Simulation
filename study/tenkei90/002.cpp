#include <bits/stdc++.h>
using namespace std;

void make_pattern(int num_left, int num_right, string s, int N, vector<string>& reached_s) {
    // 終了条件：左右ともにN/2個使った
    if (num_left == N / 2 && num_right == N / 2) {
        reached_s.push_back(s);
        return;
    }
    
    // ( を追加できる条件：まだN/2個未満
    if (num_left < N / 2) {
        make_pattern(num_left + 1, num_right, s + '(', N, reached_s);
    }
    
    // ) を追加できる条件：左カッコより少ない（バランスを保つため）
    if (num_right < num_left) {
        make_pattern(num_left, num_right + 1, s + ')', N, reached_s);
    }
}

int main() {
    // Encyclopedia of Parentheses (★3)
    int N;
    cin >> N;
    // 奇数の場合は全て成り立たないので無視。
    if (N % 2) {
        cout << "" << endl;
        return 0;
    }

    // 順列問題 10で42
    // N / 2個の ( に ) を挿入する問題 → 樹形図みたいに発展していく
    vector<char> braket = { '(',')' };
    
    string s = ""; // ここに()をぶち込んでいく
    vector<string> reached_s;
    make_pattern(0, 0, s, N, reached_s);
    for (int i = 0; i < reached_s.size(); i++) {
        cout << reached_s.at(i) << endl;
    }
}

