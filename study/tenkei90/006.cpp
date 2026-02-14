#include <bits/stdc++.h>
using namespace std;

void make_substring(string S, vector<string> s, int N, int K) {
    // ベースケース

    // 終了条件
    int cnt_notend = 0;
    for (int i = 0; i < K; i++) {
        if (s.at(i).size() != N) {
            cnt_notend++;
        }
    }   
    if (cnt_notend == 0) {
        return;
    }
    // 再帰
    for (int i = 0; i < K; i++) {
        if (s.at(i).size() != N) {
            make_substring(S, s, N, K);
        }
    }
}

int main() {
    int N, K;
    cin >> N >> K;
    string S;
    cin >> S;

    // 文字列Sの部分列sを生成する
    vector<string> s;
    make_substring(S, s, N, K);
    // その中でもっとも辞書順で小さいものを求める
    string min_s = s.at(0);
    for (int i = 1; i < K; i++) {
        if (s.at(i) < min_s) {
            min_s = s.at(i);
        }
    }
    cout << min_s << endl;
}
// 貪欲法(こっちが正解)
// #include <bits/stdc++.h>
// using namespace std;

// int main() {
//     int N, K;
//     cin >> N >> K;
//     string S;
//     cin >> S;
    
//     string ans = "";
//     int pos = 0;  // 現在の探索開始位置
    
//     for (int i = 0; i < K; i++) {
//         // あと (K - i) 文字選ぶ必要がある
//         // なので、S[pos] ~ S[N - (K - i)] の範囲から最小を選ぶ
//         int last = N - (K - i);
//         char min_char = S[pos];
//         int min_pos = pos;
        
//         for (int j = pos; j <= last; j++) {
//             if (S[j] < min_char) {
//                 min_char = S[j];
//                 min_pos = j;
//             }
//         }
        
//         ans += min_char;
//         pos = min_pos + 1;  // 選んだ位置の次から探す
//     }
    
//     cout << ans << endl;
// }