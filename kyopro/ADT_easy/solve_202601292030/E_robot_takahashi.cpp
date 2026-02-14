#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    string s;
    cin >> s;
    vector<int> W(n);
    for (int i = 0; i < n; i++) {
        cin >> W.at(i);
    }

    // sortして、比較を効率化(どこに線を引くかで1,0が決定する)
    vector<pair<int, int>> p(n);
    for (int i = 0; i < n; i++) {
        p.at(i) = {W.at(i), s.at(i) - '0'};
    }
    sort(p.begin(), p.end());

    // 境界値を1つずつずらす際に、前回の結果を利用して計算量を削減
    // 初期値は全て大人判定
    int count = 0;
    int count_child = 0;
    int max_count = 0;
    for (int i = 0; i < n; i++) {
        if (p.at(i).second == 1) {
            count++;
        } else {
            count_child++;
        }
    }
    max_count = max(count, count_child);
    
    // 境界値を1つずつずらす、ずらす値によってcountを変更
    // 同じ値は一度で飛び越えてもらわないといけないので、while処理を挟む
    // 全部子供の場合の処理を考えるのが面倒なので、例外処理で対応。
    // 0 を飛び越えたらcountは+1, 1を飛び越えたらcountは-1
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && p.at(j).first == p.at(i).first) {
            if (p.at(j).second == 1) {
                count -= 1;
            } else {
                count += 1;
            }
            j++;
        }
        i = j; // カウンタの取り扱いに注意 whileの条件式でjが増加しているのでiに代入するだけで良い
        max_count = max(max_count, count);
    }
    cout << max_count << endl;
}

// あるいはこういうふうに書くことでスキップできる
// for (int i = 0; i < n; i++) {
// 		if (a[i].second == '1')x--;
// 		else x++;
// 		if (i < (n - 1)) {
// 			if (a[i].first != a[i + 1].first)ans = max(ans, x); // ここでスキップ
// 		}
// 		else ans = max(ans, x);
// 	}