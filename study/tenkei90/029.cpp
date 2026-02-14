#include <bits/stdc++.h>
#include <atcoder/lazysegtree>
using namespace std;
using namespace atcoder;

// S: 扱うデータ型 (区間の最大値を持つため int)
using S = int;
// F: 遅延させるデータ型 (更新する値を持つため int)
using F = int;

// op: 区間の取得クエリで返す値 (最大値を求めたいので max)
S op(S a, S b) { return std::max(a, b); }

// e: opの単位元 (maxなので、影響しない最小の値 0)
S e() { return 0; }

// mapping: 遅延データを値に反映させる (fで上書き)
// fが0(ID)なら何もしない
S mapping(F f, S x) { return (f == 0 ? x : f); }

// composition: 遅延データを合成する (後から来たfで上書き)
// fが0(ID)ならgのまま
F composition(F f, F g) { return (f == 0 ? g : f); }

// id: 遅延データの単位元 (何もしないを表す 0)
F id() { return 0; }

int main() {
    int W, N;
    cin >> W >> N;
    
    // やりたいこと（愚直にやると）
    // Li~Riの範囲の値の最大値を取得
    // その範囲の値を取得した最大値+1に更新する
    // これをN回繰り返す
    // O(N*W)
    // vector<int> h(W,0);
    // for (int i = 0; i < N; i++) {
    //     int L, R;
    //     cin >> L >> R;
    //     L--, R--;
    //     int max_height = 0;
    //     for (int j = L; j <= R; j++) {
    //         max_height = max(max_height, h.at(j));
    //     }
    //     for (int j = L; j <= R; j++) {
    //         h.at(j) = max_height + 1;
    //     }
    //     cout << h.at(L) << "\n";
    // }

    
    // 効率化
    // 遅延評価セグメント木を使う

    // セグメント木の宣言 (サイズ W, 初期値 e()=0)
    lazy_segtree<S, op, e, F, mapping, composition, id> seg(W);

    for (int i = 0; i < N; ++i) {
        int L, R;
        cin >> L >> R;
        // 0-indexed に変換
        L--;
        // 区間 [L, R) の最大値を取得 (Rは半開区間なのでそのままRでOK)
        // 入力のRは閉区間なので、半開区間 [L, R+1) にするなら seg.prod(L, R) ではなく seg.prod(L, R) -> seg.prod(L-1, R)
        // 問題文: [Li, Ri] (1-indexed 閉区間)
        // ACL: [l, r) (0-indexed 半開区間)
        // なので、L -> L-1, R -> R (そのまま)
        // 例: 1〜3 -> 0〜2 (index 0, 1, 2) -> 半開区間 [0, 3) 
        // つまり L--, R はそのままでOK
        
        int height = seg.prod(L, R);
        
        // 区間 [L, R) を height + 1 で更新
        seg.apply(L, R, height + 1);
        
        cout << height + 1 << endl;
    }
}