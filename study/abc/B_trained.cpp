#include <bits/stdc++.h>
using namespace std;

//isloopだと計算時間が多くなってしまう;判定法を変える
// bool isloop(vector<int> path, int i, int count) {
//     for (int j = 0; j < count; j++) {
//         if (path.at(j) == i) {
//             return true;
//         }
//     }
//     return false;
// }

int main() {
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a.at(i);
        a.at(i)--; //入力値を-1してインデックスにする
    }
    //数列の変化を追って、a.at(i)=1になるまでを求める
    //訪問履歴を記憶する（この場合の走査計算はO(n)）
    vector<bool> visited(n, false);
    int count = 0;
    int i = 0;
    bool end = false;
    while (end == false) {
        if (a.at(i) == 1) { //1に到達するボタンだったら、押して終わり
            count++;    
            end = true;
        } else if (visited.at(i) == true) { //訪問履歴にあったら、ループなので終わり
            count = -1;
            end = true;
        } else { //それ以外は、訪問履歴に記憶して、次のボタンを押す
            visited.at(i) = true;
            i = a.at(i);
            count++;
        }
    }
    cout << count << endl;
}