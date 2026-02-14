#include <bits/stdc++.h>
using namespace std;
const double PI = acos(-1);

void calc_angle(double t, int T, int L, int X, int Y) {
    double theta = 2 * PI * t / T;
    double dx = (double)X;
    double dy = (double)Y + (double)L / 2 * sin(theta);
    double dz = (double)L / 2 - (double)L / 2 * cos(theta);
    double angle = atan2(dz, sqrt(dx * dx + dy * dy));
    double angle_deg = angle * 180 / PI;
    cout << fixed << setprecision(10) << angle_deg << "\n";
    return;
}

int main() {
    int T;
    cin >> T;
    int L, X, Y;
    cin >> L >> X >> Y;

    // T(min)で一周する観覧車から座標(X,Y,0)に存在する像を見る角度を求める
    // 観覧車の中心を(0,0,L/2)とする

    int Q;
    cin >> Q;
    for (int i = 0; i < Q; i++) {
        double t;
        cin >> t;
        calc_angle(t, T, L, X, Y);
    }
}