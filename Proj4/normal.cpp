#pragma GCC optimize("Ofast")

#include <iostream>
#include <ctime>

using namespace std;
#define ll long long
#define For(i, x, y) for(ll i = (x);i <= (y);++i)

inline ll read() {
    ll x = 0, f = 1;
    char ch = getchar();
    while (ch < '0' || ch > '9') {
        if (ch == '-')f = -1;
        ch = getchar();
    }
    while (ch >= '0' && ch <= '9') {
        x = x * 10 + ch - '0';
        ch = getchar();
    }
    return x * f;
}

inline void write(ll x) {
    if (x < 0) putchar('-'), x = -x;
    if (x > 9) write(x / 10);
    putchar(x % 10 + 48);
}

const ll mo = (ll) 1000000007;
const ll N = 505;

inline ll add(ll p, ll q) {
    return (p += q) >= mo ? p - mo : p;
}

ll dp[N][N][2], Log[N];

int main() {
    ll n = read();
    For(i, 2, n + 3) Log[i] = Log[i >> 1] + 1;
    dp[0][0][1] = 1;
    dp[1][0][0] = 1;
    dp[1][1][1] = 1;

    clock_t start, finish;
    start = clock();
    For(loop, 0, 500)
        For(i, 2, n)
            For(j, 1, Log[i + 3]) {
                For(k, 0, i - 1) {
                    dp[i][j][1] = add(dp[i][j][1], 1ll * (dp[k][j - 1][0] + dp[k][j - 1][1]) *
                                                   (dp[i - k - 1][j - 1][0] + dp[i - k - 1][j - 1][1]) % mo);
                    dp[i][j][0] = add(dp[i][j][0], 1ll * dp[k][j][1] * dp[i - k - 1][j][1] % mo);
                }
            }
    finish = clock();
    ll ans = 0;
    For(i, 0, n) ans = add(ans, dp[n][i][1]);//,cerr<<dp[n][i][1]<<endl;
    cout << ans % 1000000007 << endl;
    cout << finish - start << endl;
}