#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

#define ll long long
#define ld double
#define For(i, x, y) for(ll i = (x);i <= (y);++i)

inline ll read() //the quick read function
{
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

inline void write(ll x)  //the quick write function
{
    if (x < 0) putchar('-'), x = -x;
    if (x > 9) write(x / 10);
    putchar(x % 10 + 48);
}


//the constant to be used
const ll N = 20505;
const ld PI = acosl(-1);
const ll mo = (ll) 1e9 + 7;

//the structure for the complex
struct Complex {
    ld r, i;

    Complex() {}

    Complex(ld _r, ld _i) : r(_r), i(_i) {}

    Complex conj() const { return Complex(r, -i); }//�����

    //reload the calculation for the complex
    Complex operator+(const Complex &p) const { return Complex(r + p.r, i + p.i); }

    Complex operator-(const Complex &p) const { return Complex(r - p.r, i - p.i); }

    Complex operator*(const Complex &p) const { return Complex(r * p.r - i * p.i, r * p.i + i * p.r); }

    Complex operator/(const ld p) const { return Complex(r / p, i / p); }
};

inline ll add(ll p, ll q) { return (p += q) >= mo ? p - mo : p; }

inline ll sub(ll p, ll q) { return (p -= q) < 0 ? p + mo : p; }


/*class FFT
* the class for the fast fourier transformation(FFT)
* which can accelerate the convolution in the state transition in dp
*/
class FFT {
private:
    static const ll N = 131072;
    const ll P = 1e9 + 7;
    const ll base = 32768;
    Complex omega[N + 1], omegaInv[N + 1];
    Complex A[N + 1], B[N + 1];
    ll a0[N + 1], b0[N + 1], a1[N + 1], b1[N + 1];
    ll m;

    void init() {
        For(i, 0, N - 1) {
            omega[i] = Complex(cos(2 * PI / N * i), sin(2 * PI / N * i));
            omegaInv[i] = omega[i].conj();
        }
    }

    void reverse(Complex *a, ll n) {
        for (ll i = 0, j = 0; i < n; ++i) {
            if (i < j) swap(a[i], a[j]);
            for (ll l = n >> 1; (j ^= l) < l; l >>= 1) {}
        }
    }

    void transform(Complex *a, ll n, Complex *omega) {
        reverse(a, n);
        for (ll l = 2; l <= n; l <<= 1) {
            ll hl = l >> 1;
            for (Complex *x = a; x != a + n; x += l) {
                for (ll i = 0; i < hl; ++i) {
                    Complex t = omega[N / l * i] * x[i + hl];
                    x[i + hl] = x[i] - t;
                    x[i] = x[i] + t;
                }
            }
        }
    }

public:
    FFT() { init(); }

    int extend(int n) {
        int res = 1;
        while (res < n) res <<= 1;
        return m = res;
    }

    void dft(Complex *a, int n) { transform(a, n, omega); }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        For(i, 0, n - 1) a[i] = a[i] / n;
    }

    void mul(ll a[], ll b[], ll c[]) {
        For(i, 0, m - 1) A[i] = Complex(a[i], b[i]);
        dft(A, m);
        For(i, 0, m - 1) {
            ll j = (m - i) & (m - 1);
            B[i] = (A[i] * A[i] - (A[j] * A[j]).conj()) * Complex(0, -0.25);
        }
        idft(B, m);
        For(i, 0, m - 1) c[i] = (ll) (B[i].r + 0.5) % P;
    }

    void mulmod(ll a[], ll b[], ll c[]) {
        ll i;
        For(i, 0, m - 1) a0[i] = a[i] >> 15, b0[i] = b[i] >> 15;
        for (mul(a0, b0, a0), i = 0; i < m; ++i) {
            c[i] = 1ll * a0[i] * base * base % P;
            a1[i] = a[i] & (base - 1), b1[i] = b[i] & (base - 1);
        }
        for (mul(a1, b1, a1), i = 0; i < m; ++i) {
            c[i] = add(a1[i], c[i]), a0[i] = add(a0[i], a1[i]);
            a1[i] = (a[i] >> 15) + (a[i] & (base - 1)), b1[i] = (b[i] >> 15) + (b[i] & base - 1);
        }
        for (mul(a1, b1, a1), i = 0; i < m; ++i)
            c[i] = (1ll * base * sub(a1[i], a0[i]) + c[i]) % P;
    }
} fft;

ll n, m;

ll dp[N][20][2];//the 3d array for the dp
/*
* the dp[i][j][0] represents the number of red&black trees have i vetices, have the black height j, and its root vertex is red
* the dp[i][j][1] represents the number of red&black trees have i vetices, have the black height j, and its root vertex is black
*/


ll Log[N], aa[N << 1], bb[N << 1], cc[N << 1];

int main() {
    n = read();

    //the initial condition
    dp[0][0][1] = 1;//the null vertex has only 1 type
    dp[1][0][0] = 1;//the single red vertex has only 1 type

    clock_t start, finish;
    start = clock();
    For(loop, 0, 500) {
        For(i, 2, n + 3) Log[i] = Log[i >> 1] + 1;
        //to pre-calculate a table of log, so we needn't calc it in the loop and can decrease the time complexity

        For(j, 1, Log[n + 1])
            //our loop for dp, j is the black height, as a red&black tree with size n have the maximum black height of log(n+1)
            //we only enumerate it to log(n+1)
        {
            For(i, 0, n - 1) aa[i] = bb[i] = add(dp[i][j - 1][0], dp[i][j - 1][1]);
            //the FFT accelerate part:
            For(i, n, n * 2 - 1) aa[i] = bb[i] = 0;
            m = fft.extend(n << 1);
            fft.mulmod(aa, bb, cc);

            For(i, 1, n) dp[i][j][1] = cc[i - 1];//,cerr<<cc[i]<<' ';cerr<<endl;
            For(i, 0, n - 1) aa[i] = bb[i] = dp[i][j][1];
            For(i, n, n * 2 - 1) aa[i] = bb[i] = 0;
            m = fft.extend(n << 1);
            fft.mulmod(aa, bb, cc);
            For(i, 1, n) dp[i][j][0] = cc[i - 1];
        }
    }
    finish = clock();
    ll ans = 0;
    For(i, 0, n) ans = add(ans, dp[n][i][1]);//,cerr<<dp[n][i][1]<<endl;
    cout << ans % 1000000007 << endl;
    cout << finish - start << endl;
}