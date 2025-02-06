#include <algorithm>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>#include <algorithm>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<ll, ll> P;

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0以上1未満の小数をとる乱数
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = randxor() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

bool IsNG(int x, int y, int n)
{
  if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

double TL = 1.8;
int mode;
clock_t startTime, endTime;

int N, M, T;
int NN = 2 * 6 * (6 - 1);
vector<vector<int>> X(NN, vector<int>(15, 0));
int maxX[20];
int XX[100][20];
int order[1000];

vector<vector<int>> A(6, vector<int>(6, 0));

int uuu[10][6][6][15], vvv[10][6][6][15];

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp() {}

// 入力受け取り
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> N >> M >> T;
    for (int i = 0; i < NN; i++) {
      for (int j = 0; j < M; j++) {
        cin >> X[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    ifs >> N >> M >> T;
    for (int i = 0; i < NN; i++) {
      for (int j = 0; j < M; j++) {
        ifs >> X[i][j];
      }
    }
    rep(i, T)
    {
      rep(j, N)
      {
        rep(k, N - 1)
        {
          string sss;
          ifs >> sss;
          rep(l, M) { uuu[i][j][k][l] = sss[l] - '0'; }
        }
      }
      rep(j, N - 1)
      {
        rep(k, N)
        {
          string sss;
          ifs >> sss;
          rep(l, M) { vvv[i][j][k][l] = sss[l] - '0'; }
        }
      }
    }
  }

  rep(j, M) maxX[j] = 0;
  rep(i, NN)
  {
    rep(j, M) { maxX[j] = max(maxX[j], X[i][j]); }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
ll CalcScore()
{
  ll sumX = 0;
  rep(j, M) sumX += maxX[j];

  int ma = 0;
  rep(i, NN)
  {
    int tmp = 0;
    rep(j, M) { tmp += X[i][j]; }
    ma = max(ma, tmp);
  }

  ll res = round(1000000.0 * ma / sumX);
  return res;
}

// 初期解生成
void Initialize() {}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
  }
  else {
  }
}

void GetReturn(int turn)
{
  int now = 0;
  rep(i, N)
  {
    rep(j, N - 1)
    {
      rep(k, M)
      {
        if (uuu[turn % T][i][j][k] == 0) {
          X[now][k] = XX[A[i][j]][k];
        }
        else {
          X[now][k] = XX[A[i][j + 1]][k];
        }
      }
      now++;
    }
  }
  rep(i, N - 1)
  {
    rep(j, N)
    {
      rep(k, M)
      {
        if (vvv[turn % T][i][j][k] == 0) {
          X[now][k] = XX[A[i][j]][k];
        }
        else {
          X[now][k] = XX[A[i + 1][j]][k];
        }
      }
      now++;
    }
  }
}

void Convert()
{
  vector<P> v[M];
  rep(i, NN)
  {
    rep(j, M) { v[j].push_back(P(X[i][j], i)); }
  }
  rep(i, M) { sort(v[i].begin(), v[i].end(), greater<P>()); }
  int val[NN][M];
  rep(i, NN)
  {
    rep(j, M) { val[v[j][i].second][j] = i; }
  }

  int ma[M] = {};
  rep(i, NN)
  {
    rep(j, M) { ma[j] = max(ma[j], X[i][j]); }
  }

  rep(i, NN)
  {
    rep(j, M)
    {
      XX[i][j] = X[i][j];
      // if (X[i][j] == ma[j]) X[i][j] * 5;
      X[i][j] = X[i][j] * X[i][j] * X[i][j];
    }
  }
}

void Convert2()
{
  vector<P> v[M];
  rep(i, NN)
  {
    rep(j, M) { v[j].push_back(P(X[i][j], i)); }
  }
  rep(i, M) { sort(v[i].begin(), v[i].end(), greater<P>()); }
  int val[NN][M];
  rep(i, NN)
  {
    rep(j, M)
    {
      // cout << i << ' ' << j << ' ' << v[j][i].second << endl;
      val[v[j][i].second][j] = NN - i;
    }
  }

  rep(i, NN)
  {
    rep(j, M)
    {
      XX[i][j] = X[i][j];
      X[i][j] = val[i][j] * val[i][j] * val[i][j];
    }
  }
}

void Sort()
{
  vector<P> vp;
  rep(i, NN)
  {
    int sum = 0;
    rep(j, M) sum += X[i][j];
    vp.push_back(P(sum, i));
  }
  sort(vp.begin(), vp.end(), greater<P>());
  rep(i, vp.size()) { order[i] = vp[i].second; }
}

int o[6][6] = {
    {21, 20, 19, 18, 17, 16}, {22, 7, 6, 5, 16, 35},
    {23, 8, 1, 4, 15, 34},    {24, 9, 2, 3, 14, 33},
    {25, 10, 11, 12, 13, 32}, {26, 27, 28, 29, 30, 31},
};

vector<P> guruguruOrder;
void InitGuruGuru()
{
  guruguruOrder.clear();
  vector<pair<int, P>> vp;
  rep(i, 6)
  {
    rep(j, 6) { vp.push_back(make_pair(o[i][j], P(i, j))); }
  }
  sort(vp.begin(), vp.end());
  for (auto p : vp) {
    guruguruOrder.push_back(p.second);
  }
}

void GuruGuru(int t)
{
  rep(i, 36)
  {
    int x = guruguruOrder[i].first;
    int y = guruguruOrder[i].second;
    A[x][y] = order[i];
  }

  vector<P> vp1, vp2;
  rep(i, NN)
  {
    ll sum1 = 0, sum2 = 0;
    rep(j, M)
    {
      if (j % 2 == 0) {
        sum1 += X[i][j];
      }
      else {
        sum2 += X[i][j];
      }
    }
    vp1.push_back(P(sum1, i));
    vp2.push_back(P(sum2, i));
  }

  sort(vp1.begin(), vp1.end(), greater<P>());
  sort(vp2.begin(), vp2.end(), greater<P>());
  int use[NN] = {};
  int itr1 = 0, itr2 = 0;
  rep(i, N)
  {
    rep(j, N)
    {
      if ((i + j) % 2 == 0) {
        while (true) {
          int id = vp1[itr1].second;
          if (use[id]) {
            itr1++;
            continue;
          }
          else {
            A[i][j] = id;
            use[id] = 1;
            break;
          }
        }
      }
      else {
        while (true) {
          int id = vp2[itr2].second;
          if (use[id]) {
            itr2++;
            continue;
          }
          else {
            A[i][j] = id;
            use[id] = 1;
            break;
          }
        }
      }
    }
  }

  if (t < 100) {
    rep(_, 30000)
    {
      int ra1 = randxor() % 6;
      int ra2 = randxor() % 6;
      int ra3 = randxor() % 6;
      int ra4 = randxor() % 6;

      int before = 0;
      {
        int x = ra1;
        int y = ra2;
        int id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { before += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
        x = ra3;
        y = ra4;
        id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { before += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
      }

      swap(A[ra1][ra2], A[ra3][ra4]);

      int after = 0;
      {
        int x = ra1;
        int y = ra2;
        int id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { after += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
        x = ra3;
        y = ra4;
        id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { after += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
      }

      if (before > after) {
        swap(A[ra1][ra2], A[ra3][ra4]);
      }
    }
  }
}

ll InnerGuruGuru(int x, int y)
{
  ll res = 0;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];
    rep(j, M) { res += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
    // rep(j, M) {
    //   res +=
    //       X[id][j] * 100 + max(X[id][j], X[nd][j]) + abs(X[id][j] -
    //       X[nd][j]);
    // }
  }
  return res;
}

vector<ll> InnerGuruGuru3(int x, int y)
{
  vector<ll> res;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];

    ll tmp = 0;
    rep(j, M) { tmp += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
  }
  return res;
}

ll InnerGuruGuru2(int x, int y)
{
  ll res = 0;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];
    ll ma = 0;
    rep(loop, 5)
    {
      ll tmp = 0;
      rep(j, M)
      {
        if (randxor() % 2) {
          tmp += X[id][j];
        }
        else {
          tmp += X[nd][j];
        }
      }
      ma = max(ma, tmp);
    }
    res = max(res, ma);
  }
  return res;
}

void GuruGuru2(int t)
{
  rep(i, 36)
  {
    int x = guruguruOrder[i].first;
    int y = guruguruOrder[i].second;
    A[x][y] = order[i];
  }

  // vector<P> vp1, vp2;
  // rep(i, NN) {
  //   ll sum1 = 0, sum2 = 0;
  //   rep(j, M) {
  //     if (j % 2 == 0) {
  //       sum1 += X[i][j];
  //     } else {
  //       sum2 += X[i][j];
  //     }
  //   }
  //   vp1.push_back(P(sum1, i));
  //   vp2.push_back(P(sum2, i));
  // }

  // sort(vp1.begin(), vp1.end(), greater<P>());
  // sort(vp2.begin(), vp2.end(), greater<P>());
  // int use[NN] = {};
  // int itr1 = 0, itr2 = 0;
  // rep(i, N) {
  //   rep(j, N) {
  //     if ((i + j) % 2 == 0) {
  //       while (true) {
  //         int id = vp1[itr1].second;
  //         if (use[id]) {
  //           itr1++;
  //           continue;
  //         } else {
  //           A[i][j] = id;
  //           use[id] = 1;
  //           break;
  //         }
  //       }
  //     } else {
  //       while (true) {
  //         int id = vp2[itr2].second;
  //         if (use[id]) {
  //           itr2++;
  //           continue;
  //         } else {
  //           A[i][j] = id;
  //           use[id] = 1;
  //           break;
  //         }
  //       }
  //     }
  //   }
  // }

  if (t < 100) {
    rep(_, 50000)
    {
      int ra1 = randxor() % 6;
      int ra2 = randxor() % 6;
      int ra3 = randxor() % 6;
      int ra4 = randxor() % 6;
      if (t < 109) {
      }
      else {
        ra1 = randxor() % 4 + 1;
        ra2 = randxor() % 4 + 1;
        ra3 = randxor() % 4 + 1;
        ra4 = randxor() % 4 + 1;
      }

      // if ((ra1 + ra2) % 2 != (ra3 + ra4) % 2) continue;

      int before = 0;
      before += InnerGuruGuru(ra1, ra2);
      before += InnerGuruGuru(ra3, ra4);

      swap(A[ra1][ra2], A[ra3][ra4]);

      int after = 0;
      after += InnerGuruGuru(ra1, ra2);
      after += InnerGuruGuru(ra3, ra4);

      double temp = (50000 - _) / 50000.0 * 10;
      const double prob = exp((double)(after - before) / temp);

      if (prob > rand01()) {
      }
      else {
        swap(A[ra1][ra2], A[ra3][ra4]);
      }
    }
  }
}

ll Solve(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  InitGuruGuru();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  for (int t = 0; t < T; t++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        A[i][j] = i * N + j;
      }
    }

    Convert();
    Sort();
    GuruGuru2(t);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (mode == 0) {
          cout << A[i][j];

          if (j < N - 1) {
            cout << " ";
          }
          else {
            cout << endl;
          }
        }
        else {
          ofs << A[i][j];

          if (j < N - 1) {
            ofs << " ";
          }
          else {
            ofs << endl;
          }
        }
      }
    }

    if (mode == 0) {
      cout.flush();

      for (int i = 0; i < NN; i++) {
        for (int j = 0; j < M; j++) {
          cin >> X[i][j];
        }
      }
    }
    else {
      GetReturn(t);
    }
  }

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }

  return score;

  // 初期解生成
  Initialize();

  // 解答を出力
  Output(ofs);

  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}

#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<ll, ll> P;

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0以上1未満の小数をとる乱数
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = randxor() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

bool IsNG(int x, int y, int n)
{
  if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

double TL = 1.8;
int mode;
clock_t startTime, endTime;

const int N = 6;
const int M = 15;
const int T = 10;
const int NN = 2 * 6 * (6 - 1);
vector<vector<int>> X(NN, vector<int>(15, 0));
int maxX[20];
int XX[100][20];
int order[1000];

vector<vector<int>> A(6, vector<int>(6, 0));

int uuu[10][6][6][15], vvv[10][6][6][15];

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp() {}

// 入力受け取り
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    int n, m, t;
    cin >> n >> m >> t;
    for (int i = 0; i < NN; i++) {
      for (int j = 0; j < M; j++) {
        cin >> X[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    int n, m, t;
    ifs >> n >> m >> t;
    for (int i = 0; i < NN; i++) {
      for (int j = 0; j < M; j++) {
        ifs >> X[i][j];
      }
    }
    rep(i, T)
    {
      rep(j, N)
      {
        rep(k, N - 1)
        {
          string sss;
          ifs >> sss;
          rep(l, M) { uuu[i][j][k][l] = sss[l] - '0'; }
        }
      }
      rep(j, N - 1)
      {
        rep(k, N)
        {
          string sss;
          ifs >> sss;
          rep(l, M) { vvv[i][j][k][l] = sss[l] - '0'; }
        }
      }
    }
  }

  rep(j, M) maxX[j] = 0;
  rep(i, NN)
  {
    rep(j, M) { maxX[j] = max(maxX[j], X[i][j]); }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
ll CalcScore()
{
  ll sumX = 0;
  rep(j, M) sumX += maxX[j];

  int ma = 0;
  rep(i, NN)
  {
    int tmp = 0;
    rep(j, M) { tmp += X[i][j]; }
    ma = max(ma, tmp);
  }

  ll res = round(1000000.0 * ma / sumX);
  return res;
}

// 初期解生成
void Initialize() {}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
  }
  else {
  }
}

void GetReturn(int turn)
{
  int now = 0;
  rep(i, N)
  {
    rep(j, N - 1)
    {
      rep(k, M)
      {
        if (uuu[turn % T][i][j][k] == 0) {
          X[now][k] = XX[A[i][j]][k];
        }
        else {
          X[now][k] = XX[A[i][j + 1]][k];
        }
      }
      now++;
    }
  }
  rep(i, N - 1)
  {
    rep(j, N)
    {
      rep(k, M)
      {
        if (vvv[turn % T][i][j][k] == 0) {
          X[now][k] = XX[A[i][j]][k];
        }
        else {
          X[now][k] = XX[A[i + 1][j]][k];
        }
      }
      now++;
    }
  }
}

void Convert()
{
  vector<P> v[M];
  rep(i, NN)
  {
    rep(j, M) { v[j].push_back(P(X[i][j], i)); }
  }
  rep(i, M) { sort(v[i].begin(), v[i].end(), greater<P>()); }
  int val[NN][M];
  rep(i, NN)
  {
    rep(j, M) { val[v[j][i].second][j] = i; }
  }

  int ma[M] = {};
  rep(i, NN)
  {
    rep(j, M) { ma[j] = max(ma[j], X[i][j]); }
  }

  rep(i, NN)
  {
    rep(j, M)
    {
      XX[i][j] = X[i][j];
      // if (X[i][j] == ma[j]) X[i][j] * 5;
      X[i][j] = X[i][j] * X[i][j] * X[i][j];
    }
  }
}

void Convert2()
{
  vector<P> v[M];
  rep(i, NN)
  {
    rep(j, M) { v[j].push_back(P(X[i][j], i)); }
  }
  rep(i, M) { sort(v[i].begin(), v[i].end(), greater<P>()); }
  int val[NN][M];
  rep(i, NN)
  {
    rep(j, M)
    {
      // cout << i << ' ' << j << ' ' << v[j][i].second << endl;
      val[v[j][i].second][j] = NN - i;
    }
  }

  rep(i, NN)
  {
    rep(j, M)
    {
      XX[i][j] = X[i][j];
      X[i][j] = val[i][j] * val[i][j] * val[i][j];
    }
  }
}

void Sort()
{
  vector<P> vp;
  rep(i, NN)
  {
    int sum = 0;
    rep(j, M) sum += X[i][j];
    vp.push_back(P(sum, i));
  }
  sort(vp.begin(), vp.end(), greater<P>());
  rep(i, vp.size()) { order[i] = vp[i].second; }
}

int o[6][6] = {
    {21, 20, 19, 18, 17, 16}, {22, 7, 6, 5, 16, 35},
    {23, 8, 1, 4, 15, 34},    {24, 9, 2, 3, 14, 33},
    {25, 10, 11, 12, 13, 32}, {26, 27, 28, 29, 30, 31},
};

vector<P> guruguruOrder;
void InitGuruGuru()
{
  guruguruOrder.clear();
  vector<pair<int, P>> vp;
  rep(i, 6)
  {
    rep(j, 6) { vp.push_back(make_pair(o[i][j], P(i, j))); }
  }
  sort(vp.begin(), vp.end());
  for (auto p : vp) {
    guruguruOrder.push_back(p.second);
  }
}

void GuruGuru(int t)
{
  rep(i, 36)
  {
    int x = guruguruOrder[i].first;
    int y = guruguruOrder[i].second;
    A[x][y] = order[i];
  }

  vector<P> vp1, vp2;
  rep(i, NN)
  {
    ll sum1 = 0, sum2 = 0;
    rep(j, M)
    {
      if (j % 2 == 0) {
        sum1 += X[i][j];
      }
      else {
        sum2 += X[i][j];
      }
    }
    vp1.push_back(P(sum1, i));
    vp2.push_back(P(sum2, i));
  }

  sort(vp1.begin(), vp1.end(), greater<P>());
  sort(vp2.begin(), vp2.end(), greater<P>());
  int use[NN] = {};
  int itr1 = 0, itr2 = 0;
  rep(i, N)
  {
    rep(j, N)
    {
      if ((i + j) % 2 == 0) {
        while (true) {
          int id = vp1[itr1].second;
          if (use[id]) {
            itr1++;
            continue;
          }
          else {
            A[i][j] = id;
            use[id] = 1;
            break;
          }
        }
      }
      else {
        while (true) {
          int id = vp2[itr2].second;
          if (use[id]) {
            itr2++;
            continue;
          }
          else {
            A[i][j] = id;
            use[id] = 1;
            break;
          }
        }
      }
    }
  }

  if (t < 100) {
    rep(_, 30000)
    {
      int ra1 = randxor() % 6;
      int ra2 = randxor() % 6;
      int ra3 = randxor() % 6;
      int ra4 = randxor() % 6;

      int before = 0;
      {
        int x = ra1;
        int y = ra2;
        int id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { before += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
        x = ra3;
        y = ra4;
        id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { before += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
      }

      swap(A[ra1][ra2], A[ra3][ra4]);

      int after = 0;
      {
        int x = ra1;
        int y = ra2;
        int id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { after += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
        x = ra3;
        y = ra4;
        id = A[x][y];
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny, N)) continue;
          int nd = A[nx][ny];
          rep(j, M) { after += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
        }
      }

      if (before > after) {
        swap(A[ra1][ra2], A[ra3][ra4]);
      }
    }
  }
}

ll InnerGuruGuru(int x, int y)
{
  ll res = 0;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];
    rep(j, M) { res += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
    // rep(j, M) {
    //   res +=
    //       X[id][j] * 100 + max(X[id][j], X[nd][j]) + abs(X[id][j] -
    //       X[nd][j]);
    // }
  }
  return res;
}

vector<ll> InnerGuruGuru3(int x, int y)
{
  vector<ll> res;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];

    ll tmp = 0;
    rep(j, M) { tmp += 10 * X[id][j] + 2 * abs(X[id][j] - X[nd][j]); }
  }
  return res;
}

ll InnerGuruGuru2(int x, int y)
{
  ll res = 0;
  int id = A[x][y];
  rep(i, 4)
  {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (IsNG(nx, ny, N)) continue;
    int nd = A[nx][ny];
    ll ma = 0;
    rep(loop, 5)
    {
      ll tmp = 0;
      rep(j, M)
      {
        if (randxor() % 2) {
          tmp += X[id][j];
        }
        else {
          tmp += X[nd][j];
        }
      }
      ma = max(ma, tmp);
    }
    res = max(res, ma);
  }
  return res;
}

void GuruGuru2(int t)
{
  rep(i, 36)
  {
    int x = guruguruOrder[i].first;
    int y = guruguruOrder[i].second;
    A[x][y] = order[i];
  }

  // vector<P> vp1, vp2;
  // rep(i, NN) {
  //   ll sum1 = 0, sum2 = 0;
  //   rep(j, M) {
  //     if (j % 2 == 0) {
  //       sum1 += X[i][j];
  //     } else {
  //       sum2 += X[i][j];
  //     }
  //   }
  //   vp1.push_back(P(sum1, i));
  //   vp2.push_back(P(sum2, i));
  // }

  // sort(vp1.begin(), vp1.end(), greater<P>());
  // sort(vp2.begin(), vp2.end(), greater<P>());
  // int use[NN] = {};
  // int itr1 = 0, itr2 = 0;
  // rep(i, N) {
  //   rep(j, N) {
  //     if ((i + j) % 2 == 0) {
  //       while (true) {
  //         int id = vp1[itr1].second;
  //         if (use[id]) {
  //           itr1++;
  //           continue;
  //         } else {
  //           A[i][j] = id;
  //           use[id] = 1;
  //           break;
  //         }
  //       }
  //     } else {
  //       while (true) {
  //         int id = vp2[itr2].second;
  //         if (use[id]) {
  //           itr2++;
  //           continue;
  //         } else {
  //           A[i][j] = id;
  //           use[id] = 1;
  //           break;
  //         }
  //       }
  //     }
  //   }
  // }

  if (t < 100) {
    rep(_, 50000)
    {
      int ra1 = randxor() % 6;
      int ra2 = randxor() % 6;
      int ra3 = randxor() % 6;
      int ra4 = randxor() % 6;
      if (t < 109) {
      }
      else {
        ra1 = randxor() % 4 + 1;
        ra2 = randxor() % 4 + 1;
        ra3 = randxor() % 4 + 1;
        ra4 = randxor() % 4 + 1;
      }

      // if ((ra1 + ra2) % 2 != (ra3 + ra4) % 2) continue;

      int before = 0;
      before += InnerGuruGuru(ra1, ra2);
      before += InnerGuruGuru(ra3, ra4);

      swap(A[ra1][ra2], A[ra3][ra4]);

      int after = 0;
      after += InnerGuruGuru(ra1, ra2);
      after += InnerGuruGuru(ra3, ra4);

      double temp = (50000 - _) / 50000.0 * 10;
      const double prob = exp((double)(after - before) / temp);

      if (prob > rand01()) {
      }
      else {
        swap(A[ra1][ra2], A[ra3][ra4]);
      }
    }
  }
}

ll Solve(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  InitGuruGuru();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  for (int t = 0; t < T; t++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        A[i][j] = i * N + j;
      }
    }

    Convert();
    Sort();
    GuruGuru2(t);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (mode == 0) {
          cout << A[i][j];

          if (j < N - 1) {
            cout << " ";
          }
          else {
            cout << endl;
          }
        }
        else {
          ofs << A[i][j];

          if (j < N - 1) {
            ofs << " ";
          }
          else {
            ofs << endl;
          }
        }
      }
    }

    if (mode == 0) {
      cout.flush();

      for (int i = 0; i < NN; i++) {
        for (int j = 0; j < M; j++) {
          cin >> X[i][j];
        }
      }
    }
    else {
      GetReturn(t);
    }
  }

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }

  return score;

  // 初期解生成
  Initialize();

  // 解答を出力
  Output(ofs);

  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
