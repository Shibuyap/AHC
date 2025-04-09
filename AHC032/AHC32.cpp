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
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

namespace /* 乱数ライブラリ */
{
  static uint32_t Rand()
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
  static double rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.9;
int mode;
clock_t startTime, endTime;

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int n = 9;
const int m = 20;
const int T = 81;
const int MOD = 998244353;
int a[n][n];
int initA[n][n];
int s[m][3][3];
int ans[T + 100][3];
ll ansScore;

int best_a[n][n];
ll best_ans[T][3];
ll best_ansScore;

void CopyToBest()
{
  best_ansScore = ansScore;
  rep(i, T)
  {
    rep(j, 3)
    {
      best_ans[i][j] = ans[i][j];
    }
  }
  rep(i, n)
  {
    rep(j, n)
    {
      best_a[i][j] = a[i][j];
    }
  }
}

void CopyFromBest()
{
  ansScore = best_ansScore;
  rep(i, T)
  {
    rep(j, 3)
    {
      ans[i][j] = best_ans[i][j];
    }
  }
  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = best_a[i][j];
    }
  }
}

// スコア計算
ll CalcScore()
{
  ll res = 0;

  rep(i, n)
  {
    rep(j, n)
    {
      res += a[i][j] % MOD;
    }
  }

  return res;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  rep(i, T)
  {
    rep(j, 3)
    {
      ans[i][j] = -1;
    }
  }
}

void InitializeAns()
{
  rep(i, T)
  {
    rep(j, 3)
    {
      ans[i][j] = -1;
    }
  }
  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = initA[i][j];
    }
  }
  ansScore = CalcScore();
}

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
    int _n, _m, _k;
    cin >> _n >> _m >> _k;
    rep(i, n)
    {
      rep(j, n)
      {
        cin >> a[i][j];
      }
    }
    rep(i, m)
    {
      rep(j, 3)
      {
        rep(k, 3)
        {
          cin >> s[i][j][k];
        }
      }
    }
  }
  // ファイル入力する
  else {
    int _n, _m, _k;
    ifs >> _n >> _m >> _k;
    rep(i, n)
    {
      rep(j, n)
      {
        ifs >> a[i][j];
      }
    }
    rep(i, m)
    {
      rep(j, 3)
      {
        rep(k, 3)
        {
          ifs >> s[i][j][k];
        }
      }
    }
  }

  rep(i, n)
  {
    rep(j, n)
    {
      initA[i][j] = a[i][j];
    }
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

// 初期解生成
void Initialize()
{
  rep(i, T)
  {
    ans[i][0] = -1;
  }
  ansScore = CalcScore();
  CopyToBest();
}

// 解答出力
void Output(ofstream& ofs)
{
  int L = 0;
  rep(i, T)
  {
    if (ans[i][0] == -1) continue;
    L++;
  }
  if (mode == 0) {
    cout << L << endl;
    rep(i, T)
    {
      if (ans[i][0] == -1) continue;
      rep(j, 3) cout << ans[i][j] << ' ';
      cout << endl;
    }
  }
  else {
    ofs << L << endl;
    rep(i, T)
    {
      if (ans[i][0] == -1) continue;
      rep(j, 3) ofs << ans[i][j] << ' ';
      ofs << endl;
    }
  }
}

double nowTime = 0;
double startTemp = 1001001001;
double endTemp = 0.1;
void Method1_1()
{
  int raT = Rand() % T;
  int raX = Rand() % 7;
  int raY = Rand() % 7;
  int raM = Rand() % m;
  if (Rand() % 2 == 0) { raM = -1; }
  if (ans[raT][0] == raM) return;
  ll diff = 0;
  if (ans[raT][0] != -1) {
    rep(i, 3)
    {
      rep(j, 3)
      {
        diff -= s[ans[raT][0]][i][j];
        if (s[ans[raT][0]][i][j] > a[ans[raT][1] + i][ans[raT][2] + j]) { diff += MOD; }
      }
    }
  }
  if (raM != -1) {
    rep(i, 3)
    {
      rep(j, 3)
      {
        diff += s[raM][i][j];
        if (s[raM][i][j] + a[raX + i][raY + j] > MOD) diff -= MOD;
      }
    }
  }

  double temp = (startTemp + (endTemp - startTemp) * nowTime / TL);
  const double prob = exp((double)diff / temp);
  // if (prob > rand01()) {
  if (diff >= 0) {
    if (ans[raT][0] != -1) {
      rep(i, 3)
      {
        rep(j, 3)
        {
          a[ans[raT][1] + i][ans[raT][2] + j] -= s[ans[raT][0]][i][j];
          if (a[ans[raT][1] + i][ans[raT][2] + j] < 0) { a[ans[raT][1] + i][ans[raT][2] + j] += MOD; }
        }
      }
    }
    if (raM != -1) {
      rep(i, 3)
      {
        rep(j, 3)
        {
          a[raX + i][raY + j] += s[raM][i][j];
          if (a[raX + i][raY + j] > MOD) a[raX + i][raY + j] -= MOD;
        }
      }
    }

    ans[raT][0] = raM;
    ans[raT][1] = raX;
    ans[raT][2] = raY;

    ansScore += diff;
    cout << ansScore << endl;
    if (ansScore > best_ansScore) { CopyToBest(); }
  }
}

void Method1()
{
  CopyFromBest();
  int loopCount = 0;
  nowTime = GetNowTime();
  while (true) {
    loopCount++;
    if (loopCount % 10000 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) { break; }
    }

    Method1_1();
  }

  CopyFromBest();
  if (mode != 0) { cout << "loop = " << loopCount << ", ansScore = " << CalcScore() << endl; }
}

int use[3][3];
void Rule1(int x, int y, int dir1, int dir2)
{
  rep(k, 3) rep(l, 3) use[k][l] = 0;
  if (dir1 == 0 && dir2 == 0) {
    if (x == n - 3) {
      if (y == n - 3) {
        rep(i, 3)
        {
          rep(j, 3)
          {
            use[i][j] = 1;
          }
        }
      }
      else {
        rep(i, 3) use[i][0] = 1;
      }
    }
    else {
      if (y == n - 3) {
        rep(l, 3) use[0][l] = 1;
      }
      else {
        use[0][0] = 1;
      }
    }
  }
  if (dir1 == 0 && dir2 == 1) {
    if (x == n - 3) {
      if (y == 0) {
        rep(i, 3)
        {
          rep(j, 3)
          {
            use[i][j] = 1;
          }
        }
      }
      else {
        rep(i, 3) use[i][2] = 1;
      }
    }
    else {
      if (y == 0) {
        rep(l, 3) use[0][l] = 1;
      }
      else {
        use[0][2] = 1;
      }
    }
  }
  if (dir1 == 1 && dir2 == 0) {
    if (x == 0) {
      if (y == n - 3) {
        rep(i, 3)
        {
          rep(j, 3)
          {
            use[i][j] = 1;
          }
        }
      }
      else {
        rep(i, 3) use[i][0] = 1;
      }
    }
    else {
      if (y == n - 3) {
        rep(l, 3) use[2][l] = 1;
      }
      else {
        use[2][0] = 1;
      }
    }
  }
  if (dir1 == 1 && dir2 == 1) {
    if (x == 0) {
      if (y == 0) {
        rep(i, 3)
        {
          rep(j, 3)
          {
            use[i][j] = 1;
          }
        }
      }
      else {
        rep(i, 3) use[i][2] = 1;
      }
    }
    else {
      if (y == 0) {
        rep(l, 3) use[2][l] = 1;
      }
      else {
        use[2][2] = 1;
      }
    }
  }
}

void Rule2(int i, int j, int dir1)
{
  rep(p, 3)
  {
    rep(q, 3)
    {
      use[p][q] = 0;
      if (dir1 == 0) {
        if (i == n - 3) use[p][q] = 1;
        if (p == 0) use[p][q] = 1;
      }
      else {
        if (i == 0) use[p][q] = 1;
        if (p == 2) use[p][q] = 1;
      }
    }
  }
}

ll maxSum;
ll ma[3][3];
ll maAnssArr[10];
int maAnsCount = 0;

ll now[3][3];
ll anssArr[10];
void Method2DFS(int mm, int cnt, int lim)
{
  if (cnt == lim) return;
  ll keep[3][3];
  rep(p, 3) rep(q, 3) keep[p][q] = now[p][q];
  srep(i, mm, m)
  {
    anssArr[cnt] = i;
    cnt++;
    ll tmpSum = 0;
    rep(p, 3)
    {
      rep(q, 3)
      {
        now[p][q] = (now[p][q] + s[i][p][q]) % MOD;
        if (use[p][q]) tmpSum += now[p][q];
      }
    }
    if (tmpSum > maxSum) {
      maxSum = tmpSum;
      rep(p, 3)
      {
        rep(q, 3)
        {
          ma[p][q] = now[p][q];
        }
      }
      rep(j, cnt) maAnssArr[j] = anssArr[j];
      maAnsCount = cnt;
    }

    Method2DFS(i, cnt, lim);

    cnt--;
    rep(p, 3)
    {
      rep(q, 3)
      {
        now[p][q] = keep[p][q];
      }
    }
  }
}

void Method2(double timeLimit)
{
  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) break;

    ll hosyou = Rand() % 200000000 + MOD - 200000000;

    int ng = 0;
    InitializeAns();
    int cnt = 0;
    int dir1 = Rand() % 2;
    int dir2 = Rand() % 2;
    rep(ii, n - 2)
    {
      int i = ii;
      if (dir1) i = n - 3 - ii;
      rep(jj, n - 2)
      {
        int j = jj;
        if (dir2) j = n - 3 - jj;
        maAnsCount = 0;
        maxSum = 0;
        rep(k, 3)
        {
          rep(l, 3)
          {
            ma[k][l] = a[i + k][j + l];
            now[k][l] = ma[k][l];
          }
        }
        Rule1(i, j, dir1, dir2);
        int useCount = 0;
        rep(k, 3) rep(l, 3) useCount += use[k][l];
        rep(k, 3)
        {
          rep(l, 3)
          {
            if (use[k][l]) { maxSum += ma[k][l]; }
          }
        }

        if (useCount == 1) {
          Method2DFS(0, 0, 1);
          if (maxSum < hosyou) {
            rep(p, 3)
            {
              rep(q, 3)
              {
                now[p][q] = a[i + p][j + q];
              }
            }
            Method2DFS(0, 0, 2);
          }
        }
        else if (useCount <= 3) {
          Method2DFS(0, 0, 3);
        }
        else {
          if (cnt + 6 > T) {
            ng = 1;
            break;
          }
          Method2DFS(0, 0, 6);
        }

        rep(k, maAnsCount)
        {
          int ansM = maAnssArr[k];
          ans[cnt][0] = ansM;
          ans[cnt][1] = i;
          ans[cnt][2] = j;
          cnt++;
        }
        if (cnt > T) {
          ng = 1;
          break;
        }
        rep(k, 3)
        {
          rep(l, 3)
          {
            a[i + k][j + l] = ma[k][l];
          }
        }
      }
      if (ng) break;
    }
    if (ng) continue;

    ansScore = CalcScore();
    if (ansScore > best_ansScore) {
      // cout << hosyou << ' ' << ansScore << endl;
      CopyToBest();
    }
  }

  if (mode != 0) cout << loopCount << endl;
}

void Method3(double timeLimit)
{
  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) break;

    ll hosyou = Rand() % 200000000 + MOD - 200000000;

    int aaa[49];
    rep(i, 49) aaa[i] = i;
    FisherYates(aaa, 49);
    int last[n][n];
    rep(i, n)
    {
      rep(j, n)
      {
        last[i][j] = 0;
      }
    }
    rep(ii, 49)
    {
      int x = aaa[ii] / 7;
      int y = aaa[ii] % 7;
      rep(p, 3)
      {
        rep(q, 3)
        {
          last[x + p][y + q] = ii;
        }
      }
    }

    int ng = 0;
    InitializeAns();
    int cnt = 0;

    rep(ii, 49)
    {
      int i = aaa[ii] / 7;
      int j = aaa[ii] % 7;

      maAnsCount = 0;
      maxSum = 0;
      rep(k, 3)
      {
        rep(l, 3)
        {
          ma[k][l] = a[i + k][j + l];
          now[k][l] = ma[k][l];
        }
      }

      rep(p, 3)
      {
        rep(q, 3)
        {
          use[p][q] = 0;
          if (last[i + p][j + q] == ii) { use[p][q] = 1; }
        }
      }
      int useCount = 0;
      rep(k, 3) rep(l, 3) useCount += use[k][l];
      rep(k, 3)
      {
        rep(l, 3)
        {
          if (use[k][l]) { maxSum += ma[k][l]; }
        }
      }

      if (useCount == 1) {
        Method2DFS(0, 0, 1);
        if (maxSum < hosyou) {
          rep(p, 3)
          {
            rep(q, 3)
            {
              now[p][q] = a[i + p][j + q];
            }
          }
          Method2DFS(0, 0, 2);
        }
      }
      else if (useCount <= 6) {
        Method2DFS(0, 0, 3);
      }
      else {
        if (cnt + 6 > T) {
          ng = 1;
          break;
        }
        Method2DFS(0, 0, 6);
      }

      rep(k, maAnsCount)
      {
        int ansM = maAnssArr[k];
        ans[cnt][0] = ansM;
        ans[cnt][1] = i;
        ans[cnt][2] = j;
        cnt++;
      }
      if (cnt > T) {
        ng = 1;
        break;
      }
      rep(k, 3)
      {
        rep(l, 3)
        {
          a[i + k][j + l] = ma[k][l];
        }
      }
    }

    if (ng) continue;

    ansScore = CalcScore();
    if (ansScore > best_ansScore) {
      // cout << hosyou << ' ' << ansScore << endl;
      CopyToBest();
    }
  }

  if (mode != 0) cout << loopCount << endl;
}

int keepA[n][n];
int keepAns[110][3];
int baseA[n][n];

void Method4(double timeLimit)
{
  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) break;

    ll hosyou = Rand() % 200000000 + MOD - 200000000;

    int ng = 0;
    InitializeAns();
    int cnt = 0;
    int dir1 = Rand() % 2;
    int dir2 = Rand() % 2;
    // dir1     = 0;
    dir2 = 0;
    rep(ii, n - 2)
    {
      int i = ii;
      if (dir1) i = n - 3 - ii;
      if (ii == n - 3 && cnt + 3 * 6 + 4 > T) {
        ng = 1;
        break;
      }

      int nowCnt = cnt;
      ll maPosSum = 0;
      int maCntTail = 0;
      rep(p, n)
      {
        rep(q, n)
        {
          baseA[p][q] = a[p][q];
        }
      }

      rep(jjj, n - 2)
      {
        int ng2 = 0;
        rep(jj, jjj)
        {
          int j = jj;
          dir2 = 0;
          if (dir2) j = n - 3 - jj;
          maAnsCount = 0;
          maxSum = 0;
          rep(k, 3)
          {
            rep(l, 3)
            {
              ma[k][l] = a[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          Rule1(i, j, dir1, dir2);
          int useCount = 0;
          rep(k, 3) rep(l, 3) useCount += use[k][l];
          rep(k, 3)
          {
            rep(l, 3)
            {
              if (use[k][l]) { maxSum += ma[k][l]; }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              rep(p, 3)
              {
                rep(q, 3)
                {
                  now[p][q] = a[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            Method2DFS(0, 0, 4);
          }

          rep(k, maAnsCount)
          {
            int ansM = maAnssArr[k];
            ans[cnt][0] = ansM;
            ans[cnt][1] = i;
            ans[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          rep(k, 3)
          {
            rep(l, 3)
            {
              a[i + k][j + l] = ma[k][l];
            }
          }
        }
        if (ng) break;

        for (int jj = n - 3; jj > jjj; jj--) {
          int j = jj;
          dir2 = 1;
          maAnsCount = 0;
          maxSum = 0;
          rep(k, 3)
          {
            rep(l, 3)
            {
              ma[k][l] = a[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          Rule1(i, j, dir1, dir2);
          int useCount = 0;
          rep(k, 3) rep(l, 3) useCount += use[k][l];
          rep(k, 3)
          {
            rep(l, 3)
            {
              if (use[k][l]) { maxSum += ma[k][l]; }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              rep(p, 3)
              {
                rep(q, 3)
                {
                  now[p][q] = a[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            Method2DFS(0, 0, 4);
          }

          rep(k, maAnsCount)
          {
            int ansM = maAnssArr[k];
            ans[cnt][0] = ansM;
            ans[cnt][1] = i;
            ans[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          rep(k, 3)
          {
            rep(l, 3)
            {
              a[i + k][j + l] = ma[k][l];
            }
          }
        }

        {
          int j = jjj;
          maAnsCount = 0;
          maxSum = 0;
          rep(k, 3)
          {
            rep(l, 3)
            {
              ma[k][l] = a[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          // Rule1(i, j, dir1, dir2);
          Rule2(i, j, dir1);
          int useCount = 0;
          rep(k, 3) rep(l, 3) useCount += use[k][l];
          rep(k, 3)
          {
            rep(l, 3)
            {
              if (use[k][l]) { maxSum += ma[k][l]; }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              rep(p, 3)
              {
                rep(q, 3)
                {
                  now[p][q] = a[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            int num = min(6, T - cnt);
            Method2DFS(0, 0, num);
          }

          rep(k, maAnsCount)
          {
            int ansM = maAnssArr[k];
            ans[cnt][0] = ansM;
            ans[cnt][1] = i;
            ans[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          rep(k, 3)
          {
            rep(l, 3)
            {
              a[i + k][j + l] = ma[k][l];
            }
          }
        }

        if (ng) break;

        ll tmpPosSum = 0;

        if (dir1 == 0) {
          rep(j, n)
          {
            tmpPosSum += a[i][j];
          }
          if (i == n - 3) {
            rep(j, n)
            {
              tmpPosSum += a[i + 1][j];
              tmpPosSum += a[i + 2][j];
            }
          }
        }
        else {
          rep(j, n)
          {
            tmpPosSum += a[i + 2][j];
          }
          if (i == 0) {
            rep(j, n)
            {
              tmpPosSum += a[i][j];
              tmpPosSum += a[i + 1][j];
            }
          }
        }
        if (tmpPosSum > maPosSum) {
          maPosSum = tmpPosSum;
          srep(t, nowCnt, cnt)
          {
            rep(k, 3) keepAns[t][k] = ans[t][k];
          }
          rep(p, n)
          {
            rep(q, n)
            {
              keepA[p][q] = a[p][q];
            }
          }
          maCntTail = cnt;
        }

        cnt = nowCnt;
        rep(p, n)
        {
          rep(q, n)
          {
            a[p][q] = baseA[p][q];
          }
        }
      }

      if (ng) break;

      srep(t, nowCnt, maCntTail)
      {
        rep(k, 3)
        {
          ans[t][k] = keepAns[t][k];
        }
      }
      cnt = maCntTail;
      rep(p, n)
      {
        rep(q, n)
        {
          a[p][q] = keepA[p][q];
        }
      }
    }
    if (ng) continue;

    srep(t, cnt, T) ans[t][0] = -1;
    ansScore = CalcScore();
    // cout << ansScore << endl;
    if (ansScore > best_ansScore) {
      // cout << hosyou << ' ' << ansScore << endl;
      CopyToBest();
    }
  }

  if (mode != 0) cout << loopCount << endl;
}

ll Solve(int probNum)
{
  startTime = clock();
  endTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();
  // Method2(TL);
  Method4(TL);
  // Method3(TL);
  CopyFromBest();
  // Method1();

  CopyFromBest();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) { ofs.close(); }

  ll score = 0;
  if (mode != 0) { score = CalcScore(); }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
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
