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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数
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

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);


const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
std::chrono::steady_clock::time_point startTime, endTime;

void ResetTime()
{
  startTime = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTime - startTime;
  return elapsed.count();
}


const int MAX_N = 30;

const int n = 200;
const int m = 10;
int b[m][n];
int bCount[m];
int c[n][2];
int init_b[m][n];
int init_c[n][2];
int init_bCount[m];

int ansScore;
int ans[5100][2];
int ansSize;

int best_ansScore;

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ansScore = 0;
}

// 入力受け取り
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    int _n, _m;
    cin >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> b[i][j];
        b[i][j]--;
        c[b[i][j]][0] = i;
        c[b[i][j]][1] = j;
      }
    }
  }
  // ファイル入力する
  else {
    int _n, _m;
    ifs >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        ifs >> b[i][j];
        b[i][j]--;
        c[b[i][j]][0] = i;
        c[b[i][j]][1] = j;
      }
    }
  }

  rep(i, m)
  {
    rep(j, n / m) { init_b[i][j] = b[i][j]; }
  }
  rep(i, n)
  {
    rep(j, 2) { init_c[i][j] = c[i][j]; }
  }
  rep(i, m)
  {
    bCount[i] = n / m;
    init_bCount[i] = bCount[i];
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコア計算
int CalcScore()
{
  int tmp_b[m][n];
  int tmp_c[n][2];
  int tmp_bCount[m];

  int cnt = 0;
  rep(i, m)
  {
    rep(j, n / m) {
      tmp_b[i][j] = init_b[i][j];
    }
  }
  rep(i, n)
  {
    rep(j, 2) {
      tmp_c[i][j] = init_c[i][j];
    }
  }
  rep(i, m) {
    tmp_bCount[i] = init_bCount[i];
  }

  int res = 10000;
  rep(i, ansSize)
  {
    int num = ans[i][0];
    int x = tmp_c[num][0];
    int y = tmp_c[num][1];
    int nx = ans[i][1];
    if (nx == -1) {
      tmp_bCount[x]--;
      cnt++;
    }
    else {
      res--;
      rep(j, tmp_bCount[x] - y)
      {
        int num2 = tmp_b[x][y + j];
        tmp_c[num2][0] = nx;
        tmp_c[num2][1] = tmp_bCount[nx];
        tmp_b[nx][tmp_bCount[nx]] = num2;
        tmp_bCount[nx]++;
        res--;
      }
      tmp_bCount[x] = y;
    }
  }
  if (cnt != n) return -1;
  return res;
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    rep(i, ansSize) { cout << ans[i][0] + 1 << ' ' << ans[i][1] + 1 << endl; }
  }
  else {
    // ファイル出力
    rep(i, ansSize) { ofs << ans[i][0] + 1 << ' ' << ans[i][1] + 1 << endl; }
  }
}

int GetNum(int x) {
  if (bCount[x] == 0) {
    return -1;
  }
  return b[x][bCount[x] - 1];
}

void CarryOut(int& carryOutCount) {
  int x = c[carryOutCount][0];
  bCount[x]--;

  ans[ansSize][0] = carryOutCount;
  ans[ansSize][1] = -1;
  ansSize++;
  carryOutCount++;
}

void Move(int x, int nx, int y)
{
  int num = b[x][y];
  ans[ansSize][0] = num;
  ans[ansSize][1] = nx;
  ansSize++;

  rep(j, bCount[x] - y)
  {
    int num2 = b[x][y + j];
    b[nx][bCount[nx]] = b[x][y + j];
    c[num2][0] = nx;
    c[num2][1] = bCount[nx];
    bCount[nx]++;
  }
  bCount[x] = y;
}

// 1列ずつソートしていく
void Method1()
{
  ansSize = 0;
  int carryOutCount = 0;

  while (carryOutCount < n) {
    rep(i, m)
    {
      // 一旦すべて取り出す
      while (bCount[i] > 0) {
        int num = GetNum(i);

        if (num == carryOutCount) {
          CarryOut(carryOutCount);
          continue;
        }

        int idx = -1;
        int idxNum = -1;
        rep(j, m) {
          if (j == i)continue;

          int jNum = GetNum(j);

          if (idx == -1) {
            idx = j;
            idxNum = jNum;
          }
          else {
            if (idxNum == -1) {
              if (jNum > num) {
                idx = j;
                idxNum = jNum;
              }
            }
            else if (idxNum < num) {
              if (jNum == -1 || jNum > num) {
                idx = j;
                idxNum = jNum;
              }
              else if (jNum > idxNum) {
                idx = j;
                idxNum = jNum;
              }
            }
            else {
              if (jNum > num && jNum < idxNum) {
                idx = j;
                idxNum = jNum;
              }
            }
          }
        }

        Move(i, idx, bCount[i] - 1);
      }

      // 戻す
      int now = 999;
      while (true) {
        int ma = -1;
        int idx = -1;
        rep(j, m) {
          if (j == i)continue;
          if (bCount[j] == 0)continue;

          int jNum = GetNum(j);
          if (jNum < now && ma < jNum) {
            ma = jNum;
            idx = j;
          }
        }

        if (idx == -1)break;

        int y = bCount[idx] - 1;
        while (y > 0) {
          int num = b[idx][y];
          int nextNum = b[idx][y - 1];
          if (num < nextNum && nextNum < now) {
            y--;
          }
          else {
            break;
          }
        }

        Move(idx, i, y);
        now = GetNum(i);
      }

      // 運び出せる箱があるか確認
      while (true) {
        int ok = 0;
        rep(j, m) {
          if (GetNum(j) == carryOutCount) {
            CarryOut(carryOutCount);
            ok = 1;
            break;
          }
        }
        if (ok == 0)break;
      }
    }
  }
}

ll Solve(int probNum)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Method1();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
