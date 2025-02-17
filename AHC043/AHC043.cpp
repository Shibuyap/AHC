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

// ループの簡略化マクロ
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t RandXor()
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() { return (RandXor() + 0.5) * (1.0 / UINT_MAX); }

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化（使用されていない）
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// 非常に大きな値
const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int ex[13] = { -2,-1,-1,-1,0,0,0,0,0,1,1,1,2 };
const int ey[13] = { 0,-1,0,1,-2,-1,0,1,2,-1,0,1,0 };

double TL = 2.8; // 時間制限（Time Limit）
int mode;        // 実行モード
std::chrono::steady_clock::time_point startTimeClock; // 時間計測用

// 時間計測をリセットする関数
void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}


const int MAX_M = 1600;

const int n = 50;
int m;
int k;
const int T = 800;

vector<int> sx, sy, tx, ty;
int a[n][n];
int b[n][n];

int ansScore;

int best_ansScore;

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  sx.clear();
  sy.clear();
  tx.clear();
  ty.clear();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn, tt;
    cin >> nn >> m >> k >> tt;
    sx.resize(m);
    sy.resize(m);
    tx.resize(m);
    ty.resize(m);
    rep(i, m)
    {
      cin >> sx[i] >> sy[i] >> tx[i] >> ty[i];
    }
  }
  else {
    // ファイル入力
    int nn, tt;
    ifs >> nn >> m >> k >> tt;
    sx.resize(m);
    sy.resize(m);
    tx.resize(m);
    ty.resize(m);
    rep(i, m)
    {
      ifs >> sx[i] >> sy[i] >> tx[i] >> ty[i];
    }
  }

  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = -1;
      b[i][j] = -1;
    }
  }
  rep(i, m)
  {
    a[sx[i]][sy[i]] = i;
    b[tx[i]][ty[i]] = i;
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
ll CalcScore()
{
  ll res = 0;
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
  }
  else {
    // ファイル出力
  }
}

struct EkiMap
{
  int AFlags[MAX_M] = {};
  int BFlags[MAX_M] = {};
  int Count = 0;
  int Stations[50][2];
  int StationCount = 10;
};

// ナイーブな解法
void Method1()
{
  EkiMap ekiMap;
  EkiMap best_ekiMap;

  rep(i, ekiMap.StationCount)
  {
    int x = RandXor() % n;
    int y = RandXor() % n;
    ekiMap.Stations[i][0] = x;
    ekiMap.Stations[i][1] = y;

    rep(j, 13)
    {
      int nx = x + ex[j];
      int ny = y + ey[j];
    }
  }


  double nowTime = GetNowTime();
  const double START_TEMP = 2048.0;
  const double END_TEMP = 0.1;

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL / 3) break;
    }
    loop++;



    double diffScore = 0;

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
    double prob = exp(diffScore / temp);

    if (prob > Rand01()) {
      // 採用
    }
    else {
      // 元に戻す
    }
  }
}

// 問題を解く関数
ll Solve(int problem_num)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

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
    RandXor();
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
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
