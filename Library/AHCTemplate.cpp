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
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

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

static double Rand01() {
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}

static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
}

// [l, r]
static uint32_t RandRange(uint32_t l, uint32_t r)
{
  return l + Rand() % (r - l + 1);
}


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタ
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
std::chrono::steady_clock::time_point startTimeClock;

void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

int n;

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

bool IsNG(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
  }
  else {
    // ファイル入力
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

// ナイーブな解法
void Method1()
{

}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp[10];
  double EndTemp;
  double MultipleValue;
  int Partition[10];
};

void SimulatedAnnealing(Hypers hypers)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp[0];
  const double END_TEMP = hypers.EndTemp;


  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) break;
    }

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    // 近傍解作成
    int raMode = Rand() % hypers.Partition[1];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {

    }
    else if (raMode < hypers.Partition[1]) {

    }

    // スコア計算
    double tmpScore = CalcScore();

    // 焼きなまし
    double diffScore = (tmpScore - ansScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > Rand01()) {
      // 採用
      ansScore = tmpScore;

      // Best解よりもいいか
      if (ansScore > best_ansScore) {
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {
      }
      else if (raMode < hypers.Partition[1]) {
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
  }

  CopyToAns();
}

// 問題を解く関数
ll Solve(int problem_num, Hypers hypers)
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

  // 焼きなまし
  //SimulatedAnnealing(hypers);

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
    Rand();
  }

  mode = 2;

  Hypers hypers;
  hypers.StartTemp[0] = 2048.0;
  hypers.StartTemp[1] = 2048.0;
  hypers.StartTemp[2] = 2048.0;
  hypers.StartTemp[3] = 2048.0;
  hypers.StartTemp[4] = 2048.0;
  hypers.StartTemp[5] = 2048.0;
  hypers.StartTemp[6] = 2048.0;
  hypers.StartTemp[7] = 2048.0;
  hypers.StartTemp[8] = 2048.0;
  hypers.StartTemp[9] = 2048.0;
  hypers.EndTemp = 0.0;
  hypers.MultipleValue = 12345.0;
  hypers.Partition[0] = 100;
  hypers.Partition[1] = 200;
  hypers.Partition[2] = 300;
  hypers.Partition[3] = 400;
  hypers.Partition[4] = 500;
  hypers.Partition[5] = 600;
  hypers.Partition[6] = 700;
  hypers.Partition[7] = 800;
  hypers.Partition[8] = 900;
  hypers.Partition[9] = 1000;

  if (mode == 0) {
    Solve(0, hypers);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 15)
    {
      ll score = Solve(i, hypers);
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
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers newHypers;
      newHypers.StartTemp[0] = pow(2.0, Rand01() * 20);
      newHypers.EndTemp = 0.0;
      newHypers.MultipleValue = pow(2.0, Rand01() * 20);
      newHypers.Partition[0] = Rand() % 101;

      ll sum = 0;
      srep(i, 0, 15)
      {
        ll score = Solve(i, newHypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << newHypers.StartTemp[0]
        << ", EndTemp = " << newHypers.EndTemp
        << ", MultipleValue = " << newHypers.MultipleValue
        << ", Partition1 = " << newHypers.Partition[0]
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = newHypers;
      }

      loop++;
    }
  }

  return 0;
}
