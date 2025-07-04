﻿#include <algorithm>
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

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

static uint32_t rand_xorshift()
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


static double rand_01()
{
  return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
}


static double RandRange(double l, double r)
{
  return l + (r - l) * rand_01();
}


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = rand_xorshift() % (i + 1);
    int tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// 2次元
int queueArr2[10000][2];
int queueHead2 = 0;
int queueTail2 = 0;
void ClearQueue()
{
  queueHead2 = 0;
  queueTail2 = 0;
}
int FrontX()
{
  return queueArr2[queueHead2][0];
}
int FrontY()
{
  return queueArr2[queueHead2][1];
}
void Push(int x, int y)
{
  queueArr2[queueTail2][0] = x;
  queueArr2[queueTail2][1] = y;
  queueTail2++;
}
void Pop()
{
  queueHead2++;
}
int Size()
{
  return queueTail2 - queueHead2;
}



const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;


const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;

std::chrono::steady_clock::time_point start_time_clock;

void start_timer()
{
  start_time_clock = std::chrono::steady_clock::now();
}

double get_elapsed_time()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
  return elapsed.count();
}

const int n = 50;
const int m = 100;
int c[n + 2][n + 2];
int g[m + 1][m + 1];

int ansScore;
int d[n + 2][n + 2];

int best_ansScore;
int best_d[n + 2][n + 2];

void CopyToBest()
{
  best_ansScore = ansScore;
  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      best_d[i][j] = d[i][j];
    }
  }
}

void CopyToAns()
{
  ansScore = best_ansScore;
  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      d[i][j] = best_d[i][j];
    }
  }
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
void Input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn, mm;
    cin >> nn >> mm;
    for (int i = 1; i < n + 1; ++i) {
      for (int j = 1; j < n + 1; ++j) {
        cin >> c[i][j];
      }
    }
  }
  else {
    // ファイル入力
    int nn, mm;
    ifs >> nn >> mm;
    for (int i = 1; i < n + 1; ++i) {
      for (int j = 1; j < n + 1; ++j) {
        ifs >> c[i][j];
      }
    }
  }

  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      d[i][j] = c[i][j];
    }
  }

  for (int i = 0; i < m + 1; ++i) {
    for (int j = 0; j < (m + 1); ++j) {
      g[i][j] = 0;
    }
  }
  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 1; ++j) {
      g[c[i][j]][c[i][j + 1]] = 1;
      g[c[i][j + 1]][c[i][j]] = 1;
    }
  }
  for (int i = 0; i < n + 1; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      g[c[i][j]][c[i + 1][j]] = 1;
      g[c[i + 1][j]][c[i][j]] = 1;
    }
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int case_num, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore()
{
  int res = 1;
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < n + 1; ++j) {
      if (d[i][j] == 0) {
        res++;
      }
    }
  }
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    for (int i = 1; i < n + 1; ++i) {
      for (int j = 1; j < n + 1; ++j) {
        cout << d[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    for (int i = 1; i < n + 1; ++i) {
      for (int j = 1; j < n + 1; ++j) {
        ofs << std::setw(3) << d[i][j] << ' ';
      }
      ofs << endl;
    }
  }
}

int checkG[m + 1][m + 1];
int checkVisited[n + 2][n + 2];
int checkVisited2[m + 1];

bool Check()
{
  for (int i = 0; i < m + 1; ++i) {
    for (int j = 0; j < (m + 1); ++j) {
      checkG[i][j] = 0;
    }
  }
  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 1; ++j) {
      if (d[i][j] == d[i][j + 1])continue;
      if (g[d[i][j]][d[i][j + 1]] == 0) return false;
      if (g[d[i][j + 1]][d[i][j]] == 0) return false;
      checkG[d[i][j]][d[i][j + 1]] = 1;
      checkG[d[i][j + 1]][d[i][j]] = 1;
    }
  }
  for (int i = 0; i < n + 1; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      if (d[i][j] == d[i + 1][j])continue;
      if (g[d[i][j]][d[i + 1][j]] == 0) return false;
      if (g[d[i + 1][j]][d[i][j]] == 0) return false;
      checkG[d[i][j]][d[i + 1][j]] = 1;
      checkG[d[i + 1][j]][d[i][j]] = 1;
    }
  }

  for (int i = 0; i < m + 1; ++i) {
    for (int j = i + 1; j < m + 1; ++j) {
      if (checkG[i][j] != g[i][j])return false;
    }
  }

  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      checkVisited[i][j] = 0;
    }
  }
  for (int i = 0; i < m + 1; ++i) {
    checkVisited2[i] = 0;
  }
  ClearQueue();
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < n + 1; ++j) {
      if (d[i][j] == 0 || checkVisited[i][j] == 1)continue;
      int num = d[i][j];
      if (checkVisited2[num] == 1)return false;
      checkVisited2[num] = 1;
      checkVisited[i][j] = 1;
      Push(i, j);
      while (Size() > 0) {
        int x = FrontX();
        int y = FrontY();
        Pop();
        for (int k = 0; k < 4; ++k) {
          int nx = x + dx[k];
          int ny = y + dy[k];
          if (d[nx][ny] == num && checkVisited[nx][ny] == 0) {
            checkVisited[nx][ny] = 1;
            Push(nx, ny);
          }
        }
      }
    }
  }

  return true;
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition[10];
};

int keep[n + 2][n + 2];
void KeepD()
{
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < n + 1; ++j) {
      keep[i][j] = d[i][j];
    }
  }
}

void SimulatedAnnealing(Hypers hypers)
{
  ansScore = CalcScore();

  CopyToBest();

  double nowTime = get_elapsed_time();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;


  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = get_elapsed_time();
      if (nowTime > TL) { break; }
    }

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    // 近傍解作成
    int raMode = rand_xorshift() % hypers.Partition[2];
    int ra1, ra2, ra3, ra4, ra5, raDir;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      int raDir = rand_xorshift() % 2;
      int ra1 = rand_xorshift() % n + 1;
      int ok = 0;
      while (ok == 0) {
        raDir = rand_xorshift() % 2;
        ra1 = rand_xorshift() % n + 1;
        if (raDir == 0) {
          for (int j = 1; j < n + 1; ++j) {
            if (d[ra1][j] != 0)ok = 1;
          }
        }
        else {
          for (int i = 1; i < n + 1; ++i) {
            if (d[i][ra1] != 0)ok = 1;
          }
        }
      }

      if (raDir == 0) {
        for (int j = 1; j < n + 1; ++j) {
          if (g[d[ra1 - 1][j]][d[ra1 + 1][j]] == 0) {
            ok = 0;
          }
        }
      }
      else {
        for (int i = 1; i < n + 1; ++i) {
          if (g[d[i][ra1 - 1]][d[i][ra1 + 1]] == 0) {
            ok = 0;
          }
        }
      }

      if (ok == 0)continue;
      KeepD();

      if (raDir == 0) {
        for (int i = ra1; i < n + 1; ++i) {
          for (int j = 1; j < n + 1; ++j) {
            d[i][j] = d[i + 1][j];
          }
        }
      }
      else {
        for (int j = ra1; j < n + 1; ++j) {
          for (int i = 1; i < n + 1; ++i) {
            d[i][j] = d[i][j + 1];
          }
        }
      }
    }
    else if (raMode < hypers.Partition[1]) {

      while (true) {
        ra1 = rand_xorshift() % (n - 1) + 1;
        ra2 = rand_xorshift() % (n - 1) + 1;
        raDir = rand_xorshift() % 4;
        if (d[ra1][ra2] == 0 || d[ra1][ra2 + 1] == 0 || d[ra1 + 1][ra2] == 0 || d[ra1 + 1][ra2 + 1] == 0)continue;
        break;
      }

      keep1 = d[ra1][ra2];
      keep2 = d[ra1][ra2 + 1];
      keep3 = d[ra1 + 1][ra2];
      keep4 = d[ra1 + 1][ra2 + 1];
      if (raDir == 0) {
        d[ra1][ra2] = d[ra1 + 1][ra2];
        d[ra1][ra2 + 1] = d[ra1 + 1][ra2 + 1];
      }
      if (raDir == 1) {
        d[ra1][ra2] = d[ra1][ra2 + 1];
        d[ra1 + 1][ra2] = d[ra1 + 1][ra2 + 1];
      }
      if (raDir == 2) {
        d[ra1 + 1][ra2] = d[ra1][ra2];
        d[ra1 + 1][ra2 + 1] = d[ra1][ra2 + 1];
      }
      if (raDir == 3) {
        d[ra1][ra2 + 1] = d[ra1][ra2];
        d[ra1 + 1][ra2 + 1] = d[ra1 + 1][ra2];
      }
    }
    else if (raMode < hypers.Partition[2]) {
      while (true) {
        ra1 = rand_xorshift() % n + 1;
        ra2 = rand_xorshift() % n + 1;
        int dir = rand_xorshift() % 4;
        ra3 = ra1 + dx[dir];
        ra4 = ra2 + dy[dir];
        if (d[ra1][ra2] == 0 || d[ra1][ra2] == d[ra3][ra4])continue;
        break;
      }

      keep1 = d[ra1][ra2];
      d[ra1][ra2] = d[ra3][ra4];
    }

    double tmpScore = -INF;
    if (Check()) {
      // スコア計算
      tmpScore = CalcScore();
    }

    // 焼きなまし
    double diffScore = (tmpScore - ansScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > rand_01()) {
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
        for (int i = 1; i < n + 1; ++i) {
          for (int j = 1; j < n + 1; ++j) {
            d[i][j] = keep[i][j];
          }
        }
      }
      else if (raMode < hypers.Partition[1]) {
        d[ra1][ra2] = keep1;
        d[ra1][ra2 + 1] = keep2;
        d[ra1 + 1][ra2] = keep3;
        d[ra1 + 1][ra2 + 1] = keep4;
      }
      else if (raMode < hypers.Partition[2]) {
        d[ra1][ra2] = keep1;
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
  }

  CopyToAns();
}

// 問題を解く関数
ll Solve(int case_num, Hypers hypers)
{
  start_timer();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(case_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(case_num, ofs);

  // 焼きなまし
  SimulatedAnnealing(hypers);

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

int main()
{
  mode = 2;

  Hypers HYPERS;
  HYPERS.StartTemp = 2048.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 12345.0;
  HYPERS.Partition[0] = 100;
  HYPERS.Partition[1] = 200;
  HYPERS.Partition[2] = 300;
  HYPERS.Partition[3] = 400;
  HYPERS.Partition[4] = 500;
  HYPERS.Partition[5] = 600;
  HYPERS.Partition[6] = 700;
  HYPERS.Partition[7] = 800;
  HYPERS.Partition[8] = 900;
  HYPERS.Partition[9] = 1000;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = Solve(i, HYPERS);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << get_elapsed_time() << ", ";
        cout << endl;
      }
    }
  }
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, rand_01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, rand_01() * 20);
      hypers.Partition[0] = rand_xorshift() % 101;

      ll sum = 0;
      for (int i = 0; i < 15; ++i) {
        ll score = Solve(i, hypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition[0]
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
