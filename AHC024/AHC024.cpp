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


static double Rand01()
{
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}


static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
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

std::chrono::steady_clock::time_point startTimeClock;

void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTimeClock;
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
  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
      best_d[i][j] = d[i][j];
    }
  }
}

void CopyToAns()
{
  ansScore = best_ansScore;
  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
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
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn, mm;
    cin >> nn >> mm;
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1)
      {
        cin >> c[i][j];
      }
    }
  }
  else {
    // ファイル入力
    int nn, mm;
    ifs >> nn >> mm;
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1)
      {
        ifs >> c[i][j];
      }
    }
  }

  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
      d[i][j] = c[i][j];
    }
  }

  for (int i = 0; i < (m + 1); ++i)
  {
    for (int j = 0; j < (m + 1); ++j)
    {
      g[i][j] = 0;
    }
  }
  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 1); ++j)
    {
      g[c[i][j]][c[i][j + 1]] = 1;
      g[c[i][j + 1]][c[i][j]] = 1;
    }
  }
  for (int i = 0; i < (n + 1); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
      g[c[i][j]][c[i + 1][j]] = 1;
      g[c[i + 1][j]][c[i][j]] = 1;
    }
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
int CalcScore()
{
  int res = 1;
  srep(i, 1, n + 1)
  {
    srep(j, 1, n + 1)
    {
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
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1)
      {
        cout << d[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1)
      {
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
  for (int i = 0; i < (m + 1); ++i)
  {
    for (int j = 0; j < (m + 1); ++j)
    {
      checkG[i][j] = 0;
    }
  }
  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 1); ++j)
    {
      if (d[i][j] == d[i][j + 1])continue;
      if (g[d[i][j]][d[i][j + 1]] == 0) return false;
      if (g[d[i][j + 1]][d[i][j]] == 0) return false;
      checkG[d[i][j]][d[i][j + 1]] = 1;
      checkG[d[i][j + 1]][d[i][j]] = 1;
    }
  }
  for (int i = 0; i < (n + 1); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
      if (d[i][j] == d[i + 1][j])continue;
      if (g[d[i][j]][d[i + 1][j]] == 0) return false;
      if (g[d[i + 1][j]][d[i][j]] == 0) return false;
      checkG[d[i][j]][d[i + 1][j]] = 1;
      checkG[d[i + 1][j]][d[i][j]] = 1;
    }
  }

  for (int i = 0; i < (m + 1); ++i)
  {
    srep(j, i + 1, m + 1)
    {
      if (checkG[i][j] != g[i][j])return false;
    }
  }

  for (int i = 0; i < (n + 2); ++i)
  {
    for (int j = 0; j < (n + 2); ++j)
    {
      checkVisited[i][j] = 0;
    }
  }
  for (int i = 0; i < (m + 1); ++i)
  {
    checkVisited2[i] = 0;
  }
  ClearQueue();
  srep(i, 1, n + 1)
  {
    srep(j, 1, n + 1)
    {
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
        for (int k = 0; k < (4); ++k)
        {
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
  srep(i, 1, n + 1)
  {
    srep(j, 1, n + 1)
    {
      keep[i][j] = d[i][j];
    }
  }
}

void SimulatedAnnealing(Hypers hypers)
{
  ansScore = CalcScore();

  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
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
    int raMode = Rand() % hypers.Partition[2];
    int ra1, ra2, ra3, ra4, ra5, raDir;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      int raDir = Rand() % 2;
      int ra1 = Rand() % n + 1;
      int ok = 0;
      while (ok == 0) {
        raDir = Rand() % 2;
        ra1 = Rand() % n + 1;
        if (raDir == 0) {
          srep(j, 1, n + 1)
          {
            if (d[ra1][j] != 0)ok = 1;
          }
        }
        else {
          srep(i, 1, n + 1)
          {
            if (d[i][ra1] != 0)ok = 1;
          }
        }
      }

      if (raDir == 0) {
        srep(j, 1, n + 1)
        {
          if (g[d[ra1 - 1][j]][d[ra1 + 1][j]] == 0) {
            ok = 0;
          }
        }
      }
      else {
        srep(i, 1, n + 1)
        {
          if (g[d[i][ra1 - 1]][d[i][ra1 + 1]] == 0) {
            ok = 0;
          }
        }
      }

      if (ok == 0)continue;
      KeepD();

      if (raDir == 0) {
        srep(i, ra1, n + 1)
        {
          srep(j, 1, n + 1)
          {
            d[i][j] = d[i + 1][j];
          }
        }
      }
      else {
        srep(j, ra1, n + 1)
        {
          srep(i, 1, n + 1)
          {
            d[i][j] = d[i][j + 1];
          }
        }
      }
    }
    else if (raMode < hypers.Partition[1]) {

      while (true) {
        ra1 = Rand() % (n - 1) + 1;
        ra2 = Rand() % (n - 1) + 1;
        raDir = Rand() % 4;
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
        ra1 = Rand() % n + 1;
        ra2 = Rand() % n + 1;
        int dir = Rand() % 4;
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
        srep(i, 1, n + 1)
        {
          srep(j, 1, n + 1)
          {
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
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

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
    srep(i, 0, 15)
    {
      ll score = Solve(i, HYPERS);
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
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition[0] = Rand() % 101;

      ll sum = 0;
      srep(i, 0, 15)
      {
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
