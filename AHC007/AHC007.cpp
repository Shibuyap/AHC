#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

class UnionFind
{
public:
  const int MAX_UF = 1000006;

  int N;                // 頂点数
  vector<int> Par;      // 親
  vector<int> Rank;     // 木の深さ

  // 初期化
  void Initialize()
  {
    for (int i = 0; i < N; i++) {
      Par[i] = i;
      Rank[i] = 0;
    }
  }

  UnionFind(int n)
  {
    N = n;
    Par.resize(N);
    Rank.resize(N);

    Initialize();
  }

  UnionFind()
  {
    UnionFind(MAX_UF);
  }

  // 木の根を求める
  int Find(int x)
  {
    if (Par[x] == x) {
      return x;
    }
    else {
      return Par[x] = Find(Par[x]);
    }
  }

  // xとyの属する集合を併合
  void Unite(int x, int y)
  {
    x = Find(x);
    y = Find(y);
    if (x == y) return;

    if (Rank[x] < Rank[y]) {
      Par[x] = y;
    }
    else {
      Par[y] = x;
      if (Rank[x] == Rank[y]) Rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool IsSame(int x, int y)
  {
    return Find(x) == Find(y);
  }
};

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

double minCostRatio = 1.1;
double maxCostRatio = 2.9;
int repeat = 50;

const int n = 400;
const int m = 1995;
int x[n], y[n];
int u[m], v[m];
int distances[m];
vector<vector<P>> randomDistances;

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  randomDistances.clear();
}

int GetDistance(int i, int j)
{
  return round(sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j])));
}

// 入力ファイルストリームオープン
void OpenIfs(int problemNum, ifstream& ifs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
    ifs.open(oss.str());
  }
}

// 入力受け取り
void Input(int problemNum, ifstream& ifs)
{
  if (!ifs.is_open()) {
    // 標準入力
    for (int i = 0; i < n; ++i) cin >> x[i] >> y[i];
    for (int i = 0; i < m; ++i) cin >> u[i] >> v[i];
  }
  else {
    // ファイル入力
    for (int i = 0; i < n; ++i) ifs >> x[i] >> y[i];
    for (int i = 0; i < m; ++i) ifs >> u[i] >> v[i];
  }

  for (int i = 0; i < m; ++i) {
    distances[i] = GetDistance(u[i], v[i]);
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
ll CalcScore()
{
  ll res = 0;
  return res;
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
  }
  else {
    // ファイル出力
  }
}

void Prepare()
{
  for (int _ = 0; _ < repeat; ++_) {
    vector<P> tmpVec(m);
    for (int i = 0; i < m; ++i) {
      tmpVec[i] = make_pair(round(RandRange(minCostRatio, maxCostRatio) * distances[i]), i);
    }
    sort(tmpVec.begin(), tmpVec.end());
    randomDistances.emplace_back(tmpVec);
  }
}

int Kruskal(UnionFind uf, int start, int caseNumber)
{
  // ここの時点ではu[start]とv[start]は繋がっていない

  for (auto p : randomDistances[caseNumber]) {
    int cost = p.first;
    int idx = p.second;
    if (start < idx) {
      uf.Unite(u[idx], v[idx]);
      if (uf.IsSame(u[start], v[start])) {
        return cost;
      }
    }
  }

  // 辺startを採用しないと連結にならない
  return -1;
}

bool MonteCarlo(UnionFind uf, int start)
{
  // ここに来た時点でu[start]とv[start]は繋がっていない

  int sum_costs = 0;
  for (int i = 0; i < repeat; ++i) {
    sum_costs += Kruskal(uf, start, i);
    if (i == 0 && sum_costs == -1) {
      // 辺startを採用しないと連結にならない
      return true;
    }
  }
  return repeat * distances[start] < sum_costs;
}

// ナイーブな解法
void Method1(ifstream& ifs, ofstream& ofs)
{
  Prepare();

  UnionFind uf(n);

  for (int i = 0; i < m; ++i) {
    if (mode == 0) {
      cin >> distances[i];
    }
    else {
      ifs >> distances[i];
    }

    if (!uf.IsSame(u[i], v[i]) && MonteCarlo(uf, i)) {
      uf.Unite(u[i], v[i]);
      if (mode == 0) {
        cout << 1 << endl << flush;
      }
      else {
        ofs << 1 << endl;
      }
    }
    else {
      if (mode == 0) {
        cout << 0 << endl << flush;
      }
      else {
        ofs << 0 << endl;
      }
    }
  }
}

ll Solve(int probNum)
{
  ResetTime();

  SetUp();

  ifstream ifs;
  OpenIfs(probNum, ifs);

  Input(probNum, ifs);

  ofstream ofs;
  OpenOfs(probNum, ofs);

  Method1(ifs, ofs);

  Output(ofs);

  if (ifs.is_open()) {
    ifs.close();
  }

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

  if (mode == 0) {
    Solve(0);
  }
  else {
    ll sum = 0;
    for (int i = 0; i < 100; ++i)
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
