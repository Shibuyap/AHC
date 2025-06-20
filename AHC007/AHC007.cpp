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

class Unionfind
{
public:
  const int MAX_UF = 1000006;

  int N;                // 頂点数
  vector<int> par;      // 親
  vector<int> rank;     // 木の深さ

  // 初期化
  void init()
  {
    for (int i = 0; i < N; i++) {
      par[i] = i;
      rank[i] = 0;
    }
  }

  Unionfind(int n)
  {
    N = n;
    par.resize(N);
    rank.resize(N);

    init();
  }

  Unionfind()
  {
    Unionfind(MAX_UF);
  }

  // 木の根を求める
  int find(int x)
  {
    if (par[x] == x) {
      return x;
    }
    else {
      return par[x] = find(par[x]);
    }
  }

  // xとyの属する集合を併合
  void unite(int x, int y)
  {
    x = find(x);
    y = find(y);
    if (x == y) { return; }

    if (rank[x] < rank[y]) {
      par[x] = y;
    }
    else {
      par[y] = x;
      if (rank[x] == rank[y]) rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool same(int x, int y)
  {
    return find(x) == find(y);
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

static double rand_range(double l, double r)
{
  return l + (r - l) * Rand01();
}

void shuffle(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
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

std::chrono::steady_clock::time_point start_time;

void reset_time()
{
  start_time = std::chrono::steady_clock::now();
}

double get_time()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time;
  return elapsed.count();
}

double min_ratio = 1.1;
double max_ratio = 2.9;
int repeat = 50;

const int n = 400;
const int m = 1995;
int x[n], y[n];
int u[m], v[m];
int distances[m];
vector<vector<P>> rand_dist;

// 複数ケース回すときに内部状態を初期値に戻す
void setup()
{
  rand_dist.clear();
}

int get_dist(int i, int j)
{
  return round(sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j])));
}

// 入力ファイルストリームオープン
void open_ifs(int pn, ifstream& ifs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << pn << ".txt";
    ifs.open(oss.str());
  }
}

// 入力受け取り
void input(int pn, ifstream& ifs)
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
    distances[i] = get_dist(u[i], v[i]);
  }
}

// 出力ファイルストリームオープン
void open_ofs(int pn, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << pn << ".txt";
    ofs.open(oss.str());
  }
}

// スコア計算
ll calc_score()
{
  ll res = 0;
  return res;
}

// 解答出力
void output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
  }
  else {
    // ファイル出力
  }
}

void prepare()
{
  for (int _ = 0; _ < repeat; ++_) {
    vector<P> tmp(m);
    for (int i = 0; i < m; ++i) {
      tmp[i] = make_pair(round(rand_range(min_ratio, max_ratio) * distances[i]), i);
    }
    sort(tmp.begin(), tmp.end());
    rand_dist.emplace_back(tmp);
  }
}

int Kruskal(Unionfind uf, int start, int cn)
{
  // ここの時点ではu[start]とv[start]は繋がっていない

  for (auto p : rand_dist[cn]) {
    int cost = p.first;
    int idx = p.second;
    if (start < idx) {
      uf.unite(u[idx], v[idx]);
      if (uf.same(u[start], v[start])) {
        return cost;
      }
    }
  }

  // 辺startを採用しないと連結にならない
  return -1;
}

bool monte_carlo(Unionfind uf, int start)
{
  // ここに来た時点でu[start]とv[start]は繋がっていない

  int sum = 0;
  for (int i = 0; i < repeat; ++i) {
    sum += Kruskal(uf, start, i);
    if (i == 0 && sum == -1) {
      // 辺startを採用しないと連結にならない
      return true;
    }
  }
  return repeat * distances[start] < sum;
}

// ナイーブな解法
void Method1(ifstream& ifs, ofstream& ofs)
{
  prepare();

  Unionfind uf(n);

  for (int i = 0; i < m; ++i) {
    if (mode == 0) {
      cin >> distances[i];
    }
    else {
      ifs >> distances[i];
    }

    if (!uf.same(u[i], v[i]) && monte_carlo(uf, i)) {
      uf.unite(u[i], v[i]);
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

ll Solve(int pn)
{
  reset_time();

  setup();

  ifstream ifs;
  open_ifs(pn, ifs);

  input(pn, ifs);

  ofstream ofs;
  open_ofs(pn, ofs);

  Method1(ifs, ofs);

  output(ofs);

  if (ifs.is_open()) {
    ifs.close();
  }

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calc_score();
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
    for (int i = 0; i < 100; ++i) {
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
