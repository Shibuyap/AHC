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

struct Point {
  int x;
  int y;
};

const int MAX_N = 30;

const int n = 200;
const int m = 10;

vector<int> init_b[m];
vector<Point> init_c;

struct Problem {
  vector<int> b[m];
  vector<Point> c;
  vector<P> ans;
};

void CopyToBest()
{
}

void CopyToAns()
{
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  rep(i, m) {
    init_b[i].clear();
  }
  init_c.clear();
}

// 入力受け取り
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  rep(i, m) {
    init_b[i].resize(n / m);
  }
  init_c.resize(n);

  // 標準入力する
  if (!ifs.is_open()) {
    int _n, _m;
    cin >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> init_b[i][j];
        init_b[i][j]--;
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
        ifs >> init_b[i][j];
        init_b[i][j]--;
      }
    }
  }

  rep(i, m) {
    rep(j, n / m) {
      init_c[init_b[i][j]].x = i;
      init_c[init_b[i][j]].y = j;
    }
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
int CalcScore(const vector<P>& ans)
{
  vector<int> tmp_b[m];
  vector<Point> tmp_c = init_c;

  int cnt = 0;
  rep(i, m)
  {
    tmp_b[i] = init_b[i];
  }

  int res = 10000;
  rep(i, ans.size())
  {
    int num = ans[i].first;
    int x = tmp_c[num].x;
    int y = tmp_c[num].y;
    int nx = ans[i].second;
    if (nx == -1) {
      tmp_b[x].pop_back();
      cnt++;
    }
    else {
      res--;
      rep(j, tmp_b[x].size() - y)
      {
        int num2 = tmp_b[x][y + j];
        tmp_c[num2].x = nx;
        tmp_c[num2].y = tmp_b[nx].size();
        tmp_b[nx].push_back(num2);
        res--;
      }
      tmp_b[x].resize(y);
    }
  }
  if (cnt != n) return -1;
  return res;
}

// 解答出力
void Output(ofstream& ofs, const vector<P>& ans)
{
  if (mode == 0) {
    // 標準出力
    rep(i, ans.size()) { cout << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
  else {
    // ファイル出力
    rep(i, ans.size()) { ofs << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
}

void MoveOneTurn(Problem& prob, int turn, int from, const vector<int>& idx) {
  int turnY = prob.c[turn].y;

  for (int k = prob.b[from].size() - 1; k > turnY;) {
    int to = idx[k];
    while (idx[k - 1] == to)k--;

    srep(l, k, prob.b[from].size()) {
      int num = prob.b[from][l];
      prob.c[num].x = to;
      prob.c[num].y = prob.b[to].size();
      prob.b[to].push_back(num);
    }

    prob.ans.emplace_back(prob.b[from][k], to);
    prob.b[from].resize(k);

    k--;
  }

  prob.b[from].pop_back();
  prob.ans.emplace_back(turn, -1);
}

void DecideIdx(const Problem& prob, vector<int>& idx, const vector<P> minIs, int from, int k) {
  int bb = prob.b[from][k];
  int id = 1;
  while (id + 1 < m && minIs[id].first < bb)id++;
  idx[k] = minIs[id].second;
}

int PlayOut() {
  return 0;
}

// kotatsugameさんの貪欲
Problem Method2() {
  Problem prob;
  rep(i, m)prob.b[i] = init_b[i];
  prob.c = init_c;

  rep(turn, n) {
    vector<P> minIs(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, prob.b[i].size()) {
        minI = min(minI, prob.b[i][j]);
      }
      minIs[i] = P(minI, i);
    }

    sort(minIs.begin(), minIs.end());

    int from = minIs[0].second;
    int turnY = prob.c[turn].y;

    vector<int> idx(prob.b[from].size(), -1);
    srep(k, turnY + 1, prob.b[from].size()) {
      DecideIdx(prob, idx, minIs, from, k);
    }

    MoveOneTurn(prob, turn, from, idx);
  }

  return prob;
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
  auto problem = Method2();

  // 解答を出力
  Output(ofs, problem.ans);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore(problem.ans);
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
