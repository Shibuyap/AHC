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

using namespace std;
typedef long long int ll;
typedef pair<ll, ll> P;
typedef pair<P, P> PP;

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


  static double Rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      int tmp = data[i];
      data[i] = data[j];
      data[j] = tmp;
    }
  }
}  // namespace

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;
const int MA = 1000000000;

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

const int n = 1000;
const int MAX_M = 5000;

ll a[n], b[n];
ll L;

ll ans[5500][4];
int ansCount;
ll ansCost;

ll best_ans[5500][4];
int best_ansCount;
ll best_ansCost;

void CopyToBest()
{
  best_ansCount = ansCount;
  best_ansCost = ansCost;

  for (int i = 0; i < ansCount; ++i) {
    for (int j = 0; j < 4; ++j) {
      best_ans[i][j] = ans[i][j];
    }
  }
}

void CopyToAns()
{
  ansCount = best_ansCount;
  ansCost = best_ansCost;

  for (int i = 0; i < ansCount; ++i) {
    for (int j = 0; j < 4; ++j) {
      ans[i][j] = best_ans[i][j];
    }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ansCount = 0;
  ansCost = 0;
}

// 入力受け取り
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    int nnn;
    cin >> nnn;
    for (int i = 0; i < n; ++i) {
      cin >> a[i] >> b[i];
    }
  }
  // ファイル入力する
  else {
    int nnn;
    ifs >> nnn;
    for (int i = 0; i < n; ++i) {
      ifs >> a[i] >> b[i];
    }
  }

  L = 0;
  for (int i = 0; i < n; ++i) {
    L = max(L, a[i]);
    L = max(L, b[i]);
  }

  vector<P> vp;
  for (int i = 0; i < n; ++i) {
    vp.push_back(P(a[i], b[i]));
  }
  sort(vp.begin(), vp.end());
  for (int i = 0; i < n; ++i) {
    a[i] = vp[i].first;
    b[i] = vp[i].second;
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int case_num, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// スコア計算
ll CalcScore()
{
  ll result = round(1.0 * n * L / (ansCost + 1) * 1000000.0);
  return result;
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << ansCount << endl;
    for (int i = 0; i < ansCount; ++i) {
      for (int j = 0; j < 4; ++j) {
        cout << ans[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    ofs << ansCount << endl;
    for (int i = 0; i < ansCount; ++i) {
      for (int j = 0; j < 4; ++j) {
        ofs << ans[i][j] << ' ';
      }
      ofs << endl;
    }
  }
}

void AddOneStepToAns(int x1, int y1, int x2, int y2)
{
  ans[ansCount][0] = x1;
  ans[ansCount][1] = y1;
  ans[ansCount][2] = x2;
  ans[ansCount][3] = y2;
  ansCount++;
  ansCost += (x2 - x1) + (y2 - y1);
}

// 23位解法
void Method23()
{
  vector<P> points;
  vector<int> used;
  priority_queue<PP> que;
  for (int i = 0; i < n; ++i) {
    points.emplace_back(a[i], b[i]);
    used.push_back(0);
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      PP edge;
      edge.first.first = points[i].first + points[i].second + points[j].first + points[j].second
        - abs(points[i].first - points[j].first) - abs(points[i].second - points[j].second);
      edge.second.first = i;
      edge.second.second = j;
      que.push(edge);
    }
  }

  ansCount = 0;
  ansCost = 0;

  while (que.size()) {
    PP edge = que.top();
    que.pop();

    int idx1 = edge.second.first;
    int idx2 = edge.second.second;
    if (used[idx1] || used[idx2]) { continue; }

    int newX = min(points[idx1].first, points[idx2].first);
    int newY = min(points[idx1].second, points[idx2].second);
    AddOneStepToAns(newX, newY, points[idx1].first, points[idx1].second);
    AddOneStepToAns(newX, newY, points[idx2].first, points[idx2].second);

    used[idx1] = 1;
    used[idx2] = 1;

    for (int i = 0; i < points.size(); ++i) {
      if (used[i])continue;

      PP newEdge;
      newEdge.first.first = points[i].first + points[i].second + newX + newY
        - abs(points[i].first - newX) - abs(points[i].second - newY);
      newEdge.second.first = i;
      newEdge.second.second = points.size();
      que.push(newEdge);
    }

    points.emplace_back(newX, newY);
    used.push_back(0);
  }

  for (int i = 0; i < ansCount / 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      swap(ans[i][j], ans[ansCount - 1 - i][j]);
    }
  }
}

ll Solve(int case_num)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(case_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(case_num, ofs);

  Method23();

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
  mode = 0;

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
