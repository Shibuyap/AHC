﻿#include "Haipara.h"

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
typedef pair<int, int> P;

const ll INF = 1001001001001001001;

const ll D1 = 10LL;
const ll D2 = 100LL;
const ll D3 = 1000LL;
const ll D4 = 10000LL;
const ll D5 = 100000LL;
const ll D6 = 1000000LL;
const ll D7 = 10000000LL;
const ll D8 = 100000000LL;
const ll D9 = 1000000000LL;
const ll D10 = 10000000000LL;
const ll D11 = 100000000000LL;
const ll D12 = 1000000000000LL;
const ll D13 = 10000000000000LL;
const ll D14 = 100000000000000LL;
const ll D15 = 1000000000000000LL;
const ll D16 = 10000000000000000LL;
const ll D17 = 100000000000000000LL;
const ll D18 = 1000000000000000000LL;

namespace /* 乱数ライブラリ */
{
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
}  // namespace

std::random_device seed_gen;
std::default_random_engine engine(seed_gen());
std::exponential_distribution<> dist(1e-5);
std::mt19937 engine_mt19937(seed_gen());

class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};
Timer timer;
double elapsed;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.3;
int mode;

const int MAX_Q = 3232;
int N, D, Q;
int NN, DD, QQ;
vector<int> l[MAX_Q], r[MAX_Q];
int answers[MAX_Q][110];
string comments[MAX_Q];
char C[MAX_Q];
int ans[110];

int real_ans[110];
ll real_minScore;
int real_max_D;

int hikaku[110][110];

vector<int> pseudoItems;

// ローカルテスト用入力
ll W[110];
int karusa[110];

// ハイパラ
std::vector<int> haipara[14][40];
std::vector<long long int> haipara2[14][40];

// N * Q * D = 14 * 40 * 24 = 6720
void CalcNNDDQQ()
{
  NN = (N - 30) / 5;
  NN = min(NN, 13);  // 14
  DD = D - 2;        // 24
  vector<int> qv;
  for (int i = 0; i < 40; ++i) {
    double si = 1.0 + 0.1 * (i + 1);
    qv.push_back(round(pow(2, si) * N));
    if (i == 39) {
      qv[i] *= 10;
    }
  }
  QQ = lower_bound(qv.begin(), qv.end(), Q) - qv.begin();  // 40
}

void GenerateLocalItems()
{
  for (int i = 0; i < N; ++i) {
    double num;
    while (true) {
      num = dist(engine);
      if (num < 100000.0 * N / D) {
        break;
      }
    }
    W[i] = round(num);
    W[i] = max(W[i], 1LL);
  }
}

void GenerateCase(int _n = -1, int _d = -1, int _q = -1)
{
  N = rand_xorshift() % 71 + 30;
  if (_n != -1) N = _n;
  D = rand_xorshift() % (N / 4 - 1) + 2;
  if (_d != -1) D = _d;
  Q = round(pow(2, rand_01() * 4.0 + 1.0) * N);
  if (_q != -1) Q = _q;
  GenerateLocalItems();
}

void GenerateNNDDQQ(int _nn = -1, int _qq = -1, int _dd = -1)
{
  NN = rand_xorshift() % 14;
  if (_nn != -1) {
    NN = _nn;
  }
  QQ = rand_xorshift() % 40;
  if (_qq != -1) {
    QQ = _qq;
  }
  DD = rand_xorshift() % haipara[NN][QQ].size();
  if (_dd != -1) {
    DD = _dd;
  }
}

void GeneratecaseFromNNDDQQ()
{
  if (NN == 13) {
    N = rand_xorshift() % 6 + 30 + NN * 5;
  }
  else {
    N = rand_xorshift() % 5 + 30 + NN * 5;
  }
  double dou = rand_01() * 0.1;
  Q = round(pow(2, 1.0 + 0.1 * QQ + dou) * N);
  D = DD + 2;
  GenerateLocalItems();
}

vector<int> pseudoItemsAll[101][26];
void GeneratePseudoItems()
{
  if (!pseudoItemsAll[N][D].empty()) {
    pseudoItems = pseudoItemsAll[N][D];
  }
  else {
    double pseudoSum[110] = {};
    vector<double> nums(110);
    for (int loop = 0; loop < 100; ++loop) {
      for (int i = 0; i < N; ++i) {
        while (true) {
          nums[i] = dist(engine);
          if (nums[i] < 100000.0 * N / D) {
            break;
          }
        }
        nums[i] = max(nums[i], 1.0);
      }
      sort(nums.begin(), nums.begin() + N);
      for (int i = 0; i < N; ++i) { pseudoSum[i] += nums[i]; }
    }
    pseudoItems.clear();
    for (int i = 0; i < N; ++i) { pseudoItems.push_back(round(pseudoSum[i] / 1000)); }
    pseudoItemsAll[N][D] = pseudoItems;
  }
}

// 入力受け取り
void Input(int case_num)
{
  if (mode <= 2) {
    if (mode == 0) {
      // 標準入力する
      cin >> N >> D >> Q;
    }
    else if (mode <= 2) {
      // ファイル入力する
      std::ostringstream oss;
      oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
      ifstream ifs(oss.str());

      ifs >> N >> D >> Q;
      for (int i = 0; i < N; ++i) { ifs >> W[i]; }
    }
    CalcNNDDQQ();
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int case_num, ofstream& ofs)
{
  if (mode < 1000000) {
    if (mode != 0) {
      std::ostringstream oss;
      oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
      ofs.open(oss.str());
    }
  }
}

// Forward declarations
char Query(int& turn);
void DummyQuery(int& countQ);
void Move1(int& countQ, int cutLine = 999);
void Swap1(int& countQ, int diffLine = 999);
int binarySearchGroupPosition(int gId, const vector<int>& groups, int initialLeft, int initialRight, int& countQ);
int binarySearchGroupPosition(int gId, const vector<int>& groups, int initialLeft, int initialRight, int& countQ);

// ========== Additional Common Helper Functions ==========

// Check if time limit is approaching
bool isTimeLimitApproaching(double threshold = 0.0)
{
  return elapsed > TL + threshold;
}

// Check if time limit fraction is exceeded
bool isTimeLimitFractionExceeded(double fraction)
{
  return elapsed > TL * fraction;
}

// Check if we're near query limit
bool isNearQueryLimit(int countQ, int margin = 1)
{
  return countQ >= Q - margin;
}

// Random percentage selection
int getRandomPercentage()
{
  return rand_xorshift() % 100;
}

// Initialize answer array with default pattern
void initializeAnswerWithDefault()
{
  for (int i = 0; i < N; ++i) {
    ans[i] = i % D;
  }
}

// Complete remaining queries with dummy
void completeWithDummyQueries(int& countQ)
{
  DummyQuery(countQ);
}

// Adjust hiritu based on time
int adjustHirituByTime(int baseHiritu, double timeFraction = 1.0 / 3.0)
{
  if (baseHiritu >= 100 && isTimeLimitFractionExceeded(timeFraction)) {
    return 90;
  }
  return baseHiritu;
}

// Standard Move/Swap selection pattern
void performStandardMoveSwap(int& countQ, int hiritu)
{
  if (isNearQueryLimit(countQ)) {
    if (hiritu < 10) {
      return;
    }
    return;
  }

  int qu = getRandomPercentage();
  if (qu < hiritu) {
    Move1(countQ);
  }
  else {
    if (isNearQueryLimit(countQ)) {
      return;
    }
    Swap1(countQ);
  }
}

// ========== End of Additional Common Helper Functions ==========

// Common refactored functions
void initializeAnsArray()
{
  initializeAnswerWithDefault();
}

int countItemsInGroup(int groupId)
{
  int cnt = 0;
  for (int j = 0; j < N; ++j) {
    if (ans[j] == groupId) {
      cnt++;
    }
  }
  return cnt;
}

void collectItemsFromGroup(int groupId, vector<int>& items)
{
  items.clear();
  for (int j = 0; j < N; ++j) {
    if (ans[j] == groupId) {
      items.push_back(j);
    }
  }
}

void copyAnswersArray(int turn)
{
  for (int i = 0; i < N; ++i) {
    answers[turn][i] = ans[i];
  }
}

void initializeGroups(vector<int>& groups)
{
  for (int i = 0; i < D; ++i) {
    groups.push_back(i);
  }
}

void initializeItems(vector<int>& items)
{
  for (int i = 0; i < N; ++i) {
    items.push_back(i);
  }
}

void updateKarusaArray(const vector<int>& items)
{
  for (int i = 0; i < N; ++i) {
    karusa[items[i]] = i;
  }
}

void populateAnsItems(vector<int> ansItems[])
{
  for (int i = 0; i < N; ++i) {
    ansItems[ans[i]].push_back(i);
  }
}

void moveGroupToPosition(vector<int>& groups, int from, int to)
{
  while (from > to) {
    swap(groups[from - 1], groups[from]);
    from--;
  }
}

void assignHeaviestToGroups(const vector<int>& items, int startIdx)
{
  for (int i = 0; i < D; ++i) {
    int id = items[startIdx - i];
    ans[id] = i;
  }
}

void assignToLightestGroup(const vector<int>& items, int itemIdx, int& countQ, vector<int>& groups)
{
  int id = items[itemIdx];
  int gId = groups[D - 1];
  ans[id] = gId;
  int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
  moveGroupToPosition(groups, D - 1, left);
}

template<typename T>
char compareQuerySums(const vector<int>& left_items, const vector<int>& right_items, const T& values)
{
  ll sumL = 0, sumR = 0;
  for (int j : left_items) sumL += values[j];
  for (int j : right_items) sumR += values[j];

  if (sumL < sumR) return '<';
  else if (sumL == sumR) return '=';
  else return '>';
}

// Binary search to find position of gId in groups
int binarySearchGroupPosition(int gId, const vector<int>& groups, int initialLeft, int initialRight, int& countQ)
{
  int left = initialLeft;
  int right = initialRight;
  while (left < right) {
    int mid = (left + right) / 2;
    for (int j = 0; j < N; ++j) {
      if (ans[j] == gId) {
        l[countQ].push_back(j);
      }
      else if (ans[j] == groups[mid]) {
        r[countQ].push_back(j);
      }
    }
    char c = Query(countQ);
    if (c == '=') {
      left = right = mid;
    }
    else if (c == '>') {
      right = mid;
    }
    else {
      left = mid + 1;
    }
  }
  return left;
}

// Build query groups for comparison
void buildQueryGroups(int& countQ, int x, int y, int excludeIndex = -1)
{
  for (int j = 0; j < N; ++j) {
    if (ans[j] == x && j != excludeIndex) {
      l[countQ].push_back(j);
    }
    if (ans[j] == y && (excludeIndex == -1 || j != excludeIndex)) {
      r[countQ].push_back(j);
    }
  }
}

// Check if should skip based on cutline
bool shouldSkipBasedOnCutline(int z, int cutLine)
{
  if (cutLine >= 100) return false;

  int win = 0, lose = 0;
  for (int j = 0; j < N; ++j) {
    if (hikaku[z][j] == 1) win++;
    else if (hikaku[z][j] == -1) lose++;
  }
  return (win + lose >= cutLine && win >= lose);
}

int selectRandomFromGroup(int groupId)
{
  vector<int> groupItems;
  for (int j = 0; j < N; ++j) {
    if (ans[j] == groupId) {
      groupItems.push_back(j);
    }
  }
  return groupItems.empty() ? -1 : groupItems[rand_xorshift() % groupItems.size()];
}

int selectGroupWithMinSize(int minSize)
{
  int x = rand_xorshift() % D;
  int loop = 0;
  while (true) {
    loop++;
    if (loop >= 100) return -1;
    x = rand_xorshift() % D;
    int cnt = countItemsInGroup(x);
    if (cnt >= minSize) return x;
  }
}

bool ErrorCheck()
{
  for (int i = 0; i < Q; ++i) {
    if (l[i].empty() || r[i].empty()) {
      cerr << "NG ErrorCheck1 : " << i << endl;
      return false;
    }
    if (l[i].size() + r[i].size() > N) {
      cerr << "NG ErrorCheck2 : " << i << endl;
      return false;
    }
  }

  int cnt[110] = {};
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < N; ++j) cnt[j] = 0;
    for (auto x : l[i]) cnt[x]++;
    for (auto y : r[i]) cnt[y]++;
    for (int j = 0; j < N; ++j) {
      if (cnt[j] >= 2) {
        cerr << "NG ErrorCheck3 : " << i << endl;
        return false;
      }
    }
  }

  return true;
}

// スコア計算
ll CalcScore()
{
  double t[30] = {};
  for (int i = 0; i < N; ++i) { t[ans[i]] += W[i]; }
  double t_ave = 0;
  for (int i = 0; i < D; ++i) { t_ave += t[i]; }
  t_ave /= D;
  double V = 0;
  for (int i = 0; i < D; ++i) { V += (t[i] - t_ave) * (t[i] - t_ave); }
  V /= D;
  return 1 + round(100.0 * sqrt(V));
}

char LocalQuery(int turn)
{
  return compareQuerySums(l[turn], r[turn], W);
}

char PseudoItemsQuery(int turn)
{
  vector<int> left_vals, right_vals;
  for (int j : l[turn]) left_vals.push_back(karusa[j]);
  for (int j : r[turn]) right_vals.push_back(karusa[j]);
  return compareQuerySums(left_vals, right_vals, pseudoItems);
}

char PseudoItemsQueryGroup(vector<int> vl, vector<int> vr)
{
  return compareQuerySums(vl, vr, pseudoItems);
}

map<pair<vector<int>, vector<int>>, char> calledMap;
int queryCount = 0;
char Query(int& turn)
{
  queryCount++;
  bool isUse = false;
  if (mode < 1000000) {
    if (queryCount % 100 == 0) {
      elapsed = timer.get_elapsed_time();
    }
    if (elapsed < TL) {
      isUse = true;
    }
  }
  else {
    if (queryCount < 10000) {
      isUse = true;
    }
  }

  if (isUse) {
    if (calledMap.find(make_pair(l[turn], r[turn])) != calledMap.end()) {
      char c = calledMap[make_pair(l[turn], r[turn])];
      l[turn].clear();
      r[turn].clear();
      return c;
    }
  }

  if (turn >= Q) {
    // WA対策
    l[turn].clear();
    r[turn].clear();
    return '>';
  }

  if (elapsed > TL) {
    if (l[turn].empty()) {
      l[turn].clear();
      r[turn].clear();
      l[turn].push_back(0);
      r[turn].push_back(1);
      if (mode == 0) {
        cout << l[turn].size() << ' ' << r[turn].size();
        for (int j = 0; j < l[turn].size(); ++j) { cout << ' ' << l[turn][j]; }
        for (int j = 0; j < r[turn].size(); ++j) { cout << ' ' << r[turn][j]; }
        cout << endl;
        fflush(stdout);
        cin >> C[turn];
      }
      else {
        C[turn] = LocalQuery(turn);
      }
      for (int i = 0; i < N; ++i) { answers[turn][i] = ans[i]; }

      turn++;
      if (turn > 0) {
        comments[turn] = comments[turn - 1];
      }
      return '<';
    }
    if (r[turn].empty()) {
      l[turn].clear();
      r[turn].clear();
      l[turn].push_back(0);
      r[turn].push_back(1);
      if (mode == 0) {
        cout << l[turn].size() << ' ' << r[turn].size();
        for (int j = 0; j < l[turn].size(); ++j) { cout << ' ' << l[turn][j]; }
        for (int j = 0; j < r[turn].size(); ++j) { cout << ' ' << r[turn][j]; }
        cout << endl;
        fflush(stdout);
        cin >> C[turn];
      }
      else {
        C[turn] = LocalQuery(turn);
      }
      for (int i = 0; i < N; ++i) { answers[turn][i] = ans[i]; }

      turn++;
      if (turn > 0) {
        comments[turn] = comments[turn - 1];
      }
      return '>';
    }
  }

  if (l[turn].empty()) {
    l[turn].clear();
    r[turn].clear();
    return '<';
  }
  if (r[turn].empty()) {
    l[turn].clear();
    r[turn].clear();
    return '>';
  }

  if (mode == 0) {
    cout << l[turn].size() << ' ' << r[turn].size();
    for (int j = 0; j < l[turn].size(); ++j) { cout << ' ' << l[turn][j]; }
    for (int j = 0; j < r[turn].size(); ++j) { cout << ' ' << r[turn][j]; }
    cout << endl;
    fflush(stdout);
    cin >> C[turn];
  }
  else {
    C[turn] = LocalQuery(turn);
  }
  for (int i = 0; i < N; ++i) { answers[turn][i] = ans[i]; }

  if (true) {
    calledMap[make_pair(l[turn], r[turn])] = C[turn];
    char c2 = '=';
    if (C[turn] == '<') {
      c2 = '>';
    }
    else if (C[turn] == '>') {
      c2 = '<';
    }
    calledMap[make_pair(r[turn], l[turn])] = c2;
  }

  turn++;
  if (turn > 0) {
    comments[turn] = comments[turn - 1];
  }
  return C[turn - 1];
}

char QueryMapDirectly(const vector<int>& vl, const vector<int>& vr)
{
  if (calledMap.find(make_pair(vl, vr)) != calledMap.end()) {
    return calledMap[make_pair(vl, vr)];
  }
  return '?';
}

int myQueue[11000][2];
int query1Count;
char Query1(int& turn, int lhs, int rhs)
{
  query1Count++;
  if (query1Count % 100 == 0) {
    elapsed = timer.get_elapsed_time();
  }
  if (elapsed > TL - 0.1) {
    cerr << "Assert Query1" << endl;
  }
  if (hikaku[lhs][rhs] != -2 && elapsed < TL - 0.1) {
    if (hikaku[lhs][rhs] == -1) {
      return '<';
    }
    if (hikaku[lhs][rhs] == 0) {
      return '=';
    }
    if (hikaku[lhs][rhs] == 1) {
      return '>';
    }
  }
  l[turn].push_back(lhs);
  r[turn].push_back(rhs);
  char c = Query(turn);
  {
    int head = 0;
    int tail = 0;
    myQueue[tail][0] = lhs;
    myQueue[tail][1] = rhs;
    tail++;
    if (c == '<') {
      hikaku[lhs][rhs] = -1;
      hikaku[rhs][lhs] = 1;
    }
    else if (c == '=') {
      hikaku[lhs][rhs] = 0;
      hikaku[rhs][lhs] = 0;
    }
    else {
      hikaku[lhs][rhs] = 1;
      hikaku[rhs][lhs] = -1;
    }
    while (head < tail) {
      int x = myQueue[head][0];
      int y = myQueue[head][1];
      head++;
      int z = hikaku[x][y];
      for (int i = 0; i < N; ++i) {
        if (z == -1) {
          if (hikaku[i][y] == -2) {
            if (hikaku[i][x] == -1 || hikaku[i][x] == 0) {
              hikaku[i][y] = -1;
              hikaku[y][i] = 1;
              myQueue[tail][0] = i;
              myQueue[tail][1] = y;
              tail++;
            }
          }
          else if (hikaku[x][i] == -2) {
            if (hikaku[y][i] == -1 || hikaku[y][i] == 0) {
              hikaku[x][i] = -1;
              hikaku[i][x] = 1;
              myQueue[tail][0] = i;
              myQueue[tail][1] = x;
              tail++;
            }
          }
        }
        else if (z == 0) {
          if (hikaku[i][y] == -2) {
            if (hikaku[i][x] != -2) {
              hikaku[i][y] = hikaku[i][x];
              hikaku[y][i] = hikaku[x][i];
              myQueue[tail][0] = i;
              myQueue[tail][1] = y;
              tail++;
            }
          }
          else if (hikaku[x][i] == -2) {
            if (hikaku[i][y] != -2) {
              hikaku[i][x] = hikaku[i][y];
              hikaku[x][i] = hikaku[y][i];
              myQueue[tail][0] = i;
              myQueue[tail][1] = x;
              tail++;
            }
          }
        }
        else {
          if (hikaku[i][y] == -2) {
            if (hikaku[i][x] == 1 || hikaku[i][x] == 0) {
              hikaku[i][y] = 1;
              hikaku[y][i] = -1;
              myQueue[tail][0] = i;
              myQueue[tail][1] = y;
              tail++;
            }
          }
          else if (hikaku[x][i] == -2) {
            if (hikaku[y][i] == 1 || hikaku[y][i] == 0) {
              hikaku[x][i] = 1;
              hikaku[i][x] = -1;
              myQueue[tail][0] = i;
              myQueue[tail][1] = x;
              tail++;
            }
          }
        }
      }
    }
  }
  return c;
}

char QueryGroup(int& turn, int lgId, int rgId)
{
  for (int i = 0; i < N; ++i) {
    if (ans[i] == lgId) {
      l[turn].push_back(i);
    }
    if (ans[i] == rgId) {
      r[turn].push_back(i);
    }
  }
  char c = Query(turn);
  return c;
}

void DummyQuery(int& countQ)
{
  while (countQ < Q) {
    l[countQ].push_back(0);
    r[countQ].push_back(1);
    if (mode == 0) {
      cout << l[countQ].size() << ' ' << r[countQ].size();
      for (int j = 0; j < l[countQ].size(); ++j) { cout << ' ' << l[countQ][j]; }
      for (int j = 0; j < r[countQ].size(); ++j) { cout << ' ' << r[countQ][j]; }
      cout << endl;
      fflush(stdout);
      cin >> C[countQ];
    }
    else {
      C[countQ] = LocalQuery(countQ);
    }
    for (int i = 0; i < N; ++i) { answers[countQ][i] = ans[i]; }

    countQ++;
    if (countQ > 0) {
      comments[countQ] = comments[countQ - 1];
    }
  }
}

int memory_SwapNeighbor1[110];
vector<int> DItems[30];
void ResetMemory()
{
  for (int i = 0; i < 110; ++i) { memory_SwapNeighbor1[i] = 0; }
  for (int i = 0; i < 30; ++i) {
    DItems[i].clear();
  }
}

// Base Move1 implementation
void Move1Base(int& countQ, int srcGroup, int dstGroup, int itemToMove, const string& methodName)
{
  buildQueryGroups(countQ, srcGroup, dstGroup, itemToMove);

  char c = Query(countQ);
  if (c == '>') {
    ans[itemToMove] = dstGroup;
    comments[countQ] += methodName + " ";
    ResetMemory();
  }
}

// 1個移動
void Move1(int& countQ, int cutLine)
{
  int x = selectGroupWithMinSize(2);
  if (x == -1) { return; }

  int y = rand_xorshift() % D;
  while (x == y) {
    y = rand_xorshift() % D;
  }

  int z = selectRandomFromGroup(x);
  if (z == -1) { return; }

  if (shouldSkipBasedOnCutline(z, cutLine)) {
    return;
  }

  Move1Base(countQ, x, y, z, "Move1");
}

// 戻り値：移動先グループ
int Move1_Specify(int& countQ, int groupId)
{
  int x = groupId;
  int y = rand_xorshift() % D;
  while (x == y) {
    y = rand_xorshift() % D;
  }
  vector<int> vx;
  collectItemsFromGroup(x, vx);
  if (vx.size() <= 1) {
    return -1;
  }
  int z = vx[rand_xorshift() % vx.size()];

  buildQueryGroups(countQ, x, y, z);

  char c = Query(countQ);
  if (c == '>') {
    ans[z] = y;
    comments[countQ] += "Move1 ";
    ResetMemory();
    return y;
  }
  return -1;
}

void MoveSmall1(const vector<int>&, int& countQ, int smallLine)
{
  int _N = min(smallLine, N);
  _N = max(1, _N);
  int z = rand_xorshift() % _N;
  int x = ans[z];
  int y = rand_xorshift() % D;
  while (x == y) {
    y = rand_xorshift() % D;
  }

  buildQueryGroups(countQ, x, y, z);

  char c = Query(countQ);
  if (c == '>') {
    ans[z] = y;
    comments[countQ] += "MoveSmall1 ";
    ResetMemory();
  }
}

void Move1_Two(int& countQ, int cutLine = 999)
{
  if (D < 3) {
    Move1(countQ, cutLine);
    return;
  }

  int x = selectGroupWithMinSize(2);
  if (x == -1) { return; }

  int y1 = rand_xorshift() % D;
  while (x == y1) {
    y1 = rand_xorshift() % D;
  }
  int y2 = rand_xorshift() % D;
  while (x == y2 || y1 == y2) {
    y2 = rand_xorshift() % D;
  }

  char c1 = QueryGroup(countQ, y1, y2);
  int y = (c1 == '>') ? y2 : y1;

  int z = selectRandomFromGroup(x);
  if (z == -1) { return; }

  if (shouldSkipBasedOnCutline(z, cutLine)) {
    return;
  }

  Move1Base(countQ, x, y, z, "Move1");
}

void Move1Minimum(int& countQ, int cutLine = 999)
{
  int x = selectGroupWithMinSize(2);
  if (x == -1) { return; }

  int y = rand_xorshift() % D;
  while (x == y) {
    y = rand_xorshift() % D;
  }

  vector<int> vv;
  collectItemsFromGroup(x, vv);
  if (vv.empty()) { return; }

  int z = vv[rand_xorshift() % vv.size()];
  for (int _ = 0; _ < 30; ++_) {
    int zz = vv[rand_xorshift() % vv.size()];
    if (hikaku[zz][z] == 0 || hikaku[zz][z] == -1) {
      z = zz;
    }
  }

  if (shouldSkipBasedOnCutline(z, cutLine)) {
    return;
  }

  Move1Base(countQ, x, y, z, "Move1Minimum");
}

int Move1Combo(int& countQ, int combo, int cutLine = 999)
{
  int x = combo;
  if (x == -1) {
    x = rand_xorshift() % D;
  }
  while (true) {
    int cnt = countItemsInGroup(x);
    if (cnt >= 2) {
      break;
    }
    x = rand_xorshift() % D;
  }
  int y = rand_xorshift() % D;
  while (x == y) {
    y = rand_xorshift() % D;
  }
  vector<int> vv;
  collectItemsFromGroup(x, vv);
  int z = vv[rand_xorshift() % vv.size()];

  if (cutLine < 100) {
    int win = 0;
    int lose = 0;
    for (int j = 0; j < N; ++j) {
      if (hikaku[z][j] == 1) {
        win++;
      }
      else if (hikaku[z][j] == -1) {
        lose++;
      }
    }
    if (win + lose >= cutLine && win >= lose) {
      return -1;
    }
  }

  buildQueryGroups(countQ, x, y, z);

  char c = Query(countQ);
  if (c == '>') {
    ans[z] = y;
    comments[countQ] += "Move1Combo ";
    ResetMemory();
    return x;
  }
  return -1;
}

// Helper to select two different groups with minimum size
bool selectTwoGroupsWithMinSize(int& x, int& y, int minSize = 2)
{
  x = selectGroupWithMinSize(minSize);
  if (x == -1) return false;

  int yloop = 0;
  while (true) {
    yloop++;
    if (yloop == 30) return false;

    y = rand_xorshift() % D;
    if (x == y) { continue; }

    int cnt = countItemsInGroup(y);
    if (cnt >= minSize) return true;
  }
}

// 1個交換
int arr8_2_L[110];
int arr8_2_R[110];
void Swap1(int& countQ, int diffLine)
{
  int x, y;
  if (!selectTwoGroupsWithMinSize(x, y)) { return; }

  vector<int> vx, vy;
  collectItemsFromGroup(x, vx);
  collectItemsFromGroup(y, vy);

  int lid = vx[rand_xorshift() % vx.size()];
  int rid = vy[rand_xorshift() % vy.size()];

  if (diffLine < 100) {
    if (hikaku[lid][rid] == -1) {
      int diff = 0;
      for (int i = 0; i < N; ++i) {
        if (hikaku[lid][i] == -1 && hikaku[i][rid] == -1) {
          diff++;
          if (diff >= diffLine) {
            return;
          }
        }
      }
    }
    else if (hikaku[lid][rid] == 1) {
      int diff = 0;
      for (int i = 0; i < N; ++i) {
        if (hikaku[lid][i] == 1 && hikaku[i][rid] == 1) {
          diff++;
          if (diff >= diffLine) {
            return;
          }
        }
      }
    }
  }

  char c1 = Query1(countQ, lid, rid);
  if (c1 == '=') {
    return;
  }

  for (auto j : vx) {
    if (j != lid) {
      l[countQ].push_back(j);
    }
  }
  for (auto j : vy) {
    if (j != rid) {
      r[countQ].push_back(j);
    }
  }
  char c2 = Query(countQ);
  bool isSwap = false;
  if (c2 == '=') {
    isSwap = true;
  }
  else if (c1 == c2) {
    isSwap = true;
  }

  if (isSwap) {
    swap(ans[lid], ans[rid]);
    comments[countQ] += "Swap1 ";
    ResetMemory();
  }
}

int arrSwapHalf_L[110];
int arrSwapHalf_R[110];
void SwapHalf(int& countQ)
{
  int x = selectGroupWithMinSize(2);
  if (x == -1) { return; }

  int y = rand_xorshift() % D;
  while (true) {
    y = rand_xorshift() % D;
    if (x == y) {
      continue;
    }
    int cnt = countItemsInGroup(y);
    if (cnt >= 2) {
      break;
    }
  }

  vector<int> vx, vy;
  collectItemsFromGroup(x, vx);
  collectItemsFromGroup(y, vy);

  int szL = vx.size();
  int szR = vy.size();
  for (int i = 0; i < szL; i++) {
    arrSwapHalf_L[i] = vx[i];
  }
  for (int i = 0; i < szR; i++) {
    arrSwapHalf_R[i] = vy[i];
  }

  std::shuffle(arrSwapHalf_L, arrSwapHalf_L + szL, engine_mt19937);
  std::shuffle(arrSwapHalf_R, arrSwapHalf_R + szR, engine_mt19937);
  for (int i = 0; i < szL / 2; ++i) { l[countQ].push_back(arrSwapHalf_L[i]); }
  for (int i = 0; i < szR / 2; ++i) { r[countQ].push_back(arrSwapHalf_R[i]); }

  char c1 = Query(countQ);
  if (c1 == '=') {
    for (int i = 0; i < szL / 2; ++i) { ans[arrSwapHalf_L[i]] = y; }
    for (int i = 0; i < szR / 2; ++i) { ans[arrSwapHalf_R[i]] = x; }
    comments[countQ] += "SwapHalf ";
    ResetMemory();
    return;
  }

  for (int i = szL / 2; i < szL; ++i) { l[countQ].push_back(arrSwapHalf_L[i]); }
  for (int i = szR / 2; i < szR; ++i) { r[countQ].push_back(arrSwapHalf_R[i]); }
  char c2 = Query(countQ);
  bool isSwap = false;
  if (c2 == '=') {
    isSwap = true;
  }
  else if (c1 == c2) {
    isSwap = true;
  }

  if (isSwap) {
    for (int i = 0; i < szL / 2; ++i) { ans[arrSwapHalf_L[i]] = y; }
    for (int i = 0; i < szR / 2; ++i) { ans[arrSwapHalf_R[i]] = x; }
    comments[countQ] += "SwapHalf ";
    ResetMemory();
  }
}

void Method11_1(int& countQ)
{
  int x = selectGroupWithMinSize(2);
  if (x == -1) { return; }

  int y = rand_xorshift() % D;
  while (true) {
    y = rand_xorshift() % D;
    if (x == y) {
      continue;
    }
    int cnt = countItemsInGroup(y);
    if (cnt >= 2) {
      break;
    }
  }

  vector<int> vl, vr;
  collectItemsFromGroup(x, vl);
  collectItemsFromGroup(y, vr);
  int lid = vl[rand_xorshift() % vl.size()];
  int rid = vr[rand_xorshift() % vr.size()];
  for (int i = 0; i < 50; ++i) {
    rid = vr[rand_xorshift() % vr.size()];
    if (hikaku[lid][rid] == 1) {
      break;
    }
  }

  char c1 = Query1(countQ, lid, rid);
  if (c1 == '=') {
    return;
  }

  for (auto j : vl) {
    if (j != lid) {
      l[countQ].push_back(j);
    }
  }
  for (auto j : vr) {
    if (j != rid) {
      r[countQ].push_back(j);
    }
  }
  char c2 = Query(countQ);
  bool isSwap = false;
  if (c2 == '=') {
    isSwap = true;
  }
  else if (c1 == c2) {
    isSwap = true;
  }

  if (isSwap) {
    swap(ans[lid], ans[rid]);
    comments[countQ] += "Method11_1 ";
    ResetMemory();
  }
}

bool IsAllSearched_SwapNeighbor1(const vector<int>& items)
{
  vector<int> ansItems[30];
  populateAnsItems(ansItems);
  vector<int> vl, vr;
  for (int num = 0; num < N - 1; ++num) {
    int x = items[num];
    int y = items[num + 1];
    int xg = ans[x];
    int yg = ans[y];
    if (xg == yg) {
      continue;
    }
    if (ansItems[xg].size() <= 1 || ansItems[yg].size() <= 1) {
      continue;
    }
    vl.clear();
    vr.clear();
    for (auto i : ansItems[xg]) {
      if (i != x) {
        vl.push_back(i);
      }
    }
    for (auto i : ansItems[yg]) {
      if (i != y) {
        vr.push_back(i);
      }
    }

    char c = QueryMapDirectly(vl, vr);
    if (c == '?') {
      return false;
    }

    // 更新可能か
    if (c == '<') {
      return false;
    }
  }

  return true;
}

// Helper for SwapNeighbor methods
void buildSwapNeighborQuery(int& countQ, int x, int y, int xg, int yg)
{
  if (DItems[0].empty()) {
    populateAnsItems(DItems);
  }
  for (auto i : DItems[xg]) {
    if (i != x) {
      l[countQ].push_back(i);
    }
  }
  for (auto i : DItems[yg]) {
    if (i != y) {
      r[countQ].push_back(i);
    }
  }
}

// 重さの近いものをスワップ
bool SwapNeighbor1(const vector<int>& items, int& countQ, int _m = -1)
{
  int M = _m;
  if (_m == -1) M = N;
  bool isChange = false;
  int loop = 0;
  while (true) {
    loop++;
    int num = rand_xorshift() % (M - 1);
    for (int _ = 0; _ < 30; ++_) {
      if (memory_SwapNeighbor1[num]) {
        num = rand_xorshift() % (M - 1);
      }
      else {
        break;
      }
    }
    memory_SwapNeighbor1[num] = 1;
    int x = items[num];
    int y = items[num + 1];
    int xg = ans[x];
    int yg = ans[y];
    if (xg == yg) {
      if (loop < 100) {
        continue;
      }
      else {
        while (xg == yg) {
          num = rand_xorshift() % (N - 1);
          x = items[num];
          y = items[num + 1];
          xg = ans[x];
          yg = ans[y];
        }
      }
    }

    buildSwapNeighborQuery(countQ, x, y, xg, yg);
    char c = Query(countQ);

    if (c == '<') {
      swap(ans[x], ans[y]);
      comments[countQ] += "SwapNeighbor1 ";
      ResetMemory();
      isChange = true;
    }
    break;
  }
  return isChange;
}

bool SwapNeighborSmall1(const vector<int>& items, int& countQ, int smallLine)
{
  bool isChange = false;
  while (true) {
    int _N = min(N, smallLine);
    _N = max(2, _N);
    int num = rand_xorshift() % (_N - 1);
    int x = items[num];
    int y = items[num + 1];
    int xg = ans[x];
    int yg = ans[y];
    if (xg == yg) {
      break;
    }

    buildSwapNeighborQuery(countQ, x, y, xg, yg);

    if (l[countQ].empty() || r[countQ].empty()) {
      l[countQ].clear();
      r[countQ].clear();
      return false;
    }
    char c = Query(countQ);

    if (c == '<') {
      swap(ans[x], ans[y]);
      comments[countQ] += "SwapNeighborSmall1 ";
      ResetMemory();
      isChange = true;
    }
    break;
  }
  return isChange;
}

void SwapNeighbor1Block(const vector<vector<int>>& blocks, int& countQ)
{
  int M = blocks.size();
  while (true) {
    int num = rand_xorshift() % (M - 1);
    vector<int> xb = blocks[num];
    vector<int> yb = blocks[num + 1];
    int xg = ans[xb[0]];
    int yg = ans[yb[0]];
    if (xg == yg) {
      continue;
    }
    for (int i = 0; i < N; ++i) {
      if (ans[i] == xg) {
        int ok = 1;
        for (auto x : xb) {
          if (i == x) ok = 0;
        }
        if (ok) {
          l[countQ].push_back(i);
        }
      }
      if (ans[i] == yg) {
        int ok = 1;
        for (auto y : yb) {
          if (i == y) ok = 0;
        }
        if (ok) {
          r[countQ].push_back(i);
        }
      }
    }
    if (l[countQ].empty() || r[countQ].empty()) {
      l[countQ].clear();
      r[countQ].clear();
      return;
    }
    char c = Query(countQ);

    if (c == '<') {
      for (auto x : xb) {
        ans[x] = yg;
      }
      for (auto y : yb) {
        ans[y] = xg;
      }
      comments[countQ] += "SwapNeighbor1Block ";
      ResetMemory();
    }
    break;
  }
}

void Method8(int hiritu = 60)
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    performStandardMoveSwap(countQ, hiritu);

    if (isTimeLimitApproaching()) {
      cerr << "Assert Method8" << endl;
      break;
    }
  }

  completeWithDummyQueries(countQ);
}

void Method208(int _kireme, int hiritu1, int hiritu2)
{
  initializeAnsArray();

  int kireme = round((double)Q * _kireme / 100);
  int hiritu = 100;

  int countQ = 0;
  while (countQ < Q) {
    if (countQ < kireme) {
      hiritu = hiritu1;
    }
    else {
      hiritu = hiritu2;
    }
    performStandardMoveSwap(countQ, hiritu);
  }

  completeWithDummyQueries(countQ);
}

void Method8_Two(int hiritu = 60)
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      Move1_Two(countQ);
    }
    else {
      if (Q - 1 <= countQ) {
        if (hiritu < 10) { break; }
        continue;
      }
      Swap1(countQ);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method19(int hiritu = 60)
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      Move1(countQ);
    }
    else {
      if (Q - 1 <= countQ) {
        if (hiritu < 10) { break; }
        continue;
      }
      SwapHalf(countQ);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method18(int hiritu = 60)
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      Move1Minimum(countQ);
    }
    else {
      if (Q - 1 <= countQ) {
        if (hiritu < 10) { break; }
        continue;
      }
      Swap1(countQ);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method13(int diffLine = 999)
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    int qu = rand_xorshift() % 100;
    if (qu < 10) {
      Move1(countQ);
    }
    else if (qu < 60) {
      Move1(countQ, 10);
    }
    else {
      if (Q - 1 <= countQ) {
        continue;
      }
      Swap1(countQ, diffLine);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method14()
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ < Q) {
    if (countQ < 0.75 * Q) {
      int qu = rand_xorshift() % 100;
      if (qu < 60) {
        Move1(countQ);
      }
      else {
        if (Q - 1 <= countQ) {
          continue;
        }
        Swap1(countQ, 999);
      }
    }
    else {
      if (Q - 1 <= countQ) {
        break;
      }
      Swap1(countQ, 1);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method15(int hiritu = 60)
{
  initializeAnsArray();

  int countQ = 0;
  int combo = -1;
  while (countQ < Q) {
    if (combo != -1) {
      combo = Move1Combo(countQ, combo);
    }
    else {
      int qu = rand_xorshift() % 100;
      if (qu < hiritu) {
        Move1(countQ);
      }
      else {
        if (Q - 1 <= countQ) {
          if (hiritu < 10) { break; }
          continue;
        }
        Swap1(countQ);
      }
    }
  }

  completeWithDummyQueries(countQ);
}

void Method16(int hiritu = 60)
{
  int countQ = 0;

  // 初期解を少し工夫
  for (int i = 0; i < N; ++i) {
    if (i < D) {
      ans[i] = i;
    }
    else {
      int x = rand_xorshift() % D;
      int y = rand_xorshift() % D;
      while (x == y) {
        y = rand_xorshift() % D;
      }
      for (int j = 0; j < N; ++j) {
        if (ans[j] == x) {
          l[countQ].push_back(j);
        }
        if (ans[j] == y) {
          r[countQ].push_back(j);
        }
      }

      char c = Query(countQ);
      if (c == '<') {
        ans[i] = x;
      }
      else {
        ans[i] = y;
      }
    }
  }

  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      Move1(countQ);
    }
    else {
      if (Q - 1 <= countQ) {
        if (hiritu < 10) { break; }
        continue;
      }
      Swap1(countQ);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method11()
{
  initializeAnsArray();

  int countQ = 0;
  while (countQ <= Q - 2) {
    int qu = rand_xorshift() % 100;
    if (qu < 33) {
      Move1(countQ);
    }
    else if (qu < 66) {
      if (Q - 2 <= countQ) {
        continue;
      }
      Swap1(countQ);
    }
    else {
      if (Q - 2 <= countQ) {
        continue;
      }
      Method11_1(countQ);
    }
  }

  completeWithDummyQueries(countQ);
}

void MergeDfs(vector<int>& items, int& countQ, int left, int right,
  int ikichi = 1001001)
{
  if (left + 1 == right) {
    return;
  }

  if (left + 2 == right) {
    int lhs = items[left];
    int rhs = items[left + 1];
    char c = Query1(countQ, lhs, rhs);
    if (c == '>') {
      swap(items[left], items[left + 1]);
    }
    return;
  }

  int mid = (left + right) / 2;
  MergeDfs(items, countQ, left, mid, ikichi);
  MergeDfs(items, countQ, mid, right, ikichi);

  int ite1 = left;
  int ite2 = mid;
  vector<int> tmp;
  if (right - left >= ikichi) {
    if (left == 0 && right == N) {
      while (ite1 < mid / 2 && ite2 < mid + (right - mid) / 2) {
        tmp.push_back(items[ite1]);
        ite1++;
        tmp.push_back(items[ite2]);
        ite2++;
      }
    }
  }

  while (ite1 < mid || ite2 < right) {
    if (ite1 == mid) {
      tmp.push_back(items[ite2]);
      ite2++;
    }
    else if (ite2 == right) {
      tmp.push_back(items[ite1]);
      ite1++;
    }
    else {
      int lhs = items[ite1];
      int rhs = items[ite2];
      char c = Query1(countQ, lhs, rhs);
      if (c == '>') {
        tmp.push_back(rhs);
        ite2++;
      }
      else {
        tmp.push_back(lhs);
        ite1++;
      }
    }
  }
  for (int i = 0; i < tmp.size(); ++i) { items[left + i] = tmp[i]; }
}

void MergeDfsBlock(vector<vector<int>>& blocks, int& countQ, int left,
  int right, int ikichi = 1001001)
{
  if (left + 1 == right) {
    return;
  }

  if (left + 2 == right) {
    l[countQ] = blocks[left];
    r[countQ] = blocks[left + 1];
    char c = Query(countQ);
    if (c == '>') {
      swap(blocks[left], blocks[left + 1]);
    }
    return;
  }

  int mid = (left + right) / 2;
  MergeDfsBlock(blocks, countQ, left, mid, ikichi);
  MergeDfsBlock(blocks, countQ, mid, right, ikichi);

  int ite1 = left;
  int ite2 = mid;
  vector<vector<int>> tmp;
  if (right - left >= ikichi) {
    if (left == 0 && right == N) {
      while (ite1 < mid / 2 && ite2 < mid + (right - mid) / 2) {
        tmp.push_back(blocks[ite1]);
        ite1++;
        tmp.push_back(blocks[ite2]);
        ite2++;
      }
    }
  }

  while (ite1 < mid || ite2 < right) {
    if (ite1 == mid) {
      tmp.push_back(blocks[ite2]);
      ite2++;
    }
    else if (ite2 == right) {
      tmp.push_back(blocks[ite1]);
      ite1++;
    }
    else {
      l[countQ] = blocks[ite1];
      r[countQ] = blocks[ite2];
      char c = Query(countQ);
      if (c == '>') {
        tmp.push_back(blocks[ite2]);
        ite2++;
      }
      else {
        tmp.push_back(blocks[ite1]);
        ite1++;
      }
    }
  }
  for (int i = 0; i < tmp.size(); ++i) { blocks[left + i] = tmp[i]; }
}

void MergeSort(vector<int>& items, int& countQ, int ikichi = 1001001, int _m = -1)
{
  if (_m == -1) _m = N;
  MergeDfs(items, countQ, 0, _m, ikichi);
}

void MergeSortBlock(vector<vector<int>>& blocks, int& countQ,
  int ikichi = 1001001)
{
  MergeDfsBlock(blocks, countQ, 0, blocks.size(), ikichi);
}

int CountMaxMergeSortDfs(int left, int right)
{
  if (left + 1 == right) {
    return 0;
  }

  if (left + 2 == right) {
    return 1;
  }

  int mid = (left + right) / 2;
  int cnt = 0;
  cnt += CountMaxMergeSortDfs(left, mid);
  cnt += CountMaxMergeSortDfs(mid, right);
  cnt += right - left - 1;
  return cnt;
}
int CountMaxMergeSort()
{
  return CountMaxMergeSortDfs(0, N);
}

void MergeDfs_Group(vector<int>& groups, int& countQ, int left, int right)
{
  if (left + 1 == right) {
    return;
  }

  if (left + 2 == right) {
    char c = QueryGroup(countQ, groups[left], groups[left + 1]);
    if (c == '>') {
      swap(groups[left], groups[left + 1]);
    }
    return;
  }

  int mid = (left + right) / 2;
  MergeDfs_Group(groups, countQ, left, mid);
  MergeDfs_Group(groups, countQ, mid, right);

  int ite1 = left;
  int ite2 = mid;
  vector<int> tmp;
  while (ite1 < mid || ite2 < right) {
    if (ite1 == mid) {
      tmp.push_back(groups[ite2]);
      ite2++;
    }
    else if (ite2 == right) {
      tmp.push_back(groups[ite1]);
      ite1++;
    }
    else {
      int lgId = groups[ite1];
      int rgId = groups[ite2];
      char c = QueryGroup(countQ, lgId, rgId);
      if (c == '>') {
        tmp.push_back(rgId);
        ite2++;
      }
      else {
        tmp.push_back(lgId);
        ite1++;
      }
    }
  }
  for (int i = 0; i < tmp.size(); ++i) { groups[left + i] = tmp[i]; }
}

void MergeSort_Group(vector<int>& groups, int& countQ)
{
  MergeDfs_Group(groups, countQ, 0, D);
  reverse(groups.begin(), groups.end());  // 重い順で返す
  return;
}

bool IsAllSearched_Swap2(const vector<int>& items, int minDiff)
{
  vector<int> ansItems[30];
  populateAnsItems(ansItems);
  vector<int> vl, vr;

  for (int num = 0; num < N - 1; ++num) {
    int x1 = items[num];
    int y1 = items[num + 1];
    int xg = ans[x1];
    int yg = ans[y1];
    if (xg == yg) {
      continue;
    }
    int x2 = -1;
    int y2 = -1;
    int lastX = 1001001;
    int lastI = 1001001;
    for (int ii = N - 1; ii >= 0; --ii) {
      int i = items[ii];
      if (ans[i] == xg && i != x1) {
        lastX = i;
        lastI = ii;
      }
      if (ans[i] == yg && i != y1) {
        if (lastI - ii < minDiff) {
          minDiff = lastI - ii;
          x2 = lastX;
          y2 = i;
        }
      }
    }
    if (x2 == -1) {
      continue;
    }
    vl.clear();
    vr.clear();
    for (int i = 0; i < N; ++i) {
      if (ans[i] == xg && i != x1 && i != x2) {
        vl.push_back(i);
      }
      if (ans[i] == yg && i != y1 && i != y2) {
        vr.push_back(i);
      }
    }
    if (vl.empty() || vr.empty()) {
      continue;
    }
    char c1 = QueryMapDirectly(vl, vr);

    {
      vl.clear();
      vr.clear();
      vl.push_back(x1);
      vl.push_back(x2);
      vr.push_back(y1);
      vr.push_back(y2);
    }
    char c2 = QueryMapDirectly(vl, vr);

    if (c1 == '?' || c2 == '?') {
      return false;
    }

    if (c1 != '=' && c1 == c2) {
      return false;
    }
  }

  return true;
}

bool Swap2(const vector<int>& items, int& countQ, int minDiff = 10, int _m = -1)
{
  int M = _m;
  if (_m == -1) M = N;
  M = max(M, 5);
  bool isChange = false;
  int loop = 0;
  while (true) {
    loop++;
    int num = rand_xorshift() % (M - 1);
    int x1 = items[num];
    int y1 = items[num + 1];
    int xg = ans[x1];
    int yg = ans[y1];
    if (xg == yg) {
      if (loop > 100) {
        while (xg == yg) {
          num = rand_xorshift() % (N - 1);
          x1 = items[num];
          y1 = items[num + 1];
          xg = ans[x1];
          yg = ans[y1];
        }
      }
      else {
        continue;
      }
    }
    int x2 = -1;
    int y2 = -1;
    int lastX = 1001001;
    int lastI = 1001001;
    for (int ii = N - 1; ii >= 0; --ii) {
      int i = items[ii];
      if (ans[i] == xg && i != x1) {
        lastX = i;
        lastI = ii;
      }
      if (ans[i] == yg && i != y1) {
        if (lastI - ii < minDiff) {
          minDiff = lastI - ii;
          x2 = lastX;
          y2 = i;
        }
      }
    }
    if (x2 == -1) {
      break;
    }
    for (int i = 0; i < N; ++i) {
      if (ans[i] == xg && i != x1 && i != x2) {
        l[countQ].push_back(i);
      }
      if (ans[i] == yg && i != y1 && i != y2) {
        r[countQ].push_back(i);
      }
    }
    if (l[countQ].empty() || r[countQ].empty()) {
      l[countQ].clear();
      r[countQ].clear();
      return false;
    }
    char c1 = Query(countQ);

    {
      l[countQ].push_back(x1);
      l[countQ].push_back(x2);
      r[countQ].push_back(y1);
      r[countQ].push_back(y2);
    }
    char c2 = Query(countQ);

    if (c1 == c2) {
      swap(ans[x1], ans[y1]);
      swap(ans[x2], ans[y2]);
      comments[countQ] += "Swap2 ";
      ResetMemory();
      if (c1 != '=') {
        isChange = true;
      }
    }
    break;
  }
  return isChange;
}

bool SwapN(const vector<int>& items, int& countQ, int minDiff)
{
  bool isChange = false;
  int xg = rand_xorshift() % D;
  int yg = rand_xorshift() % D;
  while (xg == yg) {
    yg = rand_xorshift() % D;
  }

  vector<P> vpx, vpy;
  {
    int lastX = 1001001;
    int lastI = 1001001;
    for (int ii = N - 1; ii >= 0; --ii) {
      int i = items[ii];
      if (ans[i] == xg) {
        lastX = i;
        lastI = ii;
      }
      if (ans[i] == yg) {
        if (lastI - ii <= minDiff) {
          vpx.push_back(P(lastX, i));
          lastX = 1001001;
          lastI = 1001001;
        }
      }
    }
  }

  {
    int lastY = 1001001;
    int lastI = 1001001;
    for (int ii = N - 1; ii >= 0; --ii) {
      int i = items[ii];
      if (ans[i] == yg) {
        lastY = i;
        lastI = ii;
      }
      if (ans[i] == xg) {
        if (lastI - ii <= minDiff) {
          vpy.push_back(P(i, lastY));
          lastY = 1001001;
          lastI = 1001001;
        }
      }
    }
  }

  if (vpx.empty() || vpy.empty()) {
    return false;
  }

  vector<P> vpx2, vpy2;
  set<int> use;
  int loop = 0;
  while (true) {
    int ok = 0;

    for (int winter = 0; winter < 10; ++winter) {
      if (loop % 2 == 0) {
        P p = vpx[rand_xorshift() % vpx.size()];
        int x = p.first;
        int y = p.second;
        if (use.find(x) == use.end() && use.find(y) == use.end()) {
          vpx2.push_back(p);
          use.insert(x);
          use.insert(y);
          ok = 1;
          break;
        }
      }
      else {
        P p = vpy[rand_xorshift() % vpy.size()];
        int x = p.first;
        int y = p.second;
        if (use.find(x) == use.end() && use.find(y) == use.end()) {
          vpy2.push_back(p);
          use.insert(x);
          use.insert(y);
          ok = 1;
          break;
        }
      }
    }

    if (ok == 0) {
      break;
    }
    loop++;
  }

  if (vpx2.empty() || vpy2.empty()) {
    return false;
  }

  std::shuffle(vpx2.begin(), vpx2.end(), engine_mt19937);
  std::shuffle(vpy2.begin(), vpy2.end(), engine_mt19937);

  int lNum = rand_xorshift() % 3 + 1;
  lNum = min(lNum, (int)vpx2.size());
  int rNum = rand_xorshift() % 3 + 1;
  rNum = min(rNum, (int)vpy2.size());

  vector<int> vQuery1X, vQuery1Y, vQuery2X, vQuery2Y;
  for (int i = 0; i < lNum; ++i) {
    vQuery1X.push_back(vpx2[i].first);
    vQuery1Y.push_back(vpx2[i].second);
  }
  for (int i = 0; i < rNum; ++i) {
    vQuery1X.push_back(vpy2[i].first);
    vQuery1Y.push_back(vpy2[i].second);
  }
  for (int i = 0; i < N; ++i) {
    if (ans[i] == xg) {
      int notUse = 1;
      for (auto xx : vQuery1X) {
        if (xx == i) {
          notUse = 0;
        }
      }
      if (notUse) {
        vQuery2X.push_back(i);
      }
    }
    if (ans[i] == yg) {
      int notUse = 1;
      for (auto yy : vQuery1Y) {
        if (yy == i) {
          notUse = 0;
        }
      }
      if (notUse) {
        vQuery2Y.push_back(i);
      }
    }
  }

  l[countQ] = vQuery1X;
  r[countQ] = vQuery1Y;
  char c1 = Query(countQ);
  l[countQ] = vQuery2X;
  r[countQ] = vQuery2Y;
  char c2 = Query(countQ);

  if (c1 == c2) {
    for (auto xx : vQuery1X) {
      ans[xx] = yg;
    }
    for (auto yy : vQuery1Y) {
      ans[yy] = xg;
    }
    comments[countQ] += "SwapN ";
    ResetMemory();
    if (c1 != '=') {
      isChange = true;
    }
  }

  return isChange;
}

void Swap2Block(const vector<vector<int>>& blocks, int& countQ)
{
  int M = blocks.size();
  while (true) {
    int num = rand_xorshift() % (M - 1);
    vector<int> xb1 = blocks[num];
    vector<int> yb1 = blocks[num + 1];
    int xg = ans[xb1[0]];
    int yg = ans[yb1[0]];

    if (xg == yg) {
      continue;
    }

    vector<int> vx, vy;
    for (int i = 0; i < M; ++i) {
      if (ans[blocks[i][0]] == xg && i != num) {
        vx.push_back(i);
      }
      if (ans[blocks[i][0]] == yg && i != num + 1) {
        vy.push_back(i);
      }
    }
    if (vx.empty() || vy.empty()) {
      break;
    }
    int x2 = vx[rand_xorshift() % vx.size()];
    int y2 = vy[rand_xorshift() % vy.size()];

    for (int i = 0; i < N; ++i) {
      if (ans[i] == xg && i != num && i != x2) {
        l[countQ].push_back(i);
      }
      if (ans[i] == yg && i != num + 1 && i != y2) {
        r[countQ].push_back(i);
      }
    }

    char c1 = Query(countQ);

    vector<int> vl, vr;
    {
      l[countQ] = xb1;
      for (auto x : blocks[x2]) l[countQ].push_back(x);
      vl = l[countQ];
      r[countQ] = yb1;
      for (auto y : blocks[y2]) r[countQ].push_back(y);
      vr = r[countQ];
    }
    char c2 = Query(countQ);

    if (c1 == c2) {
      for (auto x : vl) ans[x] = yg;
      for (auto y : vr) ans[y] = xg;
      comments[countQ] += "Swap2Block ";
      ResetMemory();
    }
    break;
  }
}

bool IsAllCandidateSearched_Method6(const vector<int>& items, int minDiff)
{
  if (!IsAllSearched_SwapNeighbor1(items)) {
    return false;
  }

  if (!IsAllSearched_Swap2(items, minDiff)) {
    return false;
  }

  return true;
}

// ========== Common Helper Functions for Method226/266 ==========

// Calculate max_D (heaviest group)
int calculateMaxD(int& countQ)
{
  int tmp_max_D = 0;
  for (int i = 1; i < D; ++i) {
    for (int j = 0; j < N; ++j) {
      if (ans[j] == tmp_max_D) {
        l[countQ].push_back(j);
      }
      if (ans[j] == i) {
        r[countQ].push_back(j);
      }
    }
    char c = Query(countQ);
    if (c == '<') {
      tmp_max_D = i;
    }
  }
  return tmp_max_D;
}

// Compare with provisional best and update if better
void compareAndUpdateBest(int tmp_max_D, int& real_max_D, int real_ans[], int& countQ)
{
  int cntBoth[110] = {};
  for (int j = 0; j < N; ++j) {
    if (ans[j] == tmp_max_D) {
      cntBoth[j]++;
    }
  }
  for (int j = 0; j < N; ++j) {
    if (real_ans[j] == real_max_D) {
      cntBoth[j]++;
    }
  }
  for (int j = 0; j < N; ++j) {
    if (ans[j] == tmp_max_D && cntBoth[j] == 1) {
      l[countQ].push_back(j);
    }
  }
  for (int j = 0; j < N; ++j) {
    if (real_ans[j] == real_max_D && cntBoth[j] == 1) {
      r[countQ].push_back(j);
    }
  }
  char c2 = Query(countQ);
  if (c2 == '<') {
    real_max_D = tmp_max_D;
    for (int j = 0; j < N; ++j) {
      real_ans[j] = ans[j];
    }
  }
}

// Common optimization loop for Method226/266
bool runOptimizationWithRestarts(vector<int>& items, int& countQ, int hiritu, int minDiff,
  int maxFailedCount, bool checkTime = true)
{
  int failedCount = 0;

  while (countQ < Q - D) {
    if (checkTime && isTimeLimitApproaching(0.1)) {
      cerr << "TLE : Method226 " << "time = " << elapsed << endl;
      break;
    }

    int qu = getRandomPercentage();
    bool isChange = false;

    if (qu < hiritu) {
      isChange = SwapNeighbor1(items, countQ);
      if (isChange) {
        failedCount = 0;
      }
    }
    else {
      if (isNearQueryLimit(countQ, D + 1)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      isChange = Swap2(items, countQ, minDiff);
      if (isChange) {
        failedCount = 0;
      }
    }

    failedCount++;
    if (failedCount > maxFailedCount) {
      break;
    }
  }

  return true;
}

// Save current state to real_ans
void saveCurrentState(int real_ans[])
{
  for (int i = 0; i < N; ++i) {
    real_ans[i] = ans[i];
  }
}

// Restore state from keepInitialAns
void restoreInitialState(int keepInitialAns[])
{
  for (int i = 0; i < N; ++i) {
    ans[i] = keepInitialAns[i];
  }
}

// Random swap for Method266
void randomSwap(vector<int>& items, int kosuu, int saidai)
{
  vector<int> ningning;
  saidai = min(saidai, N - 1);
  saidai = max(saidai, 1);
  for (int i = 0; i < saidai; ++i) {
    ningning.push_back(i);
  }
  std::shuffle(ningning.begin(), ningning.end(), engine_mt19937);
  kosuu = min(kosuu, saidai);
  kosuu = max(kosuu, 1);
  for (int i = 0; i < kosuu; ++i) {
    swap(ans[items[i]], ans[items[i + 1]]);
  }
}

// ========== End of Common Helper Functions for Method226/266 ==========

// ========== Common Helper Functions for Method106/246/112 ==========

// Initialize items with merge sort and assign heaviest items to groups
void initializeItemsWithHeaviestAssignment(vector<int>& items, int& countQ, int ikichi = 1001001)
{
  initializeItems(items);
  MergeSort(items, countQ, ikichi);
  updateKarusaArray(items);

  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    int id = items[N - 1 - i];
    ans[id] = i;
  }
}

// Assign items using PseudoItems until totyuu threshold
void assignItemsWithPseudoQuery(vector<int>& items, vector<int>& groups, int& countQ, int totyuu)
{
  for (int i = N - D - 1; i >= 0; --i) {
    if (N - i >= totyuu) {
      break;
    }
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;

    int left = 0;
    int right = D - 1;
    while (left < right) {
      int mid = (left + right) / 2;
      for (int j = 0; j < N; ++j) {
        if (ans[j] == gId) {
          l[countQ].push_back(j);
        }
        if (ans[j] == groups[mid]) {
          r[countQ].push_back(j);
        }
      }
      char c = PseudoItemsQuery(countQ);
      l[countQ].clear();
      r[countQ].clear();
      if (c == '=') {
        left = right = mid;
      }
      else if (c == '>') {
        right = mid;
      }
      else {
        left = mid + 1;
      }
    }
    moveGroupToPosition(groups, D - 1, left);
  }
}

// Find tail position for remaining items (Method106/246 style)
int findTailPosition(const vector<int>& items)
{
  int tail = N;
  while (tail > 0) {
    if (ans[items[tail - 1]] != -1) {
      tail--;
    }
    else {
      break;
    }
  }
  return tail;
}

// Assign remaining items to lightest group
void assignRemainingItems(const vector<int>& items, vector<int>& groups, int& countQ, int startPos, int endPos = 0)
{
  for (int i = startPos - 1; i >= endPos; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;
    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }
}

// Common structure for Method106/246/112
void initializeWithPseudoQueryMethod(vector<int>& items, vector<int>& groups, int& countQ, int totyuu, int ikichi = 1001001)
{
  // Initialize items and assign heaviest
  initializeItemsWithHeaviestAssignment(items, countQ, ikichi);

  // Initialize groups
  initializeGroups(groups);

  // Assign items with pseudo query until totyuu
  assignItemsWithPseudoQuery(items, groups, countQ, totyuu);

  // Merge sort groups
  MergeSort_Group(groups, countQ);
}

// ========== End of Common Helper Functions for Method106/246/112 ==========

void Method266(int hiritu, int minDiff, int kosuu, int saidai, int maxFailedCount = 1000)
{
  if (maxFailedCount < 10) {
    maxFailedCount = 10;
  }
  vector<int> items;
  initializeItems(items);
  int countQ = 0;

  // アイテムをマージソート(軽い順)
  MergeSort(items, countQ);
  updateKarusaArray(items);

  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    int id = items[N - 1 - i];
    ans[id] = i;
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);

  // 一番軽いグループに入れていく
  for (int i = N - D - 1; i >= 0; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;
    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }

  int keepInitialAns[110];
  for (int i = 0; i < N; ++i) {
    keepInitialAns[i] = ans[i];
  }

  int setCount = 0;
  while (countQ <= Q - D) {
    // Run optimization
    runOptimizationWithRestarts(items, countQ, hiritu, minDiff, maxFailedCount, false);

    if (setCount == 0) {
      saveCurrentState(real_ans);
      if (countQ <= Q - D) {
        real_max_D = calculateMaxD(countQ);
      }
      else {
        break;
      }
    }
    else {
      if (countQ <= Q - D) {
        int tmp_max_D = calculateMaxD(countQ);
        compareAndUpdateBest(tmp_max_D, real_max_D, real_ans, countQ);
      }
      else {
        break;
      }
    }

    setCount++;

    if (isTimeLimitApproaching()) {
      if (mode != 0) {
        cout << "Assert Method226 : N = " << N << ", Q = " << Q << ", setCount = " << setCount << endl;
      }
      break;
    }

    // ランダムにいくつかスワップ
    randomSwap(items, kosuu, saidai);
  }

  if (setCount > 0) {
    for (int i = 0; i < N; ++i) {
      ans[i] = real_ans[i];
    }
  }

  completeWithDummyQueries(countQ);
}

// Common initialization for methods using merge sort
void initializeWithMergeSort(vector<int>& items, vector<int>& groups, int& countQ, int ikichi = 1001001)
{
  initializeItems(items);
  MergeSort(items, countQ, ikichi);
  updateKarusaArray(items);

  assignHeaviestToGroups(items, N - 1);
  initializeGroups(groups);

  for (int i = N - D - 1; i >= 0; --i) {
    assignToLightestGroup(items, i, countQ, groups);
  }
}

void Method226(int hiritu = 100, int minDiff = 10)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeWithMergeSort(items, groups, countQ);

  int keepInitialAns[110];
  for (int i = 0; i < N; ++i) {
    keepInitialAns[i] = ans[i];
  }

  int setCount = 0;
  while (countQ <= Q - D) {
    // Run optimization with time check
    runOptimizationWithRestarts(items, countQ, hiritu, minDiff, 1000, true);

    if (setCount == 0) {
      saveCurrentState(real_ans);
      if (countQ <= Q - D) {
        real_max_D = calculateMaxD(countQ);
      }
      else {
        break;
      }
    }
    else {
      if (countQ <= Q - D) {
        int tmp_max_D = calculateMaxD(countQ);
        compareAndUpdateBest(tmp_max_D, real_max_D, real_ans, countQ);
      }
      else {
        break;
      }
    }

    setCount++;

    if (isTimeLimitApproaching()) {
      if (mode != 0) {
        cout << "Assert Method226 : N = " << N << ", Q = " << Q << ", setCount = " << setCount << ", time = " << elapsed << endl;
      }
      break;
    }

    restoreInitialState(keepInitialAns);
  }

  if (setCount > 0) {
    for (int i = 0; i < N; ++i) {
      ans[i] = real_ans[i];
    }
  }

  completeWithDummyQueries(countQ);
}

void Method6(int hiritu = 100, int minDiff = 10)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeWithMergeSort(items, groups, countQ);

  int failedCount = 0;
  while (countQ < Q) {
    hiritu = adjustHirituByTime(hiritu);
    if (hiritu == 90) {
      minDiff = 999;
    }
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      bool isChange = SwapNeighbor1(items, countQ);
      if (isChange) {
        failedCount = 0;
      }
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }

      bool isChange = Swap2(items, countQ, minDiff);
      if (isChange) {
        failedCount = 0;
      }
    }
    failedCount++;
  }
  completeWithDummyQueries(countQ);
}

void Method706(int hiritu1, int minDiff, int hiritu2)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeWithMergeSort(items, groups, countQ);

  int failedCount = 0;
  if (hiritu1 < 10) {
    hiritu1 = 10;
  }
  while (countQ < Q) {
    if (elapsed > TL / 3) {
      hiritu1 = 33;
      hiritu2 = 66;
      minDiff = 999;
    }
    int qu = rand_xorshift() % 100;
    if (qu < hiritu1) {
      bool isChange = SwapNeighbor1(items, countQ);
      if (isChange) {
        failedCount = 0;
      }
    }
    else if (qu < hiritu2) {
      bool isChange = Swap2(items, countQ, minDiff);
      if (isChange) {
        failedCount = 0;
      }
    }
    else {
      if (failedCount > 500) {
        int minDiff2 = minDiff;
        minDiff2 = max(1, minDiff2);
        bool isChange = SwapN(items, countQ, minDiff2);
        if (isChange) {
          failedCount = 0;
        }
      }
    }
    failedCount++;
  }
  completeWithDummyQueries(countQ);
}

// Common initialization for block-based methods
bool initializeBlocks(vector<vector<int>>& blocks, vector<int>& groups, int& countQ, int blockSize)
{
  for (int i = 0; i < N; ++i) {
    if (i % blockSize == 0) {
      blocks.push_back({});
    }
    blocks.back().push_back(i);
  }
  if (blocks.size() < D) {
    for (int i = 0; i < N; ++i) { ans[i] = i % D; }
    completeWithDummyQueries(countQ);
    return false;
  }

  int M = blocks.size();
  MergeSortBlock(blocks, countQ);

  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < blocks[M - 1 - i].size(); ++j) {
      int id = blocks[M - 1 - i][j];
      ans[id] = i;
    }
  }

  initializeGroups(groups);

  // 一番軽いグループに入れていく
  for (int i = M - D - 1; i >= 0; --i) {
    int gId = groups[D - 1];
    for (int j = 0; j < blocks[i].size(); ++j) {
      int id = blocks[i][j];
      ans[id] = gId;
    }

    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }
  return true;
}

void Method206(int hiritu1, int hiritu2, int timing, int blockSize)
{
  int countQ = 0;
  vector<vector<int>> blocks;
  vector<int> groups;

  if (!initializeBlocks(blocks, groups, countQ, blockSize)) { return; }

  int loopCount = 0;
  while (countQ < Q) {
    loopCount++;
    int isKouhan = false;

    if (countQ >= round((double)Q * timing / 100)) {
      isKouhan = true;
    }
    if (loopCount > 1000) {
      isKouhan = true;
    }

    if (isKouhan) {
      int qu = rand_xorshift() % 100;
      int hiritu = hiritu1;
      if (qu < hiritu) {
        Move1(countQ);
      }
      else {
        if (Q - 1 <= countQ) {
          if (hiritu < 10) { break; }
          continue;
        }
        Swap1(countQ);
      }
    }
    else {
      int hiritu = hiritu2;
      if (hiritu >= 100 && elapsed > TL / 3) {
        hiritu = 90;
      }
      int qu = rand_xorshift() % 100;
      if (qu < hiritu) {
        SwapNeighbor1Block(blocks, countQ);
      }
      else {
        if (countQ >= Q - 1) {
          if (hiritu < 10) {
            break;
          }
          continue;
        }

        Swap2Block(blocks, countQ);
      }
    }
  }
  completeWithDummyQueries(countQ);
}

// ========== Common Helper Functions for Method216/316/916/516 ==========

// Initialize blocks from items
bool initializeBlocks(vector<vector<int>>& blocks, int blockSize, int& countQ)
{
  for (int i = 0; i < N; ++i) {
    if (i % blockSize == 0) {
      blocks.push_back({});
    }
    blocks.back().push_back(i);
  }
  if (blocks.size() < D) {
    for (int i = 0; i < N; ++i) { ans[i] = i % D; }
    completeWithDummyQueries(countQ);
    return false;  // Early exit
  }
  return true;
}

// Sort blocks and handle edge cases
bool sortAndValidateBlocks(vector<vector<int>>& blocks, int destroySize, int& countQ)
{
  int M = blocks.size();
  MergeSortBlock(blocks, countQ);

  if (destroySize >= M) {
    for (int i = 0; i < N; ++i) { ans[i] = i % D; }
    completeWithDummyQueries(countQ);
    return false;  // Early exit
  }
  return true;
}

// Split large blocks into individual items
void splitLargeBlocks(vector<vector<int>>& blocks, int destroySize,
  vector<vector<int>>& tmp1, vector<vector<int>>& tmp2)
{
  int M = blocks.size();
  for (int i = 0; i < M; ++i) {
    if (i < M - destroySize) {
      tmp2.push_back(blocks[i]);
    }
    else {
      for (int j = 0; j < blocks[i].size(); ++j) {
        tmp1.push_back({ blocks[i][j] });
      }
    }
  }
}

// Insert split items back into blocks using binary search
void insertSplitItems(vector<vector<int>>& blocks, const vector<vector<int>>& tmp1,
  const vector<vector<int>>& tmp2, int& countQ)
{
  blocks = tmp2;
  for (auto tmp : tmp1) {
    int M = blocks.size();
    int left = 0;
    int right = M;
    while (left < right) {
      int mid = (left + right) / 2;
      l[countQ] = tmp;
      r[countQ] = blocks[mid];
      char c = Query(countQ);
      if (c == '=') {
        left = right = mid;
      }
      else if (c == '<') {
        right = mid;
      }
      else {
        left = mid + 1;
      }
    }
    blocks.insert(blocks.begin() + left, tmp);
  }
}

// Create items array from single-item blocks
int createItemsArray(const vector<vector<int>>& blocks, vector<int>& items)
{
  int used[110] = {};
  for (auto block : blocks) {
    if (block.size() == 1) {
      items.push_back(block[0]);
      used[block[0]] = 1;
    }
  }
  int M2 = items.size();
  for (int i = 0; i < N; ++i) {
    if (used[i] == 0) {
      items.push_back(i);
    }
  }
  return M2;
}

// Check if should switch to later phase optimization
bool shouldSwitchToLaterPhase(int countQ, int loopCount, int timing)
{
  if (countQ >= round((double)Q * timing / 100)) {
    return true;
  }
  if (loopCount > 1000) {
    return true;
  }
  return false;
}

// Main optimization loop handler
void runOptimizationLoop(vector<vector<int>>& blocks, vector<int>& items, int M2,
  int& countQ, int hiritu1, int hiritu2, int timing,
  bool useMethod216Style)
{
  int loopCount = 0;
  while (countQ < Q) {
    loopCount++;
    bool isKouhan = shouldSwitchToLaterPhase(countQ, loopCount, timing);

    if (isKouhan) {
      if (useMethod216Style) {
        int qu = rand_xorshift() % 100;
        int hiritu = hiritu1;
        if (qu < hiritu) {
          Move1(countQ);
        }
        else {
          if (Q - 1 <= countQ) {
            if (hiritu < 10) { break; }
            continue;
          }
          Swap1(countQ);
        }
      }
      else {
        SwapNeighbor1(items, countQ, M2);
      }
    }
    else {
      int hiritu = hiritu2;
      if (hiritu >= 100 && elapsed > TL / 3) {
        hiritu = 90;
      }
      int qu = rand_xorshift() % 100;
      if (qu < hiritu) {
        SwapNeighbor1Block(blocks, countQ);
      }
      else {
        if (countQ >= Q - 1) {
          if (hiritu < 10) {
            break;
          }
          continue;
        }
        Swap2Block(blocks, countQ);
      }
    }
  }
  completeWithDummyQueries(countQ);
}

// ========== End of Common Helper Functions ==========

void Method216(int hiritu1, int hiritu2, int timing, int blockSize, int destroySize)
{
  int countQ = 0;
  vector<vector<int>> blocks;

  // Initialize blocks
  if (!initializeBlocks(blocks, blockSize, countQ)) {
    return;
  }

  // Sort and validate
  if (!sortAndValidateBlocks(blocks, destroySize, countQ)) {
    return;
  }

  // Split large blocks
  vector<vector<int>> tmp1, tmp2;
  splitLargeBlocks(blocks, destroySize, tmp1, tmp2);

  // Insert split items
  insertSplitItems(blocks, tmp1, tmp2, countQ);

  int M = blocks.size();
  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < blocks[M - 1 - i].size(); ++j) {
      int id = blocks[M - 1 - i][j];
      ans[id] = i;
    }
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);

  // 一番軽いグループに入れていく
  for (int i = M - D - 1; i >= 0; --i) {
    int gId = groups[D - 1];
    for (int j = 0; j < blocks[i].size(); ++j) {
      int id = blocks[i][j];
      ans[id] = gId;
    }

    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }

  // Run optimization loop (Method216 style)
  vector<int> dummyItems;  // Not used in Method216
  runOptimizationLoop(blocks, dummyItems, 0, countQ, hiritu1, hiritu2, timing, true);
}

void Method316(int hiritu1, int hiritu2, int timing, int blockSize, int destroySize)
{
  int countQ = 0;
  vector<vector<int>> blocks;

  // Initialize blocks
  if (!initializeBlocks(blocks, blockSize, countQ)) {
    return;
  }

  // Sort and validate
  if (!sortAndValidateBlocks(blocks, destroySize, countQ)) {
    return;
  }

  // Split large blocks
  vector<vector<int>> tmp1, tmp2;
  splitLargeBlocks(blocks, destroySize, tmp1, tmp2);

  // Insert split items
  insertSplitItems(blocks, tmp1, tmp2, countQ);

  // Create items array
  vector<int> items;
  int M2 = createItemsArray(blocks, items);

  if (M2 < 5) {
    for (int i = 0; i < N; ++i) {
      ans[i] = i % D;
    }
    completeWithDummyQueries(countQ);
    return;
  }

  int M = blocks.size();
  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < blocks[M - 1 - i].size(); ++j) {
      int id = blocks[M - 1 - i][j];
      ans[id] = i;
    }
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);

  // 一番軽いグループに入れていく
  for (int i = M - D - 1; i >= 0; --i) {
    int gId = groups[D - 1];
    for (int j = 0; j < blocks[i].size(); ++j) {
      int id = blocks[i][j];
      ans[id] = gId;
    }

    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }

  // Run optimization loop (Method316 style - uses items)
  runOptimizationLoop(blocks, items, M2, countQ, hiritu1, hiritu2, timing, false);
}

void Method916(int totyuu, int hiritu1, int hiritu2, int timing, int blockSize, int destroySize)
{
  int countQ = 0;
  vector<vector<int>> blocks;

  // Initialize blocks
  if (!initializeBlocks(blocks, blockSize, countQ)) {
    return;
  }

  // Sort and validate
  if (!sortAndValidateBlocks(blocks, destroySize, countQ)) {
    return;
  }

  // Split large blocks
  vector<vector<int>> tmp1, tmp2;
  splitLargeBlocks(blocks, destroySize, tmp1, tmp2);

  // Insert split items
  insertSplitItems(blocks, tmp1, tmp2, countQ);

  // Create items array
  vector<int> items;
  int M2 = createItemsArray(blocks, items);

  if (M2 < 5) {
    for (int i = 0; i < N; ++i) {
      ans[i] = i % D;
    }
    completeWithDummyQueries(countQ);
    return;
  }

  int M = blocks.size();
  vector<int> groupNakami[30];
  int flag[110] = {};
  // 各グループに1個ずつ入れる
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < blocks[M - 1 - i].size(); ++j) {
      int id = blocks[M - 1 - i][j];
      ans[id] = i;
    }
    groupNakami[i].push_back(M - 1 - i);
    flag[M - 1 - i] = 1;
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする

  for (int i = 0; i < D; ++i) {
    groups.push_back(i);
  }

  // 途中までpseudoItemsを用いて一番軽いグループに入れていく
  for (int i = M - D - 1; i >= 0; --i) {
    if (M - i >= totyuu) {
      break;
    }
    int gId = groups[D - 1];
    for (int j = 0; j < blocks[i].size(); ++j) {
      int id = blocks[i][j];
      ans[id] = gId;
    }
    groupNakami[gId].push_back(i);
    flag[i] = 1;

    int left = 0;
    int right = D - 1;
    while (left < right) {
      int mid = (left + right) / 2;
      int midId = groups[mid];
      char c = PseudoItemsQueryGroup(groupNakami[gId], groupNakami[midId]);
      if (c == '=') {
        left = right = mid;
      }
      else if (c == '>') {
        right = mid;
      }
      else {
        left = mid + 1;
      }
    }
    moveGroupToPosition(groups, D - 1, left);
  }

  // ここでグループを一度マージソートして正しい順序に並び替える
  MergeSort_Group(groups, countQ);

  // 一番軽いグループに入れていく
  int tail = M;
  while (tail > 0) {
    if (flag[tail - 1] == 1) {
      tail--;
    }
    else {
      break;
    }
  }
  for (int i = tail - 1; i >= 0; --i) {
    int gId = groups[D - 1];
    for (int j = 0; j < blocks[i].size(); ++j) {
      int id = blocks[i][j];
      ans[id] = gId;
    }

    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }

  // Run optimization loop (Method916 style - uses items)
  runOptimizationLoop(blocks, items, M2, countQ, hiritu1, hiritu2, timing, false);
}

void Method516(int hiritu1, int hiritu2, int timing, int blockSize, int destroySize)
{
  int countQ = 0;
  vector<vector<int>> blocks;

  // Initialize blocks
  if (!initializeBlocks(blocks, blockSize, countQ)) {
    return;
  }

  // Sort and validate
  if (!sortAndValidateBlocks(blocks, destroySize, countQ)) {
    return;
  }

  // Split large blocks
  vector<vector<int>> tmp1, tmp2;
  splitLargeBlocks(blocks, destroySize, tmp1, tmp2);

  // Insert split items
  insertSplitItems(blocks, tmp1, tmp2, countQ);

  // Create items array
  vector<int> items;
  int M2 = createItemsArray(blocks, items);

  if (M2 < 5) {
    for (int i = 0; i < N; ++i) {
      ans[i] = i % D;
    }
    completeWithDummyQueries(countQ);
    return;
  }

  int M = blocks.size();
  // つづら折りに入れていく
  int junban[60];
  for (int i = 0; i < D; ++i) {
    junban[i] = i;
    junban[D * 2 - 1 - i] = i;
  }
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < blocks[M - 1 - i].size(); ++j) {
      int id = blocks[M - 1 - i][j];
      ans[id] = junban[i % (D * 2)];
    }
  }

  // Run optimization loop (Method516 style - uses items)
  runOptimizationLoop(blocks, items, M2, countQ, hiritu1, hiritu2, timing, false);
}

void Method12(int ikichi = N)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeWithMergeSort(items, groups, countQ, ikichi);

  while (countQ < Q) {
    SwapNeighbor1(items, countQ);
  }

  completeWithDummyQueries(countQ);
}

void Method106(int hiritu = 100, int minDiff = 10, int totyuu = 999)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  // Common initialization with pseudo query method
  initializeWithPseudoQueryMethod(items, groups, countQ, totyuu);

  // Find tail position and assign remaining items
  int tail = findTailPosition(items);
  assignRemainingItems(items, groups, countQ, tail);

  // Optimization loop
  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      SwapNeighbor1(items, countQ);
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      Swap2(items, countQ, minDiff);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method306(int hiritu = 100, int minDiff = 10, int totyuu = 999, int _m = 30, int aketoku = 0)
{
  int M = _m;
  M = min(M, N);
  M = max(M, 10);
  vector<int> items;
  initializeItems(items);
  int countQ = 0;

  // 一部アイテムをマージソート(軽い順)
  MergeSort(items, countQ, 1001001, M);

  // ソートしていないアイテムを入れる
  int _D = D - aketoku;
  _D = max(_D, 2);
  for (int i = M; i < N; ++i) {
    int id = items[i];
    ans[id] = (i - M) % _D;
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);
  // ここでグループを一度マージソートして正しい順序に並び替える
  MergeSort_Group(groups, countQ);
  // 一番軽いグループに入れていく
  int tail = N;
  while (tail > 0) {
    if (ans[items[tail - 1]] != -1) {
      tail--;
    }
    else {
      break;
    }
  }
  for (int i = tail - 1; i >= 0; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;
    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }
  int loop = 0;
  while (countQ < Q) {
    loop++;

    if (loop > 100000) {
      break;
    }
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      SwapNeighbor1(items, countQ, M);
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      Swap2(items, countQ, minDiff, M);
    }
  }
  completeWithDummyQueries(countQ);
}

void Method806(int hiritu, int minDiff, int _m)
{
  int M = _m;
  M = min(M, N);
  M = max(M, 10);
  vector<int> items;
  initializeItems(items);
  int countQ = 0;

  // 1回だけ比べる
  for (int i = 0; i < M; ++i) {
    int j = N - 1 - i;
    if (j <= i) { break; }
    char c = Query1(countQ, i, j);
    if (c == '>') {
      swap(items[i], items[j]);
    }
  }

  // 一部アイテムをマージソート(軽い順)
  MergeSort(items, countQ, 1001001, M);

  // ソートしていないアイテムを入れる
  for (int i = M; i < N; ++i) {
    int id = items[i];
    ans[id] = (i - M) % D;
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);
  // ここでグループを一度マージソートして正しい順序に並び替える
  MergeSort_Group(groups, countQ);
  // 一番軽いグループに入れていく
  int tail = N;
  while (tail > 0) {
    if (ans[items[tail - 1]] != -1) {
      tail--;
    }
    else {
      break;
    }
  }
  for (int i = tail - 1; i >= 0; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;
    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }
  int loop = 0;
  while (countQ < Q) {
    loop++;

    if (loop > 100000) {
      break;
    }
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      SwapNeighbor1(items, countQ, M);
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      Swap2(items, countQ, minDiff, M);
    }
  }
  completeWithDummyQueries(countQ);
}

void Method606(int hiritu = 100, int minDiff = 10, int totyuu = 999, int _m = 30, int aketoku = 0)
{
  int M = _m;
  M = min(M, N);
  M = max(M, 10);
  vector<int> items;
  initializeItems(items);
  int countQ = 0;

  // 一部アイテムをマージソート(軽い順)
  MergeSort(items, countQ, 1001001, M);

  // ソートしていないアイテムを入れる
  int _D = D - aketoku;
  _D = max(_D, 2);
  for (int i = M; i < N; ++i) {
    int id = items[i];
    ans[id] = (i - M) % _D;
  }

  vector<int> groups;  // 常に重い順に並んでいるようにする
  initializeGroups(groups);
  // ここでグループを一度マージソートして正しい順序に並び替える
  MergeSort_Group(groups, countQ);
  // 一番軽いグループに入れていく
  int tail = N;
  while (tail > 0) {
    if (ans[items[tail - 1]] != -1) {
      tail--;
    }
    else {
      break;
    }
  }
  for (int i = tail - 1; i >= 0; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;

    // 初回だけ特別処理
    if (i == tail - 1) {
      for (int karina = 0; karina < 10; ++karina) {
        // 今のグループに数回Move1試す
        int dstId = Move1_Specify(countQ, gId);
        if (dstId != -1) {
          int dstPos = 0;
          for (int j = 0; j < D; ++j) {
            if (groups[j] == dstId) {
              dstPos = j;
            }
          }
          // 整列を更新
          int left = binarySearchGroupPosition(gId, groups, 0, dstPos, countQ);
          moveGroupToPosition(groups, dstPos, left);
        }
      }
    }

    int left = binarySearchGroupPosition(gId, groups, 0, D - 1, countQ);
    moveGroupToPosition(groups, D - 1, left);
  }
  int loop = 0;
  while (countQ < Q) {
    loop++;

    if (loop > 100000) {
      break;
    }
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      SwapNeighbor1(items, countQ, M);
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      Swap2(items, countQ, minDiff, M);
    }
  }
  completeWithDummyQueries(countQ);
}

void Method246(int hiritu, int totyuu, int small1, int small2)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  // Common initialization with pseudo query method
  initializeWithPseudoQueryMethod(items, groups, countQ, totyuu);

  // Find tail position and assign remaining items
  int tail = findTailPosition(items);
  assignRemainingItems(items, groups, countQ, tail);

  // Optimization loop with special small parameters
  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      int sma = small1;
      if (hiritu == 100) {
        if (rand_xorshift() % 20 == 0) {
          sma = 999;
        }
      }
      SwapNeighborSmall1(items, countQ, sma);
    }
    else {
      MoveSmall1(items, countQ, small2);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method112(int ikichi = N, int totyuu = 999)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  // Common initialization with pseudo query method (with custom ikichi)
  initializeWithPseudoQueryMethod(items, groups, countQ, totyuu, ikichi);

  // Assign ALL remaining items (Method112 doesn't use tail finding)
  assignRemainingItems(items, groups, countQ, N - D);

  // Optimization loop - simple SwapNeighbor1 only
  while (countQ < Q) {
    SwapNeighbor1(items, countQ);
  }

  completeWithDummyQueries(countQ);
}

void Method17(int ikichi, int hiritu)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeWithMergeSort(items, groups, countQ, ikichi);

  while (countQ < Q) {
    int qu = getRandomPercentage();
    if (qu < hiritu) {
      SwapNeighbor1(items, countQ);
    }
    else {
      if (isNearQueryLimit(countQ)) {
        if (hiritu < 10) {
          break;
        }
        continue;
      }
      Swap2(items, countQ, 999);
    }
  }

  completeWithDummyQueries(countQ);
}

void Method10(int hiritu = 70, int minDiff = 10, bool isMethod9 = false)
{
  vector<int> items;
  vector<int> groups;
  int countQ = 0;

  initializeItems(items);
  MergeSort(items, countQ);
  updateKarusaArray(items);

  assignHeaviestToGroups(items, N - 1);
  initializeGroups(groups);

  // pseudoItemsを用いて一番軽いグループに入れていく
  for (int i = N - D - 1; i >= 0; --i) {
    int id = items[i];
    int gId = groups[D - 1];
    ans[id] = gId;
    int left = 0;
    int right = D - 1;
    while (left < right) {
      int mid = (left + right) / 2;
      for (int j = 0; j < N; ++j) {
        if (ans[j] == gId) {
          l[countQ].push_back(j);
        }
        if (ans[j] == groups[mid]) {
          r[countQ].push_back(j);
        }
      }
      char c = PseudoItemsQuery(countQ);
      l[countQ].clear();
      r[countQ].clear();
      if (c == '=') {
        left = right = mid;
      }
      else if (c == '>') {
        right = mid;
      }
      else {
        left = mid + 1;
      }
    }
    moveGroupToPosition(groups, D - 1, left);
  }

  if (isMethod9) {
    while (countQ < Q) {
      int qu = rand_xorshift() % 100;
      if (qu < 100) {
        SwapNeighbor1(items, countQ);
      }
    }
  }
  else {
    while (countQ <= Q - 2) {
      int qu = rand_xorshift() % 100;
      if (qu < hiritu) {
        SwapNeighbor1(items, countQ);
      }
      else {
        Swap2(items, countQ, minDiff);
      }
    }
  }

  completeWithDummyQueries(countQ);
}

void PrintAns(ofstream& ofs)
{
  if (mode == 0) {
    for (int i = 0; i < N; ++i) { cout << ans[i] << ' '; }
    cout << endl;
  }
  else if (mode == 1) {
    ofs << "# " << haipara[NN][QQ][DD] << endl;
    for (int i = 0; i < Q; ++i) {
      ofs << "#c ";
      for (int j = 0; j < N; ++j) { ofs << answers[i][j] << ' '; }
      ofs << endl;
      ofs << comments[i] << endl;
      ofs << l[i].size() << ' ' << r[i].size();
      for (int j = 0; j < l[i].size(); ++j) { ofs << ' ' << l[i][j]; }
      for (int j = 0; j < r[i].size(); ++j) { ofs << ' ' << r[i][j]; }
      ofs << endl;
    }
    for (int i = 0; i < N; ++i) { ofs << ans[i] << ' '; }
    ofs << endl;

    for (int i = 0; i < N; ++i) {
      ofs << "# ";
      for (int j = 0; j < N; ++j) {
        if (hikaku[i][j] == -1) {
          ofs << "< ";
        }
        if (hikaku[i][j] == 0) {
          ofs << "= ";
        }
        if (hikaku[i][j] == 1) {
          ofs << "> ";
        }
        if (hikaku[i][j] == -2) {
          ofs << "? ";
        }
      }
      ofs << endl;
    }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  timer.start();
  elapsed = timer.get_elapsed_time();

  for (int i = 0; i < MAX_Q; ++i) {
    l[i].clear();
    r[i].clear();
    comments[i].clear();
  }
  comments[0] = "# ";
  for (int i = 0; i < 110; ++i) { ans[i] = -1; }
  for (int i = 0; i < 110; ++i) {
    for (int j = 0; j < 110; ++j) { hikaku[i][j] = -2; }
    hikaku[i][i] = 0;
  }
  calledMap.clear();
  queryCount = 0;
  query1Count = 0;
  ResetMemory();

  real_max_D = -1;
  real_minScore = INF;
}

ll Solve(int case_num, ll hai2 = D18)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(case_num);

  GeneratePseudoItems();

  int hai = 8;
  if (mode < 1000000) {
    hai = haipara[NN][QQ][DD];
    hai2 = haipara2[NN][QQ][DD];
  }
  else {
    hai = case_num;
  }
  hai %= 1000000;

  if (hai == 200208) {
    Method208(hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200206) {
    Method206(hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200216) {
    Method216(hai2 % D10 / D8, hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200226) {
    Method226(hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200227) {
    Method226(70, 10);
  }
  else if (hai == 200236) {
    Method216(hai2 % D10 / D8, hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200246) {
    Method246(hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200256) {
    Method246(100, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200266) {
    Method266(70, 10, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200366) {
    Method266(70, 10, hai2 % D4 / D2, hai2 % D2, hai2 % D7 / D4);
  }
  else if (hai == 200306) {
    Method306(hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200806) {
    Method806(hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200606) {
    Method606(hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200376) {
    Method306(hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2, hai2 % D10 / D8);
  }
  else if (hai == 200316) {
    Method316(hai2 % D10 / D8, hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200916) {
    Method916(hai2 % D12 / D10, hai2 % D10 / D8, hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200516) {
    Method516(hai2 % D10 / D8, hai2 % D8 / D6, hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (hai == 200706) {
    Method706(hai2 % D6 / D4, hai2 % D4 / D2, hai2 % D2);
  }
  else if (8000 <= hai && hai <= 8100) {
    Method8(hai - 8000);
  }
  else if (70000 <= hai && hai <= 79999) {
    Method6((hai - 70000) / 100, hai % 100);
  }
  else if (80000 <= hai && hai <= 89999) {
    Method17((hai - 80000) / 100, hai % 100);
  }
  else if (1200 <= hai && hai <= 1299) {
    Method12(hai - 1200);
  }
  else if (1900 <= hai && hai <= 1999) {
    Method19(hai - 1900);
  }
  else if (90000 <= hai && hai <= 99999) {
    Method106((hai - 90000) / 100, 10, hai % 100);
  }
  else if (100000 <= hai && hai <= 109999) {
    Method112((hai - 100000) / 100, hai % 100);
  }
  else if (110000 <= hai && hai <= 110100) {
    Method8_Two(hai - 110000);
  }
  else {
    switch (hai) {
      case 6:
        Method6();
        break;
      case 610:
        Method6(70, 10);
        break;
      case 61003:
        Method6(50, 3);
        break;
      case 8:
        Method8();
        break;
      case 9:
        Method10(70, 10, true);
        break;
      case 10:
        Method10();
        break;
      case 1003:
        Method10(50, 3);
        break;
      case 11:
        Method11();
        break;
      case 12:
        Method12();
        break;
      case 13:
        Method13();
        break;
      case 14:
        Method14();
        break;
      case 15:
        Method15();
        break;
      case 16:
        Method16();
        break;
      case 1220:
        Method12(20);
        break;
      default:
        cerr << "NG hai : " << hai << endl;
        Method8();
        break;
    }
  }

  if (!ErrorCheck()) {
    cerr << "ErrorCheck :  haipara = " << haipara[NN][QQ][DD] << ", haiapara2 = " << haipara2[NN][QQ][DD] << endl;
  }

  if (mode != 0) {
    elapsed = timer.get_elapsed_time();
    if (elapsed > 1.95) {
      cerr << "!!!TLE!!! << endl";
      cerr << "NN = " << NN << ", QQ = " << QQ << ", DD = " << DD << endl;
      cerr << "haipara = " << haipara[NN][QQ][DD] << ", haiapara2 = " << haipara2[NN][QQ][DD] << endl;
    }
  }

  // 出力ファイルストリームオープン
  ofstream ofs;
  if (mode == 0 || mode == 1) {
    OpenOfs(case_num, ofs);
    PrintAns(ofs);
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

void PrintHaipara(int loop)
{
  string str = "./haipara/haipara" + to_string(loop) + ".txt";
  std::ofstream file(str);
  file << "vector<int> haipara[14][40] = {" << endl;
  for (int i = 0; i < 14; ++i) {
    file << "{" << endl;
    int ii = i * 5 + 30 + 4;
    if (i == 13) ii = 100;
    for (int j = 0; j < 40; ++j) {
      file << "{";
      for (int k = 0; k < haipara[i][j].size(); ++k) {
        file << setw(7) << haipara[i][j][k];
        if (k < haipara[i][j].size() - 1) {
          file << ",";
        }
      }
      if (j < 39) {
        file << "}," << endl;
      }
      else {
        file << "}" << endl;
      }
    }
    if (i < 13) {
      file << "}," << endl;
    }
    else {
      file << "}" << endl;
    }
  }
  file << "};" << endl;
  file << endl;
  file << "vector<ll> haipara2[14][40] = {" << endl;
  for (int i = 0; i < 14; ++i) {
    file << "{" << endl;
    int ii = i * 5 + 30 + 4;
    if (i == 13) ii = 100;
    for (int j = 0; j < 40; ++j) {
      file << "{";
      for (int k = 0; k < haipara2[i][j].size(); ++k) {
        file << setw(19) << haipara2[i][j][k];
        if (k < haipara2[i][j].size() - 1) {
          file << ",";
        }
      }
      if (j < 39) {
        file << "}," << endl;
      }
      else {
        file << "}" << endl;
      }
    }
    if (i < 13) {
      file << "}," << endl;
    }
    else {
      file << "}" << endl;
    }
  }
  file << "};" << endl;
  file.close();
}

int main()
{
  for (int i = 0; i < 14; i++) {
    for (int j = 0; j < 40; j++) {
      haipara[i][j].clear();
      haipara[i][j] = Haipara[i][j];
      haipara2[i][j].clear();
      haipara2[i][j] = Haipara2[i][j];
    }
  }

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int j = 0; j < 1; ++j) {
      for (int i = 0; i < 100; ++i) {
        ll score = Solve(i);
        sum += score;
        if (score == 0) { continue; }
        cout << "num = " << setw(2) << i << ", ";
        cout << "N = " << N << ", Q = " << Q << ", D = " << D << ", queryCount = " << queryCount << ", ";
        cout << "nowTime = " << elapsed << ", ";
        cout << "haipara = " << haipara[NN][QQ][DD] << ", ";
        cout << "score = " << setw(7) << score << ", ";
        cout << "sum = " << setw(9) << sum << endl;
      }
    }
  }
  else if (mode == 2) {
    ll sum = 0;
    for (int ii = 0; ii < 100; ++ii) {
      int i = 209;
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }
  else if (mode == 3) {
    int loop = 0;
    while (true) {
      loop++;
      GenerateNNDDQQ();

      ll hai2 = haipara2[NN][QQ][DD];

      GeneratecaseFromNNDDQQ();

      Solve(2, hai2);
      elapsed = timer.get_elapsed_time();

      if (loop % 100 == 0) cout << "loop = " << loop << endl;
    }
  }
  else if (mode == 1000000) {
    int loop = 0;
    int winCount = 0;
    queue<P> winQueue;
    while (true) {
      int hai = 0;
      int newHai = 0;
      ll hai2 = D18;
      ll newHai2 = D18;
      if (winQueue.size()) {
        int newNNQQDD = winQueue.front().first;
        int challengeNNQQDD = winQueue.front().second;
        winQueue.pop();

        NN = challengeNNQQDD / 10000;
        QQ = challengeNNQQDD % 10000 / 100;
        DD = challengeNNQQDD % 100;
        hai = haipara[NN][QQ][DD];
        hai2 = haipara2[NN][QQ][DD];

        int NNN = newNNQQDD / 10000;
        int QQQ = newNNQQDD % 10000 / 100;
        int DDD = newNNQQDD % 100;
        newHai = haipara[NNN][QQQ][DDD];
        newHai2 = haipara2[NNN][QQQ][DDD];
      }
      else {
        GenerateNNDDQQ();

        hai = haipara[NN][QQ][DD];
        hai2 = haipara2[NN][QQ][DD];

        if (hai != 200316 && hai != 200916) {
          continue;
        }

        newHai = 200916;
        newHai2 = hai2 % D10 + rand_xorshift() % D2 * D10;
        if (hai == 200916 && rand_xorshift() % 2 == 0) {
          newHai2 = hai2 + D10;
        }

        if (false) {
          if (rand_xorshift() % 2 == 0) {
            int NNN = NN + rand_xorshift() % 3 - 1;
            NNN = max(NNN, 0);
            NNN = min(NNN, 13);
            int QQQ = QQ + rand_xorshift() % 3 - 1;
            QQQ = max(QQQ, 0);
            QQQ = min(QQQ, 39);
            int DDD = DD + rand_xorshift() % 3 - 1;
            DDD = max(0, DDD);
            DDD = min(DDD, (int)haipara[NNN][QQQ].size() - 1);
            newHai = haipara[NNN][QQQ][DDD];
            newHai2 = haipara2[NNN][QQQ][DDD];
            if (newHai != 200916) {
              continue;
            }
          }
        }
      }
      GeneratecaseFromNNDDQQ();

      if (newHai == hai && newHai2 == hai2) {
        continue;
      }
      bool isWin = true;
      int winSum = 0;
      int loseSum = 0;
      int drawSum = 0;
      for (int aespa = 0; aespa < 1; ++aespa) {
        int win = 0;
        int lose = 0;
        int draw = 0;

        int WinKijun = 0;
        if (WinKijun == 0) {
          for (int _ = 0; _ < 200; ++_) {
            GeneratecaseFromNNDDQQ();
            ll oldScore = Solve(hai, hai2);
            ll newScore = Solve(newHai, newHai2);
            if (newScore < oldScore) {
              win++;
            }
            else if (newScore == oldScore) {
              draw++;
            }
            else {
              lose++;
            }
            if (win - lose <= -10) { break; }
            if (win >= 120) { break; }
            if (win >= 60 && lose <= 30) { break; }
            if (lose >= 80) { break; }
            if (lose >= 20 && win <= lose) { break; }
            if (win >= 20 && lose * 11 / 10 >= win) { break; }
            if (win - lose >= 40) { break; }
          }
          winSum += win;
          loseSum += lose;
          drawSum += draw;
          if (win == 0 || (double)win / (win + lose) < 0.6) {
            isWin = false;
            break;
          }
          if (draw >= 100 && win <= 10) {
            isWin = false;
            break;
          }
        }
        else {
          for (int _ = 0; _ < 1000; ++_) {
            GeneratecaseFromNNDDQQ();
            ll oldScore = Solve(hai, hai2);
            ll newScore = Solve(newHai, newHai2);
            if (newScore < oldScore) {
              win++;
            }
            else if (newScore == oldScore) {
              draw++;
            }
            else {
              lose++;
            }
            if (win - lose <= -10) { break; }
            if (win >= 550) { break; }
            if (lose >= 20 && win <= lose) { break; }
            if (win - lose >= 100) { break; }
          }
          winSum += win;
          loseSum += lose;
          drawSum += draw;
          if (win <= 30 || (double)win / (win + lose) < 0.55) {
            isWin = false;
            break;
          }
        }
      }


      if (isWin) {
        if ((double)winSum / (winSum + loseSum) < 0.55) {
          isWin = false;
        }
      }

      if (isWin) {
        winCount++;
        haipara[NN][QQ][DD] = newHai;
        haipara2[NN][QQ][DD] = newHai2;
        cout << "loop = " << setw(5) << loop << ", ";
        cout << "N = " << setw(3) << N << ", ";
        cout << "D = " << setw(2) << D << ", ";
        cout << "Q = " << setw(4) << Q << ", ";
        cout << "Q / N = " << fixed << setprecision(2) << setw(5) << (double)Q / N << ", ";
        cout << "old = " << setw(7) << hai << ", new = " << setw(7) << newHai << ", ";
        cout << "oldHai2 = " << setw(19) << hai2 << ", ";
        cout << "newHai2 = " << setw(19) << newHai2;
        cout << " : " << winSum << "/" << loseSum << "/" << drawSum << endl;

        // winQueueに入れる
        for (int _ = 0; _ < 6; ++_) {
          int NNN = NN;
          int QQQ = QQ;
          int DDD = DD;
          if (_ == 0) {
            NNN = NN - 1;
            if (NNN < 0) { continue; }
          }
          else if (_ == 1) {
            NNN = NN + 1;
            if (NNN >= 14) { continue; }
          }
          else if (_ == 2) {
            QQQ = QQ - 1;
            if (QQQ < 0) { continue; }
          }
          else if (_ == 3) {
            QQQ = QQ + 1;
            if (QQQ >= 40) { continue; }
          }
          else if (_ == 4) {
            DDD = DD - 1;
            if (DDD < 0) { continue; }
          }
          else if (_ == 5) {
            DDD = DD + 1;
          }
          if (DDD >= haipara[NNN][QQQ].size()) { continue; }
          winQueue.push(
            P(NN * 10000 + QQ * 100 + DD, NNN * 10000 + QQQ * 100 + DDD));
        }
      }

      loop++;
      if (loop % 100 == 0) {
        PrintHaipara(loop);
        cout << "PrintHaipara : " << loop << endl;
        winCount = 0;
      }
    }
  }

  return 0;
}
