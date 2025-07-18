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

const int dx[5] = { 1, 0, -1, 0, 0 };
const int dy[5] = { 0, -1, 0, 1, 0 };
const char dirChar[5] = { 'D', 'L', 'U', 'R', '.' };
const char rotChar[3] = { 'L', '.', 'R' };
const char tipChar[2] = { '.', 'P' };
const int BASE_DIR = 3;
const int ddx[9] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };
const int ddy[9] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };

int order[5] = { 0,1,2,3,4 };
vector<int> ordVec[6] = { {0,-1,1},{0,1,-1},{-1,0,1},{-1,1,0},{1,-1,0},{1,0,-1} };

const double TL = 2.8;
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

const int MAX_N = 30;
const int MAX_M = 450;
const int MAX_V = 15;
const int MAX_T = 100005;

const double ACTION_RATIO_A = 50005.0;
const double ACTION_RATIO_B = 50004.0;
const double ACTION_RATIO_BASE = 50000.0;
const double POSITION_RATIO = 1.01;

int n, m, v;
int init_a[MAX_N][MAX_N];
int init_b[MAX_N][MAX_N];

int aRow[MAX_N];
int aColumn[MAX_N];
int bRow[MAX_N];
int bColumn[MAX_N];

int V;
int pa[MAX_V];
int le[MAX_V];
int ansCount;
int dir[MAX_T];
int rot[MAX_T][MAX_V];
int tip[MAX_T][MAX_V];
int sx;
int sy;
int leafs1;
int Method;
int ArmLengthMethod;
vector<int> PCount[100];

int best_V;
int best_pa[MAX_V];
int best_le[MAX_V];
int best_ansCount[100];
int best_total_ansCount;
int best_dir[MAX_T];
int best_rot[MAX_T][MAX_V];
int best_tip[MAX_T][MAX_V];
int best_sx;
int best_sy;
int best_leafs1;
int best_Method;
int best_ArmLengthMethod;
vector<int> best_PCount[100];

int route[MAX_N * MAX_N * 3][2];

double centerX;
double centerY;


void CopyToBest()
{
  best_V = V;
  for (int i = 0; i < V; ++i) {
    best_pa[i] = pa[i];
    best_le[i] = le[i];
  }

  best_total_ansCount = ansCount;
  for (int i = 0; i < ansCount; ++i) {
    best_dir[i] = dir[i];
    for (int j = 0; j < V; ++j) {
      best_rot[i][j] = rot[i][j];
      best_tip[i][j] = tip[i][j];
    }
  }
  best_sx = sx;
  best_sy = sy;
  best_leafs1 = leafs1;
  best_Method = Method;
  best_ArmLengthMethod = ArmLengthMethod;

  best_PCount[Method] = PCount[Method];
  best_ansCount[Method] = ansCount;
}

void CopyFromBest()
{
  V = best_V;
  for (int i = 0; i < V; ++i) {
    pa[i] = best_pa[i];
    le[i] = best_le[i];
  }
  ansCount = best_total_ansCount;
  for (int i = 0; i < ansCount; ++i) {
    dir[i] = best_dir[i];
    for (int j = 0; j < V; ++j) {
      rot[i][j] = best_rot[i][j];
      tip[i][j] = best_tip[i][j];
    }
  }
  sx = best_sx;
  sy = best_sy;
  leafs1 = best_leafs1;
  Method = best_Method;
  ArmLengthMethod = best_ArmLengthMethod;

  PCount[Method] = best_PCount[Method];
  //ansCount = best_ansCount[Method];
}

void CopyToBestPCount()
{
  best_ansCount[Method] = ansCount;
  best_PCount[Method] = PCount[Method];
}

bool IsNG(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  V = 0;
  ansCount = 0;
  for (int i = 0; i < 100; ++i) {
    PCount[i].clear();
    best_PCount[i].clear();
    best_ansCount[i] = 999;
  }
}

// 入力受け取り
void Input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> n >> m >> v;
    for (int i = 0; i < n; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < n; ++j) {
        init_a[i][j] = s[j] - '0';
      }
    }
    for (int i = 0; i < n; ++i) {
      string t;
      cin >> t;
      for (int j = 0; j < n; ++j) {
        init_b[i][j] = t[j] - '0';
      }
    }
  }
  else {
    // ファイル入力する
    ifs >> n >> m >> v;
    for (int i = 0; i < n; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < n; ++j) {
        init_a[i][j] = s[j] - '0';
      }
    }
    for (int i = 0; i < n; ++i) {
      string t;
      ifs >> t;
      for (int j = 0; j < n; ++j) {
        init_b[i][j] = t[j] - '0';
      }
    }
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
int CalcScore()
{
  return ansCount;
}

// NGチェック
bool IsValidAnswer()
{
  if (mode == 0) {
    return true;
  }

  vector<vector<int>> a, b;
  a.resize(n, vector<int>(n));
  b.resize(n, vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = init_a[i][j];
      b[i][j] = init_b[i][j];
    }
  }

  vector<int> nowRot(V);
  vector<int> nowTip(V);
  int x = sx;
  int y = sy;

  for (int t = 0; t < ansCount; ++t) {
    x = x + dx[dir[t]];
    y = y + dy[dir[t]];
    if (IsNG(x, y)) {
      cout << "Out of Range: turn = " << t << ", x = " << x << ", y = " << y << endl;
      return false;
    }

    for (int i = 1; i < V; ++i) {
      nowRot[i] = (nowRot[i] + rot[t][i]) % 4;
    }

    for (int i = 0; i < V; ++i) {
      if (tip[t][i] == 1) {
        // ここでアーム位置を計算
      }
    }
  }
  return true;
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << V << endl;
    for (int i = 1; i < V; ++i) {
      cout << pa[i] << ' ' << le[i] << endl;
    }
    cout << sx << ' ' << sy << endl;
    for (int i = 0; i < ansCount; ++i) {
      cout << dirChar[dir[i]];
      for (int j = 1; j < v; ++j)cout << rotChar[rot[i][j] + 1]; // rotは-1〜1のため帳尻合わせ
      for (int j = 0; j < V; ++j)cout << tipChar[tip[i][j]];
      cout << endl;
    }
  }
  else {
    ofs << V << endl;
    for (int i = 1; i < V; ++i) {
      ofs << pa[i] << ' ' << le[i] << endl;
    }
    ofs << sx << ' ' << sy << endl;
    for (int i = 0; i < ansCount; ++i) {
      ofs << dirChar[dir[i]];
      for (int j = 1; j < v; ++j)ofs << rotChar[rot[i][j] + 1]; // rotは-1〜1のため帳尻合わせ
      for (int j = 0; j < V; ++j)ofs << tipChar[tip[i][j]];
      ofs << endl;
    }
  }
}

void CopyAB(vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  mCount = 0;
  centerX = 0;
  centerY = 0;
  a.resize(n, vector<int>(n));
  b.resize(n, vector<int>(n));
  for (int i = 0; i < n; ++i) {
    aRow[i] = 0;
    aColumn[i] = 0;
    bRow[i] = 0;
    bColumn[i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j] = init_a[i][j];
      b[i][j] = init_b[i][j];
      if (a[i][j] == 1 && b[i][j] == 1) {
        mCount++;
        a[i][j] = 0;
        b[i][j] = 0;
      }
      aRow[i] += a[i][j];
      aColumn[j] += a[i][j];
      bRow[i] += b[i][j];
      bColumn[j] += b[i][j];

      centerX += a[i][j] * i + b[i][j] * i;
      centerY += a[i][j] * j + b[i][j] * j;
    }
  }

  centerX /= 2.0 * (m - mCount);
  centerY /= 2.0 * (m - mCount);
}

class RotTip
{
private:
  int m_V;

public:
  vector<int> Rot;
  vector<int> Tip;
  vector<int> NowRot;
  vector<int> NowTip;

  RotTip()
  {
    m_V = v;
    Rot.resize(v);
    Tip.resize(v);
    NowRot.resize(v);
    NowTip.resize(v);
  }

  RotTip(int v)
  {
    m_V = v;
    Rot.resize(v);
    Tip.resize(v);
    NowRot.resize(v);
    NowTip.resize(v);
  }

  void Initialize(const vector<int>& nowRot, const vector<int>& nowTip)
  {
    for (int i = 0; i < m_V; ++i) {
      Rot[i] = 0;
      Tip[i] = 0;
      NowRot[i] = nowRot[i];
      NowTip[i] = nowTip[i];
    }
  }

  void Reflect(vector<int>& nowRot, vector<int>& nowTip, const int _t) const
  {
    for (int i = 0; i < V; ++i) {
      rot[_t][i] = Rot[i];
      tip[_t][i] = Tip[i];
      nowRot[i] = NowRot[i];
      nowTip[i] = NowTip[i];
    }
  }
};

class KeepAB
{
public:
  int KeepA[MAX_V][3];
  int KeepB[MAX_V][3];
  int KeepACount;
  int KeepBCount;

  KeepAB()
  {
    KeepACount = 0;
    KeepBCount = 0;
  }

  void AddA(int x, int y, int val)
  {
    KeepA[KeepACount][0] = x;
    KeepA[KeepACount][1] = y;
    KeepA[KeepACount][2] = val;
    KeepACount++;
  }

  void AddB(int x, int y, int val)
  {
    KeepB[KeepBCount][0] = x;
    KeepB[KeepBCount][1] = y;
    KeepB[KeepBCount][2] = val;
    KeepBCount++;
  }

  void Clear()
  {
    KeepACount = 0;
    KeepBCount = 0;
  }

  void Copy(const KeepAB& src)
  {
    KeepACount = src.KeepACount;
    for (int i = 0; i < KeepACount; ++i) {
      for (int j = 0; j < 3; ++j) {
        KeepA[i][j] = src.KeepA[i][j];
      }
    }
    KeepBCount = src.KeepBCount;
    for (int i = 0; i < KeepBCount; ++i) {
      for (int j = 0; j < 3; ++j) {
        KeepB[i][j] = src.KeepB[i][j];
      }
    }
  }
};

// 配列の要素を復元して行・列のカウントを更新する共通処理
void RestoreArrayWithCount(vector<vector<int>>& arr, int row[], int column[],
  const int keep[][3], int count)
{
  for (int i = 0; i < count; ++i) {
    arr[keep[i][0]][keep[i][1]] = keep[i][2];
    row[keep[i][0]]++;
    column[keep[i][1]]++;
  }
}

void RollBackFromKeepAB(const KeepAB& keepAB, vector<vector<int>>& a, vector<vector<int>>& b)
{
  RestoreArrayWithCount(a, aRow, aColumn, keepAB.KeepA, keepAB.KeepACount);
  RestoreArrayWithCount(b, bRow, bColumn, keepAB.KeepB, keepAB.KeepBCount);
}

// 配列の要素をクリアして行・列のカウントを更新する共通処理
void ClearArrayWithCount(vector<vector<int>>& arr, int row[], int column[],
  const int keep[][3], int count)
{
  for (int i = 0; i < count; ++i) {
    if (arr[keep[i][0]][keep[i][1]] == 0) {
      assert(false);
    }
    arr[keep[i][0]][keep[i][1]] = 0;
    row[keep[i][0]]--;
    column[keep[i][1]]--;
  }
}

void ReflectFromMaxAB(const KeepAB& maxAB, vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  mCount += maxAB.KeepBCount;
  ClearArrayWithCount(a, aRow, aColumn, maxAB.KeepA, maxAB.KeepACount);
  ClearArrayWithCount(b, bRow, bColumn, maxAB.KeepB, maxAB.KeepBCount);
}

// 配列の要素を1に設定して行・列のカウントを更新する共通処理
void SetArrayToOneWithCount(vector<vector<int>>& arr, int row[], int column[],
  const int keep[][3], int count)
{
  for (int i = 0; i < count; ++i) {
    if (arr[keep[i][0]][keep[i][1]] == 1) {
      assert(false);
    }
    arr[keep[i][0]][keep[i][1]] = 1;
    row[keep[i][0]]++;
    column[keep[i][1]]++;
  }
}

void RollBackFromMaxAB(const KeepAB& maxAB, vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  mCount -= maxAB.KeepBCount;
  SetArrayToOneWithCount(a, aRow, aColumn, maxAB.KeepA, maxAB.KeepACount);
  SetArrayToOneWithCount(b, bRow, bColumn, maxAB.KeepB, maxAB.KeepBCount);
}

class MaxCandidate
{
public:
  double maxActionScore;
  int finishTurn;
  int maxMarginCount;
  int maxDir;
  RotTip maxRT;
  KeepAB maxAB;

  MaxCandidate()
  {
    maxActionScore = -1;
    finishTurn = 999;
    maxMarginCount = 999;
    maxDir = -1;
  }
};

// 周囲のセルをチェックしてスコアを計算する共通処理
void CalculatePositionBonus(double& actionScore, int newX, int newY,
  const vector<vector<int>>& arr)
{
  for (int k = 0; k < 4; ++k) {
    for (int l = 1; l < 2; ++l) {
      int nkrx = newX + dx[k] * l;
      int nkry = newY + dy[k] * l;
      if (IsNG(nkrx, nkry)) {
        actionScore += POSITION_RATIO;
      }
      else if (arr[nkrx][nkry] == 0) {
        actionScore += POSITION_RATIO;
      }
      else {
        break;
      }
    }
  }
}

bool CanCatch(const vector<int>& nowRot, const vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b,
  RotTip& tmpRT, KeepAB& keepAB, double& actionScore,
  const int i, const int j, const int newX, const int newY)
{
  actionScore = 0;

  if (nowTip[i] == 0) {
    if (a[newX][newY] == 1) {
      keepAB.AddA(newX, newY, a[newX][newY]);
      a[newX][newY] = 0;
      aRow[newX]--;
      aColumn[newY]--;
      tmpRT.Rot[i] = j;
      tmpRT.NowRot[i] = (nowRot[i] + j) % 4;
      tmpRT.Tip[i] = 1;
      tmpRT.NowTip[i] = 1;
      actionScore = ACTION_RATIO_A;
      actionScore += abs(centerX - newX);
      actionScore += abs(centerY - newY);
      //actionScore += 10 * (20 - le[i]);

      CalculatePositionBonus(actionScore, newX, newY, a);

      //for (int l = 1; l < 10; ++l)
      //{
      //  int nkrx = newX + l;
      //  if (IsNG(nkrx, newY))break;
      //  if (a[nkrx][newY] == 0)break;
      //  nkrx = newX - l;
      //  if (IsNG(nkrx, newY))break;
      //  if (a[nkrx][newY] == 0)break;
      //  actionScore -= POSITION_RATIO * 3453;
      //}
      //for (int l = 1; l < 10; ++l)
      //{
      //  int nkry = newY + l;
      //  if (IsNG(newX, nkry))break;
      //  if (a[newX][nkry] == 0)break;
      //  nkry = newY - l;
      //  if (IsNG(newX, nkry))break;
      //  if (a[newX][nkry] == 0)break;
      //  actionScore -= POSITION_RATIO * 3453;
      //}

      //if (aRow[newX] == 0)actionScore += POSITION_RATIO * 100;
      //if (aColumn[newY] == 0)actionScore += POSITION_RATIO * 100;
      return true;
    }
  }
  else {
    if (b[newX][newY] == 1) {
      keepAB.AddB(newX, newY, b[newX][newY]);
      b[newX][newY] = 0;
      bRow[newX]--;
      bColumn[newY]--;
      tmpRT.Rot[i] = j;
      tmpRT.NowRot[i] = (nowRot[i] + j) % 4;
      tmpRT.Tip[i] = 1;
      tmpRT.NowTip[i] = 0;
      actionScore = ACTION_RATIO_B;
      actionScore += abs(centerX - newX);
      actionScore += abs(centerY - newY);
      //actionScore += 10 * (20 - le[i]);

      CalculatePositionBonus(actionScore, newX, newY, b);

      //if (bRow[newX] == 0)actionScore += POSITION_RATIO * 100;
      //if (bColumn[newY] == 0)actionScore += POSITION_RATIO * 100;
      return true;
    }
  }

  return false;
}

int MakeLength(int rand)
{
  int length = 1;
  if (rand < 30) {
    ArmLengthMethod = 0;
    length = Rand() % ((n - 1) * 1 / 3) + 1;
  }
  else if (rand < 60) {
    ArmLengthMethod = 1;
    length = Rand() % ((n - 1) * 1 / 2) + 1;
  }
  else if (rand < 80) {
    ArmLengthMethod = 2;
    length = Rand() % ((n - 1) * 2 / 3) + 1;
  }
  else if (rand < 100) {
    ArmLengthMethod = 3;
    length = Rand() % ((n - 1) * 3 / 4) + 1;
  }
  else {
    ArmLengthMethod = 4;
    length = Rand() % (n - 1) + 1;
  }
  return length;
}

// 木構造のパレント設定共通処理
void SetTreeParents(int* parentRules, int ruleCount)
{
  for (int i = 1; i < V && i <= ruleCount; i++) {
    pa[i] = parentRules[i - 1];
  }
  for (int i = ruleCount + 1; i < V; i++) {
    pa[i] = parentRules[ruleCount - 1];
  }
}

// 木構造の長さ設定共通処理
void SetTreeLengths(int rand, int startIdx, int baseOffset, int range)
{
  for (int i = 1; i < V; ++i) {
    le[i] = MakeLength(rand);
  }

  if (Rand() % 2) {
    int st = Rand() % range + 1;
    for (int i = startIdx; i < V; ++i)le[i] = i - startIdx + st;
  }
}

void MakeTree1()
{
  int rand = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1 };
  SetTreeParents(parentRules, 2);
  SetTreeLengths(rand, 2, 2, 2);
}

void MakeTree2()
{
  int rand = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2 };
  SetTreeParents(parentRules, 3);
  SetTreeLengths(rand, 3, 3, 3);
}

void MakeTree6()
{
  int rand = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3 };
  SetTreeParents(parentRules, 4);
  SetTreeLengths(rand, 4, 4, 4);
}

void MakeTree4()
{
  int rand = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3, 4 };
  SetTreeParents(parentRules, 5);
  SetTreeLengths(rand, 5, 5, 5);
}

void MakeTree5()
{
  int rand = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3, 4, 5 };
  SetTreeParents(parentRules, 6);
  SetTreeLengths(rand, 6, 6, 5);
}

void MakeTree22()
{
  int rand = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 6) + 1;
  for (int i = 1; i < V; ++i) {
    if (i == 1) {
      pa[i] = 0;
    }
    else if (i == 2) {
      pa[i] = 1;
    }
    else if (i == 3) {
      pa[i] = 1;
    }
    else if (i < 4 + leafs1) {
      pa[i] = 2;
    }
    else {
      pa[i] = 3;
    }
    le[i] = MakeLength(rand);
  }

  if (Rand() % 2) {
    int st = Rand() % 4 + 1;
    for (int i = 4; i < 4 + leafs1; ++i)le[i] = i - 4 + st;
    st = Rand() % 4 + 1;
    for (int i = 4 + leafs1; i < V; ++i)le[i] = i - (4 + leafs1) + st;
  }
}

void MakeTree32()
{
  int rand = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 7) + 1;
  for (int i = 1; i < V; ++i) {
    if (i == 1) {
      pa[i] = 0;
    }
    else if (i == 2) {
      pa[i] = 1;
    }
    else if (i == 3) {
      pa[i] = 2;
    }
    else if (i == 4) {
      pa[i] = 2;
    }
    else if (i < 5 + leafs1) {
      pa[i] = 3;
    }
    else {
      pa[i] = 4;
    }
    le[i] = MakeLength(rand);
  }

  if (Rand() % 2) {
    int st = Rand() % 5 + 1;
    for (int i = 5; i < 5 + leafs1; ++i)le[i] = i - 5 + st;
    st = Rand() % 5 + 1;
    for (int i = 5 + leafs1; i < V; ++i)le[i] = i - (5 + leafs1) + st;
  }
}

void MakeTree42()
{
  int rand = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 8) + 1;
  for (int i = 1; i < V; ++i) {
    if (i == 1) {
      pa[i] = 0;
    }
    else if (i == 2) {
      pa[i] = 1;
    }
    else if (i == 3) {
      pa[i] = 2;
    }
    else if (i == 4) {
      pa[i] = 3;
    }
    else if (i == 5) {
      pa[i] = 3;
    }
    else if (i < 6 + leafs1) {
      pa[i] = 4;
    }
    else {
      pa[i] = 5;
    }
    le[i] = MakeLength(rand);
  }

  if (Rand() % 2) {
    int st = Rand() % 5 + 1;
    for (int i = 6; i < 6 + leafs1; ++i)le[i] = i - 6 + st;
    st = Rand() % 5 + 1;
    for (int i = 6 + leafs1; i < V; ++i)le[i] = i - (6 + leafs1) + st;
  }
}

int CalcNeedLength(const int baseX, const int baseY)
{
  int needLength = 999;
  if (baseX < 0) {
    needLength = 0 - baseX;
  }
  else if (baseX < n) {
    needLength = 0;
  }
  else {
    needLength = baseX - (n - 1);
  }
  if (baseY < 0) {
    needLength = min(needLength, 0 - baseY);
  }
  else if (baseX < n) {
    needLength = 0;
  }
  else {
    needLength = min(needLength, baseY - (n - 1));
  }

  return needLength;
}

int CalcMarginCount(const int leafCount, const RotTip& maxRT)
{
  int marginCount = leafCount;
  for (int i = 0; i < V; ++i) {
    if (maxRT.Tip[i] == 1) {
      marginCount--;
    }
  }
  return marginCount;
}

// CallDecideBestの前方宣言
void CallDecideBest(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b);

int doMarginCount = 0;
int doOneSetCount = 0;
int doRowColumnCount = 0;
double DoOneSet(RotTip& tmpRT, KeepAB& keepAB,
  const int baseX, const int baseY, const int needLength, const int prenRot, const int startLeaf,
  const vector<int>& nowRot, const vector<int>& nowTip,
  const MaxCandidate& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{

  doOneSetCount++;

  if (baseX < 0 || n <= baseX || (aRow[baseX] == 0 && bRow[baseX] == 0)) {
    if (baseY < 0 || n <= baseY || (aColumn[baseY] == 0 && bColumn[baseY] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  int margin = candidates.maxMarginCount;
  double actionScore = 0;
  double actionScoreSum = 0;
  keepAB.Clear();

  for (int i = startLeaf; i < V; ++i) {
    if (le[i] < needLength)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j) {
      int newRot = (prenRot + nowRot[i] + j + 4) % 4;
      int newX = baseX + le[i] * dx[newRot];
      int newY = baseY + le[i] * dy[newRot];

      if (IsNG(newX, newY))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, newX, newY);
      if (isCatch)break;
    }

    actionScoreSum += actionScore;

    if (!isCatch) {
      margin--;
    }

    if (margin < 0) {
      doMarginCount++;
      break;
    }
  }

  RollBackFromKeepAB(keepAB, a, b);

  return actionScoreSum;
}

double DoOneSet2(RotTip& tmpRT, KeepAB& keepAB,
  const int baseX1, const int baseY1, const int needLength1, const int prenRot1,
  const int baseX2, const int baseY2, const int needLength2, const int prenRot2,
  const int startLeaf, const vector<int>& nowRot, const vector<int>& nowTip,
  const MaxCandidate& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{

  doOneSetCount++;

  if (baseX1 < 0 || n <= baseX1 || (aRow[baseX1] == 0 && bRow[baseX1] == 0)) {
    if (baseY1 < 0 || n <= baseY1 || (aColumn[baseY1] == 0 && bColumn[baseY1] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  if (baseX2 < 0 || n <= baseX2 || (aRow[baseX2] == 0 && bRow[baseX2] == 0)) {
    if (baseY2 < 0 || n <= baseY2 || (aColumn[baseY2] == 0 && bColumn[baseY2] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  int margin = candidates.maxMarginCount;
  double actionScore = 0;
  double actionScoreSum = 0;
  keepAB.Clear();

  for (int i = startLeaf; i < startLeaf + leafs1; ++i) {
    if (le[i] < needLength1)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j) {
      int newRot = (prenRot1 + nowRot[i] + j + 4) % 4;
      int newX = baseX1 + le[i] * dx[newRot];
      int newY = baseY1 + le[i] * dy[newRot];

      if (IsNG(newX, newY))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, newX, newY);
      if (isCatch)break;
    }

    actionScoreSum += actionScore;

    if (!isCatch) {
      margin--;
    }

    if (margin < 0) {
      doMarginCount++;
      break;
    }
  }

  for (int i = startLeaf + leafs1; i < V; ++i) {
    if (margin < 0) {
      doMarginCount++;
      break;
    }

    if (le[i] < needLength2)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j) {
      int newRot = (prenRot2 + nowRot[i] + j + 4) % 4;
      int newX = baseX2 + le[i] * dx[newRot];
      int newY = baseY2 + le[i] * dy[newRot];

      if (IsNG(newX, newY))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, newX, newY);
      if (isCatch)break;
    }

    actionScoreSum += actionScore;

    if (!isCatch) {
      margin--;
    }
  }

  RollBackFromKeepAB(keepAB, a, b);

  return actionScoreSum;
}

void SetCandidate(vector<MaxCandidate>& candidates, const int ordDir, const double actionScoreSum, const RotTip& tmpRT, const KeepAB& keepAB, const int startLeaf)
{
  int beamWidth = candidates.size();

  candidates[beamWidth - 1].maxDir = ordDir;
  candidates[beamWidth - 1].maxActionScore = actionScoreSum;
  candidates[beamWidth - 1].maxRT = tmpRT;
  candidates[beamWidth - 1].maxAB.Copy(keepAB);
  candidates[beamWidth - 1].maxMarginCount = CalcMarginCount(V - startLeaf, candidates[beamWidth - 1].maxRT);

  int idx = beamWidth - 1;
  while (idx > 0 && candidates[idx].maxActionScore > candidates[idx - 1].maxActionScore) {
    swap(candidates[idx], candidates[idx - 1]);
    idx--;
  }
}

void DecideBest42(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      tmpRT.Initialize(nowRot, nowTip);
      tmpRT.Rot[1] = ii;
      tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;

      int baseX = nx + le[1] * dx[newRot1];
      int baseY = ny + le[1] * dy[newRot1];

      int needLength = CalcNeedLength(baseX, baseY);

      double actionScoreSum = DoOneSet(tmpRT, keepAB, baseX, baseY, needLength, newRot1,
        2, nowRot, nowTip, candidates[beamWidth - 1], a, b);

      if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
        SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 2);
      }
    }
  }
}

void DecideBest52(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[2]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        tmpRT.Initialize(nowRot, nowTip);
        tmpRT.Rot[1] = ii;
        tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
        tmpRT.Rot[2] = iii;
        tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;

        int baseX = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2];
        int baseY = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2];

        int needLength = CalcNeedLength(baseX, baseY);

        double actionScoreSum = DoOneSet(tmpRT, keepAB, baseX, baseY, needLength, newRot2,
          3, nowRot, nowTip, candidates[beamWidth - 1], a, b);

        if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
          SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 3);
        }
      }
    }
  }
}

void DecideBest62(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[2]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[4]]) {
          int newRot3 = (newRot2 + nowRot[3] + iiii + 4) % 4;

          tmpRT.Initialize(nowRot, nowTip);
          tmpRT.Rot[1] = ii;
          tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
          tmpRT.Rot[2] = iii;
          tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
          tmpRT.Rot[3] = iiii;
          tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;

          int baseX = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3];
          int baseY = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3];

          int needLength = CalcNeedLength(baseX, baseY);

          double actionScoreSum = DoOneSet(tmpRT, keepAB, baseX, baseY, needLength, newRot3,
            4, nowRot, nowTip, candidates[beamWidth - 1], a, b);

          if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
            SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 4);
          }
        }
      }
    }
  }
}

void DecideBest72(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[2]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[4]]) {
          int newRot3 = (newRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[rand[5]]) {
            int newRot4 = (newRot3 + nowRot[4] + iiiii + 4) % 4;

            tmpRT.Initialize(nowRot, nowTip);
            tmpRT.Rot[1] = ii;
            tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
            tmpRT.Rot[2] = iii;
            tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
            tmpRT.Rot[3] = iiii;
            tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
            tmpRT.Rot[4] = iiiii;
            tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;

            int baseX = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3] + le[4] * dx[newRot4];
            int baseY = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3] + le[4] * dy[newRot4];

            int needLength = CalcNeedLength(baseX, baseY);

            double actionScoreSum = DoOneSet(tmpRT, keepAB, baseX, baseY, needLength, newRot4,
              5, nowRot, nowTip, candidates[beamWidth - 1], a, b);

            if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
              SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 5);
            }
          }
        }
      }
    }
  }
}

void DecideBest82(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[1]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[2]]) {
          int newRot3 = (newRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[rand[3]]) {
            int newRot4 = (newRot3 + nowRot[4] + iiiii + 4) % 4;

            for (const auto iiiiii : ordVec[rand[4]]) {
              int newRot5 = (newRot4 + nowRot[5] + iiiiii + 4) % 4;

              tmpRT.Initialize(nowRot, nowTip);
              tmpRT.Rot[1] = ii;
              tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
              tmpRT.Rot[2] = iii;
              tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
              tmpRT.Rot[3] = iiii;
              tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
              tmpRT.Rot[4] = iiiii;
              tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;
              tmpRT.Rot[5] = iiiiii;
              tmpRT.NowRot[5] = (nowRot[5] + iiiiii) % 4;

              int baseX = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3] + le[4] * dx[newRot4] + le[5] * dx[newRot5];
              int baseY = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3] + le[4] * dy[newRot4] + le[5] * dy[newRot5];

              int needLength = CalcNeedLength(baseX, baseY);

              double actionScoreSum = DoOneSet(tmpRT, keepAB, baseX, baseY, needLength, newRot5,
                6, nowRot, nowTip, candidates[beamWidth - 1], a, b);

              if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
                SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 6);
              }
            }
          }
        }
      }
    }
  }
}


void DecideBest53(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[2]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[4]]) {
          int newRot3 = (newRot1 + nowRot[3] + iiii + 4) % 4;

          tmpRT.Initialize(nowRot, nowTip);
          tmpRT.Rot[1] = ii;
          tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
          tmpRT.Rot[2] = iii;
          tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
          tmpRT.Rot[3] = iiii;
          tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;

          int baseX1 = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2];
          int baseY1 = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2];
          int baseX2 = nx + le[1] * dx[newRot1] + le[3] * dx[newRot3];
          int baseY2 = ny + le[1] * dy[newRot1] + le[3] * dy[newRot3];

          int needLength1 = CalcNeedLength(baseX1, baseY1);
          int needLength2 = CalcNeedLength(baseX2, baseY2);

          double actionScoreSum = DoOneSet2(tmpRT, keepAB,
            baseX1, baseY1, needLength1, newRot2,
            baseX2, baseY2, needLength2, newRot3,
            4, nowRot, nowTip, candidates[beamWidth - 1], a, b);

          if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
            SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 4);
          }
        }
      }
    }
  }
}

void DecideBest63(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[2]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[4]]) {
          int newRot3 = (newRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[rand[5]]) {
            int newRot4 = (newRot2 + nowRot[4] + iiiii + 4) % 4;

            tmpRT.Initialize(nowRot, nowTip);
            tmpRT.Rot[1] = ii;
            tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
            tmpRT.Rot[2] = iii;
            tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
            tmpRT.Rot[3] = iiii;
            tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
            tmpRT.Rot[4] = iiiii;
            tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;

            int baseX1 = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3];
            int baseY1 = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3];
            int baseX2 = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[4] * dx[newRot4];
            int baseY2 = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[4] * dy[newRot4];

            int needLength1 = CalcNeedLength(baseX1, baseY1);
            int needLength2 = CalcNeedLength(baseX2, baseY2);

            double actionScoreSum = DoOneSet2(tmpRT, keepAB,
              baseX1, baseY1, needLength1, newRot3,
              baseX2, baseY2, needLength2, newRot4,
              5, nowRot, nowTip, candidates[beamWidth - 1], a, b);

            if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
              SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 5);
            }
          }
        }
      }
    }
  }
}

void DecideBest73(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = candidates.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord) {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int rand[6];
    for (int i = 0; i < 6; ++i)rand[i] = Rand() % 6;
    for (const auto ii : ordVec[rand[0]]) {
      int newRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[rand[1]]) {
        int newRot2 = (newRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[rand[2]]) {
          int newRot3 = (newRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[rand[3]]) {
            int newRot4 = (newRot3 + nowRot[4] + iiiii + 4) % 4;

            for (const auto iiiiii : ordVec[rand[4]]) {
              int newRot5 = (newRot3 + nowRot[5] + iiiiii + 4) % 4;

              tmpRT.Initialize(nowRot, nowTip);
              tmpRT.Rot[1] = ii;
              tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
              tmpRT.Rot[2] = iii;
              tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
              tmpRT.Rot[3] = iiii;
              tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
              tmpRT.Rot[4] = iiiii;
              tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;
              tmpRT.Rot[5] = iiiiii;
              tmpRT.NowRot[5] = (nowRot[5] + iiiiii) % 4;

              int baseX1 = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3] + le[4] * dx[newRot4];
              int baseY1 = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3] + le[4] * dy[newRot4];
              int baseX2 = nx + le[1] * dx[newRot1] + le[2] * dx[newRot2] + le[3] * dx[newRot3] + le[5] * dx[newRot5];
              int baseY2 = ny + le[1] * dy[newRot1] + le[2] * dy[newRot2] + le[3] * dy[newRot3] + le[5] * dy[newRot5];

              int needLength1 = CalcNeedLength(baseX1, baseY1);
              int needLength2 = CalcNeedLength(baseX2, baseY2);

              double actionScoreSum = DoOneSet2(tmpRT, keepAB,
                baseX1, baseY1, needLength1, newRot4,
                baseX2, baseY2, needLength2, newRot5,
                6, nowRot, nowTip, candidates[beamWidth - 1], a, b);

              if (actionScoreSum > candidates[beamWidth - 1].maxActionScore) {
                SetCandidate(candidates, order[ord], actionScoreSum, tmpRT, keepAB, 6);
              }
            }
          }
        }
      }
    }
  }
}

bool UpdateTurn(const MaxCandidate& candidates, vector<int>& nowRot, vector<int>& nowTip, int& x, int& y, int& _t,
  vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  bool isOk = true;

  candidates.maxRT.Reflect(nowRot, nowTip, _t);

  dir[_t] = candidates.maxDir;

  ReflectFromMaxAB(candidates.maxAB, a, b, mCount);
  PCount[Method][_t] += candidates.maxAB.KeepACount + candidates.maxAB.KeepBCount;

  if (Method == best_Method) {
    if (best_PCount[Method].size() > 0 && _t >= best_PCount[Method].size()) {
      isOk = false;
    }
    else if (PCount[Method][_t] < best_PCount[Method][_t] - 5) {
      isOk = false;
    }
  }

  x += dx[candidates.maxDir];
  y += dy[candidates.maxDir];
  _t++;

  PCount[Method].push_back(PCount[Method][_t - 1]);

  return isOk;
}

void RollBackTurn(const MaxCandidate& candidates, int& x, int& y, int& _t,
  vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  RollBackFromMaxAB(candidates.maxAB, a, b, mCount);
  PCount[Method].pop_back();

  x -= dx[candidates.maxDir];
  y -= dy[candidates.maxDir];
  _t--;
}

MaxCandidate Beam(int& _t, int& x, int& y, vector<int>& nowRot, vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b, const int beamDepth, int& mCount)
{
  FisherYates(order, 5);

  int BEAM_WIDTH = 3;
  if (beamDepth == 0) {
    BEAM_WIDTH = 1;
  }

  vector<MaxCandidate> candidates(BEAM_WIDTH);
  for (int i = 0; i < BEAM_WIDTH; ++i) {
    candidates[i].maxRT.Initialize(nowRot, nowTip);
  }

  CallDecideBest(x, y, nowRot, nowTip, candidates, a, b);

  if (mCount + candidates[0].maxAB.KeepBCount == m) {
    candidates[0].finishTurn = 0;
    return candidates[0];
  }

  if (beamDepth == 0) {
    return candidates[0];
  }

  auto keepNowRot2 = nowRot;
  auto keepNowTip2 = nowTip;


  int minFinishTurn = 999;
  double maxActionScore = -1;
  int bestIndex = -1;
  for (int idx = 0; idx < BEAM_WIDTH; ++idx) {
    auto keepNowRot = nowRot;
    auto keepNowTip = nowTip;
    UpdateTurn(candidates[idx], nowRot, nowTip, x, y, _t, a, b, mCount);

    MaxCandidate tmpMaxCand = Beam(_t, x, y, nowRot, nowTip, a, b, beamDepth - 1, mCount);

    nowRot = keepNowRot;
    nowTip = keepNowTip;
    RollBackTurn(candidates[idx], x, y, _t, a, b, mCount);

    candidates[idx].finishTurn = min(tmpMaxCand.finishTurn + 1, 999);
    candidates[idx].maxActionScore += tmpMaxCand.maxActionScore;

    if (candidates[idx].finishTurn < minFinishTurn) {
      minFinishTurn = candidates[idx].finishTurn;
      maxActionScore = candidates[idx].maxActionScore;
      bestIndex = idx;
    }
    else if (candidates[idx].finishTurn == minFinishTurn && candidates[idx].maxActionScore > maxActionScore) {
      minFinishTurn = candidates[idx].finishTurn;
      maxActionScore = candidates[idx].maxActionScore;
      bestIndex = idx;
    }
  }

  return candidates[bestIndex];
}

// DecideBest関数を呼び出す共通処理
void CallDecideBest(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& candidates, vector<vector<int>>& a, vector<vector<int>>& b)
{
  if (Method == 42) {
    DecideBest42(x, y, nowRot, nowTip, candidates, a, b);
  }
  else if (Method == 52) {
    DecideBest52(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 62) {
    DecideBest62(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 72) {
    DecideBest72(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 82) {
    DecideBest82(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 53) {
    DecideBest53(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 63) {
    DecideBest63(x, y, nowRot, nowTip, candidates, a, b);
  }
  else  if (Method == 73) {
    DecideBest73(x, y, nowRot, nowTip, candidates, a, b);
  }
  else {
    if (mode != 0) {
      cout << "NG" << endl;
    }
  }
}

void Beam2_Internal(int& _t, int& x, int& y, vector<int>& nowRot, vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b, const int beamDepth, int& mCount, const int beamWidth, vector<MaxCandidate>& candidates)
{
  FisherYates(order, 5);
  CallDecideBest(x, y, nowRot, nowTip, candidates, a, b);
}

MaxCandidate Beam2(int& _t, int& x, int& y, vector<int>& nowRot, vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b, const int beamDepth, int& mCount)
{
  int BEAM_WIDTH = 30;

  vector<MaxCandidate> candidates(BEAM_WIDTH);
  for (int i = 0; i < BEAM_WIDTH; ++i) {
    candidates[i].maxRT.Initialize(nowRot, nowTip);
  }

  for (int depth = 0; depth < beamDepth; ++depth) {
    if (depth == 0) {
      Beam2_Internal(_t, x, y, nowRot, nowTip, a, b, 0, mCount, BEAM_WIDTH, candidates);
    }
    else {
      vector<MaxCandidate> candidates2(BEAM_WIDTH);
      for (int i = 0; i < BEAM_WIDTH; ++i) {
        candidates2[i].maxRT.Initialize(nowRot, nowTip);
      }

      for (int idx = 0; idx < BEAM_WIDTH; ++idx) {
        auto keepNowRot = nowRot;
        auto keepNowTip = nowTip;
        UpdateTurn(candidates[idx], nowRot, nowTip, x, y, _t, a, b, mCount);

        vector<MaxCandidate> tmpCand(BEAM_WIDTH);
        for (int i = 0; i < BEAM_WIDTH; ++i) {
          tmpCand[i].maxRT.Initialize(nowRot, nowTip);
        }

        Beam2_Internal(_t, x, y, nowRot, nowTip, a, b, 0, mCount, BEAM_WIDTH, tmpCand);

        nowRot = keepNowRot;
        nowTip = keepNowTip;
        RollBackTurn(candidates[idx], x, y, _t, a, b, mCount);

        for (int i = 0; i < BEAM_WIDTH; ++i) {
          if (tmpCand[i].finishTurn < candidates2[BEAM_WIDTH - 1].finishTurn) {
            candidates2[BEAM_WIDTH - 1] = tmpCand[i];
          }
          else if (tmpCand[i].finishTurn == candidates2[BEAM_WIDTH - 1].finishTurn && tmpCand[i].maxActionScore > candidates2[BEAM_WIDTH - 1].maxActionScore) {
            candidates2[BEAM_WIDTH - 1] = tmpCand[i];
          }

          for (int j = BEAM_WIDTH - 1; j >= 1; --j) {
            if (candidates2[j].finishTurn < candidates2[j - 1].finishTurn) {
              swap(candidates2[j], candidates2[j - 1]);
            }
            else if (candidates2[j].finishTurn == candidates2[j - 1].finishTurn && candidates2[j].maxActionScore > candidates2[j - 1].maxActionScore) {
              swap(candidates2[j], candidates2[j - 1]);
            }
            else {
              break;
            }
          }
        }
      }

      candidates = candidates2;
    }
  }

  return candidates[0];
}

// メソッド選択の共通処理
int SelectMethod(int v)
{
  if (v < 7) {
    int rand = Rand() % 100;
    if (rand < 20) return 42;
    else if (rand < 60) return 52;
    else return 62;
  }
  else {
    int rand = Rand() % 100;
    if (rand < 40) return 52;
    else if (rand < 60) return 62;
    else if (rand < 80) return 72;
    else return 82;
  }
}

// メソッドに基づいて木を作成する共通処理
void CreateTreeByMethod(int method)
{
  switch (method) {
    case 42: MakeTree1(); break;
    case 52: MakeTree2(); break;
    case 62: MakeTree6(); break;
    case 72: MakeTree4(); break;
    case 82: MakeTree5(); break;
    case 53: MakeTree22(); break;
    case 63: MakeTree32(); break;
    case 73: MakeTree42(); break;
  }
}

void Method100(double timeLimit, int case_num, ofstream& ofs)
{
  ResetTime();

  int loop[100] = {};
  int nn2 = n * n * 2;
  while (true) {
    double time = GetNowTime();
    if (time > timeLimit)break;

    // メソッド決定
    Method = SelectMethod(v);


    //if (v >= 8) {
    //  Method = 73;
    //}

    if (time > timeLimit * 0.5) {
      Method = best_Method;
    }


    loop[Method]++;

    // 木作成
    CreateTreeByMethod(Method);

    // 初期位置作成
    sx = Rand() % n;
    sy = Rand() % n;

    while (sx < 3 || n - 3 <= sx) {
      sx = Rand() % n;
    }
    while (sy < 3 || n - 3 <= sy) {
      sy = Rand() % n;
    }

    if (Rand() % 2 && Method == best_Method && best_ansCount[Method] < 999) {
      for (int i = 0; i < v; ++i) {
        pa[i] = best_pa[i];
        le[i] = best_le[i];
      }
    }

    // シミュレーション開始
    int x = sx;
    int y = sy;
    int _t = 0;

    vector<int> nowRot(v);
    vector<int> nowTip(v);

    vector<vector<int>> a;
    vector<vector<int>> b;
    int mCount = 0;
    CopyAB(a, b, mCount);
    PCount[Method].clear();
    PCount[Method].push_back(mCount * 2);

    int walkCount = 0;
    while (mCount < m && _t < best_ansCount[Method] && _t < best_total_ansCount + 2) {
      int xx = x;
      int yy = y;
      int tt = _t;
      int mCount2 = mCount;

      MaxCandidate candidates = Beam(_t, x, y, nowRot, nowTip, a, b, 0, mCount);
      //MaxCandidate candidates = Beam2(_t, x, y, nowRot, nowTip, a, b, 3, mCount);

      //FisherYates(order, 5);

      //int BEAM_WIDTH = 1;
      //vector<MaxCandidate> candidates(BEAM_WIDTH);
      //for (int i = 0; i < BEAM_WIDTH; ++i)
      //{
      //  candidates[i].maxRT.Initialize(nowRot, nowTip);
      //}

      //if (Method == 42) {
      //  DecideBest42(x, y, nowRot, nowTip, candidates, a, b);
      //}
      //else if (Method == 52) {
      //  DecideBest52(x, y, nowRot, nowTip, candidates, a, b);
      //}
      //else  if (Method == 62) {
      //  DecideBest62(x, y, nowRot, nowTip, candidates, a, b);
      //}
      //else  if (Method == 72) {
      //  DecideBest72(x, y, nowRot, nowTip, candidates, a, b);
      //}
      //else {
      //  if (mode != 0) {
      //    cout << "NG" << endl;
      //  }
      //}

      if (candidates.maxAB.KeepACount == 0 && candidates.maxAB.KeepBCount == 0) {
        walkCount++;
      }
      else {
        walkCount = 0;
      }
      if (walkCount == 10) {
        break;
      }

      bool isOk = UpdateTurn(candidates, nowRot, nowTip, x, y, _t, a, b, mCount);
      if (!isOk)break;
    }

    ansCount = _t;
    if (mCount == m && ansCount <= best_total_ansCount) {
      CopyToBest();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << best_ansCount[Method];
        cout << ", sx = " << sx << ", sy = " << sy;
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;

        if (ofs.is_open()) {
          ofs.close();
        }
        OpenOfs(case_num, ofs);
        Output(ofs);
      }
    }
    else if (mCount == m && ansCount <= best_ansCount[Method]) {
      CopyToBestPCount();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << best_ansCount[Method];
        cout << ", sx = " << sx << ", sy = " << sy;
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;
      }
    }
  }

  if (mode == 2) {
    cout << "Method42 loop = " << loop[42] << " " << endl;
    cout << "Method52 loop = " << loop[52] << " " << endl;
    cout << "Method53 loop = " << loop[53] << " " << endl;
    cout << "Method62 loop = " << loop[62] << " " << endl;
    cout << "Method63 loop = " << loop[63] << " " << endl;
    cout << "Method72 loop = " << loop[72] << " " << endl;
    cout << "Method73 loop = " << loop[73] << " " << endl;
    cout << "Method82 loop = " << loop[82] << " " << endl;
    cout << doRowColumnCount << " " << doMarginCount << " / " << doOneSetCount << endl;
    doRowColumnCount = 0;
    doMarginCount = 0;
    doOneSetCount = 0;
  }

  //　ベストな解に対してビームサーチ
  ResetTime();
  int iter = 0;
  while (true) {
    iter++;

    CopyFromBest();


    if (mode != 0) {
      if (Method == 0) {
        cout << "Miss" << endl;
        return;
      }
    }

    // シミュレーション開始
    int x = sx;
    int y = sy;
    int _t = 0;

    vector<int> nowRot(v);
    vector<int> nowTip(v);

    vector<vector<int>> a;
    vector<vector<int>> b;
    int mCount = 0;
    CopyAB(a, b, mCount);
    PCount[Method].clear();
    PCount[Method].push_back(mCount * 2);

    while (mCount < m && _t < best_total_ansCount) {
      MaxCandidate candidates = Beam(_t, x, y, nowRot, nowTip, a, b, iter, mCount);
      bool isOk = UpdateTurn(candidates, nowRot, nowTip, x, y, _t, a, b, mCount);
      if (!isOk)break;

      double time = GetNowTime();
      if (time > 0.1)break;
    }

    double time = GetNowTime();
    if (time > 0.1)break;
    if (mode == 2) {
      cout << iter << ' ' << _t << ' ' << best_total_ansCount << endl;
    }
    ansCount = _t;
    if (mCount == m && ansCount < best_total_ansCount) {
      if (mode != 0) {
        cout << "OK" << iter << ' ' << best_total_ansCount << ' ' << ansCount << endl;
      }

      CopyToBest();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << best_ansCount[Method];
        cout << ", sx = " << sx << ", sy = " << sy;
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;
      }
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

  Method = 0;
  ansCount = 99999;

  CopyToBest();

  Method100(TL, case_num, ofs);

  CopyFromBest();

  // 解答を出力
  if (ofs.is_open()) {
    ofs.close();
  }
  OpenOfs(case_num, ofs);
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
 TODO
 ・取り順と置き順すごい大事そう

 ・腕の長さをランダムにしない

 ・枝刈り高速化
   ・縦と横に分ける
     ・盤外確定は枝刈り

 ・ビームサーチ化

 ・近傍探索

 ・アドホックな手法考える

落ちているたこ焼きの取り順
 ・孤立してるやつ
 ・壁に近いやつ
 ・腕は一直線に取った方がお得？

一旦置くのアリ？

リファクタ

 ・関数の共通化
   ・ひとまず共通化する

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode <= 100) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(i);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "N = " << setw(2) << n << ", ";
        cout << "M = " << setw(3) << m << ", ";
        cout << "V = " << setw(2) << v << ", ";
        cout << "Method = " << Method << ", ";
        cout << "ArmLengthMethod = " << ArmLengthMethod << ", ";
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;
      }
    }
  }

  return 0;
}
