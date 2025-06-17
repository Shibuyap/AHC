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
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
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

int real_V;
int real_pa[MAX_V];
int real_le[MAX_V];
int real_ansCount[100];
int real_real_ansCount;
int real_dir[MAX_T];
int real_rot[MAX_T][MAX_V];
int real_tip[MAX_T][MAX_V];
int real_sx;
int real_sy;
int real_leafs1;
int real_Method;
int real_ArmLengthMethod;
vector<int> real_PCount[100];

int route[MAX_N * MAX_N * 3][2];

double heavyX;
double heavyY;


void CopyToReal()
{
  real_V = V;
  for (int i = 0; i < V; ++i)
  {
    real_pa[i] = pa[i];
    real_le[i] = le[i];
  }

  real_real_ansCount = ansCount;
  for (int i = 0; i < ansCount; ++i)
  {
    real_dir[i] = dir[i];
    for (int j = 0; j < V; ++j)
    {
      real_rot[i][j] = rot[i][j];
      real_tip[i][j] = tip[i][j];
    }
  }
  real_sx = sx;
  real_sy = sy;
  real_leafs1 = leafs1;
  real_Method = Method;
  real_ArmLengthMethod = ArmLengthMethod;

  real_PCount[Method] = PCount[Method];
  real_ansCount[Method] = ansCount;
}

void CopyToAns()
{
  V = real_V;
  for (int i = 0; i < V; ++i)
  {
    pa[i] = real_pa[i];
    le[i] = real_le[i];
  }
  ansCount = real_real_ansCount;
  for (int i = 0; i < ansCount; ++i)
  {
    dir[i] = real_dir[i];
    for (int j = 0; j < V; ++j)
    {
      rot[i][j] = real_rot[i][j];
      tip[i][j] = real_tip[i][j];
    }
  }
  sx = real_sx;
  sy = real_sy;
  leafs1 = real_leafs1;
  Method = real_Method;
  ArmLengthMethod = real_ArmLengthMethod;

  PCount[Method] = real_PCount[Method];
  //ansCount = real_ansCount[Method];
}

void CopyToRealPCount()
{
  real_ansCount[Method] = ansCount;
  real_PCount[Method] = PCount[Method];
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
  for (int i = 0; i < 100; ++i)
  {
    PCount[i].clear();
    real_PCount[i].clear();
    real_ansCount[i] = 999;
  }
}

// 入力受け取り
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> n >> m >> v;
    for (int i = 0; i < n; ++i)
    {
      string s;
      cin >> s;
      for (int j = 0; j < n; ++j)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    for (int i = 0; i < n; ++i)
    {
      string t;
      cin >> t;
      for (int j = 0; j < n; ++j)
      {
        init_b[i][j] = t[j] - '0';
      }
    }
  }
  else {
    // ファイル入力する
    ifs >> n >> m >> v;
    for (int i = 0; i < n; ++i)
    {
      string s;
      ifs >> s;
      for (int j = 0; j < n; ++j)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    for (int i = 0; i < n; ++i)
    {
      string t;
      ifs >> t;
      for (int j = 0; j < n; ++j)
      {
        init_b[i][j] = t[j] - '0';
      }
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
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      a[i][j] = init_a[i][j];
      b[i][j] = init_b[i][j];
    }
  }

  vector<int> nowRot(V);
  vector<int> nowTip(V);
  int x = sx;
  int y = sy;

  for (int t = 0; t < ansCount; ++t)
  {
    x = x + dx[dir[t]];
    y = y + dy[dir[t]];
    if (IsNG(x, y)) {
      cout << "Out of Range: turn = " << t << ", x = " << x << ", y = " << y << endl;
      return false;
    }

    for (int i = 1; i < V; ++i)
    {
      nowRot[i] = (nowRot[i] + rot[t][i]) % 4;
    }

    for (int i = 0; i < V; ++i)
    {
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
    for (int i = 1; i < V; ++i)
    {
      cout << pa[i] << ' ' << le[i] << endl;
    }
    cout << sx << ' ' << sy << endl;
    for (int i = 0; i < ansCount; ++i)
    {
      cout << dirChar[dir[i]];
      for (int j = 1; j < v; ++j)cout << rotChar[rot[i][j] + 1]; // rotは-1〜1のため帳尻合わせ
      for (int j = 0; j < V; ++j)cout << tipChar[tip[i][j]];
      cout << endl;
    }
  }
  else {
    ofs << V << endl;
    for (int i = 1; i < V; ++i)
    {
      ofs << pa[i] << ' ' << le[i] << endl;
    }
    ofs << sx << ' ' << sy << endl;
    for (int i = 0; i < ansCount; ++i)
    {
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
  heavyX = 0;
  heavyY = 0;
  a.resize(n, vector<int>(n));
  b.resize(n, vector<int>(n));
  for (int i = 0; i < n; ++i)
  {
    aRow[i] = 0;
    aColumn[i] = 0;
    bRow[i] = 0;
    bColumn[i] = 0;
  }
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
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

      heavyX += a[i][j] * i + b[i][j] * i;
      heavyY += a[i][j] * j + b[i][j] * j;
    }
  }

  heavyX /= 2.0 * (m - mCount);
  heavyY /= 2.0 * (m - mCount);
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
    for (int i = 0; i < m_V; ++i)
    {
      Rot[i] = 0;
      Tip[i] = 0;
      NowRot[i] = nowRot[i];
      NowTip[i] = nowTip[i];
    }
  }

  void Reflect(vector<int>& nowRot, vector<int>& nowTip, const int _t) const
  {
    for (int i = 0; i < V; ++i)
    {
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
    for (int i = 0; i < KeepACount; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        KeepA[i][j] = src.KeepA[i][j];
      }
    }
    KeepBCount = src.KeepBCount;
    for (int i = 0; i < KeepBCount; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        KeepB[i][j] = src.KeepB[i][j];
      }
    }
  }
};

// 配列の要素を復元して行・列のカウントを更新する共通処理
void RestoreArrayWithCount(vector<vector<int>>& arr, int row[], int column[],
  const int keep[][3], int count)
{
  for (int i = 0; i < count; ++i)
  {
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
  for (int i = 0; i < count; ++i)
  {
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
  for (int i = 0; i < count; ++i)
  {
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
void CalculatePositionBonus(double& actionScore, int nrx, int nry,
  const vector<vector<int>>& arr)
{
  for (int k = 0; k < 4; ++k)
  {
    for (int l = 1; l < 2; ++l)
    {
      int nkrx = nrx + dx[k] * l;
      int nkry = nry + dy[k] * l;
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
  const int i, const int j, const int nrx, const int nry)
{
  actionScore = 0;

  if (nowTip[i] == 0) {
    if (a[nrx][nry] == 1) {
      keepAB.AddA(nrx, nry, a[nrx][nry]);
      a[nrx][nry] = 0;
      aRow[nrx]--;
      aColumn[nry]--;
      tmpRT.Rot[i] = j;
      tmpRT.NowRot[i] = (nowRot[i] + j) % 4;
      tmpRT.Tip[i] = 1;
      tmpRT.NowTip[i] = 1;
      actionScore = ACTION_RATIO_A;
      actionScore += abs(heavyX - nrx);
      actionScore += abs(heavyY - nry);
      //actionScore += 10 * (20 - le[i]);

      CalculatePositionBonus(actionScore, nrx, nry, a);

      //for (int l = 1; l < 10; ++l)
      //{
      //  int nkrx = nrx + l;
      //  if (IsNG(nkrx, nry))break;
      //  if (a[nkrx][nry] == 0)break;
      //  nkrx = nrx - l;
      //  if (IsNG(nkrx, nry))break;
      //  if (a[nkrx][nry] == 0)break;
      //  actionScore -= POSITION_RATIO * 3453;
      //}
      //for (int l = 1; l < 10; ++l)
      //{
      //  int nkry = nry + l;
      //  if (IsNG(nrx, nkry))break;
      //  if (a[nrx][nkry] == 0)break;
      //  nkry = nry - l;
      //  if (IsNG(nrx, nkry))break;
      //  if (a[nrx][nkry] == 0)break;
      //  actionScore -= POSITION_RATIO * 3453;
      //}

      //if (aRow[nrx] == 0)actionScore += POSITION_RATIO * 100;
      //if (aColumn[nry] == 0)actionScore += POSITION_RATIO * 100;
      return true;
    }
  }
  else {
    if (b[nrx][nry] == 1) {
      keepAB.AddB(nrx, nry, b[nrx][nry]);
      b[nrx][nry] = 0;
      bRow[nrx]--;
      bColumn[nry]--;
      tmpRT.Rot[i] = j;
      tmpRT.NowRot[i] = (nowRot[i] + j) % 4;
      tmpRT.Tip[i] = 1;
      tmpRT.NowTip[i] = 0;
      actionScore = ACTION_RATIO_B;
      actionScore += abs(heavyX - nrx);
      actionScore += abs(heavyY - nry);
      //actionScore += 10 * (20 - le[i]);

      CalculatePositionBonus(actionScore, nrx, nry, b);

      //if (bRow[nrx] == 0)actionScore += POSITION_RATIO * 100;
      //if (bColumn[nry] == 0)actionScore += POSITION_RATIO * 100;
      return true;
    }
  }

  return false;
}

int MakeLength(int ra)
{
  int length = 1;
  if (ra < 30) {
    ArmLengthMethod = 0;
    length = Rand() % ((n - 1) * 1 / 3) + 1;
  }
  else if (ra < 60) {
    ArmLengthMethod = 1;
    length = Rand() % ((n - 1) * 1 / 2) + 1;
  }
  else if (ra < 80) {
    ArmLengthMethod = 2;
    length = Rand() % ((n - 1) * 2 / 3) + 1;
  }
  else if (ra < 100) {
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
void SetTreeLengths(int ra, int startIdx, int baseOffset, int range)
{
  for (int i = 1; i < V; ++i)
  {
    le[i] = MakeLength(ra);
  }

  if (Rand() % 2) {
    int st = Rand() % range + 1;
    for (int i = startIdx; i < V; ++i)le[i] = i - startIdx + st;
  }
}

void MakeTree1()
{
  int ra = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1 };
  SetTreeParents(parentRules, 2);
  SetTreeLengths(ra, 2, 2, 2);
}

void MakeTree2()
{
  int ra = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2 };
  SetTreeParents(parentRules, 3);
  SetTreeLengths(ra, 3, 3, 3);
}

void MakeTree6()
{
  int ra = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3 };
  SetTreeParents(parentRules, 4);
  SetTreeLengths(ra, 4, 4, 4);
}

void MakeTree4()
{
  int ra = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3, 4 };
  SetTreeParents(parentRules, 5);
  SetTreeLengths(ra, 5, 5, 5);
}

void MakeTree5()
{
  int ra = Rand() % 100;
  V = v;

  int parentRules[] = { 0, 1, 2, 3, 4, 5 };
  SetTreeParents(parentRules, 6);
  SetTreeLengths(ra, 6, 6, 5);
}

void MakeTree22()
{
  int ra = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 6) + 1;
  for (int i = 1; i < V; ++i)
  {
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
    le[i] = MakeLength(ra);
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
  int ra = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 7) + 1;
  for (int i = 1; i < V; ++i)
  {
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
    le[i] = MakeLength(ra);
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
  int ra = Rand() % 100;
  V = v;
  leafs1 = Rand() % (V - 8) + 1;
  for (int i = 1; i < V; ++i)
  {
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
    le[i] = MakeLength(ra);
  }

  if (Rand() % 2) {
    int st = Rand() % 5 + 1;
    for (int i = 6; i < 6 + leafs1; ++i)le[i] = i - 6 + st;
    st = Rand() % 5 + 1;
    for (int i = 6 + leafs1; i < V; ++i)le[i] = i - (6 + leafs1) + st;
  }
}

int CalcNeedLength(const int prenrx, const int prenry)
{
  int needLength = 999;
  if (prenrx < 0) {
    needLength = 0 - prenrx;
  }
  else if (prenrx < n) {
    needLength = 0;
  }
  else {
    needLength = prenrx - (n - 1);
  }
  if (prenry < 0) {
    needLength = min(needLength, 0 - prenry);
  }
  else if (prenrx < n) {
    needLength = 0;
  }
  else {
    needLength = min(needLength, prenry - (n - 1));
  }

  return needLength;
}

int CalcMarginCount(const int leafCount, const RotTip& maxRT)
{
  int marginCount = leafCount;
  for (int i = 0; i < V; ++i)
  {
    if (maxRT.Tip[i] == 1) {
      marginCount--;
    }
  }
  return marginCount;
}

// CallDecideBestの前方宣言
void CallDecideBest(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b);

int doMarginCount = 0;
int doOneSetCount = 0;
int doRowColumnCount = 0;
double DoOneSet(RotTip& tmpRT, KeepAB& keepAB,
  const int prenrx, const int prenry, const int needLength, const int prenRot, const int startLeaf,
  const vector<int>& nowRot, const vector<int>& nowTip,
  const MaxCandidate& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{

  doOneSetCount++;

  if (prenrx < 0 || n <= prenrx || (aRow[prenrx] == 0 && bRow[prenrx] == 0)) {
    if (prenry < 0 || n <= prenry || (aColumn[prenry] == 0 && bColumn[prenry] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  int margin = maxCand.maxMarginCount;
  double actionScore = 0;
  double actionScoreSum = 0;
  keepAB.Clear();

  for (int i = startLeaf; i < V; ++i)
  {
    if (le[i] < needLength)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j)
    {
      int nRot = (prenRot + nowRot[i] + j + 4) % 4;
      int nrx = prenrx + le[i] * dx[nRot];
      int nry = prenry + le[i] * dy[nRot];

      if (IsNG(nrx, nry))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, nrx, nry);
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
  const int prenrx1, const int prenry1, const int needLength1, const int prenRot1,
  const int prenrx2, const int prenry2, const int needLength2, const int prenRot2,
  const int startLeaf, const vector<int>& nowRot, const vector<int>& nowTip,
  const MaxCandidate& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{

  doOneSetCount++;

  if (prenrx1 < 0 || n <= prenrx1 || (aRow[prenrx1] == 0 && bRow[prenrx1] == 0)) {
    if (prenry1 < 0 || n <= prenry1 || (aColumn[prenry1] == 0 && bColumn[prenry1] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  if (prenrx2 < 0 || n <= prenrx2 || (aRow[prenrx2] == 0 && bRow[prenrx2] == 0)) {
    if (prenry2 < 0 || n <= prenry2 || (aColumn[prenry2] == 0 && bColumn[prenry2] == 0)) {
      doRowColumnCount++;
      return 0;
    }
  }

  int margin = maxCand.maxMarginCount;
  double actionScore = 0;
  double actionScoreSum = 0;
  keepAB.Clear();

  for (int i = startLeaf; i < startLeaf + leafs1; ++i)
  {
    if (le[i] < needLength1)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j)
    {
      int nRot = (prenRot1 + nowRot[i] + j + 4) % 4;
      int nrx = prenrx1 + le[i] * dx[nRot];
      int nry = prenry1 + le[i] * dy[nRot];

      if (IsNG(nrx, nry))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, nrx, nry);
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

  for (int i = startLeaf + leafs1; i < V; ++i)
  {
    if (margin < 0) {
      doMarginCount++;
      break;
    }

    if (le[i] < needLength2)continue;

    bool isCatch = false;
    for (int j = -1; j < 2; ++j)
    {
      int nRot = (prenRot2 + nowRot[i] + j + 4) % 4;
      int nrx = prenrx2 + le[i] * dx[nRot];
      int nry = prenry2 + le[i] * dy[nRot];

      if (IsNG(nrx, nry))continue;

      isCatch = CanCatch(nowRot, nowTip, a, b, tmpRT, keepAB, actionScore, i, j, nrx, nry);
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

void SetMaxCand(vector<MaxCandidate>& maxCand, const int ordDir, const double actionScoreSum, const RotTip& tmpRT, const KeepAB& keepAB, const int startLeaf)
{
  int beamWidth = maxCand.size();

  maxCand[beamWidth - 1].maxDir = ordDir;
  maxCand[beamWidth - 1].maxActionScore = actionScoreSum;
  maxCand[beamWidth - 1].maxRT = tmpRT;
  maxCand[beamWidth - 1].maxAB.Copy(keepAB);
  maxCand[beamWidth - 1].maxMarginCount = CalcMarginCount(V - startLeaf, maxCand[beamWidth - 1].maxRT);

  int idx = beamWidth - 1;
  while (idx > 0 && maxCand[idx].maxActionScore > maxCand[idx - 1].maxActionScore) {
    swap(maxCand[idx], maxCand[idx - 1]);
    idx--;
  }
}

void DecideBest42(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      tmpRT.Initialize(nowRot, nowTip);
      tmpRT.Rot[1] = ii;
      tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;

      int prenrx = nx + le[1] * dx[nRot1];
      int prenry = ny + le[1] * dy[nRot1];

      int needLength = CalcNeedLength(prenrx, prenry);

      double actionScoreSum = DoOneSet(tmpRT, keepAB, prenrx, prenry, needLength, nRot1,
        2, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

      if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
        SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 2);
      }
    }
  }
}

void DecideBest52(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[2]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        tmpRT.Initialize(nowRot, nowTip);
        tmpRT.Rot[1] = ii;
        tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
        tmpRT.Rot[2] = iii;
        tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;

        int prenrx = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2];
        int prenry = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2];

        int needLength = CalcNeedLength(prenrx, prenry);

        double actionScoreSum = DoOneSet(tmpRT, keepAB, prenrx, prenry, needLength, nRot2,
          3, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

        if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
          SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 3);
        }
      }
    }
  }
}

void DecideBest62(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[2]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[4]]) {
          int nRot3 = (nRot2 + nowRot[3] + iiii + 4) % 4;

          tmpRT.Initialize(nowRot, nowTip);
          tmpRT.Rot[1] = ii;
          tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
          tmpRT.Rot[2] = iii;
          tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
          tmpRT.Rot[3] = iiii;
          tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;

          int prenrx = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3];
          int prenry = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3];

          int needLength = CalcNeedLength(prenrx, prenry);

          double actionScoreSum = DoOneSet(tmpRT, keepAB, prenrx, prenry, needLength, nRot3,
            4, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

          if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
            SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 4);
          }
        }
      }
    }
  }
}

void DecideBest72(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[2]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[4]]) {
          int nRot3 = (nRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[ra[5]]) {
            int nRot4 = (nRot3 + nowRot[4] + iiiii + 4) % 4;

            tmpRT.Initialize(nowRot, nowTip);
            tmpRT.Rot[1] = ii;
            tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
            tmpRT.Rot[2] = iii;
            tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
            tmpRT.Rot[3] = iiii;
            tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
            tmpRT.Rot[4] = iiiii;
            tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;

            int prenrx = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3] + le[4] * dx[nRot4];
            int prenry = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3] + le[4] * dy[nRot4];

            int needLength = CalcNeedLength(prenrx, prenry);

            double actionScoreSum = DoOneSet(tmpRT, keepAB, prenrx, prenry, needLength, nRot4,
              5, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

            if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
              SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 5);
            }
          }
        }
      }
    }
  }
}

void DecideBest82(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[1]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[2]]) {
          int nRot3 = (nRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[ra[3]]) {
            int nRot4 = (nRot3 + nowRot[4] + iiiii + 4) % 4;

            for (const auto iiiiii : ordVec[ra[4]]) {
              int nRot5 = (nRot4 + nowRot[5] + iiiiii + 4) % 4;

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

              int prenrx = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3] + le[4] * dx[nRot4] + le[5] * dx[nRot5];
              int prenry = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3] + le[4] * dy[nRot4] + le[5] * dy[nRot5];

              int needLength = CalcNeedLength(prenrx, prenry);

              double actionScoreSum = DoOneSet(tmpRT, keepAB, prenrx, prenry, needLength, nRot5,
                6, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

              if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
                SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 6);
              }
            }
          }
        }
      }
    }
  }
}


void DecideBest53(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[2]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[4]]) {
          int nRot3 = (nRot1 + nowRot[3] + iiii + 4) % 4;

          tmpRT.Initialize(nowRot, nowTip);
          tmpRT.Rot[1] = ii;
          tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
          tmpRT.Rot[2] = iii;
          tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
          tmpRT.Rot[3] = iiii;
          tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;

          int prenrx1 = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2];
          int prenry1 = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2];
          int prenrx2 = nx + le[1] * dx[nRot1] + le[3] * dx[nRot3];
          int prenry2 = ny + le[1] * dy[nRot1] + le[3] * dy[nRot3];

          int needLength1 = CalcNeedLength(prenrx1, prenry1);
          int needLength2 = CalcNeedLength(prenrx2, prenry2);

          double actionScoreSum = DoOneSet2(tmpRT, keepAB,
            prenrx1, prenry1, needLength1, nRot2,
            prenrx2, prenry2, needLength2, nRot3,
            4, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

          if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
            SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 4);
          }
        }
      }
    }
  }
}

void DecideBest63(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[2]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[4]]) {
          int nRot3 = (nRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[ra[5]]) {
            int nRot4 = (nRot2 + nowRot[4] + iiiii + 4) % 4;

            tmpRT.Initialize(nowRot, nowTip);
            tmpRT.Rot[1] = ii;
            tmpRT.NowRot[1] = (nowRot[1] + ii) % 4;
            tmpRT.Rot[2] = iii;
            tmpRT.NowRot[2] = (nowRot[2] + iii) % 4;
            tmpRT.Rot[3] = iiii;
            tmpRT.NowRot[3] = (nowRot[3] + iiii) % 4;
            tmpRT.Rot[4] = iiiii;
            tmpRT.NowRot[4] = (nowRot[4] + iiiii) % 4;

            int prenrx1 = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3];
            int prenry1 = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3];
            int prenrx2 = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[4] * dx[nRot4];
            int prenry2 = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[4] * dy[nRot4];

            int needLength1 = CalcNeedLength(prenrx1, prenry1);
            int needLength2 = CalcNeedLength(prenrx2, prenry2);

            double actionScoreSum = DoOneSet2(tmpRT, keepAB,
              prenrx1, prenry1, needLength1, nRot3,
              prenrx2, prenry2, needLength2, nRot4,
              5, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

            if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
              SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 5);
            }
          }
        }
      }
    }
  }
}

void DecideBest73(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  int beamWidth = maxCand.size();

  RotTip tmpRT(v);
  KeepAB keepAB;

  for (int ord = 0; ord < 5; ++ord)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    for (int i = 0; i < 6; ++i)ra[i] = Rand() % 6;
    for (const auto ii : ordVec[ra[0]]) {
      int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;

      for (const auto iii : ordVec[ra[1]]) {
        int nRot2 = (nRot1 + nowRot[2] + iii + 4) % 4;

        for (const auto iiii : ordVec[ra[2]]) {
          int nRot3 = (nRot2 + nowRot[3] + iiii + 4) % 4;

          for (const auto iiiii : ordVec[ra[3]]) {
            int nRot4 = (nRot3 + nowRot[4] + iiiii + 4) % 4;

            for (const auto iiiiii : ordVec[ra[4]]) {
              int nRot5 = (nRot3 + nowRot[5] + iiiiii + 4) % 4;

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

              int prenrx1 = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3] + le[4] * dx[nRot4];
              int prenry1 = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3] + le[4] * dy[nRot4];
              int prenrx2 = nx + le[1] * dx[nRot1] + le[2] * dx[nRot2] + le[3] * dx[nRot3] + le[5] * dx[nRot5];
              int prenry2 = ny + le[1] * dy[nRot1] + le[2] * dy[nRot2] + le[3] * dy[nRot3] + le[5] * dy[nRot5];

              int needLength1 = CalcNeedLength(prenrx1, prenry1);
              int needLength2 = CalcNeedLength(prenrx2, prenry2);

              double actionScoreSum = DoOneSet2(tmpRT, keepAB,
                prenrx1, prenry1, needLength1, nRot4,
                prenrx2, prenry2, needLength2, nRot5,
                6, nowRot, nowTip, maxCand[beamWidth - 1], a, b);

              if (actionScoreSum > maxCand[beamWidth - 1].maxActionScore) {
                SetMaxCand(maxCand, order[ord], actionScoreSum, tmpRT, keepAB, 6);
              }
            }
          }
        }
      }
    }
  }
}

bool UpdateTurn(const MaxCandidate& maxCand, vector<int>& nowRot, vector<int>& nowTip, int& x, int& y, int& _t,
  vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  bool isOk = true;

  maxCand.maxRT.Reflect(nowRot, nowTip, _t);

  dir[_t] = maxCand.maxDir;

  ReflectFromMaxAB(maxCand.maxAB, a, b, mCount);
  PCount[Method][_t] += maxCand.maxAB.KeepACount + maxCand.maxAB.KeepBCount;

  if (Method == real_Method) {
    if (real_PCount[Method].size() > 0 && _t >= real_PCount[Method].size()) {
      isOk = false;
    }
    else if (PCount[Method][_t] < real_PCount[Method][_t] - 5) {
      isOk = false;
    }
  }

  x += dx[maxCand.maxDir];
  y += dy[maxCand.maxDir];
  _t++;

  PCount[Method].push_back(PCount[Method][_t - 1]);

  return isOk;
}

void RollBackTurn(const MaxCandidate& maxCand, int& x, int& y, int& _t,
  vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  RollBackFromMaxAB(maxCand.maxAB, a, b, mCount);
  PCount[Method].pop_back();

  x -= dx[maxCand.maxDir];
  y -= dy[maxCand.maxDir];
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

  vector<MaxCandidate> maxCand(BEAM_WIDTH);
  for (int i = 0; i < BEAM_WIDTH; ++i)
  {
    maxCand[i].maxRT.Initialize(nowRot, nowTip);
  }

  CallDecideBest(x, y, nowRot, nowTip, maxCand, a, b);

  if (mCount + maxCand[0].maxAB.KeepBCount == m) {
    maxCand[0].finishTurn = 0;
    return maxCand[0];
  }

  if (beamDepth == 0) {
    return maxCand[0];
  }

  auto keepNowRot2 = nowRot;
  auto keepNowTip2 = nowTip;


  int minFinishTurn = 999;
  double maxActionScore = -1;
  int bestIndex = -1;
  for (int aespa = 0; aespa < BEAM_WIDTH; ++aespa)
  {
    auto keepNowRot = nowRot;
    auto keepNowTip = nowTip;
    UpdateTurn(maxCand[aespa], nowRot, nowTip, x, y, _t, a, b, mCount);

    MaxCandidate tmpMaxCand = Beam(_t, x, y, nowRot, nowTip, a, b, beamDepth - 1, mCount);

    nowRot = keepNowRot;
    nowTip = keepNowTip;
    RollBackTurn(maxCand[aespa], x, y, _t, a, b, mCount);

    maxCand[aespa].finishTurn = min(tmpMaxCand.finishTurn + 1, 999);
    maxCand[aespa].maxActionScore += tmpMaxCand.maxActionScore;

    if (maxCand[aespa].finishTurn < minFinishTurn) {
      minFinishTurn = maxCand[aespa].finishTurn;
      maxActionScore = maxCand[aespa].maxActionScore;
      bestIndex = aespa;
    }
    else if (maxCand[aespa].finishTurn == minFinishTurn && maxCand[aespa].maxActionScore > maxActionScore) {
      minFinishTurn = maxCand[aespa].finishTurn;
      maxActionScore = maxCand[aespa].maxActionScore;
      bestIndex = aespa;
    }
  }

  return maxCand[bestIndex];
}

// DecideBest関数を呼び出す共通処理
void CallDecideBest(const int x, const int y, const vector<int>& nowRot, const vector<int>& nowTip,
  vector<MaxCandidate>& maxCand, vector<vector<int>>& a, vector<vector<int>>& b)
{
  if (Method == 42) {
    DecideBest42(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else if (Method == 52) {
    DecideBest52(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 62) {
    DecideBest62(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 72) {
    DecideBest72(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 82) {
    DecideBest82(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 53) {
    DecideBest53(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 63) {
    DecideBest63(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else  if (Method == 73) {
    DecideBest73(x, y, nowRot, nowTip, maxCand, a, b);
  }
  else {
    if (mode != 0) {
      cout << "NG" << endl;
    }
  }
}

void Beam2_Internal(int& _t, int& x, int& y, vector<int>& nowRot, vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b, const int beamDepth, int& mCount, const int beamWidth, vector<MaxCandidate>& maxCand)
{
  FisherYates(order, 5);
  CallDecideBest(x, y, nowRot, nowTip, maxCand, a, b);
}

MaxCandidate Beam2(int& _t, int& x, int& y, vector<int>& nowRot, vector<int>& nowTip,
  vector<vector<int>>& a, vector<vector<int>>& b, const int beamDepth, int& mCount)
{
  int BEAM_WIDTH = 30;

  vector<MaxCandidate> maxCand(BEAM_WIDTH);
  for (int i = 0; i < BEAM_WIDTH; ++i)
  {
    maxCand[i].maxRT.Initialize(nowRot, nowTip);
  }

  for (int winter = 0; winter < beamDepth; ++winter)
  {
    if (winter == 0) {
      Beam2_Internal(_t, x, y, nowRot, nowTip, a, b, 0, mCount, BEAM_WIDTH, maxCand);
    }
    else {
      vector<MaxCandidate> maxCand2(BEAM_WIDTH);
      for (int i = 0; i < BEAM_WIDTH; ++i)
      {
        maxCand2[i].maxRT.Initialize(nowRot, nowTip);
      }

      for (int aespa = 0; aespa < BEAM_WIDTH; ++aespa)
      {
        auto keepNowRot = nowRot;
        auto keepNowTip = nowTip;
        UpdateTurn(maxCand[aespa], nowRot, nowTip, x, y, _t, a, b, mCount);

        vector<MaxCandidate> tmpCand(BEAM_WIDTH);
        for (int i = 0; i < BEAM_WIDTH; ++i)
        {
          tmpCand[i].maxRT.Initialize(nowRot, nowTip);
        }

        Beam2_Internal(_t, x, y, nowRot, nowTip, a, b, 0, mCount, BEAM_WIDTH, tmpCand);

        nowRot = keepNowRot;
        nowTip = keepNowTip;
        RollBackTurn(maxCand[aespa], x, y, _t, a, b, mCount);

        for (int i = 0; i < BEAM_WIDTH; ++i)
        {
          if (tmpCand[i].finishTurn < maxCand2[BEAM_WIDTH - 1].finishTurn) {
            maxCand2[BEAM_WIDTH - 1] = tmpCand[i];
          }
          else if (tmpCand[i].finishTurn == maxCand2[BEAM_WIDTH - 1].finishTurn && tmpCand[i].maxActionScore > maxCand2[BEAM_WIDTH - 1].maxActionScore) {
            maxCand2[BEAM_WIDTH - 1] = tmpCand[i];
          }

          for (int j = BEAM_WIDTH - 1; j >= 1; --j)
          {
            if (maxCand2[j].finishTurn < maxCand2[j - 1].finishTurn) {
              swap(maxCand2[j], maxCand2[j - 1]);
            }
            else if (maxCand2[j].finishTurn == maxCand2[j - 1].finishTurn && maxCand2[j].maxActionScore > maxCand2[j - 1].maxActionScore) {
              swap(maxCand2[j], maxCand2[j - 1]);
            }
            else {
              break;
            }
          }
        }
      }

      maxCand = maxCand2;
    }
  }

  return maxCand[0];
}

// メソッド選択の共通処理
int SelectMethod(int v)
{
  if (v < 7) {
    int ra = Rand() % 100;
    if (ra < 20) return 42;
    else if (ra < 60) return 52;
    else return 62;
  }
  else {
    int ra = Rand() % 100;
    if (ra < 40) return 52;
    else if (ra < 60) return 62;
    else if (ra < 80) return 72;
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

void Method100(double timeLimit, int probNum, ofstream& ofs)
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
      Method = real_Method;
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

    if (Rand() % 2 && Method == real_Method && real_ansCount[Method] < 999) {
      for (int i = 0; i < v; ++i)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
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
    while (mCount < m && _t < real_ansCount[Method] && _t < real_real_ansCount + 2) {
      int xx = x;
      int yy = y;
      int tt = _t;
      int mCount2 = mCount;

      MaxCandidate maxCand = Beam(_t, x, y, nowRot, nowTip, a, b, 0, mCount);
      //MaxCandidate maxCand = Beam2(_t, x, y, nowRot, nowTip, a, b, 3, mCount);

      //FisherYates(order, 5);

      //int BEAM_WIDTH = 1;
      //vector<MaxCandidate> maxCand(BEAM_WIDTH);
      //for (int i = 0; i < BEAM_WIDTH; ++i)
      //{
      //  maxCand[i].maxRT.Initialize(nowRot, nowTip);
      //}

      //if (Method == 42) {
      //  DecideBest42(x, y, nowRot, nowTip, maxCand, a, b);
      //}
      //else if (Method == 52) {
      //  DecideBest52(x, y, nowRot, nowTip, maxCand, a, b);
      //}
      //else  if (Method == 62) {
      //  DecideBest62(x, y, nowRot, nowTip, maxCand, a, b);
      //}
      //else  if (Method == 72) {
      //  DecideBest72(x, y, nowRot, nowTip, maxCand, a, b);
      //}
      //else {
      //  if (mode != 0) {
      //    cout << "NG" << endl;
      //  }
      //}

      if (maxCand.maxAB.KeepACount == 0 && maxCand.maxAB.KeepBCount == 0) {
        walkCount++;
      }
      else {
        walkCount = 0;
      }
      if (walkCount == 10) {
        break;
      }

      bool isOk = UpdateTurn(maxCand, nowRot, nowTip, x, y, _t, a, b, mCount);
      if (!isOk)break;
    }

    ansCount = _t;
    if (mCount == m && ansCount <= real_real_ansCount) {
      CopyToReal();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << real_ansCount[Method];
        cout << ", sx = " << sx << ", sy = " << sy;
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;

        if (ofs.is_open()) {
          ofs.close();
        }
        OpenOfs(probNum, ofs);
        Output(ofs);
      }
    }
    else if (mCount == m && ansCount <= real_ansCount[Method]) {
      CopyToRealPCount();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << real_ansCount[Method];
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
  int karina = 0;
  while (true) {
    karina++;

    CopyToAns();


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

    while (mCount < m && _t < real_real_ansCount) {
      MaxCandidate maxCand = Beam(_t, x, y, nowRot, nowTip, a, b, karina, mCount);
      bool isOk = UpdateTurn(maxCand, nowRot, nowTip, x, y, _t, a, b, mCount);
      if (!isOk)break;

      double time = GetNowTime();
      if (time > 0.1)break;
    }

    double time = GetNowTime();
    if (time > 0.1)break;
    if (mode == 2) {
      cout << karina << ' ' << _t << ' ' << real_real_ansCount << endl;
    }
    ansCount = _t;
    if (mCount == m && ansCount < real_real_ansCount) {
      if (mode != 0) {
        cout << "OK" << karina << ' ' << real_real_ansCount << ' ' << ansCount << endl;
      }

      CopyToReal();
      if (mode == 2) {
        cout << "Method" << Method << ", " << "loop = " << loop[Method];
        cout << ", score = " << real_ansCount[Method];
        cout << ", sx = " << sx << ", sy = " << sy;
        for (int i = 1; i < V; ++i)cout << ", " << le[i];
        cout << endl;
      }
    }
  }
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

  Method = 0;
  ansCount = 99999;

  CopyToReal();

  Method100(TL, probNum, ofs);

  CopyToAns();

  // 解答を出力
  if (ofs.is_open()) {
    ofs.close();
  }
  OpenOfs(probNum, ofs);
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
    for (int i = 0; i < 10; ++i)
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
