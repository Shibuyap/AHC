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
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

namespace /* 乱数ライブラリ */
{
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
clock_t startTime, endTime;

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
  rep(i, V)
  {
    real_pa[i] = pa[i];
    real_le[i] = le[i];
  }

  real_real_ansCount = ansCount;
  rep(i, ansCount)
  {
    real_dir[i] = dir[i];
    rep(j, V)
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
  rep(i, V)
  {
    pa[i] = real_pa[i];
    le[i] = real_le[i];
  }
  ansCount = real_real_ansCount;
  rep(i, ansCount)
  {
    dir[i] = real_dir[i];
    rep(j, V)
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

void CopyToRealPCount() {
  real_ansCount[Method] = ansCount;
  real_PCount[Method] = PCount[Method];
}

void ResetTime()
{
  startTime = clock();
}

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
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
  rep(i, 100)
  {
    PCount[i].clear();
    real_PCount[i].clear();
    real_ansCount[i] = 999;
  }
}

// 入力受け取り
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> n >> m >> v;
    rep(i, n)
    {
      string s;
      cin >> s;
      rep(j, n)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    rep(i, n)
    {
      string t;
      cin >> t;
      rep(j, n)
      {
        init_b[i][j] = t[j] - '0';
      }
    }
  }
  else {
    // ファイル入力する
    ifs >> n >> m >> v;
    rep(i, n)
    {
      string s;
      ifs >> s;
      rep(j, n)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    rep(i, n)
    {
      string t;
      ifs >> t;
      rep(j, n)
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
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
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
  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = init_a[i][j];
      b[i][j] = init_b[i][j];
    }
  }

  vector<int> nowRot(V);
  vector<int> nowTip(V);
  int x = sx;
  int y = sy;

  rep(t, ansCount)
  {
    x = x + dx[dir[t]];
    y = y + dy[dir[t]];
    if (IsNG(x, y)) {
      cout << "Out of Range: turn = " << t << ", x = " << x << ", y = " << y << endl;
      return false;
    }

    srep(i, 1, V)
    {
      nowRot[i] = (nowRot[i] + rot[t][i]) % 4;
    }

    rep(i, V)
    {
      if (tip[t][i] == 1) {
        // ここでアーム位置を計算
      }
    }
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    //if (ansCount < 100) {
    //  assert(false);
    //}
    cout << V << endl;
    srep(i, 1, V)
    {
      cout << pa[i] << ' ' << le[i] << endl;
    }
    cout << sx << ' ' << sy << endl;
    rep(i, ansCount)
    {
      cout << dirChar[dir[i]];
      srep(j, 1, v)cout << rotChar[rot[i][j] + 1]; // rotは-1～1のため帳尻合わせ
      rep(j, V)cout << tipChar[tip[i][j]];
      cout << endl;
    }
  }
  else {
    ofs << V << endl;
    srep(i, 1, V)
    {
      ofs << pa[i] << ' ' << le[i] << endl;
    }
    ofs << sx << ' ' << sy << endl;
    rep(i, ansCount)
    {
      ofs << dirChar[dir[i]];
      srep(j, 1, v)ofs << rotChar[rot[i][j] + 1]; // rotは-1～1のため帳尻合わせ
      rep(j, V)ofs << tipChar[tip[i][j]];
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
  rep(i, n)
  {
    aRow[i] = 0;
    aColumn[i] = 0;
    bRow[i] = 0;
    bColumn[i] = 0;
  }
  rep(i, n)
  {
    rep(j, n)
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
    rep(i, m_V)
    {
      Rot[i] = 0;
      Tip[i] = 0;
      NowRot[i] = nowRot[i];
      NowTip[i] = nowTip[i];
    }
  }

  void Reflect(vector<int>& nowRot, vector<int>& nowTip, const int _t) const
  {
    rep(i, V)
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
    rep(i, KeepACount)
    {
      rep(j, 3)
      {
        KeepA[i][j] = src.KeepA[i][j];
      }
    }
    KeepBCount = src.KeepBCount;
    rep(i, KeepBCount)
    {
      rep(j, 3)
      {
        KeepB[i][j] = src.KeepB[i][j];
      }
    }
  }
};

void RollBackFromKeepAB(const KeepAB& keepAB, vector<vector<int>>& a, vector<vector<int>>& b)
{
  rep(i, keepAB.KeepACount)
  {
    a[keepAB.KeepA[i][0]][keepAB.KeepA[i][1]] = keepAB.KeepA[i][2];
    aRow[keepAB.KeepA[i][0]]++;
    aColumn[keepAB.KeepA[i][1]]++;
  }
  rep(i, keepAB.KeepBCount)
  {
    b[keepAB.KeepB[i][0]][keepAB.KeepB[i][1]] = keepAB.KeepB[i][2];
    bRow[keepAB.KeepB[i][0]]++;
    bColumn[keepAB.KeepB[i][1]]++;
  }
}

void ReflectFromMaxAB(const KeepAB& maxAB, vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  mCount += maxAB.KeepBCount;
  rep(i, maxAB.KeepACount)
  {
    if (a[maxAB.KeepA[i][0]][maxAB.KeepA[i][1]] == 0) {
      assert(false);
    }
    a[maxAB.KeepA[i][0]][maxAB.KeepA[i][1]] = 0;
    aRow[maxAB.KeepA[i][0]]--;
    aColumn[maxAB.KeepA[i][1]]--;
  }
  rep(i, maxAB.KeepBCount)
  {
    if (b[maxAB.KeepB[i][0]][maxAB.KeepB[i][1]] == 0) {
      assert(false);
    }
    b[maxAB.KeepB[i][0]][maxAB.KeepB[i][1]] = 0;
    bRow[maxAB.KeepB[i][0]]--;
    bColumn[maxAB.KeepB[i][1]]--;
  }
}

void RollBackFromMaxAB(const KeepAB& maxAB, vector<vector<int>>& a, vector<vector<int>>& b, int& mCount)
{
  mCount -= maxAB.KeepBCount;
  rep(i, maxAB.KeepACount)
  {
    if (a[maxAB.KeepA[i][0]][maxAB.KeepA[i][1]] == 1) {
      assert(false);
    }
    a[maxAB.KeepA[i][0]][maxAB.KeepA[i][1]] = 1;
    aRow[maxAB.KeepA[i][0]]++;
    aColumn[maxAB.KeepA[i][1]]++;
  }
  rep(i, maxAB.KeepBCount)
  {
    if (b[maxAB.KeepB[i][0]][maxAB.KeepB[i][1]] == 1) {
      assert(false);
    }
    b[maxAB.KeepB[i][0]][maxAB.KeepB[i][1]] = 1;
    bRow[maxAB.KeepB[i][0]]++;
    bColumn[maxAB.KeepB[i][1]]++;
  }
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

      rep(k, 4)
      {
        srep(l, 1, 2)
        {
          int nkrx = nrx + dx[k] * l;
          int nkry = nry + dy[k] * l;
          if (IsNG(nkrx, nkry)) {
            actionScore += POSITION_RATIO;
          }
          else if (a[nkrx][nkry] == 0) {
            actionScore += POSITION_RATIO;
          }
          else {
            break;
          }
        }
      }

      //srep(l, 1, 10)
      //{
      //  int nkrx = nrx + l;
      //  if (IsNG(nkrx, nry))break;
      //  if (a[nkrx][nry] == 0)break;
      //  nkrx = nrx - l;
      //  if (IsNG(nkrx, nry))break;
      //  if (a[nkrx][nry] == 0)break;
      //  actionScore -= POSITION_RATIO * 3453;
      //}
      //srep(l, 1, 10)
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

      rep(k, 4)
      {
        srep(l, 1, 2)
        {
          int nkrx = nrx + dx[k] * l;
          int nkry = nry + dy[k] * l;
          if (IsNG(nkrx, nkry)) {
            actionScore += POSITION_RATIO;
          }
          else if (b[nkrx][nkry] == 0) {
            actionScore += POSITION_RATIO;
          }
          else {
            break;
          }
        }
      }

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
    length = randxor() % ((n - 1) * 1 / 3) + 1;
  }
  else if (ra < 60) {
    ArmLengthMethod = 1;
    length = randxor() % ((n - 1) * 1 / 2) + 1;
  }
  else if (ra < 80) {
    ArmLengthMethod = 2;
    length = randxor() % ((n - 1) * 2 / 3) + 1;
  }
  else if (ra < 100) {
    ArmLengthMethod = 3;
    length = randxor() % ((n - 1) * 3 / 4) + 1;
  }
  else {
    ArmLengthMethod = 4;
    length = randxor() % (n - 1) + 1;
  }
  return length;
}

void MakeTree1()
{
  int ra = randxor() % 100;
  V = v;
  srep(i, 1, V)
  {
    if (i == 1) {
      pa[i] = 0;
    }
    else {
      pa[i] = 1;
    }
    le[i] = MakeLength(ra);
  }

  if (randxor() % 2) {
    int st = randxor() % 2 + 1;
    srep(i, 2, V)le[i] = i - 2 + st;
  }
}

void MakeTree2()
{
  int ra = randxor() % 100;
  V = v;
  srep(i, 1, V)
  {
    if (i == 1) {
      pa[i] = 0;
    }
    else if (i == 2) {
      pa[i] = 1;
    }
    else {
      pa[i] = 2;
    }
    le[i] = MakeLength(ra);
  }

  if (randxor() % 2) {
    int st = randxor() % 3 + 1;
    srep(i, 3, V)le[i] = i - 3 + st;
    //if (le[V - 1] < n - 1 && V >= 5 && randxor() % 2 == 0) {
    //  int st2 = randxor() % (V - 4) + 4;
    //  srep(i, st2, V)le[i]++;
    //}
  }
}

void MakeTree6()
{
  int ra = randxor() % 100;
  V = v;
  srep(i, 1, V)
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
    else {
      pa[i] = 3;
    }
    le[i] = MakeLength(ra);
  }

  if (randxor() % 2) {
    int st = randxor() % 4 + 1;
    srep(i, 4, V)le[i] = i - 4 + st;
    //if (le[V - 1] < n - 1 && V >= 6 && randxor() % 2 == 0) {
    //  int st2 = randxor() % (V - 5) + 5;
    //  srep(i, st2, V)le[i]++;
    //}
  }
}

void MakeTree4()
{
  int ra = randxor() % 100;
  V = v;
  srep(i, 1, V)
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
    else {
      pa[i] = 4;
    }
    le[i] = MakeLength(ra);
  }

  if (randxor() % 2) {
    int st = randxor() % 5 + 1;
    srep(i, 5, V)le[i] = i - 5 + st;
    //if (le[V - 1] < n - 1 && V >= 7 && randxor() % 2 == 0) {
    //  int st2 = randxor() % (V - 6) + 6;
    //  srep(i, st2, V)le[i]++;
    //}
  }
}

void MakeTree5()
{
  int ra = randxor() % 100;
  V = v;
  srep(i, 1, V)
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
      pa[i] = 4;
    }
    else {
      pa[i] = 5;
    }
    le[i] = MakeLength(ra);
  }

  if (randxor() % 2) {
    int st = randxor() % 5 + 1;
    srep(i, 6, V)le[i] = i - 6 + st;
    //if (le[V - 1] < n - 1 && V >= 8 && randxor() % 2 == 0) {
    //  int st2 = randxor() % (V - 7) + 7;
    //  srep(i, st2, V)le[i]++;
    //}
  }
}

void MakeTree22()
{
  int ra = randxor() % 100;
  V = v;
  leafs1 = randxor() % (V - 6) + 1;
  srep(i, 1, V)
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

  if (randxor() % 2) {
    int st = randxor() % 4 + 1;
    srep(i, 4, 4 + leafs1)le[i] = i - 4 + st;
    st = randxor() % 4 + 1;
    srep(i, 4 + leafs1, V)le[i] = i - (4 + leafs1) + st;
  }
}

void MakeTree32()
{
  int ra = randxor() % 100;
  V = v;
  leafs1 = randxor() % (V - 7) + 1;
  srep(i, 1, V)
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

  if (randxor() % 2) {
    int st = randxor() % 5 + 1;
    srep(i, 5, 5 + leafs1)le[i] = i - 5 + st;
    st = randxor() % 5 + 1;
    srep(i, 5 + leafs1, V)le[i] = i - (5 + leafs1) + st;
  }
}

void MakeTree42()
{
  int ra = randxor() % 100;
  V = v;
  leafs1 = randxor() % (V - 8) + 1;
  srep(i, 1, V)
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

  if (randxor() % 2) {
    int st = randxor() % 5 + 1;
    srep(i, 6, 6 + leafs1)le[i] = i - 6 + st;
    st = randxor() % 5 + 1;
    srep(i, 6 + leafs1, V)le[i] = i - (6 + leafs1) + st;
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
  rep(i, V)
  {
    if (maxRT.Tip[i] == 1) {
      marginCount--;
    }
  }
  return marginCount;
}

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

  srep(i, startLeaf, V)
  {
    if (le[i] < needLength)continue;

    bool isCatch = false;
    srep(j, -1, 2)
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

  srep(i, startLeaf, startLeaf + leafs1)
  {
    if (le[i] < needLength1)continue;

    bool isCatch = false;
    srep(j, -1, 2)
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

  srep(i, startLeaf + leafs1, V)
  {
    if (margin < 0) {
      doMarginCount++;
      break;
    }

    if (le[i] < needLength2)continue;

    bool isCatch = false;
    srep(j, -1, 2)
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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

  rep(ord, 5)
  {
    // dir
    int nx = x + dx[order[ord]];
    int ny = y + dy[order[ord]];

    if (IsNG(nx, ny))continue;

    int ra[6];
    rep(i, 6)ra[i] = randxor() % 6;
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
  rep(i, BEAM_WIDTH)
  {
    maxCand[i].maxRT.Initialize(nowRot, nowTip);
  }

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
  rep(aespa, BEAM_WIDTH)
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

void Method100(double timeLimit, int probNum, ofstream& ofs)
{
  ResetTime();

  int loop[100] = {};
  int nn2 = n * n * 2;
  while (true) {
    double time = GetNowTime();
    if (time > timeLimit)break;

    // メソッド決定
    if (v < 7) {
      int ra = randxor() % 100;
      if (ra < 20) {
        Method = 42;
      }
      else if (ra < 60) {
        Method = 52;
      }
      else {
        Method = 62;
      }
    }
    else {
      int ra = randxor() % 100;
      if (ra < 40) {
        Method = 52;
      }
      else if (ra < 60) {
        Method = 62;
      }
      else if (ra < 80) {
        Method = 72;
      }
      else {
        Method = 82;
      }
    }


    //if (v >= 8) {
    //  Method = 73;
    //}

    if (time > timeLimit * 0.5) {
      Method = real_Method;
    }


    loop[Method]++;

    // 木作成
    switch (Method) {
    case 42:
      MakeTree1();
      break;
    case 52:
      MakeTree2();
      break;
    case 62:
      MakeTree6();
      break;
    case 72:
      MakeTree4();
      break;
    case 82:
      MakeTree5();
      break;
    case 53:
      MakeTree22();
      break;
    case 63:
      MakeTree32();
      break;
    case 73:
      MakeTree42();
      break;
    }

    // 初期位置作成
    sx = randxor() % n;
    sy = randxor() % n;

    while (sx < 3 || n - 3 <= sx) {
      sx = randxor() % n;
    }
    while (sy < 3 || n - 3 <= sy) {
      sy = randxor() % n;
    }

    if (randxor() % 2 && Method == real_Method && real_ansCount[Method] < 999) {
      rep(i, v) {
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

      //FisherYates(order, 5);

      //int BEAM_WIDTH = 1;
      //vector<MaxCandidate> maxCand(BEAM_WIDTH);
      //rep(i, BEAM_WIDTH)
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
        srep(i, 1, V)cout << ", " << le[i];
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
        srep(i, 1, V)cout << ", " << le[i];
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
  while (true)
  {
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
        srep(i, 1, V)cout << ", " << le[i];
        cout << endl;
      }
    }
  }
}

ll Solve(int probNum)
{
  startTime = clock();

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
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode <= 100) {
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
        cout << "N = " << setw(2) << n << ", ";
        cout << "M = " << setw(3) << m << ", ";
        cout << "V = " << setw(2) << v << ", ";
        cout << "Method = " << Method << ", ";
        cout << "ArmLengthMethod = " << ArmLengthMethod << ", ";
        srep(i, 1, V)cout << ", " << le[i];
        cout << endl;
      }
    }
  }

  return 0;
}
