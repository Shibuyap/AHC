#include <algorithm>
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
#include <array>
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef pair<int, int> P;
#define MAX_N 200005
#define INF 1001001001

const array<int, 4> dx = { -1, 0, 1, 0 };
const array<int, 4> dy = { 0, -1, 0, 1 };

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
}  // namespace

class GameState
{
public:
  int maxScore;
  int ope1, ope2;
  vector<vector<int>> ans1, ans2;
  array<array<int, 100>, 100> a;
  vector<int> x, y;
  vector<int> R, D;
  int viewOrder = 0;
  vector<vector<int>> moves;
  vector<int> moveCnt;
  array<array<int, 100>, 100> cellUse;
  vector<vector<int>> udlr;
  vector<int> parent;
  vector<int> unionSize;
  set<P> vp;
};

namespace /* 変数 */
{
  // 入力用変数
  int n, K;
  array<string, 100> s;

  // 解答用変数
  int maxScore;
  int ope1, ope2;
  vector<vector<int>> ans1, ans2;
  array<array<int, 100>, 100> a;
  vector<int> x, y;
  vector<int> R, D;
  int viewOrder = 0;
  vector<vector<int>> moves;
  vector<int> moveCnt;
  array<array<int, 100>, 100> cellUse;
  vector<vector<int>> udlr;
  vector<int> parent;
  vector<int> unionSize;
  set<P> vp;

  GameState real_GameState;

  int seed_maxScore;
  int seed_GameState.ope1, seed_ope2;
  vector<vector<int>> seed_GameState.ans1, seed_ans2;
  array<array<int, 100>, 100> seed_a;
  vector<int> seed_GameState.x, seed_y;
  vector<int> seed_GameState.R, seed_D;
  int seed_viewOrder = 0;
  vector<vector<int>> seed_moves;
  vector<int> seed_moveCnt;
  array<array<int, 100>, 100> seed_cellUse;
  vector<vector<int>> seed_udlr;
  vector<int> seed_parent;
  vector<int> seed_unionSize;
  set<P> seed_vp;

  GameState seed_GameState;
  int outer_maxScore;
  int outer_GameState.ope1, outer_ope2;
  vector<vector<int>> outer_GameState.ans1, outer_ans2;
  array<array<int, 100>, 100> outer_a;
  vector<int> outer_GameState.x, outer_y;
  vector<int> outer_GameState.R, outer_D;
  int outer_viewOrder = 0;
  vector<vector<int>> outer_moves;
  vector<int> outer_moveCnt;
  array<array<int, 100>, 100> outer_cellUse;
  vector<vector<int>> outer_udlr;
  vector<int> outer_parent;
  vector<int> outer_unionSize;
  set<P> outer_vp;

  GameState outer_GameState;
  // その他
  int K100;
  int methodCount[20][2];
  int methodSum[2];
  int outer_Split = 1;
  vector<int> visited;
  int visitedCnt;
  vector<int> que;

  vector<vector<int>> keep_vp;
  int vpcnt;
  vector<vector<int>> keep_parent;
  int parentcnt;
  vector<vector<int>> keep_udlr;
  int udlrcnt;
  vector<vector<int>> keepA;
  int acnt;
  vector<vector<int>> vv;

}  // namespace

inline bool IsNG(int xx, int yy)
{
  if (xx < 0 || n <= xx || yy < 0 || n <= yy) {
    return true;
  }
  return false;
}

inline bool HasServer(int xx, int yy)
{
  if (0 <= a[xx][yy] && a[xx][yy] < K100) {
    return true;
  }
  return false;
}

inline int MakeAValue(int ite1, int ite2)
{
  int num = -1 * (ite1 * 1000 + ite2);
  if (ite2 < ite1) {
    num = -1 * (ite2 * 1000 + ite1);
  }
  return num;
}

int GetIte(int xx, int yy, char cc)
{
  if (cc == 'U') {
    xx--;
    while (xx >= 0 && !HasServer(xx, yy)) {
      xx--;
    }
    if (xx < 0) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'D') {
    xx++;
    while (xx < n && !HasServer(xx, yy)) {
      xx++;
    }
    if (xx >= n) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'L') {
    yy--;
    while (yy >= 0 && !HasServer(xx, yy)) {
      yy--;
    }
    if (yy < 0) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'R') {
    yy++;
    while (yy < n && !HasServer(xx, yy)) {
      yy++;
    }
    if (yy >= n) {
      return -1;
    }
    return a[xx][yy];
  }

  return -1;
}

// サーバーまでたどり着けるか
int GetIte2(int xx, int yy, char cc)
{
  if (cc == 'U') {
    xx--;
    while (xx >= 0 && !HasServer(xx, yy)) {
      if (a[xx][yy] != INF) {
        return -1;
      }
      xx--;
    }
    if (xx < 0) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'D') {
    xx++;
    while (xx < n && !HasServer(xx, yy)) {
      if (a[xx][yy] != INF) {
        return -1;
      }
      xx++;
    }
    if (xx >= n) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'L') {
    yy--;
    while (yy >= 0 && !HasServer(xx, yy)) {
      if (a[xx][yy] != INF) {
        return -1;
      }
      yy--;
    }
    if (yy < 0) {
      return -1;
    }
    return a[xx][yy];
  }

  if (cc == 'R') {
    yy++;
    while (yy < n && !HasServer(xx, yy)) {
      if (a[xx][yy] != INF) {
        return -1;
      }
      yy++;
    }
    if (yy >= n) {
      return -1;
    }
    return a[xx][yy];
  }

  return -1;
}

void MethodCountReset()
{
  rep(i, 20)
  {
    rep(j, 2) { methodCount[i][j] = 0; }
  }
  methodSum[0] = 0;
  methodSum[1] = 0;
}

void UpdateR(int i)
{
  int now = -1;
  rep(j, n)
  {
    if (0 <= a[i][j] && a[i][j] < K100) {
      R[a[i][j]] = -1;
      if (now != -1) {
        R[now] = a[i][j];
      }
      now = a[i][j];
    }
  }
}

void UpdateD(int j)
{
  int now = -1;
  rep(i, n)
  {
    if (0 <= a[i][j] && a[i][j] < K100) {
      D[a[i][j]] = -1;
      if (now != -1) {
        D[now] = a[i][j];
      }
      now = a[i][j];
    }
  }
}

// スコア計算
// 要素数の大きい順に使う
int CalcScore(int times, bool makeAns = false)
{
  int res = 0;

  {
    int nokori = times;

    for (auto&& p : vp) {
      if (nokori == 0) {
        break;
      }
      int countSize = -p.first;
      int ite = p.second;
      if (nokori >= countSize - 1) {
        res += countSize * (countSize - 1) / 2;
        nokori -= countSize - 1;
      }
      else {
        res += (nokori + 1) * nokori / 2;
        nokori = 0;
      }

      if (nokori == 0) {
        break;
      }
    }
  }

  if (makeAns) {
    visitedCnt++;

    ope2 = 0;

    for (auto&& p : vp) {
      if (times == 0) {
        break;
      }
      int countSize = -p.first;
      int ite = p.second;

      visited[ite] = visitedCnt;
      int queL = 0;
      int queR = 0;
      que[queR] = ite;
      queR++;
      while (queL < queR) {
        int ite = que[queL];
        queL++;
        rep(j, 4)
        {
          int nxt = udlr[ite][j];
          if (nxt != -1 && visited[nxt] != visitedCnt) {
            que[queR] = nxt;
            queR++;
            visited[nxt] = visitedCnt;

            ans2[ope2][0] = x[ite];
            ans2[ope2][1] = y[ite];
            ans2[ope2][2] = x[nxt];
            ans2[ope2][3] = y[nxt];

            ope2++;
            if (ope2 == times) {
              break;
            }
          }
        }
        if (ope2 == times) {
          break;
        }
      }
      if (ope2 == times) {
        break;
      }
    }
  }

  return res;
}

void Init()
{
  visitedCnt = 0;

  // vectorの初期化
  // 最大K=5なので、5*100=500が最大コンピュータ数
  const int maxComputers = 500;
  const int maxOperations = 500;

  x.resize(maxComputers);
  y.resize(maxComputers);
  R.resize(maxComputers);
  D.resize(maxComputers);
  parent.resize(maxComputers);
  unionSize.resize(maxComputers);
  moveCnt.resize(maxComputers);
  visited.resize(maxComputers);
  que.resize(maxComputers);
  ans1.resize(maxOperations, vector<int>(5));
  ans2.resize(maxOperations, vector<int>(4));
  moves.resize(maxComputers, vector<int>(maxOperations));
  udlr.resize(maxComputers, vector<int>(4));

  real_GameState.x.resize(maxComputers);
  real_GameState.y.resize(maxComputers);
  real_GameState.R.resize(maxComputers);
  real_GameState.D.resize(maxComputers);
  real_GameState.parent.resize(maxComputers);
  real_GameState.unionSize.resize(maxComputers);
  real_GameState.moveCnt.resize(maxComputers);
  real_GameState.ans1.resize(maxOperations, vector<int>(5));
  real_GameState.ans2.resize(maxOperations, vector<int>(4));
  real_GameState.moves.resize(maxComputers, vector<int>(maxOperations));
  real_GameState.udlr.resize(maxComputers, vector<int>(4));

  seed_x.resize(maxComputers);
  seed_y.resize(maxComputers);
  seed_R.resize(maxComputers);
  seed_D.resize(maxComputers);
  seed_parent.resize(maxComputers);
  seed_unionSize.resize(maxComputers);
  seed_moveCnt.resize(maxComputers);
  seed_ans1.resize(maxOperations, vector<int>(5));
  seed_ans2.resize(maxOperations, vector<int>(4));
  seed_moves.resize(maxComputers, vector<int>(maxOperations));
  seed_udlr.resize(maxComputers, vector<int>(4));

  outer_x.resize(maxComputers);
  outer_y.resize(maxComputers);
  outer_R.resize(maxComputers);
  outer_D.resize(maxComputers);
  outer_parent.resize(maxComputers);
  outer_unionSize.resize(maxComputers);
  outer_moveCnt.resize(maxComputers);
  outer_ans1.resize(maxOperations, vector<int>(5));
  outer_ans2.resize(maxOperations, vector<int>(4));
  outer_moves.resize(maxComputers, vector<int>(maxOperations));
  outer_udlr.resize(maxComputers, vector<int>(4));

  keep_vp.resize(10000, vector<int>(3));
  keep_parent.resize(10000, vector<int>(2));
  keep_udlr.resize(10000, vector<int>(3));
  keepA.resize(10000, vector<int>(3));
  vv.resize(maxComputers);

  // a,x,y
  int cnt[5] = {};
  rep(i, n)
  {
    rep(j, n)
    {
      int val = s[i][j] - '0' - 1;
      if (val != -1) {
        x[val * 100 + cnt[val]] = i;
        y[val * 100 + cnt[val]] = j;
        a[i][j] = val * 100 + cnt[val];
        cnt[val]++;
      }
      else {
        a[i][j] = INF;
      }
    }
  }

  // cellUse
  rep(i, n)
  {
    rep(j, n)
    {
      cellUse[i][j] = 0;
      if (a[i][j] != INF) {
        cellUse[i][j] = 1;
      }
    }
  }

  // R,D
  rep(i, n) { UpdateR(i); }
  rep(j, n) { UpdateD(j); }

  // udlr
  rep(i, K100)
  {
    rep(j, 4) { udlr[i][j] = -1; }
  }

  if (viewOrder == 0) {
    // 横縦の順
    rep(i, K100)
    {
      if (R[i] == -1) {
        continue;
      }
      if (i / 100 == R[i] / 100) {
        udlr[i][3] = R[i];
        udlr[R[i]][2] = i;
        int aVal = MakeAValue(i, R[i]);
        srep(k, y[i] + 1, y[R[i]]) { a[x[i]][k] = aVal; }
      }
    }

    rep(i, K100)
    {
      if (D[i] == -1) {
        continue;
      }
      if (i / 100 == D[i] / 100) {
        int ok = 1;
        srep(k, x[i] + 1, x[D[i]])
        {
          if (a[k][y[i]] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          udlr[i][1] = D[i];
          udlr[D[i]][0] = i;
          int aVal = MakeAValue(i, D[i]);
          srep(k, x[i] + 1, x[D[i]]) { a[k][y[i]] = aVal; }
        }
      }
    }
  }
  else {
    // 縦横の順
    rep(i, K100)
    {
      if (D[i] == -1) {
        continue;
      }
      if (i / 100 == D[i] / 100) {
        udlr[i][1] = D[i];
        udlr[D[i]][0] = i;
        int aVal = MakeAValue(i, D[i]);
        srep(k, x[i] + 1, x[D[i]]) { a[k][y[i]] = aVal; }
      }
    }

    rep(i, K100)
    {
      if (R[i] == -1) {
        continue;
      }
      if (i / 100 == R[i] / 100) {
        int ok = 1;
        srep(k, y[i] + 1, y[R[i]])
        {
          if (a[x[i]][k] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          udlr[i][3] = R[i];
          udlr[R[i]][2] = i;
          int aVal = MakeAValue(i, R[i]);
          srep(k, y[i] + 1, y[R[i]]) { a[x[i]][k] = aVal; }
        }
      }
    }
  }

  ope1 = 0;
  ope2 = 0;
  maxScore = 0;
  rep(i, K * 100) { moveCnt[i] = 0; }

  // parent, vp
  vp.clear();
  visitedCnt++;
  rep(i, K100) { visited[i] = -1; }

  rep(i, K100) { unionSize[i] = 0; }

  rep(i, K100)
  {
    if (visited[i] == visitedCnt) continue;

    int queL = 0;
    int queR = 0;

    que[queR] = i;
    queR++;
    visited[i] = visitedCnt;
    parent[i] = i;

    int countSize = 1;
    while (queL < queR) {
      int ite = que[queL];
      queL++;
      rep(j, 4)
      {
        int nxt = udlr[ite][j];
        if (nxt != -1 && visited[nxt] != visitedCnt) {
          que[queR] = nxt;
          queR++;
          visited[nxt] = visitedCnt;
          parent[nxt] = i;
          countSize++;
        }
      }
    }

    if (countSize >= 2) {
      vp.insert(P(-countSize, i));
      unionSize[i] = countSize;
    }
  }
}

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

  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    cin >> n >> K;
    rep(i, n) cin >> s[i];
  }
  else {  // ファイル入力する
    ifs >> n >> K;
    rep(i, n) ifs >> s[i];
  }

  K100 = K * 100;

  Init();
}

void Output(int mode, int problemNum)
{
  if (mode == 0) {
    cout << ope1 << endl;
    rep(i, ope1)
    {
      rep(j, 4) { cout << ans1[i][j] << ' '; }
      cout << endl;
    }

    cout << ope2 << endl;
    rep(i, ope2)
    {
      rep(j, 4) { cout << ans2[i][j] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    ofs << ope1 << endl;
    rep(i, ope1)
    {
      rep(j, 4) { ofs << ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << ope2 << endl;
    rep(i, ope2)
    {
      rep(j, 4) { ofs << ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs << endl;

    rep(i, ope1 + 1)
    {
      ofs << i << endl;
      rep(j, i)
      {
        rep(k, 4) { ofs << ans1[j][k] << ' '; }
        ofs << endl;
      }
      ofs << 0 << endl;
      ofs << endl;
    }

    ofs << ope1 << endl;
    rep(i, ope1)
    {
      rep(j, 4) { ofs << ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << ope2 << endl;
    rep(i, ope2)
    {
      rep(j, 4) { ofs << ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs.close();
  }
}

void NormalClear()
{
  viewOrder = 0;
  maxScore = 0;
  ope1 = 0;
  ope2 = 0;
}

void RealClear()
{
  real_GameState.viewOrder = 0;
  real_GameState.maxScore = 0;
  real_GameState.ope1 = 0;
  real_GameState.ope2 = 0;
}

void SeedClear()
{
  seed_viewOrder = 0;
  seed_maxScore = 0;
  seed_ope1 = 0;
  seed_ope2 = 0;
}

void OuterClear()
{
  outer_viewOrder = 0;
  outer_maxScore = 0;
  outer_ope1 = 0;
  outer_ope2 = 0;
}

// maxとreal_maxを初期化
void AllClear()
{
  NormalClear();
  RealClear();
  MethodCountReset();
}

void AllClear_seed()
{
  NormalClear();
  RealClear();
  SeedClear();
  MethodCountReset();
}

void AllClear_outer()
{
  NormalClear();
  RealClear();
  SeedClear();
  OuterClear();
  MethodCountReset();
}

void CopyToReal()
{
  real_GameState.viewOrder = viewOrder;
  real_GameState.maxScore = maxScore;
  real_GameState.ope1 = ope1;
  real_GameState.ope2 = ope2;
  rep(i, ope1)
  {
    rep(j, 5) { real_GameState.ans1[i][j] = ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { real_GameState.ans2[i][j] = ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { real_GameState.a[i][j] = a[i][j]; }
  }
  rep(i, K100)
  {
    real_GameState.x[i] = x[i];
    real_GameState.y[i] = y[i];
    real_GameState.R[i] = R[i];
    real_GameState.D[i] = D[i];
  }

  rep(i, n)
  {
    rep(j, n) { real_GameState.cellUse[i][j] = cellUse[i][j]; }
  }

  rep(i, K100)
  {
    real_GameState.moveCnt[i] = moveCnt[i];
    rep(j, moveCnt[i]) { real_GameState.moves[i][j] = moves[i][j]; }
    rep(j, 4) { real_GameState.udlr[i][j] = udlr[i][j]; }
  }

  rep(i, K100)
  {
    real_GameState.parent[i] = parent[i];
    real_GameState.unionSize[i] = unionSize[i];
  }
  real_GameState.vp = vp;
}

void CopyToSeed()
{
  seed_GameState.viewOrder = viewOrder;
  seed_GameState.maxScore = maxScore;
  seed_GameState.ope1 = ope1;
  seed_GameState.ope2 = ope2;
  rep(i, ope1)
  {
    rep(j, 5) { seed_GameState.ans1[i][j] = ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { seed_GameState.ans2[i][j] = ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { seed_GameState.a[i][j] = a[i][j]; }
  }
  rep(i, K100)
  {
    seed_GameState.x[i] = x[i];
    seed_GameState.y[i] = y[i];
    seed_GameState.R[i] = R[i];
    seed_GameState.D[i] = D[i];
  }

  rep(i, n)
  {
    rep(j, n) { seed_GameState.cellUse[i][j] = cellUse[i][j]; }
  }

  rep(i, K100)
  {
    seed_GameState.moveCnt[i] = moveCnt[i];
    rep(j, moveCnt[i]) { seed_GameState.moves[i][j] = moves[i][j]; }
    rep(j, 4) { seed_GameState.udlr[i][j] = udlr[i][j]; }
  }

  rep(i, K100)
  {
    seed_GameState.parent[i] = parent[i];
    seed_GameState.unionSize[i] = unionSize[i];
  }
  seed_GameState.vp = vp;
}

void CopyToOuter()
{
  outer_GameState.viewOrder = viewOrder;
  outer_GameState.maxScore = maxScore;
  outer_GameState.ope1 = ope1;
  outer_GameState.ope2 = ope2;
  rep(i, ope1)
  {
    rep(j, 5) { outer_GameState.ans1[i][j] = ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { outer_GameState.ans2[i][j] = ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { outer_GameState.a[i][j] = a[i][j]; }
  }
  rep(i, K100)
  {
    outer_GameState.x[i] = x[i];
    outer_GameState.y[i] = y[i];
    outer_GameState.R[i] = R[i];
    outer_GameState.D[i] = D[i];
  }

  rep(i, n)
  {
    rep(j, n) { outer_GameState.cellUse[i][j] = cellUse[i][j]; }
  }

  rep(i, K100)
  {
    outer_GameState.moveCnt[i] = moveCnt[i];
    rep(j, moveCnt[i]) { outer_GameState.moves[i][j] = moves[i][j]; }
    rep(j, 4) { outer_GameState.udlr[i][j] = udlr[i][j]; }
  }

  rep(i, K100)
  {
    outer_GameState.parent[i] = parent[i];
    outer_GameState.unionSize[i] = unionSize[i];
  }
  outer_GameState.vp = vp;
}

void RollBackFromReal()
{
  viewOrder = real_GameState.viewOrder;
  maxScore = real_GameState.maxScore;
  ope1 = real_GameState.ope1;
  ope2 = real_GameState.ope2;
  rep(i, ope1)
  {
    rep(j, 5) { ans1[i][j] = real_GameState.ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { ans2[i][j] = real_GameState.ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { a[i][j] = real_GameState.a[i][j]; }
  }
  rep(i, K100)
  {
    x[i] = real_GameState.x[i];
    y[i] = real_GameState.y[i];
    R[i] = real_GameState.R[i];
    D[i] = real_GameState.D[i];
  }

  rep(i, n)
  {
    rep(j, n) { cellUse[i][j] = real_GameState.cellUse[i][j]; }
  }

  rep(i, K100)
  {
    moveCnt[i] = real_GameState.moveCnt[i];
    rep(j, moveCnt[i]) { moves[i][j] = real_GameState.moves[i][j]; }
    rep(j, 4) { udlr[i][j] = real_GameState.udlr[i][j]; }
  }

  rep(i, K100)
  {
    parent[i] = real_GameState.parent[i];
    unionSize[i] = real_GameState.unionSize[i];
  }
  vp = real_GameState.vp;
}

void RollBackFromSeed()
{
  viewOrder = seed_viewOrder;
  maxScore = seed_maxScore;
  ope1 = seed_ope1;
  ope2 = seed_ope2;
  rep(i, ope1)
  {
    rep(j, 5) { ans1[i][j] = seed_GameState.ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { ans2[i][j] = seed_GameState.ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { a[i][j] = seed_GameState.a[i][j]; }
  }
  rep(i, K100)
  {
    x[i] = seed_GameState.x[i];
    y[i] = seed_GameState.y[i];
    R[i] = seed_GameState.R[i];
    D[i] = seed_GameState.D[i];
  }

  rep(i, n)
  {
    rep(j, n) { cellUse[i][j] = seed_GameState.cellUse[i][j]; }
  }

  rep(i, K100)
  {
    moveCnt[i] = seed_GameState.moveCnt[i];
    rep(j, moveCnt[i]) { moves[i][j] = seed_GameState.moves[i][j]; }
    rep(j, 4) { udlr[i][j] = seed_GameState.udlr[i][j]; }
  }

  rep(i, K100)
  {
    parent[i] = seed_GameState.parent[i];
    unionSize[i] = seed_GameState.unionSize[i];
  }
  vp = seed_vp;
}

void RollBackFromOuter()
{
  viewOrder = outer_viewOrder;
  maxScore = outer_maxScore;
  ope1 = outer_ope1;
  ope2 = outer_ope2;
  rep(i, ope1)
  {
    rep(j, 5) { ans1[i][j] = outer_GameState.ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { ans2[i][j] = outer_GameState.ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { a[i][j] = outer_GameState.a[i][j]; }
  }
  rep(i, K100)
  {
    x[i] = outer_GameState.x[i];
    y[i] = outer_GameState.y[i];
    R[i] = outer_GameState.R[i];
    D[i] = outer_GameState.D[i];
  }

  rep(i, n)
  {
    rep(j, n) { cellUse[i][j] = outer_GameState.cellUse[i][j]; }
  }

  rep(i, K100)
  {
    moveCnt[i] = outer_GameState.moveCnt[i];
    rep(j, moveCnt[i]) { moves[i][j] = outer_GameState.moves[i][j]; }
    rep(j, 4) { udlr[i][j] = outer_GameState.udlr[i][j]; }
  }

  rep(i, K100)
  {
    parent[i] = outer_GameState.parent[i];
    unionSize[i] = outer_GameState.unionSize[i];
  }
  vp = outer_vp;
}

// コンピュータをランダムに1マス移動
void PushACnt(int xx, int yy)
{
  keepA[acnt][0] = xx;
  keepA[acnt][1] = yy;
  keepA[acnt][2] = a[xx][yy];
  acnt++;
}
void BackA()
{
  drep(i, acnt) { a[keepA[i][0]][keepA[i][1]] = keepA[i][2]; }
  acnt = 0;
}

void Push_udlr(int ite, int dir)
{
  keep_udlr[udlrcnt][0] = ite;
  keep_udlr[udlrcnt][1] = dir;
  keep_udlr[udlrcnt][2] = udlr[ite][dir];
  udlrcnt++;
}
void Update_udlr(int ite, int dir, int val)
{
  Push_udlr(ite, dir);
  udlr[ite][dir] = val;
}
void Back_udlr()
{
  drep(i, udlrcnt) { udlr[keep_udlr[i][0]][keep_udlr[i][1]] = keep_udlr[i][2]; }
  udlrcnt = 0;
}

void PushParent(int ite)
{
  keep_parent[parentcnt][0] = ite;
  keep_parent[parentcnt][1] = parent[ite];
  parentcnt++;
}
void BackParent()
{
  drep(i, parentcnt) { parent[keep_parent[i][0]] = keep_parent[i][1]; }
  parentcnt = 0;
}

void PushVp(int val, int ite, int pushpop)
{
  keep_vp[vpcnt][0] = val;
  keep_vp[vpcnt][1] = ite;
  keep_vp[vpcnt][2] = pushpop;
  vpcnt++;
}
void BackVp()
{
  drep(i, vpcnt)
  {
    if (keep_vp[i][2] == 0) {
      // pushされたのでpopする
      vp.erase(P(keep_vp[i][0], keep_vp[i][1]));
      unionSize[keep_vp[i][1]] = 0;
    }
    else {
      vp.insert(P(keep_vp[i][0], keep_vp[i][1]));
      unionSize[keep_vp[i][1]] = -keep_vp[i][0];
    }
  }
  vpcnt = 0;
}

inline void EraseUnion(int ite)
{
  if (unionSize[parent[ite]] != 0) {
    vp.erase(P(-unionSize[parent[ite]], parent[ite]));
    PushVp(-unionSize[parent[ite]], parent[ite], 1);
    unionSize[parent[ite]] = 0;
  }
}

// 戻り値：更新したかどうか
int InnerMethod(double start_temp, double end_temp, double now_progress,
  int ite, int dir, bool forceDo = false, int MethodeMode = 0)
{
  acnt = 0;
  udlrcnt = 0;
  parentcnt = 0;
  vpcnt = 0;

  int nx = x[ite] + dx[dir];
  int ny = y[ite] + dy[dir];

  int xx = x[ite];
  int yy = y[ite];

  int na = a[nx][ny];

  PushACnt(xx, yy);
  PushACnt(nx, ny);

  set<int> se;     // 集合を再構成したい頂点を保持
  vector<P> beam;  // 新たに繋げるか見る頂点

  /*
    左右に動くときは左右との結合関係は保持
    次に動いた2マスの縦を優先
    最後に上下の行の横
  */
  if (nx == xx) {
    // 右
    if (ny == yy + 1) {
      // 元のマスの左のつながり
      // 繋がっている
      if (udlr[ite][2] != -1) {
        if (HasServer(xx, yy - 1)) {
          a[xx][yy] = MakeAValue(ite, a[xx][yy - 1]);
        }
        else {
          a[xx][yy] = a[xx][yy - 1];
        }
      }
      else {
        a[xx][yy] = INF;
      }

      // 元のマスの右のつながり
      a[xx][ny] = ite;

      // 元のマスの上下のつながり
      // 左と繋がっている場合、切る
      if (udlr[ite][2] != -1) {
        if (udlr[ite][0] != -1) {
          int iteU = udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }

        if (udlr[ite][1] != -1) {
          int iteD = udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
      }
      else {
        // 上下繋がっていた
        if (udlr[ite][0] != -1 && udlr[ite][1] != -1) {
          // 繋ぎなおす
          int iteU = udlr[ite][0];
          int iteD = udlr[ite][1];
          Update_udlr(iteU, 1, iteD);
          Update_udlr(iteD, 0, iteU);
          Update_udlr(ite, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          int aVal = MakeAValue(iteU, iteD);
          srep(i, x[iteU] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = aVal;
          }
        }
        // 上繋がっていた
        else if (udlr[ite][0] != -1) {
          // 切る
          int iteU = udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 下繋がっていた
        else if (udlr[ite][1] != -1) {
          // 切る
          int iteD = udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 繋がりなし
        else {
          // 上と下が同じ色の場合
          int iteU = GetIte2(xx, yy, 'U');
          int iteD = GetIte2(xx, yy, 'D');
          if (iteU != -1 && iteD != -1 && iteU / 100 == iteD / 100) {
            // 繋ぐ
            Update_udlr(iteU, 1, iteD);
            Update_udlr(iteD, 0, iteU);

            EraseUnion(iteU);
            EraseUnion(iteD);
            se.insert(iteU);

            int aVal = MakeAValue(iteU, iteD);
            srep(i, x[iteU] + 1, x[iteD])
            {
              PushACnt(i, yy);
              a[i][yy] = aVal;
            }
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && udlr[ite][3] == -1) {
        int iteU = -na / 1000;
        int iteD = -na % 1000;
        if (x[iteU] > x[iteD]) {
          swap(iteU, iteD);
        }

        // 同じ色
        if (iteU / 100 == ite / 100) {
          // 繋ぎなおす
          Update_udlr(iteU, 1, ite);
          Update_udlr(iteD, 0, ite);
          Update_udlr(ite, 0, iteU);
          Update_udlr(ite, 1, iteD);

          EraseUnion(ite);
          EraseUnion(iteU);

          se.insert(ite);

          int aVal = MakeAValue(iteU, ite);
          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }

          aVal = MakeAValue(ite, iteD);
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }
        // 違う色
        else {
          // 切る
          Update_udlr(iteU, 1, -1);
          Update_udlr(iteD, 0, -1);

          EraseUnion(iteU);
          se.insert(iteU);
          se.insert(iteD);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
        }
      }
      // 繋がっていない
      else {
        // 上と繋げられるかどうか
        int iteU = GetIte2(xx, ny, 'U');
        if (iteU != -1 && iteU / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(iteU, 1, ite);
          Update_udlr(ite, 0, iteU);

          EraseUnion(iteU);
          EraseUnion(ite);
          se.insert(iteU);

          int aVal = MakeAValue(iteU, ite);
          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }

        // 下と繋げられるかどうか
        int iteD = GetIte2(xx, ny, 'D');
        if (iteD != -1 && iteD / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(ite, 1, iteD);
          Update_udlr(iteD, 0, ite);

          EraseUnion(ite);
          EraseUnion(iteD);
          se.insert(ite);

          int aVal = MakeAValue(ite, iteD);
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }
      }
    }
    // 左
    else {
      // 元のマスの右のつながり
      // 繋がっている
      if (udlr[ite][3] != -1) {
        if (HasServer(xx, yy + 1)) {
          a[xx][yy] = MakeAValue(ite, a[xx][yy + 1]);
        }
        else {
          a[xx][yy] = a[xx][yy + 1];
        }
      }
      else {
        a[xx][yy] = INF;
      }

      // 元のマスの左のつながり
      a[xx][ny] = ite;

      // 元のマスの上下のつながり
      // 右と繋がっている場合、切る
      if (udlr[ite][3] != -1) {
        if (udlr[ite][0] != -1) {
          int iteU = udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }

        if (udlr[ite][1] != -1) {
          int iteD = udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
      }
      else {
        // 上下繋がっていた
        if (udlr[ite][0] != -1 && udlr[ite][1] != -1) {
          // 繋ぎなおす
          int iteU = udlr[ite][0];
          int iteD = udlr[ite][1];
          Update_udlr(iteU, 1, iteD);
          Update_udlr(iteD, 0, iteU);
          Update_udlr(ite, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          int aVal = MakeAValue(iteU, iteD);
          srep(i, x[iteU] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = aVal;
          }
        }
        // 上繋がっていた
        else if (udlr[ite][0] != -1) {
          // 切る
          int iteU = udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(iteU);
          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 下繋がっていた
        else if (udlr[ite][1] != -1) {
          // 切る
          int iteD = udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(iteD);
          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, yy);
            a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 繋がりなし
        else {
          // 上と下が同じ色の場合
          int iteU = GetIte2(xx, yy, 'U');
          int iteD = GetIte2(xx, yy, 'D');
          if (iteU != -1 && iteD != -1 && iteU / 100 == iteD / 100) {
            // 繋ぐ
            Update_udlr(iteU, 1, iteD);
            Update_udlr(iteD, 0, iteU);

            EraseUnion(iteU);
            EraseUnion(iteD);
            se.insert(iteU);

            int aVal = MakeAValue(iteU, iteD);
            srep(i, x[iteU] + 1, x[iteD])
            {
              PushACnt(i, yy);
              a[i][yy] = aVal;
            }
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && udlr[ite][2] == -1) {
        int iteU = -na / 1000;
        int iteD = -na % 1000;
        if (x[iteU] > x[iteD]) {
          swap(iteU, iteD);
        }

        // 同じ色
        if (iteU / 100 == ite / 100) {
          // 繋ぎなおす
          Update_udlr(iteU, 1, ite);
          Update_udlr(iteD, 0, ite);
          Update_udlr(ite, 0, iteU);
          Update_udlr(ite, 1, iteD);

          EraseUnion(ite);
          EraseUnion(iteU);

          se.insert(ite);

          int aVal = MakeAValue(iteU, ite);
          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }

          aVal = MakeAValue(ite, iteD);
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }
        // 違う色
        else {
          // 切る
          Update_udlr(iteU, 1, -1);
          Update_udlr(iteD, 0, -1);

          EraseUnion(iteU);
          se.insert(iteU);
          se.insert(iteD);

          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
        }
      }
      // 繋がっていない
      else {
        // 上と繋げられるかどうか
        int iteU = GetIte2(xx, ny, 'U');
        if (iteU != -1 && iteU / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(iteU, 1, ite);
          Update_udlr(ite, 0, iteU);

          EraseUnion(iteU);
          EraseUnion(ite);
          se.insert(iteU);

          int aVal = MakeAValue(iteU, ite);
          srep(i, x[iteU] + 1, x[ite])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }

        // 下と繋げられるかどうか
        int iteD = GetIte2(xx, ny, 'D');
        if (iteD != -1 && iteD / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(ite, 1, iteD);
          Update_udlr(iteD, 0, ite);

          EraseUnion(ite);
          EraseUnion(iteD);
          se.insert(ite);

          int aVal = MakeAValue(ite, iteD);
          srep(i, x[ite] + 1, x[iteD])
          {
            PushACnt(i, ny);
            a[i][ny] = aVal;
          }
        }
      }
    }
  }
  else {
    // 下
    if (nx == xx + 1) {
      // 元のマスの上のつながり
      // 繋がっている
      if (udlr[ite][0] != -1) {
        if (HasServer(xx - 1, yy)) {
          a[xx][yy] = MakeAValue(ite, a[xx - 1][yy]);
        }
        else {
          a[xx][yy] = a[xx - 1][yy];
        }
      }
      else {
        a[xx][yy] = INF;
      }

      // 元のマスの下のつながり
      a[nx][yy] = ite;

      // 元のマスの左右のつながり
      // 上と繋がっている場合、切る
      if (udlr[ite][0] != -1) {
        if (udlr[ite][2] != -1) {
          int iteL = udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }

        if (udlr[ite][3] != -1) {
          int iteR = udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
      }
      else {
        // 左右繋がっていた
        if (udlr[ite][2] != -1 && udlr[ite][3] != -1) {
          // 繋ぎなおす
          int iteL = udlr[ite][2];
          int iteR = udlr[ite][3];
          Update_udlr(iteL, 3, iteR);
          Update_udlr(iteR, 2, iteL);
          Update_udlr(ite, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          int aVal = MakeAValue(iteL, iteR);
          srep(i, y[iteL] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = aVal;
          }
        }
        // 左繋がっていた
        else if (udlr[ite][2] != -1) {
          // 切る
          int iteL = udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 右繋がっていた
        else if (udlr[ite][3] != -1) {
          // 切る
          int iteR = udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(iteR);
          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 繋がりなし
        else {
          // 左と右が同じ色の場合
          int iteL = GetIte2(xx, yy, 'L');
          int iteR = GetIte2(xx, yy, 'R');
          if (iteL != -1 && iteR != -1 && iteL / 100 == iteR / 100) {
            // 繋ぐ
            Update_udlr(iteL, 3, iteR);
            Update_udlr(iteR, 2, iteL);

            EraseUnion(iteL);
            EraseUnion(iteR);
            se.insert(iteL);

            int aVal = MakeAValue(iteL, iteR);
            srep(i, y[iteL] + 1, y[iteR])
            {
              PushACnt(xx, i);
              a[xx][i] = aVal;
            }
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && udlr[ite][1] == -1) {
        int iteL = -na / 1000;
        int iteR = -na % 1000;
        if (y[iteL] > y[iteR]) {
          swap(iteL, iteR);
        }

        // 同じ色
        if (iteL / 100 == ite / 100) {
          // 繋ぎなおす
          Update_udlr(iteL, 3, ite);
          Update_udlr(iteR, 2, ite);
          Update_udlr(ite, 2, iteL);
          Update_udlr(ite, 3, iteR);

          EraseUnion(ite);
          EraseUnion(iteL);

          se.insert(ite);

          int aVal = MakeAValue(iteL, ite);
          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }

          aVal = MakeAValue(ite, iteR);
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }
        // 違う色
        else {
          // 切る
          Update_udlr(iteL, 3, -1);
          Update_udlr(iteR, 2, -1);

          EraseUnion(iteL);
          se.insert(iteL);
          se.insert(iteR);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
        }
      }
      // 繋がっていない
      else {
        // 左と繋げられるかどうか
        int iteL = GetIte2(nx, yy, 'L');
        if (iteL != -1 && iteL / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(iteL, 3, ite);
          Update_udlr(ite, 2, iteL);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);

          int aVal = MakeAValue(iteL, ite);
          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }

        // 右と繋げられるかどうか
        int iteR = GetIte2(nx, yy, 'R');
        if (iteR != -1 && iteR / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(ite, 3, iteR);
          Update_udlr(iteR, 2, ite);

          EraseUnion(ite);
          EraseUnion(iteR);
          se.insert(ite);

          int aVal = MakeAValue(ite, iteR);
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }
      }
    }
    // 上
    else {
      // 元のマスの下のつながり
      // 繋がっている
      if (udlr[ite][1] != -1) {
        if (HasServer(xx + 1, yy)) {
          a[xx][yy] = MakeAValue(ite, a[xx + 1][yy]);
        }
        else {
          a[xx][yy] = a[xx + 1][yy];
        }
      }
      else {
        a[xx][yy] = INF;
      }

      // 元のマスの上のつながり
      a[nx][yy] = ite;

      // 元のマスの左右のつながり
      // 下と繋がっている場合、切る
      if (udlr[ite][1] != -1) {
        if (udlr[ite][2] != -1) {
          int iteL = udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }

        if (udlr[ite][3] != -1) {
          int iteR = udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
      }
      else {
        // 左右繋がっていた
        if (udlr[ite][2] != -1 && udlr[ite][3] != -1) {
          // 繋ぎなおす
          int iteL = udlr[ite][2];
          int iteR = udlr[ite][3];
          Update_udlr(iteL, 3, iteR);
          Update_udlr(iteR, 2, iteL);
          Update_udlr(ite, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          int aVal = MakeAValue(iteL, iteR);
          srep(i, y[iteL] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = aVal;
          }
        }
        // 左繋がっていた
        else if (udlr[ite][2] != -1) {
          // 切る
          int iteL = udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 右繋がっていた
        else if (udlr[ite][3] != -1) {
          // 切る
          int iteR = udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(iteR);
          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(xx, i);
            a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 繋がりなし
        else {
          // 左と右が同じ色の場合
          int iteL = GetIte2(xx, yy, 'L');
          int iteR = GetIte2(xx, yy, 'R');
          if (iteL != -1 && iteR != -1 && iteL / 100 == iteR / 100) {
            // 繋ぐ
            Update_udlr(iteL, 3, iteR);
            Update_udlr(iteR, 2, iteL);

            EraseUnion(iteL);
            EraseUnion(iteR);
            se.insert(iteL);

            int aVal = MakeAValue(iteL, iteR);
            srep(i, y[iteL] + 1, y[iteR])
            {
              PushACnt(xx, i);
              a[xx][i] = aVal;
            }
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && udlr[ite][0] == -1) {
        int iteL = -na / 1000;
        int iteR = -na % 1000;
        if (y[iteL] > y[iteR]) {
          swap(iteL, iteR);
        }

        // 同じ色
        if (iteL / 100 == ite / 100) {
          // 繋ぎなおす
          Update_udlr(iteL, 3, ite);
          Update_udlr(iteR, 2, ite);
          Update_udlr(ite, 2, iteL);
          Update_udlr(ite, 3, iteR);

          EraseUnion(ite);
          EraseUnion(iteL);

          se.insert(ite);

          int aVal = MakeAValue(iteL, ite);
          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }

          aVal = MakeAValue(ite, iteR);
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }
        // 違う色
        else {
          // 切る
          Update_udlr(iteL, 3, -1);
          Update_udlr(iteR, 2, -1);

          EraseUnion(iteL);
          se.insert(iteL);
          se.insert(iteR);

          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
        }
      }
      // 繋がっていない
      else {
        // 左と繋げられるかどうか
        int iteL = GetIte2(nx, yy, 'L');
        if (iteL != -1 && iteL / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(iteL, 3, ite);
          Update_udlr(ite, 2, iteL);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);

          int aVal = MakeAValue(iteL, ite);
          srep(i, y[iteL] + 1, y[ite])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }

        // 右と繋げられるかどうか
        int iteR = GetIte2(nx, yy, 'R');
        if (iteR != -1 && iteR / 100 == ite / 100) {
          // 繋ぐ
          Update_udlr(ite, 3, iteR);
          Update_udlr(iteR, 2, ite);

          EraseUnion(ite);
          EraseUnion(iteR);
          se.insert(ite);

          int aVal = MakeAValue(ite, iteR);
          srep(i, y[ite] + 1, y[iteR])
          {
            PushACnt(nx, i);
            a[nx][i] = aVal;
          }
        }
      }
    }
  }

  // 新たな辺の作成
  // 左右
  if (nx == xx) {
    for (auto&& p : beam) {
      int px = p.first;
      int py = p.second;
      if (a[px][py] != INF) continue;
      int iteL = GetIte2(px, py, 'L');
      int iteR = GetIte2(px, py, 'R');
      if (iteL != -1 && iteR != -1 && iteL / 100 == iteR / 100) {
        Update_udlr(iteL, 3, iteR);
        Update_udlr(iteR, 2, iteL);

        EraseUnion(iteL);
        EraseUnion(iteR);
        se.insert(iteL);

        int aVal = MakeAValue(iteL, iteR);
        srep(i, y[iteL] + 1, y[iteR])
        {
          PushACnt(px, i);
          a[px][i] = aVal;
        }
      }
    }
  }
  // 上下
  else {
    for (auto&& p : beam) {
      int px = p.first;
      int py = p.second;
      if (a[px][py] != INF) continue;
      int iteU = GetIte2(px, py, 'U');
      int iteD = GetIte2(px, py, 'D');
      if (iteU != -1 && iteD != -1 && iteU / 100 == iteD / 100) {
        Update_udlr(iteU, 1, iteD);
        Update_udlr(iteD, 0, iteU);

        EraseUnion(iteU);
        EraseUnion(iteD);
        se.insert(iteU);

        int aVal = MakeAValue(iteU, iteD);
        srep(i, x[iteU] + 1, x[iteD])
        {
          PushACnt(i, py);
          a[i][py] = aVal;
        }
      }
    }
  }

  x[ite] = nx;
  y[ite] = ny;



  if (nx != xx) {
    UpdateR(xx);
    UpdateR(nx);
  }


  if (ny != yy) {
    UpdateD(yy);
    UpdateD(ny);
  }


  // parent, vpの更新
  visitedCnt++;
  for (auto&& vpIte : se) {
    if (visited[vpIte] == visitedCnt) {
      continue;
    }

    int queL = 0;
    int queR = 0;

    que[queR] = vpIte;
    queR++;
    visited[vpIte] = visitedCnt;
    PushParent(vpIte);
    parent[vpIte] = vpIte;

    int countSize = 1;
    while (queL < queR) {
      int nowIte = que[queL];
      queL++;
      rep(j, 4)
      {
        int nxt = udlr[nowIte][j];
        if (nxt != -1 && visited[nxt] != visitedCnt) {
          que[queR] = nxt;
          queR++;
          visited[nxt] = visitedCnt;
          PushParent(nxt);
          parent[nxt] = vpIte;
          countSize++;
        }
      }
    }

    if (countSize >= 2) {
      vp.insert(P(-countSize, vpIte));
      PushVp(-countSize, vpIte, 0);
      unionSize[vpIte] = countSize;
    }
  }

  int tmpScore = 0;
  if (MethodeMode == 5) {
    tmpScore = CalcScore(K100 - (ope1 - 1));
  }
  else {
    tmpScore = CalcScore(K100 - (ope1 + 1));
  }

  methodCount[1][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  int isDo = 0;
  if (forceDo || prob > Rand01()) {
    isDo = 1;
    maxScore += diffScore;

    methodCount[1][0]++;
    methodSum[0]++;

    ans1[ope1][0] = xx;
    ans1[ope1][1] = yy;
    ans1[ope1][2] = nx;
    ans1[ope1][3] = ny;
    ans1[ope1][4] = ite;
    ope1++;

    if (maxScore > real_GameState.maxScore) {
      if (MethodeMode == 5) {
        isDo = 5;
      }
      CopyToReal();
    }
  }
  else {
    // 元に戻す
    x[ite] = xx;
    y[ite] = yy;

    BackA();
    Back_udlr();
    BackParent();
    BackVp();

    if (nx != xx) {
      UpdateR(xx);
      UpdateR(nx);
    }


    if (ny != yy) {
      UpdateD(yy);
      UpdateD(ny);
    }
  }

  return isDo;
}

void Method1(double start_temp, double end_temp, double now_progress)
{
  int ite, dir, nx, ny;
  int randCnt = 0;
  while (true) {
    randCnt++;
    if (randCnt == 100) {
      return;
    }
    ite = Rand() % K100;  // 1マス動かすコンピュータ
    dir = Rand() % 4;

    nx = x[ite] + dx[dir];
    ny = y[ite] + dy[dir];
    int ok = 0;
    if (0 <= nx && nx < n && 0 <= ny && ny < n && !HasServer(nx, ny)) {
      break;
    }
  }

  InnerMethod(start_temp, end_temp, now_progress, ite, dir);
}


// 空白を2マス動かす
void Method3(double start_temp, double end_temp, double now_progress)
{
  int xx, yy, dir1, dir2;
  while (true) {
    xx = Rand() % n;
    yy = Rand() % n;
    if (a[xx][yy] == -1) {
      break;
    }
  }

  dir1 = Rand() % 4;
  dir2 = Rand() % 4;
  int nx1 = xx + dx[dir1];
  int ny1 = yy + dy[dir1];
  int nx2 = nx1 + dx[dir2];
  int ny2 = ny1 + dy[dir2];
  if (nx1 < 0 || n <= nx1 || ny1 < 0 || n <= ny1 || a[nx1][ny1] == -1) {
    return;
  }
  if (nx2 < 0 || n <= nx2 || ny2 < 0 || n <= ny2 || a[nx2][ny2] == -1) {
    return;
  }

  int ite1 = a[nx1][ny1];
  std::swap(a[nx1][ny1], a[xx][yy]);
  x[ite1] = xx;
  y[ite1] = yy;

  int ite2 = a[nx2][ny2];
  std::swap(a[nx2][ny2], a[nx1][ny1]);
  x[ite2] = nx1;
  y[ite2] = ny1;

  UpdateR(xx);
  UpdateR(nx1);
  UpdateR(nx2);
  UpdateD(yy);
  UpdateD(ny1);
  UpdateD(ny2);

  int tmpScore = CalcScore(K100 - (ope1 + 2));
  methodCount[3][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    maxScore += diffScore;

    methodCount[3][0]++;
    methodSum[0]++;

    ans1[ope1][0] = nx1;
    ans1[ope1][1] = ny1;
    ans1[ope1][2] = xx;
    ans1[ope1][3] = yy;
    ans1[ope1][4] = ite1;
    ope1++;

    ans1[ope1][0] = nx2;
    ans1[ope1][1] = ny2;
    ans1[ope1][2] = nx1;
    ans1[ope1][3] = ny1;
    ans1[ope1][4] = ite2;
    ope1++;

    if (maxScore > real_GameState.maxScore) {
      CopyToReal();
    }
  }
  else {
    // 元に戻す
    x[ite2] = nx2;
    y[ite2] = ny2;
    std::swap(a[nx2][ny2], a[nx1][ny1]);

    x[ite1] = nx1;
    y[ite1] = ny1;
    std::swap(a[nx1][ny1], a[xx][yy]);

    UpdateR(xx);
    UpdateR(nx1);
    UpdateR(nx2);
    UpdateD(yy);
    UpdateD(ny1);
    UpdateD(ny2);
  }
}

// コンピュータをランダムに2マス移動
void Method4(double start_temp, double end_temp, double now_progress)
{
  int ite = Rand() % K100;  // 1マス動かすコンピュータ
  int dir1 = Rand() % 4;
  int dir2 = Rand() % 4;

  int xx = x[ite];
  int yy = y[ite];

  int nx1 = xx + dx[dir1];
  int ny1 = yy + dy[dir1];
  int nx2 = nx1 + dx[dir2];
  int ny2 = ny1 + dy[dir2];

  if (IsNG(nx1, ny1) || IsNG(nx2, ny2) || a[nx1][ny1] != -1 ||
    a[nx2][ny2] != -1) {
    return;
  }

  swap(a[xx][yy], a[nx2][ny2]);
  x[ite] = nx2;
  y[ite] = ny2;

  UpdateR(xx);
  UpdateR(nx1);
  UpdateR(nx2);
  UpdateD(yy);
  UpdateD(ny1);
  UpdateD(ny2);

  int tmpScore = CalcScore(K100 - (ope1 + 2));
  methodCount[4][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    maxScore += diffScore;

    methodCount[4][0]++;
    methodSum[0]++;

    ans1[ope1][0] = xx;
    ans1[ope1][1] = yy;
    ans1[ope1][2] = nx1;
    ans1[ope1][3] = ny1;
    ans1[ope1][4] = ite;
    ope1++;

    ans1[ope1][0] = nx1;
    ans1[ope1][1] = ny1;
    ans1[ope1][2] = nx2;
    ans1[ope1][3] = ny2;
    ans1[ope1][4] = ite;
    ope1++;

    if (maxScore > real_GameState.maxScore) {
      CopyToReal();
    }
  }
  else {
    // 元に戻す
    swap(a[xx][yy], a[nx2][ny2]);
    x[ite] = xx;
    y[ite] = yy;

    UpdateR(xx);
    UpdateR(nx1);
    UpdateR(nx2);
    UpdateD(yy);
    UpdateD(ny1);
    UpdateD(ny2);
  }
}

// 移動をランダムに1つ削除
void Method5(double start_temp, double end_temp, double now_progress)
{
  if (ope1 == 0) return;
  int ite = Rand() % ope1;

  // NGチェック
  // ite以降の操作で、操作元が移動後のマス、操作後が移動前のマス、の操作が出てこなければOK
  srep(i, ite + 1, ope1)
  {
    if (ans1[i][4] == ans1[ite][4]) {
      return;
    }
    if (ans1[i][0] == ans1[ite][2] && ans1[i][1] == ans1[ite][3]) {
      return;
    }
    if (ans1[i][2] == ans1[ite][0] && ans1[i][3] == ans1[ite][1]) {
      return;
    }
  }

  // 消したい移動の逆の操作を足してInnerMethodを実行
  int reverseDir = 0;
  rep(i, 4)
  {
    if (ans1[ite][2] == ans1[ite][0] + dx[i] &&
      ans1[ite][3] == ans1[ite][1] + dy[i]) {
      reverseDir = (i + 2) % 4;
      break;
    }
  }

  int isDo = InnerMethod(start_temp, end_temp, now_progress, ans1[ite][4],
    reverseDir, false, 5);

  // 実行した場合、2つ消す
  if (isDo) {
    /*
    修正するもの
    int ope1;
    int ans1[1000][5];
    */
    // 消すのはiteとope1-1
    ope1--;
    srep(i, ite, ope1 - 1)
    {
      rep(j, 5) { ans1[i][j] = ans1[i + 1][j]; }
    }
    ope1--;
    if (isDo == 5) {
      CopyToReal();
    }
  }
  else {
    // 何もしない
  }
}

// スワップしてるだけの2つの移動を削除
void Method6(double start_temp, double end_temp, double now_progress)
{
  methodCount[6][1]++;
  methodSum[1]++;

  while (true) {
    int ite = -1;
    rep(i, ope1 - 1)
    {
      if (ans1[i][0] == ans1[i + 1][2] && ans1[i][1] == ans1[i + 1][3] &&
        ans1[i][2] == ans1[i + 1][0] && ans1[i][3] == ans1[i + 1][1]) {
        ite = i;
        break;
      }
    }

    if (ite == -1) {
      return;
    }

    maxScore = CalcScore(K100 - (ope1 - 2));

    methodCount[6][0]++;
    methodSum[0]++;

    srep(i, ite, ope1 - 2)
    {
      rep(j, 5) { ans1[i][j] = ans1[i + 2][j]; }
    }

    ope1 -= 2;

    if (maxScore > real_GameState.maxScore) {
      CopyToReal();
    }
  }
}

// スワップしてるだけの2つの移動を削除
void Method7(double start_temp, double end_temp, double now_progress)
{
  methodCount[7][1]++;
  methodSum[1]++;


  while (true) {
    rep(i, K * 100)
    {
      vv[i].clear();
    }
    rep(i, ope1)
    {
      int ite = ans1[i][4];
      vv[ite].push_back(i);
    }

    int it1 = -1;
    int it2 = -1;
    rep(i, K * 100)
    {
      int sz = vv[i].size();
      drep(j, sz - 1)
      {
        int ite1 = vv[i][j];
        int ite2 = vv[i][j + 1];
        if (ans1[ite1][0] == ans1[ite2][2] && ans1[ite1][1] == ans1[ite2][3] &&
          ans1[ite1][2] == ans1[ite2][0] && ans1[ite1][3] == ans1[ite2][1]) {
          // ngチェック
          int ng = 0;
          srep(k, ite1 + 1, ite2)
          {
            if (ans1[k][2] == ans1[ite1][0] && ans1[k][3] == ans1[ite1][1]) {
              ng = 1;
              break;
            }
          }
          if (ng == 0) {
            it1 = ite1;
            it2 = ite2;
            break;
          }

        }
      }
      if (it1 != -1) {
        break;
      }
    }

    if (it1 == -1) {
      break;
    }

    if (it1 != -1) {
      srep(i, it2, ope1 - 1)
      {
        rep(j, 5)
        {
          ans1[i][j] = ans1[i + 1][j];
        }
      }
      ope1--;
      srep(i, it1, ope1 - 1)
      {
        rep(j, 5)
        {
          ans1[i][j] = ans1[i + 1][j];
        }
      }
      ope1--;

      maxScore = CalcScore(K100 - ope1);
      methodCount[7][0]++;
      methodSum[0]++;

      if (maxScore > real_GameState.maxScore) {
        CopyToReal();
      }

    }
  }
}

int Solve(int mode, int problemNum = 0)
{
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();

  Init();

  // 愚直解
  maxScore = CalcScore(K100, true);
  CopyToReal();
  CopyToSeed();

  // シード作り
  int seedCount = 10;
  rep(tei, seedCount)
  {
    start_time = clock();

    Init();
    viewOrder = tei % 2;
    maxScore = CalcScore(K100, true);

    // 焼きなまし
    end_time = clock();
    double now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double TL = 1.0 / seedCount;
    double now_progress = now_time / TL;
    double start_temp = 10.0 + 10.0 * K100 / (n * n - K100);
    double end_temp = 0;
    int loop = 0;
    int rollbackCount = 0;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        end_time = clock();
        now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        now_progress = now_time / TL;
      }
      if (now_progress > 1.0) break;

      // 現在のスコアが悪いときは元に戻す
      if (maxScore * 1.2 < real_GameState.maxScore || Rand() % 123456 == 0) {
        RollBackFromReal();
        rollbackCount++;
      }

      int me = 1;

      if (Rand() % 2 == 0) {
        me = 5;
      }

      if (Rand() % 203 == 0) {
        me = 6;
      }

      // コンピュータをランダムに1マス移動
      if (me == 1) {
        if (ope1 >= K100) {
          continue;
        }
        Method1(start_temp, end_temp, now_progress);
      }

      // 空白を2マス動かす
      if (me == 3) {
        if (ope1 >= K100 - 1) {
          continue;
        }
        Method3(start_temp, end_temp, now_progress);
      }

      // コンピュータをランダムに2マス移動
      if (me == 4) {
        if (ope1 >= K100 - 1) {
          continue;
        }
        Method4(start_temp, end_temp, now_progress);
      }

      // 移動をランダムに1つ削除
      if (me == 5) {
        if (ope1 == 0) {
          continue;
        }
        Method5(start_temp, end_temp, now_progress);
      }

      // スワップしてるだけの2つの移動を削除
      if (me == 6) {
        Method6(start_temp, end_temp, now_progress);
      }
    }

    // スコアが良ければシードを更新
    RollBackFromReal();
    if (maxScore > seed_GameState.maxScore) {
      CopyToSeed();
    }

    AllClear();
  }

  // シードから戻す
  RollBackFromSeed();
  CopyToReal();

  // 焼きなまし
  start_time = clock();
  end_time = clock();
  double now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  double TL = 1.9 / outer_Split;
  double now_progress = now_time / TL;
  double start_temp = 20.0 + 10.0 * K100 / (n * n - K100);
  double end_temp = 0;
  int loop = 0;
  int rollbackCount = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      end_time = clock();
      now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
      now_progress = now_time / TL;
    }
    if (now_progress > 1.0) break;

    // 現在のスコアが悪いときは元に戻す
    if (maxScore * 1.2 < real_GameState.maxScore || Rand() % 123456 == 0) {
      RollBackFromReal();
      rollbackCount++;
    }

    int me = 1;

    if (Rand() % 2 == 0) {
      me = 5;
    }

    if (Rand() % 203 == 0) {
      me = 6;
    }

    if (Rand() % 1011 == 0) {
      me = 7;
    }

    // コンピュータをランダムに1マス移動
    if (me == 1) {
      if (ope1 >= K100) {
        continue;
      }
      Method1(start_temp, end_temp, now_progress);
    }

    // 空白を2マス動かす
    if (me == 3) {
      if (ope1 >= K100 - 1) {
        continue;
      }
      Method3(start_temp, end_temp, now_progress);
    }

    // コンピュータをランダムに2マス移動
    if (me == 4) {
      if (ope1 >= K100 - 1) {
        continue;
      }
      Method4(start_temp, end_temp, now_progress);
    }

    // 移動をランダムに1つ削除
    if (me == 5) {
      if (ope1 == 0) {
        continue;
      }
      Method5(start_temp, end_temp, now_progress);
    }

    // スワップしてるだけの2つの移動を削除
    if (me == 6) {
      Method6(start_temp, end_temp, now_progress);
    }


    // スワップしてるだけの2つの移動を削除
    if (me == 7) {
      Method7(start_temp, end_temp, now_progress);
    }
  }

  RollBackFromReal();

  int cal = CalcScore(K100 - ope1, true);

  if (true) {
    cerr << "problemNum = " << problemNum << ", N = " << n << ", K = " << K << endl;
    cerr << "start_temp = " << start_temp << ", viewOrder = " << viewOrder << endl;
    cerr << "maxScore = " << maxScore << ", ope1 = " << ope1 << ", ope2 = " << ope2 << endl;
    cerr << "cal = " << cal << endl;
    cerr << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    srep(i, 1, 8)
    {
      cerr << "Method" << i << " = " << methodCount[i][0] << " / "
        << methodCount[i][1] << endl;
    }
    cerr << "MethodSum = " << methodSum[0] << " / " << methodSum[1] << endl;
    cerr << endl;
  }

  return maxScore;
}

int SolveOuter(int mode, int problemNum)
{
  // 入力部
  Input(problemNum);

  rep(_, outer_Split)
  {
    viewOrder = _ % 2;
    int score = Solve(mode, problemNum);
    if (score >= outer_GameState.maxScore) {
      CopyToOuter();
    }
    AllClear();
  }
  RollBackFromOuter();

  // 解の出力
  Output(mode, problemNum);

  return maxScore;
}

int main()
{
  int mode = 2;

  if (mode == 0) {
    SolveOuter(mode, 0);
  }
  else if (mode == 1) {
    SolveOuter(mode, 1);
  }
  else if (mode == 2) {
    int sum = 0;
    srep(i, 0, 10)
    {
      sum += SolveOuter(mode, i);
      AllClear_outer();
    }
    cout << "sum = " << sum << endl;
  }
  else if (mode == 3) {
    int problemNum = 0;
    cin >> problemNum;
    SolveOuter(mode, problemNum);
  }

  return 0;
}
