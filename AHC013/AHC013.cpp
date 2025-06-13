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

  // Member functions
  void clear();
  void copyTo(GameState& target) const;
  void copyFrom(const GameState& source);
  void updateSingleR(int i);  // Update R for a single computer
  void updateSingleD(int j);  // Update D for a single computer
  void updateRowR(int i);     // Update all R values in row i
  void updateColD(int j);     // Update all D values in column j
  int calcScore(int times, bool makeAns = false);
};


namespace /* 変数 */
{
  // 入力用変数
  int n, K;
  array<string, 100> s;

  // 解答用変数
  GameState gameState;
  GameState real_GameState;
  GameState seed_GameState;
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
  if (0 <= gameState.a[xx][yy] && gameState.a[xx][yy] < K100) {
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
    return gameState.a[xx][yy];
  }

  if (cc == 'D') {
    xx++;
    while (xx < n && !HasServer(xx, yy)) {
      xx++;
    }
    if (xx >= n) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  if (cc == 'L') {
    yy--;
    while (yy >= 0 && !HasServer(xx, yy)) {
      yy--;
    }
    if (yy < 0) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  if (cc == 'R') {
    yy++;
    while (yy < n && !HasServer(xx, yy)) {
      yy++;
    }
    if (yy >= n) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  return -1;
}

// サーバーまでたどり着けるか
int GetIte2(int xx, int yy, char cc)
{
  if (cc == 'U') {
    xx--;
    while (xx >= 0 && !HasServer(xx, yy)) {
      if (gameState.a[xx][yy] != INF) {
        return -1;
      }
      xx--;
    }
    if (xx < 0) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  if (cc == 'D') {
    xx++;
    while (xx < n && !HasServer(xx, yy)) {
      if (gameState.a[xx][yy] != INF) {
        return -1;
      }
      xx++;
    }
    if (xx >= n) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  if (cc == 'L') {
    yy--;
    while (yy >= 0 && !HasServer(xx, yy)) {
      if (gameState.a[xx][yy] != INF) {
        return -1;
      }
      yy--;
    }
    if (yy < 0) {
      return -1;
    }
    return gameState.a[xx][yy];
  }

  if (cc == 'R') {
    yy++;
    while (yy < n && !HasServer(xx, yy)) {
      if (gameState.a[xx][yy] != INF) {
        return -1;
      }
      yy++;
    }
    if (yy >= n) {
      return -1;
    }
    return gameState.a[xx][yy];
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


// Forward declaration
bool HasServer(int xx, int yy);

// GameState member function implementations
void GameState::clear()
{
  viewOrder = 0;
  maxScore = 0;
  ope1 = 0;
  ope2 = 0;
}

void GameState::copyTo(GameState& target) const
{
  target.viewOrder = viewOrder;
  target.maxScore = maxScore;
  target.ope1 = ope1;
  target.ope2 = ope2;

  rep(i, ope1)
  {
    rep(j, 5) { target.ans1[i][j] = ans1[i][j]; }
  }
  rep(i, ope2)
  {
    rep(j, 4) { target.ans2[i][j] = ans2[i][j]; }
  }
  rep(i, n)
  {
    rep(j, n) { target.a[i][j] = a[i][j]; }
  }
  rep(i, K100)
  {
    target.x[i] = x[i];
    target.y[i] = y[i];
    target.R[i] = R[i];
    target.D[i] = D[i];
  }

  rep(i, n)
  {
    rep(j, n) { target.cellUse[i][j] = cellUse[i][j]; }
  }

  rep(i, K100)
  {
    target.moveCnt[i] = moveCnt[i];
    rep(j, moveCnt[i]) { target.moves[i][j] = moves[i][j]; }
    rep(j, 4) { target.udlr[i][j] = udlr[i][j]; }
  }

  rep(i, K100)
  {
    target.parent[i] = parent[i];
    target.unionSize[i] = unionSize[i];
  }
  target.vp = vp;
}

void GameState::copyFrom(const GameState& source)
{
  source.copyTo(*this);
}

void GameState::updateSingleR(int i)
{
  int xx = x[i];
  int yy = y[i];
  yy++;
  while (yy < n && !HasServer(xx, yy)) {
    yy++;
  }
  if (yy < n) {
    R[i] = a[xx][yy];
  }
  else {
    R[i] = -1;
  }
}

void GameState::updateSingleD(int j)
{
  int xx = x[j];
  int yy = y[j];
  xx++;
  while (xx < n && !HasServer(xx, yy)) {
    xx++;
  }
  if (xx < n) {
    D[j] = a[xx][yy];
  }
  else {
    D[j] = -1;
  }
}

void GameState::updateRowR(int i)
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

void GameState::updateColD(int j)
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

int GameState::calcScore(int times, bool makeAns)
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

  gameState.x.resize(maxComputers);
  gameState.y.resize(maxComputers);
  gameState.R.resize(maxComputers);
  gameState.D.resize(maxComputers);
  gameState.parent.resize(maxComputers);
  gameState.unionSize.resize(maxComputers);
  gameState.moveCnt.resize(maxComputers);
  visited.resize(maxComputers);
  que.resize(maxComputers);
  gameState.ans1.resize(maxOperations, vector<int>(5));
  gameState.ans2.resize(maxOperations, vector<int>(4));
  gameState.moves.resize(maxComputers, vector<int>(maxOperations));
  gameState.udlr.resize(maxComputers, vector<int>(4));

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

  seed_GameState.x.resize(maxComputers);
  seed_GameState.y.resize(maxComputers);
  seed_GameState.R.resize(maxComputers);
  seed_GameState.D.resize(maxComputers);
  seed_GameState.parent.resize(maxComputers);
  seed_GameState.unionSize.resize(maxComputers);
  seed_GameState.moveCnt.resize(maxComputers);
  seed_GameState.ans1.resize(maxOperations, vector<int>(5));
  seed_GameState.ans2.resize(maxOperations, vector<int>(4));
  seed_GameState.moves.resize(maxComputers, vector<int>(maxOperations));
  seed_GameState.udlr.resize(maxComputers, vector<int>(4));

  outer_GameState.x.resize(maxComputers);
  outer_GameState.y.resize(maxComputers);
  outer_GameState.R.resize(maxComputers);
  outer_GameState.D.resize(maxComputers);
  outer_GameState.parent.resize(maxComputers);
  outer_GameState.unionSize.resize(maxComputers);
  outer_GameState.moveCnt.resize(maxComputers);
  outer_GameState.ans1.resize(maxOperations, vector<int>(5));
  outer_GameState.ans2.resize(maxOperations, vector<int>(4));
  outer_GameState.moves.resize(maxComputers, vector<int>(maxOperations));
  outer_GameState.udlr.resize(maxComputers, vector<int>(4));

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
        gameState.x[val * 100 + cnt[val]] = i;
        gameState.y[val * 100 + cnt[val]] = j;
        gameState.a[i][j] = val * 100 + cnt[val];
        cnt[val]++;
      }
      else {
        gameState.a[i][j] = INF;
      }
    }
  }

  // gameState.cellUse
  rep(i, n)
  {
    rep(j, n)
    {
      gameState.cellUse[i][j] = 0;
      if (gameState.a[i][j] != INF) {
        gameState.cellUse[i][j] = 1;
      }
    }
  }

  // R,D
  rep(i, n) { gameState.updateRowR(i); }
  rep(j, n) { gameState.updateColD(j); }

  // gameState.udlr
  rep(i, K100)
  {
    rep(j, 4) { gameState.udlr[i][j] = -1; }
  }

  if (gameState.viewOrder == 0) {
    // 横縦の順
    rep(i, K100)
    {
      if (gameState.R[i] == -1) {
        continue;
      }
      if (i / 100 == gameState.R[i] / 100) {
        gameState.udlr[i][3] = gameState.R[i];
        gameState.udlr[gameState.R[i]][2] = i;
        int aVal = MakeAValue(i, gameState.R[i]);
        srep(k, gameState.y[i] + 1, gameState.y[gameState.R[i]]) { gameState.a[gameState.x[i]][k] = aVal; }
      }
    }

    rep(i, K100)
    {
      if (gameState.D[i] == -1) {
        continue;
      }
      if (i / 100 == gameState.D[i] / 100) {
        int ok = 1;
        srep(k, gameState.x[i] + 1, gameState.x[gameState.D[i]])
        {
          if (gameState.a[k][gameState.y[i]] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          gameState.udlr[i][1] = gameState.D[i];
          gameState.udlr[gameState.D[i]][0] = i;
          int aVal = MakeAValue(i, gameState.D[i]);
          srep(k, gameState.x[i] + 1, gameState.x[gameState.D[i]]) { gameState.a[k][gameState.y[i]] = aVal; }
        }
      }
    }
  }
  else {
    // 縦横の順
    rep(i, K100)
    {
      if (gameState.D[i] == -1) {
        continue;
      }
      if (i / 100 == gameState.D[i] / 100) {
        gameState.udlr[i][1] = gameState.D[i];
        gameState.udlr[gameState.D[i]][0] = i;
        int aVal = MakeAValue(i, gameState.D[i]);
        srep(k, gameState.x[i] + 1, gameState.x[gameState.D[i]]) { gameState.a[k][gameState.y[i]] = aVal; }
      }
    }

    rep(i, K100)
    {
      if (gameState.R[i] == -1) {
        continue;
      }
      if (i / 100 == gameState.R[i] / 100) {
        int ok = 1;
        srep(k, gameState.y[i] + 1, gameState.y[gameState.R[i]])
        {
          if (gameState.a[gameState.x[i]][k] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          gameState.udlr[i][3] = gameState.R[i];
          gameState.udlr[gameState.R[i]][2] = i;
          int aVal = MakeAValue(i, gameState.R[i]);
          srep(k, gameState.y[i] + 1, gameState.y[gameState.R[i]]) { gameState.a[gameState.x[i]][k] = aVal; }
        }
      }
    }
  }

  gameState.ope1 = 0;
  gameState.ope2 = 0;
  gameState.maxScore = 0;
  rep(i, K * 100) { gameState.moveCnt[i] = 0; }

  // gameState.parent, gameState.vp
  gameState.vp.clear();
  visitedCnt++;
  rep(i, K100) { visited[i] = -1; }

  rep(i, K100) { gameState.unionSize[i] = 0; }

  rep(i, K100)
  {
    if (visited[i] == visitedCnt) continue;

    int queL = 0;
    int queR = 0;

    que[queR] = i;
    queR++;
    visited[i] = visitedCnt;
    gameState.parent[i] = i;

    int countSize = 1;
    while (queL < queR) {
      int ite = que[queL];
      queL++;
      rep(j, 4)
      {
        int nxt = gameState.udlr[ite][j];
        if (nxt != -1 && visited[nxt] != visitedCnt) {
          que[queR] = nxt;
          queR++;
          visited[nxt] = visitedCnt;
          gameState.parent[nxt] = i;
          countSize++;
        }
      }
    }

    if (countSize >= 2) {
      gameState.vp.insert(P(-countSize, i));
      gameState.unionSize[i] = countSize;
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
    cout << gameState.ope1 << endl;
    rep(i, gameState.ope1)
    {
      rep(j, 4) { cout << gameState.ans1[i][j] << ' '; }
      cout << endl;
    }

    cout << gameState.ope2 << endl;
    rep(i, gameState.ope2)
    {
      rep(j, 4) { cout << gameState.ans2[i][j] << ' '; }
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

    ofs << gameState.ope1 << endl;
    rep(i, gameState.ope1)
    {
      rep(j, 4) { ofs << gameState.ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << gameState.ope2 << endl;
    rep(i, gameState.ope2)
    {
      rep(j, 4) { ofs << gameState.ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs << endl;

    rep(i, gameState.ope1 + 1)
    {
      ofs << i << endl;
      rep(j, i)
      {
        rep(k, 4) { ofs << gameState.ans1[j][k] << ' '; }
        ofs << endl;
      }
      ofs << 0 << endl;
      ofs << endl;
    }

    ofs << gameState.ope1 << endl;
    rep(i, gameState.ope1)
    {
      rep(j, 4) { ofs << gameState.ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << gameState.ope2 << endl;
    rep(i, gameState.ope2)
    {
      rep(j, 4) { ofs << gameState.ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs.close();
  }
}


// maxとreal_maxを初期化
void AllClear()
{
  gameState.clear();
  real_GameState.clear();
  MethodCountReset();
}

void AllClear_seed()
{
  gameState.clear();
  real_GameState.clear();
  seed_GameState.clear();
  MethodCountReset();
}

void AllClear_outer()
{
  gameState.clear();
  real_GameState.clear();
  seed_GameState.clear();
  outer_GameState.clear();
  MethodCountReset();
}







// コンピュータをランダムに1マス移動
void PushACnt(int xx, int yy)
{
  keepA[acnt][0] = xx;
  keepA[acnt][1] = yy;
  keepA[acnt][2] = gameState.a[xx][yy];
  acnt++;
}
void BackA()
{
  drep(i, acnt) { gameState.a[keepA[i][0]][keepA[i][1]] = keepA[i][2]; }
  acnt = 0;
}

void Push_udlr(int ite, int dir)
{
  keep_udlr[udlrcnt][0] = ite;
  keep_udlr[udlrcnt][1] = dir;
  keep_udlr[udlrcnt][2] = gameState.udlr[ite][dir];
  udlrcnt++;
}
void Update_udlr(int ite, int dir, int val)
{
  Push_udlr(ite, dir);
  gameState.udlr[ite][dir] = val;
}
void Back_udlr()
{
  drep(i, udlrcnt) { gameState.udlr[keep_udlr[i][0]][keep_udlr[i][1]] = keep_udlr[i][2]; }
  udlrcnt = 0;
}

void PushParent(int ite)
{
  keep_parent[parentcnt][0] = ite;
  keep_parent[parentcnt][1] = gameState.parent[ite];
  parentcnt++;
}
void BackParent()
{
  drep(i, parentcnt) { gameState.parent[keep_parent[i][0]] = keep_parent[i][1]; }
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
      gameState.vp.erase(P(keep_vp[i][0], keep_vp[i][1]));
      gameState.unionSize[keep_vp[i][1]] = 0;
    }
    else {
      gameState.vp.insert(P(keep_vp[i][0], keep_vp[i][1]));
      gameState.unionSize[keep_vp[i][1]] = -keep_vp[i][0];
    }
  }
  vpcnt = 0;
}

inline void EraseUnion(int ite)
{
  if (gameState.unionSize[gameState.parent[ite]] != 0) {
    gameState.vp.erase(P(-gameState.unionSize[gameState.parent[ite]], gameState.parent[ite]));
    PushVp(-gameState.unionSize[gameState.parent[ite]], gameState.parent[ite], 1);
    gameState.unionSize[gameState.parent[ite]] = 0;
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

  int nx = gameState.x[ite] + dx[dir];
  int ny = gameState.y[ite] + dy[dir];

  int xx = gameState.x[ite];
  int yy = gameState.y[ite];

  int na = gameState.a[nx][ny];

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
      if (gameState.udlr[ite][2] != -1) {
        if (HasServer(xx, yy - 1)) {
          gameState.a[xx][yy] = MakeAValue(ite, gameState.a[xx][yy - 1]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx][yy - 1];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの右のつながり
      gameState.a[xx][ny] = ite;

      // 元のマスの上下のつながり
      // 左と繋がっている場合、切る
      if (gameState.udlr[ite][2] != -1) {
        if (gameState.udlr[ite][0] != -1) {
          int iteU = gameState.udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }

        if (gameState.udlr[ite][1] != -1) {
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
      }
      else {
        // 上下繋がっていた
        if (gameState.udlr[ite][0] != -1 && gameState.udlr[ite][1] != -1) {
          // 繋ぎなおす
          int iteU = gameState.udlr[ite][0];
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteU, 1, iteD);
          Update_udlr(iteD, 0, iteU);
          Update_udlr(ite, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          int aVal = MakeAValue(iteU, iteD);
          srep(i, gameState.x[iteU] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = aVal;
          }
        }
        // 上繋がっていた
        else if (gameState.udlr[ite][0] != -1) {
          // 切る
          int iteU = gameState.udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 下繋がっていた
        else if (gameState.udlr[ite][1] != -1) {
          // 切る
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
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
            srep(i, gameState.x[iteU] + 1, gameState.x[iteD])
            {
              PushACnt(i, yy);
              gameState.a[i][yy] = aVal;
            }
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[ite][3] == -1) {
        int iteU = -na / 1000;
        int iteD = -na % 1000;
        if (gameState.x[iteU] > gameState.x[iteD]) {
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
          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
          }

          aVal = MakeAValue(ite, iteD);
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
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

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = INF;
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
          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
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
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
          }
        }
      }
    }
    // 左
    else {
      // 元のマスの右のつながり
      // 繋がっている
      if (gameState.udlr[ite][3] != -1) {
        if (HasServer(xx, yy + 1)) {
          gameState.a[xx][yy] = MakeAValue(ite, gameState.a[xx][yy + 1]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx][yy + 1];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの左のつながり
      gameState.a[xx][ny] = ite;

      // 元のマスの上下のつながり
      // 右と繋がっている場合、切る
      if (gameState.udlr[ite][3] != -1) {
        if (gameState.udlr[ite][0] != -1) {
          int iteU = gameState.udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }

        if (gameState.udlr[ite][1] != -1) {
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
      }
      else {
        // 上下繋がっていた
        if (gameState.udlr[ite][0] != -1 && gameState.udlr[ite][1] != -1) {
          // 繋ぎなおす
          int iteU = gameState.udlr[ite][0];
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteU, 1, iteD);
          Update_udlr(iteD, 0, iteU);
          Update_udlr(ite, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          int aVal = MakeAValue(iteU, iteD);
          srep(i, gameState.x[iteU] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = aVal;
          }
        }
        // 上繋がっていた
        else if (gameState.udlr[ite][0] != -1) {
          // 切る
          int iteU = gameState.udlr[ite][0];
          Update_udlr(iteU, 1, -1);
          Update_udlr(ite, 0, -1);

          EraseUnion(iteU);
          EraseUnion(ite);
          se.insert(iteU);
          se.insert(ite);

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
            beam.push_back(P(i, yy));
          }
        }
        // 下繋がっていた
        else if (gameState.udlr[ite][1] != -1) {
          // 切る
          int iteD = gameState.udlr[ite][1];
          Update_udlr(iteD, 0, -1);
          Update_udlr(ite, 1, -1);

          EraseUnion(iteD);
          EraseUnion(ite);
          se.insert(iteD);
          se.insert(ite);

          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, yy);
            gameState.a[i][yy] = INF;
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
            srep(i, gameState.x[iteU] + 1, gameState.x[iteD])
            {
              PushACnt(i, yy);
              gameState.a[i][yy] = aVal;
            }
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[ite][2] == -1) {
        int iteU = -na / 1000;
        int iteD = -na % 1000;
        if (gameState.x[iteU] > gameState.x[iteD]) {
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
          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
          }

          aVal = MakeAValue(ite, iteD);
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
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

          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = INF;
            beam.push_back(P(i, ny));
          }
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = INF;
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
          srep(i, gameState.x[iteU] + 1, gameState.x[ite])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
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
          srep(i, gameState.x[ite] + 1, gameState.x[iteD])
          {
            PushACnt(i, ny);
            gameState.a[i][ny] = aVal;
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
      if (gameState.udlr[ite][0] != -1) {
        if (HasServer(xx - 1, yy)) {
          gameState.a[xx][yy] = MakeAValue(ite, gameState.a[xx - 1][yy]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx - 1][yy];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの下のつながり
      gameState.a[nx][yy] = ite;

      // 元のマスの左右のつながり
      // 上と繋がっている場合、切る
      if (gameState.udlr[ite][0] != -1) {
        if (gameState.udlr[ite][2] != -1) {
          int iteL = gameState.udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }

        if (gameState.udlr[ite][3] != -1) {
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
      }
      else {
        // 左右繋がっていた
        if (gameState.udlr[ite][2] != -1 && gameState.udlr[ite][3] != -1) {
          // 繋ぎなおす
          int iteL = gameState.udlr[ite][2];
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteL, 3, iteR);
          Update_udlr(iteR, 2, iteL);
          Update_udlr(ite, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          int aVal = MakeAValue(iteL, iteR);
          srep(i, gameState.y[iteL] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = aVal;
          }
        }
        // 左繋がっていた
        else if (gameState.udlr[ite][2] != -1) {
          // 切る
          int iteL = gameState.udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 右繋がっていた
        else if (gameState.udlr[ite][3] != -1) {
          // 切る
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(iteR);
          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
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
            srep(i, gameState.y[iteL] + 1, gameState.y[iteR])
            {
              PushACnt(xx, i);
              gameState.a[xx][i] = aVal;
            }
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[ite][1] == -1) {
        int iteL = -na / 1000;
        int iteR = -na % 1000;
        if (gameState.y[iteL] > gameState.y[iteR]) {
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
          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
          }

          aVal = MakeAValue(ite, iteR);
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
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

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = INF;
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
          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
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
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
          }
        }
      }
    }
    // 上
    else {
      // 元のマスの下のつながり
      // 繋がっている
      if (gameState.udlr[ite][1] != -1) {
        if (HasServer(xx + 1, yy)) {
          gameState.a[xx][yy] = MakeAValue(ite, gameState.a[xx + 1][yy]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx + 1][yy];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの上のつながり
      gameState.a[nx][yy] = ite;

      // 元のマスの左右のつながり
      // 下と繋がっている場合、切る
      if (gameState.udlr[ite][1] != -1) {
        if (gameState.udlr[ite][2] != -1) {
          int iteL = gameState.udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }

        if (gameState.udlr[ite][3] != -1) {
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
      }
      else {
        // 左右繋がっていた
        if (gameState.udlr[ite][2] != -1 && gameState.udlr[ite][3] != -1) {
          // 繋ぎなおす
          int iteL = gameState.udlr[ite][2];
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteL, 3, iteR);
          Update_udlr(iteR, 2, iteL);
          Update_udlr(ite, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          int aVal = MakeAValue(iteL, iteR);
          srep(i, gameState.y[iteL] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = aVal;
          }
        }
        // 左繋がっていた
        else if (gameState.udlr[ite][2] != -1) {
          // 切る
          int iteL = gameState.udlr[ite][2];
          Update_udlr(iteL, 3, -1);
          Update_udlr(ite, 2, -1);

          EraseUnion(iteL);
          EraseUnion(ite);
          se.insert(iteL);
          se.insert(ite);

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
            beam.push_back(P(xx, i));
          }
        }
        // 右繋がっていた
        else if (gameState.udlr[ite][3] != -1) {
          // 切る
          int iteR = gameState.udlr[ite][3];
          Update_udlr(iteR, 2, -1);
          Update_udlr(ite, 3, -1);

          EraseUnion(iteR);
          EraseUnion(ite);
          se.insert(iteR);
          se.insert(ite);

          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(xx, i);
            gameState.a[xx][i] = INF;
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
            srep(i, gameState.y[iteL] + 1, gameState.y[iteR])
            {
              PushACnt(xx, i);
              gameState.a[xx][i] = aVal;
            }
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[ite][0] == -1) {
        int iteL = -na / 1000;
        int iteR = -na % 1000;
        if (gameState.y[iteL] > gameState.y[iteR]) {
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
          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
          }

          aVal = MakeAValue(ite, iteR);
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
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

          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = INF;
            beam.push_back(P(nx, i));
          }
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = INF;
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
          srep(i, gameState.y[iteL] + 1, gameState.y[ite])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
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
          srep(i, gameState.y[ite] + 1, gameState.y[iteR])
          {
            PushACnt(nx, i);
            gameState.a[nx][i] = aVal;
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
      if (gameState.a[px][py] != INF) continue;
      int iteL = GetIte2(px, py, 'L');
      int iteR = GetIte2(px, py, 'R');
      if (iteL != -1 && iteR != -1 && iteL / 100 == iteR / 100) {
        Update_udlr(iteL, 3, iteR);
        Update_udlr(iteR, 2, iteL);

        EraseUnion(iteL);
        EraseUnion(iteR);
        se.insert(iteL);

        int aVal = MakeAValue(iteL, iteR);
        srep(i, gameState.y[iteL] + 1, gameState.y[iteR])
        {
          PushACnt(px, i);
          gameState.a[px][i] = aVal;
        }
      }
    }
  }
  // 上下
  else {
    for (auto&& p : beam) {
      int px = p.first;
      int py = p.second;
      if (gameState.a[px][py] != INF) continue;
      int iteU = GetIte2(px, py, 'U');
      int iteD = GetIte2(px, py, 'D');
      if (iteU != -1 && iteD != -1 && iteU / 100 == iteD / 100) {
        Update_udlr(iteU, 1, iteD);
        Update_udlr(iteD, 0, iteU);

        EraseUnion(iteU);
        EraseUnion(iteD);
        se.insert(iteU);

        int aVal = MakeAValue(iteU, iteD);
        srep(i, gameState.x[iteU] + 1, gameState.x[iteD])
        {
          PushACnt(i, py);
          gameState.a[i][py] = aVal;
        }
      }
    }
  }

  gameState.x[ite] = nx;
  gameState.y[ite] = ny;



  if (nx != xx) {
    gameState.updateRowR(xx);
    gameState.updateRowR(nx);
  }


  if (ny != yy) {
    gameState.updateColD(yy);
    gameState.updateColD(ny);
  }


  // gameState.parent, vpの更新
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
    gameState.parent[vpIte] = vpIte;

    int countSize = 1;
    while (queL < queR) {
      int nowIte = que[queL];
      queL++;
      rep(j, 4)
      {
        int nxt = gameState.udlr[nowIte][j];
        if (nxt != -1 && visited[nxt] != visitedCnt) {
          que[queR] = nxt;
          queR++;
          visited[nxt] = visitedCnt;
          PushParent(nxt);
          gameState.parent[nxt] = vpIte;
          countSize++;
        }
      }
    }

    if (countSize >= 2) {
      gameState.vp.insert(P(-countSize, vpIte));
      PushVp(-countSize, vpIte, 0);
      gameState.unionSize[vpIte] = countSize;
    }
  }

  int tmpScore = 0;
  if (MethodeMode == 5) {
    tmpScore = gameState.calcScore(K100 - (gameState.ope1 - 1));
  }
  else {
    tmpScore = gameState.calcScore(K100 - (gameState.ope1 + 1));
  }

  methodCount[1][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - gameState.maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  int isDo = 0;
  if (forceDo || prob > Rand01()) {
    isDo = 1;
    gameState.maxScore += diffScore;

    methodCount[1][0]++;
    methodSum[0]++;

    gameState.ans1[gameState.ope1][0] = xx;
    gameState.ans1[gameState.ope1][1] = yy;
    gameState.ans1[gameState.ope1][2] = nx;
    gameState.ans1[gameState.ope1][3] = ny;
    gameState.ans1[gameState.ope1][4] = ite;
    gameState.ope1++;

    if (gameState.maxScore > real_GameState.maxScore) {
      if (MethodeMode == 5) {
        isDo = 5;
      }
      gameState.copyTo(real_GameState);
    }
  }
  else {
    // 元に戻す
    gameState.x[ite] = xx;
    gameState.y[ite] = yy;

    BackA();
    Back_udlr();
    BackParent();
    BackVp();

    if (nx != xx) {
      gameState.updateRowR(xx);
      gameState.updateRowR(nx);
    }


    if (ny != yy) {
      gameState.updateColD(yy);
      gameState.updateColD(ny);
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

    nx = gameState.x[ite] + dx[dir];
    ny = gameState.y[ite] + dy[dir];
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
    if (gameState.a[xx][yy] == -1) {
      break;
    }
  }

  dir1 = Rand() % 4;
  dir2 = Rand() % 4;
  int nx1 = xx + dx[dir1];
  int ny1 = yy + dy[dir1];
  int nx2 = nx1 + dx[dir2];
  int ny2 = ny1 + dy[dir2];
  if (nx1 < 0 || n <= nx1 || ny1 < 0 || n <= ny1 || gameState.a[nx1][ny1] == -1) {
    return;
  }
  if (nx2 < 0 || n <= nx2 || ny2 < 0 || n <= ny2 || gameState.a[nx2][ny2] == -1) {
    return;
  }

  int ite1 = gameState.a[nx1][ny1];
  std::swap(gameState.a[nx1][ny1], gameState.a[xx][yy]);
  gameState.x[ite1] = xx;
  gameState.y[ite1] = yy;

  int ite2 = gameState.a[nx2][ny2];
  std::swap(gameState.a[nx2][ny2], gameState.a[nx1][ny1]);
  gameState.x[ite2] = nx1;
  gameState.y[ite2] = ny1;

  gameState.updateRowR(xx);
  gameState.updateRowR(nx1);
  gameState.updateRowR(nx2);
  gameState.updateColD(yy);
  gameState.updateColD(ny1);
  gameState.updateColD(ny2);

  int tmpScore = gameState.calcScore(K100 - (gameState.ope1 + 2));
  methodCount[3][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - gameState.maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    gameState.maxScore += diffScore;

    methodCount[3][0]++;
    methodSum[0]++;

    gameState.ans1[gameState.ope1][0] = nx1;
    gameState.ans1[gameState.ope1][1] = ny1;
    gameState.ans1[gameState.ope1][2] = xx;
    gameState.ans1[gameState.ope1][3] = yy;
    gameState.ans1[gameState.ope1][4] = ite1;
    gameState.ope1++;

    gameState.ans1[gameState.ope1][0] = nx2;
    gameState.ans1[gameState.ope1][1] = ny2;
    gameState.ans1[gameState.ope1][2] = nx1;
    gameState.ans1[gameState.ope1][3] = ny1;
    gameState.ans1[gameState.ope1][4] = ite2;
    gameState.ope1++;

    if (gameState.maxScore > real_GameState.maxScore) {
      gameState.copyTo(real_GameState);
    }
  }
  else {
    // 元に戻す
    gameState.x[ite2] = nx2;
    gameState.y[ite2] = ny2;
    std::swap(gameState.a[nx2][ny2], gameState.a[nx1][ny1]);

    gameState.x[ite1] = nx1;
    gameState.y[ite1] = ny1;
    std::swap(gameState.a[nx1][ny1], gameState.a[xx][yy]);

    gameState.updateRowR(xx);
    gameState.updateRowR(nx1);
    gameState.updateRowR(nx2);
    gameState.updateColD(yy);
    gameState.updateColD(ny1);
    gameState.updateColD(ny2);
  }
}

// コンピュータをランダムに2マス移動
void Method4(double start_temp, double end_temp, double now_progress)
{
  int ite = Rand() % K100;  // 1マス動かすコンピュータ
  int dir1 = Rand() % 4;
  int dir2 = Rand() % 4;

  int xx = gameState.x[ite];
  int yy = gameState.y[ite];

  int nx1 = xx + dx[dir1];
  int ny1 = yy + dy[dir1];
  int nx2 = nx1 + dx[dir2];
  int ny2 = ny1 + dy[dir2];

  if (IsNG(nx1, ny1) || IsNG(nx2, ny2) || gameState.a[nx1][ny1] != -1 ||
    gameState.a[nx2][ny2] != -1) {
    return;
  }

  swap(gameState.a[xx][yy], gameState.a[nx2][ny2]);
  gameState.x[ite] = nx2;
  gameState.y[ite] = ny2;

  gameState.updateRowR(xx);
  gameState.updateRowR(nx1);
  gameState.updateRowR(nx2);
  gameState.updateColD(yy);
  gameState.updateColD(ny1);
  gameState.updateColD(ny2);

  int tmpScore = gameState.calcScore(K100 - (gameState.ope1 + 2));
  methodCount[4][1]++;
  methodSum[1]++;

  int diffScore = tmpScore - gameState.maxScore;

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    gameState.maxScore += diffScore;

    methodCount[4][0]++;
    methodSum[0]++;

    gameState.ans1[gameState.ope1][0] = xx;
    gameState.ans1[gameState.ope1][1] = yy;
    gameState.ans1[gameState.ope1][2] = nx1;
    gameState.ans1[gameState.ope1][3] = ny1;
    gameState.ans1[gameState.ope1][4] = ite;
    gameState.ope1++;

    gameState.ans1[gameState.ope1][0] = nx1;
    gameState.ans1[gameState.ope1][1] = ny1;
    gameState.ans1[gameState.ope1][2] = nx2;
    gameState.ans1[gameState.ope1][3] = ny2;
    gameState.ans1[gameState.ope1][4] = ite;
    gameState.ope1++;

    if (gameState.maxScore > real_GameState.maxScore) {
      gameState.copyTo(real_GameState);
    }
  }
  else {
    // 元に戻す
    swap(gameState.a[xx][yy], gameState.a[nx2][ny2]);
    gameState.x[ite] = xx;
    gameState.y[ite] = yy;

    gameState.updateRowR(xx);
    gameState.updateRowR(nx1);
    gameState.updateRowR(nx2);
    gameState.updateColD(yy);
    gameState.updateColD(ny1);
    gameState.updateColD(ny2);
  }
}

// 移動をランダムに1つ削除
void Method5(double start_temp, double end_temp, double now_progress)
{
  if (gameState.ope1 == 0) return;
  int ite = Rand() % gameState.ope1;

  // NGチェック
  // ite以降の操作で、操作元が移動後のマス、操作後が移動前のマス、の操作が出てこなければOK
  srep(i, ite + 1, gameState.ope1)
  {
    if (gameState.ans1[i][4] == gameState.ans1[ite][4]) {
      return;
    }
    if (gameState.ans1[i][0] == gameState.ans1[ite][2] && gameState.ans1[i][1] == gameState.ans1[ite][3]) {
      return;
    }
    if (gameState.ans1[i][2] == gameState.ans1[ite][0] && gameState.ans1[i][3] == gameState.ans1[ite][1]) {
      return;
    }
  }

  // 消したい移動の逆の操作を足してInnerMethodを実行
  int reverseDir = 0;
  rep(i, 4)
  {
    if (gameState.ans1[ite][2] == gameState.ans1[ite][0] + dx[i] &&
      gameState.ans1[ite][3] == gameState.ans1[ite][1] + dy[i]) {
      reverseDir = (i + 2) % 4;
      break;
    }
  }

  int isDo = InnerMethod(start_temp, end_temp, now_progress, gameState.ans1[ite][4],
    reverseDir, false, 5);

  // 実行した場合、2つ消す
  if (isDo) {
    /*
    修正するもの
    int ope1;
    int gameState.ans1[1000][5];
    */
    // 消すのはiteとope1-1
    gameState.ope1--;
    srep(i, ite, gameState.ope1 - 1)
    {
      rep(j, 5) { gameState.ans1[i][j] = gameState.ans1[i + 1][j]; }
    }
    gameState.ope1--;
    if (isDo == 5) {
      gameState.copyTo(real_GameState);
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
    rep(i, gameState.ope1 - 1)
    {
      if (gameState.ans1[i][0] == gameState.ans1[i + 1][2] && gameState.ans1[i][1] == gameState.ans1[i + 1][3] &&
        gameState.ans1[i][2] == gameState.ans1[i + 1][0] && gameState.ans1[i][3] == gameState.ans1[i + 1][1]) {
        ite = i;
        break;
      }
    }

    if (ite == -1) {
      return;
    }

    gameState.maxScore = gameState.calcScore(K100 - (gameState.ope1 - 2));

    methodCount[6][0]++;
    methodSum[0]++;

    srep(i, ite, gameState.ope1 - 2)
    {
      rep(j, 5) { gameState.ans1[i][j] = gameState.ans1[i + 2][j]; }
    }

    gameState.ope1 -= 2;

    if (gameState.maxScore > real_GameState.maxScore) {
      gameState.copyTo(real_GameState);
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
    rep(i, gameState.ope1)
    {
      int ite = gameState.ans1[i][4];
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
        if (gameState.ans1[ite1][0] == gameState.ans1[ite2][2] && gameState.ans1[ite1][1] == gameState.ans1[ite2][3] &&
          gameState.ans1[ite1][2] == gameState.ans1[ite2][0] && gameState.ans1[ite1][3] == gameState.ans1[ite2][1]) {
          // ngチェック
          int ng = 0;
          srep(k, ite1 + 1, ite2)
          {
            if (gameState.ans1[k][2] == gameState.ans1[ite1][0] && gameState.ans1[k][3] == gameState.ans1[ite1][1]) {
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
      srep(i, it2, gameState.ope1 - 1)
      {
        rep(j, 5)
        {
          gameState.ans1[i][j] = gameState.ans1[i + 1][j];
        }
      }
      gameState.ope1--;
      srep(i, it1, gameState.ope1 - 1)
      {
        rep(j, 5)
        {
          gameState.ans1[i][j] = gameState.ans1[i + 1][j];
        }
      }
      gameState.ope1--;

      gameState.maxScore = gameState.calcScore(K100 - gameState.ope1);
      methodCount[7][0]++;
      methodSum[0]++;

      if (gameState.maxScore > real_GameState.maxScore) {
        gameState.copyTo(real_GameState);
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
  gameState.maxScore = gameState.calcScore(K100, true);
  gameState.copyTo(real_GameState);
  gameState.copyTo(seed_GameState);

  // シード作り
  int seedCount = 10;
  rep(tei, seedCount)
  {
    start_time = clock();

    Init();
    gameState.viewOrder = tei % 2;
    gameState.maxScore = gameState.calcScore(K100, true);

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
      if (gameState.maxScore * 1.2 < real_GameState.maxScore || Rand() % 123456 == 0) {
        gameState.copyFrom(real_GameState);
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
        if (gameState.ope1 >= K100) {
          continue;
        }
        Method1(start_temp, end_temp, now_progress);
      }

      // 空白を2マス動かす
      if (me == 3) {
        if (gameState.ope1 >= K100 - 1) {
          continue;
        }
        Method3(start_temp, end_temp, now_progress);
      }

      // コンピュータをランダムに2マス移動
      if (me == 4) {
        if (gameState.ope1 >= K100 - 1) {
          continue;
        }
        Method4(start_temp, end_temp, now_progress);
      }

      // 移動をランダムに1つ削除
      if (me == 5) {
        if (gameState.ope1 == 0) {
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
    gameState.copyFrom(real_GameState);
    if (gameState.maxScore > seed_GameState.maxScore) {
      gameState.copyTo(seed_GameState);
    }

    AllClear();
  }

  // シードから戻す
  gameState.copyFrom(seed_GameState);
  gameState.copyTo(real_GameState);

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
    if (gameState.maxScore * 1.2 < real_GameState.maxScore || Rand() % 123456 == 0) {
      gameState.copyFrom(real_GameState);
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
      if (gameState.ope1 >= K100) {
        continue;
      }
      Method1(start_temp, end_temp, now_progress);
    }

    // 空白を2マス動かす
    if (me == 3) {
      if (gameState.ope1 >= K100 - 1) {
        continue;
      }
      Method3(start_temp, end_temp, now_progress);
    }

    // コンピュータをランダムに2マス移動
    if (me == 4) {
      if (gameState.ope1 >= K100 - 1) {
        continue;
      }
      Method4(start_temp, end_temp, now_progress);
    }

    // 移動をランダムに1つ削除
    if (me == 5) {
      if (gameState.ope1 == 0) {
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

  gameState.copyFrom(real_GameState);

  int cal = gameState.calcScore(K100 - gameState.ope1, true);

  if (true) {
    cerr << "problemNum = " << problemNum << ", N = " << n << ", K = " << K << endl;
    cerr << "start_temp = " << start_temp << ", gameState.viewOrder = " << gameState.viewOrder << endl;
    cerr << "gameState.maxScore = " << gameState.maxScore << ", gameState.ope1 = " << gameState.ope1 << ", gameState.ope2 = " << gameState.ope2 << endl;
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

  return gameState.maxScore;
}

int SolveOuter(int mode, int problemNum)
{
  // 入力部
  Input(problemNum);

  rep(_, outer_Split)
  {
    gameState.viewOrder = _ % 2;
    int score = Solve(mode, problemNum);
    if (score >= outer_GameState.maxScore) {
      gameState.copyTo(outer_GameState);
    }
    AllClear();
  }
  gameState.copyFrom(outer_GameState);

  // 解の出力
  Output(mode, problemNum);

  return gameState.maxScore;
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
