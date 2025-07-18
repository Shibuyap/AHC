﻿#include <algorithm>
#include <algorithm>
#include <array>
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
typedef pair<int, int> P;
#define MAX_N 200005
#define INF 1001001001

// Constants
const int SERVERS_PER_TYPE = 100;
const int ENCODING_MULTIPLIER = 1000;
const int MAX_COMPUTERS = 500;
const int MAX_OPERATIONS = 500;
const int MAX_ITERATIONS = 10000;

// Direction constants
enum Direction { UP = 0, LEFT = 1, DOWN = 2, RIGHT = 3 };
const array<int, 4> dx = { -1, 0, 1, 0 };
const array<int, 4> dy = { 0, -1, 0, 1 };
const array<char, 4> DIR_CHARS = { 'U', 'L', 'D', 'R' };

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

namespace /* 乱数ライブラリ */
{
  static uint32_t rand32()
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
    return (rand32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

// Helper functions for server ID encoding/decoding
inline int encodeServerId(int type, int index)
{
  return type * SERVERS_PER_TYPE + index;
}

inline int getServerType(int serverId)
{
  return serverId / SERVERS_PER_TYPE;
}

inline int getServerIndex(int serverId)
{
  return serverId % SERVERS_PER_TYPE;
}

// Helper function for bounds checking
inline bool isInBounds(int x, int y, int n)
{
  return x >= 0 && x < n && y >= 0 && y < n;
}

// Helper function to check if two servers are same type
inline bool areSameType(int server1, int server2)
{
  return server1 != -1 && server2 != -1 &&
    getServerType(server1) == getServerType(server2);
}

// Helper function to convert char to Direction
Direction charToDir(char c)
{
  switch (c) {
    case 'U': return UP;
    case 'D': return DOWN;
    case 'L': return LEFT;
    case 'R': return RIGHT;
  }
  return UP; // Default
}

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
  void initialize(int maxComputers, int maxOperations);
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

inline bool is_out_of_bounds(int xx, int yy)
{
  return !isInBounds(xx, yy, n);
}

inline bool has_server(int xx, int yy)
{
  if (0 <= gameState.a[xx][yy] && gameState.a[xx][yy] < K100) {
    return true;
  }
  return false;
}

inline int make_edge_value(int server1, int server2)
{
  int num = -1 * (server1 * ENCODING_MULTIPLIER + server2);
  if (server2 < server1) {
    num = -1 * (server2 * ENCODING_MULTIPLIER + server1);
  }
  return num;
}

// Helper function to find server in a direction
int find_server_in_direction(int xx, int yy, int deltaX, int deltaY, bool checkINF = false)
{
  xx += deltaX;
  yy += deltaY;

  while (isInBounds(xx, yy, n) && !has_server(xx, yy)) {
    if (checkINF && gameState.a[xx][yy] != INF) {
      return -1;
    }
    xx += deltaX;
    yy += deltaY;
  }

  if (!isInBounds(xx, yy, n)) {
    return -1;
  }

  return gameState.a[xx][yy];
}

int get_server(int xx, int yy, char cc)
{
  Direction dir = charToDir(cc);
  return find_server_in_direction(xx, yy, dx[dir], dy[dir], false);
}

// サーバーまでたどり着けるか
int get_server_through_inf(int xx, int yy, char cc)
{
  Direction dir = charToDir(cc);
  return find_server_in_direction(xx, yy, dx[dir], dy[dir], true);
}

void reset_method_count()
{
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 2; ++j) { methodCount[i][j] = 0; }
  }
  methodSum[0] = 0;
  methodSum[1] = 0;
}


// Forward declaration
bool has_server(int xx, int yy);

// GameState member function implementations
void GameState::clear()
{
  viewOrder = 0;
  maxScore = 0;
  ope1 = 0;
  ope2 = 0;
}

void GameState::initialize(int maxComputers, int maxOperations)
{
  x.resize(maxComputers);
  y.resize(maxComputers);
  R.resize(maxComputers);
  D.resize(maxComputers);
  parent.resize(maxComputers);
  unionSize.resize(maxComputers);
  moveCnt.resize(maxComputers);
  ans1.resize(maxOperations, vector<int>(5));
  ans2.resize(maxOperations, vector<int>(4));
  moves.resize(maxComputers, vector<int>(maxOperations));
  udlr.resize(maxComputers, vector<int>(4));
}

void GameState::copyTo(GameState& target) const
{
  target.viewOrder = viewOrder;
  target.maxScore = maxScore;
  target.ope1 = ope1;
  target.ope2 = ope2;

  for (int i = 0; i < ope1; ++i) {
    for (int j = 0; j < 5; ++j) { target.ans1[i][j] = ans1[i][j]; }
  }
  for (int i = 0; i < ope2; ++i) {
    for (int j = 0; j < 4; ++j) { target.ans2[i][j] = ans2[i][j]; }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { target.a[i][j] = a[i][j]; }
  }
  for (int i = 0; i < K100; ++i) {
    target.x[i] = x[i];
    target.y[i] = y[i];
    target.R[i] = R[i];
    target.D[i] = D[i];
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { target.cellUse[i][j] = cellUse[i][j]; }
  }

  for (int i = 0; i < K100; ++i) {
    target.moveCnt[i] = moveCnt[i];
    for (int j = 0; j < moveCnt[i]; ++j) { target.moves[i][j] = moves[i][j]; }
    for (int j = 0; j < 4; ++j) { target.udlr[i][j] = udlr[i][j]; }
  }

  for (int i = 0; i < K100; ++i) {
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
  while (yy < n && !has_server(xx, yy)) {
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
  while (xx < n && !has_server(xx, yy)) {
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
  for (int j = 0; j < n; ++j) {
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
  for (int i = 0; i < n; ++i) {
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
        for (int j = 0; j < 4; ++j) {
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

void init()
{
  visitedCnt = 0;

  // Initialize GameState objects
  gameState.initialize(MAX_COMPUTERS, MAX_OPERATIONS);
  real_GameState.initialize(MAX_COMPUTERS, MAX_OPERATIONS);
  seed_GameState.initialize(MAX_COMPUTERS, MAX_OPERATIONS);
  outer_GameState.initialize(MAX_COMPUTERS, MAX_OPERATIONS);

  // Initialize global vectors
  visited.resize(MAX_COMPUTERS);
  que.resize(MAX_COMPUTERS);
  keep_vp.resize(MAX_ITERATIONS, vector<int>(3));
  keep_parent.resize(MAX_ITERATIONS, vector<int>(2));
  keep_udlr.resize(MAX_ITERATIONS, vector<int>(3));
  keepA.resize(MAX_ITERATIONS, vector<int>(3));
  vv.resize(MAX_COMPUTERS);

  // a,x,y
  int cnt[5] = {};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int val = s[i][j] - '0' - 1;
      if (val != -1) {
        gameState.x[encodeServerId(val, cnt[val])] = i;
        gameState.y[encodeServerId(val, cnt[val])] = j;
        gameState.a[i][j] = encodeServerId(val, cnt[val]);
        cnt[val]++;
      }
      else {
        gameState.a[i][j] = INF;
      }
    }
  }

  // gameState.cellUse
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      gameState.cellUse[i][j] = 0;
      if (gameState.a[i][j] != INF) {
        gameState.cellUse[i][j] = 1;
      }
    }
  }

  // R,D
  for (int i = 0; i < n; ++i) { gameState.updateRowR(i); }
  for (int j = 0; j < n; ++j) { gameState.updateColD(j); }

  // gameState.udlr
  for (int i = 0; i < K100; ++i) {
    for (int j = 0; j < 4; ++j) { gameState.udlr[i][j] = -1; }
  }

  if (gameState.viewOrder == 0) {
    // 横縦の順
    for (int i = 0; i < K100; ++i) {
      if (gameState.R[i] == -1) {
        continue;
      }
      if (getServerType(i) == gameState.R[i] / 100) {
        gameState.udlr[i][3] = gameState.R[i];
        gameState.udlr[gameState.R[i]][2] = i;
        int aVal = make_edge_value(i, gameState.R[i]);
        for (int k = gameState.y[i] + 1; k < gameState.y[gameState.R[i]]; ++k) { gameState.a[gameState.x[i]][k] = aVal; }
      }
    }

    for (int i = 0; i < K100; ++i) {
      if (gameState.D[i] == -1) {
        continue;
      }
      if (getServerType(i) == gameState.D[i] / 100) {
        int ok = 1;
        for (int k = gameState.x[i] + 1; k < gameState.x[gameState.D[i]]; ++k) {
          if (gameState.a[k][gameState.y[i]] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          gameState.udlr[i][1] = gameState.D[i];
          gameState.udlr[gameState.D[i]][0] = i;
          int aVal = make_edge_value(i, gameState.D[i]);
          for (int k = gameState.x[i] + 1; k < gameState.x[gameState.D[i]]; ++k) { gameState.a[k][gameState.y[i]] = aVal; }
        }
      }
    }
  }
  else {
    // 縦横の順
    for (int i = 0; i < K100; ++i) {
      if (gameState.D[i] == -1) {
        continue;
      }
      if (getServerType(i) == gameState.D[i] / 100) {
        gameState.udlr[i][1] = gameState.D[i];
        gameState.udlr[gameState.D[i]][0] = i;
        int aVal = make_edge_value(i, gameState.D[i]);
        for (int k = gameState.x[i] + 1; k < gameState.x[gameState.D[i]]; ++k) { gameState.a[k][gameState.y[i]] = aVal; }
      }
    }

    for (int i = 0; i < K100; ++i) {
      if (gameState.R[i] == -1) {
        continue;
      }
      if (getServerType(i) == gameState.R[i] / 100) {
        int ok = 1;
        for (int k = gameState.y[i] + 1; k < gameState.y[gameState.R[i]]; ++k) {
          if (gameState.a[gameState.x[i]][k] < 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          gameState.udlr[i][3] = gameState.R[i];
          gameState.udlr[gameState.R[i]][2] = i;
          int aVal = make_edge_value(i, gameState.R[i]);
          for (int k = gameState.y[i] + 1; k < gameState.y[gameState.R[i]]; ++k) { gameState.a[gameState.x[i]][k] = aVal; }
        }
      }
    }
  }

  gameState.ope1 = 0;
  gameState.ope2 = 0;
  gameState.maxScore = 0;
  for (int i = 0; i < K * 100; ++i) { gameState.moveCnt[i] = 0; }

  // gameState.parent, gameState.vp
  gameState.vp.clear();
  visitedCnt++;
  for (int i = 0; i < K100; ++i) { visited[i] = -1; }

  for (int i = 0; i < K100; ++i) { gameState.unionSize[i] = 0; }

  for (int i = 0; i < K100; ++i) {
    if (visited[i] == visitedCnt) { continue; }

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
      for (int j = 0; j < 4; ++j) {
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

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {  // 標準入力する
    cin >> n >> K;
    for (int i = 0; i < n; ++i) cin >> s[i];
  }
  else {  // ファイル入力する
    ifs >> n >> K;
    for (int i = 0; i < n; ++i) ifs >> s[i];
  }

  K100 = K * 100;

  init();
}

void output_data(int mode, int case_num)
{
  if (mode == 0) {
    cout << gameState.ope1 << endl;
    for (int i = 0; i < gameState.ope1; ++i) {
      for (int j = 0; j < 4; ++j) { cout << gameState.ans1[i][j] << ' '; }
      cout << endl;
    }

    cout << gameState.ope2 << endl;
    for (int i = 0; i < gameState.ope2; ++i) {
      for (int j = 0; j < 4; ++j) { cout << gameState.ans2[i][j] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << gameState.ope1 << endl;
    for (int i = 0; i < gameState.ope1; ++i) {
      for (int j = 0; j < 4; ++j) { ofs << gameState.ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << gameState.ope2 << endl;
    for (int i = 0; i < gameState.ope2; ++i) {
      for (int j = 0; j < 4; ++j) { ofs << gameState.ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs << endl;

    for (int i = 0; i < gameState.ope1 + 1; ++i) {
      ofs << i << endl;
      for (int j = 0; j < i; ++j) {
        for (int k = 0; k < 4; ++k) { ofs << gameState.ans1[j][k] << ' '; }
        ofs << endl;
      }
      ofs << 0 << endl;
      ofs << endl;
    }

    ofs << gameState.ope1 << endl;
    for (int i = 0; i < gameState.ope1; ++i) {
      for (int j = 0; j < 4; ++j) { ofs << gameState.ans1[i][j] << ' '; }
      ofs << endl;
    }

    ofs << gameState.ope2 << endl;
    for (int i = 0; i < gameState.ope2; ++i) {
      for (int j = 0; j < 4; ++j) { ofs << gameState.ans2[i][j] << ' '; }
      ofs << endl;
    }
    ofs.close();
  }
}


// maxとreal_maxを初期化
void clear_all()
{
  gameState.clear();
  real_GameState.clear();
  reset_method_count();
}

void clear_seed_all()
{
  gameState.clear();
  real_GameState.clear();
  seed_GameState.clear();
  reset_method_count();
}

void clear_outer_all()
{
  gameState.clear();
  real_GameState.clear();
  seed_GameState.clear();
  outer_GameState.clear();
  reset_method_count();
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
  for (int i = acnt - 1; i >= 0; --i) { gameState.a[keepA[i][0]][keepA[i][1]] = keepA[i][2]; }
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
  for (int i = udlrcnt - 1; i >= 0; --i) { gameState.udlr[keep_udlr[i][0]][keep_udlr[i][1]] = keep_udlr[i][2]; }
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
  for (int i = parentcnt - 1; i >= 0; --i) { gameState.parent[keep_parent[i][0]] = keep_parent[i][1]; }
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
  for (int i = vpcnt - 1; i >= 0; --i) {
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

// サーバー間の接続を切る共通処理
void DisconnectServers(int server1, int server2, int udlr_idx1, int udlr_idx2, set<int>& se)
{
  // 切る
  Update_udlr(server1, udlr_idx1, -1);
  Update_udlr(server2, udlr_idx2, -1);

  EraseUnion(server1);
  se.insert(server2);
}

// サーバー間の接続を切る共通処理（両方のサーバーをEraseUnion/insertするパターン）
void DisconnectServersBoth(int server1, int server2, int udlr_idx1, int udlr_idx2, set<int>& se)
{
  // 切る
  Update_udlr(server1, udlr_idx1, -1);
  Update_udlr(server2, udlr_idx2, -1);

  EraseUnion(server1);
  EraseUnion(server2);
  se.insert(server1);
  se.insert(server2);
}

// サーバー間の経路を設定する共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向)
// useBeam: trueの場合beam.push_backも実行
void SetServerPath(int server1, int server2, int coord, int value, bool isVertical, bool useBeam, vector<P>& beam)
{
  if (isVertical) {
    // 縦方向
    for (int i = gameState.x[server1] + 1; i < gameState.x[server2]; ++i) {
      PushACnt(i, coord);
      gameState.a[i][coord] = value;
      if (useBeam && value == INF) {
        beam.push_back(P(i, coord));
      }
    }
  }
  else {
    // 横方向
    for (int i = gameState.y[server1] + 1; i < gameState.y[server2]; ++i) {
      PushACnt(coord, i);
      gameState.a[coord][i] = value;
      if (useBeam && value == INF) {
        beam.push_back(P(coord, i));
      }
    }
  }
}

// サーバー間を接続する共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向)
void ConnectServers(int server1, int server2, int udlr_idx1, int udlr_idx2, int coord, bool isVertical, set<int>& se)
{
  // 繋ぐ
  Update_udlr(server1, udlr_idx1, server2);
  Update_udlr(server2, udlr_idx2, server1);

  EraseUnion(server1);
  EraseUnion(server2);
  se.insert(server1);

  int aVal = make_edge_value(server1, server2);
  vector<P> dummy; // ConnectServersではbeamを使わない
  SetServerPath(server1, server2, coord, aVal, isVertical, false, dummy);
}

// 中間のサーバーを除外して両端をつなぎなおす共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向) 
void ReconnectServers(int server1, int server2, int middle, int udlr_idx1, int udlr_idx2, int udlr_idx_m1, int udlr_idx_m2, int coord, bool isVertical, set<int>& se)
{
  // 繋ぎなおす
  Update_udlr(server1, udlr_idx1, server2);
  Update_udlr(server2, udlr_idx2, server1);
  Update_udlr(middle, udlr_idx_m1, -1);
  Update_udlr(middle, udlr_idx_m2, -1);

  EraseUnion(middle);
  se.insert(server1);
  se.insert(middle);

  int aVal = make_edge_value(server1, server2);
  vector<P> dummy;
  SetServerPath(server1, server2, coord, aVal, isVertical, false, dummy);
}

// 既存の接続を切って新しく接続しなおす共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向)
void ReconnectServersWithCut(int serverU, int serverD, int newServer, int udlr_U, int udlr_D, int udlr_new1, int udlr_new2, int coord, bool isVertical, set<int>& se)
{
  // 繋ぎなおす
  Update_udlr(serverU, udlr_U, newServer);
  Update_udlr(serverD, udlr_D, newServer);
  Update_udlr(newServer, udlr_new1, serverU);
  Update_udlr(newServer, udlr_new2, serverD);

  EraseUnion(newServer);
  EraseUnion(serverU);

  se.insert(newServer);

  int aVal1 = make_edge_value(serverU, newServer);
  int aVal2 = make_edge_value(newServer, serverD);
  vector<P> dummy;
  SetServerPath(serverU, newServer, coord, aVal1, isVertical, false, dummy);
  SetServerPath(newServer, serverD, coord, aVal2, isVertical, false, dummy);
}

// 既存の2つのサーバー間の接続を切って、間に新しくサーバーを挟む共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向)
void DisconnectAndInsertServer(int server1, int server2, int newServer, int udlr_idx1, int udlr_idx2, int coord, bool isVertical, set<int>& se, vector<P>& beam)
{
  // 切る
  Update_udlr(server1, udlr_idx1, -1);
  Update_udlr(server2, udlr_idx2, -1);

  EraseUnion(server1);
  se.insert(server1);
  se.insert(server2);

  // 新しい経路を設定
  SetServerPath(server1, newServer, coord, INF, isVertical, true, beam);
  SetServerPath(newServer, server2, coord, INF, isVertical, true, beam);
}

// サーバー移動時に既存の接続を切ってINF経路を設定する共通処理
// isVertical: true=縦方向(y方向), false=横方向(x方向)
// movingServer: 移動するサーバー、connectedServer: 接続されていたサーバー
void DisconnectWithPath(int movingServer, int connectedServer, int udlr_idx_moving, int udlr_idx_connected, int coord, bool isVertical, set<int>& se, vector<P>& beam)
{
  // 切る
  Update_udlr(movingServer, udlr_idx_moving, -1);
  Update_udlr(connectedServer, udlr_idx_connected, -1);

  EraseUnion(movingServer);
  se.insert(connectedServer);
  se.insert(movingServer);

  // 経路の向きを考慮：上から下、左から右の方向で設定
  if (isVertical) {
    if (gameState.x[connectedServer] < gameState.x[movingServer]) {
      SetServerPath(connectedServer, movingServer, coord, INF, isVertical, true, beam);
    }
    else {
      SetServerPath(movingServer, connectedServer, coord, INF, isVertical, true, beam);
    }
  }
  else {
    if (gameState.y[connectedServer] < gameState.y[movingServer]) {
      SetServerPath(connectedServer, movingServer, coord, INF, isVertical, true, beam);
    }
    else {
      SetServerPath(movingServer, connectedServer, coord, INF, isVertical, true, beam);
    }
  }
}

// 戻り値：更新したかどうか
int apply_move(double temp, int server_id, int dir, bool forceDo = false, int MethodeMode = 0)
{
  acnt = 0;
  udlrcnt = 0;
  parentcnt = 0;
  vpcnt = 0;

  int nx = gameState.x[server_id] + dx[dir];
  int ny = gameState.y[server_id] + dy[dir];

  int xx = gameState.x[server_id];
  int yy = gameState.y[server_id];

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
      if (gameState.udlr[server_id][2] != -1) {
        if (has_server(xx, yy - 1)) {
          gameState.a[xx][yy] = make_edge_value(server_id, gameState.a[xx][yy - 1]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx][yy - 1];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの右のつながり
      gameState.a[xx][ny] = server_id;

      // 元のマスの上下のつながり
      // 左と繋がっている場合、切る
      if (gameState.udlr[server_id][2] != -1) {
        if (gameState.udlr[server_id][0] != -1) {
          int iteU = gameState.udlr[server_id][0];
          DisconnectWithPath(server_id, iteU, 0, 1, yy, true, se, beam);
        }

        if (gameState.udlr[server_id][1] != -1) {
          int iteD = gameState.udlr[server_id][1];
          DisconnectWithPath(server_id, iteD, 1, 0, yy, true, se, beam);
        }
      }
      else {
        // 上下繋がっていた
        if (gameState.udlr[server_id][0] != -1 && gameState.udlr[server_id][1] != -1) {
          // 繋ぎなおす
          int iteU = gameState.udlr[server_id][0];
          int iteD = gameState.udlr[server_id][1];
          ReconnectServers(iteU, iteD, server_id, 1, 0, 0, 1, yy, true, se);
        }
        // 上繋がっていた
        else if (gameState.udlr[server_id][0] != -1) {
          // 切る
          int iteU = gameState.udlr[server_id][0];
          DisconnectServers(server_id, iteU, 0, 1, se);
          se.insert(server_id);

          SetServerPath(iteU, server_id, yy, INF, true, true, beam);
        }
        // 下繋がっていた
        else if (gameState.udlr[server_id][1] != -1) {
          // 切る
          int iteD = gameState.udlr[server_id][1];
          DisconnectServers(server_id, iteD, 1, 0, se);
          se.insert(server_id);

          SetServerPath(server_id, iteD, yy, INF, true, true, beam);
        }
        // 繋がりなし
        else {
          // 上と下が同じ色の場合
          int iteU = get_server_through_inf(xx, yy, 'U');
          int iteD = get_server_through_inf(xx, yy, 'D');
          if (iteU != -1 && iteD != -1 && getServerType(iteU) == getServerType(iteD)) {
            // 繋ぐ
            ConnectServers(iteU, iteD, 1, 0, yy, true, se);
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[server_id][3] == -1) {
        int iteU = -na / ENCODING_MULTIPLIER;
        int iteD = -na % ENCODING_MULTIPLIER;
        if (gameState.x[iteU] > gameState.x[iteD]) {
          swap(iteU, iteD);
        }

        // 同じ色
        if (getServerType(iteU) == getServerType(server_id)) {
          // 繋ぎなおす
          ReconnectServersWithCut(iteU, iteD, server_id, 1, 0, 0, 1, ny, true, se);
        }
        // 違う色
        else {
          // 切る
          DisconnectAndInsertServer(iteU, iteD, server_id, 1, 0, ny, true, se, beam);
        }
      }
      // 繋がっていない
      else {
        // 上と繋げられるかどうか
        int iteU = get_server_through_inf(xx, ny, 'U');
        if (iteU != -1 && getServerType(iteU) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(iteU, server_id, 1, 0, ny, true, se);
        }

        // 下と繋げられるかどうか
        int iteD = get_server_through_inf(xx, ny, 'D');
        if (iteD != -1 && getServerType(iteD) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(server_id, iteD, 1, 0, ny, true, se);
        }
      }
    }
    // 左
    else {
      // 元のマスの右のつながり
      // 繋がっている
      if (gameState.udlr[server_id][3] != -1) {
        if (has_server(xx, yy + 1)) {
          gameState.a[xx][yy] = make_edge_value(server_id, gameState.a[xx][yy + 1]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx][yy + 1];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの左のつながり
      gameState.a[xx][ny] = server_id;

      // 元のマスの上下のつながり
      // 右と繋がっている場合、切る
      if (gameState.udlr[server_id][3] != -1) {
        if (gameState.udlr[server_id][0] != -1) {
          int iteU = gameState.udlr[server_id][0];
          DisconnectWithPath(server_id, iteU, 0, 1, yy, true, se, beam);
        }

        if (gameState.udlr[server_id][1] != -1) {
          int iteD = gameState.udlr[server_id][1];
          DisconnectWithPath(server_id, iteD, 1, 0, yy, true, se, beam);
        }
      }
      else {
        // 上下繋がっていた
        if (gameState.udlr[server_id][0] != -1 && gameState.udlr[server_id][1] != -1) {
          // 繋ぎなおす
          int iteU = gameState.udlr[server_id][0];
          int iteD = gameState.udlr[server_id][1];
          ReconnectServers(iteU, iteD, server_id, 1, 0, 0, 1, yy, true, se);
        }
        // 上繋がっていた
        else if (gameState.udlr[server_id][0] != -1) {
          // 切る
          int iteU = gameState.udlr[server_id][0];
          DisconnectServersBoth(iteU, server_id, 1, 0, se);

          SetServerPath(iteU, server_id, yy, INF, true, true, beam);
        }
        // 下繋がっていた
        else if (gameState.udlr[server_id][1] != -1) {
          // 切る
          int iteD = gameState.udlr[server_id][1];
          DisconnectServersBoth(iteD, server_id, 0, 1, se);

          SetServerPath(server_id, iteD, yy, INF, true, true, beam);
        }
        // 繋がりなし
        else {
          // 上と下が同じ色の場合
          int iteU = get_server_through_inf(xx, yy, 'U');
          int iteD = get_server_through_inf(xx, yy, 'D');
          if (iteU != -1 && iteD != -1 && getServerType(iteU) == getServerType(iteD)) {
            // 繋ぐ
            ConnectServers(iteU, iteD, 1, 0, yy, true, se);
          }
        }
      }

      // 先のマスの上下のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[server_id][2] == -1) {
        int iteU = -na / ENCODING_MULTIPLIER;
        int iteD = -na % ENCODING_MULTIPLIER;
        if (gameState.x[iteU] > gameState.x[iteD]) {
          swap(iteU, iteD);
        }

        // 同じ色
        if (getServerType(iteU) == getServerType(server_id)) {
          // 繋ぎなおす
          ReconnectServersWithCut(iteU, iteD, server_id, 1, 0, 0, 1, ny, true, se);
        }
        // 違う色
        else {
          // 切る
          DisconnectAndInsertServer(iteU, iteD, server_id, 1, 0, ny, true, se, beam);
        }
      }
      // 繋がっていない
      else {
        // 上と繋げられるかどうか
        int iteU = get_server_through_inf(xx, ny, 'U');
        if (iteU != -1 && getServerType(iteU) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(iteU, server_id, 1, 0, ny, true, se);
        }

        // 下と繋げられるかどうか
        int iteD = get_server_through_inf(xx, ny, 'D');
        if (iteD != -1 && getServerType(iteD) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(server_id, iteD, 1, 0, ny, true, se);
        }
      }
    }
  }
  else {
    // 下
    if (nx == xx + 1) {
      // 元のマスの上のつながり
      // 繋がっている
      if (gameState.udlr[server_id][0] != -1) {
        if (has_server(xx - 1, yy)) {
          gameState.a[xx][yy] = make_edge_value(server_id, gameState.a[xx - 1][yy]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx - 1][yy];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの下のつながり
      gameState.a[nx][yy] = server_id;

      // 元のマスの左右のつながり
      // 上と繋がっている場合、切る
      if (gameState.udlr[server_id][0] != -1) {
        if (gameState.udlr[server_id][2] != -1) {
          int iteL = gameState.udlr[server_id][2];
          DisconnectWithPath(server_id, iteL, 2, 3, xx, false, se, beam);
        }

        if (gameState.udlr[server_id][3] != -1) {
          int iteR = gameState.udlr[server_id][3];
          DisconnectWithPath(server_id, iteR, 3, 2, xx, false, se, beam);
        }
      }
      else {
        // 左右繋がっていた
        if (gameState.udlr[server_id][2] != -1 && gameState.udlr[server_id][3] != -1) {
          // 繋ぎなおす
          int iteL = gameState.udlr[server_id][2];
          int iteR = gameState.udlr[server_id][3];
          ReconnectServers(iteL, iteR, server_id, 3, 2, 2, 3, xx, false, se);
        }
        // 左繋がっていた
        else if (gameState.udlr[server_id][2] != -1) {
          // 切る
          int iteL = gameState.udlr[server_id][2];
          DisconnectServersBoth(iteL, server_id, 3, 2, se);

          SetServerPath(iteL, server_id, xx, INF, false, true, beam);
        }
        // 右繋がっていた
        else if (gameState.udlr[server_id][3] != -1) {
          // 切る
          int iteR = gameState.udlr[server_id][3];
          DisconnectServersBoth(iteR, server_id, 2, 3, se);

          SetServerPath(server_id, iteR, xx, INF, false, true, beam);
        }
        // 繋がりなし
        else {
          // 左と右が同じ色の場合
          int iteL = get_server_through_inf(xx, yy, 'L');
          int iteR = get_server_through_inf(xx, yy, 'R');
          if (iteL != -1 && iteR != -1 && getServerType(iteL) == getServerType(iteR)) {
            // 繋ぐ
            ConnectServers(iteL, iteR, 3, 2, xx, false, se);
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[server_id][1] == -1) {
        int iteL = -na / ENCODING_MULTIPLIER;
        int iteR = -na % ENCODING_MULTIPLIER;
        if (gameState.y[iteL] > gameState.y[iteR]) {
          swap(iteL, iteR);
        }

        // 同じ色
        if (getServerType(iteL) == getServerType(server_id)) {
          // 繋ぎなおす
          ReconnectServersWithCut(iteL, iteR, server_id, 3, 2, 2, 3, nx, false, se);
        }
        // 違う色
        else {
          // 切る
          DisconnectAndInsertServer(iteL, iteR, server_id, 3, 2, nx, false, se, beam);
        }
      }
      // 繋がっていない
      else {
        // 左と繋げられるかどうか
        int iteL = get_server_through_inf(nx, yy, 'L');
        if (iteL != -1 && getServerType(iteL) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(iteL, server_id, 3, 2, nx, false, se);
        }

        // 右と繋げられるかどうか
        int iteR = get_server_through_inf(nx, yy, 'R');
        if (iteR != -1 && getServerType(iteR) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(server_id, iteR, 3, 2, nx, false, se);
        }
      }
    }
    // 上
    else {
      // 元のマスの下のつながり
      // 繋がっている
      if (gameState.udlr[server_id][1] != -1) {
        if (has_server(xx + 1, yy)) {
          gameState.a[xx][yy] = make_edge_value(server_id, gameState.a[xx + 1][yy]);
        }
        else {
          gameState.a[xx][yy] = gameState.a[xx + 1][yy];
        }
      }
      else {
        gameState.a[xx][yy] = INF;
      }

      // 元のマスの上のつながり
      gameState.a[nx][yy] = server_id;

      // 元のマスの左右のつながり
      // 下と繋がっている場合、切る
      if (gameState.udlr[server_id][1] != -1) {
        if (gameState.udlr[server_id][2] != -1) {
          int iteL = gameState.udlr[server_id][2];
          DisconnectWithPath(server_id, iteL, 2, 3, xx, false, se, beam);
        }

        if (gameState.udlr[server_id][3] != -1) {
          int iteR = gameState.udlr[server_id][3];
          DisconnectWithPath(server_id, iteR, 3, 2, xx, false, se, beam);
        }
      }
      else {
        // 左右繋がっていた
        if (gameState.udlr[server_id][2] != -1 && gameState.udlr[server_id][3] != -1) {
          // 繋ぎなおす
          int iteL = gameState.udlr[server_id][2];
          int iteR = gameState.udlr[server_id][3];
          ReconnectServers(iteL, iteR, server_id, 3, 2, 2, 3, xx, false, se);
        }
        // 左繋がっていた
        else if (gameState.udlr[server_id][2] != -1) {
          // 切る
          int iteL = gameState.udlr[server_id][2];
          DisconnectServersBoth(iteL, server_id, 3, 2, se);

          SetServerPath(iteL, server_id, xx, INF, false, true, beam);
        }
        // 右繋がっていた
        else if (gameState.udlr[server_id][3] != -1) {
          // 切る
          int iteR = gameState.udlr[server_id][3];
          DisconnectServersBoth(iteR, server_id, 2, 3, se);

          SetServerPath(server_id, iteR, xx, INF, false, true, beam);
        }
        // 繋がりなし
        else {
          // 左と右が同じ色の場合
          int iteL = get_server_through_inf(xx, yy, 'L');
          int iteR = get_server_through_inf(xx, yy, 'R');
          if (iteL != -1 && iteR != -1 && getServerType(iteL) == getServerType(iteR)) {
            // 繋ぐ
            ConnectServers(iteL, iteR, 3, 2, xx, false, se);
          }
        }
      }

      // 先のマスの左右のつながり
      // 元々繋がっている
      if (na != INF && gameState.udlr[server_id][0] == -1) {
        int iteL = -na / ENCODING_MULTIPLIER;
        int iteR = -na % ENCODING_MULTIPLIER;
        if (gameState.y[iteL] > gameState.y[iteR]) {
          swap(iteL, iteR);
        }

        // 同じ色
        if (getServerType(iteL) == getServerType(server_id)) {
          // 繋ぎなおす
          ReconnectServersWithCut(iteL, iteR, server_id, 3, 2, 2, 3, nx, false, se);
        }
        // 違う色
        else {
          // 切る
          DisconnectAndInsertServer(iteL, iteR, server_id, 3, 2, nx, false, se, beam);
        }
      }
      // 繋がっていない
      else {
        // 左と繋げられるかどうか
        int iteL = get_server_through_inf(nx, yy, 'L');
        if (iteL != -1 && getServerType(iteL) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(iteL, server_id, 3, 2, nx, false, se);
        }

        // 右と繋げられるかどうか
        int iteR = get_server_through_inf(nx, yy, 'R');
        if (iteR != -1 && getServerType(iteR) == getServerType(server_id)) {
          // 繋ぐ
          ConnectServers(server_id, iteR, 3, 2, nx, false, se);
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
      if (gameState.a[px][py] != INF) { continue; }
      int iteL = get_server_through_inf(px, py, 'L');
      int iteR = get_server_through_inf(px, py, 'R');
      if (iteL != -1 && iteR != -1 && getServerType(iteL) == getServerType(iteR)) {
        ConnectServers(iteL, iteR, 3, 2, px, false, se);
      }
    }
  }
  // 上下
  else {
    for (auto&& p : beam) {
      int px = p.first;
      int py = p.second;
      if (gameState.a[px][py] != INF) { continue; }
      int iteU = get_server_through_inf(px, py, 'U');
      int iteD = get_server_through_inf(px, py, 'D');
      if (iteU != -1 && iteD != -1 && getServerType(iteU) == getServerType(iteD)) {
        ConnectServers(iteU, iteD, 1, 0, py, true, se);
      }
    }
  }

  gameState.x[server_id] = nx;
  gameState.y[server_id] = ny;



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
      for (int j = 0; j < 4; ++j) {
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

  double prob = exp((double)diffScore / temp);
  int isDo = 0;
  if (forceDo || prob > rand_01()) {
    isDo = 1;
    gameState.maxScore += diffScore;

    methodCount[1][0]++;
    methodSum[0]++;

    gameState.ans1[gameState.ope1][0] = xx;
    gameState.ans1[gameState.ope1][1] = yy;
    gameState.ans1[gameState.ope1][2] = nx;
    gameState.ans1[gameState.ope1][3] = ny;
    gameState.ans1[gameState.ope1][4] = server_id;
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
    gameState.x[server_id] = xx;
    gameState.y[server_id] = yy;

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
    ite = rand32() % K100;  // 1マス動かすコンピュータ
    dir = rand32() % 4;

    nx = gameState.x[ite] + dx[dir];
    ny = gameState.y[ite] + dy[dir];
    int ok = 0;
    if (0 <= nx && nx < n && 0 <= ny && ny < n && !has_server(nx, ny)) {
      break;
    }
  }

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  apply_move(temp, ite, dir);
}


// 空白を2マス動かす
void Method3(double start_temp, double end_temp, double now_progress)
{
  int xx, yy, dir1, dir2;
  while (true) {
    xx = rand32() % n;
    yy = rand32() % n;
    if (gameState.a[xx][yy] == -1) {
      break;
    }
  }

  dir1 = rand32() % 4;
  dir2 = rand32() % 4;
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
  if (prob > rand_01()) {
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
  int ite = rand32() % K100;  // 1マス動かすコンピュータ
  int dir1 = rand32() % 4;
  int dir2 = rand32() % 4;

  int xx = gameState.x[ite];
  int yy = gameState.y[ite];

  int nx1 = xx + dx[dir1];
  int ny1 = yy + dy[dir1];
  int nx2 = nx1 + dx[dir2];
  int ny2 = ny1 + dy[dir2];

  if (is_out_of_bounds(nx1, ny1) || is_out_of_bounds(nx2, ny2) || gameState.a[nx1][ny1] != -1 ||
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
  if (prob > rand_01()) {
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
  if (gameState.ope1 == 0) { return; }
  int ite = rand32() % gameState.ope1;

  // NGチェック
  // ite以降の操作で、操作元が移動後のマス、操作後が移動前のマス、の操作が出てこなければOK
  for (int i = ite + 1; i < gameState.ope1; ++i) {
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
  for (int i = 0; i < 4; ++i) {
    if (gameState.ans1[ite][2] == gameState.ans1[ite][0] + dx[i] &&
      gameState.ans1[ite][3] == gameState.ans1[ite][1] + dy[i]) {
      reverseDir = (i + 2) % 4;
      break;
    }
  }

  double temp = start_temp + (end_temp - start_temp) * now_progress;
  int isDo = apply_move(temp, gameState.ans1[ite][4], reverseDir, false, 5);

  // 実行した場合、2つ消す
  if (isDo) {
    /*
    修正するもの
    int ope1;
    int gameState.ans1[1000][5];
    */
    // 消すのはiteとope1-1
    gameState.ope1--;
    for (int i = ite; i < gameState.ope1 - 1; ++i) {
      for (int j = 0; j < 5; ++j) { gameState.ans1[i][j] = gameState.ans1[i + 1][j]; }
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
    for (int i = 0; i < gameState.ope1 - 1; ++i) {
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

    for (int i = ite; i < gameState.ope1 - 2; ++i) {
      for (int j = 0; j < 5; ++j) { gameState.ans1[i][j] = gameState.ans1[i + 2][j]; }
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
    for (int i = 0; i < K * 100; ++i) {
      vv[i].clear();
    }
    for (int i = 0; i < gameState.ope1; ++i) {
      int ite = gameState.ans1[i][4];
      vv[ite].push_back(i);
    }

    int it1 = -1;
    int it2 = -1;
    for (int i = 0; i < K * 100; ++i) {
      int sz = vv[i].size();
      for (int j = sz - 2; j >= 0; --j) {
        int ite1 = vv[i][j];
        int ite2 = vv[i][j + 1];
        if (gameState.ans1[ite1][0] == gameState.ans1[ite2][2] && gameState.ans1[ite1][1] == gameState.ans1[ite2][3] &&
          gameState.ans1[ite1][2] == gameState.ans1[ite2][0] && gameState.ans1[ite1][3] == gameState.ans1[ite2][1]) {
          // ngチェック
          int ng = 0;
          for (int k = ite1 + 1; k < ite2; ++k) {
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
      for (int i = it2; i < gameState.ope1 - 1; ++i) {
        for (int j = 0; j < 5; ++j) {
          gameState.ans1[i][j] = gameState.ans1[i + 1][j];
        }
      }
      gameState.ope1--;
      for (int i = it1; i < gameState.ope1 - 1; ++i) {
        for (int j = 0; j < 5; ++j) {
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

int solve(int mode, int problemNum = 0)
{
  Timer timer;
  timer.start();

  init();

  // 愚直解
  gameState.maxScore = gameState.calcScore(K100, true);
  gameState.copyTo(real_GameState);
  gameState.copyTo(seed_GameState);

  // シード作り
  int seedCount = 10;
  for (int tei = 0; tei < seedCount; ++tei) {
    timer.start();

    init();
    gameState.viewOrder = tei % 2;
    gameState.maxScore = gameState.calcScore(K100, true);

    // 焼きなまし
    double now_time = timer.get_elapsed_time();
    double TL = 1.0 / seedCount;
    double now_progress = now_time / TL;
    double start_temp = 10.0 + 10.0 * K100 / (n * n - K100);
    double end_temp = 0;
    int loop = 0;
    int rollbackCount = 0;
    while (true) {
      loop++;
      if (getServerIndex(loop) == 1) {
        now_time = timer.get_elapsed_time();
        now_progress = now_time / TL;
      }
      if (now_progress > 1.0) {
        break;
      }

      // 現在のスコアが悪いときは元に戻す
      if (gameState.maxScore * 1.2 < real_GameState.maxScore || rand32() % 123456 == 0) {
        gameState.copyFrom(real_GameState);
        rollbackCount++;
      }

      int me = 1;

      if (rand32() % 2 == 0) {
        me = 5;
      }

      if (rand32() % 203 == 0) {
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

    clear_all();
  }

  // シードから戻す
  gameState.copyFrom(seed_GameState);
  gameState.copyTo(real_GameState);

  // 焼きなまし
  timer.start();
  double now_time = timer.get_elapsed_time();
  double TL = 1.9 / outer_Split;
  double now_progress = now_time / TL;
  double start_temp = 20.0 + 10.0 * K100 / (n * n - K100);
  double end_temp = 0;
  int loop = 0;
  int rollbackCount = 0;
  while (true) {
    loop++;
    if (getServerIndex(loop) == 1) {
      now_time = timer.get_elapsed_time();
      now_progress = now_time / TL;
    }
    if (now_progress > 1.0) {
      break;
    }

    // 現在のスコアが悪いときは元に戻す
    if (gameState.maxScore * 1.2 < real_GameState.maxScore || rand32() % 123456 == 0) {
      gameState.copyFrom(real_GameState);
      rollbackCount++;
    }

    int me = 1;

    if (rand32() % 2 == 0) {
      me = 5;
    }

    if (rand32() % 203 == 0) {
      me = 6;
    }

    if (rand32() % 1011 == 0) {
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
    for (int i = 1; i < 8; ++i) {
      cerr << "Method" << i << " = " << methodCount[i][0] << " / "
        << methodCount[i][1] << endl;
    }
    cerr << "MethodSum = " << methodSum[0] << " / " << methodSum[1] << endl;
    cerr << endl;
  }

  return gameState.maxScore;
}

int solve_outer(int mode, int problemNum)
{
  // 入力部
  input_data(problemNum);

  for (int _ = 0; _ < outer_Split; ++_) {
    gameState.viewOrder = _ % 2;
    int score = solve(mode, problemNum);
    if (score >= outer_GameState.maxScore) {
      gameState.copyTo(outer_GameState);
    }
    clear_all();
  }
  gameState.copyFrom(outer_GameState);

  // 解の出力
  output_data(mode, problemNum);

  return gameState.maxScore;
}

int main()
{
  int mode = 2;

  if (mode == 0) {
    solve_outer(mode, 0);
  }
  else if (mode == 1) {
    solve_outer(mode, 1);
  }
  else if (mode == 2) {
    int sum = 0;
    for (int i = 0; i < 10; ++i) {
      sum += solve_outer(mode, i);
      clear_outer_all();
    }
    cout << "sum = " << sum << endl;
  }
  else if (mode == 3) {
    int problemNum = 0;
    cin >> problemNum;
    solve_outer(mode, problemNum);
  }

  return 0;
}
