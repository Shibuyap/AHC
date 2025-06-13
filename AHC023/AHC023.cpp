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
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
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
}  // namespace

// 定数
namespace Constants
{
  constexpr int T = 100;
  constexpr int H = 20;
  constexpr int W = 20;
  constexpr int HW = H * W;
  constexpr int HWT = HW * T;
  constexpr int SY = 0;
  constexpr int INF = 1001001;
  constexpr double TL = 1.8;
  constexpr int MAX_K = 44000;
  constexpr int GRID_SIZE = 25;
  constexpr int QUEUE_SIZE = 10000;
  constexpr int DIJKSTRA_QUEUE_SIZE = 110;
  constexpr ll SCORE_BASE = 1000000LL;
}

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
using namespace Constants;

// 座標変換ヘルパー関数
inline int toIndex(int x, int y)
{
  return x * H + y;
}

inline P toCoord(int index)
{
  return P(index / H, index % H);
}

std::mt19937 engine;
int mode;

/////////////////////////////////////////////////////////////////////////////////
// 構造体・比較関数

struct Crop
{
  int id;
  int s;
  int d;
};

bool compD(const Crop& a, const Crop& b) { return a.d < b.d; }
bool compS(const Crop& a, const Crop& b) { return a.s < b.s; }

// グローバル変数を名前空間に整理
namespace GlobalVars
{
  // 入力
  int SX;
  int ho[GRID_SIZE][GRID_SIZE];
  int ve[GRID_SIZE][GRID_SIZE];
  int K;
  int S[MAX_K], D[MAX_K];
  int wall[GRID_SIZE][GRID_SIZE][4];
  vector<Crop> SVec[T], DVec[T];
  int Used[MAX_K];

  // 出力
  int M;
  int ansK[MAX_K], ansX[MAX_K], ansY[MAX_K], ansS[MAX_K];

  // キュー
  int bfsQueue[QUEUE_SIZE][2];
  int dijkstraQueue[DIJKSTRA_QUEUE_SIZE][HW + 10][2];
}

using namespace GlobalVars;

/////////////////////////////////////////////////////////////////////////////////

namespace /* 関節点ライブラリ */
{
  struct LowLinkEdge
  {
    int to;
  };
  using LowLinkGraph = vector<vector<LowLinkEdge>>;

  /* Lowlink: グラフの関節点・橋を列挙する構造体
      作成: O(E+V)
      関節点の集合: vector<int> aps
      橋の集合: vector<P> bridges
  */
  struct LowLink
  {
    const LowLinkGraph& G;
    vector<int> used, ord, low;
    vector<int> aps;  // articulation points
    vector<P> bridges;

    LowLink(const LowLinkGraph& G_) : G(G_)
    {
      used.assign(G.size(), 0);
      ord.assign(G.size(), 0);
      low.assign(G.size(), 0);
      int k = 0;
      for (int i = 0; i < (int)G.size(); i++) {
        if (!used[i]) { k = dfs(i, k, -1); }
      }
      sort(aps.begin(), aps.end());          // 必要ならソートする
      sort(bridges.begin(), bridges.end());  // 必要ならソートする
    }

    // id:探索中の頂点, k:dfsで何番目に探索するか, par:idの親
    int dfs(int id, int k, int par)
    {
      used[id] = true;
      ord[id] = k++;
      low[id] = ord[id];
      bool is_aps = false;
      int count = 0;  // 子の数
      for (auto& e : G[id]) {
        if (!used[e.to]) {
          count++;
          k = dfs(e.to, k, id);
          low[id] = min(low[id], low[e.to]);
          if (par != -1 && ord[id] <= low[e.to]) {
            is_aps = true;  // 条件2を満たすので関節点
          }
          if (ord[id] < low[e.to]) {
            bridges.emplace_back(min(id, e.to), max(id, e.to));
          }
        }
        else if (e.to != par) {  // eが後退辺の時
          low[id] = min(low[id], ord[e.to]);
        }
      }
      if (par == -1 && count >= 2) { is_aps = true; }  // 条件1を満たすので関節点
      if (is_aps) { aps.push_back(id); }
      return k;
    }
  };
}  // namespace

// 複数ケース回すときに内部状態を初期値に戻す
void resetGlobalState()
{
  M = 0;
  for (int i = 0; i < T; ++i) {
    SVec[i].clear();
    DVec[i].clear();
  }
  for (int i = 0; i < MAX_K; ++i) {
    Used[i] = 0;
  }
}

// 数値をゼロ埋め4桁の文字列に変換
string formatProblemNumber(int problemNum)
{
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  return strNum;
}

// 入力データを読み込む共通処理
template<typename Stream>
void readInputData(Stream& stream)
{
  int TTT, HHH, WWW;
  stream >> TTT >> HHH >> WWW >> SX;

  for (int i = 0; i < H - 1; ++i) {
    string str;
    stream >> str;
    for (int j = 0; j < W; ++j) {
      ho[i][j] = str[j] - '0';
    }
  }

  for (int i = 0; i < H; ++i) {
    string str;
    stream >> str;
    for (int j = 0; j < W - 1; ++j) {
      ve[i][j] = str[j] - '0';
    }
  }

  stream >> K;
  for (int i = 0; i < K; ++i) {
    stream >> S[i] >> D[i];
  }
}

// 入力受け取り
void readInput(int problemNum)
{
  string fileNameIfs = "./in/" + formatProblemNumber(problemNum) + ".txt";
  ifstream ifs(fileNameIfs);

  if (!ifs.is_open()) {
    readInputData(cin);
  }
  else {
    readInputData(ifs);
  }

  // 壁情報を初期化
  for (int i = 0; i < H; ++i) {
    for (int j = 0; j < W; ++j) {
      for (int k = 0; k < 4; ++k) {
        wall[i][j][k] = 0;
      }

      // 上の壁
      wall[i][j][0] = (i == 0) || (ho[i - 1][j] != 0);

      // 左の壁
      wall[i][j][1] = (j == 0) || (ve[i][j - 1] != 0);

      // 下の壁
      wall[i][j][2] = (i == H - 1) || (ho[i][j] != 0);

      // 右の壁
      wall[i][j][3] = (j == W - 1) || (ve[i][j] != 0);
    }
  }

  for (int i = 0; i < (K); ++i) {
    S[i]--;
    D[i]--;
    Crop crop;
    crop.id = i + 1;
    crop.s = S[i];
    crop.d = D[i];
    for (int l = 0; l < (1); ++l) {
      if (S[i] - l < 0) { break; }
      SVec[S[i] - l].push_back(crop);
    }
    DVec[D[i]].push_back(crop);
  }

  for (int i = 0; i < (T); ++i) {
    sort(SVec[i].begin(), SVec[i].end(), compD);
    sort(DVec[i].begin(), DVec[i].end(), compS);
  }

  // for (int i = 0; i < (T); ++i) {
  //   cout << i << " : ";
  //   for (auto x : SVec[i]) {
  //     cout << x.d << ' ';
  //   }
  //   cout << endl;
  // }
}

// 出力ファイルストリームオープン
void openOutputFile(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/" + formatProblemNumber(probNum) + ".txt";
    ofs.open(fileNameOfs);
  }
}

// スコア計算
ll calculateScore()
{
  ll sum = 0;
  for (int i = 0; i < (M); ++i) {
    sum += D[ansK[i] - 1] - S[ansK[i] - 1] + 1;
  }
  ll res = SCORE_BASE * sum / HWT;
  return res;
}

// 初期解生成
void initializeSolution() {}

bool isValidPlacement(const int x, const int y, const vector<vector<int>>& use)
{
  if (x == SX && y == SY) {
    return false;
  }
  int blankCount = 0;
  for (int i = 0; i < (4); ++i) {
    if (wall[x][y][i]) { continue; }
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (use[nx][ny] == -1) {
      blankCount++;
    }
    else if (use[nx][ny] < use[x][y]) {
      return false;
    }
  }
  if (blankCount == 1) {
    return true;
  }
  return false;
}

bool isConnectedBFS(const vector<vector<int>>& use)
{
  int f[H][W];
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      f[i][j] = 0;
    }
  }
  int cnt = 0;
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      if (use[i][j] == -1) {
        cnt++;
      }
    }
  }

  int head = 0;
  int tail = 0;
  if (use[SX][SY] == -1) {
    bfsQueue[tail][0] = SX;
    bfsQueue[tail][1] = SY;
    tail++;
    f[SX][SY] = 1;
    cnt--;
  }
  while (head < tail) {
    int x = bfsQueue[head][0];
    int y = bfsQueue[head][1];
    head++;
    for (int i = 0; i < (4); ++i) {
      if (wall[x][y][i]) {
        continue;
      }
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (f[nx][ny] == 0 && use[nx][ny] == -1) {
        f[nx][ny] = 1;
        cnt--;
        bfsQueue[tail][0] = nx;
        bfsQueue[tail][1] = ny;
        tail++;
      }
    }
  }
  if (cnt == 0) {
    return true;
  }
  return false;
}

bool isReachableFromStart(const vector<vector<int>>& use)
{
  int f[H][W];
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      f[i][j] = INF;
    }
  }

  int cnt = 0;

  int heads[DIJKSTRA_QUEUE_SIZE] = {};
  int tails[DIJKSTRA_QUEUE_SIZE] = {};
  f[SX][SY] = use[SX][SY];
  cnt++;
  int hash = f[SX][SY] + 1;
  dijkstraQueue[hash][tails[hash]][0] = SX;
  dijkstraQueue[hash][tails[hash]][1] = SY;
  tails[hash]++;

  for (int turn = 0; turn < (T); ++turn) {
    while (heads[turn] < tails[turn]) {
      int x = dijkstraQueue[turn][heads[turn]][0];
      int y = dijkstraQueue[turn][heads[turn]][1];
      int val = turn - 1;
      heads[turn]++;
      for (int i = 0; i < (4); ++i) {
        if (wall[x][y][i]) {
          continue;
        }
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (f[nx][ny] == INF && val <= use[nx][ny]) {
          f[nx][ny] = use[nx][ny];
          cnt++;
          hash = f[nx][ny] + 1;
          dijkstraQueue[hash][tails[hash]][0] = nx;
          dijkstraQueue[hash][tails[hash]][1] = ny;
          tails[hash]++;
        }
      }
    }
  }

  if (cnt == HW) {
    return true;
  }
  return false;
}

double Score_1_Dijkstra(const int sx, const int sy, const int d, const vector<vector<int>>& use, const int walkCount)
{
  int f[H][W];
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      f[i][j] = INF;
    }
  }


  int heads[DIJKSTRA_QUEUE_SIZE] = {};
  int tails[DIJKSTRA_QUEUE_SIZE] = {};
  f[sx][sy] = use[sx][sy];

  int hash = f[sx][sy] + 1;
  dijkstraQueue[hash][tails[hash]][0] = sx;
  dijkstraQueue[hash][tails[hash]][1] = sy;
  tails[hash]++;

  vector<int> days;
  int kind[DIJKSTRA_QUEUE_SIZE] = {};
  for (int turn = 0; turn < (T); ++turn) {
    while (heads[turn] < tails[turn]) {
      kind[turn] = 1;
      int x = dijkstraQueue[turn][heads[turn]][0];
      int y = dijkstraQueue[turn][heads[turn]][1];
      int val = turn - 1;
      heads[turn]++;
      for (int i = 0; i < (4); ++i) {
        if (wall[x][y][i]) {
          continue;
        }
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (f[nx][ny] == INF && val <= use[nx][ny]) {
          f[nx][ny] = use[nx][ny];
          if (use[nx][ny] != -1) {
            days.push_back(use[nx][ny]);
          }
          hash = f[nx][ny] + 1;
          dijkstraQueue[hash][tails[hash]][0] = nx;
          dijkstraQueue[hash][tails[hash]][1] = ny;
          tails[hash]++;
        }
      }
    }
  }

  double score = 1e9;
  for (int i = 0; i < (days.size()); ++i) {
    int day = days[i];
    if (d == day) {
      score += 1e5 / ((i + 1) * (i + 1));
    }
    else {
      score -= (double)abs(d - day) / ((i + 1) * (i + 1));
    }
  }
  for (int i = 0; i < DIJKSTRA_QUEUE_SIZE; ++i) {
    score += kind[i] * 1e6;
  }
  score += walkCount * 1e5;
  return score;
}

// 見えてる種類は多く、見えてる数は少なく
double Score_2(const int sx, const int sy, const int d, vector<vector<int>>& use)
{
  double score = 100;
  use[sx][sy] = d;

  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      if (use[i][j] != -1) {
        for (int k = 0; k < (4); ++k) {
          if (wall[i][j][k]) { continue; }
          if (use[i + dx[k]][j + dy[k]] != -1) {
            if (use[i + dx[k]][j + dy[k]] == use[i][j]) {
              score += 1e5;
            }

            if (abs(use[i + dx[k]][j + dy[k]] - use[i][j]) >= 1) {
              score -= 1e4;
            }
          }
        }
      }
    }
  }
  use[sx][sy] = -1;
  return score;
}


void greedyPlacement()
{
  vector<vector<int>> use(20, vector<int>(20));
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      use[i][j] = -1;
    }
  }

  for (int turn = 0; turn < (T); ++turn) {
    // 設置
    drep(i, SVec[turn].size())
    {
      Crop crop = SVec[turn][i];
      if (Used[crop.id]) {
        continue;
      }
      int d = crop.d;

      // 関節点じゃない空白マスを列挙
      LowLinkGraph Graph;
      map<int, P> mp;     // {頂点番号,座標}
      map<P, int> mpInv;  // {座標,頂点番号}
      // 空のグラフ作成
      for (int i = 0; i < HW; ++i) {
        auto [x, y] = toCoord(i);
        if (use[x][y] == -1) {
          int num = mp.size();
          mp[num] = P(x, y);
          mpInv[P(x, y)] = num;
          Graph.push_back(vector<LowLinkEdge>());
        }
      }
      // エッジ
      for (int num = 0; num < mp.size(); ++num) {
        int x = mp[num].first;
        int y = mp[num].second;
        for (int j = 0; j < (4); ++j) {
          if (wall[x][y][j]) { continue; }
          int nx = x + dx[j];
          int ny = y + dy[j];
          if (use[nx][ny] == -1) {
            LowLinkEdge e;
            e.to = mpInv[P(nx, ny)];
            Graph[num].push_back(e);
          }
        }
      }
      // 関節点列挙
      LowLink lowLink(Graph);
      set<int> aps;
      for (auto ap : lowLink.aps) {
        aps.insert(ap);
      }
      vector<P> blanks;
      for (int i = 0; i < mp.size(); ++i) {
        if (aps.find(i) != aps.end()) {
          continue;
        }
        blanks.push_back(mp[i]);
      }
      std::shuffle(blanks.begin(), blanks.end(), engine);

      int walkCount[H][W];
      for (int j = 0; j < (H); ++j) {
        for (int k = 0; k < (W); ++k) {
          walkCount[j][k] = INF;
        }
      }
      walkCount[SX][SY] = 0;
      queue<P> que;
      que.push(P(SX, SY));
      while (que.size()) {
        int x = que.front().first;
        int y = que.front().second;
        que.pop();
        for (int j = 0; j < (4); ++j) {
          if (wall[x][y][j]) { continue; }
          int nx = x + dx[j];
          int ny = y + dy[j];
          if (use[nx][ny] == -1 && walkCount[x][y] + 1 < walkCount[nx][ny]) {
            walkCount[nx][ny] = walkCount[x][y] + 1;
            que.push(P(nx, ny));
          }
        }
      }

      // 配置決め
      vector<P> OKs;
      for (auto ap : blanks) {
        int x = ap.first;
        int y = ap.second;
        if (use[x][y] != -1) {
          assert(false);
        }

        use[x][y] = d;
        bool OK = false;
        if (isValidPlacement(x, y, use)) {
          OK = true;
        }
        if (!OK) {
          if (isConnectedBFS(use) && isReachableFromStart(use)) {
            OK = true;
          }
        }
        if (OK) {
          OKs.push_back(P(x, y));
        }
        use[x][y] = -1;
      }
      if (!OKs.empty()) {
        double ma = 0;
        P best;
        best.first = -1;
        for (auto p : OKs) {
          double tmpScore = Score_1_Dijkstra(p.first, p.second, d, use, walkCount[p.first][p.second]);
          //double tmpScore = Score_2(p.first, p.second, d, use);
          if (ma < tmpScore) {
            ma = tmpScore;
            best = p;
          }
        }
        if (best.first != -1) {
          int x = best.first;
          int y = best.second;
          use[x][y] = d;
          ansK[M] = crop.id;
          ansX[M] = x;
          ansY[M] = y;
          ansS[M] = turn;
          Used[crop.id] = 1;
          M++;
        }
      }
    }

    // 収穫
    for (int i = 0; i < (H); ++i) {
      for (int j = 0; j < (W); ++j) {
        if (use[i][j] == turn) {
          use[i][j] = -1;
        }
      }
    }
  }
}

// 解答出力
void outputSolution(ofstream& ofs)
{
  if (mode == 0) {
    cout << M << endl;
    for (int i = 0; i < (M); ++i) {
      cout << ansK[i] << ' ' << ansX[i] << ' ' << ansY[i] << ' ' << ansS[i] + 1
        << endl;
    }
  }
  else {
    ofs << M << endl;
    for (int i = 0; i < (M); ++i) {
      ofs << ansK[i] << ' ' << ansX[i] << ' ' << ansY[i] << ' ' << ansS[i] + 1
        << endl;
    }
  }
}

ll solveProblem(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  resetGlobalState();

  // 入力受け取り
  readInput(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  openOutputFile(probNum, ofs);

  // 初期解生成
  initializeSolution();

  // 貪欲1
  greedyPlacement();

  // 解答を出力
  outputSolution(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calculateScore();
  }
  return score;
}

int main()
{
  std::random_device rnd;
  engine.seed(rnd());

  mode = 1;

  if (mode == 0) {
    solveProblem(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 10)
    {
      for (int j = 0; j < (1); ++j) {
        ll score = solveProblem(i);
        sum += score;
        cout << "num = " << i << ", ";
        cout << "score = " << score << ", ";
        cout << "sum = " << sum << endl;
      }
    }
  }

  return 0;
}
