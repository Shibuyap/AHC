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

struct CropType
{
  int id;
  int plantingMonth;
  int harvestMonth;
};

bool compareByHarvestMonth(const CropType& a, const CropType& b) { return a.harvestMonth < b.harvestMonth; }
bool compareByPlantingMonth(const CropType& a, const CropType& b) { return a.plantingMonth < b.plantingMonth; }

// グローバル変数を名前空間に整理
namespace GlobalVars
{
  // 入力
  int entranceX;
  int horizontalCanals[GRID_SIZE][GRID_SIZE];
  int verticalCanals[GRID_SIZE][GRID_SIZE];
  int totalCropTypes;
  int plantingMonth[MAX_K], harvestMonth[MAX_K];
  int canalsAround[GRID_SIZE][GRID_SIZE][4];
  vector<CropType> plantableCrops[T], harvestableCrops[T];
  int isPlanted[MAX_K];

  // 出力
  int totalPlantedCrops;
  int plantedCropId[MAX_K], plantedX[MAX_K], plantedY[MAX_K], plantedMonth[MAX_K];

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
  totalPlantedCrops = 0;
  for (int i = 0; i < T; ++i) {
    plantableCrops[i].clear();
    harvestableCrops[i].clear();
  }
  for (int i = 0; i < MAX_K; ++i) {
    isPlanted[i] = 0;
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
  stream >> TTT >> HHH >> WWW >> entranceX;

  for (int i = 0; i < H - 1; ++i) {
    string str;
    stream >> str;
    for (int j = 0; j < W; ++j) {
      horizontalCanals[i][j] = str[j] - '0';
    }
  }

  for (int i = 0; i < H; ++i) {
    string str;
    stream >> str;
    for (int j = 0; j < W - 1; ++j) {
      verticalCanals[i][j] = str[j] - '0';
    }
  }

  stream >> totalCropTypes;
  for (int i = 0; i < totalCropTypes; ++i) {
    stream >> plantingMonth[i] >> harvestMonth[i];
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
        canalsAround[i][j][k] = 0;
      }

      // 上の壁
      canalsAround[i][j][0] = (i == 0) || (horizontalCanals[i - 1][j] != 0);

      // 左の壁
      canalsAround[i][j][1] = (j == 0) || (verticalCanals[i][j - 1] != 0);

      // 下の壁
      canalsAround[i][j][2] = (i == H - 1) || (horizontalCanals[i][j] != 0);

      // 右の壁
      canalsAround[i][j][3] = (j == W - 1) || (verticalCanals[i][j] != 0);
    }
  }

  for (int i = 0; i < totalCropTypes; ++i) {
    plantingMonth[i]--;
    harvestMonth[i]--;
    CropType crop;
    crop.id = i + 1;
    crop.plantingMonth = plantingMonth[i];
    crop.harvestMonth = harvestMonth[i];
    for (int l = 0; l < (1); ++l) {
      if (plantingMonth[i] - l < 0) { break; }
      plantableCrops[plantingMonth[i] - l].push_back(crop);
    }
    harvestableCrops[harvestMonth[i]].push_back(crop);
  }

  for (int i = 0; i < (T); ++i) {
    sort(plantableCrops[i].begin(), plantableCrops[i].end(), compareByHarvestMonth);
    sort(harvestableCrops[i].begin(), harvestableCrops[i].end(), compareByPlantingMonth);
  }
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
  for (int i = 0; i < (totalPlantedCrops); ++i) {
    sum += harvestMonth[plantedCropId[i] - 1] - plantingMonth[plantedCropId[i] - 1] + 1;
  }
  ll res = SCORE_BASE * sum / HWT;
  return res;
}

// 初期解生成
void initializeSolution() {}

bool isValidPlotForCrop(const int x, const int y, const vector<vector<int>>& fieldStatus)
{
  if (x == entranceX && y == SY) {
    return false;
  }
  int blankCount = 0;
  for (int i = 0; i < (4); ++i) {
    if (canalsAround[x][y][i]) { continue; }
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (fieldStatus[nx][ny] == -1) {
      blankCount++;
    }
    else if (fieldStatus[nx][ny] < fieldStatus[x][y]) {
      return false;
    }
  }
  if (blankCount == 1) {
    return true;
  }
  return false;
}

bool isFieldAccessible(const vector<vector<int>>& fieldStatus)
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
      if (fieldStatus[i][j] == -1) {
        cnt++;
      }
    }
  }

  int head = 0;
  int tail = 0;
  if (fieldStatus[entranceX][SY] == -1) {
    bfsQueue[tail][0] = entranceX;
    bfsQueue[tail][1] = SY;
    tail++;
    f[entranceX][SY] = 1;
    cnt--;
  }
  while (head < tail) {
    int x = bfsQueue[head][0];
    int y = bfsQueue[head][1];
    head++;
    for (int i = 0; i < (4); ++i) {
      if (canalsAround[x][y][i]) {
        continue;
      }
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (f[nx][ny] == 0 && fieldStatus[nx][ny] == -1) {
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

bool canReachAllPlotsFromEntrance(const vector<vector<int>>& fieldStatus)
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
  f[entranceX][SY] = fieldStatus[entranceX][SY];
  cnt++;
  int hash = f[entranceX][SY] + 1;
  dijkstraQueue[hash][tails[hash]][0] = entranceX;
  dijkstraQueue[hash][tails[hash]][1] = SY;
  tails[hash]++;

  for (int turn = 0; turn < (T); ++turn) {
    while (heads[turn] < tails[turn]) {
      int x = dijkstraQueue[turn][heads[turn]][0];
      int y = dijkstraQueue[turn][heads[turn]][1];
      int val = turn - 1;
      heads[turn]++;
      for (int i = 0; i < (4); ++i) {
        if (canalsAround[x][y][i]) {
          continue;
        }
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (f[nx][ny] == INF && val <= fieldStatus[nx][ny]) {
          f[nx][ny] = fieldStatus[nx][ny];
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

double Score_1_Dijkstra(const int sx, const int sy, const int d, const vector<vector<int>>& fieldStatus, const int walkCount)
{
  int f[H][W];
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      f[i][j] = INF;
    }
  }


  int heads[DIJKSTRA_QUEUE_SIZE] = {};
  int tails[DIJKSTRA_QUEUE_SIZE] = {};
  f[sx][sy] = fieldStatus[sx][sy];

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
        if (canalsAround[x][y][i]) {
          continue;
        }
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (f[nx][ny] == INF && val <= fieldStatus[nx][ny]) {
          f[nx][ny] = fieldStatus[nx][ny];
          if (fieldStatus[nx][ny] != -1) {
            days.push_back(fieldStatus[nx][ny]);
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
double Score_2(const int sx, const int sy, const int d, vector<vector<int>>& fieldStatus)
{
  double score = 100;
  fieldStatus[sx][sy] = d;

  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      if (fieldStatus[i][j] != -1) {
        for (int k = 0; k < (4); ++k) {
          if (canalsAround[i][j][k]) { continue; }
          if (fieldStatus[i + dx[k]][j + dy[k]] != -1) {
            if (fieldStatus[i + dx[k]][j + dy[k]] == fieldStatus[i][j]) {
              score += 1e5;
            }

            if (abs(fieldStatus[i + dx[k]][j + dy[k]] - fieldStatus[i][j]) >= 1) {
              score -= 1e4;
            }
          }
        }
      }
    }
  }
  fieldStatus[sx][sy] = -1;
  return score;
}


void planCropPlacement()
{
  vector<vector<int>> fieldStatus(20, vector<int>(20));
  for (int i = 0; i < (H); ++i) {
    for (int j = 0; j < (W); ++j) {
      fieldStatus[i][j] = -1;
    }
  }

  for (int turn = 0; turn < (T); ++turn) {
    // 設置
    drep(i, plantableCrops[turn].size())
    {
      CropType crop = plantableCrops[turn][i];
      if (isPlanted[crop.id]) {
        continue;
      }
      int d = crop.harvestMonth;

      // 関節点じゃない空白マスを列挙
      LowLinkGraph Graph;
      map<int, P> mp;     // {頂点番号,座標}
      map<P, int> mpInv;  // {座標,頂点番号}
      // 空のグラフ作成
      for (int i = 0; i < HW; ++i) {
        auto [x, y] = toCoord(i);
        if (fieldStatus[x][y] == -1) {
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
          if (canalsAround[x][y][j]) { continue; }
          int nx = x + dx[j];
          int ny = y + dy[j];
          if (fieldStatus[nx][ny] == -1) {
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
      walkCount[entranceX][SY] = 0;
      queue<P> que;
      que.push(P(entranceX, SY));
      while (que.size()) {
        int x = que.front().first;
        int y = que.front().second;
        que.pop();
        for (int j = 0; j < (4); ++j) {
          if (canalsAround[x][y][j]) { continue; }
          int nx = x + dx[j];
          int ny = y + dy[j];
          if (fieldStatus[nx][ny] == -1 && walkCount[x][y] + 1 < walkCount[nx][ny]) {
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
        if (fieldStatus[x][y] != -1) {
          assert(false);
        }

        fieldStatus[x][y] = d;
        bool OK = false;
        if (isValidPlotForCrop(x, y, fieldStatus)) {
          OK = true;
        }
        if (!OK) {
          if (isFieldAccessible(fieldStatus) && canReachAllPlotsFromEntrance(fieldStatus)) {
            OK = true;
          }
        }
        if (OK) {
          OKs.push_back(P(x, y));
        }
        fieldStatus[x][y] = -1;
      }
      if (!OKs.empty()) {
        double ma = 0;
        P best;
        best.first = -1;
        for (auto p : OKs) {
          double tmpScore = Score_1_Dijkstra(p.first, p.second, d, fieldStatus, walkCount[p.first][p.second]);
          //double tmpScore = Score_2(p.first, p.second, d, fieldStatus);
          if (ma < tmpScore) {
            ma = tmpScore;
            best = p;
          }
        }
        if (best.first != -1) {
          int x = best.first;
          int y = best.second;
          fieldStatus[x][y] = d;
          plantedCropId[totalPlantedCrops] = crop.id;
          plantedX[totalPlantedCrops] = x;
          plantedY[totalPlantedCrops] = y;
          plantedMonth[totalPlantedCrops] = turn;
          isPlanted[crop.id] = 1;
          totalPlantedCrops++;
        }
      }
    }

    // 収穫
    for (int i = 0; i < (H); ++i) {
      for (int j = 0; j < (W); ++j) {
        if (fieldStatus[i][j] == turn) {
          fieldStatus[i][j] = -1;
        }
      }
    }
  }
}

// 解答出力
void outputSolution(ofstream& ofs)
{
  if (mode == 0) {
    cout << totalPlantedCrops << endl;
    for (int i = 0; i < (totalPlantedCrops); ++i) {
      cout << plantedCropId[i] << ' ' << plantedX[i] << ' ' << plantedY[i] << ' ' << plantedMonth[i] + 1
        << endl;
    }
  }
  else {
    ofs << totalPlantedCrops << endl;
    for (int i = 0; i < (totalPlantedCrops); ++i) {
      ofs << plantedCropId[i] << ' ' << plantedX[i] << ' ' << plantedY[i] << ' ' << plantedMonth[i] + 1
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
  planCropPlacement();

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
