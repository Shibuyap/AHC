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

// ループの簡略化マクロ
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

class UnionFind
{
public:
  const int MAX_UF = 1000006;

  int N;              // 頂点数
  vector<int> Par;    // 親
  vector<int> Rank;   // 木の深さ
  vector<int> Count;  // 属する頂点の個数(親のみ正しい)

  // 初期化
  void Initialize()
  {
    for (int i = 0; i < N; i++) {
      Par[i]   = i;
      Rank[i]  = 0;
      Count[i] = 1;
    }
  }

  UnionFind(int n)
  {
    N = n;
    Par.resize(N);
    Rank.resize(N);
    Count.resize(N);

    Initialize();
  }

  UnionFind()
  {
    UnionFind(MAX_UF);
  }

  // 木の根を求める
  int Find(int x)
  {
    if (Par[x] == x) {
      return x;
    }
    else {
      return Par[x] = Find(Par[x]);
    }
  }

  // xとyの属する集合を併合
  void Unite(int x, int y)
  {
    x = Find(x);
    y = Find(y);
    if (x == y) return;

    if (Rank[x] < Rank[y]) {
      Par[x] = y;
      Count[y] += Count[x];
    }
    else {
      Par[y] = x;
      Count[x] += Count[y];
      if (Rank[x] == Rank[y]) Rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool IsSame(int x, int y)
  {
    return Find(x) == Find(y);
  }

  // xの属する集合のサイズ
  int GetCount(int x)
  {
    return Count[Find(x)];
  }
};

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t RandXor()
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() { return (RandXor() + 0.5) * (1.0 / UINT_MAX); }

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化（使用されていない）
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// 非常に大きな値
const int INF = 1001001001;

// 移動方向の配列(上,左,下,右)
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int ex[13] = { -2,-1,-1,-1,0,0,0,0,0,1,1,1,2 };
const int ey[13] = { 0,-1,0,1,-2,-1,0,1,2,-1,0,1,0 };

double TL = 2.8; // 時間制限（Time Limit）
int mode;        // 実行モード
std::chrono::steady_clock::time_point startTimeClock; // 時間計測用

// 時間計測をリセットする関数
void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

const int MAX_M = 1600;

const int n = 50;
int m;
int k;
const int T = 800;

vector<int> sx, sy, tx, ty;
vector<int> a[n][n];
vector<int> b[n][n];

struct Turn
{
  int num;
  int x;
  int y;
};

struct Ans
{
  vector<int> order;
  Turn turns[T];
  string comment[T];
  int score;
  int stationCount;
};

Ans ans;
Ans best_ans;

int board[n][n];
int turn = 0;
int money = 0;
void ClearBoard()
{
  turn = 0;
  rep(i, n)
  {
    rep(j, n)
    {
      board[i][j] = -1;
    }
  }
}

void Action(int num, int x, int y, int income)
{
  ans.comment[turn] = "income = " + to_string(income);
  ans.turns[turn].num = num;
  ans.turns[turn].x = x;
  ans.turns[turn].y = y;
  if (num != -1) {
    board[x][y] = num;
  }
  if (num == 0) {
    ans.stationCount++;
  }
  turn++;
}

void CopyToBest()
{
  best_ans = ans;
}

void CopyToAns()
{
  ans = best_ans;
}

bool IsNG(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

int GetRailNum(int x1, int y1, int x2, int y2, int x3, int y3)
{
  int dir1 = 0, dir2 = 0;
  rep(i, 4)
  {
    int nx = x2 + dx[i];
    int ny = y2 + dy[i];
    if (nx == x1 && ny == y1) dir1 = i;
    if (nx == x3 && ny == y3) dir2 = i;
  }

  if (dir1 > dir2)swap(dir1, dir2);

  if (dir1 == 0 && dir2 == 1) {
    return 4;
  }
  else if (dir1 == 0 && dir2 == 2) {
    return 2;
  }
  else if (dir1 == 0 && dir2 == 3) {
    return 5;
  }
  else if (dir1 == 1 && dir2 == 2) {
    return 3;
  }
  else if (dir1 == 1 && dir2 == 3) {
    return 1;
  }
  else if (dir1 == 2 && dir2 == 3) {
    return 6;
  }

  cerr << "NG:Rail" << endl;
  return 1;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j].clear();
      b[i][j].clear();
    }
  }

  sx.clear();
  sy.clear();
  tx.clear();
  ty.clear();

  ans.order.clear();
  ans.score =0;
  ans.stationCount = 0;
  rep(i, T)
  {
    ans.turns[i].num = -1;
  }

  CopyToBest();

  ClearBoard();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn, tt;
    cin >> nn >> m >> k >> tt;
    sx.resize(m);
    sy.resize(m);
    tx.resize(m);
    ty.resize(m);
    rep(i, m)
    {
      cin >> sx[i] >> sy[i] >> tx[i] >> ty[i];
    }
  }
  else {
    // ファイル入力
    int nn, tt;
    ifs >> nn >> m >> k >> tt;
    sx.resize(m);
    sy.resize(m);
    tx.resize(m);
    ty.resize(m);
    rep(i, m)
    {
      ifs >> sx[i] >> sy[i] >> tx[i] >> ty[i];
    }
  }

  rep(i, m)
  {
    a[sx[i]][sy[i]].push_back(i);
    b[tx[i]][ty[i]].push_back(i);
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore()
{
  int res = ans.score;
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    rep(i, T)
    {
      if (ans.turns[i].num == -1) {
        cout << ans.turns[i].num << endl;
      }
      else {
        cout << ans.turns[i].num << ' ' << ans.turns[i].x << ' ' << ans.turns[i].y << endl;
      }
    }
  }
  else {
    // ファイル出力
    rep(i, T)
    {
      ofs << "# " << ans.comment[i] << endl;
      if (ans.turns[i].num == -1) {
        ofs << ans.turns[i].num << endl;
      }
      else {
        ofs << ans.turns[i].num << ' ' << ans.turns[i].x << ' ' << ans.turns[i].y << endl;
      }
    }
  }
}

const int MAX_STATION = 200;
struct EkiMap
{
  int ACount[MAX_M] = {};
  int BCount[MAX_M] = {};
  int Count = 0;
  int Stations[MAX_STATION][2];
  int StationCount = 50;

  EkiMap(int stationCount)
  {
    StationCount = stationCount;
  }
};

int MoveOneStation(EkiMap& ekiMap, int index, int x, int y)
{
  int beforeCount = ekiMap.Count;

  int bx = ekiMap.Stations[index][0];
  int by = ekiMap.Stations[index][1];

  rep(i, 13)
  {
    int nx = bx + ex[i];
    int ny = by + ey[i];
    if (IsNG(nx, ny))continue;
    for (auto num : a[nx][ny]) {
      ekiMap.ACount[num]--;
      if (ekiMap.ACount[num] == 0 && ekiMap.BCount[num] > 0) {
        ekiMap.Count--;
      }
    }
    for (auto num : b[nx][ny]) {
      ekiMap.BCount[num]--;
      if (ekiMap.BCount[num] == 0 && ekiMap.ACount[num] > 0) {
        ekiMap.Count--;
      }
    }
  }

  ekiMap.Stations[index][0] = x;
  ekiMap.Stations[index][1] = y;

  rep(i, 13)
  {
    int nx = x + ex[i];
    int ny = y + ey[i];
    if (IsNG(nx, ny))continue;
    for (auto num : a[nx][ny]) {
      ekiMap.ACount[num]++;
      if (ekiMap.ACount[num] == 1 && ekiMap.BCount[num] > 0) {
        ekiMap.Count++;
      }
    }
    for (auto num : b[nx][ny]) {
      ekiMap.BCount[num]++;
      if (ekiMap.BCount[num] == 1 && ekiMap.ACount[num] > 0) {
        ekiMap.Count++;
      }
    }
  }

  return ekiMap.Count - beforeCount;
}

struct Edge
{
  int from;
  int to;
  int manhattan;
  vector<P> route;

  bool operator<(const Edge& other) const
  {
    return manhattan < other.manhattan;
  }
};

vector<Edge> G[MAX_STATION];
vector<Edge> RootG[MAX_STATION];
int Parent[MAX_STATION];
Edge ParentEdge[MAX_STATION];

void MakeRootG(int root)
{
  rep(i, MAX_STATION)
  {
    RootG[i].clear();
    Parent[i] = -1;
  }

  queue<int> que;
  que.push(root);
  while (que.size()) {
    int x = que.front();
    que.pop();
    for (auto e : G[x]) {
      int y = e.to;
      if (Parent[y] != -1 || y == root)continue;
      Parent[y] = x;
      ParentEdge[y] = e;
      swap(ParentEdge[y].from, ParentEdge[y].to);
      reverse(ParentEdge[y].route.begin(), ParentEdge[y].route.end());
      RootG[x].push_back(e);
      que.push(y);
    }
  }
}

int UpdateIncome(const vector<int>& people, const vector<P>& stations)
{
  int income = 0;

  for (auto id : people) {
    bool okA = false;
    int ssx = sx[id];
    int ssy = sy[id];

    for (auto p : stations) {
      int x = p.first;
      int y = p.second;
      if (abs(ssx - x) + abs(ssy - y) <= 2) {
        okA = true;
        break;
      }
    }

    if (okA) {
      bool okB = false;
      int ttx = tx[id];
      int tty = ty[id];

      for (auto p : stations) {
        int x = p.first;
        int y = p.second;
        if (abs(ttx - x) + abs(tty - y) <= 2) {
          okB = true;
          break;
        }
      }

      if (okB) {
        income += abs(ssx - ttx) + abs(ssy - tty);
      }
    }
  }

  return income;
}

int Simulate(const EkiMap& ekiMap, const vector<int>& people)
{
  ClearBoard();
  ans.stationCount = 0;

  int money = k;
  int income = 0;
  int isBuilt[MAX_STATION] = {};

  int now = 0;
  int root = -1;
  while (turn < T) {
    //{
    //  vector<P> stations;
    //  rep(i, ekiMap.StationCount)
    //  {
    //    if (isBuilt[i]) {
    //      stations.emplace_back(ekiMap.Stations[i][0], ekiMap.Stations[i][1]);
    //    }
    //  }
    //  income = UpdateIncome(people, stations);
    //}
    if (now == 0) {
      int num = ans.order[now];
      Action(0, ekiMap.Stations[num][0], ekiMap.Stations[num][1], income);
      isBuilt[num] = 1;
      money -= 5000;
      MakeRootG(num);
      root = num;
      now++;
    }
    else if (now < ekiMap.StationCount) {
      int num = ans.order[now];
      int num2 = num;
      int preNum2 = -1;
      while (num2 != root) {
        int x = ekiMap.Stations[num2][0];
        int y = ekiMap.Stations[num2][1];
        if (board[x][y] == 0) {
          break;
        }

        if (num2 == num) {
          // 駅作る
          while (money < 5000 && turn < T) {
            Action(-1, 0, 0, income);
            money += income;
          }
          if (turn == T)break;
          Action(0, ekiMap.Stations[num2][0], ekiMap.Stations[num2][1], income);
          money -= 5000;
          isBuilt[num2] = 1;

          if (board[ParentEdge[num2].route[1].first][ParentEdge[num2].route[1].second] == -1) {
            money += income;
          }
          else {
            vector<P> stations;
            rep(i, ekiMap.StationCount)
            {
              if (isBuilt[i]) {
                stations.emplace_back(ekiMap.Stations[i][0], ekiMap.Stations[i][1]);
              }
            }
            income = UpdateIncome(people, stations);
            money += income;
            break;
          }

          // 道作る
          srep(i, 1, ParentEdge[num2].route.size() - 1)
          {
            while (money < 100 && turn < T) {
              Action(-1, 0, 0, income);
              money += income;
            }
            if (turn == T)break;
            int railNum = GetRailNum(
              ParentEdge[num2].route[i - 1].first, ParentEdge[num2].route[i - 1].second,
              ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second,
              ParentEdge[num2].route[i + 1].first, ParentEdge[num2].route[i + 1].second
            );
            Action(railNum, ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second, income);
            money -= 100;

            if (i == ParentEdge[num2].route.size() - 2) {
              int num3 = Parent[num2];
              int nx = ekiMap.Stations[num3][0];
              int ny = ekiMap.Stations[num3][1];
              if (board[nx][ny] == 0) {
                vector<P> stations;
                rep(i, ekiMap.StationCount)
                {
                  if (isBuilt[i]) {
                    stations.emplace_back(ekiMap.Stations[i][0], ekiMap.Stations[i][1]);
                  }
                }
                income = UpdateIncome(people, stations);
              }
            }

            money += income;
            if (turn == T)break;
          }
          if (turn == T)break;

          preNum2 = num2;
          num2 = Parent[num2];
        }
        else {
          if (board[x][y] == -1) {
            // 一時的に線路作る
            while (money < 100 && turn < T) {
              Action(-1, 0, 0, income);
              money += income;
            }
            if (turn == T)break;
            int railNum = GetRailNum(
              ParentEdge[preNum2].route[ParentEdge[preNum2].route.size() - 2].first, ParentEdge[preNum2].route[ParentEdge[preNum2].route.size() - 2].second,
              x, y,
              ParentEdge[num2].route[1].first, ParentEdge[num2].route[1].second
            );
            Action(railNum, x, y, income);
            money -= 100;
            money += income;

            // 道作る
            srep(i, 1, ParentEdge[num2].route.size() - 1)
            {
              while (money < 100 && turn < T) {
                Action(-1, 0, 0, income);
                money += income;
              }
              if (turn == T)break;
              int railNum = GetRailNum(
                ParentEdge[num2].route[i - 1].first, ParentEdge[num2].route[i - 1].second,
                ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second,
                ParentEdge[num2].route[i + 1].first, ParentEdge[num2].route[i + 1].second
              );
              Action(railNum, ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second, income);
              money -= 100;

              if (i == ParentEdge[num2].route.size() - 2) {
                int num3 = Parent[num2];
                int nx = ekiMap.Stations[num3][0];
                int ny = ekiMap.Stations[num3][1];
                if (board[nx][ny] == 0) {
                  vector<P> stations;
                  rep(i, ekiMap.StationCount)
                  {
                    if (isBuilt[i]) {
                      stations.emplace_back(ekiMap.Stations[i][0], ekiMap.Stations[i][1]);
                    }
                  }
                  income = UpdateIncome(people, stations);
                }
              }

              money += income;
              if (turn == T)break;
            }
            if (turn == T)break;

            preNum2 = num2;
            num2 = Parent[num2];
          }
          else {
            // 駅作る
            while (money < 5000 && turn < T) {
              Action(-1, 0, 0, income);
              money += income;
            }
            if (turn == T)break;
            Action(0, ekiMap.Stations[num2][0], ekiMap.Stations[num2][1], income);
            money -= 5000;
            isBuilt[num2] = 1;

            vector<P> stations;
            rep(i, ekiMap.StationCount)
            {
              if (isBuilt[i]) {
                stations.emplace_back(ekiMap.Stations[i][0], ekiMap.Stations[i][1]);
              }
            }
            income = UpdateIncome(people, stations);
            money += income;
            break;
          }
        }
      }
      now++;
    }
    else {
      Action(-1, 0, 0, income);
      money += income;
    }
  }

  ans.score = money;
  return money;
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition;
  int StationCount;
};

// ナイーブな解法
void Method1(Hypers hypers)
{
  //EkiMap ekiMap(hypers.StationCount);
  //EkiMap best_ekiMap(hypers.StationCount);
  EkiMap ekiMap(10 + m / 100 * 5);
  EkiMap best_ekiMap(10 + m / 100 * 5);

  // 初期解
  rep(i, ekiMap.StationCount)
  {
    int x = RandXor() % (n - 2) + 1;
    int y = RandXor() % (n - 2) + 1;
    ekiMap.Stations[i][0] = x;
    ekiMap.Stations[i][1] = y;

    rep(j, 13)
    {
      int nx = x + ex[j];
      int ny = y + ey[j];
      if (IsNG(nx, ny))continue;
      for (auto num : a[nx][ny]) {
        ekiMap.ACount[num]++;
        if (ekiMap.ACount[num] == 1 && ekiMap.BCount[num] > 0) {
          ekiMap.Count++;
        }
      }
      for (auto num : b[nx][ny]) {
        ekiMap.BCount[num]++;
        if (ekiMap.BCount[num] == 1 && ekiMap.ACount[num] > 0) {
          ekiMap.Count++;
        }
      }
    }
  }
  best_ekiMap = ekiMap;

  {
    double nowTime = GetNowTime();
    const double START_TEMP = 200.0;
    const double END_TEMP = 0.1;

    int loop = 0;
    while (true) {
      if (loop % 100 == 0) {
        nowTime = GetNowTime();
        if (nowTime > TL / 3) break;
      }
      loop++;

      int raIndex = RandXor() % ekiMap.StationCount;
      int rax = RandXor() % (n - 2) + 1;
      int ray = RandXor() % (n - 2) + 1;

      int keepx = ekiMap.Stations[raIndex][0];
      int keepy = ekiMap.Stations[raIndex][1];

      double diffScore = MoveOneStation(ekiMap, raIndex, rax, ray) * 12345.5;

      double progressRatio = nowTime / TL;
      double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
      double prob = exp(diffScore / temp);

      if (prob > Rand01()) {
        // 採用
        if (ekiMap.Count > best_ekiMap.Count) {
          best_ekiMap = ekiMap;
        }
      }
      else {
        // 元に戻す
        MoveOneStation(ekiMap, raIndex, keepx, keepy);
      }
    }

    ekiMap = best_ekiMap;

    if (mode != 0) {
      cout << "loop = " << loop << ", Count = " << ekiMap.Count << endl;
    }
  }

  vector<int> people;
  rep(i, m)
  {
    if (ekiMap.ACount[i] > 0 && ekiMap.BCount[i]) {
      people.push_back(i);
    }
  }

  // クラスカル法
  rep(i, MAX_STATION)
  {
    G[i].clear();
  }
  UnionFind uf(ekiMap.StationCount);
  vector<Edge> es;
  rep(i, ekiMap.StationCount)
  {
    srep(j, i + 1, ekiMap.StationCount)
    {
      Edge e;
      e.from = i;
      e.to = j;
      e.manhattan = max(0, abs(ekiMap.Stations[i][0] - ekiMap.Stations[j][0]) + abs(ekiMap.Stations[i][1] - ekiMap.Stations[j][1]) - 1);
      es.push_back(e);
    }
  }
  sort(es.begin(), es.end());

  int fff[n][n];
  rep(i, n)
  {
    rep(j, n)
    {
      fff[i][j] = 0;
    }
  }
  rep(i, ekiMap.StationCount)
  {
    fff[ekiMap.Stations[i][0]][ekiMap.Stations[i][1]] = 2;
  }

  for (auto e : es) {
    if (!uf.IsSame(e.from, e.to)) {
      int ssx = ekiMap.Stations[e.from][0];
      int ssy = ekiMap.Stations[e.from][1];
      int ttx = ekiMap.Stations[e.to][0];
      int tty = ekiMap.Stations[e.to][1];

      // BFSで線路つなぐ
      int f[n][n];
      int f2[n][n];
      rep(i, n)
      {
        rep(j, n)
        {
          f[i][j] = INF;
          f2[i][j] = -1;
        }
      }
      queue<P> que;
      que.push(P(ssx, ssy));
      f[ssx][ssy] = 0;
      while (que.size()) {
        int x = que.front().first;
        int y = que.front().second;
        que.pop();
        if (fff[x][y] == 2 && (x != ssx || y != ssy)) {
          int ok = 0;
          rep(i, ekiMap.StationCount)
          {
            int xx = ekiMap.Stations[i][0];
            int yy = ekiMap.Stations[i][1];
            if (x == xx && y == yy && !uf.IsSame(e.from, i)) {
              e.to = i;
              ok = 1;
              break;
            }
          }
          if (ok) {
            break;
          }
          else {
            continue;
          }
        }
        rep(i, 4)
        {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (IsNG(nx, ny))continue;
          if (fff[nx][ny] == 1)continue;
          if (f[nx][ny] > f[x][y] + 1) {
            if (fff[nx][ny] == 2) {
              f[nx][ny] = f[x][y] + 1;
              f2[nx][ny] = (i + 2) % 4;
              que.push(P(nx, ny));
            }
            else {
              f[nx][ny] = f[x][y] + 1;
              f2[nx][ny] = (i + 2) % 4;
              que.push(P(nx, ny));
            }
          }
        }
      }

      ttx = ekiMap.Stations[e.to][0];
      tty = ekiMap.Stations[e.to][1];
      if (f[ttx][tty] == INF) {
        cerr << "NG:BFS" << endl;
        continue;
        fff[ssx][ssy] = 3;
        fff[ttx][tty] = 4;
        rep(i, n)
        {
          rep(j, n)
          {
            cerr << fff[i][j];
          }
          cerr << endl;
        }
      }
      vector<P> route;
      int x = ttx;
      int y = tty;
      while (x != ssx || y != ssy) {
        route.push_back(P(x, y));
        int nx = x + dx[f2[x][y]];
        int ny = y + dy[f2[x][y]];
        x = nx;
        y = ny;
      }
      route.push_back(P(x, y));

      reverse(route.begin(), route.end());
      e.route = route;
      G[e.from].push_back(e);

      reverse(e.route.begin(), e.route.end());
      swap(e.from, e.to);
      G[e.from].push_back(e);

      uf.Unite(e.from, e.to);

      srep(i, 1, route.size() - 1)
      {
        fff[route[i].first][route[i].second] = 1;
      }
    }
  }

  //rep(i, n)
  //{
  //  rep(j, n)
  //  {
  //    cerr << fff[i][j];
  //  }
  //  cerr << endl;
  //}

  {
    ans.order.clear();
    rep(i, ekiMap.StationCount)
    {
      ans.order.push_back(i);
    }

    double nowTime = GetNowTime();
    const double START_TEMP = 2.0;
    const double END_TEMP = 0.1;

    int loop = 0;
    while (true) {
      if (loop % 100 == 0) {
        nowTime = GetNowTime();
        if (nowTime > TL) break;
      }
      loop++;

      int ra1 = RandXor() % ekiMap.StationCount;
      int ra2 = RandXor() % ekiMap.StationCount;
      while (ra1 == ra2) {
        ra2 = RandXor() % ekiMap.StationCount;
      }

      swap(ans.order[ra1], ans.order[ra2]);

      if (ans.score <= k) {
        std::shuffle(ans.order.begin(), ans.order.end(), engine);
      }

      int beforeScore = ans.score;
      int beforeStationCount = ans.stationCount;
      int afterScore = Simulate(ekiMap, people);
      int afterStationCount = ans.stationCount;
      double RATIO = 10000;

      double diffScore = ((double)afterScore + afterStationCount * RATIO - beforeScore - beforeStationCount * RATIO) * 1234.5;

      double progressRatio = nowTime / TL;
      double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
      double prob = exp(diffScore / temp);

      if (prob > Rand01()) {
        // 採用
        if (afterScore + afterStationCount * RATIO > best_ans.score + best_ans.stationCount * RATIO) {
          best_ans = ans;
        }
      }
      else {
        // 元に戻す
        swap(ans.order[ra1], ans.order[ra2]);
        ans.score = beforeScore;
        ans.stationCount =beforeStationCount;
      }
    }

    ans = best_ans;

    if (mode != 0) {
      cout << "loop = " << loop << ", Count = " << ekiMap.Count << endl;
    }
  }
}

void SimulatedAnnealing(Hypers hypers)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;


  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    int raMode = RandXor() % 100;
    if (raMode < hypers.Partition) {
      // 近傍解作成

      // スコア計算
      double tmpScore = CalcScore();

      // 焼きなまし
      double diffScore = (tmpScore - ans.score) * hypers.MultipleValue;
      double prob = exp(diffScore / temp);
      if (prob > Rand01()) {
        // 採用
        ans.score = tmpScore;

        // Best解よりもいいか
        if (ans.score > best_ans.score) {
          CopyToBest();
        }
      }
      else {
        // 元に戻す
      }
    }
    else if (raMode < 100) {

    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
  }

  CopyToAns();
}


// 問題を解く関数
ll Solve(int problem_num, Hypers hypers)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // 初期解生成
  Method1(hypers);

  // 焼きなまし
  //SimulatedAnnealing(hypers);

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

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;

  Hypers HYPERS;
  HYPERS.StartTemp = 2048.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 1.0;
  HYPERS.Partition = 50;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i, HYPERS);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "M = " << setw(4) << m << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
  }
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition = RandXor() % 101;
      hypers.StationCount = RandXor() % 151 + 10;

      ll sum = 0;
      srep(i, 25, 26)
      {
        ll score = Solve(i, hypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        //<< ", StartTemp = " << hypers.StartTemp
        //<< ", EndTemp = " << hypers.EndTemp
        //<< ", MultipleValue = " << hypers.MultipleValue
        //<< ", Partition1 = " << hypers.Partition
        << ", StationCount = " << hypers.StationCount
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
