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
      Par[i] = i;
      Rank[i] = 0;
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
    if (x == y) { return; }

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


static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
}


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const int INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int station_dx[13] = { -2,-1,-1,-1,0,0,0,0,0,1,1,1,2 };
const int station_dy[13] = { 0,-1,0,1,-2,-1,0,1,2,-1,0,1,0 };

double TL = 2.8;
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
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
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

bool OutOfBounds(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

int CalcRailType(int x1, int y1, int x2, int y2, int x3, int y3)
{
  int dir1 = 0, dir2 = 0;
  for (int i = 0; i < 4; ++i) {
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
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i][j].clear();
      b[i][j].clear();
    }
  }

  sx.clear();
  sy.clear();
  tx.clear();
  ty.clear();

  ans.order.clear();
  ans.score = 0;
  ans.stationCount = 0;
  for (int i = 0; i < T; ++i) {
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
    for (int i = 0; i < m; ++i) {
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
    for (int i = 0; i < m; ++i) {
      ifs >> sx[i] >> sy[i] >> tx[i] >> ty[i];
    }
  }

  for (int i = 0; i < m; ++i) {
    a[sx[i]][sy[i]].push_back(i);
    b[tx[i]][ty[i]].push_back(i);
  }
}

// 出力ファイルストリームを開く関数
void OpenOutput(int probNum, ofstream& ofs)
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
    for (int i = 0; i < T; ++i) {
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
    for (int i = 0; i < T; ++i) {
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

int MoveOneStation(EkiMap& emap, int index, int x, int y)
{
  int beforeCount = emap.Count;

  int bx = emap.Stations[index][0];
  int by = emap.Stations[index][1];

  for (int i = 0; i < 13; ++i) {
    int nx = bx + station_dx[i];
    int ny = by + station_dy[i];
    if (OutOfBounds(nx, ny))continue;
    for (auto num : a[nx][ny]) {
      emap.ACount[num]--;
      if (emap.ACount[num] == 0 && emap.BCount[num] > 0) {
        emap.Count -= abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
    for (auto num : b[nx][ny]) {
      emap.BCount[num]--;
      if (emap.BCount[num] == 0 && emap.ACount[num] > 0) {
        emap.Count -= abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
  }

  emap.Stations[index][0] = x;
  emap.Stations[index][1] = y;

  for (int i = 0; i < 13; ++i) {
    int nx = x + station_dx[i];
    int ny = y + station_dy[i];
    if (OutOfBounds(nx, ny))continue;
    for (auto num : a[nx][ny]) {
      emap.ACount[num]++;
      if (emap.ACount[num] == 1 && emap.BCount[num] > 0) {
        emap.Count += abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
    for (auto num : b[nx][ny]) {
      emap.BCount[num]++;
      if (emap.BCount[num] == 1 && emap.ACount[num] > 0) {
        emap.Count += abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
  }

  return emap.Count - beforeCount;
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
  for (int i = 0; i < MAX_STATION; ++i) {
    RootG[i].clear();
    Parent[i] = -1;
  }

  queue<int> q;
  q.push(root);
  while (q.size()) {
    int x = q.front();
    q.pop();
    for (auto e : G[x]) {
      int y = e.to;
      if (Parent[y] != -1 || y == root)continue;
      Parent[y] = x;
      ParentEdge[y] = e;
      swap(ParentEdge[y].from, ParentEdge[y].to);
      reverse(ParentEdge[y].route.begin(), ParentEdge[y].route.end());
      RootG[x].push_back(e);
      q.push(y);
    }
  }
}

int okA[MAX_M];
int okB[MAX_M];
vector<vector<int>> numsA;
vector<vector<int>> numsB;
int UpdateIncome(const vector<int>& people, const vector<int>& stations, int income)
{
  for (auto st : stations) {
    for (auto num : numsA[st]) {
      okA[num]++;
      if (okA[num] == 1 && okB[num] > 0) {
        income += abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
    for (auto num : numsB[st]) {
      okB[num]++;
      if (okB[num] == 1 && okA[num] > 0) {
        income += abs(sx[num] - tx[num]) + abs(sy[num] - ty[num]);
      }
    }
  }

  return income;
}

int Simulate(int simMode, const EkiMap& emap, const vector<int>& people, bool isSkip = true)
{
  int TT = T;
  if (simMode != 0) {
    TT = 100 * simMode;
  }

  ClearBoard();
  ans.stationCount = 0;
  for (auto id : people) {
    okA[id] = 0;
    okB[id] = 0;
  }

  int money = k;
  int income = 0;
  int isBuilt[MAX_STATION] = {};

  int now = 0;
  int root = -1;
  vector<int> stations;
  while (turn < TT) {
    if (now == 0) {
      int num = ans.order[now];
      Action(0, emap.Stations[num][0], emap.Stations[num][1], income);
      isBuilt[num] = 1;
      stations.push_back(num);
      money -= 5000;
      MakeRootG(num);
      root = num;
      now++;
    }
    else if (now < emap.StationCount) {
      int num = ans.order[now];
      int num2 = num;
      int preNum2 = -1;
      while (num2 != root) {
        int x = emap.Stations[num2][0];
        int y = emap.Stations[num2][1];
        if (board[x][y] == 0) {
          break;
        }

        if (num2 == num) {
          // 駅作る
          while (money < 5000 && turn < TT) {
            Action(-1, 0, 0, income);
            money += income;
            if (isSkip) {
              if (income == 0) {
                turn = TT;
              }
            }
          }
          if (turn == TT)break;
          Action(0, emap.Stations[num2][0], emap.Stations[num2][1], income);
          money -= 5000;
          isBuilt[num2] = 1;
          stations.push_back(num2);

          if (board[ParentEdge[num2].route[1].first][ParentEdge[num2].route[1].second] == -1) {
            money += income;
          }
          else {
            income = UpdateIncome(people, stations, income);
            stations.clear();
            money += income;
            break;
          }

          // 道作る
          for (int i = 1; i < ParentEdge[num2].route.size() - 1; ++i) {
            while (money < 100 && turn < TT) {
              Action(-1, 0, 0, income);
              money += income;
              if (isSkip) {
                if (income == 0) {
                  turn = TT;
                }
              }
            }
            if (turn == TT)break;
            int railNum = CalcRailType(
              ParentEdge[num2].route[i - 1].first, ParentEdge[num2].route[i - 1].second,
              ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second,
              ParentEdge[num2].route[i + 1].first, ParentEdge[num2].route[i + 1].second
            );
            Action(railNum, ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second, income);
            money -= 100;

            if (i == ParentEdge[num2].route.size() - 2) {
              int num3 = Parent[num2];
              int nx = emap.Stations[num3][0];
              int ny = emap.Stations[num3][1];
              if (board[nx][ny] == 0) {
                income = UpdateIncome(people, stations, income);
                stations.clear();
              }
            }

            money += income;
            if (turn == TT)break;
          }
          if (turn == TT)break;

          preNum2 = num2;
          num2 = Parent[num2];
        }
        else {
          if (board[x][y] == -1) {
            // 一時的に線路作る
            while (money < 100 && turn < TT) {
              Action(-1, 0, 0, income);
              money += income;
              if (isSkip) {
                if (income == 0) {
                  turn = TT;
                }
              }
            }
            if (turn == TT)break;
            int railNum = CalcRailType(
              ParentEdge[preNum2].route[ParentEdge[preNum2].route.size() - 2].first, ParentEdge[preNum2].route[ParentEdge[preNum2].route.size() - 2].second,
              x, y,
              ParentEdge[num2].route[1].first, ParentEdge[num2].route[1].second
            );
            Action(railNum, x, y, income);
            money -= 100;
            money += income;

            // 道作る
            for (int i = 1; i < ParentEdge[num2].route.size() - 1; ++i) {
              while (money < 100 && turn < TT) {
                Action(-1, 0, 0, income);
                money += income;
                if (isSkip) {
                  if (income == 0) {
                    turn = TT;
                  }
                }
              }
              if (turn == TT)break;
              int railNum = CalcRailType(
                ParentEdge[num2].route[i - 1].first, ParentEdge[num2].route[i - 1].second,
                ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second,
                ParentEdge[num2].route[i + 1].first, ParentEdge[num2].route[i + 1].second
              );
              Action(railNum, ParentEdge[num2].route[i].first, ParentEdge[num2].route[i].second, income);
              money -= 100;

              if (i == ParentEdge[num2].route.size() - 2) {
                int num3 = Parent[num2];
                int nx = emap.Stations[num3][0];
                int ny = emap.Stations[num3][1];
                if (board[nx][ny] == 0) {
                  income = UpdateIncome(people, stations, income);
                  stations.clear();
                }
              }

              money += income;
              if (turn == TT)break;
            }
            if (turn == TT)break;

            preNum2 = num2;
            num2 = Parent[num2];
          }
          else {
            // 駅作る
            while (money < 5000 && turn < TT) {
              Action(-1, 0, 0, income);
              money += income;
              if (isSkip) {
                if (income == 0) {
                  turn = TT;
                }
              }
            }
            if (turn == TT)break;
            Action(0, emap.Stations[num2][0], emap.Stations[num2][1], income);
            money -= 5000;
            isBuilt[num2] = 1;
            stations.push_back(num2);

            income = UpdateIncome(people, stations, income);
            stations.clear();
            money += income;
            break;
          }
        }
      }
      now++;
    }
    else {
      if (isSkip) {
        money += income * (TT - turn);
        turn = TT;
      }
      else {
        Action(-1, 0, 0, income);
        money += income;
      }
    }
  }

  money += (T - TT) * income;

  ans.score = money;
  return money;
}

// ハイパーパラメータ
struct Params
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition;
  int StationCount;
};

// ナイーブな解法
void StationPlacement(Params params)
{
  int stationCount = 40 + m / 100 * 5;
  if (mode == 3) {
    stationCount = params.StationCount;
  }
  EkiMap emap(stationCount);
  EkiMap best_emap(stationCount);

  // 初期解
  for (int i = 0; i < emap.StationCount; ++i) {
    int x = Rand() % n;
    int y = Rand() % n;
    emap.Stations[i][0] = x;
    emap.Stations[i][1] = y;

    for (int j = 0; j < 13; ++j) {
      int nx = x + station_dx[j];
      int ny = y + station_dy[j];
      if (OutOfBounds(nx, ny))continue;
      for (auto num : a[nx][ny]) {
        emap.ACount[num]++;
        if (emap.ACount[num] == 1 && emap.BCount[num] > 0) {
          emap.Count++;
        }
      }
      for (auto num : b[nx][ny]) {
        emap.BCount[num]++;
        if (emap.BCount[num] == 1 && emap.ACount[num] > 0) {
          emap.Count++;
        }
      }
    }
  }
  best_emap = emap;

  {
    double nowTime = GetNowTime();
    const double START_TEMP = 200.0;
    const double END_TEMP = 0.1;

    int loop = 0;
    while (true) {
      if (loop % 100 == 0) {
        nowTime = GetNowTime();
        if (nowTime > TL / 3) { break; }
      }
      loop++;

      int r_idx = Rand() % emap.StationCount;

      int r_x = Rand() % n;
      int r_y = Rand() % n;

      int keep_x = emap.Stations[r_idx][0];
      int keep_y = emap.Stations[r_idx][1];

      double diffScore = MoveOneStation(emap, r_idx, r_x, r_y) * 12345.5;

      double progressRatio = nowTime / TL;
      double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
      double prob = exp(diffScore / temp);

      if (prob > Rand01()) {
        // 採用
        if (emap.Count > best_emap.Count) {
          best_emap = emap;
        }
      }
      else {
        // 元に戻す
        MoveOneStation(emap, r_idx, keep_x, keep_y);
      }
    }

    emap = best_emap;

    if (mode != 0) {
      cout << "loop = " << loop << ", Count = " << emap.Count << endl;
    }
  }

  vector<int> people;
  numsA.resize(emap.StationCount);
  numsB.resize(emap.StationCount);
  for (int i = 0; i < emap.StationCount; ++i) {
    numsA[i].clear();
    numsB[i].clear();
  }
  for (int i = 0; i < m; ++i) {
    if (emap.ACount[i] > 0 && emap.BCount[i] > 0) {
      people.push_back(i);
      for (int j = 0; j < emap.StationCount; ++j) {
        int x = emap.Stations[j][0];
        int y = emap.Stations[j][1];
        if (abs(x - sx[i]) + abs(y - sy[i]) <= 2) {
          numsA[j].push_back(i);
        }
        if (abs(x - tx[i]) + abs(y - ty[i]) <= 2) {
          numsB[j].push_back(i);
        }
      }
    }
  }

  // クラスカル法
  for (int i = 0; i < MAX_STATION; ++i) {
    G[i].clear();
  }
  UnionFind uf(emap.StationCount);
  vector<Edge> es;
  for (int i = 0; i < emap.StationCount; ++i) {
    for (int j = i + 1; j < emap.StationCount; ++j) {
      Edge e;
      e.from = i;
      e.to = j;
      e.manhattan = max(0, abs(emap.Stations[i][0] - emap.Stations[j][0]) + abs(emap.Stations[i][1] - emap.Stations[j][1]) - 1);
      es.push_back(e);
    }
  }
  sort(es.begin(), es.end());

  int rail_map[n][n];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      rail_map[i][j] = 0;
    }
  }
  for (int i = 0; i < emap.StationCount; ++i) {
    rail_map[emap.Stations[i][0]][emap.Stations[i][1]] = 2;
  }

  for (auto e : es) {
    if (!uf.IsSame(e.from, e.to)) {
      int ssx = emap.Stations[e.from][0];
      int ssy = emap.Stations[e.from][1];
      int ttx = emap.Stations[e.to][0];
      int tty = emap.Stations[e.to][1];

      // BFSで線路つなぐ
      int f[n][n];
      int prev_dir[n][n];
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          f[i][j] = INF;
          prev_dir[i][j] = -1;
        }
      }
      queue<P> q;
      q.push(P(ssx, ssy));
      f[ssx][ssy] = 0;
      while (q.size()) {
        int x = q.front().first;
        int y = q.front().second;
        q.pop();
        if (rail_map[x][y] == 2 && (x != ssx || y != ssy)) {
          int ok = 0;
          for (int i = 0; i < emap.StationCount; ++i) {
            int xx = emap.Stations[i][0];
            int yy = emap.Stations[i][1];
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
        for (int i = 0; i < 4; ++i) {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (OutOfBounds(nx, ny))continue;
          if (rail_map[nx][ny] == 1)continue;
          if (f[nx][ny] > f[x][y] + 1) {
            if (rail_map[nx][ny] == 2) {
              f[nx][ny] = f[x][y] + 1;
              prev_dir[nx][ny] = (i + 2) % 4;
              q.push(P(nx, ny));
            }
            else {
              f[nx][ny] = f[x][y] + 1;
              prev_dir[nx][ny] = (i + 2) % 4;
              q.push(P(nx, ny));
            }
          }
        }
      }

      ttx = emap.Stations[e.to][0];
      tty = emap.Stations[e.to][1];
      if (f[ttx][tty] == INF) {
        cerr << "NG:BFS" << endl;
        continue;
        rail_map[ssx][ssy] = 3;
        rail_map[ttx][tty] = 4;
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            cerr << rail_map[i][j];
          }
          cerr << endl;
        }
      }
      vector<P> route;
      int x = ttx;
      int y = tty;
      while (x != ssx || y != ssy) {
        route.push_back(P(x, y));
        int nx = x + dx[prev_dir[x][y]];
        int ny = y + dy[prev_dir[x][y]];
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

      for (int i = 1; i < route.size() - 1; ++i) {
        rail_map[route[i].first][route[i].second] = 1;
      }
    }
  }

  {
    ans.order.clear();
    for (int i = 0; i < emap.StationCount; ++i) {
      ans.order.push_back(i);
    }

    double startTime = GetNowTime();
    double nowTime = GetNowTime();
    const double START_TEMP = 2.0;
    const double END_TEMP = 0.1;

    int loop = 0;
    while (true) {
      if (loop % 100 == 0) {
        nowTime = GetNowTime();
        if (nowTime > TL) { break; }
      }
      loop++;
      double progressRatio = (nowTime - startTime) / (TL - startTime);

      int r1 = Rand() % emap.StationCount;
      int r2 = Rand() % emap.StationCount;
      while (r1 == r2) {
        r2 = Rand() % emap.StationCount;
      }

      swap(ans.order[r1], ans.order[r2]);

      if (ans.score <= k) {
        std::shuffle(ans.order.begin(), ans.order.end(), engine);
      }

      int simMode = 0;
      if (progressRatio < 0.5 && Rand() % 2 == 0) {
        simMode = progressRatio * 10 + 1;
      }

      int beforeScore = ans.score;
      int beforeStationCount = ans.stationCount;
      int afterScore = Simulate(simMode, emap, people);
      int afterStationCount = ans.stationCount;
      double RATIO = 0;

      double diffScore = ((double)afterScore + afterStationCount * RATIO - beforeScore - beforeStationCount * RATIO) * 1234.5;


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
        swap(ans.order[r1], ans.order[r2]);
        ans.score = beforeScore;
        ans.stationCount = beforeStationCount;
      }
    }

    ans = best_ans;
    Simulate(0, emap, people, false);

    if (mode != 0) {
      cout << "loop = " << loop << ", Count = " << emap.Count << endl;
    }
  }
}

// 貪欲
void GreedyPlacement()
{
  ClearBoard();

  // 2点決める

}

void SA(Params params)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = params.StartTemp;
  const double END_TEMP = params.EndTemp;


  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) { break; }
    }

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    int r_mode = Rand() % 100;
    if (r_mode < params.Partition) {
      // 近傍解作成

      // スコア計算
      double tmpScore = CalcScore();

      // 焼きなまし
      double diffScore = (tmpScore - ans.score) * params.MultipleValue;
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
    else if (r_mode < 100) {

    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
  }

  CopyToAns();
}


// 問題を解く関数
ll Solve(int problem_num, Params params)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOutput(problem_num, ofs);

  // 初期解生成
  StationPlacement(params);
  //GreedyPlacement();

  // 焼きなまし
  //SA(params);

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

int main()
{
  mode = 2;

  Params PARAMS;
  PARAMS.StartTemp = 2048.0;
  PARAMS.EndTemp = 0.0;
  PARAMS.MultipleValue = 1.0;
  PARAMS.Partition = 50;

  if (mode == 0) {
    Solve(0, PARAMS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(i, PARAMS);
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
    Params bestParams;
    ll bestSumScore = 0;

    while (true) {
      Params params;
      params.StartTemp = pow(2.0, Rand01() * 20);
      params.EndTemp = 0.0;
      params.MultipleValue = pow(2.0, Rand01() * 20);
      params.Partition = Rand() % 101;
      params.StationCount = Rand() % 1 + 30;
      cout << "StationCount = " << params.StationCount << endl;

      ll sum = 0;
      for (int i = 0; i < 10; ++i) {
        ll score = Solve(i, params);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StationCount = " << params.StationCount
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestParams = params;
      }

      loop++;
      break;
    }
  }

  return 0;
}
