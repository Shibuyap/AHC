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

// 移動方向の配列
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
int a[n][n];
int b[n][n];

struct Turn
{
  int num;
  int x;
  int y;
};

Turn ans[T];

int ansScore;

int best_ansScore;

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

void Action(int num, int x, int y)
{
  ans[turn].num =num;
  ans[turn].x = x;
  ans[turn].y = y;
  if (num != -1) {
    board[x][y] = num;
  }
  turn++;
}

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

bool IsNG(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  sx.clear();
  sy.clear();
  tx.clear();
  ty.clear();

  rep(i, T)
  {
    ans[i].num = -1;
  }

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

  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = -1;
      b[i][j] = -1;
    }
  }
  rep(i, m)
  {
    a[sx[i]][sy[i]] = i;
    b[tx[i]][ty[i]] = i;
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
  int res = ansScore;
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    rep(i, T)
    {
      if (ans[i].num == -1) {
        cout << ans[i].num << endl;
      }
      else {
        cout << ans[i].num << ' ' << ans[i].x << ' ' << ans[i].y << endl;
      }
    }
  }
  else {
    // ファイル出力
    rep(i, T)
    {
      if (ans[i].num == -1) {
        ofs << ans[i].num << endl;
      }
      else {
        ofs << ans[i].num << ' ' << ans[i].x << ' ' << ans[i].y << endl;
      }
    }
  }
}

struct EkiMap
{
  int ACount[MAX_M] = {};
  int BCount[MAX_M] = {};
  int Count = 0;
  int Stations[50][2];
  int StationCount = 10;
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
    if (a[nx][ny] != -1) {
      int num = a[nx][ny];
      ekiMap.ACount[num]--;
      if (ekiMap.ACount[num] == 0 && ekiMap.BCount[num] > 0) {
        ekiMap.Count--;
      }
    }
    if (b[nx][ny] != -1) {
      int num = b[nx][ny];
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
    if (a[nx][ny] != -1) {
      int num = a[nx][ny];
      ekiMap.ACount[num]++;
      if (ekiMap.ACount[num] == 1 && ekiMap.BCount[num] > 0) {
        ekiMap.Count++;
      }
    }
    if (b[nx][ny] != -1) {
      int num = b[nx][ny];
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
  int cost;

  bool operator<(const Edge& other) const
  {
    return cost < other.cost;
  }
};

// ナイーブな解法
void Method1()
{
  EkiMap ekiMap;
  EkiMap best_ekiMap;

  // 初期解
  rep(i, ekiMap.StationCount)
  {
    int x = RandXor() % n;
    int y = RandXor() % n;
    ekiMap.Stations[i][0] = x;
    ekiMap.Stations[i][1] = y;

    rep(j, 13)
    {
      int nx = x + ex[j];
      int ny = y + ey[j];
      if (IsNG(nx, ny))continue;
      if (a[nx][ny] != -1) {
        int num = a[nx][ny];
        ekiMap.ACount[num]++;
        if (ekiMap.ACount[num] == 1 && ekiMap.BCount[num] > 0) {
          ekiMap.Count++;
        }
      }
      if (b[nx][ny] != -1) {
        int num = b[nx][ny];
        ekiMap.BCount[num]++;
        if (ekiMap.BCount[num] == 1 && ekiMap.ACount[num] > 0) {
          ekiMap.Count++;
        }
      }
    }
  }
  best_ekiMap = ekiMap;

  double nowTime = GetNowTime();
  const double START_TEMP = 2.0;
  const double END_TEMP = 0.1;

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL / 3) break;
    }
    loop++;

    int raIndex = RandXor() % ekiMap.StationCount;
    int rax = RandXor() % n;
    int ray = RandXor() % n;

    int keepx = ekiMap.Stations[raIndex][0];
    int keepy = ekiMap.Stations[raIndex][1];

    double diffScore = MoveOneStation(ekiMap, raIndex, rax, ray) * 1234.5;

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

  vector<int> nums, oks;
  vector<vector<int>> as, bs;
  rep(i, m)
  {
    if (ekiMap.ACount[i] > 0 && ekiMap.BCount[i]) {
      nums.push_back(i);
      oks.push_back(0);
      vector<int> aa, bb;
      rep(j, ekiMap.StationCount)
      {
        if (abs(sx[i] - ekiMap.Stations[j][0]) + abs(sy[i] - ekiMap.Stations[j][1])) {
          aa.push_back(j);
        }
        if (abs(tx[i] - ekiMap.Stations[j][0]) + abs(ty[i] - ekiMap.Stations[j][1])) {
          bb.push_back(j);
        }
      }
      as.push_back(aa);
      bs.push_back(bb);
    }
  }

  // クラスカル法
  int money = k;
  int income = 0;
  UnionFind uf(ekiMap.StationCount);
  vector<Edge> es;
  rep(i, ekiMap.StationCount)
  {
    srep(j, i + 1, ekiMap.StationCount)
    {
      Edge e;
      e.from = i;
      e.to = j;
      e.cost = max(0, abs(ekiMap.Stations[i][0] - ekiMap.Stations[j][0]) + abs(ekiMap.Stations[i][1] - ekiMap.Stations[j][1]) - 1);
      es.push_back(e);
    }
  }
  sort(es.begin(), es.end());
  for (auto e : es) {
    if (!uf.IsSame(e.from, e.to)) {

      int ssx = ekiMap.Stations[e.from][0];
      int ssy = ekiMap.Stations[e.from][1];
      int ttx = ekiMap.Stations[e.to][0];
      int tty = ekiMap.Stations[e.to][1];

      // BFSで線路つなぐ
      {
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
          if (x == ttx && y == tty)break;
          rep(i, 4)
          {
            int nx = x + dx[i];
            int ny = y + dy[i];
            if (IsNG(nx, ny))continue;
            if (f[nx][ny] > f[x][y] + 1) {
              f[nx][ny] = f[x][y] + 1;
              f2[nx][ny] = (i + 2) % 4;
              que.push(P(nx, ny));
            }
          }
        }
        if (f[ttx][tty] == INF) {
          cerr << "NG:BFS" << endl;
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

        srep(i, 1, route.size() - 1)
        {

        }
      }

      if (board[ssx][ssy] != 0) {
        while (money < 5000 && turn < T) {
          Action(-1, 0, 0);
          money +=income;
        }
        if (turn >= T) break;
        Action(0, ssx, ssy);
        money -= 5000;
        money += income;
      }

      if (board[ttx][tty] != 0) {
        while (money < 5000 && turn < T) {
          Action(-1, 0, 0);
          money +=income;
        }
        if (turn >= T) break;
        Action(0, ttx, tty);
        money -= 5000;
        money += income;
      }

      uf.Unite(e.from, e.to);
      rep(i, nums.size())
      {
        if (oks[i])continue;
        rep(j, as[i].size())
        {
          rep(k, bs[i].size())
          {
            if (uf.IsSame(as[i][j], bs[i][k])) {
              oks[i] = 1;
              income += abs(sx[nums[i]] - tx[nums[i]]) + abs(sy[nums[i]] - ty[nums[i]]);
              break;
            }
          }
          if (oks[i])break;
        }
      }
    }
    if (turn >= T) break;
  }

  ansScore = money;
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition;
};

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
      double diffScore = (tmpScore - ansScore) * hypers.MultipleValue;
      double prob = exp(diffScore / temp);
      if (prob > Rand01()) {
        // 採用
        ansScore = tmpScore;

        // Best解よりもいいか
        if (ansScore > best_ansScore) {
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
  Method1();

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

      ll sum = 0;
      srep(i, 0, 15)
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
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition
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
