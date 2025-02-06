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
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
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
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const int INF = 1001001;

std::mt19937 engine;
double TL = 1.8;
int mode;

/////////////////////////////////////////////////////////////////////////////////
// 変数

struct Crop
{
  int id;
  int s;
  int d;
};

bool compD(const Crop& a, const Crop& b) { return a.d < b.d; }
bool compS(const Crop& a, const Crop& b) { return a.s < b.s; }

// 入力
const int T = 100;
const int H = 20;
const int W = 20;
const int HW = 400;
const int HWT = 40000;
const int SY = 0;
int SX;
int ho[25][25];
int ve[25][25];
int K;
int S[44000], D[44000];
int wall[25][25][4];
vector<Crop> SVec[100], DVec[100];
int Used[44000];

// 出力
int M;
int ansK[44000], ansX[44000], ansY[44000], ansS[44000];

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
        if (!used[i]) k = dfs(i, k, -1);
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
          if (par != -1 && ord[id] <= low[e.to])
            is_aps = true;  // 条件2を満たすので関節点
          if (ord[id] < low[e.to])
            bridges.emplace_back(min(id, e.to), max(id, e.to));
        }
        else if (e.to != par) {  // eが後退辺の時
          low[id] = min(low[id], ord[e.to]);
        }
      }
      if (par == -1 && count >= 2) is_aps = true;  // 条件1を満たすので関節点
      if (is_aps) aps.push_back(id);
      return k;
    }
  };
}  // namespace

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  M = 0;
  rep(i, 100)
  {
    SVec[i].clear();
    DVec[i].clear();
  }
  rep(i, 44000)
  {
    Used[i] = 0;
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
    int TTT, HHH, WWW;
    cin >> TTT >> HHH >> WWW >> SX;
    rep(i, H - 1)
    {
      string str;
      cin >> str;
      rep(j, W) { ho[i][j] = str[j] - '0'; }
    }
    rep(i, H)
    {
      string str;
      cin >> str;
      rep(j, W - 1) { ve[i][j] = str[j] - '0'; }
    }
    cin >> K;
    rep(i, K) { cin >> S[i] >> D[i]; }
  }
  // ファイル入力する
  else {
    int TTT, HHH, WWW;
    ifs >> TTT >> HHH >> WWW >> SX;
    rep(i, H - 1)
    {
      string str;
      ifs >> str;
      rep(j, W) { ho[i][j] = str[j] - '0'; }
    }
    rep(i, H)
    {
      string str;
      ifs >> str;
      rep(j, W - 1) { ve[i][j] = str[j] - '0'; }
    }
    ifs >> K;
    rep(i, K) { ifs >> S[i] >> D[i]; }
  }

  rep(i, H)
  {
    rep(j, W)
    {
      rep(k, 4) { wall[i][j][k] = 0; }
      if (i == 0) {
        wall[i][j][0] = 1;
      }
      else {
        if (ho[i - 1][j]) {
          wall[i][j][0] = 1;
        }
      }
      if (j == 0) {
        wall[i][j][1] = 1;
      }
      else {
        if (ve[i][j - 1]) {
          wall[i][j][1] = 1;
        }
      }
      if (i == H - 1) {
        wall[i][j][2] = 1;
      }
      else {
        if (ho[i][j]) {
          wall[i][j][2] = 1;
        }
      }
      if (j == W - 1) {
        wall[i][j][3] = 1;
      }
      else {
        if (ve[i][j]) {
          wall[i][j][3] = 1;
        }
      }
    }
  }

  rep(i, K)
  {
    S[i]--;
    D[i]--;
    Crop crop;
    crop.id = i + 1;
    crop.s = S[i];
    crop.d = D[i];
    rep(l, 1)
    {
      if (S[i] - l < 0)break;
      SVec[S[i] - l].push_back(crop);
    }
    DVec[D[i]].push_back(crop);
  }

  rep(i, T)
  {
    sort(SVec[i].begin(), SVec[i].end(), compD);
    sort(DVec[i].begin(), DVec[i].end(), compS);
  }

  // rep(i, T) {
  //   cout << i << " : ";
  //   for (auto x : SVec[i]) {
  //     cout << x.d << ' ';
  //   }
  //   cout << endl;
  // }
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
ll CalcScore()
{
  ll sum = 0;
  rep(i, M) { sum += D[ansK[i] - 1] - S[ansK[i] - 1] + 1; }
  ll res = 1000000LL * sum / HWT;
  return res;
}

// 初期解生成
void Initialize() {}

bool OKCheck(const int x, const int y, const vector<vector<int>>& use)
{
  if (x == SX && y == SY) {
    return false;
  }
  int blankCount = 0;
  rep(i, 4)
  {
    if (wall[x][y][i]) continue;
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

int bfsQueue[10000][2];
bool NGCheck_BFS(const vector<vector<int>>& use)
{
  int f[H][W];
  rep(i, H)
  {
    rep(j, W) { f[i][j] = 0; }
  }
  int cnt = 0;
  rep(i, H)
  {
    rep(j, W)
    {
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
    rep(i, 4)
    {
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

int dijkstraQueue[110][410][2];
bool NGCheck_Dijkstra(const vector<vector<int>>& use)
{
  int f[H][W];
  rep(i, H)
  {
    rep(j, W) { f[i][j] = INF; }
  }

  int cnt = 0;

  int heads[110] = {};
  int tails[110] = {};
  f[SX][SY] = use[SX][SY];
  cnt++;
  int hash = f[SX][SY] + 1;
  dijkstraQueue[hash][tails[hash]][0] = SX;
  dijkstraQueue[hash][tails[hash]][1] = SY;
  tails[hash]++;

  rep(turn, T)
  {
    while (heads[turn] < tails[turn]) {
      int x = dijkstraQueue[turn][heads[turn]][0];
      int y = dijkstraQueue[turn][heads[turn]][1];
      int val = turn - 1;
      heads[turn]++;
      rep(i, 4)
      {
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
  rep(i, H)
  {
    rep(j, W) { f[i][j] = INF; }
  }


  int heads[110] = {};
  int tails[110] = {};
  f[sx][sy] = use[sx][sy];

  int hash = f[sx][sy] + 1;
  dijkstraQueue[hash][tails[hash]][0] = sx;
  dijkstraQueue[hash][tails[hash]][1] = sy;
  tails[hash]++;

  vector<int> days;
  int kind[110] = {};
  rep(turn, T)
  {
    while (heads[turn] < tails[turn]) {
      kind[turn] = 1;
      int x = dijkstraQueue[turn][heads[turn]][0];
      int y = dijkstraQueue[turn][heads[turn]][1];
      int val = turn - 1;
      heads[turn]++;
      rep(i, 4)
      {
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
  rep(i, days.size())
  {
    int day = days[i];
    if (d == day) {
      score += 1e5 / ((i + 1) * (i + 1));
    }
    else {
      score -= (double)abs(d - day) / ((i + 1) * (i + 1));
    }
  }
  rep(i, 110)
  {
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

  rep(i, H)
  {
    rep(j, W)
    {
      if (use[i][j] != -1) {
        rep(k, 4)
        {
          if (wall[i][j][k])continue;
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


void Method1()
{
  vector<vector<int>> use(20, vector<int>(20));
  rep(i, H)
  {
    rep(j, W) { use[i][j] = -1; }
  }

  rep(turn, T)
  {
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
      rep(i, HW)
      {
        int x = i / H;
        int y = i % H;
        if (use[x][y] == -1) {
          int num = mp.size();
          mp[num] = P(x, y);
          mpInv[P(x, y)] = num;
          Graph.push_back(vector<LowLinkEdge>());
        }
      }
      // エッジ
      rep(num, mp.size())
      {
        int x = mp[num].first;
        int y = mp[num].second;
        rep(j, 4)
        {
          if (wall[x][y][j]) continue;
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
      rep(i, mp.size())
      {
        if (aps.find(i) != aps.end()) {
          continue;
        }
        blanks.push_back(mp[i]);
      }
      std::shuffle(blanks.begin(), blanks.end(), engine);

      int walkCount[H][W];
      rep(j, H)
      {
        rep(k, W)
        {
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
        rep(j, 4)
        {
          if (wall[x][y][j])continue;
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
        if (OKCheck(x, y, use)) {
          OK = true;
        }
        if (!OK) {
          if (NGCheck_BFS(use) && NGCheck_Dijkstra(use)) {
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
    rep(i, H)
    {
      rep(j, W)
      {
        if (use[i][j] == turn) {
          use[i][j] = -1;
        }
      }
    }
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << M << endl;
    rep(i, M)
    {
      cout << ansK[i] << ' ' << ansX[i] << ' ' << ansY[i] << ' ' << ansS[i] + 1
        << endl;
    }
  }
  else {
    ofs << M << endl;
    rep(i, M)
    {
      ofs << ansK[i] << ' ' << ansX[i] << ' ' << ansY[i] << ' ' << ansS[i] + 1
        << endl;
    }
  }
}

ll Solve(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();

  // 貪欲1
  Method1();

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
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }
  std::random_device rnd;
  engine.seed(rnd());

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 20)
    {
      rep(j, 1)
      {
        ll score = Solve(i);
        sum += score;
        cout << "num = " << i << ", ";
        cout << "score = " << score << ", ";
        cout << "sum = " << sum << endl;
      }
    }
  }

  return 0;
}
