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

int run_mode;

namespace /* 乱数ライブラリ */
{
  static uint32_t rand_u32()
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

  static double rand_unit() {
    return (rand_u32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

const ll INF = 1001001001001001001;
struct edge
{
  int u, v, cost;
  int id;
};

const int MAX_N = 100;
const int MAX_M = 300;
const int MAX_K = 5000;

// 定数
vector<edge> edges;
int node_cnt, edge_cnt, res_cnt;
ll node_x[MAX_N], node_y[MAX_N];
ll u[MAX_M], v[MAX_M], w[MAX_M];
ll res_x[MAX_K], res_y[MAX_K];
ll distn[MAX_N][MAX_K];
ll distk[MAX_K][MAX_N];
ll ndist[MAX_N][MAX_N];
ll ncost[MAX_N][MAX_N];
vector<P> nears[MAX_K];
vector<P> nearsn[MAX_N];

// 更新する変数
ll maxScore;
ll real_maxScore;
ll p[MAX_N];
ll B[MAX_M];
ll real_p[MAX_N];
ll real_B[MAX_M];
int maxPowers[MAX_N] = {};
int real_maxPowers[MAX_N] = {};

// Union-Find
namespace
{
  void init(vector<int>& par, vector<int>& rank, int n)
  {
    for (int i = 0; i < n; i++) {
      par[i] = i;
      rank[i] = 0;
    }
  }

  // 木の根を求める
  int find(vector<int>& par, int x)
  {
    if (par[x] == x) {
      return x;
    }
    else {
      return par[x] = find(par, par[x]);
    }
  }

  // xとyの属する集合を併合
  void unite(vector<int>& par, vector<int>& rank, int x, int y)
  {
    x = find(par, x);
    y = find(par, y);
    if (x == y) return;

    if (rank[x] < rank[y]) {
      par[x] = y;
    }
    else {
      par[y] = x;
      if (rank[x] == rank[y]) rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool same(vector<int>& par, int x, int y)
  {
    return find(par, x) == find(par, y);
  }
}  // namespace

bool comp(const edge& e1, const edge& e2) { return e1.cost < e2.cost; }

int isNeed[MAX_N];
int real_isNeed[MAX_N];
int isUse[MAX_K];
long long int kruscal(int V)
{
  int needCount = 0;
  rep(i, V) { needCount += isNeed[i]; }
  int E = edges.size();
  rep(i, E) isUse[i] = 0;

  vector<int> par(V + 10);   // 親
  vector<int> rank(V + 10);  // 木の深さ
  init(par, rank, V + 10);   // Union-Findの初期化
  long long int res = 0;
  int uniteCount = 0;
  for (int i = 0; i < E; i++) {
    edge e = edges[i];
    if (!isNeed[e.u] || !isNeed[e.v]) {
      continue;
    }
    if (!same(par, e.u, e.v)) {
      isUse[e.id] = 1;
      unite(par, rank, e.u, e.v);
      res += e.cost;
      uniteCount++;
    }
  }

  if (uniteCount != needCount - 1) res = INF;
  return res;
}

ll outer_kruskal() { return kruscal(node_cnt); }

ll outer_kruskal_lns()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  rep(i, node_cnt) { isNeed[i] = 1; }
  ll mi = outer_kruskal();
  int loop = 0;
  double TL = 0.1;
  while (true) {
    loop++;
    {
      endTime = clock();
      double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      double nowProgress = nowTime / TL;
      if (nowProgress > 1.0) break;
    }

    int num = rand_u32() % (node_cnt - 1) + 1;
    if (p[num] != 0) {
      continue;
    }

    isNeed[num] = 1 - isNeed[num];
    ll tmp = outer_kruskal();
    if (tmp <= mi) {
      mi = tmp;
    }
    else {
      isNeed[num] = 1 - isNeed[num];
    }
  }

  return mi;
}

bool check_coverage()
{
  rep(i, res_cnt)
  {
    int ok = 0;
    int sz = nears[i].size();
    rep(j, sz)
    {
      int jj = nears[i][j].second;
      if (distk[i][jj] <= p[jj]) {
        ok = 1;
        break;
      }
    }
    if (!ok) {
      return false;
    }
  }
  return true;
}

bool check_coverage_node(int nn)
{
  for (auto ii : nearsn[nn]) {
    int i = ii.second;
    int ok = 0;
    int sz = nears[i].size();
    rep(j, sz)
    {
      int jj = nears[i][j].second;
      if (distk[i][jj] <= p[jj]) {
        ok = 1;
        break;
      }
    }
    if (!ok) {
      return false;
    }
  }
  return true;
}

ll keepW;
ll calc_score(bool isALL = true)
{
  if (!check_coverage()) {
    return -1;
  }

  double S = 0;
  rep(i, node_cnt) { S += p[i] * p[i]; }
  if (isALL) {
    keepW = 0;
    rep(i, edge_cnt)
    {
      if (res_y[i]) {
        S += w[i];
        keepW += w[i];
      }
    }
  }
  else {
    S += keepW;
  }


  ll point = round(1e6 * (1.0 + 1e8 / (S + 1e7)));
  return point;
}

ll calc_score_partial(int nn, bool isALL = true)
{
  if (!check_coverage_node(nn)) {
    return -1;
  }

  double S = 0;
  rep(i, node_cnt) { S += p[i] * p[i]; }
  if (isALL) {
    keepW = 0;
    rep(i, edge_cnt)
    {
      if (res_y[i]) {
        S += w[i];
        keepW += w[i];
      }
    }
  }
  else {
    S += keepW;
  }


  ll point = round(1e6 * (1.0 + 1e8 / (S + 1e7)));
  return point;
}

ll sq_dist(ll x1, ll y1, ll x2, ll y2)
{
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

ll min_radius_le_5000(ll x1, ll y1, ll x2, ll y2)
{
  if (sq_dist(x1, y1, x2, y2) > 25000000) {
    return INF;
  }
  if (sq_dist(x1, y1, x2, y2) == 0) {
    return 0;
  }
  int ng = 0, ok = 5000;
  while (ng + 1 < ok) {
    int mid = (ok + ng) / 2;
    if (sq_dist(x1, y1, x2, y2) <= mid * mid) {
      ok = mid;
    }
    else {
      ng = mid;
    }
  }
  return ok;
}

void init_state()
{
  rep(i, node_cnt) { p[i] = 0; }
  rep(i, edge_cnt) { res_y[i] = 0; }
  rep(i, node_cnt) { isNeed[i] = 1; }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void read_input(int problemNum)
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
    cin >> node_cnt >> edge_cnt >> res_cnt;
    rep(i, node_cnt) { cin >> node_x[i] >> node_y[i]; }
    rep(i, edge_cnt) { cin >> u[i] >> v[i] >> w[i]; }
    rep(i, res_cnt) { cin >> res_x[i] >> res_y[i]; }

  }
  // ファイル入力する
  else {
    ifs >> node_cnt >> edge_cnt >> res_cnt;
    rep(i, node_cnt) { ifs >> node_x[i] >> node_y[i]; }
    rep(i, edge_cnt) { ifs >> u[i] >> v[i] >> w[i]; }
    rep(i, res_cnt) { ifs >> res_x[i] >> res_y[i]; }
  }

  rep(i, node_cnt)
  {
    rep(j, node_cnt)
    {
      ndist[i][j] = INF;
      ncost[i][j] = INF;
    }
  }

  rep(i, edge_cnt)
  {
    u[i]--;
    v[i]--;
    ncost[u[i]][v[i]] = w[i];
    ncost[v[i]][u[i]] = w[i];
  }

  rep(i, res_cnt) { nears[i].clear(); }

  rep(i, node_cnt)
  {
    rep(j, res_cnt)
    {
      distn[i][j] = min_radius_le_5000(node_x[i], node_y[i], res_x[j], res_y[j]);
      distk[j][i] = distn[i][j];
      if (distk[j][i] <= 5000) {
        nears[j].push_back(P(distk[j][i], i));
      }
    }
  }

  rep(i, res_cnt) { sort(nears[i].begin(), nears[i].end()); }

  rep(i, edge_cnt)
  {
    edge e;
    e.u = u[i];
    e.v = v[i];
    e.cost = w[i];
    e.id = i;
    edges.push_back(e);
  }

  sort(edges.begin(), edges.end(), comp);  // edge.costが小さい順にソートする

  rep(i, node_cnt)
  {
    rep(j, res_cnt)
    {
      if (distn[i][j] <= 5000) {
        nearsn[i].push_back(P(distn[i][j], j));
      }
    }
  }

  rep(i, node_cnt)
  {
    sort(nearsn[i].begin(), nearsn[i].end());
  }

  init_state();
}

// 解答出力
void write_output(int mode, int problemNum)
{
  if (mode == 0) {
    rep(i, node_cnt) { cout << p[i] << ' '; }
    cout << endl;
    rep(i, edge_cnt) { cout << res_y[i] << ' '; }
    cout << endl;

  }
  else {
    // ファイル出力
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

    rep(i, node_cnt) { ofs << p[i] << ' '; }
    ofs << endl;
    rep(i, edge_cnt) { ofs << res_y[i] << ' '; }
    ofs << endl;

    ofs.close();
  }
}

void snapshot_best()
{
  real_maxScore = maxScore;
  rep(i, node_cnt) { real_p[i] = p[i]; }
  rep(i, edge_cnt) { real_B[i] = B[i]; }
  rep(i, node_cnt) { real_maxPowers[i] = maxPowers[i]; }
  rep(i, node_cnt)
  {
    real_isNeed[i] = isNeed[i];
  }
}

// ランダムに1つ拡大縮小する
void sa_single_power_perturb(double temperature)
{
  int num = rand_u32() % node_cnt;
  int pre = p[num];

  p[num] += rand_u32() % 51 - 25;
  p[num] = max(0LL, p[num]);
  p[num] = min(5000LL, p[num]);

  ll tmpScore = calc_score();

  ll diffScore = tmpScore - maxScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > rand_unit()) {
    maxScore += diffScore;
    if (maxScore > real_maxScore) {
      snapshot_best();
    }
  }
  else {
    // 元に戻す
    p[num] = pre;
  }
}

void solve_all_max()
{
  // 解1
  rep(i, node_cnt) { p[i] = 5000; }
  rep(i, edge_cnt) { res_y[i] = 1; }
}

void solve_greedy_sa()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  outer_kruskal();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }

  int flag[MAX_K] = {};
  vector<P> farest;
  rep(i, res_cnt) { farest.push_back(nears[i][0]); }
  sort(farest.begin(), farest.end());

  drep(i, res_cnt)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, res_cnt)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  // 焼きなまし
  maxScore = calc_score();
  real_maxScore = maxScore;
  rep(i, node_cnt) { real_p[i] = p[i]; }
  rep(i, edge_cnt) { real_B[i] = B[i]; }

  double TL = 1.5;
  int loop = 0;
  double startTemperature = 200048;
  double endTemperature = 0;
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  double nowProgress = nowTime / TL;
  while (true) {
    loop++;
    if (loop % 100 == 0) {
      endTime = clock();
      nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
      if (nowProgress > 1.0) break;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    sa_single_power_perturb(temperature);
  }

  if (run_mode != 0) {
    cout << loop << endl;
  }

  maxScore = real_maxScore;
  rep(i, node_cnt) { p[i] = real_p[i]; }
  rep(i, edge_cnt) { B[i] = real_B[i]; }

  outer_kruskal_lns();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }
}

void solve_layered_sa()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  outer_kruskal();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }

  int flag[MAX_K] = {};

  rep(i, res_cnt)
  {
    maxPowers[nears[i][0].second] =
      max(maxPowers[nears[i][0].second], nears[i][0].first);
  }
  vector<P> farest;
  rep(i, node_cnt) { farest.push_back(P(maxPowers[i], i)); }
  sort(farest.begin(), farest.end());

  drep(i, node_cnt)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, res_cnt)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  // 焼きなまし
  maxScore = calc_score();
  real_maxScore = maxScore;
  rep(i, node_cnt) { real_p[i] = p[i]; }
  rep(i, edge_cnt) { real_B[i] = B[i]; }

  rep(i, node_cnt) { real_maxPowers[i] = maxPowers[i]; }

  double TL_ALL = 1.5;
  int SET_COUNT = 2;
  rep(haibara, SET_COUNT)
  {
    startTime = clock();
    double TL = TL_ALL / SET_COUNT;

    int loop = 0;
    double startTemperature = 1000;
    double endTemperature = 0;
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
    double nowProgress = nowTime / TL;
    while (true) {
      loop++;
      if (loop % 10 == 0) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
        if (nowProgress > 1.0) break;
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      int pre_maxPowers[MAX_N];
      rep(i, node_cnt)
      {
        p[i] = 0;
        pre_maxPowers[i] = maxPowers[i];
      }

      int num = rand_u32() % node_cnt;


      maxPowers[num] += rand_u32() % 101 - 50;
      if (rand_u32() % 10 == 0) {
        maxPowers[num] += rand_u32() % 1001 - 500;
      }
      if (rand_u32() % 100 == 0) {
        maxPowers[num] += rand_u32() % 10001 - 5000;
      }
      maxPowers[num] = max(0, maxPowers[num]);
      maxPowers[num] = min(5000, maxPowers[num]);

      if (rand() % 10 == 0) {
        int num2 = rand_u32() % node_cnt;
        maxPowers[num2] += rand_u32() % 101 - 50;
        if (rand_u32() % 10 == 0) {
          maxPowers[num2] += rand_u32() % 1001 - 500;
        }
        if (rand_u32() % 100 == 0) {
          maxPowers[num2] += rand_u32() % 10001 - 5000;
        }
        maxPowers[num2] = max(0, maxPowers[num2]);
        maxPowers[num2] = min(5000, maxPowers[num2]);
      }


      rep(i, res_cnt) { flag[i] = 0; }

      farest.clear();
      rep(i, node_cnt) { farest.push_back(P(maxPowers[i], i)); }
      sort(farest.begin(), farest.end());

      drep(i, node_cnt)
      {
        int num = farest[i].second;
        int power = farest[i].first;
        if (p[num] != 0) {
          continue;
        }
        ll need = 0;
        for (const auto& pa : nearsn[num]) {
          int j = pa.second;
          ll kyori = pa.first;
          if (kyori > power) {
            break;
          }
          if (!flag[j]) {
            flag[j] = 1;
            need = max(need, kyori);
          }
        }
        p[num] = need;
      }

      ll tmpScore = calc_score(false);

      ll diffScore = tmpScore - maxScore;

      double prob = exp((double)diffScore / temperature);
      if (prob > rand_unit()) {
        maxScore += diffScore;
        if (maxScore > real_maxScore) {
          if (run_mode != 0) {
          }
          snapshot_best();
        }
      }
      else {
        // 元に戻す
        rep(i, node_cnt) { maxPowers[i] = pre_maxPowers[i]; }
      }
    }

    if (run_mode != 0) {
      cout << loop << endl;
    }

    maxScore = real_maxScore;
    rep(i, node_cnt) { p[i] = real_p[i]; }
    rep(i, edge_cnt) { B[i] = real_B[i]; }
    rep(i, node_cnt)
    {
      maxPowers[i] = real_maxPowers[i];
    }
  }

  outer_kruskal_lns();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }
}


ll outer_kruskal_with_score()
{

  ll tmp = outer_kruskal();

  ll res = 0;
  if (tmp != INF) {

    res = calc_score();
    cout << res << endl;
  }


  return res;
}

void solve_sa_edge_toggle()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  outer_kruskal();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }

  int flag[MAX_K] = {};

  rep(i, res_cnt)
  {
    maxPowers[nears[i][0].second] =
      max(maxPowers[nears[i][0].second], nears[i][0].first);
  }
  vector<P> farest;
  rep(i, node_cnt) { farest.push_back(P(maxPowers[i], i)); }
  sort(farest.begin(), farest.end());

  drep(i, node_cnt)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, res_cnt)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  outer_kruskal_lns();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }

  // 焼きなまし
  maxScore = calc_score();
  real_maxScore = maxScore;
  rep(i, node_cnt) { real_p[i] = p[i]; }
  rep(i, edge_cnt) { real_B[i] = B[i]; }
  rep(i, node_cnt) { real_maxPowers[i] = maxPowers[i]; }
  rep(i, node_cnt)
  {
    real_isNeed[i] = isNeed[i];
  }

  double TL_ALL = 1.5;
  int SET_COUNT = 2;
  rep(haibara, SET_COUNT)
  {
    startTime = clock();
    double TL = TL_ALL / SET_COUNT;

    int loop = 0;
    double startTemperature = 1000;
    double endTemperature = 0;
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
    double nowProgress = nowTime / TL;
    while (true) {
      loop++;
      if (loop % 10 == 0) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
        if (nowProgress > 1.0) break;
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      if (rand_u32() % 10 != 0) {
        int pre_maxPowers[MAX_N];
        int pre_p[MAX_N];
        rep(i, node_cnt)
        {
          pre_p[i] = p[i];
          p[i] = 0;
          pre_maxPowers[i] = maxPowers[i];
        }

        int num = rand_u32() % node_cnt;


        maxPowers[num] += rand_u32() % 101 - 50;
        maxPowers[num] = max(0, maxPowers[num]);
        maxPowers[num] = min(5000, maxPowers[num]);

        rep(i, res_cnt) { flag[i] = 0; }

        farest.clear();
        rep(i, node_cnt) { farest.push_back(P(maxPowers[i], i)); }
        sort(farest.begin(), farest.end());

        drep(i, node_cnt)
        {
          int num = farest[i].second;
          int power = farest[i].first;
          if (p[num] != 0) {
            continue;
          }
          ll need = 0;
          for (const auto& pa : nearsn[num]) {
            int j = pa.second;
            ll kyori = pa.first;
            if (kyori > power) {
              break;
            }
            if (!flag[j]) {
              flag[j] = 1;
              need = max(need, kyori);
            }
          }
          p[num] = need;
        }

        int isUpdate = 0;
        vector<int> isNeedKeep;
        rep(i, node_cnt)
        {
          if (pre_p[i] == 0 && p[i] > 0 && isNeed[i] == 0) {
            isNeed[i] = 1;
            isNeedKeep.push_back(i);
            isUpdate = 1;
            break;
          }
        }
        ll tmpScore = 0;
        if (isUpdate) {
          tmpScore = outer_kruskal_with_score();
        }
        else {
          tmpScore = calc_score(false);
        }


        ll diffScore = tmpScore - maxScore;

        double prob = exp((double)diffScore / temperature);
        if (prob > rand_unit()) {
          maxScore += diffScore;
          if (maxScore > real_maxScore) {
            snapshot_best();
          }
        }
        else {
          // 元に戻す
          rep(i, node_cnt) { maxPowers[i] = pre_maxPowers[i]; }
          for (auto ite : isNeedKeep) {
            isNeed[ite] = 0;
          }
        }
      }
      else {
        int num = rand_u32() % (node_cnt - 1) + 1;
        if (p[num] != 0) {
          continue;
        }
        isNeed[num] = 1 - isNeed[num];
        ll tmpScore = outer_kruskal_with_score();

        ll diffScore = tmpScore - maxScore;

        double prob = exp((double)diffScore / temperature);
        if (prob > rand_unit()) {
          maxScore += diffScore;
          if (maxScore > real_maxScore) {
            if (run_mode != 0) {
              cout << maxScore << ' ' << real_maxScore << endl;
            }
            snapshot_best();
          }
        }
        else {
          // 元に戻す
          isNeed[num] = 1 - isNeed[num];
        }
      }
    }

    if (run_mode != 0) {
      cout << loop << endl;
    }

    maxScore = real_maxScore;
    rep(i, node_cnt) { p[i] = real_p[i]; }
    rep(i, edge_cnt) { B[i] = real_B[i]; }
    rep(i, node_cnt)
    {
      maxPowers[i] = real_maxPowers[i];
    }
    rep(i, node_cnt)
    {
      isNeed[i] = real_isNeed[i];
    }
  }

  outer_kruskal_lns();
  rep(i, edge_cnt) { res_y[i] = isUse[i]; }
}


double solve_case(int mode, int problemNum = 0)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // solve_all_max();
  // solve_greedy_sa();
  solve_layered_sa();
  //solve_sa_edge_toggle();

  check_coverage();

  if (mode != 0) {
    cout << calc_score() << endl;
  }
  return calc_score();
}

double run_with_io(int mode, int problemNum = 0)
{
  // 入力受け取り
  read_input(problemNum);

  double score = solve_case(mode, problemNum);

  // 解答の出力
  write_output(mode, problemNum);

  return score;
}

int main()
{
  run_mode = 2;

  // 提出用
  if (run_mode == 0) {
    run_with_io(run_mode);
  }
  // 1ケース試す
  else if (run_mode == 1) {
    run_with_io(run_mode, 9);
  }
  // 複数ケース試す
  else if (run_mode == 2) {
    rep(i, 10) { run_with_io(run_mode, i); }
  }

  return 0;
}
