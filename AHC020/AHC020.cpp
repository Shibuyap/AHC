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

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

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

  static double rand_unit()
  {
    return (rand_u32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

const ll INF = 1001001001001001001;

struct Point
{
  ll x;
  ll y;
};

struct Edge
{
  int u, v, cost;
  int id;
};

const int MAX_N = 100;
const int MAX_M = 300;
const int MAX_K = 5000;

// 定数
vector<Edge> edges;
int node_cnt, edge_cnt, res_cnt;
Point nodes[MAX_N];
ll edge_u[MAX_M], edge_v[MAX_M], edge_cost_raw[MAX_M];
Point ress[MAX_K];
ll node_res_dist[MAX_N][MAX_K];
ll res_node_dist[MAX_K][MAX_N];
ll node_node_dist[MAX_N][MAX_N];
ll node_node_cost[MAX_N][MAX_N];
vector<P> res_near_nodes[MAX_K];
vector<P> node_near_res[MAX_N];

// 更新する変数
ll best_score_cur;
ll best_score;
ll power_rad[MAX_N];
ll edge_sel[MAX_M];
ll best_power_rad[MAX_N];
ll best_edge_sel[MAX_M];
int power_cap[MAX_N] = {};
int best_power_cap[MAX_N] = {};

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

bool comp(const Edge& e1, const Edge& e2)
{
  return e1.cost < e2.cost;
}

int mst_node_req[MAX_N];
int best_mst_node_req[MAX_N];
int mst_edge_use[MAX_K];
/*----------  Kruskal core  ----------*/
long long kruskal(int vertex_cnt)
{
  int need_cnt = 0;
  rep(v_idx, vertex_cnt)
  {
    need_cnt += mst_node_req[v_idx];
  }

  int edge_cnt = edges.size();
  rep(e_idx, edge_cnt) mst_edge_use[e_idx] = 0;

  /* Union-Find work arrays */
  vector<int> uf_parent(vertex_cnt + 10);
  vector<int> uf_rank(vertex_cnt + 10);
  init(uf_parent, uf_rank, vertex_cnt + 10);

  long long mst_cost = 0;
  int unite_cnt = 0;

  for (int e_idx = 0; e_idx < edge_cnt; ++e_idx) {
    const Edge& ed = edges[e_idx];
    if (!mst_node_req[ed.u] || !mst_node_req[ed.v]) continue;

    if (!same(uf_parent, ed.u, ed.v)) {
      mst_edge_use[ed.id] = 1;
      unite(uf_parent, uf_rank, ed.u, ed.v);
      mst_cost += ed.cost;
      ++unite_cnt;
    }
  }

  if (unite_cnt != need_cnt - 1) mst_cost = INF;
  return mst_cost;
}

ll outer_kruskal()
{
  return kruskal(node_cnt);
}

/*----------  LNS-style toggle search  ----------*/
ll outer_kruskal_lns()
{
  start_timer();

  rep(v_idx, node_cnt) mst_node_req[v_idx] = 1;

  ll best_cost = outer_kruskal();
  int iter_cnt = 0;
  const double time_limit = 0.1;

  while (true) {
    ++iter_cnt;

    /* タイムリミット監視 */
    double sec_elapsed = get_elapsed_time();
    double ratio_elapsed = sec_elapsed / time_limit;
    if (ratio_elapsed > 1.0) break;

    /* ランダム頂点をトグル */
    int node_id = rand_u32() % (node_cnt - 1) + 1;
    if (power_rad[node_id] != 0) continue;

    mst_node_req[node_id] ^= 1;
    ll cur_cost = outer_kruskal();

    if (cur_cost <= best_cost) {
      best_cost = cur_cost;
    }
    else {
      mst_node_req[node_id] ^= 1;
    }
  }
  return best_cost;
}

bool check_coverage()
{
  rep(res_idx, res_cnt)
  {
    bool is_covered = false;
    int near_sz = res_near_nodes[res_idx].size();
    rep(near_idx, near_sz)
    {
      int node_id = res_near_nodes[res_idx][near_idx].second;
      if (res_node_dist[res_idx][node_id] <= power_rad[node_id]) {
        is_covered = true;
        break;
      }
    }
    if (!is_covered) return false;
  }
  return true;
}

bool check_coverage_node(int target_node)
{
  for (auto near_pair : node_near_res[target_node]) {
    int res_idx = near_pair.second;
    bool is_covered = false;
    int near_sz = res_near_nodes[res_idx].size();
    rep(near_idx, near_sz)
    {
      int node_id = res_near_nodes[res_idx][near_idx].second;
      if (res_node_dist[res_idx][node_id] <= power_rad[node_id]) {
        is_covered = true;
        break;
      }
    }
    if (!is_covered) return false;
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
  rep(i, node_cnt)
  {
    S += power_rad[i] * power_rad[i];
  }
  if (isALL) {
    keepW = 0;
    rep(i, edge_cnt)
    {
      if (ress[i].y) {
        S += edge_cost_raw[i];
        keepW += edge_cost_raw[i];
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
  rep(i, node_cnt)
  {
    S += power_rad[i] * power_rad[i];
  }
  if (isALL) {
    keepW = 0;
    rep(i, edge_cnt)
    {
      if (ress[i].y) {
        S += edge_cost_raw[i];
        keepW += edge_cost_raw[i];
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
  rep(i, node_cnt)
  {
    power_rad[i] = 0;
  }
  rep(i, edge_cnt)
  {
    ress[i].y = 0;
  }
  rep(i, node_cnt)
  {
    mst_node_req[i] = 1;
  }
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
    rep(i, node_cnt)
    {
      cin >> nodes[i].x >> nodes[i].y;
    }
    rep(i, edge_cnt)
    {
      cin >> edge_u[i] >> edge_v[i] >> edge_cost_raw[i];
    }
    rep(i, res_cnt)
    {
      cin >> ress[i].x >> ress[i].y;
    }

  }
  // ファイル入力する
  else {
    ifs >> node_cnt >> edge_cnt >> res_cnt;
    rep(i, node_cnt)
    {
      ifs >> nodes[i].x >> nodes[i].y;
    }
    rep(i, edge_cnt)
    {
      ifs >> edge_u[i] >> edge_v[i] >> edge_cost_raw[i];
    }
    rep(i, res_cnt)
    {
      ifs >> ress[i].x >> ress[i].y;
    }
  }

  rep(i, node_cnt)
  {
    rep(j, node_cnt)
    {
      node_node_dist[i][j] = INF;
      node_node_cost[i][j] = INF;
    }
  }

  rep(i, edge_cnt)
  {
    edge_u[i]--;
    edge_v[i]--;
    node_node_cost[edge_u[i]][edge_v[i]] = edge_cost_raw[i];
    node_node_cost[edge_v[i]][edge_u[i]] = edge_cost_raw[i];
  }

  rep(i, res_cnt)
  {
    res_near_nodes[i].clear();
  }

  rep(i, node_cnt)
  {
    rep(j, res_cnt)
    {
      node_res_dist[i][j] = min_radius_le_5000(nodes[i].x, nodes[i].y, ress[j].x, ress[j].y);
      res_node_dist[j][i] = node_res_dist[i][j];
      if (res_node_dist[j][i] <= 5000) {
        res_near_nodes[j].push_back(P(res_node_dist[j][i], i));
      }
    }
  }

  rep(i, res_cnt)
  {
    sort(res_near_nodes[i].begin(), res_near_nodes[i].end());
  }

  rep(i, edge_cnt)
  {
    Edge e;
    e.u = edge_u[i];
    e.v = edge_v[i];
    e.cost = edge_cost_raw[i];
    e.id = i;
    edges.push_back(e);
  }

  sort(edges.begin(), edges.end(), comp);  // Edge.costが小さい順にソートする

  rep(i, node_cnt)
  {
    rep(j, res_cnt)
    {
      if (node_res_dist[i][j] <= 5000) {
        node_near_res[i].push_back(P(node_res_dist[i][j], j));
      }
    }
  }

  rep(i, node_cnt)
  {
    sort(node_near_res[i].begin(), node_near_res[i].end());
  }

  init_state();
}

// 解答出力
void write_output(int mode, int problemNum)
{
  if (mode == 0) {
    rep(i, node_cnt)
    {
      cout << power_rad[i] << ' ';
    }
    cout << endl;
    rep(i, edge_cnt)
    {
      cout << ress[i].y << ' ';
    }
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

    rep(i, node_cnt)
    {
      ofs << power_rad[i] << ' ';
    }
    ofs << endl;
    rep(i, edge_cnt)
    {
      ofs << ress[i].y << ' ';
    }
    ofs << endl;

    ofs.close();
  }
}

void snapshot_best()
{
  best_score = best_score_cur;
  rep(i, node_cnt)
  {
    best_power_rad[i] = power_rad[i];
  }
  rep(i, edge_cnt)
  {
    best_edge_sel[i] = edge_sel[i];
  }
  rep(i, node_cnt)
  {
    best_power_cap[i] = power_cap[i];
  }
  rep(i, node_cnt)
  {
    best_mst_node_req[i] = mst_node_req[i];
  }
}

// ランダムに1つ拡大縮小する
void sa_single_power_perturb(double temperature)
{
  int num = rand_u32() % node_cnt;
  int pre = power_rad[num];

  power_rad[num] += rand_u32() % 51 - 25;
  power_rad[num] = max(0LL, power_rad[num]);
  power_rad[num] = min(5000LL, power_rad[num]);

  ll tmpScore = calc_score();

  ll diffScore = tmpScore - best_score_cur;

  double prob = exp((double)diffScore / temperature);
  if (prob > rand_unit()) {
    best_score_cur += diffScore;
    if (best_score_cur > best_score) {
      snapshot_best();
    }
  }
  else {
    // 元に戻す
    power_rad[num] = pre;
  }
}

void solve_layered_sa()
{
  /*--- 初期 MST & 半径決定 ---*/
  outer_kruskal();
  rep(edge_idx, edge_cnt) { ress[edge_idx].y = mst_edge_use[edge_idx]; }

  int res_cover_flag[MAX_K] = {};

  /* 各住民の最も近い頂点で power_cap を下限設定 */
  rep(res_idx, res_cnt)
  {
    int node_id = res_near_nodes[res_idx][0].second;
    power_cap[node_id] =
      max(power_cap[node_id], res_near_nodes[res_idx][0].first);
  }

  vector<P> farthest;
  rep(node_idx, node_cnt) { farthest.emplace_back(power_cap[node_idx], node_idx); }
  sort(farthest.begin(), farthest.end());

  /* power_cap が大きい頂点から順に半径を割り当て */
  drep(f_idx, node_cnt)
  {
    int node_id = farthest[f_idx].second;
    int cap_value = farthest[f_idx].first;
    if (power_rad[node_id] != 0) continue;

    ll need_radius = 0;
    rep(res_idx, res_cnt)
    {
      if (!res_cover_flag[res_idx] && res_node_dist[res_idx][node_id] <= cap_value) {
        res_cover_flag[res_idx] = 1;
        need_radius = max(need_radius, res_node_dist[res_idx][node_id]);
      }
    }
    power_rad[node_id] = need_radius;
  }

  /*--- SA 用のベストスナップショット保存 ---*/
  best_score_cur = calc_score();
  best_score = best_score_cur;
  rep(node_idx, node_cnt) { best_power_rad[node_idx] = power_rad[node_idx]; }
  rep(edge_idx, edge_cnt) { best_edge_sel[edge_idx] = edge_sel[edge_idx]; }
  rep(node_idx, node_cnt) { best_power_cap[node_idx] = power_cap[node_idx]; }

  /*--- 多段 SA ---*/
  const double tl_total = 1.5;
  const int    stage_cnt = 2;

  rep(stage_idx, stage_cnt)
  {
    start_timer();
    double tl_stage = tl_total / stage_cnt;

    int    iter_cnt = 0;
    double temp_ini = 1000.0;
    double temp_fin = 0.0;

    /* -------- SA 本体 -------- */
    while (true) {
      ++iter_cnt;
      if (iter_cnt % 10 == 0 && get_elapsed_time() / tl_stage > 1.0) break;

      double progress = get_elapsed_time() / tl_stage;
      double temperature = temp_ini + (temp_fin - temp_ini) * progress;

      int prev_power_cap[MAX_N];
      rep(node_idx, node_cnt)
      {
        power_rad[node_idx] = 0;
        prev_power_cap[node_idx] = power_cap[node_idx];
      }

      /*----- ランダムに power_cap をいじる -----*/
      int pick_node = rand_u32() % node_cnt;
      auto mutate_cap = [&](int id) {
        power_cap[id] += rand_u32() % 101 - 50;
        if (rand_u32() % 10 == 0) power_cap[id] += rand_u32() % 1001 - 500;
        if (rand_u32() % 100 == 0) power_cap[id] += rand_u32() % 10001 - 5000;
        power_cap[id] = clamp(power_cap[id], 0, 5000);
        };
      mutate_cap(pick_node);
      if (rand_u32() % 10 == 0) mutate_cap(rand_u32() % node_cnt);

      /*----- 新しい半径割り当て -----*/
      rep(res_idx, res_cnt) res_cover_flag[res_idx] = 0;

      farthest.clear();
      rep(node_idx, node_cnt) farthest.emplace_back(power_cap[node_idx], node_idx);
      sort(farthest.begin(), farthest.end());

      drep(f_idx, node_cnt)
      {
        int node_id = farthest[f_idx].second;
        int cap_value = farthest[f_idx].first;
        if (power_rad[node_id] != 0) continue;

        ll need_radius = 0;
        for (const auto& pr : node_near_res[node_id]) {
          int res_idx = pr.second;
          ll dist_val = pr.first;
          if (dist_val > cap_value) break;
          if (!res_cover_flag[res_idx]) {
            res_cover_flag[res_idx] = 1;
            need_radius = max(need_radius, dist_val);
          }
        }
        power_rad[node_id] = need_radius;
      }

      /*----- 受理判定 -----*/
      ll new_score = calc_score(false);
      ll delta_score = new_score - best_score_cur;
      double prob = exp(static_cast<double>(delta_score) / temperature);

      if (prob > rand_unit()) {
        best_score_cur += delta_score;
        if (best_score_cur > best_score) snapshot_best();
      }
      else {
        rep(node_idx, node_cnt) power_cap[node_idx] = prev_power_cap[node_idx];
      }
    }

    if (run_mode != 0) cout << iter_cnt << '\n';

    /* ステージ終了：状態をリセット */
    best_score_cur = best_score;
    rep(node_idx, node_cnt)
    {
      power_rad[node_idx] = best_power_rad[node_idx];
      power_cap[node_idx] = best_power_cap[node_idx];
    }
    rep(edge_idx, edge_cnt) edge_sel[edge_idx] = best_edge_sel[edge_idx];
  }

  /*--- LNS で接続を再最適化 ---*/
  outer_kruskal_lns();
  rep(edge_idx, edge_cnt)
  {
    ress[edge_idx].y = mst_edge_use[edge_idx];
  }
}

double solve_case(int mode, int problemNum = 0)
{
  start_timer();

  solve_layered_sa();

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
