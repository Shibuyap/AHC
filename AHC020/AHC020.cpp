#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

int mode;

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

const ll INF = 1001001001001001001;

struct Point
{
  ll x;
  ll y;
};

struct Edge
{
  int u;
  int v;
  ll cost;
  int id;
};

const int MAX_N = 100;
const int MAX_M = 300;
const int MAX_K = 5000;

// 定数
vector<Edge> edges;
int n, m, k;
Point nodes[MAX_N];

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
int best_mst_node_req[MAX_N];
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
    if (x == y) { return; }

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
int mst_edge_use[MAX_K];
/*----------  Kruskal core  ----------*/
long long kruskal(int vertex_cnt)
{
  int need_cnt = 0;
  for (int v_idx = 0; v_idx < vertex_cnt; ++v_idx) {
    need_cnt += mst_node_req[v_idx];
  }

  int edge_cnt = edges.size();
  for (int e_idx = 0; e_idx < edge_cnt; ++e_idx) {
    mst_edge_use[e_idx] = 0;
  }

  /* Union-Find work arrays */
  vector<int> uf_parent(vertex_cnt + 10);
  vector<int> uf_rank(vertex_cnt + 10);
  init(uf_parent, uf_rank, vertex_cnt + 10);

  long long mst_cost = 0;
  int unite_cnt = 0;

  for (int e_idx = 0; e_idx < edge_cnt; ++e_idx) {
    const Edge& ed = edges[e_idx];
    if (!mst_node_req[ed.u] || !mst_node_req[ed.v]) {
      continue;
    }

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
  return kruskal(n);
}

/*----------  LNS-style toggle search  ----------*/
ll outer_kruskal_lns()
{
  start_timer();

  for (int v_idx = 0; v_idx < n; ++v_idx) {
    mst_node_req[v_idx] = 1;
  }

  ll best_cost = outer_kruskal();
  int iter_cnt = 0;
  const double time_limit = 0.1;

  while (true) {
    ++iter_cnt;

    /* タイムリミット監視 */
    double sec_elapsed = get_elapsed_time();
    double ratio_elapsed = sec_elapsed / time_limit;
    if (ratio_elapsed > 1.0) {
      break;
    }

    /* ランダム頂点をトグル */
    int node_id = rand32() % (n - 1) + 1;
    if (power_rad[node_id] != 0) {
      continue;
    }

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
  for (int res_idx = 0; res_idx < k; ++res_idx) {
    bool is_covered = false;
    int near_sz = res_near_nodes[res_idx].size();
    for (int near_idx = 0; near_idx < near_sz; ++near_idx) {
      int node_id = res_near_nodes[res_idx][near_idx].second;
      if (res_node_dist[res_idx][node_id] <= power_rad[node_id]) {
        is_covered = true;
        break;
      }
    }
    if (!is_covered) {
      return false;
    }
  }
  return true;
}

bool check_coverage_node(int target_node)
{
  for (auto near_pair : node_near_res[target_node]) {
    int res_idx = near_pair.second;
    bool is_covered = false;
    int near_sz = res_near_nodes[res_idx].size();
    for (int near_idx = 0; near_idx < near_sz; ++near_idx) {
      int node_id = res_near_nodes[res_idx][near_idx].second;
      if (res_node_dist[res_idx][node_id] <= power_rad[node_id]) {
        is_covered = true;
        break;
      }
    }
    if (!is_covered) {
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
  for (int i = 0; i < n; ++i) {
    S += power_rad[i] * power_rad[i];
  }
  if (isALL) {
    keepW = 0;
    for (int i = 0; i < m; ++i) {
      if (ress[i].y) {
        S += edges[i].cost;
        keepW += edges[i].cost;
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
  for (int i = 0; i < n; ++i) {
    S += power_rad[i] * power_rad[i];
  }
  if (isALL) {
    keepW = 0;
    for (int i = 0; i < m; ++i) {
      if (ress[i].y) {
        S += edges[i].cost;
        keepW += edges[i].cost;
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
  for (int i = 0; i < n; ++i) {
    power_rad[i] = 0;
  }
  for (int i = 0; i < m; ++i) {
    ress[i].y = 0;
  }
  for (int i = 0; i < n; ++i) {
    mst_node_req[i] = 1;
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void read_input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  edges.clear();

  if (!ifs.is_open()) {
    // 標準入力する
    cin >> n >> m >> k;
    for (int i = 0; i < n; ++i) {
      cin >> nodes[i].x >> nodes[i].y;
    }
    for (int i = 0; i < m; ++i) {
      Edge e;
      cin >> e.u >> e.v >> e.cost;
      e.u--;
      e.v--;
      e.id = i;
      edges.push_back(e);
    }
    for (int i = 0; i < k; ++i) {
      cin >> ress[i].x >> ress[i].y;
    }

  }
  else {
    // ファイル入力する
    ifs >> n >> m >> k;
    for (int i = 0; i < n; ++i) {
      ifs >> nodes[i].x >> nodes[i].y;
    }
    for (int i = 0; i < m; ++i) {
      Edge e;
      ifs >> e.u >> e.v >> e.cost;
      e.u--;
      e.v--;
      e.id = i;
      edges.push_back(e);
    }
    for (int i = 0; i < k; ++i) {
      ifs >> ress[i].x >> ress[i].y;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      node_node_dist[i][j] = INF;
      node_node_cost[i][j] = INF;
    }
  }

  for (int i = 0; i < m; ++i) {
    node_node_cost[edges[i].u][edges[i].v] = edges[i].cost;
    node_node_cost[edges[i].v][edges[i].u] = edges[i].cost;
  }

  for (int i = 0; i < k; ++i) {
    res_near_nodes[i].clear();
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      node_res_dist[i][j] = min_radius_le_5000(nodes[i].x, nodes[i].y, ress[j].x, ress[j].y);
      res_node_dist[j][i] = node_res_dist[i][j];
      if (res_node_dist[j][i] <= 5000) {
        res_near_nodes[j].push_back(P(res_node_dist[j][i], i));
      }
    }
  }

  for (int i = 0; i < k; ++i) {
    sort(res_near_nodes[i].begin(), res_near_nodes[i].end());
  }

  sort(edges.begin(), edges.end(), comp);  // Edge.costが小さい順にソートする

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      if (node_res_dist[i][j] <= 5000) {
        node_near_res[i].push_back(P(node_res_dist[i][j], j));
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    sort(node_near_res[i].begin(), node_near_res[i].end());
  }

  init_state();
}

// 解答出力
void write_output(int mode, int problemNum)
{
  if (mode == 0) {
    for (int i = 0; i < n; ++i) {
      cout << power_rad[i] << ' ';
    }
    cout << endl;
    for (int i = 0; i < m; ++i) {
      cout << ress[i].y << ' ';
    }
    cout << endl;
  }
  else {
    // ファイル出力
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    for (int i = 0; i < n; ++i) {
      ofs << power_rad[i] << ' ';
    }
    ofs << endl;
    for (int i = 0; i < m; ++i) {
      ofs << ress[i].y << ' ';
    }
    ofs << endl;

    ofs.close();
  }
}

void snapshot_best()
{
  best_score = best_score_cur;
  for (int i = 0; i < n; ++i) {
    best_power_rad[i] = power_rad[i];
  }
  for (int i = 0; i < m; ++i) {
    best_edge_sel[i] = edge_sel[i];
  }
  for (int i = 0; i < n; ++i) {
    best_power_cap[i] = power_cap[i];
  }
  for (int i = 0; i < n; ++i) {
    best_mst_node_req[i] = mst_node_req[i];
  }
}

// ランダムに1つ拡大縮小する
void sa_single_power_perturb(double temperature)
{
  int num = rand32() % n;
  int pre = power_rad[num];

  power_rad[num] += rand32() % 51 - 25;
  power_rad[num] = max(0LL, power_rad[num]);
  power_rad[num] = min(5000LL, power_rad[num]);

  ll tmpScore = calc_score();

  ll diffScore = tmpScore - best_score_cur;

  double prob = exp((double)diffScore / temperature);
  if (prob > rand_01()) {
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
  for (int edge_idx = 0; edge_idx < m; ++edge_idx) {
    ress[edge_idx].y = mst_edge_use[edge_idx];
  }

  int res_cover_flag[MAX_K] = {};

  /* 各住民の最も近い頂点で power_cap を下限設定 */
  for (int res_idx = 0; res_idx < k; ++res_idx) {
    int node_id = res_near_nodes[res_idx][0].second;
    power_cap[node_id] = max(power_cap[node_id], res_near_nodes[res_idx][0].first);
  }

  vector<P> farthest;
  for (int node_idx = 0; node_idx < n; ++node_idx) {
    farthest.emplace_back(power_cap[node_idx], node_idx);
  }
  sort(farthest.begin(), farthest.end());

  /* power_cap が大きい頂点から順に半径を割り当て */
  for (int f_idx = n - 1; f_idx >= 0; --f_idx) {
    int node_id = farthest[f_idx].second;
    int cap_value = farthest[f_idx].first;
    if (power_rad[node_id] != 0) { continue; }

    ll need_radius = 0;
    for (int res_idx = 0; res_idx < k; ++res_idx) {
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
  for (int node_idx = 0; node_idx < n; ++node_idx) {
    best_power_rad[node_idx] = power_rad[node_idx];
  }
  for (int edge_idx = 0; edge_idx < m; ++edge_idx) {
    best_edge_sel[edge_idx] = edge_sel[edge_idx];
  }
  for (int node_idx = 0; node_idx < n; ++node_idx) {
    best_power_cap[node_idx] = power_cap[node_idx];
  }

  /*--- 多段 SA ---*/
  const double tl_total = 1.5;
  const int    stage_cnt = 2;

  for (int stage_idx = 0; stage_idx < stage_cnt; ++stage_idx) {
    start_timer();
    double tl_stage = tl_total / stage_cnt;

    int    iter_cnt = 0;
    double temp_ini = 1000.0;
    double temp_fin = 0.0;

    /* -------- SA 本体 -------- */
    int prev_power_cap[MAX_N];
    while (true) {
      ++iter_cnt;
      if (iter_cnt % 10 == 0 && get_elapsed_time() / tl_stage > 1.0) {
        break;
      }

      double progress = get_elapsed_time() / tl_stage;
      double temperature = temp_ini + (temp_fin - temp_ini) * progress;

      for (int node_idx = 0; node_idx < n; ++node_idx) {
        power_rad[node_idx] = 0;
        prev_power_cap[node_idx] = power_cap[node_idx];
      }

      /*----- ランダムに power_cap をいじる -----*/
      int pick_node = rand32() % n;
      auto mutate_cap = [&](int id) {
        power_cap[id] += rand32() % 101 - 50;
        if (rand32() % 10 == 0) {
          power_cap[id] += rand32() % 1001 - 500;
        }
        if (rand32() % 100 == 0) {
          power_cap[id] += rand32() % 10001 - 5000;
        }
        power_cap[id] = clamp(power_cap[id], 0, 5000);
        };
      mutate_cap(pick_node);
      if (rand32() % 10 == 0) mutate_cap(rand32() % n);

      /*----- 新しい半径割り当て -----*/
      for (int res_idx = 0; res_idx < k; ++res_idx) {
        res_cover_flag[res_idx] = 0;
      }

      farthest.clear();
      for (int node_idx = 0; node_idx < n; ++node_idx) {
        farthest.emplace_back(power_cap[node_idx], node_idx);
      }
      sort(farthest.begin(), farthest.end());

      for (int f_idx = n - 1; f_idx >= 0; --f_idx) {
        int node_id = farthest[f_idx].second;
        int cap_value = farthest[f_idx].first;
        if (power_rad[node_id] != 0) {
          continue;
        }

        ll need_radius = 0;
        for (const auto& pr : node_near_res[node_id]) {
          int res_idx = pr.second;
          ll dist_val = pr.first;
          if (dist_val > cap_value) {
            break;
          }
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

      if (prob > rand_01()) {
        best_score_cur += delta_score;
        if (best_score_cur > best_score) {
          snapshot_best();
        }
      }
      else {
        for (int node_idx = 0; node_idx < n; ++node_idx) {
          power_cap[node_idx] = prev_power_cap[node_idx];
        }
      }
    }

    if (mode != 0) cout << iter_cnt << '\n';

    /* ステージ終了：状態をリセット */
    best_score_cur = best_score;
    for (int node_idx = 0; node_idx < n; ++node_idx) {
      power_rad[node_idx] = best_power_rad[node_idx];
      power_cap[node_idx] = best_power_cap[node_idx];
    }
    for (int edge_idx = 0; edge_idx < m; ++edge_idx) {
      edge_sel[edge_idx] = best_edge_sel[edge_idx];
    }
  }

  /*--- LNS で接続を再最適化 ---*/
  outer_kruskal_lns();
  for (int edge_idx = 0; edge_idx < m; ++edge_idx) {
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
  read_input(problemNum);

  double score = solve_case(mode, problemNum);

  write_output(mode, problemNum);

  return score;
}

int main()
{
  mode = 2;

  // 提出用
  if (mode == 0) {
    run_with_io(mode);
  }
  // 1ケース試す
  else if (mode == 1) {
    run_with_io(mode, 9);
  }
  // 複数ケース試す
  else if (mode == 2) {
    for (int i = 0; i < 10; ++i) { run_with_io(mode, i); }
  }

  return 0;
}
