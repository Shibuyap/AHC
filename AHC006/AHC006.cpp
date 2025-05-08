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
#define MAX_N 200005

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

int exec_mode = 0;

struct Point
{
  int x;
  int y;

  Point()
  {

  }

  Point(int _x, int _y)
  {
    x = _x;
    y = _y;
  }
};

namespace /* 変数 */
{
  const int INIT_N = 1000;

  int n = 1000;
  const int m = 55;
  vector<Point> start_points(INIT_N), goal_points(INIT_N);
  vector<int> orig_indices;

  ll curr_score;
  ll curr_length;
  vector<Point> current_path;

  vector<int> sel_pair_idx(m * 2);
  vector<int> reduced_pair_idx(m * 2);
  vector<int> is_used(INIT_N);
  vector<int> partner_idx(INIT_N * 2);

  ll best_score;
  ll best_length;
  vector<Point> best_path;

  vector<int> best_sel_pair_idx;
  vector<int> best_reduced_pair_idx;
  vector<int> best_is_used;
  vector<int> best_partner_idx;

  void ResetParam()
  {
    n = INIT_N;

    start_points.resize(INIT_N);
    goal_points.resize(INIT_N);
    orig_indices.clear();

    curr_score = 0;
    curr_length = 0;
    current_path.clear();

    sel_pair_idx.clear();
    reduced_pair_idx.clear();
    is_used.clear();
    partner_idx.clear();

    sel_pair_idx.resize(m * 2);
    reduced_pair_idx.resize(m * 2);
    is_used.resize(INIT_N);
    partner_idx.resize(INIT_N * 2);

    best_score = 0;
    best_length = 0;
    best_path.clear();

    //best_sel_pair_idx.clear();
    //best_reduced_pair_idx.clear();
    //best_is_used.clear();
    //best_partner_idx.clear();

    //best_sel_pair_idx.resize(m * 2);
    //best_reduced_pair_idx.resize(m * 2);
    //best_is_used.resize(INIT_N);
    //best_partner_idx.resize(INIT_N * 2);
  }
}  // namespace

// スコア計算
ll CalcScore(ll time)
{
  ll res = round(100000000.0 / (1000.0 + time));
  return res;
}

inline int manhattan(const Point& a, const Point& b)
{
  return abs(a.x - b.x) + abs(a.y - b.y);
}

ll compute_path_time()
{
  ll timeSum = 0;
  srep(i, 1, current_path.size())
  {
    timeSum += manhattan(current_path[i], current_path[i - 1]);
  }
  return timeSum;
}

inline bool is_inside_rect(int ite, int L, int R, int U, int D)
{
  if (start_points[ite].x < L || R < start_points[ite].x) return false;
  if (start_points[ite].y < U || D < start_points[ite].y) return false;
  if (goal_points[ite].x < L || R < goal_points[ite].x) return false;
  if (goal_points[ite].y < U || D < goal_points[ite].y) return false;
  return true;
}

void filter_input_points()
{
  // --- 収集済みの最良データを保持 ---
  vector<Point> best_start_points, best_goal_points;
  int         best_metric = 1'001'001;   // 今のところ最小の「矩形スコア」
  int         iter_count = 0;

  // --- 改善ループ ---
  while (true) {
    ++iter_count;
    if (iter_count % 100 == 1 && get_elapsed_time() > 0.2) break;

    // 候補データを一時的に保持
    vector<Point> cand_start_points, cand_goal_points;
    vector<int> cand_indices;

    // ランダムに矩形を生成
    int rect_left = Rand() % 801;
    int rect_right = Rand() % 801;
    int rect_upper = Rand() % 801;
    int rect_lower = Rand() % 801;
    if (rect_left > rect_right)  swap(rect_left, rect_right);
    if (rect_upper > rect_lower)  swap(rect_upper, rect_lower);

    // “広さ” の代わりに幅＋高さを評価指標にしている
    const int rect_metric = (rect_right - rect_left) + (rect_lower - rect_upper);
    if (rect_metric >= best_metric) continue;

    // 矩形内に収まる点を抽出
    for (int i = 0; i < n; ++i) {
      if (is_inside_rect(i, rect_left, rect_right, rect_upper, rect_lower)) {
        cand_start_points.push_back(start_points[i]);
        cand_goal_points.push_back(goal_points[i]);
        cand_indices.push_back(i);
      }
    }

    // 十分な点数があり、しかも指標が更新されたら採用
    if (cand_start_points.size() >= m + 1) {
      best_start_points = std::move(cand_start_points);
      best_goal_points = std::move(cand_goal_points);
      orig_indices = std::move(cand_indices);
      best_metric = rect_metric;
    }
  }

  // ベスト結果でグローバルを更新
  if (best_start_points.size() >= m) {
    start_points.swap(best_start_points);
    goal_points.swap(best_goal_points);
    n = static_cast<int>(start_points.size());
  }
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    for (int i = 0; i < n; ++i) {
      cin >> start_points[i].x >> start_points[i].y >> goal_points[i].x >> goal_points[i].y;
    }
  }
  else {
    // ファイル入力
    for (int i = 0; i < n; ++i) {
      ifs >> start_points[i].x >> start_points[i].y >> goal_points[i].x >> goal_points[i].y;
    }
  }
}

void output_data(int case_num)
{
  if (exec_mode == 0) {
    // 標準出力
    cout << 50;
    for (int i = 0; i < 50; ++i) {
      if (sel_pair_idx[i] < INIT_N) {
        cout << " " << orig_indices[sel_pair_idx[i]] + 1;
      }
    }
    cout << endl;
    cout << current_path.size();
    for (int i = 0; i < current_path.size(); ++i) {
      cout << " " << current_path[i].x << " " << current_path[i].y;
    }
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << 50;
    for (int i = 0; i < 50; ++i) {
      if (sel_pair_idx[i] < INIT_N) {
        ofs << " " << orig_indices[sel_pair_idx[i]] + 1;
      }
    }
    ofs << endl;
    ofs << current_path.size();
    for (int i = 0; i < current_path.size(); ++i) {
      ofs << " " << current_path[i].x << " " << current_path[i].y;
    }
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

void build_initial_path()
{
  // 愚直解
  current_path.push_back(Point(400, 400));
  for (int i = 0; i < m; ++i) {
    current_path.push_back(Point(start_points[i].x, start_points[i].y));
    sel_pair_idx[i] = i;
    is_used[i] = 1;
  }
  for (int i = 0; i < m; ++i) {
    current_path.push_back(Point(goal_points[i].x, goal_points[i].y));
  }
  current_path.push_back(Point(400, 400));
  curr_length = compute_path_time();
  curr_score = CalcScore(curr_length);

  best_path = current_path;
  best_length = curr_length;
  best_score = curr_score;
}

void sa_point_swap()
{
  // 山登り解、焼きなまし解
  double time_limit = 0.4;
  double start_temp = 4800;
  double end_temp = 1000;
  int loop = 0;
  double now_time = get_elapsed_time();
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) {
        break;
      }
    }

    int cand_idx = Rand() % n;
    while (is_used[cand_idx]) {
      cand_idx = Rand() % n;
    }

    int path_pos = Rand() % m;

    int keepNum = sel_pair_idx[path_pos];
    int diff = 0;
    diff += manhattan(start_points[cand_idx], current_path[path_pos]);
    diff += manhattan(start_points[cand_idx], current_path[path_pos + 2]);
    diff += manhattan(goal_points[cand_idx], current_path[path_pos + m]);
    diff += manhattan(goal_points[cand_idx], current_path[path_pos + m + 2]);
    diff -= manhattan(current_path[path_pos + 1], current_path[path_pos]);
    diff -= manhattan(current_path[path_pos + 1], current_path[path_pos + 2]);
    diff -= manhattan(current_path[path_pos + m + 1], current_path[path_pos + m]);
    diff -= manhattan(current_path[path_pos + m + 1], current_path[path_pos + m + 2]);

    int new_length = curr_length + diff;
    int diffScore = -diff;

    double temp = start_temp + (end_temp - start_temp) * now_time / time_limit;
    double prob = exp((double)diffScore / temp);
    if (prob > Rand01()) {
      curr_length = new_length;
      is_used[sel_pair_idx[path_pos]] = 0;
      is_used[cand_idx] = 1;
      sel_pair_idx[path_pos] = cand_idx;
      current_path[path_pos + 1] = Point(start_points[cand_idx].x, start_points[cand_idx].y);
      current_path[path_pos + m + 1] = Point(goal_points[cand_idx].x, goal_points[cand_idx].y);
      curr_score = CalcScore(curr_length);
      if (curr_score > best_score) {
        best_length = curr_length;
        best_score = curr_score;
        best_path = current_path;
        best_sel_pair_idx = sel_pair_idx;
        best_is_used = is_used;
      }
    }
    else {
      // 元に戻す
      ;
    }
  }

  // 最高スコアを戻す
  current_path = best_path;
  curr_length = best_length;
  curr_score = best_score;
  sel_pair_idx = best_sel_pair_idx;
  is_used = best_is_used;

  for (int i = 0; i < m; ++i) {
    partner_idx[sel_pair_idx[i]] = i + m;
    partner_idx[sel_pair_idx[i] + INIT_N] = i;
  }
  for (int i = 0; i < m; ++i) {
    sel_pair_idx[i + m] = sel_pair_idx[i] + INIT_N;
  }
  best_sel_pair_idx = sel_pair_idx;
  best_partner_idx = partner_idx;

  reduced_pair_idx = sel_pair_idx;
  best_reduced_pair_idx = reduced_pair_idx;
}

void sa_two_opt_path()
{
  double start_temp = 48;
  double end_temp = 0.0001;

  double time_limit = 1.8;

  int loop = 0;
  double now_time = get_elapsed_time();

  // それぞれをTSP
  for (int ui_tei = 0; ui_tei < 10; ++ui_tei) {
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
      }
      if (now_time > (time_limit / 15) * (6 + ui_tei)) break;

      int ite1 = Rand() % (m * 2);
      int ite2 = Rand() % (m * 2);
      while (ite1 == ite2) {
        ite2 = Rand() % (m * 2);
      }
      if (ite1 > ite2) swap(ite1, ite2);
      if (ite2 - ite1 == 1) continue;

      if (sel_pair_idx[ite1] + INIT_N == sel_pair_idx[ite2]) continue;
      if (sel_pair_idx[ite1] < INIT_N) {
        int pa = partner_idx[sel_pair_idx[ite1]];
        if (pa < ite2) {
          continue;
        }
      }
      if (sel_pair_idx[ite2] >= INIT_N) {
        int pa = partner_idx[sel_pair_idx[ite2]];
        if (ite1 < pa) {
          continue;
        }
      }

      int diff = 0;

      if (true || loop % 2 == 0) {
        diff += manhattan(current_path[ite1 + 1], current_path[ite2]);
        diff += manhattan(current_path[ite1 + 1], current_path[ite2 + 2]);
        diff -= manhattan(current_path[ite1 + 1], current_path[ite1]);
        diff -= manhattan(current_path[ite1 + 1], current_path[ite1 + 2]);
        diff += manhattan(current_path[ite2 + 1], current_path[ite1]);
        diff += manhattan(current_path[ite2 + 1], current_path[ite1 + 2]);
        diff -= manhattan(current_path[ite2 + 1], current_path[ite2]);
        diff -= manhattan(current_path[ite2 + 1], current_path[ite2 + 2]);

        int new_length = curr_length + diff;
        int diffScore = -diff;

        double temp = start_temp + (end_temp - start_temp) * now_time / time_limit;
        double prob = exp((double)diffScore / temp);
        if (prob > Rand01()) {
          curr_length = new_length;
          if (sel_pair_idx[ite1] < INIT_N) {
            partner_idx[sel_pair_idx[ite1] + INIT_N] = ite2;
          }
          else {
            partner_idx[sel_pair_idx[ite1] - INIT_N] = ite2;
          }
          if (sel_pair_idx[ite2] < INIT_N) {
            partner_idx[sel_pair_idx[ite2] + INIT_N] = ite1;
          }
          else {
            partner_idx[sel_pair_idx[ite2] - INIT_N] = ite1;
          }
          swap(sel_pair_idx[ite1], sel_pair_idx[ite2]);
          swap(current_path[ite1 + 1], current_path[ite2 + 1]);
          curr_score = CalcScore(curr_length);
          if (curr_score > best_score) {
            best_length = curr_length;
            best_score = curr_score;
            best_path = current_path;
            best_sel_pair_idx = sel_pair_idx;
            best_is_used = is_used;
            best_partner_idx = partner_idx;
          }
        }
        else {
          // 元に戻す
          ;
        }
      }
    }

    // 最高スコアを戻す
    current_path = best_path;
    curr_length = best_length;
    curr_score = best_score;
    sel_pair_idx = best_sel_pair_idx;
    reduced_pair_idx = best_reduced_pair_idx;
    is_used = best_is_used;
    partner_idx = best_partner_idx;
  }
}

void sa_path_pruning()
{
  double time_limit = 1.95;

  const int INF = 1001001001;
  int mi = INF;
  reduced_pair_idx = sel_pair_idx;
  vector<Point> ans2 = current_path;

  random_device seed_gen;
  mt19937 engine(seed_gen());

  vector<int> rand_order;
  for (int i = 0; i < m; ++i) {
    rand_order.push_back(i);
  }

  int loop = 0;
  double now_time = get_elapsed_time();
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
    }
    if (now_time > time_limit) break;

    std::shuffle(rand_order.begin(), rand_order.end(), engine);
    vector<int> pick_order;
    for (int i = 0; i < 50; ++i) pick_order.push_back(rand_order[i]);
    sort(pick_order.begin(), pick_order.end());

    set<int> pick_set;
    int cnt = 0;
    int now = 0;
    for (int i = 0; i < m * 2; ++i) {
      if (sel_pair_idx[i] < INIT_N) {
        if (cnt == pick_order[now]) {
          pick_set.insert(sel_pair_idx[i]);
          now++;
        }
        cnt++;
      }
    }

    int new_length = 0;
    int cur_x = 400;
    int cur_y = 400;
    vector<Point> ans3;
    ans3.push_back(Point(400, 400));
    for (int i = 0; i < m * 2; ++i) {
      int ite = sel_pair_idx[i];
      if (ite >= INIT_N) ite -= INIT_N;
      if (pick_set.find(ite) != pick_set.end()) {
        ans3.push_back(current_path[i + 1]);
        new_length += abs(cur_x - current_path[i + 1].x) + abs(cur_y - current_path[i + 1].y);
        cur_x = current_path[i + 1].x;
        cur_y = current_path[i + 1].y;
      }
    }
    ans3.push_back(Point(400, 400));

    if (new_length < mi) {
      mi = new_length;
      ans2 = ans3;
      reduced_pair_idx.clear();

      for (auto ite : pick_set) {
        reduced_pair_idx.push_back(ite);
      }
    }
  }

  current_path = ans2;
  sel_pair_idx = reduced_pair_idx;

  curr_length = compute_path_time();
  curr_score = CalcScore(curr_length);
}

int Solve(int case_num)
{
  start_timer();

  ResetParam();

  input_data(case_num);

  filter_input_points();

  build_initial_path();

  sa_point_swap();

  sa_two_opt_path();

  sa_path_pruning();

  output_data(case_num);

  cerr << curr_score << endl;
  cerr << get_elapsed_time() << "sec." << endl;

  return curr_score;
}

int main()
{
  exec_mode = 1;
  if (exec_mode == 0) {
    Solve(0);
  }
  else if (exec_mode == 1) {
    srep(i, 0, 10)
    {
      Solve(i);
    }
  }

  return 0;
}
