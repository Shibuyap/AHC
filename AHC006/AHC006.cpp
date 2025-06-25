#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <random>
#include <set>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

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

int mode = 0;

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

inline int manhattan(const Point& a, const Point& b)
{
  return abs(a.x - b.x) + abs(a.y - b.y);
}


const int INIT_N = 1000;

int n = 1000;
const int m = 51;
vector<Point> start_points(INIT_N), goal_points(INIT_N);
vector<int> orig_idx;

inline bool inside_rect(int ite, int L, int R, int U, int D)
{
  if (start_points[ite].x < L || R < start_points[ite].x) return false;
  if (start_points[ite].y < U || D < start_points[ite].y) return false;
  if (goal_points[ite].x < L || R < goal_points[ite].x) return false;
  if (goal_points[ite].y < U || D < goal_points[ite].y) return false;
  return true;
}

class Answer
{
public:
  int score;
  int length;
  vector<Point> path;

  vector<int> sel_idx;
  vector<int> red_idx;
  vector<int> used;
  vector<int> partner;

  Answer()
    : score(0), length(0)
  {
    sel_idx.resize(m * 2);
    red_idx.resize(m * 2);
    used.resize(INIT_N);
    partner.resize(INIT_N * 2);
  }

  void clear()
  {
    score = 0;
    length = 0;
    path.clear();
    fill(sel_idx.begin(), sel_idx.end(), INIT_N);
    fill(red_idx.begin(), red_idx.end(), INIT_N);
    fill(used.begin(), used.end(), 0);
    fill(partner.begin(), partner.end(), INIT_N);
  }

  int calc_score()
  {
    score = round(100000000.0 / (1000.0 + length));
    return score;
  }

  int compute_path_time()
  {
    length = 0;
    for (int i = 1; i < path.size(); ++i) {
      length += manhattan(path[i], path[i - 1]);
    }
    return length;
  }
};

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

void output_data(const Answer& answer, int case_num)
{
  if (mode == 0) {
    // 標準出力
    cout << 50;
    for (int i = 0; i < 50; ++i) {
      if (answer.sel_idx[i] < INIT_N) {
        cout << " " << orig_idx[answer.sel_idx[i]] + 1;
      }
    }
    cout << endl;
    cout << answer.path.size();
    for (int i = 0; i < answer.path.size(); ++i) {
      cout << " " << answer.path[i].x << " " << answer.path[i].y;
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
      if (answer.sel_idx[i] < INIT_N) {
        ofs << " " << orig_idx[answer.sel_idx[i]] + 1;
      }
    }
    ofs << endl;
    ofs << answer.path.size();
    for (int i = 0; i < answer.path.size(); ++i) {
      ofs << " " << answer.path[i].x << " " << answer.path[i].y;
    }
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

void reset_param()
{
  n = INIT_N;

  start_points.resize(INIT_N);
  goal_points.resize(INIT_N);
  orig_idx.clear();
}

void filter_input_points(double time_limit)
{
  // --- 収集済みの最良データを保持 ---
  vector<Point> best_s, best_g;
  int         best = 1'001'001;   // 今のところ最小の「矩形スコア」
  int         iter = 0;

  // --- 改善ループ ---
  while (true) {
    ++iter;
    if (iter % 100 == 1 && get_elapsed_time() > time_limit) { break; }

    // 候補データを一時的に保持
    vector<Point> cand_s, cand_g;
    vector<int> cand_idx;

    // ランダムに矩形を生成
    int l = Rand() % 801;
    int r = Rand() % 801;
    int u = Rand() % 801;
    int d = Rand() % 801;
    if (l > r)  swap(l, r);
    if (u > d)  swap(u, d);

    // “広さ” の代わりに幅＋高さを評価指標にしている
    const int metric = (r - l) + (d - u);
    if (metric >= best) { continue; }

    // 矩形内に収まる点を抽出
    for (int i = 0; i < n; ++i) {
      if (inside_rect(i, l, r, u, d)) {
        cand_s.push_back(start_points[i]);
        cand_g.push_back(goal_points[i]);
        cand_idx.push_back(i);
      }
    }

    // 十分な点数があり、しかも指標が更新されたら採用
    if (cand_s.size() >= m + 1) {
      best_s = std::move(cand_s);
      best_g = std::move(cand_g);
      orig_idx = std::move(cand_idx);
      best = metric;
    }
  }

  // ベスト結果でグローバルを更新
  if (best_s.size() >= m) {
    start_points.swap(best_s);
    goal_points.swap(best_g);
    n = static_cast<int>(start_points.size());
  }
}

void build_initial_path(Answer& answer)
{
  answer.clear();
  answer.path.push_back(Point(400, 400));
  for (int i = 0; i < m; ++i) {
    answer.path.push_back(Point(start_points[i].x, start_points[i].y));
    answer.sel_idx[i] = i;
    answer.used[i] = 1;
  }
  for (int i = 0; i < m; ++i) {
    answer.path.push_back(Point(goal_points[i].x, goal_points[i].y));
  }
  answer.path.push_back(Point(400, 400));
  answer.compute_path_time();
  answer.calc_score();
}

void sa_point_swap(Answer& answer, double time_limit)
{
  Answer best_answer = answer;

  // 山登り解、焼きなまし解
  double s_temp = 4800;
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
    while (answer.used[cand_idx]) {
      cand_idx = Rand() % n;
    }

    int path_pos = Rand() % m;

    int diff = 0;
    diff += manhattan(start_points[cand_idx], answer.path[path_pos]);
    diff += manhattan(start_points[cand_idx], answer.path[path_pos + 2]);
    diff += manhattan(goal_points[cand_idx], answer.path[path_pos + m]);
    diff += manhattan(goal_points[cand_idx], answer.path[path_pos + m + 2]);
    diff -= manhattan(answer.path[path_pos + 1], answer.path[path_pos]);
    diff -= manhattan(answer.path[path_pos + 1], answer.path[path_pos + 2]);
    diff -= manhattan(answer.path[path_pos + m + 1], answer.path[path_pos + m]);
    diff -= manhattan(answer.path[path_pos + m + 1], answer.path[path_pos + m + 2]);

    int new_length = answer.length + diff;
    int diffScore = -diff;

    double temp = s_temp + (end_temp - s_temp) * now_time / time_limit;
    double prob = exp((double)diffScore / temp);
    if (prob > Rand01()) {
      answer.length = new_length;
      answer.used[answer.sel_idx[path_pos]] = 0;
      answer.used[cand_idx] = 1;
      answer.sel_idx[path_pos] = cand_idx;
      answer.path[path_pos + 1] = Point(start_points[cand_idx].x, start_points[cand_idx].y);
      answer.path[path_pos + m + 1] = Point(goal_points[cand_idx].x, goal_points[cand_idx].y);
      answer.calc_score();
      if (answer.score > best_answer.score) {
        best_answer = answer;
      }
    }
    else {
      // 元に戻す
      ;
    }
  }

  // 最高スコアを戻す
  answer = best_answer;

  for (int i = 0; i < m; ++i) {
    answer.partner[answer.sel_idx[i]] = i + m;
    answer.partner[answer.sel_idx[i] + INIT_N] = i;
  }
  for (int i = 0; i < m; ++i) {
    answer.sel_idx[i + m] = answer.sel_idx[i] + INIT_N;
  }
  answer.red_idx = answer.sel_idx;

  best_answer = answer;

  cerr << "Point Swap iteration: " << loop << endl;
}

void sa_two_opt_path(Answer& answer, double time_limit)
{
  Answer best_answer = answer;

  double start_time = get_elapsed_time();
  double s_temp = 48;
  double end_temp = 0.0001;

  // それぞれをTSP
  for (int ui_tei = 0; ui_tei < 15; ++ui_tei) {
    double this_start_time = get_elapsed_time();
    double this_time_limit = start_time + (time_limit - start_time) / 15.0 * (ui_tei + 1);
    double now_time = get_elapsed_time();
    int loop = 0;

    while (true) {
      loop++;
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
      }
      if (now_time > this_time_limit) {
        break;
      }

      int ite1 = Rand() % (m * 2);
      int ite2 = Rand() % (m * 2);
      while (ite1 == ite2) {
        ite2 = Rand() % (m * 2);
      }
      if (ite1 > ite2) {
        swap(ite1, ite2);
      }
      if (ite2 - ite1 == 1) {
        continue;
      }

      if (answer.sel_idx[ite1] + INIT_N == answer.sel_idx[ite2]) {
        continue;
      }

      if (answer.sel_idx[ite1] < INIT_N) {
        int pa = answer.partner[answer.sel_idx[ite1]];
        if (pa < ite2) {
          continue;
        }
      }
      if (answer.sel_idx[ite2] >= INIT_N) {
        int pa = answer.partner[answer.sel_idx[ite2]];
        if (ite1 < pa) {
          continue;
        }
      }

      int diff = 0;
      diff += manhattan(answer.path[ite1 + 1], answer.path[ite2]);
      diff += manhattan(answer.path[ite1 + 1], answer.path[ite2 + 2]);
      diff -= manhattan(answer.path[ite1 + 1], answer.path[ite1]);
      diff -= manhattan(answer.path[ite1 + 1], answer.path[ite1 + 2]);
      diff += manhattan(answer.path[ite2 + 1], answer.path[ite1]);
      diff += manhattan(answer.path[ite2 + 1], answer.path[ite1 + 2]);
      diff -= manhattan(answer.path[ite2 + 1], answer.path[ite2]);
      diff -= manhattan(answer.path[ite2 + 1], answer.path[ite2 + 2]);

      int new_length = answer.length + diff;
      int diffScore = -diff;

      double temp = s_temp + (end_temp - s_temp) * (now_time - this_start_time) / (this_time_limit - this_start_time);
      double prob = exp((double)diffScore / temp);
      if (prob > Rand01()) {
        answer.length = new_length;
        if (answer.sel_idx[ite1] < INIT_N) {
          answer.partner[answer.sel_idx[ite1] + INIT_N] = ite2;
        }
        else {
          answer.partner[answer.sel_idx[ite1] - INIT_N] = ite2;
        }
        if (answer.sel_idx[ite2] < INIT_N) {
          answer.partner[answer.sel_idx[ite2] + INIT_N] = ite1;
        }
        else {
          answer.partner[answer.sel_idx[ite2] - INIT_N] = ite1;
        }
        swap(answer.sel_idx[ite1], answer.sel_idx[ite2]);
        swap(answer.path[ite1 + 1], answer.path[ite2 + 1]);
        answer.calc_score();
        if (answer.score > best_answer.score) {
          best_answer = answer;
        }
      }
      else {
        // 元に戻す
      }
    }

    cerr << "Set " << ui_tei << " iteration: " << loop << endl;

    // 最高スコアを戻す
    answer = best_answer;
  }
}

void sa_path_pruning(Answer& answer, double time_limit)
{
  const int INF = 1001001001;
  int mi = INF;
  answer.red_idx = answer.sel_idx;
  vector<Point> ans2 = answer.path;

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
    if (now_time > time_limit) { break; }

    std::shuffle(rand_order.begin(), rand_order.end(), engine);
    vector<int> pick_order;
    for (int i = 0; i < 50; ++i) pick_order.push_back(rand_order[i]);
    sort(pick_order.begin(), pick_order.end());

    set<int> pick_set;
    int cnt = 0;
    int now = 0;
    for (int i = 0; i < m * 2; ++i) {
      if (answer.sel_idx[i] < INIT_N) {
        if (cnt == pick_order[now]) {
          pick_set.insert(answer.sel_idx[i]);
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
      int ite = answer.sel_idx[i];
      if (ite >= INIT_N) ite -= INIT_N;
      if (pick_set.find(ite) != pick_set.end()) {
        ans3.push_back(answer.path[i + 1]);
        new_length += abs(cur_x - answer.path[i + 1].x) + abs(cur_y - answer.path[i + 1].y);
        cur_x = answer.path[i + 1].x;
        cur_y = answer.path[i + 1].y;
      }
    }
    ans3.push_back(Point(400, 400));

    if (new_length < mi) {
      mi = new_length;
      ans2 = ans3;
      answer.red_idx.clear();

      for (auto ite : pick_set) {
        answer.red_idx.push_back(ite);
      }
    }
  }

  answer.path = ans2;
  answer.sel_idx = answer.red_idx;

  answer.compute_path_time();
  answer.calc_score();

  cerr << "Path Pruning iteration: " << loop << endl;
}

int solve_case(int cn)
{
  start_timer();

  reset_param();

  input_data(cn);

  filter_input_points(0.2);

  Answer answer;

  build_initial_path(answer);

  sa_point_swap(answer, 0.4);

  sa_two_opt_path(answer, 1.8);

  sa_path_pruning(answer, 1.95);

  output_data(answer, cn);

  cerr << answer.score << endl;
  cerr << get_elapsed_time() << "sec." << endl;

  return answer.score;
}

int main()
{
  mode = 1;
  if (mode == 0) {
    solve_case(0);
  }
  else if (mode == 1) {
    for (int i = 0; i < 10; ++i) {
      solve_case(i);
    }
  }

  return 0;
}
