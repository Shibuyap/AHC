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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
#define MAX_N 200005

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

namespace /* 乱数ライブラリ */
{
  static uint32_t Rand() {
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


  static double Rand01() {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

int exec_mode = 0;

struct Point
{
  int x;
  int y;

  Point(int _x, int _y) {
    x = _x;
    y = _y;
  }
};

namespace /* 変数 */
{
  const int INIT_N = 1000;

  // 入力用変数
  int n = 1000;
  int m = 60;
  vector<int> start_x(INIT_N), start_y(INIT_N), goal_x(INIT_N), goal_y(INIT_N);
  vector<int> argA;

  // 解答用変数
  ll maxScore;
  ll minTime;
  vector<Point> current_path;
  vector<int> argAns(140);
  vector<int> argAns2(140);
  vector<int> use(1100);
  vector<int> pair_(2100);

  // 焼きなまし用変数
  ll real_maxScore;
  ll real_minTime;
  vector<Point> best_path;
  vector<int> real_argAns(140);
  vector<int> real_argAns2(140);
  vector<int> real_use(1100);
  vector<int> real_pair_(2100);

  void ResetParam() {
    argA.clear();
    maxScore = 0;
    minTime = 0;
    current_path.clear();
    argAns.resize(140);
    argAns2.resize(140);
    use.resize(1100);
    pair_.resize(2100);
    real_maxScore = 0;
    real_minTime = 0;
    best_path.clear();
  }

}  // namespace

// スコア計算
ll CalcScore(ll time) {
  ll res = round(100000000.0 / (1000.0 + time));
  return res;
}

ll compute_path_time() {
  ll timeSum = 0;
  srep(i, 1, current_path.size()) {
    timeSum += abs(current_path[i].x - current_path[i - 1].x) + abs(current_path[i].y - current_path[i - 1].y);
  }
  return timeSum;
}

inline bool is_inside_rect(int ite, int L, int R, int U, int D) {
  if (start_x[ite] < L || R < start_x[ite]) return false;
  if (start_y[ite] < U || D < start_y[ite]) return false;
  if (goal_x[ite] < L || R < goal_x[ite]) return false;
  if (goal_y[ite] < U || D < goal_y[ite]) return false;
  return true;
}

void filter_input_points() {
  // --- 収集済みの最良データを保持 ---
  vector<int> best_start_x, best_start_y, best_goal_x, best_goal_y;
  int         best_metric = 1'001'001;   // 今のところ最小の「矩形スコア」
  int         iter_count = 0;

  // --- 改善ループ ---
  while (true) {
    ++iter_count;
    if (iter_count % 100 == 1 && get_elapsed_time() > 0.2) break;

    // 候補データを一時的に保持
    vector<int> cand_start_x, cand_start_y, cand_goal_x, cand_goal_y, cand_indices;

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
    rep(i, n) {
      if (is_inside_rect(i, rect_left, rect_right, rect_upper, rect_lower)) {
        cand_start_x.push_back(start_x[i]);
        cand_start_y.push_back(start_y[i]);
        cand_goal_x.push_back(goal_x[i]);
        cand_goal_y.push_back(goal_y[i]);
        cand_indices.push_back(i);
      }
    }

    // 十分な点数があり、しかも指標が更新されたら採用
    if (cand_start_x.size() >= m + 1) {
      best_start_x = std::move(cand_start_x);
      best_start_y = std::move(cand_start_y);
      best_goal_x = std::move(cand_goal_x);
      best_goal_y = std::move(cand_goal_y);
      argA = std::move(cand_indices);
      best_metric = rect_metric;
    }
  }

  // ベスト結果でグローバルを更新
  if (best_start_x.size() >= m) {
    start_x.swap(best_start_x);
    start_y.swap(best_start_y);
    goal_x.swap(best_goal_x);
    goal_y.swap(best_goal_y);
    n = static_cast<int>(start_x.size());
  }
}

void input_data(int case_num) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    rep(i, n) {
      cin >> start_x[i] >> start_y[i] >> goal_x[i] >> goal_y[i];
    }
  }
  else {
    // ファイル入力
    rep(i, n) {
      ifs >> start_x[i] >> start_y[i] >> goal_x[i] >> goal_y[i];
    }
  }
}

void output_data(int case_num) {
  if (exec_mode == 0) {
    // 標準出力
    cout << 50;
    rep(i, 50) {
      if (argAns[i] < INIT_N) {
        cout << " " << argA[argAns[i]] + 1;
      }
    }
    cout << endl;
    cout << current_path.size();
    rep(i, current_path.size()) {
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
    rep(i, 50) {
      if (argAns[i] < INIT_N) {
        ofs << " " << argA[argAns[i]] + 1;
      }
    }
    ofs << endl;
    ofs << current_path.size();
    rep(i, current_path.size()) {
      ofs << " " << current_path[i].x << " " << current_path[i].y;
    }
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

inline int manhattan(const Point& a, const Point& b) {
  return abs(a.x - b.x) + abs(a.y - b.y);
}

int Solve(int mode) {
  start_timer();

  input_data(0);

  filter_input_points();

  // 愚直解
  current_path.push_back(Point(400, 400));
  rep(i, m) {
    current_path.push_back(Point(start_x[i], start_y[i]));
    argAns[i] = i;
    use[i] = 1;
  }
  rep(i, m) {
    current_path.push_back(Point(goal_x[i], goal_y[i]));
  }
  current_path.push_back(Point(400, 400));
  minTime = compute_path_time();
  maxScore = CalcScore(minTime);

  best_path = current_path;
  real_minTime = minTime;
  real_maxScore = maxScore;

  // 山登り解、焼きなまし解
  double TL = 1.8;
  double start_temp = 48;
  double end_temp = 0.0001;
  int loop = 0;
  double now_time = get_elapsed_time();
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
      if (now_time > 0.4) {
        break;
      }
    }

    int x = Rand() % n;
    while (use[x]) {
      x = Rand() % n;
    }

    int ite = Rand() % m;

    int keepNum = argAns[ite];
    int diff = 0;
    diff += abs(start_x[x] - current_path[ite].x) + abs(start_y[x] - current_path[ite].y);
    diff += abs(start_x[x] - current_path[ite + 2].x) + abs(start_y[x] - current_path[ite + 2].y);
    diff += abs(goal_x[x] - current_path[ite + m].x) + abs(goal_y[x] - current_path[ite + m].y);
    diff += abs(goal_x[x] - current_path[ite + m + 2].x) + abs(goal_y[x] - current_path[ite + m + 2].y);

    diff -= abs(current_path[ite + 1].x - current_path[ite].x) + abs(current_path[ite + 1].y - current_path[ite].y);
    diff -= abs(current_path[ite + 1].x - current_path[ite + 2].x) + abs(current_path[ite + 1].y - current_path[ite + 2].y);
    diff -= abs(current_path[ite + m + 1].x - current_path[ite + m].x) + abs(current_path[ite + m + 1].y - current_path[ite + m].y);
    diff -= abs(current_path[ite + m + 1].x - current_path[ite + m + 2].x) + abs(current_path[ite + m + 1].y - current_path[ite + m + 2].y);

    int tmpTime = minTime + diff;
    int diffScore = -diff;

    double temp = start_temp + (end_temp - start_temp) * now_time / TL;
    double prob = exp((double)diffScore / temp);
    if (prob > Rand01()) {
      minTime = tmpTime;
      use[argAns[ite]] = 0;
      use[x] = 1;
      argAns[ite] = x;
      current_path[ite + 1] = Point(start_x[x], start_y[x]);
      current_path[ite + m + 1] = Point(goal_x[x], goal_y[x]);
      maxScore = CalcScore(minTime);
      if (maxScore > real_maxScore) {
        real_minTime = minTime;
        real_maxScore = maxScore;
        best_path = current_path;
        real_argAns = argAns;
        real_use = use;
      }
    }
    else {
      // 元に戻す
      ;
    }
  }

  // 最高スコアを戻す
  current_path = best_path;
  minTime = real_minTime;
  maxScore = real_maxScore;
  argAns = real_argAns;
  use = real_use;

  rep(i, m) {
    pair_[argAns[i]] = i + m;
    pair_[argAns[i] + INIT_N] = i;
  }
  rep(i, m) {
    argAns[i + m] = argAns[i] + INIT_N;
  }
  real_argAns = argAns;
  real_pair_ = pair_;

  argAns2 = argAns;
  real_argAns2 = argAns2;

  // それぞれをTSP
  rep(ui_tei, 10) {
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
      }
      if (now_time > (TL / 15) * (6 + ui_tei)) break;

      int ite1 = Rand() % (m * 2);
      int ite2 = Rand() % (m * 2);
      while (ite1 == ite2) {
        ite2 = Rand() % (m * 2);
      }
      if (ite1 > ite2) swap(ite1, ite2);
      if (ite2 - ite1 == 1) continue;

      if (argAns[ite1] + INIT_N == argAns[ite2]) continue;
      if (argAns[ite1] < INIT_N) {
        int pa = pair_[argAns[ite1]];
        if (pa < ite2) {
          continue;
        }
      }
      if (argAns[ite2] >= INIT_N) {
        int pa = pair_[argAns[ite2]];
        if (ite1 < pa) {
          continue;
        }
      }

      int diff = 0;

      if (true || loop % 2 == 0) {
        diff += abs(current_path[ite1 + 1].x - current_path[ite2].x) + abs(current_path[ite1 + 1].y - current_path[ite2].y);
        diff += abs(current_path[ite1 + 1].x - current_path[ite2 + 2].x) + abs(current_path[ite1 + 1].y - current_path[ite2 + 2].y);
        diff -= abs(current_path[ite1 + 1].x - current_path[ite1].x) + abs(current_path[ite1 + 1].y - current_path[ite1].y);
        diff -= abs(current_path[ite1 + 1].x - current_path[ite1 + 2].x) + abs(current_path[ite1 + 1].y - current_path[ite1 + 2].y);

        diff += abs(current_path[ite2 + 1].x - current_path[ite1].x) + abs(current_path[ite2 + 1].y - current_path[ite1].y);
        diff += abs(current_path[ite2 + 1].x - current_path[ite1 + 2].x) + abs(current_path[ite2 + 1].y - current_path[ite1 + 2].y);
        diff -= abs(current_path[ite2 + 1].x - current_path[ite2].x) + abs(current_path[ite2 + 1].y - current_path[ite2].y);
        diff -= abs(current_path[ite2 + 1].x - current_path[ite2 + 2].x) + abs(current_path[ite2 + 1].y - current_path[ite2 + 2].y);

        int tmpTime = minTime + diff;
        int diffScore = -diff;

        double temp = start_temp + (end_temp - start_temp) * now_time / TL;
        double prob = exp((double)diffScore / temp);
        if (prob > Rand01()) {
          minTime = tmpTime;
          if (argAns[ite1] < INIT_N) {
            pair_[argAns[ite1] + INIT_N] = ite2;
          }
          else {
            pair_[argAns[ite1] - INIT_N] = ite2;
          }
          if (argAns[ite2] < INIT_N) {
            pair_[argAns[ite2] + INIT_N] = ite1;
          }
          else {
            pair_[argAns[ite2] - INIT_N] = ite1;
          }
          swap(argAns[ite1], argAns[ite2]);
          swap(current_path[ite1 + 1], current_path[ite2 + 1]);
          maxScore = CalcScore(minTime);
          if (maxScore > real_maxScore) {
            real_minTime = minTime;
            real_maxScore = maxScore;
            best_path = current_path;
            real_argAns = argAns;
            real_use = use;
            real_pair_ = pair_;
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
    minTime = real_minTime;
    maxScore = real_maxScore;
    argAns = real_argAns;
    argAns2 = real_argAns2;
    use = real_use;
    pair_ = real_pair_;
  }

  TL = 1.95;

  const int INF = 1001001001;
  int mi = INF;
  argAns2 = argAns;
  vector<Point> ans2 = current_path;

  random_device seed_gen;
  mt19937 engine(seed_gen());

  vector<int> v;
  rep(i, m) {
    v.push_back(i);
  }

  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
    }
    if (now_time > TL) break;

    std::shuffle(v.begin(), v.end(), engine);
    vector<int> vv;
    rep(i, 50) vv.push_back(v[i]);
    sort(vv.begin(), vv.end());

    set<int> useIte;
    int cnt = 0;
    int now = 0;
    rep(i, m * 2) {
      if (argAns[i] < INIT_N) {
        if (cnt == vv[now]) {
          useIte.insert(argAns[i]);
          now++;
        }
        cnt++;
      }
    }

    int tmpTime = 0;
    int nowx = 400;
    int nowy = 400;
    vector<Point> ans3;
    ans3.push_back(Point(400, 400));
    rep(i, m * 2) {
      int ite = argAns[i];
      if (ite >= INIT_N) ite -= INIT_N;
      if (useIte.find(ite) != useIte.end()) {
        ans3.push_back(current_path[i + 1]);
        tmpTime += abs(nowx - current_path[i + 1].x) + abs(nowy - current_path[i + 1].y);
        nowx = current_path[i + 1].x;
        nowy = current_path[i + 1].y;
      }
    }
    ans3.push_back(Point(400, 400));

    if (tmpTime < mi) {
      mi = tmpTime;
      ans2 = ans3;
      argAns2.clear();

      for (auto ite : useIte) {
        argAns2.push_back(ite);
      }
    }
  }

  current_path = ans2;
  argAns = argAns2;

  minTime = compute_path_time();
  maxScore = CalcScore(minTime);

  output_data(0);

  // デバッグ用
  if (mode != 0) {
    cout << "loop = " << loop << endl;
    cout << maxScore << endl;
    cout << get_elapsed_time() << "sec." << endl;
  }

  return maxScore;
}

int main() {
  exec_mode = 1;
  m = 55;
  if (exec_mode == 0) {
    Solve(exec_mode);
  }
  else if (exec_mode == 1) {
    Solve(exec_mode);
  }

  return 0;
}
