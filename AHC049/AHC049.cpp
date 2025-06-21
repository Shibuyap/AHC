#include <algorithm>
#include <array>
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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

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

// 乱数
namespace
{
  static uint32_t rand_xorshift()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
  }

  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U', 'L', 'D', 'R' };

const double TIME_LIMIT = 1.9;
int exec_mode;

const int MAX_N = 20;

struct Box
{
  int w;  // 重さ
  int d;  // 耐久力
  int initial_d;  // 初期耐久力
  bool exist;  // 箱がまだ存在するか

  Box() : w(0), d(0), initial_d(0), exist(false) {}
  Box(int weight, int durability) : w(weight), d(durability), initial_d(durability), exist(true) {}
};

struct Input
{
  int n;
  Box grid[MAX_N][MAX_N];  // グリッド上の箱

  void read(istream& is)
  {
    is >> n;

    // 重さを読み込み
    rep(i, n)
    {
      rep(j, n)
      {
        int w;
        is >> w;
        grid[i][j].w = w;
        grid[i][j].exist = (w > 0);
      }
    }

    // 耐久力を読み込み
    rep(i, n)
    {
      rep(j, n)
      {
        int d;
        is >> d;
        grid[i][j].d = d;
        grid[i][j].initial_d = d;
      }
    }
  }
};

vector<P> pick_order;
vector<P> best_pick_order;
ll best_pick_order_score;

struct State
{
  int x, y;  // 現在位置
  vector<Box> hand;  // 手に持っている箱（下から順）
  Box grid[MAX_N][MAX_N];  // グリッド上の箱
  int remaining_boxes;  // 残っている箱の数
  int move_count;  // 移動回数
  vector<char> actions;  // 実行した操作の履歴

  State() : x(0), y(0), remaining_boxes(0), move_count(0) {}

  void init(const Input& input)
  {
    x = 0;
    y = 0;
    hand.clear();
    remaining_boxes = 0;
    move_count = 0;
    actions.clear();

    rep(i, input.n)
    {
      rep(j, input.n)
      {
        grid[i][j] = input.grid[i][j];
        if (grid[i][j].exist && !(i == 0 && j == 0)) {
          remaining_boxes++;
        }
      }
    }
  }

  bool pickup()
  {
    if (!grid[x][y].exist) return false;
    if (x == 0 && y == 0) return false;  // 出入り口では拾えない

    hand.push_back(grid[x][y]);
    grid[x][y].exist = false;
    remaining_boxes--;
    actions.push_back('1');
    return true;
  }

  bool put()
  {
    if (hand.empty()) return false;
    if (grid[x][y].exist) return false;
    if (x == 0 && y == 0) return false;  // 出入り口には置けない

    grid[x][y] = hand.back();
    hand.pop_back();
    remaining_boxes++;
    actions.push_back('2');
    return true;
  }

  bool move(int dir)
  {
    int nx = x + DX[dir];
    int ny = y + DY[dir];

    if (nx < 0 || nx >= MAX_N || ny < 0 || ny >= MAX_N) return false;

    x = nx;
    y = ny;
    move_count++;

    // 手に持っている箱の耐久力を減少
    int total_weight = 0;
    for (int i = hand.size() - 1; i >= 0; i--) {
      hand[i].d -= total_weight;
      total_weight += hand[i].w;
    }

    actions.push_back(DIR_CHAR[dir]);

    // 出入り口に到着したら箱を全て運び出す
    if (x == 0 && y == 0 && !hand.empty()) {
      hand.clear();
    }

    return true;
  }

  ll calc_score(int n) const
  {
    if (remaining_boxes > 0) {
      return n * n - remaining_boxes;
    }
    else {
      return n * n + 2 * n * n * n - move_count;
    }
  }

  // 現在位置から(0,0)への長方形内で安全に運べる箱の一覧を取得
  vector<P> get_safe_pickable_boxes(int min_d) const
  {
    vector<P> pickable;

    // (0,0)と(x,y)で作る長方形の範囲を決定
    int min_x = min(0, x);
    int max_x = max(0, x);
    int min_y = min(0, y);
    int max_y = max(0, y);

    // 長方形内の各マスをチェック
    for (int i = min_x; i <= max_x; i++) {
      for (int j = min_y; j <= max_y; j++) {
        if (!grid[i][j].exist) continue;
        if (i == 0 && j == 0) continue;  // 出入り口はスキップ

        // その箱を拾った場合の(0,0)までの最短距離
        int dist_to_home = i + j;

        // その箱が耐えられるかチェック
        // 箱の耐久力 > (0,0)までの距離 × 現在の総重量
        if (min_d > grid[i][j].w * dist_to_home) {
          pickable.push_back(P(i, j));
        }
      }
    }

    return pickable;
  }
};

Input input;
State current_state;
State best_state;

void store_best_state()
{
  best_state = current_state;
}

void restore_best_state()
{
  current_state = best_state;
}

bool is_out_of_range(int x, int y)
{
  if (x < 0 || MAX_N <= x || y < 0 || MAX_N <= y) return true;
  return false;
}

void initialize_state()
{
  current_state.init(input);
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    input.read(cin);
  }
  else {
    // ファイル入力
    input.read(ifs);
  }
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

void output_data(ofstream& ofs)
{
  const vector<char>& actions = best_state.actions;

  if (exec_mode == 0) {
    // 標準出力
    for (char c : actions) {
      cout << c << endl;
    }
  }
  else {
    // ファイル出力
    for (char c : actions) {
      ofs << c << endl;
    }
  }
}

ll calculate_score()
{
  return best_state.calc_score(input.n);
}

void build_initial_solution()
{
  current_state.init(input);

  // 複数箱を同時運搬する貪欲解法
  while (current_state.remaining_boxes > 0) {
    // 一番重い箱
    int max_w = -1;
    int target_x = -1, target_y = -1;

    for (int i = 0; i < input.n; i++) {
      for (int j = 0; j < input.n; j++) {
        if (current_state.grid[i][j].exist) {
          if (current_state.grid[i][j].w > max_w) {
            max_w = current_state.grid[i][j].w;
            target_x = i;
            target_y = j;
          }
        }
      }
    }

    if (target_x == -1) break;

    // 目標位置への移動経路上で箱を拾う
    // まずx方向に移動
    while (current_state.x != target_x) {
      if (current_state.x < target_x) {
        current_state.move(2);  // Down
      }
      else {
        current_state.move(0);  // Up
      }
    }

    // 次にy方向に移動
    while (current_state.y != target_y) {
      if (current_state.y < target_y) {
        current_state.move(3);  // Right
      }
      else {
        current_state.move(1);  // Left
      }
    }

    int min_d = INT_MAX;

    // 目標位置の箱を拾う
    if (current_state.grid[current_state.x][current_state.y].exist) {
      min_d = current_state.grid[current_state.x][current_state.y].d;
      current_state.pickup();
      pick_order.push_back(P(current_state.x, current_state.y));
    }

    // 出入り口へ戻る（途中で箱を拾いながら）
    while (current_state.x != 0 || current_state.y != 0) {
      auto pickable_boxes = current_state.get_safe_pickable_boxes(min_d);

      if (!pickable_boxes.empty()) {
        // 安全に拾える箱がある場合、最も近い箱を拾う
        P closest_box = pickable_boxes[0];
        int min_dist = closest_box.first + closest_box.second;
        for (const auto& box : pickable_boxes) {
          int dist = box.first + box.second;
          if (dist > min_dist) {
            closest_box = box;
            min_dist = dist;
          }
        }
        while (current_state.x != closest_box.first || current_state.y != closest_box.second) {
          if (current_state.x < closest_box.first) {
            current_state.move(2);  // Down
          }
          else if (current_state.x > closest_box.first) {
            current_state.move(0);  // Up
          }
          else if (current_state.y < closest_box.second) {
            current_state.move(3);  // Right
          }
          else if (current_state.y > closest_box.second) {
            current_state.move(1);  // Left
          }
        }
        current_state.pickup();
        pick_order.push_back(P(current_state.x, current_state.y));
        min_d -= (closest_box.first + closest_box.second) * current_state.grid[current_state.x][current_state.y].w;
        min_d = min(min_d, current_state.grid[current_state.x][current_state.y].d);
      }
      else {
        while (current_state.x != 0 || current_state.y != 0) {
          // 安全に拾える箱がない場合、出入り口へ戻る
          if (current_state.x > 0) {
            current_state.move(0);  // Up
          }
          else if (current_state.y > 0) {
            current_state.move(1);  // Left
          }
        }
      }
    }
  }

  best_state = current_state;
}

int calc_score_with_pick_order()
{
  int move_count = 0;
  int min_d = INT_MAX;
  int x = 0, y = 0;

  for (auto& p : pick_order) {
    int target_x = p.first;
    int target_y = p.second;

    if (target_x > x || target_y > y) {
      // まず(0,0)に行く
      move_count += x + y;  // (0,0)までの距離
      x = 0;
      y = 0;
      min_d = INT_MAX;
      // 次に目的地へ移動
      move_count += target_x + target_y;  // (target_x, target_y)までの距離
      x = target_x;
      y = target_y;
    }
    else {
      if (min_d > (target_x + target_y) * current_state.grid[target_x][target_y].w) {
        // 目的地へ直接移動
        move_count += (x - target_x) + (y - target_y);
        x = target_x;
        y = target_y;
      }
      else {
        // まず(0,0)に行く
        move_count += x + y;  // (0,0)までの距離
        x = 0;
        y = 0;
        min_d = INT_MAX;
        // 次に目的地へ移動
        move_count += target_x + target_y;  // (target_x, target_y)までの距離
        x = target_x;
        y = target_y;
      }
    }

    min_d -= (target_x + target_y) * current_state.grid[x][y].w;
    min_d = min(min_d, current_state.grid[x][y].d);
  }

  // 0,0に移動
  move_count += x + y;

  int n = 20;
  int score = n * n + 2 * n * n * n - move_count;
  return score;
}

vector<vector<P>> calc_score_with_pick_order2()
{
  vector<vector<P>> res;
  int move_count = 0;
  int min_d = INT_MAX;
  int x = 0, y = 0;

  vector<P> now;
  for (auto& p : pick_order) {
    int target_x = p.first;
    int target_y = p.second;

    if (target_x > x || target_y > y) {
      // まず(0,0)に行く
      move_count += x + y;  // (0,0)までの距離
      x = 0;
      y = 0;
      min_d = INT_MAX;
      res.push_back(now);
      now.clear();
      // 次に目的地へ移動
      move_count += target_x + target_y;  // (target_x, target_y)までの距離
      x = target_x;
      y = target_y;
    }
    else {
      if (min_d > (target_x + target_y) * current_state.grid[target_x][target_y].w) {
        // 目的地へ直接移動
        move_count += (x - target_x) + (y - target_y);
        x = target_x;
        y = target_y;
      }
      else {
        // まず(0,0)に行く
        move_count += x + y;  // (0,0)までの距離
        x = 0;
        y = 0;
        min_d = INT_MAX;
        res.push_back(now);
        now.clear();
        // 次に目的地へ移動
        move_count += target_x + target_y;  // (target_x, target_y)までの距離
        x = target_x;
        y = target_y;
      }
    }

    min_d -= (target_x + target_y) * current_state.grid[x][y].w;
    min_d = min(min_d, current_state.grid[x][y].d);
    now.push_back(P(target_x, target_y));
  }

  // 0,0に移動
  move_count += x + y;
  res.push_back(now);
  now.clear();

  int n = 20;
  int score = n * n + 2 * n * n * n - move_count;
  return res;
}

vector<int> modoru;
void remake_state(bool modoraaa = false)
{
  current_state.init(input);
  int min_d = INT_MAX;
  rep(i, pick_order.size())
  {
    auto p = pick_order[i];
    int target_x = p.first;
    int target_y = p.second;

    if ((target_x > current_state.x || target_y > current_state.y) && (!modoraaa || !modoru[i])) {
      // まず(0,0)に行く
      while (current_state.x != 0 || current_state.y != 0) {
        if (current_state.x > 0) {
          current_state.move(0);  // Up
        }
        else if (current_state.y > 0) {
          current_state.move(1);  // Left
        }
      }
      min_d = INT_MAX;
      // 次に目的地へ移動
      while (current_state.x != target_x || current_state.y != target_y) {
        if (current_state.x < target_x) {
          current_state.move(2);  // Down
        }
        else if (current_state.x > target_x) {
          current_state.move(0);  // Up
        }
        else if (current_state.y < target_y) {
          current_state.move(3);  // Right
        }
        else if (current_state.y > target_y) {
          current_state.move(1);  // Left
        }
      }
    }
    else {
      if (min_d > (target_x + target_y) * current_state.grid[target_x][target_y].w) {
        // 目的地へ直接移動
        while (current_state.x != target_x || current_state.y != target_y) {
          if (current_state.x < target_x) {
            current_state.move(2);  // Down
          }
          else if (current_state.x > target_x) {
            current_state.move(0);  // Up
          }
          else if (current_state.y < target_y) {
            current_state.move(3);  // Right
          }
          else if (current_state.y > target_y) {
            current_state.move(1);  // Left
          }
        }
      }
      else {
        // まず(0,0)に行く
        while (current_state.x != 0 || current_state.y != 0) {
          if (current_state.x > 0) {
            current_state.move(0);  // Up
          }
          else if (current_state.y > 0) {
            current_state.move(1);  // Left
          }
        }
        min_d = INT_MAX;
        // 次に目的地へ移動
        while (current_state.x != target_x || current_state.y != target_y) {
          if (current_state.x < target_x) {
            current_state.move(2);  // Down
          }
          else if (current_state.x > target_x) {
            current_state.move(0);  // Up
          }
          else if (current_state.y < target_y) {
            current_state.move(3);  // Right
          }
          else if (current_state.y > target_y) {
            current_state.move(1);  // Left
          }
        }
      }
    }

    current_state.pickup();
    min_d -= (target_x + target_y) * current_state.grid[current_state.x][current_state.y].w;
    min_d = min(min_d, current_state.grid[current_state.x][current_state.y].d);
  }

  while (current_state.x != 0 || current_state.y != 0) {
    if (current_state.x > 0) {
      current_state.move(0);  // Up
    }
    else if (current_state.y > 0) {
      current_state.move(1);  // Left
    }
  }
}

void next_move(int& x, int& y, int& min_d, int& move_count, int target_x, int target_y, bool aaa = true)
{
  if (aaa && (target_x > x || target_y > y)) {
    // まず(0,0)に行く
    move_count += x + y;  // (0,0)までの距離
    x = 0;
    y = 0;
    min_d = INT_MAX;
    // 次に目的地へ移動
    move_count += target_x + target_y;  // (target_x, target_y)までの距離
    x = target_x;
    y = target_y;
  }
  else {
    if (min_d > (target_x + target_y) * current_state.grid[target_x][target_y].w) {
      // 目的地へ直接移動
      move_count += (x - target_x) + (y - target_y);
      x = target_x;
      y = target_y;
    }
    else {
      // まず(0,0)に行く
      move_count += x + y;  // (0,0)までの距離
      x = 0;
      y = 0;
      min_d = INT_MAX;
      // 次に目的地へ移動
      move_count += target_x + target_y;  // (target_x, target_y)までの距離
      x = target_x;
      y = target_y;
    }
  }

  min_d -= (target_x + target_y) * current_state.grid[x][y].w;
  min_d = min(min_d, current_state.grid[x][y].d);
}

int b[MAX_N][MAX_N];
vector<P> get_safe_pickable_boxes(int x, int y, int min_d)
{
  vector<P> pickable;

  // (0,0)と(x,y)で作る長方形の範囲を決定
  int min_x = min(0, x);
  int max_x = max(0, x);
  int min_y = min(0, y);
  int max_y = max(0, y);

  // 長方形内の各マスをチェック
  for (int i = min_x; i <= max_x; i++) {
    for (int j = min_y; j <= max_y; j++) {
      if (b[i][j] == 0) continue;
      if (i == 0 && j == 0) continue;  // 出入り口はスキップ

      // その箱を拾った場合の(0,0)までの最短距離
      int dist_to_home = i + j;

      // その箱が耐えられるかチェック
      // 箱の耐久力 > (0,0)までの距離 × 現在の総重量
      if (min_d > current_state.grid[i][j].w * dist_to_home) {
        pickable.push_back(P(i, j));
      }
    }
  }

  return pickable;
}

inline int Min(int a, int b)
{
  if (a <= b)return a;
  return b;
}

vector<P> after_method()
{
  auto vv = calc_score_with_pick_order2();

  vector<P> new_pick_order;

  int move_count = 0;
  int min_d = INT_MAX;
  int x = 0, y = 0;

  rep(i, MAX_N)
  {
    rep(j, MAX_N)
    {
      b[i][j] = 1;
    }
  }
  b[0][0] = 0;

  for (int i = 0; i < vv.size(); i++) {
    auto& v = vv[i];
    for (int j = 0; j < v.size(); j++) {
      if (b[v[j].first][v[j].second] == 0) continue;  // 既に拾われた箱はスキップ
      next_move(x, y, min_d, move_count, v[j].first, v[j].second);
      new_pick_order.push_back(P(x, y));
      b[x][y] = 0;

      bool isLast = true;
      srep(k, j + 1, v.size())
      {
        if (b[v[k].first][v[k].second] != 0) {
          isLast = false;
          break;
        }
      }

      if (isLast) {
        // 出入り口へ戻る（途中で箱を拾いながら）
        while (true) {
          auto pickable_boxes = get_safe_pickable_boxes(x, y, min_d);

          if (!pickable_boxes.empty()) {
            // 安全に拾える箱がある場合、最も近い箱を拾う
            P closest_box = pickable_boxes[0];

            int after_mind = 0;
            int min_dist = closest_box.first + closest_box.second;
            for (const auto& box : pickable_boxes) {
              int dist = box.first + box.second;
              int tmp_min_d = min(min_d - (box.first + box.second) * current_state.grid[box.first][box.second].w, current_state.grid[box.first][box.second].d);
              //if (tmp_min_d > after_mind) {
              if (dist > min_dist) {
                closest_box = box;
                after_mind = tmp_min_d;
                min_dist = dist;
              }
            }

            next_move(x, y, min_d, move_count, closest_box.first, closest_box.second);
            new_pick_order.push_back(P(x, y));
            b[x][y] = 0;
          }
          else {
            break;
          }
        }
      }
    }
  }

  return new_pick_order;
}

vector<P> bbb;
int simulate_need(vector<P>& v, int x, int y, int move_count)
{
  int need_d = 0;
  bbb.clear();
  int min_d = INT_MAX;
  for (int j = 0; j < v.size(); j++) {
    if (b[v[j].first][v[j].second] == 0) continue;  // 既に拾われた箱はスキップ
    next_move(x, y, min_d, move_count, v[j].first, v[j].second);
    b[x][y] = 0;
    bbb.push_back(P(x, y));
    need_d += (x + y) * current_state.grid[x][y].w;  // (x,y)までの距離 × 重さ

    bool isLast = true;
    srep(k, j + 1, v.size())
    {
      if (b[v[k].first][v[k].second] != 0) {
        isLast = false;
        break;
      }
    }

    if (isLast) {
      // 出入り口へ戻る（途中で箱を拾いながら）
      while (true) {
        auto pickable_boxes = get_safe_pickable_boxes(x, y, min_d);

        if (!pickable_boxes.empty()) {
          // 安全に拾える箱がある場合、最も近い箱を拾う
          P closest_box = pickable_boxes[0];

          int after_mind = 0;
          int min_dist = closest_box.first + closest_box.second;
          for (const auto& box : pickable_boxes) {
            int dist = box.first + box.second;
            int tmp_min_d = min(min_d - (box.first + box.second) * current_state.grid[box.first][box.second].w, current_state.grid[box.first][box.second].d);
            //if (tmp_min_d > after_mind) {
            if (dist > min_dist) {
              closest_box = box;
              after_mind = tmp_min_d;
              min_dist = dist;
            }
          }

          next_move(x, y, min_d, move_count, closest_box.first, closest_box.second);
          b[x][y] = 0;
          bbb.push_back(P(x, y));
          need_d += (x + y) * current_state.grid[x][y].w;  // (x,y)までの距離 × 重さ
        }
        else {
          break;
        }
      }
    }
  }

  for(auto& p : bbb) {
    b[p.first][p.second] = 1;  // 出入り口は拾えないので、ここで拾う
  }

  return need_d;
}

P can_pick(int gx, int gy, int need_d)
{
  P p(-1, -1);
  int max_mul = 0;
  rep(i, gx)
  {
    rep(j, gy)
    {
      if (b[i][j] == 0) continue;  // 既に拾われた箱はスキップ
      if (current_state.grid[i][j].d > need_d) {
        int tmp = (i + j) * current_state.grid[i][j].w;  // (i,j)までの距離 × 重さ
        if (tmp > max_mul) {
          max_mul = tmp;
          p = P(i, j);
        }
      }
    }
  }
  return p;
}

vector<P> after_method2(bool aaa = false)
{
  auto vv = calc_score_with_pick_order2();

  vector<P> new_pick_order;

  int move_count = 0;
  int min_d = INT_MAX;
  int x = 0, y = 0;

  rep(i, MAX_N)
  {
    rep(j, MAX_N)
    {
      b[i][j] = 1;
    }
  }
  b[0][0] = 0;

  for (int i = 0; i < vv.size(); i++) {
    auto& v = vv[i];
    int need_d = simulate_need(v, x, y, move_count);

    int first_x = -1;
    int first_y = -1;  
    for (int j = 0; j < v.size(); j++) {
      if (b[v[j].first][v[j].second] != 0) {
        first_x = v[j].first;
        first_y = v[j].second;
        break;
      }
    }

    P p(-1, -1);
    if(first_x != -1 && first_y != -1) {
      p = can_pick(first_x, first_y, need_d);
    }
    if(p.first != -1) {
      next_move(x, y, min_d, move_count, p.first, p.second);
      modoru.push_back(1);
      new_pick_order.push_back(P(x, y));
      b[x][y] = 0;
    }

    for (int j = 0; j < v.size(); j++) {
      if (b[v[j].first][v[j].second] == 0) continue;  // 既に拾われた箱はスキップ
      if (p.first != -1) {
        next_move(x, y, min_d, move_count, v[j].first, v[j].second, false);
        modoru.push_back(0);
        p.first = -1;
      }
      else {
        next_move(x, y, min_d, move_count, v[j].first, v[j].second);
        modoru.push_back(1);
      }
      
      new_pick_order.push_back(P(x, y));
      b[x][y] = 0;

      bool isLast = true;
      srep(k, j + 1, v.size())
      {
        if (b[v[k].first][v[k].second] != 0) {
          isLast = false;
          break;
        }
      }

      if (isLast) {
        // 出入り口へ戻る（途中で箱を拾いながら）
        while (true) {
          auto pickable_boxes = get_safe_pickable_boxes(x, y, min_d);

          if (!pickable_boxes.empty()) {
            // 安全に拾える箱がある場合、最も近い箱を拾う
            P closest_box = pickable_boxes[0];

            int after_mind = 0;
            int min_dist = closest_box.first + closest_box.second;
            for (const auto& box : pickable_boxes) {
              int dist = box.first + box.second;
              int tmp_min_d = min(min_d - (box.first + box.second) * current_state.grid[box.first][box.second].w, current_state.grid[box.first][box.second].d);
              //if (tmp_min_d > after_mind) {
              if (dist > min_dist) {
                closest_box = box;
                after_mind = tmp_min_d;
                min_dist = dist;
              }
            }

            next_move(x, y, min_d, move_count, closest_box.first, closest_box.second);
            modoru.push_back(1);
            new_pick_order.push_back(P(x, y));
            b[x][y] = 0;
          }
          else {
            break;
          }
        }
      }
    }

    next_move(x, y, min_d, move_count, 0, 0);
    min_d = INT_MAX;
  }

  return new_pick_order;
}


struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams)
{
  store_best_state();

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;

  current_state.init(input);
  ll old_score = calc_score_with_pick_order();
  best_pick_order_score = old_score;
  best_pick_order = pick_order;
  vector<P> keep_order = pick_order;  // 巻き戻し用の保持
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) { break; }
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[3];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % pick_order.size();
      ra2 = rand_xorshift() % pick_order.size();
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % pick_order.size();
      }
      swap(pick_order[ra1], pick_order[ra2]);
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
      auto vv = calc_score_with_pick_order2();
      // vvをシャッフル
      std::shuffle(vv.begin(), vv.end(), std::mt19937(std::random_device()()));
      pick_order.clear();
      for (const auto& v : vv) {
        for (const auto& p : v) {
          pick_order.push_back(p);
        }
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
      ra1 = rand_xorshift() % pick_order.size();
      ra2 = rand_xorshift() % pick_order.size();
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % pick_order.size();
      }
      if (ra1 > ra2) {
        swap(ra1, ra2);
      }
      for (int i = ra1; i < ra2; ++i) {
        swap(pick_order[i], pick_order[i + 1]);
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[3]) {
      keep_order = pick_order;  // 巻き戻し用の保持
      pick_order = after_method();
    }

    // スコア計算
    ll tmp_score = calc_score_with_pick_order();

    // 焼きなましで採用判定
    double diff_score = (tmp_score - old_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      old_score = tmp_score;

      // ベスト更新
      if (tmp_score > best_pick_order_score) {
        //cerr << "Best score updated: " << tmp_score << endl;
        best_pick_order_score = tmp_score;
        best_pick_order = pick_order;
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        swap(pick_order[ra1], pick_order[ra2]);
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
        // 近傍操作3 の巻き戻し
        for (int i = ra2; i > ra1; --i) {
          swap(pick_order[i], pick_order[i - 1]);
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[3]) {
        pick_order = keep_order;  // 巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << "Final score: " << best_pick_order_score << ", loop count: " << loop << endl;
  }

  pick_order = best_pick_order;
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  input_data(case_num);

  initialize_state();

  ofstream ofs;
  open_ofs(case_num, ofs);

  pick_order.clear();
  build_initial_solution();

  // 焼きなまし実行
  run_simulated_annealing(annealingParams);

  modoru.clear();
  pick_order = after_method();

  remake_state(true);
  // 最終的な状態をベスト状態に更新
  store_best_state();

  // 解答を出力
  output_data(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main()
{
  exec_mode = 2;  // 0: 標準入出力, 1-2: ファイル入出力

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 100048.0;
  annealingParams.start_temperature[1] = 2048.0;
  annealingParams.start_temperature[2] = 2048.0;
  annealingParams.start_temperature[3] = 2048.0;
  annealingParams.start_temperature[4] = 2048.0;
  annealingParams.start_temperature[5] = 2048.0;
  annealingParams.start_temperature[6] = 2048.0;
  annealingParams.start_temperature[7] = 2048.0;
  annealingParams.start_temperature[8] = 2048.0;
  annealingParams.start_temperature[9] = 2048.0;
  annealingParams.end_temperature = 0.0;
  annealingParams.score_scale = 12345.0;
  annealingParams.operation_thresholds[0] = 10000;
  annealingParams.operation_thresholds[1] = 10110;
  annealingParams.operation_thresholds[2] = 11000;
  annealingParams.operation_thresholds[3] = 11100;
  annealingParams.operation_thresholds[4] = 500;
  annealingParams.operation_thresholds[5] = 600;
  annealingParams.operation_thresholds[6] = 700;
  annealingParams.operation_thresholds[7] = 800;
  annealingParams.operation_thresholds[8] = 900;
  annealingParams.operation_thresholds[9] = 1000;

  if (exec_mode == 0) {
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 150; ++i) {
      ll score = solve_case(i, annealingParams);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(3) << i << ", "
          << "score = " << setw(6) << score << ", "
          << "sum = " << setw(8) << sum_score << ", "
          << "time = " << setw(5) << get_elapsed_time() << ", "
          << endl;
      }
    }
  }
  else if (exec_mode == 3) {
    int loop_count = 0;
    AnnealingParams best_annealingParams;
    ll best_sum_score = 0;

    while (true) {
      AnnealingParams new_annealingParams;
      new_annealingParams.start_temperature[0] = pow(2.0, rand_01() * 20);
      new_annealingParams.end_temperature = 0.0;
      new_annealingParams.score_scale = pow(2.0, rand_01() * 20);
      new_annealingParams.operation_thresholds[0] = rand() % 101;

      ll sum_score = 0;
      for (int i = 0; i < 150; ++i) {
        ll score = solve_case(i, new_annealingParams);
        sum_score += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cerr << "loop_count = " << loop_count
        << ", sum_score = " << sum_score
        << ", start_temperature = " << new_annealingParams.start_temperature[0]
        << ", end_temperature = " << new_annealingParams.end_temperature
        << ", score_scale = " << new_annealingParams.score_scale
        << ", operation_thresholds = " << new_annealingParams.operation_thresholds[0]
        << endl;

      if (sum_score > best_sum_score) {
        best_sum_score = sum_score;
        best_annealingParams = new_annealingParams;
      }

      loop_count++;
    }
  }

  return 0;
}
