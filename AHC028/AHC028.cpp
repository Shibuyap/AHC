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

namespace Timer
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double getElapsedTime()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

namespace
{
  namespace Random
  {
    uint32_t xorshift()
    {
      static uint32_t x = 123456789;
      static uint32_t y = 362436069;
      static uint32_t z = 521288629;
      static uint32_t w = 88675123;
      uint32_t t = x ^ (x << 11);
      x = y;
      y = z;
      z = w;
      return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    double rand01()
    {
      return (xorshift() + 0.5) * (1.0 / UINT_MAX);
    }
  }
}

namespace Constants
{
  constexpr int COLOR_COUNT = 26;
  constexpr int TASK_COUNT = 200;
  constexpr int PATH_LENGTH = 5;
  constexpr int BOARD_SIZE = 15;
  constexpr int INF = 1001001001;
  constexpr double TIME_LIMIT_SEC = 1.8;
  constexpr int MAX_SCORE = 10000;
  constexpr int TASK_PENALTY = 5;
  constexpr int DP_SIZE = 40;
  constexpr int GREEDY_ITERATIONS = 400;
  constexpr int SWAP_RATIO = 50;
  constexpr double START_TEMPERATURE = 2.0;
  constexpr double END_TEMPERATURE = 0.0;
}

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

using namespace Constants;
int mode;

struct Path
{
  int start_i;
  int start_j;
  int goal_i;
  int goal_j;
  int path_cost;
  vector<P> cells;

  bool operator<(const Path& right) const { return path_cost < right.path_cost; }
};

namespace GlobalData
{
  // 入力データ
  int start_i, start_j;
  string board_chars[BOARD_SIZE];
  int board_color[BOARD_SIZE][BOARD_SIZE];
  string task_strings[TASK_COUNT];
  int task_colors[TASK_COUNT][PATH_LENGTH];
  vector<P> cells_by_color[COLOR_COUNT];
  int cell_count_by_color[COLOR_COUNT];

  // 解の情報
  int task_order[TASK_COUNT];
  int path_index[TASK_COUNT];
  vector<P> final_cells;
  int best_task_order[TASK_COUNT];
  int best_path_index[TASK_COUNT];
  int best_score;

  // パス情報
  vector<Path> task_paths[TASK_COUNT];
  int min_path_cost[TASK_COUNT];
}

using namespace GlobalData;

inline int manhattanDistance(int i1, int j1, int i2, int j2)
{
  return abs(i1 - i2) + abs(j1 - j2);
}

void resetGlobalState()
{
  for (int i = 0; i < COLOR_COUNT; ++i) { cells_by_color[i].clear(); }
  for (int i = 0; i < TASK_COUNT; ++i) { task_paths[i].clear(); }
}

void readInput(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  int _n, _m;
  if (!ifs.is_open()) {
    cin >> _n >> _m;
    cin >> start_i >> start_j;
    for (int i = 0; i < BOARD_SIZE; ++i) cin >> board_chars[i];
    for (int i = 0; i < TASK_COUNT; ++i) cin >> task_strings[i];
  }
  else {
    ifs >> _n >> _m;
    ifs >> start_i >> start_j;
    for (int i = 0; i < BOARD_SIZE; ++i) ifs >> board_chars[i];
    for (int i = 0; i < TASK_COUNT; ++i) ifs >> task_strings[i];
  }

  for (int i = 0; i < BOARD_SIZE; ++i) {
    for (int j = 0; j < BOARD_SIZE; ++j) {
      board_color[i][j] = board_chars[i][j] - 'A';
    }
  }
  for (int i = 0; i < TASK_COUNT; ++i) {
    for (int j = 0; j < PATH_LENGTH; ++j) {
      task_colors[i][j] = task_strings[i][j] - 'A';
    }
  }

  for (int i = 0; i < BOARD_SIZE; ++i) {
    for (int j = 0; j < BOARD_SIZE; ++j) {
      cells_by_color[board_color[i][j]].push_back(P(i, j));
    }
  }
  for (int i = 0; i < COLOR_COUNT; ++i) {
    cell_count_by_color[i] = cells_by_color[i].size();
  }
}

void openOutputFile(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

int score_when_share_last3(const Path& x, const Path& y)
{
  if (x.cells[2] == y.cells[0] && x.cells[3] == y.cells[1] && x.cells[4] == y.cells[2]) {
    int res = manhattanDistance(x.cells[4].first, x.cells[4].second, y.cells[0].first, y.cells[0].second)
      + manhattanDistance(y.cells[0].first, y.cells[0].second, y.cells[1].first, y.cells[1].second)
      + manhattanDistance(y.cells[1].first, y.cells[1].second, y.cells[2].first, y.cells[2].second) + 3;
    return res;
  }
  return 0;
}

int score_when_share_last2(const Path& x, const Path& y)
{
  if (x.cells[3] == y.cells[0] && x.cells[4] == y.cells[1]) {
    int res = manhattanDistance(x.cells[4].first, x.cells[4].second, y.cells[0].first, y.cells[0].second)
      + manhattanDistance(y.cells[0].first, y.cells[0].second, y.cells[1].first, y.cells[1].second) + 2;
    return res;
  }
  return 0;
}

ll calculateTotalScore()
{
  int score = MAX_SCORE - TASK_COUNT * TASK_PENALTY;
  for (int i = 0; i < TASK_COUNT; ++i) {
    if (i == 0) {
      score -= manhattanDistance(task_paths[task_order[i]][path_index[i]].start_i, task_paths[task_order[i]][path_index[i]].start_j, start_i, start_j);
    }
    else {
      int dist =
        manhattanDistance(task_paths[task_order[i]][path_index[i]].start_i, task_paths[task_order[i]][path_index[i]].start_j,
          task_paths[task_order[i - 1]][path_index[i - 1]].goal_i, task_paths[task_order[i - 1]][path_index[i - 1]].goal_j);
      score -= dist;
      int same3 = score_when_share_last3(task_paths[task_order[i - 1]][path_index[i - 1]], task_paths[task_order[i]][path_index[i]]);
      int same2 = score_when_share_last2(task_paths[task_order[i - 1]][path_index[i - 1]], task_paths[task_order[i]][path_index[i]]);
      if (same3 > 0) {
        score += same3;
      }
      else if (same2 > 0) {
        score += same2;
      }
      else {
        if (dist == 0) score++;
      }
    }

    score -= task_paths[task_order[i]][path_index[i]].path_cost;
  }
  return score;
}

// タスクとパスの最適な組み合わせを選択
struct TaskPathSelection
{
  int task_id;
  int path_id;
  int delta;
};

// デルタ値を計算
int calculateDelta(int task_id, int path_id, int pos, int cur_i, int cur_j)
{
  const Path& path = task_paths[task_id][path_id];
  int dist_start = manhattanDistance(path.start_i, path.start_j, cur_i, cur_j);

  int bonus3 = 0, bonus2 = 0;
  if (pos > 0) {
    bonus3 = score_when_share_last3(
      task_paths[task_order[pos - 1]][path_index[pos - 1]],
      path);
    bonus2 = score_when_share_last2(
      task_paths[task_order[pos - 1]][path_index[pos - 1]],
      path);
  }

  int delta;
  if (bonus3 > 0) {
    delta = dist_start - bonus3 + path.path_cost - min_path_cost[task_id];
  }
  else if (bonus2 > 0) {
    delta = dist_start - bonus2 + path.path_cost - min_path_cost[task_id];
  }
  else {
    int same_cell = (pos > 0 && dist_start == 0) ? 1 : 0;
    delta = dist_start - same_cell + path.path_cost - min_path_cost[task_id];
  }

  return delta;
}

// 最適なタスクとパスを選択
TaskPathSelection selectBestTaskPath(const vector<int>& shuffled_tasks, const int* used,
  int pos, int cur_i, int cur_j)
{
  TaskPathSelection best = { -1, -1, INF };

  for (int task_id : shuffled_tasks) {
    if (used[task_id]) { continue; }

    for (int path_id = 0; path_id < min((int)task_paths[task_id].size(), BOARD_SIZE); ++path_id) {
      int delta = calculateDelta(task_id, path_id, pos, cur_i, cur_j);

      if (delta < best.delta || (delta == best.delta && Random::xorshift() % 2)) {
        best.task_id = task_id;
        best.path_id = path_id;
        best.delta = delta;
      }
    }
  }

  return best;
}

// 貪欲法による解の構築
void greedyConstruction(std::mt19937& rng, int start_pos, int end_pos, const int* initial_used = nullptr)
{
  int used[TASK_COUNT]{};
  if (initial_used) {
    for (int i = 0; i < TASK_COUNT; ++i) {
      used[i] = initial_used[i];
    }
  }

  int cur_i = start_i, cur_j = start_j;
  if (start_pos > 0) {
    const Path& prev_path = task_paths[task_order[start_pos - 1]][path_index[start_pos - 1]];
    cur_i = prev_path.goal_i;
    cur_j = prev_path.goal_j;
  }

  vector<int> shuffled_tasks;
  for (int t = 0; t < TASK_COUNT; ++t) {
    shuffled_tasks.push_back(t);
  }

  for (int pos = start_pos; pos < end_pos; ++pos) {
    std::shuffle(shuffled_tasks.begin(), shuffled_tasks.end(), rng);

    TaskPathSelection selection = selectBestTaskPath(shuffled_tasks, used, pos, cur_i, cur_j);

    task_order[pos] = selection.task_id;
    path_index[pos] = selection.path_id;
    used[selection.task_id] = 1;

    const Path& selected_path = task_paths[selection.task_id][selection.path_id];
    cur_i = selected_path.goal_i;
    cur_j = selected_path.goal_j;
  }
}

// ベスト解を保存
void saveBestSolution()
{
  for (int i = 0; i < TASK_COUNT; ++i) {
    best_task_order[i] = task_order[i];
    best_path_index[i] = path_index[i];
  }
}

// ベスト解を復元
void restoreBestSolution()
{
  for (int i = 0; i < TASK_COUNT; ++i) {
    task_order[i] = best_task_order[i];
    path_index[i] = best_path_index[i];
  }
}

void buildInitialSolution()
{
  std::random_device seed_gen;
  std::mt19937 rng(seed_gen());

  best_score = -1;

  // 初期解の構築（完全ランダム）
  for (int attempt = 0; attempt < GREEDY_ITERATIONS; ++attempt) {
    greedyConstruction(rng, 0, TASK_COUNT);

    int score = calculateTotalScore();
    if (score > best_score) {
      best_score = score;
      saveBestSolution();
    }
  }

  restoreBestSolution();

  // フェーズベースの改善
  const struct
  {
    int iterations;
    int prefix_len;
  } phases[] = {
    {GREEDY_ITERATIONS, 150},
    {1000, 180},
    {1000, 190}
  };

  for (auto [iter_cnt, prefix_len] : phases) {
    for (int attempt = 0; attempt < iter_cnt; ++attempt) {
      // 既存の解の一部を保持して再構築
      int used[TASK_COUNT]{};
      for (int pos = 0; pos < prefix_len; ++pos) {
        used[task_order[pos]] = 1;
      }

      greedyConstruction(rng, prefix_len, TASK_COUNT, used);

      int score = calculateTotalScore();
      if (score > best_score) {
        best_score = score;
        saveBestSolution();
      }
    }
    restoreBestSolution();
  }
}

void writeAnswer(ofstream& ofs)
{
  final_cells.clear();
  for (int i = 0; i < TASK_COUNT; ++i) {
    Path path = task_paths[task_order[i]][path_index[i]];
    for (int j = 0; j < PATH_LENGTH; ++j) {
      if (j <= 2 && i > 0 && score_when_share_last3(task_paths[task_order[i - 1]][path_index[i - 1]], path) > 0) {
        continue;
      }
      else if (j <= 1 && i > 0 && score_when_share_last2(task_paths[task_order[i - 1]][path_index[i - 1]], path) > 0) {
        continue;
      }
      else if (j == 0 && !final_cells.empty()) {
        if (final_cells.back() == path.cells[0]) {
          continue;
        }
      }
      final_cells.push_back(path.cells[j]);
    }
  }
  if (mode == 0) {
    for (auto p : final_cells) {
      cout << p.first << ' ' << p.second << endl;
    }
  }
  else {
    for (auto p : final_cells) {
      ofs << p.first << ' ' << p.second << endl;
    }
  }
}

// ヘルパー関数：タスクのパスを取得
inline const Path& getTaskPath(int idx)
{
  return task_paths[task_order[idx]][path_index[idx]];
}

// ヘルパー関数：2点間の距離を計算
inline int calcDistance(int idx1, int idx2, bool useGoalForFirst, bool useStartForSecond)
{
  const Path& path1 = getTaskPath(idx1);
  const Path& path2 = getTaskPath(idx2);

  int i1 = useGoalForFirst ? path1.goal_i : path1.start_i;
  int j1 = useGoalForFirst ? path1.goal_j : path1.start_j;
  int i2 = useStartForSecond ? path2.start_i : path2.goal_i;
  int j2 = useStartForSecond ? path2.start_j : path2.goal_j;

  return manhattanDistance(i1, j1, i2, j2);
}

// ヘルパー関数：スタート地点からの距離を計算
inline int calcDistanceFromStart(int idx)
{
  const Path& path = getTaskPath(idx);
  return manhattanDistance(path.start_i, path.start_j, start_i, start_j);
}

// ヘルパー関数：コストを計算（距離の合計 - 同一地点ボーナス）
inline int calcCost(const vector<int>& distances)
{
  int total = 0;
  int sameCount = 0;
  for (int d : distances) {
    total += d;
    if (d == 0) sameCount++;
  }
  return total - sameCount;
}

void twoSwap(double temperature)
{
  // 2つのタスクインデックスをランダムに選択
  int idx_a = Random::xorshift() % TASK_COUNT;
  int idx_b = Random::xorshift() % TASK_COUNT;
  if (idx_a > idx_b) std::swap(idx_a, idx_b);

  // 隣接要素の入れ替えは無効
  if (idx_b - idx_a <= 1) return;

  int prev_a = idx_a - 1;
  int next_a = idx_a + 1;
  int prev_b = idx_b - 1;
  int next_b = idx_b + 1;

  int cost_before = 0;
  int cost_after = 0;

  // ケース1: 両方とも中間要素
  if (idx_a != 0 && idx_b != TASK_COUNT - 1) {
    // 現在のコスト
    vector<int> distances_before = {
      calcDistance(prev_a, idx_a, true, true),   // prev_a -> idx_a
      calcDistance(idx_a, next_a, true, true),   // idx_a -> next_a
      calcDistance(prev_b, idx_b, true, true),   // prev_b -> idx_b
      calcDistance(idx_b, next_b, true, true)    // idx_b -> next_b
    };
    cost_before = calcCost(distances_before);

    // 入れ替え後のコスト
    vector<int> distances_after = {
      calcDistance(prev_a, idx_b, true, true),   // prev_a -> idx_b
      calcDistance(idx_b, next_a, true, true),   // idx_b -> next_a
      calcDistance(prev_b, idx_a, true, true),   // prev_b -> idx_a
      calcDistance(idx_a, next_b, true, true)    // idx_a -> next_b
    };
    cost_after = calcCost(distances_after);
  }
  // ケース2: idx_aが先頭
  else if (idx_a == 0 && idx_b != TASK_COUNT - 1) {
    // 現在のコスト
    vector<int> distances_before = {
      calcDistanceFromStart(idx_a),              // start -> idx_a
      calcDistance(idx_a, next_a, true, true),   // idx_a -> next_a
      calcDistance(prev_b, idx_b, true, true),   // prev_b -> idx_b
      calcDistance(idx_b, next_b, true, true)    // idx_b -> next_b
    };
    cost_before = calcCost(distances_before);

    // 入れ替え後のコスト（d1がカウントされないバグを修正）
    int d1 = calcDistanceFromStart(idx_b);
    int d2 = calcDistance(idx_b, next_a, true, true);
    int d3 = calcDistance(prev_b, idx_a, true, true);
    int d4 = calcDistance(idx_a, next_b, true, true);
    int same_cnt = (d1 == 0) + (d2 == 0) + (d4 == 0);
    cost_after = d1 + d2 + d3 + d4 - same_cnt;
  }
  // ケース3: idx_bが末尾
  else if (idx_a != 0 && idx_b == TASK_COUNT - 1) {
    // 現在のコスト
    vector<int> distances_before = {
      calcDistance(prev_a, idx_a, true, true),   // prev_a -> idx_a
      calcDistance(idx_a, next_a, true, true),   // idx_a -> next_a
      calcDistance(prev_b, idx_b, true, true)    // prev_b -> idx_b
    };
    cost_before = calcCost(distances_before);

    // 入れ替え後のコスト
    vector<int> distances_after = {
      calcDistance(prev_a, idx_b, true, true),   // prev_a -> idx_b
      calcDistance(idx_b, next_a, true, true),   // idx_b -> next_a
      calcDistance(prev_b, idx_a, true, true)    // prev_b -> idx_a
    };
    cost_after = calcCost(distances_after);
  }
  else {
    return;
  }

  // 受理判定
  int delta_cost = cost_before - cost_after;
  if (delta_cost >= 0) {
    std::swap(task_order[idx_a], task_order[idx_b]);
    std::swap(path_index[idx_a], path_index[idx_b]);
  }
}

// パス変更時のコストを計算
int calculatePathChangeScore(int task_idx, int path_id)
{
  int prev_idx = task_idx - 1;
  int next_idx = task_idx + 1;

  const Path& prev_path = task_paths[task_order[prev_idx]][path_index[prev_idx]];
  const Path& curr_path = task_paths[task_order[task_idx]][path_id];
  const Path& next_path = task_paths[task_order[next_idx]][path_index[next_idx]];

  int dist1 = manhattanDistance(curr_path.start_i, curr_path.start_j, prev_path.goal_i, prev_path.goal_j);
  int dist2 = manhattanDistance(curr_path.goal_i, curr_path.goal_j, next_path.start_i, next_path.start_j);

  int same_count = (dist1 == 0) + (dist2 == 0);
  return dist1 + dist2 - same_count + curr_path.path_cost;
}

void changePathId(double temperature)
{
  // ランダムにタスクと新しいパスを選択
  int task_idx = Random::xorshift() % TASK_COUNT;
  int new_path_id = Random::xorshift() % 10;

  // 無効な操作をスキップ
  if (task_paths[task_order[task_idx]].size() <= new_path_id) return;
  if (task_idx == 0 || task_idx == TASK_COUNT - 1) return;
  if (new_path_id == path_index[task_idx]) return;

  // 現在のコストを計算
  int before_score = calculatePathChangeScore(task_idx, path_index[task_idx]);

  // 新しいパスでのコストを計算
  int after_score = calculatePathChangeScore(task_idx, new_path_id);

  // 差分計算と受理判定
  int diff_score = before_score - after_score;
  if (diff_score >= 0) {
    path_index[task_idx] = new_path_id;
  }
}


ll solveSingleCase(int probNum)
{
  Timer::start();

  resetGlobalState();

  readInput(probNum);

  ofstream ofs;
  openOutputFile(probNum, ofs);

  {
    int dp[PATH_LENGTH][DP_SIZE][DP_SIZE];
    int dp2[PATH_LENGTH][DP_SIZE][DP_SIZE];
    for (int i = 0; i < TASK_COUNT; ++i) {
      for (int j = 0; j < PATH_LENGTH; ++j) {
        int x = task_colors[i][0];
        if (j == 0) {
          for (int k = 0; k < cell_count_by_color[x]; ++k) {
            for (int l = 0; l < cell_count_by_color[x]; ++l) { dp[j][k][l] = INF; }
          }
        }
        else {
          int y = task_colors[i][j];
          for (int k = 0; k < cell_count_by_color[x]; ++k) {
            for (int l = 0; l < cell_count_by_color[y]; ++l) { dp[j][k][l] = INF; }
          }
        }
      }

      int x = task_colors[i][0];
      for (int j = 0; j < PATH_LENGTH; ++j) {
        if (j == 0) {
          for (int k = 0; k < cell_count_by_color[x]; ++k) {
            for (int l = 0; l < cell_count_by_color[x]; ++l) {
              if (k == l)
                dp[j][k][l] = 0;
              else
                dp[j][k][l] = INF;
            }
          }
        }
        else {
          int y = task_colors[i][j - 1];
          int z = task_colors[i][j];
          for (int k = 0; k < cell_count_by_color[x]; ++k) {
            for (int l = 0; l < cell_count_by_color[y]; ++l) {
              for (int o = 0; o < cell_count_by_color[z]; ++o) {
                int dist = abs(cells_by_color[z][o].first - cells_by_color[y][l].first) +
                  abs(cells_by_color[z][o].second - cells_by_color[y][l].second);
                if (dp[j][k][o] > dp[j - 1][k][l] + dist) {
                  dp[j][k][o] = dp[j - 1][k][l] + dist;
                  dp2[j][k][o] = l;
                }
              }
            }
          }
        }
      }

      int z = task_colors[i][4];
      for (int j = 0; j < cell_count_by_color[x]; ++j) {
        for (int k = 0; k < cell_count_by_color[z]; ++k) {
          Path path;
          path.start_i = cells_by_color[x][j].first;
          path.start_j = cells_by_color[x][j].second;
          path.goal_i = cells_by_color[z][k].first;
          path.goal_j = cells_by_color[z][k].second;
          path.path_cost = dp[4][j][k];
          vector<P> tmp;
          int now = k;
          for (int l = PATH_LENGTH - 1; l >= 0; --l)
          {
            int y = task_colors[i][l];
            tmp.push_back(cells_by_color[y][now]);
            if (l == 0) { break; }
            now = dp2[l][j][now];
          }
          reverse(tmp.begin(), tmp.end());
          path.cells = tmp;
          task_paths[i].push_back(path);
        }
      }

      sort(task_paths[i].begin(), task_paths[i].end());
      min_path_cost[i] = task_paths[i][0].path_cost;
    }
  }

  buildInitialSolution();

  writeAnswer(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calculateTotalScore();
    cerr << Timer::getElapsedTime() << " sec" << endl;
  }
  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    solveSingleCase(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i)
    {
      ll score = solveSingleCase(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
