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

// �^�C�}�[
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

namespace /* �������C�u���� */
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
}  // namespace

// 定数
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

// パス情報を表す構造体
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

// グローバル変数を名前空間に整理
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

// �����P�[�X�񂷂Ƃ��ɓ����Ԃ�����l�ɖ߂�
void resetGlobalState()
{
  rep(i, COLOR_COUNT) { cells_by_color[i].clear(); }
  rep(i, TASK_COUNT) { task_paths[i].clear(); }
}

// ���͎󂯎��
void readInput(int problemNum)
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

  // �W�����͂���
  int _n, _m;
  if (!ifs.is_open()) {
    cin >> _n >> _m;
    cin >> start_i >> start_j;
    rep(i, BOARD_SIZE) cin >> board_chars[i];
    rep(i, TASK_COUNT) cin >> task_strings[i];
  }
  // �t�@�C�����͂���
  else {
    ifs >> _n >> _m;
    ifs >> start_i >> start_j;
    rep(i, BOARD_SIZE) ifs >> board_chars[i];
    rep(i, TASK_COUNT) ifs >> task_strings[i];
  }

  rep(i, BOARD_SIZE)
  {
    rep(j, BOARD_SIZE)
    {
      board_color[i][j] = board_chars[i][j] - 'A';
    }
  }
  rep(i, TASK_COUNT)
  {
    rep(j, PATH_LENGTH)
    {
      task_colors[i][j] = task_strings[i][j] - 'A';
    }
  }

  rep(i, BOARD_SIZE)
  {
    rep(j, BOARD_SIZE)
    {
      cells_by_color[board_color[i][j]].push_back(P(i, j));
    }
  }
  rep(i, COLOR_COUNT)
  {
    cell_count_by_color[i] = cells_by_color[i].size();
  }
}

// �o�̓t�@�C���X�g���[���I�[�v��
void openOutputFile(int probNum, ofstream& ofs)
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

// �X�R�A�v�Z
ll calculateTotalScore()
{
  int score = MAX_SCORE - TASK_COUNT * TASK_PENALTY;
  rep(i, TASK_COUNT)
  {
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

// �����𐶐�
void buildInitialSolution()
{
  std::random_device seed_gen;
  std::mt19937   rng(seed_gen());

  best_score = -1;

  /* ---------- 1 st  pass (�S�^�X�N�Ώ�) ---------- */
  rep(attempt_idx, GREEDY_ITERATIONS)
  {
    int   used[TASK_COUNT]{};          // ���łɑI�΂ꂽ�^�X�N
    int   cur_i = start_i, cur_j = start_j;

    std::vector<int> shuffled_tasks;
    rep(t, TASK_COUNT) shuffled_tasks.push_back(t);

    rep(pos, TASK_COUNT)
    {
      int best_delta = INF;
      int best_task = -1;
      int best_path_id = -1;

      std::shuffle(shuffled_tasks.begin(), shuffled_tasks.end(), rng);

      rep(task_idx, TASK_COUNT)
      {
        int task_id = shuffled_tasks[task_idx];
        if (used[task_id]) continue;

        rep(path_id, BOARD_SIZE)
        {
          if (task_paths[task_id].size() <= path_id) break;

          int dist_start = manhattanDistance(
            task_paths[task_id][path_id].start_i,
            task_paths[task_id][path_id].start_j,
            cur_i, cur_j);

          int bonus3 = 0, bonus2 = 0;
          if (pos > 0) {
            bonus3 = score_when_share_last3(
              task_paths[task_order[pos - 1]][path_index[pos - 1]],
              task_paths[task_id][path_id]);
            bonus2 = score_when_share_last2(
              task_paths[task_order[pos - 1]][path_index[pos - 1]],
              task_paths[task_id][path_id]);
          }

          int delta;
          if (bonus3 > 0) {
            delta = dist_start - bonus3 + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
          }
          else if (bonus2 > 0) {
            delta = dist_start - bonus2 + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
          }
          else {
            int same_cell = (pos > 0 && dist_start == 0) ? 1 : 0;
            delta = dist_start - same_cell + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
          }

          if (delta < best_delta || (delta == best_delta && Random::xorshift() % 2)) {
            best_delta = delta;
            best_task = task_id;
            best_path_id = path_id;
          }
        }
      }

      task_order[pos] = best_task;
      path_index[pos] = best_path_id;
      used[best_task] = 1;
      cur_i = task_paths[best_task][best_path_id].goal_i;
      cur_j = task_paths[best_task][best_path_id].goal_j;
    }

    int score = calculateTotalScore();
    if (score > best_score) {
      best_score = score;
      rep(i, TASK_COUNT)
      {
        best_task_order[i] = task_order[i];
        best_path_index[i] = path_index[i];
      }
    }
  }

  /* ---- �ȍ~�Aprefix ��Œ肵�Ȃ��� 3 �i�K�ŉ��� ---- */
  auto restore_best = [&]() {
    rep(i, TASK_COUNT)
    {
      task_order[i] = best_task_order[i];
      path_index[i] = best_path_index[i];
    }
    };
  restore_best();

  const struct
  {
    int iterations;
    int prefix_len;
  } phases[] = {
      { GREEDY_ITERATIONS , 150 },
      {1000 , 180 },
      {1000 , 190 }
  };

  for (auto [iter_cnt, prefix_len] : phases) {
    rep(attempt_idx, iter_cnt)
    {
      int   used[TASK_COUNT]{};
      int   cur_i = start_i, cur_j = start_j;

      std::vector<int> shuffled_tasks;
      rep(t, TASK_COUNT) shuffled_tasks.push_back(t);

      /* prefix �͑O�̃x�X�g����̂܂܌Œ� */
      rep(pos, prefix_len) used[task_order[pos]] = 1;

      srep(pos, prefix_len, TASK_COUNT)
      {
        int best_delta = INF;
        int best_task = -1;
        int best_path_id = -1;

        std::shuffle(shuffled_tasks.begin(), shuffled_tasks.end(), rng);

        rep(task_idx, TASK_COUNT)
        {
          int task_id = shuffled_tasks[task_idx];
          if (used[task_id]) continue;

          rep(path_id, BOARD_SIZE)
          {
            if (task_paths[task_id].size() <= path_id) break;

            int dist_start = manhattanDistance(
              task_paths[task_id][path_id].start_i,
              task_paths[task_id][path_id].start_j,
              cur_i, cur_j);

            int bonus3 = 0, bonus2 = 0;
            if (pos > 0) {
              bonus3 = score_when_share_last3(
                task_paths[task_order[pos - 1]][path_index[pos - 1]],
                task_paths[task_id][path_id]);
              bonus2 = score_when_share_last2(
                task_paths[task_order[pos - 1]][path_index[pos - 1]],
                task_paths[task_id][path_id]);
            }

            int delta;
            if (bonus3 > 0) {
              delta = dist_start - bonus3 + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
            }
            else if (bonus2 > 0) {
              delta = dist_start - bonus2 + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
            }
            else {
              int same_cell = (pos > 0 && dist_start == 0) ? 1 : 0;
              delta = dist_start - same_cell + task_paths[task_id][path_id].path_cost - min_path_cost[task_id];
            }

            if (delta < best_delta || (delta == best_delta && Random::xorshift() % 2)) {
              best_delta = delta;
              best_task = task_id;
              best_path_id = path_id;
            }
          }
        }

        task_order[pos] = best_task;
        path_index[pos] = best_path_id;
        used[best_task] = 1;
        cur_i = task_paths[best_task][best_path_id].goal_i;
        cur_j = task_paths[best_task][best_path_id].goal_j;
      }

      int score = calculateTotalScore();
      if (score > best_score) {
        best_score = score;
        rep(i, TASK_COUNT)
        {
          best_task_order[i] = task_order[i];
          best_path_index[i] = path_index[i];
        }
      }
    }
    restore_best();
  }
}

// �𓚏o��
void writeAnswer(ofstream& ofs)
{
  final_cells.clear();
  rep(i, TASK_COUNT)
  {
    Path path = task_paths[task_order[i]][path_index[i]];
    rep(j, PATH_LENGTH)
    {
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

// �C���f�b�N�X 2 �_�����ւ���ߖT����
void twoSwap(double temperature)
{
  /* --- �X���b�v�ΏۃC���f�b�N�X�𖳍�ׂɌ��� --- */
  int idx_a = Random::xorshift() % TASK_COUNT;
  int idx_b = Random::xorshift() % TASK_COUNT;
  if (idx_a > idx_b) std::swap(idx_a, idx_b);
  if (idx_b - idx_a <= 1) return;                 // �A���v�f�͖���

  int prev_a = idx_a - 1;
  int next_a = idx_a + 1;
  int prev_b = idx_b - 1;
  int next_b = idx_b + 1;

  int cost_before = 0;
  int cost_after = 0;

  /* ---------- �P�[�X�@: �ǂ����[�v�f�ł͂Ȃ� ---------- */
  if (idx_a != 0 && idx_b != TASK_COUNT - 1) {
    { // ����O�R�X�g
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_i,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_j);
      int d2 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].goal_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int d4 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].goal_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].goal_j,
        task_paths[task_order[next_b]][path_index[next_b]].start_i,
        task_paths[task_order[next_b]][path_index[next_b]].start_j);
      int same_cnt = (d1 == 0) + (d2 == 0) + (d3 == 0) + (d4 == 0);
      cost_before = d1 + d2 + d3 + d4 - same_cnt;
    }
    { // �����R�X�g
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int d2 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].goal_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].goal_j,
        task_paths[task_order[next_b]][path_index[next_b]].start_i,
        task_paths[task_order[next_b]][path_index[next_b]].start_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_i,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_j);
      int d4 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].goal_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int same_cnt = (d1 == 0) + (d2 == 0) + (d3 == 0) + (d4 == 0);
      cost_after = d1 + d2 + d3 + d4 - same_cnt;
    }
  }
  /* ---------- �P�[�X�A: idx_a ���擪 ---------- */
  else if (idx_a == 0 && idx_b != TASK_COUNT - 1) {
    { // ����O
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        start_i, start_j);
      int d2 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].goal_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int d4 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].goal_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].goal_j,
        task_paths[task_order[next_b]][path_index[next_b]].start_i,
        task_paths[task_order[next_b]][path_index[next_b]].start_j);
      int same_cnt = (d2 == 0) + (d3 == 0) + (d4 == 0);
      cost_before = d1 + d2 + d3 + d4 - same_cnt;
    }
    { // �����
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int d2 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].goal_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].goal_j,
        task_paths[task_order[next_b]][path_index[next_b]].start_i,
        task_paths[task_order[next_b]][path_index[next_b]].start_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        start_i, start_j);
      int d4 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].goal_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int same_cnt = (d1 == 0) + (d2 == 0) + (d4 == 0);
      cost_after = d1 + d2 + d3 + d4 - same_cnt;
    }
  }
  /* ---------- �P�[�X�B: idx_b ������ ---------- */
  else if (idx_a != 0 && idx_b == TASK_COUNT - 1) {
    { // ����O
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_i,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_j);
      int d2 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].goal_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int same_cnt = (d1 == 0) + (d2 == 0) + (d3 == 0);
      cost_before = d1 + d2 + d3 - same_cnt;
    }
    { // �����
      int d1 = manhattanDistance(task_paths[task_order[idx_a]][path_index[idx_a]].start_i,
        task_paths[task_order[idx_a]][path_index[idx_a]].start_j,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_i,
        task_paths[task_order[prev_b]][path_index[prev_b]].goal_j);
      int d3 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].start_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].start_j,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_i,
        task_paths[task_order[prev_a]][path_index[prev_a]].goal_j);
      int d4 = manhattanDistance(task_paths[task_order[idx_b]][path_index[idx_b]].goal_i,
        task_paths[task_order[idx_b]][path_index[idx_b]].goal_j,
        task_paths[task_order[next_a]][path_index[next_a]].start_i,
        task_paths[task_order[next_a]][path_index[next_a]].start_j);
      int same_cnt = (d1 == 0) + (d3 == 0) + (d4 == 0);
      cost_after = d1 + d3 + d4 - same_cnt;
    }
  }
  else {
    return;
  }

  int    delta_cost = cost_before - cost_after;
  double accept_prob = std::exp(static_cast<double>(delta_cost) / temperature);

  // �Ă��Ȃ܂��F���P or ����P�ł�m���Ŏ�
  // if (Random::rand01() < accept_prob)
  if (delta_cost >= 0) {
    std::swap(task_order[idx_a], task_order[idx_b]);
    std::swap(path_index[idx_a], path_index[idx_b]);
  }
}

// ID�ύX
void changePathId(double temperature)
{
  int x = Random::xorshift() % TASK_COUNT;
  int y = Random::xorshift() % 10;
  if (task_paths[task_order[x]].size() <= y) return;
  if (x == 0 || x == TASK_COUNT - 1) return;
  if (y == path_index[x]) return;
  int x11 = x - 1;
  int x12 = x + 1;

  int beforeScore = 0;
  int afterScore = 0;
  {
    int dist1 = manhattanDistance(task_paths[task_order[x]][path_index[x]].start_i, task_paths[task_order[x]][path_index[x]].start_j,
      task_paths[task_order[x11]][path_index[x11]].goal_i, task_paths[task_order[x11]][path_index[x11]].goal_j);
    int dist2 = manhattanDistance(task_paths[task_order[x]][path_index[x]].goal_i, task_paths[task_order[x]][path_index[x]].goal_j,
      task_paths[task_order[x12]][path_index[x12]].start_i, task_paths[task_order[x12]][path_index[x12]].start_j);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    beforeScore = dist1 + dist2 - same + task_paths[task_order[x]][path_index[x]].path_cost;
  }
  {
    int dist1 = manhattanDistance(task_paths[task_order[x]][y].start_i, task_paths[task_order[x]][y].start_j,
      task_paths[task_order[x11]][path_index[x11]].goal_i, task_paths[task_order[x11]][path_index[x11]].goal_j);
    int dist2 = manhattanDistance(task_paths[task_order[x]][y].goal_i, task_paths[task_order[x]][y].goal_j,
      task_paths[task_order[x12]][path_index[x12]].start_i, task_paths[task_order[x12]][path_index[x12]].start_j);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    afterScore = dist1 + dist2 - same + task_paths[task_order[x]][y].path_cost;
  }

  int diffScore = beforeScore - afterScore;
  double prob = exp((double)diffScore / temperature);
  //if (prob > Random::rand01()) {
  if (diffScore >= 0) {
    path_index[x] = y;
  }
}

void simulatedAnnealing()
{
  int loop = 0;
  double startTemperature = START_TEMPERATURE;
  double endTemperature = END_TEMPERATURE;
  double nowProgress = 0;
  while (false) {
    loop++;
    if (loop % 100 == 0) {
      double nowTime = Timer::getElapsedTime();
      nowProgress = nowTime / TIME_LIMIT_SEC;
      if (nowProgress > 1.0) break;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    int ra = Random::xorshift() % 100;
    if (ra < SWAP_RATIO) {
      // 2�_�X���b�v
      twoSwap(temperature);
    }
    else {
      changePathId(temperature);
    }
  }

  if (mode != 0) {
    cout << "loop = " << loop << endl;
  }
}

ll solveSingleCase(int probNum)
{
  Timer::start();

  // �����P�[�X�񂷂Ƃ��ɓ����Ԃ�����l�ɖ߂�
  resetGlobalState();

  // ���͎󂯎��
  readInput(probNum);

  // �o�̓t�@�C���X�g���[���I�[�v��
  ofstream ofs;
  openOutputFile(probNum, ofs);

  // dp
  {
    int dp[PATH_LENGTH][DP_SIZE][DP_SIZE];
    int dp2[PATH_LENGTH][DP_SIZE][DP_SIZE];
    rep(i, TASK_COUNT)
    {
      rep(j, PATH_LENGTH)
      {
        int x = task_colors[i][0];
        if (j == 0) {
          rep(k, cell_count_by_color[x])
          {
            rep(l, cell_count_by_color[x]) { dp[j][k][l] = INF; }
          }
        }
        else {
          int y = task_colors[i][j];
          rep(k, cell_count_by_color[x])
          {
            rep(l, cell_count_by_color[y]) { dp[j][k][l] = INF; }
          }
        }
      }

      int x = task_colors[i][0];
      rep(j, PATH_LENGTH)
      {
        if (j == 0) {
          rep(k, cell_count_by_color[x])
          {
            rep(l, cell_count_by_color[x])
            {
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
          rep(k, cell_count_by_color[x])
          {
            rep(l, cell_count_by_color[y])
            {
              rep(o, cell_count_by_color[z])
              {
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
      rep(j, cell_count_by_color[x])
      {
        rep(k, cell_count_by_color[z])
        {
          Path path;
          path.start_i = cells_by_color[x][j].first;
          path.start_j = cells_by_color[x][j].second;
          path.goal_i = cells_by_color[z][k].first;
          path.goal_j = cells_by_color[z][k].second;
          path.path_cost = dp[4][j][k];
          vector<P> tmp;
          int now = k;
          drep(l, PATH_LENGTH)
          {
            int y = task_colors[i][l];
            tmp.push_back(cells_by_color[y][now]);
            if (l == 0) break;
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

  // �����𐶐�
  buildInitialSolution();

  simulatedAnnealing();

  // �𓚂�o��
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
    srep(i, 0, 10)
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
