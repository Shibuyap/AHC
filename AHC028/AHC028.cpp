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
  static uint32_t rand_uint32()
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
    return (rand_uint32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int COLOR_COUNT = 26;
const int MAX_TASK_COUNT = 210;
const int MAX_PATH_LENGTH = 5;
const int MAX_BOARD_SIZE = 20;

const int INF = 1001001001;
double time_limit_sec = 1.8;
int mode;

int board_size, task_count;
int start_i, start_j;
string board_chars[MAX_BOARD_SIZE];
int board_color[MAX_BOARD_SIZE][MAX_BOARD_SIZE];
string task_strings[MAX_TASK_COUNT];
int task_colors[MAX_TASK_COUNT][5];
vector<P> cells_by_color[COLOR_COUNT];
int cell_count_by_color[COLOR_COUNT];

int task_order[MAX_TASK_COUNT];
int path_index[MAX_TASK_COUNT];
vector<P> final_cells;
int best_task_order[MAX_TASK_COUNT];
int best_path_index[MAX_TASK_COUNT];
int best_score;

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
vector<Path> V[MAX_TASK_COUNT];
int baseVal[MAX_TASK_COUNT];

int manhattan_distance(int i1, int j1, int i2, int j2)
{
  return abs(i1 - i2) + abs(j1 - j2);
}

// 複数ケース回すときに内部状態を初期値に戻す
void reset_global_state()
{
  rep(i, COLOR_COUNT) { cells_by_color[i].clear(); }
  rep(i, MAX_TASK_COUNT) { V[i].clear(); }
}

// 入力受け取り
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
    cin >> board_size >> task_count;
    cin >> start_i >> start_j;
    rep(i, board_size) cin >> board_chars[i];
    rep(i, task_count) cin >> task_strings[i];
  }
  // ファイル入力する
  else {
    ifs >> board_size >> task_count;
    ifs >> start_i >> start_j;
    rep(i, board_size) ifs >> board_chars[i];
    rep(i, task_count) ifs >> task_strings[i];
  }

  rep(i, board_size)
  {
    rep(j, board_size) { board_color[i][j] = board_chars[i][j] - 'A'; }
  }
  rep(i, task_count)
  {
    rep(j, 5) { task_colors[i][j] = task_strings[i][j] - 'A'; }
  }

  rep(i, board_size)
  {
    rep(j, board_size) { cells_by_color[board_color[i][j]].push_back(P(i, j)); }
  }
  rep(i, COLOR_COUNT) { cell_count_by_color[i] = cells_by_color[i].size(); }
}

// 出力ファイルストリームオープン
void open_output_file(int probNum, ofstream& ofs)
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
    int res = manhattan_distance(x.cells[4].first, x.cells[4].second, y.cells[0].first, y.cells[0].second)
      + manhattan_distance(y.cells[0].first, y.cells[0].second, y.cells[1].first, y.cells[1].second)
      + manhattan_distance(y.cells[1].first, y.cells[1].second, y.cells[2].first, y.cells[2].second) + 3;
    return res;
  }
  return 0;
}

int score_when_share_last2(const Path& x, const Path& y)
{
  if (x.cells[3] == y.cells[0] && x.cells[4] == y.cells[1]) {
    int res = manhattan_distance(x.cells[4].first, x.cells[4].second, y.cells[0].first, y.cells[0].second)
      + manhattan_distance(y.cells[0].first, y.cells[0].second, y.cells[1].first, y.cells[1].second) + 2;
    return res;
  }
  return 0;
}

// スコア計算
ll calc_total_score()
{
  int score = 10000 - task_count * 5;
  rep(i, task_count)
  {
    if (i == 0) {
      score -= manhattan_distance(V[task_order[i]][path_index[i]].start_i, V[task_order[i]][path_index[i]].start_j, start_i, start_j);
    }
    else {
      int dist =
        manhattan_distance(V[task_order[i]][path_index[i]].start_i, V[task_order[i]][path_index[i]].start_j,
          V[task_order[i - 1]][path_index[i - 1]].goal_i, V[task_order[i - 1]][path_index[i - 1]].goal_j);
      score -= dist;
      int same3 = score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], V[task_order[i]][path_index[i]]);
      int same2 = score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], V[task_order[i]][path_index[i]]);
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

    score -= V[task_order[i]][path_index[i]].path_cost;
  }
  return score;
}

// 初期解生成
void build_initial_solution()
{
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());

  best_score = -1;

  rep(aespa, 400)
  {
    // 貪欲に作る
    int f[MAX_TASK_COUNT] = {};
    int ti = start_i;
    int tj = start_j;

    vector<int> nums;
    rep(j, task_count)nums.push_back(j);
    rep(i, task_count)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, task_count)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, MAX_BOARD_SIZE)
        {
          if (V[j].size() <= k)break;
          int dist1 = manhattan_distance(V[j][k].start_i, V[j][k].start_j, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
            same2 = score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].path_cost - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].path_cost - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].path_cost - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && rand_uint32() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      task_order[i] = miAns;
      path_index[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].goal_i;
      tj = V[miAns][miId].goal_j;
    }

    int score = calc_total_score();
    if (score > best_score) {
      best_score = score;
      rep(i, task_count)
      {
        best_task_order[i] = task_order[i];
        best_path_index[i] = path_index[i];
      }
    }
  }

  rep(i, task_count)
  {
    task_order[i] = best_task_order[i];
    path_index[i] = best_path_index[i];
  }


  rep(aespa, 400)
  {
    // 貪欲に作る
    int f[MAX_TASK_COUNT] = {};
    int ti = start_i;
    int tj = start_j;

    vector<int> nums;
    rep(j, task_count)nums.push_back(j);
    rep(i, 150)
    {
      f[task_order[i]] = 1;
    }
    srep(i, 150, task_count)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, task_count)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, MAX_BOARD_SIZE)
        {
          if (V[j].size() <= k)break;
          int dist1 = manhattan_distance(V[j][k].start_i, V[j][k].start_j, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
            same2 = score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].path_cost - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].path_cost - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].path_cost - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && rand_uint32() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      task_order[i] = miAns;
      path_index[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].goal_i;
      tj = V[miAns][miId].goal_j;
    }

    int score = calc_total_score();
    if (score > best_score) {
      best_score = score;
      rep(i, task_count)
      {
        best_task_order[i] = task_order[i];
        best_path_index[i] = path_index[i];
      }
    }
  }

  rep(i, task_count)
  {
    task_order[i] = best_task_order[i];
    path_index[i] = best_path_index[i];
  }

  rep(aespa, 1000)
  {
    // 貪欲に作る
    int f[MAX_TASK_COUNT] = {};
    int ti = start_i;
    int tj = start_j;

    vector<int> nums;
    rep(j, task_count)nums.push_back(j);
    rep(i, 180)
    {
      f[task_order[i]] = 1;
    }
    srep(i, 180, task_count)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, task_count)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, MAX_BOARD_SIZE)
        {
          if (V[j].size() <= k)break;
          int dist1 = manhattan_distance(V[j][k].start_i, V[j][k].start_j, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
            same2 = score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].path_cost - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].path_cost - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].path_cost - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && rand_uint32() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      task_order[i] = miAns;
      path_index[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].goal_i;
      tj = V[miAns][miId].goal_j;
    }

    int score = calc_total_score();
    if (score > best_score) {
      best_score = score;
      rep(i, task_count)
      {
        best_task_order[i] = task_order[i];
        best_path_index[i] = path_index[i];
      }
    }
  }

  rep(i, task_count)
  {
    task_order[i] = best_task_order[i];
    path_index[i] = best_path_index[i];
  }

  rep(aespa, 1000)
  {
    // 貪欲に作る
    int f[MAX_TASK_COUNT] = {};
    int ti = start_i;
    int tj = start_j;

    vector<int> nums;
    rep(j, task_count)nums.push_back(j);
    rep(i, 190)
    {
      f[task_order[i]] = 1;
    }
    srep(i, 190, task_count)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, task_count)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, MAX_BOARD_SIZE)
        {
          if (V[j].size() <= k)break;
          int dist1 = manhattan_distance(V[j][k].start_i, V[j][k].start_j, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
            same2 = score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].path_cost - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].path_cost - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].path_cost - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && rand_uint32() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      task_order[i] = miAns;
      path_index[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].goal_i;
      tj = V[miAns][miId].goal_j;
    }

    int score = calc_total_score();
    if (score > best_score) {
      best_score = score;
      rep(i, task_count)
      {
        best_task_order[i] = task_order[i];
        best_path_index[i] = path_index[i];
      }
    }
  }

  rep(i, task_count)
  {
    task_order[i] = best_task_order[i];
    path_index[i] = best_path_index[i];
  }
}

// 解答出力
void write_answer(ofstream& ofs)
{
  final_cells.clear();
  rep(i, task_count)
  {
    Path path = V[task_order[i]][path_index[i]];
    rep(j, 5)
    {
      if (j <= 2 && i > 0 && score_when_share_last3(V[task_order[i - 1]][path_index[i - 1]], path) > 0) {
        continue;
      }
      else if (j <= 1 && i > 0 && score_when_share_last2(V[task_order[i - 1]][path_index[i - 1]], path) > 0) {
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

void two_swap(double temperature)
{
  // 2点スワップ
  int x1 = rand_uint32() % task_count;
  int x2 = rand_uint32() % task_count;
  if (x1 > x2) swap(x1, x2);
  if (x2 - x1 <= 1) return;
  int x11 = x1 - 1;
  int x12 = x1 + 1;
  int x21 = x2 - 1;
  int x22 = x2 + 1;

  int beforeScore = 0;
  int afterScore = 0;
  if (x1 != 0 && x2 != task_count - 1) {
    {
      int dist1 = manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j,
        V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
      int dist2 = manhattan_distance(V[task_order[x1]][path_index[x1]].goal_i, V[task_order[x1]][path_index[x1]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);
      int dist3 = manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int dist4 = manhattan_distance(V[task_order[x2]][path_index[x2]].goal_i, V[task_order[x2]][path_index[x2]].goal_j,
        V[task_order[x22]][path_index[x22]].start_i, V[task_order[x22]][path_index[x22]].start_j);
      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 + dist4 - same;
    }
    {
      int dist1 = manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int dist2 = manhattan_distance(V[task_order[x1]][path_index[x1]].goal_i, V[task_order[x1]][path_index[x1]].goal_j,
        V[task_order[x22]][path_index[x22]].start_i, V[task_order[x22]][path_index[x22]].start_j);
      int dist3 = manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j,
        V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
      int dist4 = manhattan_distance(V[task_order[x2]][path_index[x2]].goal_i, V[task_order[x2]][path_index[x2]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist2 + dist3 + dist4 - same;
    }
  }
  else if (x1 == 0 && x2 != task_count - 1) {
    {
      int dist1 =
        manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j, start_i, start_j);
      int dist2 = manhattan_distance(V[task_order[x1]][path_index[x1]].goal_i, V[task_order[x1]][path_index[x1]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);
      int dist3 = manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int dist4 = manhattan_distance(V[task_order[x2]][path_index[x2]].goal_i, V[task_order[x2]][path_index[x2]].goal_j,
        V[task_order[x22]][path_index[x22]].start_i, V[task_order[x22]][path_index[x22]].start_j);
      int same = 0;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 + dist4 - same;
    }
    {
      int dist1 = manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int dist2 = manhattan_distance(V[task_order[x1]][path_index[x1]].goal_i, V[task_order[x1]][path_index[x1]].goal_j,
        V[task_order[x22]][path_index[x22]].start_i, V[task_order[x22]][path_index[x22]].start_j);
      int dist3 =
        manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j, start_i, start_j);
      int dist4 = manhattan_distance(V[task_order[x2]][path_index[x2]].goal_i, V[task_order[x2]][path_index[x2]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist2 + dist3 + dist4 - same;
    }
  }
  else if (x1 != 0 && x2 == task_count - 1) {
    {
      int dist1 = manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j,
        V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
      int dist2 = manhattan_distance(V[task_order[x1]][path_index[x1]].goal_i, V[task_order[x1]][path_index[x1]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);
      int dist3 = manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 - same;
    }
    {
      int dist1 = manhattan_distance(V[task_order[x1]][path_index[x1]].start_i, V[task_order[x1]][path_index[x1]].start_j,
        V[task_order[x21]][path_index[x21]].goal_i, V[task_order[x21]][path_index[x21]].goal_j);
      int dist3 = manhattan_distance(V[task_order[x2]][path_index[x2]].start_i, V[task_order[x2]][path_index[x2]].start_j,
        V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
      int dist4 = manhattan_distance(V[task_order[x2]][path_index[x2]].goal_i, V[task_order[x2]][path_index[x2]].goal_j,
        V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist3 + dist4 - same;
    }
  }
  else {
    return;
  }

  int diffScore = beforeScore - afterScore;
  double prob = exp((double)diffScore / temperature);
  //if (prob > rand_01()) {
  if (diffScore >= 0) {
    swap(task_order[x1], task_order[x2]);
    swap(path_index[x1], path_index[x2]);
  }
}

// ID変更
void change_path_id(double temperature)
{
  int x = rand_uint32() % task_count;
  int y = rand_uint32() % 10;
  if (V[task_order[x]].size() <= y) return;
  if (x == 0 || x == task_count - 1) return;
  if (y == path_index[x]) return;
  int x11 = x - 1;
  int x12 = x + 1;

  int beforeScore = 0;
  int afterScore = 0;
  {
    int dist1 = manhattan_distance(V[task_order[x]][path_index[x]].start_i, V[task_order[x]][path_index[x]].start_j,
      V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
    int dist2 = manhattan_distance(V[task_order[x]][path_index[x]].goal_i, V[task_order[x]][path_index[x]].goal_j,
      V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    beforeScore = dist1 + dist2 - same + V[task_order[x]][path_index[x]].path_cost;
  }
  {
    int dist1 = manhattan_distance(V[task_order[x]][y].start_i, V[task_order[x]][y].start_j,
      V[task_order[x11]][path_index[x11]].goal_i, V[task_order[x11]][path_index[x11]].goal_j);
    int dist2 = manhattan_distance(V[task_order[x]][y].goal_i, V[task_order[x]][y].goal_j,
      V[task_order[x12]][path_index[x12]].start_i, V[task_order[x12]][path_index[x12]].start_j);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    afterScore = dist1 + dist2 - same + V[task_order[x]][y].path_cost;
  }

  int diffScore = beforeScore - afterScore;
  double prob = exp((double)diffScore / temperature);
  //if (prob > rand_01()) {
  if (diffScore >= 0) {
    path_index[x] = y;
  }
}

void simulated_annealing()
{
  int loop = 0;
  double startTemperature = 2;
  double endTemperature = 0;
  double nowProgress = 0;
  while (false) {
    loop++;
    if (loop % 100 == 0) {
      double nowTime = get_elapsed_time();
      nowProgress = nowTime / time_limit_sec;
      if (nowProgress > 1.0) break;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    int ra = rand_uint32() % 100;
    if (ra < 50) {
      // 2点スワップ
      two_swap(temperature);
    }
    else {
      change_path_id(temperature);
    }
  }

  if (mode != 0) {
    cout << "loop = " << loop << endl;
  }
}

ll solve_single_case(int probNum)
{
  start_timer();

  // 複数ケース回すときに内部状態を初期値に戻す
  reset_global_state();

  // 入力受け取り
  read_input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  open_output_file(probNum, ofs);

  // dp
  {
    int dp[5][40][40];
    int dp2[5][40][40];
    rep(i, task_count)
    {
      rep(j, 5)
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
      rep(j, 5)
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
          drep(l, 5)
          {
            int y = task_colors[i][l];
            tmp.push_back(cells_by_color[y][now]);
            if (l == 0) break;
            now = dp2[l][j][now];
          }
          reverse(tmp.begin(), tmp.end());
          path.cells = tmp;
          V[i].push_back(path);
        }
      }

      sort(V[i].begin(), V[i].end());
      baseVal[i] = V[i][0].path_cost;
    }
  }

  // 初期解生成
  build_initial_solution();

  simulated_annealing();

  // 解答を出力
  write_answer(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calc_total_score();
  }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    rand_uint32();
  }

  mode = 1;

  if (mode == 0) {
    solve_single_case(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = solve_single_case(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
