#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <utility>

#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 200005

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

static uint32_t rand_xorshift32()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

static double rand_unit_double()
{
  return (rand_xorshift32() + 0.5) * (1.0 / UINT_MAX);
}

void shuffle_array(int* arr, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = rand_xorshift32() % (i + 1);
    int swa = arr[i];
    arr[i] = arr[j];
    arr[j] = swa;
  }
}

int next_directions[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

int dx[4] = { -1, 0, 1, 0 };
int dy[4] = { 0, -1, 0, 1 };
char next_char[4] = { 'U','L','D','R' };

int exec_mode = 0;
bool show_log = true;

const int GRID_SIZE = 50;
const int MAX_PATH_LENGTH = GRID_SIZE * GRID_SIZE;
int si, sj;
int tile_id[GRID_SIZE][GRID_SIZE];
int cell_value[GRID_SIZE][GRID_SIZE];

int visited[GRID_SIZE * GRID_SIZE];
int visited_version;

class Path
{
public:
  Path() : length(0), score(0)
  {
  }

  int length;
  int direction[MAX_PATH_LENGTH];
  int x[MAX_PATH_LENGTH];
  int y[MAX_PATH_LENGTH];
  int score;

  void init(int start_x, int start_y)
  {
    x[0] = start_x;
    y[0] = start_y;
    score = cell_value[x[0]][y[0]];
    length = 1;
  }

  void add(int d)
  {
    x[length] = x[length - 1] + dx[d];
    y[length] = y[length - 1] + dy[d];
    direction[length - 1] = d;
    score += cell_value[x[length]][y[length]];
    length++;
  }

  void pop()
  {
    length--;
    score -= cell_value[x[length]][y[length]];
  }

  void reverse()
  {
    for (int i = 0; i < (length / 2); ++i) {
      swap(x[i], x[length - 1 - i]);
      swap(y[i], y[length - 1 - i]);
    }
    for (int i = 0; i < (length - 1); ++i) {
      for (int j = 0; j < (4); ++j) {
        int nx = x[i] + dx[j];
        int ny = y[i] + dy[j];
        if (nx == x[i + 1] && ny == y[i + 1]) {
          direction[i] = j;
          break;
        }
      }
    }
  }

  void copy(const Path& a)
  {
    length = a.length;
    score = a.score;
    for (int i = 0; i < (length); ++i) {
      direction[i] = a.direction[i];
      x[i] = a.x[i];
      y[i] = a.y[i];
    }
  }

  Path& operator=(const Path& a)
  {
    if (this == &a) return *this;
    length = a.length;
    score = a.score;
    for (int i = 0; i < (length); ++i) {
      direction[i] = a.direction[i];
      x[i] = a.x[i];
      y[i] = a.y[i];
    }
    return *this;
  }

  bool operator<(const Path& other) const
  {
    return score < other.score;
  }

  bool operator>(const Path& other) const
  {
    return score > other.score;
  }
};

Path best_path;

bool is_out_of_bounds(int x, int y)
{
  if (x < 0 || GRID_SIZE <= x || y < 0 || GRID_SIZE <= y) return true;
  return false;
}

static void read_case_input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) { // 標準入力する
    cin >> si >> sj;
    for (int i = 0; i < (GRID_SIZE); ++i) {
      for (int j = 0; j < (GRID_SIZE); ++j) {
        cin >> tile_id[i][j];
      }
    }
    for (int i = 0; i < (GRID_SIZE); ++i) {
      for (int j = 0; j < (GRID_SIZE); ++j) {
        cin >> cell_value[i][j];
      }
    }
  }
  else { // ファイル入力する
    ifs >> si >> sj;
    for (int i = 0; i < (GRID_SIZE); ++i) {
      for (int j = 0; j < (GRID_SIZE); ++j) {
        ifs >> tile_id[i][j];
      }
    }
    for (int i = 0; i < (GRID_SIZE); ++i) {
      for (int j = 0; j < (GRID_SIZE); ++j) {
        ifs >> cell_value[i][j];
      }
    }
  }
}

void write_case_output(int case_num)
{
  string best_string;
  for (int i = 0; i < (best_path.length - 1); ++i) {
    best_string += next_char[best_path.direction[i]];
  }

  if (exec_mode == 0) {
    // 標準出力
    cout << best_string << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << best_string << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

void reset_visited_all()
{
  visited_version++;
  visited[tile_id[si][sj]] = visited_version;
}

void mark_path_as_visited(const Path& path)
{
  visited_version++;
  for (int i = 0; i < (path.length); ++i) {
    visited[tile_id[path.x[i]][path.y[i]]] = visited_version;
  }
}

void extend_path_randomly(Path& path)
{
  int x = path.x[path.length - 1];
  int y = path.y[path.length - 1];
  while (true) {
    int ra = rand_xorshift32() % 24;
    bool ok = false;
    for (int i = 0; i < (4); ++i) {
      int nx = x + dx[next_directions[ra][i]];
      int ny = y + dy[next_directions[ra][i]];
      if (!is_out_of_bounds(nx, ny) && visited[tile_id[nx][ny]] != visited_version) {
        visited[tile_id[nx][ny]] = visited_version;
        path.add(next_directions[ra][i]);
        x = nx;
        y = ny;
        ok = true;
        break;
      }
    }
    if (!ok) break;
  }
}

bool try_connect_path(Path& path, int sx, int sy, int gx, int gy)
{
  path.init(sx, sy);
  int x = sx;
  int y = sy;
  while (x != gx || y != gy) {
    int ra = rand_xorshift32() % 24;
    bool ok = false;
    for (int i = 0; i < (4); ++i) {
      int nx = x + dx[next_directions[ra][i]];
      int ny = y + dy[next_directions[ra][i]];
      if ((nx == gx && ny == gy) || (!is_out_of_bounds(nx, ny) && visited[tile_id[nx][ny]] != visited_version)) {
        x = nx;
        y = ny;
        path.add(next_directions[ra][i]);
        visited[tile_id[x][y]] = visited_version;
        ok = true;
        break;
      }
    }
    if (!ok) break;
  }
  return x == gx && y == gy;
}

// 初期解
void build_from_best_prefix(double time_limit_1, double time_limit_2)
{
  Path current_path;
  int iter = 0;
  while (true) {
    iter++;
    reset_visited_all();
    current_path.init(si, sj);

    if (iter > 100000 || get_elapsed_time() > time_limit_1) {
      int prefix_len = rand_xorshift32() % best_path.length;
      for (int i = 0; i < (prefix_len); ++i) {
        current_path.add(best_path.direction[i]);
        int x = current_path.x[current_path.length - 1];
        int y = current_path.y[current_path.length - 1];
        visited[tile_id[x][y]] = visited_version;
      }
    }

    extend_path_randomly(current_path);

    if (current_path.score > best_path.score) {
      best_path.copy(current_path);
    }

    if (get_elapsed_time() > time_limit_2) {
      break;
    }
  }
}

struct AnnealParam
{
  double time_limit;
  double start_temp;
  double end_temp;
  double score_scale;
  int min_segment_len;
  int range_segment_len;
};

class DfsSolver
{
public:
  bool connected;
  Path best_path;
  Path path;
  int goal_x;
  int goal_y;
  int remaining_search_cnt;

  void start(int sx, int sy, int gx, int gy)
  {
    goal_x = gx;
    goal_y = gy;
    connected = false;
    path.init(sx, sy);
    best_path.init(sx, sy);
    remaining_search_cnt = 500;
    dfs(sx, sy);
  }

  void dfs(int x, int y)
  {
    int legal_dir[4];
    if (remaining_search_cnt == 0) {
      return;
    }
    remaining_search_cnt--;

    int cnt = 0;
    for (int i = 0; i < (4); ++i) {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (nx == goal_x && ny == goal_y) {
        connected = true;
        path.add(i);
        if (path.score > best_path.score) {
          best_path.copy(path);
        }
        path.pop();
      }
      else  if (!is_out_of_bounds(nx, ny) && visited[tile_id[nx][ny]] != visited_version) {
        legal_dir[cnt] = i;
        cnt++;
      }
    }

    shuffle_array(legal_dir, cnt);
    for (int ii = 0; ii < (cnt); ++ii) {
      int i = legal_dir[ii];
      int nx = x + dx[i];
      int ny = y + dy[i];

      path.add(i);
      visited[tile_id[nx][ny]] = visited_version;
      dfs(nx, ny);
      visited[tile_id[nx][ny]] = -1;
      path.pop();
    }
  }
};

bool search_best_path_2(int seg_start_x, int seg_start_y, int seg_goal_x, int seg_goal_y, Path& keep_path, const AnnealParam& param)
{
  DfsSolver dfsSolver;
  int ra = rand_xorshift32() % 2;
  if (ra == 0) {
    dfsSolver.start(seg_start_x, seg_start_y, seg_goal_x, seg_goal_y);
  }
  else {
    dfsSolver.start(seg_goal_x, seg_goal_y, seg_start_x, seg_start_y);
  }

  if (dfsSolver.connected) {
    keep_path.copy(dfsSolver.best_path);
    if (ra == 1) {
      keep_path.reverse();
    }
    return true;
  }

  return false;
}

void anneal_path_segment(const AnnealParam& param)
{
  Path current_path;
  current_path.copy(best_path);

  // ── 区間 [seg_left, seg_right) を焼き鈍して再生成 ──
  Path before_keep_path, keep_path, after_keep_path;

  int iteration_cnt = 0;

  double anneal_start_time = get_elapsed_time();

  while (true) {
    ++iteration_cnt;

    reset_visited_all();

    int path_len = current_path.length;
    int segment_len = rand_xorshift32() % param.range_segment_len + param.min_segment_len;
    int seg_left = rand_xorshift32() % (path_len - segment_len);
    int seg_right = seg_left + segment_len;

    int seg_start_x = current_path.x[seg_left];
    int seg_start_y = current_path.y[seg_left];
    int seg_goal_x = current_path.x[seg_right];
    int seg_goal_y = current_path.y[seg_right];

    /* ─ before_keep_path ─ */
    before_keep_path.init(si, sj);
    for (int i = 0; i < (seg_left); ++i) {
      before_keep_path.add(current_path.direction[i]);
      int cx = before_keep_path.x[before_keep_path.length - 1];
      int cy = before_keep_path.y[before_keep_path.length - 1];
      visited[tile_id[cx][cy]] = visited_version;
    }

    /* ─ keep_path (旧区間) ─ */
    keep_path.init(seg_start_x, seg_start_y);
    srep(i, seg_left, seg_right)
    {
      keep_path.add(current_path.direction[i]);
      int cx = keep_path.x[keep_path.length - 1];
      int cy = keep_path.y[keep_path.length - 1];
      visited[tile_id[cx][cy]] = -1;
    }

    visited[tile_id[seg_goal_x][seg_goal_y]] = visited_version;

    /* ─ after_keep_path ─ */
    after_keep_path.init(seg_goal_x, seg_goal_y);
    srep(i, seg_right, path_len - 1)
    {
      after_keep_path.add(current_path.direction[i]);
      int cx = after_keep_path.x[after_keep_path.length - 1];
      int cy = after_keep_path.y[after_keep_path.length - 1];
      visited[tile_id[cx][cy]] = visited_version;
    }

    bool res = search_best_path_2(seg_start_x, seg_start_y, seg_goal_x, seg_goal_y, keep_path, param);
    if (!res) {
      continue;
    }

    double now_time = get_elapsed_time();
    double progress_ratio = (now_time - anneal_start_time) / (param.time_limit - anneal_start_time);
    double temp = param.start_temp + (param.end_temp - param.start_temp) * progress_ratio;

    double new_score = before_keep_path.score + keep_path.score + after_keep_path.score
      - (cell_value[seg_start_x][seg_start_y] + cell_value[seg_goal_x][seg_goal_y]);
    double diff_score = (new_score - current_path.score) * param.score_scale;
    double accept_prob = exp(diff_score / temp);

    if (accept_prob > rand_unit_double()) {
      /* accept */
      current_path.copy(before_keep_path);
      for (int i = 0; i < (keep_path.length - 1); ++i) current_path.add(keep_path.direction[i]);
      for (int i = 0; i < (after_keep_path.length - 1); ++i) current_path.add(after_keep_path.direction[i]);

      //for (int i = 0; i < (current_path.length); ++i) {
      //  int x = current_path.x[i];
      //  int y = current_path.y[i];
      //  if (is_out_of_bounds(x, y)) {
      //    cerr << "out of bounds " << seg_left << ' ' << seg_right << ' ' << current_path.length << endl;
      //  }
      //}

      if (current_path.score > best_path.score) {
        best_path.copy(current_path);
      }
    }

    if (now_time > param.time_limit) break;
  }
  if (show_log) {
    cerr << "anneal_path_segment : iter = " << iteration_cnt << endl;
  }
}

int solve_case(int case_num, const AnnealParam& param)
{
  start_timer();

  read_case_input(case_num);

  best_path.init(si, sj);

  const double TIME_LIMIT_STAGE1 = 0.05;
  const double TIME_LIMIT_STAGE2 = 0.1;
  build_from_best_prefix(TIME_LIMIT_STAGE1, TIME_LIMIT_STAGE2);

  anneal_path_segment(param);

  build_from_best_prefix(0.1, 1.95);

  if (show_log) {
    cerr << "best_score = " << best_path.score << ", time = " << get_elapsed_time() << "sec." << endl;
  }

  write_case_output(case_num);

  return best_path.score;
}

int main()
{
  exec_mode = 1;
  show_log = true;

  AnnealParam param;
  param.time_limit = 1.8;
  param.start_temp = 20000048.0;
  param.end_temp = 0.1;
  param.score_scale = 12345.6;
  param.range_segment_len = 50;
  param.min_segment_len = 3;

  if (exec_mode == 0) {
    solve_case(0, param);
  }
  else if (exec_mode == 1) {
    int sum = 0;
    for (int i = 0; i < (10); ++i) {
      sum += solve_case(i, param);
      cout << i << ' ' << sum << endl;
    }
    cout << "sum = " << sum << endl;
  }
  else if (exec_mode == 2) {
    show_log = false;

    int best_score = 0;
    AnnealParam best_param = param;

    int iter = 0;
    while (true) {
      iter++;
      int new_score = 0;
      AnnealParam new_param = best_param;

      if (iter > 10) {
        int ra = rand_xorshift32() % 5;
        switch (ra) {
          case 0:
            new_param.start_temp *= rand_unit_double() * 2;
            break;
          case 1:
            new_param.score_scale *= rand_unit_double() * 2;
            break;
          case 3:
            new_param.range_segment_len *= rand_unit_double() * 2;
            new_param.range_segment_len = max(new_param.range_segment_len, 5);
            break;
          case 4:
            new_param.min_segment_len *= rand_unit_double() * 2;
            new_param.min_segment_len = max(new_param.range_segment_len, 0);
            break;
        }
      }

      for (int i = 0; i < (1); ++i) {
        new_score += solve_case(i, new_param);
      }

      if (new_score >= best_score) {
        best_score = new_score;
        best_param = new_param;

        cout << "iter = " << iter << ", score = " << best_score << endl;
        cout << "param.time_limit = " << best_param.time_limit << endl;
        cout << setprecision(10) << "param.start_temp = " << best_param.start_temp << endl;
        cout << "param.end_temp = " << best_param.end_temp << endl;
        cout << "param.score_scale = " << best_param.score_scale << endl;
        cout << "param.range_segment_len = " << best_param.range_segment_len << endl;
        cout << "param.min_segment_len = " << best_param.min_segment_len << endl;
        cout << endl;
      }
    }
  }
}
