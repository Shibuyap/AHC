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

#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 200005

std::chrono::steady_clock::time_point start_time_clock;

void start_timer() {
  start_time_clock = std::chrono::steady_clock::now();
}

double get_elapsed_time() {
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
  return elapsed.count();
}

static uint32_t rand_xorshift32() {
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

static double rand_unit_double() {
  return (rand_xorshift32() + 0.5) * (1.0 / UINT_MAX);
}

int next_directions[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

int dx[4] = { -1, 0, 1, 0 };
int dy[4] = { 0, -1, 0, 1 };
char next_char[4] = { 'U','L','D','R' };

int exec_mode = 0;

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
  Path() : length(0), score(0) {
  }

  int length;
  int direction[MAX_PATH_LENGTH];
  int x[MAX_PATH_LENGTH];
  int y[MAX_PATH_LENGTH];
  int score;

  void init(int start_x, int start_y) {
    x[0] = start_x;
    y[0] = start_y;
    score = cell_value[x[0]][y[0]];
    length = 1;
  }

  void add(int d) {
    x[length] = x[length - 1] + dx[d];
    y[length] = y[length - 1] + dy[d];
    direction[length - 1] = d;
    score += cell_value[x[length]][y[length]];
    length++;
  }

  void reverse() {
    rep(i, length / 2) {
      swap(x[i], x[length - 1 - i]);
      swap(y[i], y[length - 1 - i]);
    }
    rep(i, length - 1) {
      rep(j, 4) {
        int nx = x[i] + dx[j];
        int ny = y[i] + dy[j];
        if (nx == x[i + 1] && ny == y[i + 1]) {
          direction[i] = j;
          break;
        }
      }
    }
  }

  void copy(const Path& a) {
    length = a.length;
    score = a.score;
    rep(i, length) {
      direction[i] = a.direction[i];
      x[i] = a.x[i];
      y[i] = a.y[i];
    }
  }

  Path& operator=(const Path& a) {
    if (this == &a) return *this;
    length = a.length;
    score = a.score;
    rep(i, length) {
      direction[i] = a.direction[i];
      x[i] = a.x[i];
      y[i] = a.y[i];
    }
    return *this;
  }

  bool operator<(const Path& other) const {
    return score < other.score;
  }

  bool operator>(const Path& other) const {
    return score > other.score;
  }
};

Path best_path;

bool is_out_of_bounds(int x, int y) {
  if (x < 0 || GRID_SIZE <= x || y < 0 || GRID_SIZE <= y) return true;
  return false;
}

static void read_case_input(int case_num) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) { // 標準入力する
    cin >> si >> sj;
    rep(i, GRID_SIZE) {
      rep(j, GRID_SIZE) {
        cin >> tile_id[i][j];
      }
    }
    rep(i, GRID_SIZE) {
      rep(j, GRID_SIZE) {
        cin >> cell_value[i][j];
      }
    }
  }
  else { // ファイル入力する
    ifs >> si >> sj;
    rep(i, GRID_SIZE) {
      rep(j, GRID_SIZE) {
        ifs >> tile_id[i][j];
      }
    }
    rep(i, GRID_SIZE) {
      rep(j, GRID_SIZE) {
        ifs >> cell_value[i][j];
      }
    }
  }
}

void write_case_output(int case_num) {
  string best_string;
  rep(i, best_path.length - 1) {
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

void reset_visited_all() {
  visited_version++;
  visited[tile_id[si][sj]] = visited_version;
}

void mark_path_as_visited(const Path& path) {
  visited_version++;
  rep(i, path.length) {
    visited[tile_id[path.x[i]][path.y[i]]] = visited_version;
  }
}

void extend_path_randomly(Path& path) {
  int x = path.x[path.length - 1];
  int y = path.y[path.length - 1];
  while (true) {
    int ra = rand_xorshift32() % 24;
    bool ok = false;
    rep(i, 4) {
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

bool try_connect_path(Path& path, int sx, int sy, int gx, int gy) {
  path.init(sx, sy);
  int x = sx;
  int y = sy;
  while (x != gx || y != gy) {
    int ra = rand_xorshift32() % 24;
    bool ok = false;
    rep(i, 4) {
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
void build_from_best_prefix(double time_limit_1, double time_limit_2) {
  Path current_path;
  int iter = 0;
  while (true) {
    iter++;
    reset_visited_all();
    current_path.init(si, sj);

    if (iter > 100000 || get_elapsed_time() > time_limit_1) {
      int prefix_len = rand_xorshift32() % best_path.length;
      rep(i, prefix_len) {
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

  cerr << "build_from_best_prefix : iter = " << iter << endl;
}

void anneal_path_segment(double time_limit) {
  Path current_path;
  current_path.copy(best_path);

  // [left, right)区間の経路を再生成するため、
  // それ以前を before_keep_path、該当区間を keep_path、
  // それ以降を after_keep_path として管理。
  Path before_keep_path, keep_path, after_keep_path;
  int iter = 0;
  const double START_TEMP = 4000048.0;
  const double END_TEMP = 0.1;
  const double SCORE_SCALE = 12345.6;
  double syori3_start_time = get_elapsed_time();
  while (true) {
    iter++;

    reset_visited_all();

    int m = current_path.length;
    int len = rand_xorshift32() % 40 + 3;
    int left = rand_xorshift32() % (m - len);
    int right = left + len;

    int sx = current_path.x[left];
    int sy = current_path.y[left];
    int gx = current_path.x[right];
    int gy = current_path.y[right];

    before_keep_path.init(si, sj);
    rep(i, left) {
      before_keep_path.add(current_path.direction[i]);
      int x = before_keep_path.x[before_keep_path.length - 1];
      int y = before_keep_path.y[before_keep_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    keep_path.init(sx, sy);
    srep(i, left, right) {
      keep_path.add(current_path.direction[i]);
      int x = keep_path.x[keep_path.length - 1];
      int y = keep_path.y[keep_path.length - 1];
      visited[tile_id[x][y]] = -1;
    }

    visited[tile_id[gx][gy]] = visited_version;

    after_keep_path.init(gx, gy);
    srep(i, right, m - 1) {
      after_keep_path.add(current_path.direction[i]);
      int x = after_keep_path.x[after_keep_path.length - 1];
      int y = after_keep_path.y[after_keep_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    bool first_better_found = true;
    Path new_path;
    rep(_, 100) {
      int is_reverse = rand_xorshift32() % 2;
      bool is_connect = false;
      if (is_reverse == 0) {
        is_connect = try_connect_path(new_path, sx, sy, gx, gy);
      }
      else {
        is_connect = try_connect_path(new_path, gx, gy, sx, sy);
      }

      if (is_connect && (first_better_found || new_path.score > keep_path.score)) {
        if (is_reverse == 1) {
          new_path.reverse();
        }
        keep_path.copy(new_path);
        first_better_found = false;
      }
      rep(i, new_path.length) {
        int x = new_path.x[i];
        int y = new_path.y[i];
        if (x == sx && y == sy)continue;
        if (x == gx && y == gy)continue;
        visited[tile_id[x][y]] = -1;
      }
    }

    double now_time = get_elapsed_time();

    double progress_ratio = (now_time - syori3_start_time) / (time_limit - syori3_start_time);

    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;
    double new_score = before_keep_path.score + keep_path.score + after_keep_path.score - (cell_value[sx][sy] + cell_value[gx][gy]);
    double diff_score = (new_score - current_path.score) * SCORE_SCALE;
    double prob = exp(diff_score / temp);
    if (prob > rand_unit_double()) {
      current_path.copy(before_keep_path);
      rep(i, keep_path.length - 1) {
        current_path.add(keep_path.direction[i]);
      }
      rep(i, after_keep_path.length - 1) {
        current_path.add(after_keep_path.direction[i]);
      }

      if (current_path.score > best_path.score) {
        best_path.copy(current_path);
      }
    }

    if (now_time > time_limit) {
      break;
    }
  }

  cerr << "anneal_path_segment : iter = " << iter << endl;
}

int solve_case(int case_num) {
  start_timer();

  read_case_input(case_num);

  best_path.init(si, sj);

  const double TIME_LIMIT_STAGE1 = 0.2;
  const double TIME_LIMIT_STAGE2 = 0.5;
  const double TIME_LIMIT_FINAL = 1.9;

  build_from_best_prefix(TIME_LIMIT_STAGE1, TIME_LIMIT_STAGE2);

  anneal_path_segment(TIME_LIMIT_FINAL);

  cerr << "best_score = " << best_path.score << endl;

  write_case_output(case_num);

  return best_path.score;
}

int main() {
  exec_mode = 1;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode == 1) {
    rep(i, 5) {
      solve_case(i);
    }
  }
}
