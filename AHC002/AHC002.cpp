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

void start_timer()
{
  start_time_clock = std::chrono::steady_clock::now();
}

double get_elapsed_time()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
  return elapsed.count();
}

const double TIME_LIMIT_STAGE1 = 0.2;
const double TIME_LIMIT_STAGE2 = 0.5;
const double TIME_LIMIT_STAGE3 = 1.0;
const double TIME_LIMIT_FINAL = 1.9;

static uint32_t Rand()
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

static double Rand01()
{
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}

int next_directions[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

int dx[4] = { -1, 0, 1, 0 };
int dy[4] = { 0, -1, 0, 1 };
char next_char[4] = { 'U','L','D','R' };

const int GRID_SIZE = 50;
const int MAX_PATH_LENGTH = GRID_SIZE * GRID_SIZE;
int si, sj;
int tile_id[60][60];
int cell_value[60][60];

int visited[GRID_SIZE * GRID_SIZE];
int visited_version;

class Path
{
public:
  Path() : length(0), score(0) {}

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

  void reverse()
  {
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

  void copy(const Path& a)
  {
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

Path current_path;

const int CANDIDATE_SIZE = 1;
const int CANDIDATE_SIZE_2 = 1;
Path candidate_paths[CANDIDATE_SIZE];

Path best_path;

bool is_out_of_range(int x, int y)
{
  if (x < 0 || GRID_SIZE <= x || y < 0 || GRID_SIZE <= y) return true;
  return false;
}

static void input_data() {
  string fileNameIfs = "1120.txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs.is_open()) { // 標準入力する
    cin >> si >> sj;
    rep(i, GRID_SIZE)
    {
      rep(j, GRID_SIZE)
      {
        cin >> tile_id[i][j];
      }
    }
    rep(i, GRID_SIZE)
    {
      rep(j, GRID_SIZE)
      {
        cin >> cell_value[i][j];
      }
    }
  }
  else { // ファイル入力する
    ifs >> si >> sj;
    rep(i, GRID_SIZE)
    {
      rep(j, GRID_SIZE)
      {
        ifs >> tile_id[i][j];
      }
    }
    rep(i, GRID_SIZE)
    {
      rep(j, GRID_SIZE)
      {
        ifs >> cell_value[i][j];
      }
    }
  }

  if (ifs.is_open()) {
    ifs.close();
  }
}

void extend_random_path(Path& path) {
  int x = path.x[path.length - 1];
  int y = path.y[path.length - 1];
  while (true) {
    int ra = Rand() % 24;
    bool ok = false;
    rep(i, 4)
    {
      int nx = x + dx[next_directions[ra][i]];
      int ny = y + dy[next_directions[ra][i]];
      if (!is_out_of_range(nx, ny) && visited[tile_id[nx][ny]] != visited_version) {
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

void init_visited()
{
  visited_version++;
  visited[tile_id[si][sj]] = visited_version;
}

void init_visited(const Path& path)
{
  visited_version++;
  rep(i, path.length) {
    visited[tile_id[path.x[i]][path.y[i]]] = visited_version;
  }
}

bool attempt_connect(Path& path, int sx, int sy, int gx, int gy) {
  path.init(sx, sy);
  int x = sx;
  int y = sy;

  while (x != gx || y != gy) {
    int ra = Rand() % 24;
    bool ok = false;

    rep(i, 4)
    {
      int nx = x + dx[next_directions[ra][i]];
      int ny = y + dy[next_directions[ra][i]];

      if ((nx == gx && ny == gy) || (!is_out_of_range(nx, ny) && visited[tile_id[nx][ny]] != visited_version)) {
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

int solve()
{
  input_data();

  start_timer();

  rep(i, CANDIDATE_SIZE) {
    candidate_paths[i].init(si, sj);
  }
  best_path.init(si, sj);

  int loop1 = 0;
  // 初期解
  rep(i, 100000)
  {
    loop1++;

    init_visited();

    current_path.init(si, sj);
    extend_random_path(current_path);

    int num = i % CANDIDATE_SIZE;
    if (current_path.score > candidate_paths[num].score) {
      candidate_paths[num].copy(current_path);
    }

    if (get_elapsed_time() > TIME_LIMIT_STAGE1)
    {
      break;
    }
  }

  int loop2 = 0;
  while (true) {
    loop2++;

    init_visited();

    current_path.init(si, sj);

    int num = Rand() % CANDIDATE_SIZE;

    int m = Rand() % candidate_paths[num].length;
    rep(i, m) {
      current_path.add(candidate_paths[num].direction[i]);
      int x = current_path.x[current_path.length - 1];
      int y = current_path.y[current_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    extend_random_path(current_path);

    if (current_path.score > candidate_paths[num].score) {
      candidate_paths[num].copy(current_path);
    }

    if (get_elapsed_time() > TIME_LIMIT_STAGE2)
    {
      break;
    }
  }

  sort(candidate_paths, candidate_paths + CANDIDATE_SIZE, greater<Path>());
  best_path.copy(candidate_paths[0]);

  bool is_sorted = false;

  // [left, right)区間の経路を再生成するため、
  // それ以前を before_keep_path、該当区間を keep_path、
  // それ以降を after_keep_path として管理。
  Path before_keep_path, keep_path, after_keep_path;
  int loop3 = 0;
  const double START_TEMP = 4000048.0;
  const double END_TEMP = 0.1;
  const double SCORE_SCALE = 12345.6;
  double syori3_start_time = get_elapsed_time();
  while (true) {
    loop3++;

    int num = is_sorted ? 0 : Rand() % CANDIDATE_SIZE_2;

    init_visited();

    int m = candidate_paths[num].length;
    int len = Rand() % 40 + 3;
    int left = Rand() % (m - len);
    int right = left + len;

    int sx = candidate_paths[num].x[left];
    int sy = candidate_paths[num].y[left];
    int gx = candidate_paths[num].x[right];
    int gy = candidate_paths[num].y[right];

    before_keep_path.init(si, sj);
    rep(i, left) {
      before_keep_path.add(candidate_paths[num].direction[i]);
      int x = before_keep_path.x[before_keep_path.length - 1];
      int y = before_keep_path.y[before_keep_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    keep_path.init(sx, sy);
    srep(i, left, right)
    {
      keep_path.add(candidate_paths[num].direction[i]);
      int x = keep_path.x[keep_path.length - 1];
      int y = keep_path.y[keep_path.length - 1];
      visited[tile_id[x][y]] = -1;
    }

    visited[tile_id[gx][gy]] = visited_version;

    after_keep_path.init(gx, gy);
    srep(i, right, m - 1)
    {
      after_keep_path.add(candidate_paths[num].direction[i]);
      int x = after_keep_path.x[after_keep_path.length - 1];
      int y = after_keep_path.y[after_keep_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    bool is_first = true;
    Path new_path;
    rep(_, 100)
    {
      int is_reverse = Rand() % 2;
      bool is_connect = false;
      if (is_reverse == 0) {
        is_connect = attempt_connect(new_path, sx, sy, gx, gy);
      }
      else {
        is_connect = attempt_connect(new_path, gx, gy, sx, sy);
      }

      if (is_connect && (is_first || new_path.score > keep_path.score)) {
        if (is_reverse == 1) {
          new_path.reverse();
        }
        keep_path.copy(new_path);
        is_first = false;
      }
      rep(i, new_path.length)
      {
        int x = new_path.x[i];
        int y = new_path.y[i];
        if (x == sx && y == sy)continue;
        if (x == gx && y == gy)continue;
        visited[tile_id[x][y]] = -1;
      }
    }

    double now_time = get_elapsed_time();

    double progress_ratio = (now_time - syori3_start_time) / (TIME_LIMIT_FINAL - syori3_start_time);

    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;
    double new_score = before_keep_path.score + keep_path.score + after_keep_path.score - (cell_value[sx][sy] + cell_value[gx][gy]);
    double diff_score = (new_score - candidate_paths[num].score) * SCORE_SCALE;
    double prob = exp(diff_score / temp);
    if (prob > Rand01()) {
      candidate_paths[num].copy(before_keep_path);
      rep(i, keep_path.length - 1) {
        candidate_paths[num].add(keep_path.direction[i]);
      }
      rep(i, after_keep_path.length - 1) {
        candidate_paths[num].add(after_keep_path.direction[i]);
      }

      if (candidate_paths[num].score > best_path.score) {
        best_path.copy(candidate_paths[num]);
      }
    }

    if (!is_sorted && now_time > TIME_LIMIT_STAGE3)
    {
      sort(candidate_paths, candidate_paths + CANDIDATE_SIZE_2, greater<Path>());
      is_sorted = true;
    }

    if (now_time > TIME_LIMIT_FINAL)
    {
      break;
    }
  }

  cerr << "loop1 = " << loop1 << endl;
  cerr << "loop2 = " << loop2 << endl;
  cerr << "loop3 = " << loop3 << endl;
  cerr << "best_score = " << best_path.score << endl;

  string best_string;
  rep(i, best_path.length - 1) {
    best_string += next_char[best_path.direction[i]];
  }
  cout << best_string << endl;
  return best_path.score;
}

int main() {
  int exec_mode = 0;

  if (exec_mode == 0) {
    solve();
  }
  else if (exec_mode == 1) {
    rep(i, 5) {
      solve();
    }
  }
}
