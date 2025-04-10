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

const double TIME_LIMIT_STAGE1 = 0.5;
const double TIME_LIMIT_STAGE2 = 1.0;
const double TIME_LIMIT_STAGE3 = 1.5;
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

const int CANDIDATE_SIZE = 100;
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

int main()
{
  input_data();

  start_timer();

  rep(i, CANDIDATE_SIZE) {
    candidate_paths[i].init(si, sj);
  }

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
  while (true) {
    loop3++;

    int num = is_sorted ? 0 : Rand() % 10;

    init_visited();

    int m = candidate_paths[num].length;
    if (m <= 40) continue;
    int left = Rand() % (m - 40) + 10;
    int right = left + 1 + Rand() % 20;

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

    after_keep_path.init(gx, gy);
    srep(i, right, m - 1)
    {
      after_keep_path.add(candidate_paths[num].direction[i]);
      int x = after_keep_path.x[after_keep_path.length - 1];
      int y = after_keep_path.y[after_keep_path.length - 1];
      visited[tile_id[x][y]] = visited_version;
    }

    Path new_path;
    rep(_, 100)
    {
      new_path.init(sx, sy);
      int x = sx;
      int y = sy;

      while (x != gx || y != gy) {
        int ra = Rand() % 24;
        bool ok = false;

        rep(i, 4)
        {
          int nx = x + dx[next_directions[ra][i]];
          int ny = y + dy[next_directions[ra][i]];
          if (nx == gx && ny == gy) {
            if (visited[tile_id[nx][ny]] == visited_version) {
              ok = false;
              break;
            }
            else {
              x = nx;
              y = ny;
              new_path.add(next_directions[ra][i]);
              visited[tile_id[x][y]] = visited_version;
              ok = true;
              break;
            }
          }
          if (!is_out_of_range(nx, ny) && visited[tile_id[nx][ny]] != visited_version) {
            x = nx;
            y = ny;
            new_path.add(next_directions[ra][i]);
            visited[tile_id[x][y]] = visited_version;
            ok = true;
            break;
          }
        }

        if (!ok) break;
      }

      if (x == gx && y == gy && new_path.score > keep_path.score) {
        keep_path.copy(new_path);
      }
      srep(i, 1, new_path.length)
      {
        int x = new_path.x[i];
        int y = new_path.y[i];
        visited[tile_id[x][y]] = -1;
      }
    }

    if (before_keep_path.score + keep_path.score + after_keep_path.score > candidate_paths[num].score) {
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

    if (!is_sorted && get_elapsed_time() > TIME_LIMIT_STAGE3)
    {
      sort(candidate_paths, candidate_paths + 10, greater<Path>());
      is_sorted = true;
    }

    if (get_elapsed_time() > TIME_LIMIT_FINAL)
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
  return 0;
}
