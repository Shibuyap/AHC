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

const int n = 50;
const int nn = n * n;
int si, sj;
int f[60][60];
int value_grid[60][60];

int visited[n * n];
int visited_counter;

class Answer
{
public:
  Answer() : length(0), score(0) {}

  int length;
  int direction[nn];
  int x[nn];
  int y[nn];
  int score;

  void init(int start_x, int start_y)
  {
    x[0] = start_x;
    y[0] = start_y;
    score = value_grid[x[0]][y[0]];
    length = 1;
  }

  void add(int d)
  {
    x[length] = x[length - 1] + dx[d];
    y[length] = y[length - 1] + dy[d];
    direction[length - 1] = d;
    score += value_grid[x[length]][y[length]];
    length++;
  }

  void copy(const Answer& a)
  {
    length = a.length;
    score = a.score;
    rep(i, length) {
      direction[i] = a.direction[i];
      x[i] = a.x[i];
      y[i] = a.y[i];
    }
  }

  bool operator<(const Answer& other) const {
    return score < other.score;
  }

  bool operator>(const Answer& other) const {
    return score > other.score;
  }
};

void swap_answer(Answer& a, Answer& b)
{
  int len = max(a.length, b.length);
  for (int i = 0; i < len; i++) {
    swap(a.direction[i], b.direction[i]);
    swap(a.x[i], b.x[i]);
    swap(a.y[i], b.y[i]);
  }
  swap(a.length, b.length);
  swap(a.score, b.score);
}

Answer current_answer;

const int KOUHO_SIZE = 100;
Answer kouho[KOUHO_SIZE];

Answer best_answer;

bool is_out_of_range(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

static void input_data() {
  string fileNameIfs = "1120.txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs.is_open()) { // 標準入力する
    cin >> si >> sj;
    rep(i, n)
    {
      rep(j, n)
      {
        cin >> f[i][j];
      }
    }
    rep(i, n)
    {
      rep(j, n)
      {
        cin >> value_grid[i][j];
      }
    }
  }
  else { // ファイル入力する
    ifs >> si >> sj;
    rep(i, n)
    {
      rep(j, n)
      {
        ifs >> f[i][j];
      }
    }
    rep(i, n)
    {
      rep(j, n)
      {
        ifs >> value_grid[i][j];
      }
    }
  }
}

void random_walk(Answer& answer) {
  int x = answer.x[answer.length - 1];
  int y = answer.y[answer.length - 1];
  while (true) {
    int ra = Rand() % 24;
    bool ok = false;
    rep(i, 4)
    {
      int nx = x + dx[next_directions[ra][i]];
      int ny = y + dy[next_directions[ra][i]];
      if (!is_out_of_range(nx, ny) && visited[f[nx][ny]] != visited_counter) {
        visited[f[nx][ny]] = visited_counter;
        answer.add(next_directions[ra][i]);
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
  visited_counter++;
  visited[f[si][sj]] = visited_counter;
}

void init_visited(const Answer& answer)
{
  visited_counter++;
  rep(i, answer.length) {
    visited[f[answer.x[i]][answer.y[i]]] = visited_counter;
  }
}

int main()
{
  input_data();

  start_timer();

  rep(i, KOUHO_SIZE) {
    kouho[i].init(si, sj);
  }

  int loop1 = 0;
  // 初期解
  rep(i, 100000)
  {
    loop1++;

    init_visited();

    current_answer.init(si, sj);
    random_walk(current_answer);

    int num = i % KOUHO_SIZE;
    if (current_answer.score > kouho[num].score) {
      kouho[num].copy(current_answer);
    }

    if (get_elapsed_time() > 0.5)
    {
      break;
    }
  }

  int loop2 = 0;
  while (true) {
    loop2++;

    init_visited();

    current_answer.init(si, sj);

    int num = Rand() % KOUHO_SIZE;

    int m = rand() % kouho[num].length;
    rep(i, m) {
      current_answer.add(kouho[num].direction[i]);
      int x = current_answer.x[current_answer.length - 1];
      int y = current_answer.y[current_answer.length - 1];
      visited[f[x][y]] = visited_counter;
    }

    random_walk(current_answer);

    if (current_answer.score > kouho[num].score) {
      kouho[num].copy(current_answer);
    }

    if (get_elapsed_time() > 1.0)
    {
      break;
    }
  }

  sort(kouho, kouho + KOUHO_SIZE, greater<Answer>());
  best_answer.copy(kouho[0]);

  bool is_sorted = false;

  Answer before_keep_path, keep_path, after_keep_path;
  int loop3 = 0;
  while (true) {
    loop3++;

    int num = is_sorted ? 0 : Rand() % 10;

    init_visited();

    int m = kouho[num].length;
    int left = Rand() % (m - 40) + 10;
    int right = left + 1 + Rand() % 20;

    int sx = kouho[num].x[left];
    int sy = kouho[num].y[left];
    int gx = kouho[num].x[right];
    int gy = kouho[num].y[right];

    before_keep_path.init(si, sj);
    rep(i, left) {
      before_keep_path.add(kouho[num].direction[i]);
      int x = before_keep_path.x[before_keep_path.length - 1];
      int y = before_keep_path.y[before_keep_path.length - 1];
      visited[f[x][y]] = visited_counter;
    }

    keep_path.init(sx, sy);
    srep(i, left, right)
    {
      keep_path.add(kouho[num].direction[i]);
      int x = keep_path.x[keep_path.length - 1];
      int y = keep_path.y[keep_path.length - 1];
      visited[f[x][y]] = -1;
    }

    after_keep_path.init(gx, gy);
    srep(i, right, m - 1)
    {
      after_keep_path.add(kouho[num].direction[i]);
      int x = after_keep_path.x[after_keep_path.length - 1];
      int y = after_keep_path.y[after_keep_path.length - 1];
      visited[f[x][y]] = visited_counter;
    }

    Answer new_path;
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
            if (visited[f[nx][ny]] == visited_counter) {
              ok = false;
              break;
            }
            else {
              x = nx;
              y = ny;
              new_path.add(next_directions[ra][i]);
              visited[f[x][y]] = visited_counter;
              ok = true;
              break;
            }
          }
          if (!is_out_of_range(nx, ny) && visited[f[nx][ny]] != visited_counter) {
            x = nx;
            y = ny;
            new_path.add(next_directions[ra][i]);
            visited[f[x][y]] = visited_counter;
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
        visited[f[x][y]] = -1;
      }
    }

    if (before_keep_path.score + keep_path.score + after_keep_path.score > kouho[num].score) {

      kouho[num].copy(before_keep_path);
      rep(i, keep_path.length - 1) {
        kouho[num].add(keep_path.direction[i]);
      }
      rep(i, after_keep_path.length - 1) {
        kouho[num].add(after_keep_path.direction[i]);
      }

      if (kouho[num].score > best_answer.score) {
        best_answer.copy(kouho[num]);
      }
    }

    if (get_elapsed_time() > 1.9)
    {
      break;
    }
    if (!is_sorted && get_elapsed_time() > 1.5)
    {
      sort(kouho, kouho + 10, greater<Answer>());
      is_sorted = true;
    }
  }

  cerr << "loop1 = " << loop1 << endl;
  cerr << "loop2 = " << loop2 << endl;
  cerr << "loop3 = " << loop3 << endl;
  cerr << "best_score = " << best_answer.score << endl;

  string best_string;
  rep(i, best_answer.length - 1) {
    best_string += next_char[best_answer.direction[i]];
  }
  cout << best_string << endl;
  return 0;
}
