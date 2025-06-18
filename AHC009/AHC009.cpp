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
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char DC[4] = { 'U', 'L', 'D', 'R' };

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
  static uint32_t Rand()
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

  static double Rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
} // namespace

int mode;

namespace /* 変数 */
{
  // 入力用変数
  const int BOARD_SIZE = 20;
  const int MAX_ROUTE_LEN = 200;
  int sx, sy, tx, ty;
  double forget_prob;
  int h_wall[BOARD_SIZE][BOARD_SIZE];
  int v_wall[BOARD_SIZE][BOARD_SIZE];

  // 解答用変数
  ll score;
  vector<int> route;

  // 焼きなまし用変数
  ll best_score;
  vector<int> best_route;

} // namespace

void input_data(int cn)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << cn << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> sx >> sy >> tx >> ty >> forget_prob;
    for (int i = 0; i < BOARD_SIZE; ++i) {
      for (int j = 0; j < BOARD_SIZE - 1; ++j) {
        char c;
        cin >> c;
        h_wall[i][j] = c - '0';
      }
    }
    for (int i = 0; i < BOARD_SIZE - 1; ++i) {
      for (int j = 0; j < BOARD_SIZE; ++j) {
        char c;
        cin >> c;
        v_wall[i][j] = c - '0';
      }
    }
  }
  else {
    // ファイル入力
    ifs >> sx >> sy >> tx >> ty >> forget_prob;
    for (int i = 0; i < BOARD_SIZE; ++i) {
      for (int j = 0; j < BOARD_SIZE - 1; ++j) {
        char c;
        ifs >> c;
        h_wall[i][j] = c - '0';
      }
    }
    for (int i = 0; i < BOARD_SIZE - 1; ++i) {
      for (int j = 0; j < BOARD_SIZE; ++j) {
        char c;
        ifs >> c;
        v_wall[i][j] = c - '0';
      }
    }
  }
}

void output_data(int cn)
{
  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < static_cast<int>(min(route.size(), static_cast<size_t>(MAX_ROUTE_LEN))); ++i) {
      cout << DC[route[i]];
    }
    cout << '\n';
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << cn << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < static_cast<int>(min(route.size(), static_cast<size_t>(MAX_ROUTE_LEN))); ++i) {
      ofs << DC[route[i]];
    }
    ofs << '\n';
  }
}

// スコア計算
ll simulate_score(vector<int>& vec)
{
  const int TRIALS = 100;
  double total = 0.0;

  for (int trial = 0; trial < TRIALS; ++trial) {
    int len = vec.size();
    int row = sx, col = sy;

    for (int step = 0; step < len; ++step) {
      if (Rand01() < forget_prob) {
        // 忘却：その場に留まる
      }
      else {
        int dir = vec[step];

        if (dir == 0) {
          if (row == 0 || v_wall[row - 1][col])
            continue;
        }
        if (dir == 1) {
          if (col == 0 || h_wall[row][col - 1])
            continue;
        }
        if (dir == 2) {
          if (row == 19 || v_wall[row][col])
            continue;
        }
        if (dir == 3) {
          if (col == 19 || h_wall[row][col])
            continue;
        }

        row += dx[dir];
        col += dy[dir];
        if (row == tx && col == ty) {
          total += 400 - step;
          break;
        }
      }
    }
  }
  return (total / TRIALS) * 250000.0;
}

int Solve(int cn)
{
  start_timer();

  input_data(cn);

  route.clear();
  for (int i = 0; i < MAX_ROUTE_LEN; ++i) {
    route.push_back(Rand() % 4);
  }
  score = simulate_score(route);

  /* ──────────────────────  Dijkstra 前処理  ────────────────────── */
  int cost[BOARD_SIZE][BOARD_SIZE];
  int prev_dir[21][21];
  int seg_len[BOARD_SIZE][BOARD_SIZE];
  for (int r = 0; r < BOARD_SIZE; ++r) {
    for (int c = 0; c < BOARD_SIZE; ++c) {
      cost[r][c] = INF;
      prev_dir[20][20] = -1;
      seg_len[r][c] = -1;
    }
  }

  cost[sx][sy] = 0;
  using Node = pair<int, P>;
  priority_queue<Node, vector<Node>, greater<Node>> pq;
  pq.push({ 0, {sx, sy} });

  while (!pq.empty()) {
    auto node = pq.top();
    pq.pop();
    int row = node.second.first;
    int col = node.second.second;
    if (node.first > cost[row][col])
      continue;
    if (cost[row][col] > 250)
      continue;

    for (int d = 0; d < 4; ++d) {
      int nr = row, nc = col;
      int steps = 0;

      if (d == 0) {
        while (nr != 0 && v_wall[nr - 1][nc] == 0) {
          --nr;
          ++steps;
          if (nr == tx && nc == ty)
            break;
        }
      }
      if (d == 1) {
        while (nc != 0 && h_wall[nr][nc - 1] == 0) {
          --nc;
          ++steps;
          if (nr == tx && nc == ty)
            break;
        }
      }
      if (d == 2) {
        while (nr != 19 && v_wall[nr][nc] == 0) {
          ++nr;
          ++steps;
          if (nr == tx && nc == ty)
            break;
        }
      }
      if (d == 3) {
        while (nc != 19 && h_wall[nr][nc] == 0) {
          ++nc;
          ++steps;
          if (nr == tx && nc == ty)
            break;
        }
      }
      if (nr == row && nc == col)
        continue;

      int newCost = cost[row][col] + steps + forget_prob * 15;
      if (cost[nr][nc] > newCost) {
        cost[nr][nc] = newCost;
        seg_len[nr][nc] = steps;
        prev_dir[nr][nc] = d;
        pq.push({ newCost, {nr, nc} });
      }
    }
  }

  for (int i = 0; i < MAX_ROUTE_LEN + 10; ++i) {
    route.push_back(rand() % 4);
  }

  /* ──────────────────────  最短経路ベース経路生成  ────────────────────── */
  if (cost[tx][ty] <= MAX_ROUTE_LEN) {
    route.clear();
    int row = tx, col = ty;
    while (row != sx || col != sy) {
      int d = prev_dir[row][col];
      int steps = seg_len[row][col];
      int nr = row - dx[d] * steps;
      int nc = col - dy[d] * steps;

      for (int k = 0; k < static_cast<int>(steps / (1.0 - forget_prob) + forget_prob * 15); ++k)
        route.push_back(d);

      row = nr;
      col = nc;
    }
    reverse(route.begin(), route.end());
    while (route.size() < MAX_ROUTE_LEN)
      route.push_back(rand() % 4);
    score = simulate_score(route);
  }

  bool need_random = false;
  best_score = -1;
  if (route.size() > MAX_ROUTE_LEN) {
    bool found = false;
    for (int di = 0; di < 10; ++di) {
      for (int dj = 0; dj < 10; ++dj) {
        route.clear();
        int row = tx - di, col = ty - dj;
        if (cost[row][col] >= MAX_ROUTE_LEN) {
          continue;
        }

        while (row != sx || col != sy) {
          int d = prev_dir[row][col];
          int steps = seg_len[row][col];
          int nr = row - dx[d] * steps;
          int nc = col - dy[d] * steps;

          for (int k = 0; k < static_cast<int>(steps / (1.0 - forget_prob) + forget_prob * 15); ++k) {
            route.push_back(d);
          }

          row = nr;
          col = nc;
        }
        reverse(route.begin(), route.end());
        while (route.size() < MAX_ROUTE_LEN) {
          route.push_back(rand() % 4);
        }

        if (route.size() == MAX_ROUTE_LEN) {
          found = true;
          need_random = true;
          score = simulate_score(route);
          if (score > best_score) {
            best_score = score;
            best_route = route;
          }
        }
      }
    }
    if (!found) {
      route.clear();
      for (int i = 0; i < MAX_ROUTE_LEN + 10; ++i) route.push_back(rand() % 4);
    }
    else {
      route = best_route;
      score = best_score;
      if (score == -1)
        for (int i = 0; i < MAX_ROUTE_LEN + 10; ++i) route.push_back(rand() % 4);
    }
  }

  /* ──────────────────────  焼きなまし  ────────────────────── */
  int iter = 0;
  if (need_random || route.size() > MAX_ROUTE_LEN) {
    if (!need_random) {
      route.clear();
      for (int i = 0; i < MAX_ROUTE_LEN; ++i) route.push_back(Rand() % 4);
    }
    score = simulate_score(route);
    best_route = route;
    best_score = score;

    double elapsed = get_elapsed_time();
    const double timeLimit = 1.9;
    const double tempInit = 2048;
    const double tempFinal = 0.0001;

    while (true) {
      ++iter;
      if (iter % 100 == 1)
        elapsed = get_elapsed_time();
      if (elapsed > timeLimit)
        break;

      int pos = Rand() % MAX_ROUTE_LEN;
      int newDir = Rand() % 4;
      int backup = route[pos];
      route[pos] = newDir;

      int candScore = simulate_score(route);
      int diff = candScore - score;

      double temp = tempInit + (tempFinal - tempInit) * elapsed / timeLimit;
      double prob = exp(diff / temp);

      if (diff > 0 || prob > Rand01()) {
        score = candScore;
        if (score > best_score) {
          best_score = score;
          best_route = route;
        }
      }
      else {
        route[pos] = backup;
      }
    }
    route = best_route;
    score = best_score;
  }

  output_data(cn);

  if (mode != 0) {
    cout << "route.size() = " << route.size() << endl;
    cout << "iter = " << iter << endl;
    cout << score << endl;
    cout << get_elapsed_time() << " sec." << endl;
  }
  return 0;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    Solve(8);
  }
  else if (mode == 1) {
    for (int i = 1; i < 10; ++i)
    {
      cout << i << endl;
      Solve(i);
    }
  }

  return 0;
}
