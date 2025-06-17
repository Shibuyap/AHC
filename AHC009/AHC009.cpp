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

#define srep(i, s, t) for (int i = s; i < t; ++i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U', 'L', 'D', 'R' };

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

int exec_mode;

namespace /* 変数 */
{
  // 入力用変数
  const int BOARD_SIZE = 20;
  const int MAX_ROUTE_LEN = 200;
  int sx, sy, tx, ty;
  double forgetProb;
  int hWall[BOARD_SIZE][BOARD_SIZE];
  int vWall[BOARD_SIZE][BOARD_SIZE];

  // 解答用変数
  ll cur_score;
  vector<int> route;

  // 焼きなまし用変数
  ll best_score;
  vector<int> best_route;

} // namespace

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> sx >> sy >> tx >> ty >> forgetProb;
    for (int i = 0; i < BOARD_SIZE; ++i) {
      for (int j = 0; j < BOARD_SIZE - 1; ++j) {
        char ccc;
        cin >> ccc;
        hWall[i][j] = ccc - '0';
      }
    }
    for (int i = 0; i < BOARD_SIZE - 1; ++i) {
      for (int j = 0; j < BOARD_SIZE; ++j) {
        char ccc;
        cin >> ccc;
        vWall[i][j] = ccc - '0';
      }
    }
  }
  else {
    // ファイル入力
    ifs >> sx >> sy >> tx >> ty >> forgetProb;
    for (int i = 0; i < BOARD_SIZE; ++i) {
      for (int j = 0; j < BOARD_SIZE - 1; ++j) {
        char ccc;
        ifs >> ccc;
        hWall[i][j] = ccc - '0';
      }
    }
    for (int i = 0; i < BOARD_SIZE - 1; ++i) {
      for (int j = 0; j < BOARD_SIZE; ++j) {
        char ccc;
        ifs >> ccc;
        vWall[i][j] = ccc - '0';
      }
    }
  }
}

void output_data(int case_num)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < static_cast<int>(min(route.size(), static_cast<size_t>(MAX_ROUTE_LEN))); ++i) {
      cout << DIR_CHAR[route[i]];
    }
    cout << '¥n';
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < static_cast<int>(min(route.size(), static_cast<size_t>(MAX_ROUTE_LEN))); ++i) {
      ofs << DIR_CHAR[route[i]];
    }
    ofs << '¥n';
  }
}

// スコア計算
ll simulate_score(vector<int>& vec)
{
  const int TRIALS = 100;
  double totalReward = 0.0;

  for (int trial = 0; trial < TRIALS; ++trial) {
    int len = vec.size();
    int row = sx, col = sy;

    for (int step = 0; step < len; ++step) {
      if (Rand01() < forgetProb) {
        // 忘却：その場に留まる
      }
      else {
        int dir = vec[step];

        if (dir == 0) {
          if (row == 0 || vWall[row - 1][col])
            continue;
        }
        if (dir == 1) {
          if (col == 0 || hWall[row][col - 1])
            continue;
        }
        if (dir == 2) {
          if (row == 19 || vWall[row][col])
            continue;
        }
        if (dir == 3) {
          if (col == 19 || hWall[row][col])
            continue;
        }

        row += dx[dir];
        col += dy[dir];
        if (row == tx && col == ty) {
          totalReward += 400 - step;
          break;
        }
      }
    }
  }
  return (totalReward / TRIALS) * 250000.0;
}

int Solve(int caseId)
{
  start_timer();

  input_data(caseId);

  route.clear();
  for (int i = 0; i < MAX_ROUTE_LEN; ++i) {
    route.push_back(Rand() % 4);
  }
  cur_score = simulate_score(route);

  /* ──────────────────────  Dijkstra 前処理  ────────────────────── */
  int cost[BOARD_SIZE][BOARD_SIZE];
  int prevDir[21][21];
  int segLen[BOARD_SIZE][BOARD_SIZE];
  for (int r = 0; r < BOARD_SIZE; ++r) {
    for (int c = 0; c < BOARD_SIZE; ++c) {
      cost[r][c] = INF;
      prevDir[20][20] = -1;
      segLen[r][c] = -1;
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
      int nextRow = row, nextCol = col;
      int steps = 0;

      if (d == 0) {
        while (nextRow != 0 && vWall[nextRow - 1][nextCol] == 0) {
          --nextRow;
          ++steps;
          if (nextRow == tx && nextCol == ty)
            break;
        }
      }
      if (d == 1) {
        while (nextCol != 0 && hWall[nextRow][nextCol - 1] == 0) {
          --nextCol;
          ++steps;
          if (nextRow == tx && nextCol == ty)
            break;
        }
      }
      if (d == 2) {
        while (nextRow != 19 && vWall[nextRow][nextCol] == 0) {
          ++nextRow;
          ++steps;
          if (nextRow == tx && nextCol == ty)
            break;
        }
      }
      if (d == 3) {
        while (nextCol != 19 && hWall[nextRow][nextCol] == 0) {
          ++nextCol;
          ++steps;
          if (nextRow == tx && nextCol == ty)
            break;
        }
      }
      if (nextRow == row && nextCol == col)
        continue;

      int newCost = cost[row][col] + steps + forgetProb * 15;
      if (cost[nextRow][nextCol] > newCost) {
        cost[nextRow][nextCol] = newCost;
        segLen[nextRow][nextCol] = steps;
        prevDir[nextRow][nextCol] = d;
        pq.push({ newCost, {nextRow, nextCol} });
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
      int d = prevDir[row][col];
      int steps = segLen[row][col];
      int nextRow = row - dx[d] * steps;
      int nextCol = col - dy[d] * steps;

      for (int k = 0; k < static_cast<int>(steps / (1.0 - forgetProb) + forgetProb * 15); ++k)
        route.push_back(d);

      row = nextRow;
      col = nextCol;
    }
    reverse(route.begin(), route.end());
    while (route.size() < MAX_ROUTE_LEN)
      route.push_back(rand() % 4);
    cur_score = simulate_score(route);
  }

  bool needRandomInit = false;
  best_score = -1;
  if (route.size() > MAX_ROUTE_LEN) {
    bool foundFit = false;
    for (int di = 0; di < 10; ++di) {
      for (int dj = 0; dj < 10; ++dj) {
        route.clear();
        int row = tx - di, col = ty - dj;
        if (cost[row][col] >= MAX_ROUTE_LEN) {
          continue;
        }

        while (row != sx || col != sy) {
          int d = prevDir[row][col];
          int steps = segLen[row][col];
          int nextRow = row - dx[d] * steps;
          int nextCol = col - dy[d] * steps;

          for (int k = 0; k < static_cast<int>(steps / (1.0 - forgetProb) + forgetProb * 15); ++k) {
            route.push_back(d);
          }

          row = nextRow;
          col = nextCol;
        }
        reverse(route.begin(), route.end());
        while (route.size() < MAX_ROUTE_LEN) {
          route.push_back(rand() % 4);
        }

        if (route.size() == MAX_ROUTE_LEN) {
          foundFit = true;
          needRandomInit = true;
          cur_score = simulate_score(route);
          if (cur_score > best_score) {
            best_score = cur_score;
            best_route = route;
          }
        }
      }
    }
    if (!foundFit) {
      route.clear();
      for (int i = 0; i < MAX_ROUTE_LEN + 10; ++i) route.push_back(rand() % 4);
    }
    else {
      route = best_route;
      cur_score = best_score;
      if (cur_score == -1)
        for (int i = 0; i < MAX_ROUTE_LEN + 10; ++i) route.push_back(rand() % 4);
    }
  }

  /* ──────────────────────  焼きなまし  ────────────────────── */
  int iter = 0;
  if (needRandomInit || route.size() > MAX_ROUTE_LEN) {
    if (!needRandomInit) {
      route.clear();
      for (int i = 0; i < MAX_ROUTE_LEN; ++i) route.push_back(Rand() % 4);
    }
    cur_score = simulate_score(route);
    best_route = route;
    best_score = cur_score;

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
      int diff = candScore - cur_score;

      double temp = tempInit + (tempFinal - tempInit) * elapsed / timeLimit;
      double prob = exp(diff / temp);

      if (diff > 0 || prob > Rand01()) {
        cur_score = candScore;
        if (cur_score > best_score) {
          best_score = cur_score;
          best_route = route;
        }
      }
      else {
        route[pos] = backup;
      }
    }
    route = best_route;
    cur_score = best_score;
  }

  output_data(caseId);

  if (exec_mode != 0) {
    cout << "route.size() = " << route.size() << endl;
    cout << "iter = " << iter << endl;
    cout << cur_score << endl;
    cout << get_elapsed_time() << " sec." << endl;
  }
  return 0;
}

int main()
{
  exec_mode = 1;

  if (exec_mode == 0) {
    Solve(8);
  }
  else if (exec_mode == 1) {
    srep(i, 1, 10)
    {
      cout << i << endl;
      Solve(i);
    }
  }

  return 0;
}
