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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U','L','D','R' };

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
}  // namespace

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

}  // namespace

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> sx >> sy >> tx >> ty >> forgetProb;
    rep(i, BOARD_SIZE)
    {
      rep(j, BOARD_SIZE-1)
      {
        char ccc;
        cin >> ccc;
        hWall[i][j] = ccc - '0';
      }
    }
    rep(i, BOARD_SIZE-1)
    {
      rep(j, BOARD_SIZE)
      {
        char ccc;
        cin >> ccc;
        vWall[i][j] = ccc - '0';
      }
    }
  }
  else {
    // ファイル入力
    ifs >> sx >> sy >> tx >> ty >> forgetProb;
    rep(i, BOARD_SIZE)
    {
      rep(j, BOARD_SIZE-1)
      {
        char ccc;
        ifs >> ccc;
        hWall[i][j] = ccc - '0';
      }
    }
    rep(i, BOARD_SIZE-1)
    {
      rep(j, BOARD_SIZE)
      {
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
    rep(i, min((int)route.size(), MAX_ROUTE_LEN))
    {
      cout << DIR_CHAR[route[i]];
    }
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    rep(i, min((int)route.size(), MAX_ROUTE_LEN))
    {
      ofs << DIR_CHAR[route[i]];
    }
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

// スコア計算
ll simulate_score(vector<int>& vec)
{
  int times = 100;
  double ret = 0;
  rep(_, times)
  {
    int le = vec.size();
    int x = sx, y = sy;
    rep(i, le)
    {
      if (Rand01() < forgetProb) {
        ;
      }
      else {
        int z = vec[i];
        if (z == 0) {
          if (x == 0 || vWall[x - 1][y]) continue;
        }
        if (z == 1) {
          if (y == 0 || hWall[x][y - 1]) continue;
        }
        if (z == 2) {
          if (x == 19 || vWall[x][y]) continue;
        }
        if (z == 3) {
          if (y == 19 || hWall[x][y]) continue;
        }

        x += dx[z];
        y += dy[z];
        if (x == tx && y == ty) {
          ret += 400 - i;
          break;
        }
      }
    }
  }
  return (ret / times) * 250000.0;
}

int Solve(int num)
{
  start_timer();

  input_data(num);

  route.clear();
  rep(i, MAX_ROUTE_LEN)
  {
    route.push_back(Rand() % 4);
  }
  cur_score = simulate_score(route);


  int dp[BOARD_SIZE][BOARD_SIZE];
  int prev_dir[21][21];
  int dist[BOARD_SIZE][BOARD_SIZE];
  rep(i, BOARD_SIZE)
  {
    rep(j, BOARD_SIZE)
    {
      dp[i][j] = INF;
      prev_dir[20][20] = -1;
      dist[i][j] = -1;
    }
  }

  dp[sx][sy] = 0;
  priority_queue<pair<int, P>, vector<pair<int, P>>, greater<pair<int, P>>> que;
  pair<int, P> tmp;
  tmp.first = 0;
  tmp.second.first = sx;
  tmp.second.second = sy;
  que.push(tmp);
  while (que.size()) {
    tmp = que.top();
    que.pop();
    int x = tmp.second.first;
    int y = tmp.second.second;
    if (tmp.first > dp[x][y]) continue;
    if (dp[x][y] > 250) continue;
    rep(i, 4)
    {
      int nx = x, ny = y;
      int cnt = 0;
      if (i == 0) {
        while (nx != 0 && vWall[nx - 1][ny] == 0) {
          nx--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 1) {
        while (ny != 0 && hWall[nx][ny - 1] == 0) {
          ny--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 2) {
        while (nx != 19 && vWall[nx][ny] == 0) {
          nx++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 3) {
        while (ny != 19 && hWall[nx][ny] == 0) {
          ny++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (nx == x && ny == y) continue;
      if (dp[nx][ny] > dp[x][y] + cnt + forgetProb * 15) {
        dp[nx][ny] = dp[x][y] + cnt + forgetProb * 15;
        dist[nx][ny] = cnt;
        prev_dir[nx][ny] = i;
        tmp.first = dp[nx][ny];
        tmp.second.first = nx;
        tmp.second.second = ny;
        que.push(tmp);
      }
    }
  }

  rep(i, MAX_ROUTE_LEN + 10)
  {
    route.push_back(rand() % 4);
  }


  if (dp[tx][ty] <= MAX_ROUTE_LEN) {
    route.clear();
    int x = tx, y = ty;
    while (x != sx || y != sy) {
      int nx = x - dx[prev_dir[x][y]] * dist[x][y];
      int ny = y - dy[prev_dir[x][y]] * dist[x][y];
      rep(i, (dist[x][y] * (1.0 / (1.0 - forgetProb))) + forgetProb * 15)
      {
        route.push_back(prev_dir[x][y]);
      }
      x = nx;
      y = ny;
    }

    reverse(route.begin(), route.end());

    while (route.size() < MAX_ROUTE_LEN) {
      route.push_back(rand() % 4);
    }

    cur_score = simulate_score(route);

  }
  int aaa = 0;
  best_score = -1;
  if (route.size() > MAX_ROUTE_LEN) {
    int ok = 0;
    rep(i, 10)
    {
      rep(j, 10)
      {
        route.clear();
        int x = tx - i, y = ty - j;
        if (dp[x][y] >= MAX_ROUTE_LEN) continue;
        while (x != sx || y != sy) {
          int nx = x - dx[prev_dir[x][y]] * dist[x][y];
          int ny = y - dy[prev_dir[x][y]] * dist[x][y];
          rep(i, (dist[x][y] * (1.0 / (1.0 - forgetProb))) + forgetProb * 15)
          {
            route.push_back(prev_dir[x][y]);
          }
          x = nx;
          y = ny;
        }

        reverse(route.begin(), route.end());

        while (route.size() < MAX_ROUTE_LEN) {
          route.push_back(rand() % 4);
        }
        if (route.size() == MAX_ROUTE_LEN) {
          ok = 1;
          aaa = 1;
          cur_score = simulate_score(route);
          if (cur_score > best_score) {
            best_score = cur_score;
            best_route = route;
          }
        }
      }
    }

    if (ok == 0) {
      route.clear();
      rep(i, MAX_ROUTE_LEN + 10)
      {
        route.push_back(rand() % 4);
      }
    }
    else {
      route = best_route;
      cur_score = best_score;
      if (cur_score == -1) {
        rep(i, MAX_ROUTE_LEN + 10)
        {
          route.push_back(rand() % 4);
        }
      }
    }
  }

  int loop = 0;
  if (aaa || route.size() > MAX_ROUTE_LEN) {
    // 愚直解
    if (aaa == 0) {
      route.clear();
      rep(i, MAX_ROUTE_LEN)
      {
        route.push_back(Rand() % 4);
      }
    }

    cur_score = simulate_score(route);

    best_route = route;
    best_score = cur_score;

    // 山登り解、焼きなまし解
    double now_time = get_elapsed_time();
    double TL = 1.9;
    double start_temp = 2048;
    double end_temp = 0.0001;
    int keep[1000][3];
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
      }
      if (now_time > TL) break;

      int x = Rand() % MAX_ROUTE_LEN;
      int y = Rand() % 4;
      int keepy = route[x];
      route[x] = y;

      int tmpScore = simulate_score(route);

      int diffScore = tmpScore - cur_score;

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      double prob = exp((double)diffScore / temp);
      if (prob > Rand01()) {
        cur_score += diffScore;
        if (cur_score > best_score) {
          best_score = cur_score;
          best_route = route;
        }
      }
      else {
        // 元に戻す
        route[x] = keepy;
      }
    }

    // 最高スコアを戻す
    route = best_route;
    cur_score = best_score;
  }

  output_data(num);

  // デバッグ用
  if (exec_mode != 0) {
    cout << "route.size() = " << route.size() << endl;
    cout << "loop = " << loop << endl;
    cout << cur_score << endl;
    cout << get_elapsed_time() << "sec." << endl;
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