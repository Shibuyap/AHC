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
const char cc[4] = { 'U','L','D','R' };

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
  const int n = 20;
  const int MAX_LENGTH = 200;
  int sx, sy, tx, ty;
  double pro;
  int h[n][n];
  int v[n][n];

  // 解答用変数
  ll maxScore;
  vector<int> ans;

  // 焼きなまし用変数
  ll real_maxScore;
  vector<int> real_ans;

}  // namespace

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> sx >> sy >> tx >> ty >> pro;
    rep(i, n)
    {
      rep(j, n-1)
      {
        char ccc;
        cin >> ccc;
        h[i][j] = ccc - '0';
      }
    }
    rep(i, n-1)
    {
      rep(j, n)
      {
        char ccc;
        cin >> ccc;
        v[i][j] = ccc - '0';
      }
    }
  }
  else {
    // ファイル入力
    ifs >> sx >> sy >> tx >> ty >> pro;
    rep(i, n)
    {
      rep(j, n-1)
      {
        char ccc;
        ifs >> ccc;
        h[i][j] = ccc - '0';
      }
    }
    rep(i, n-1)
    {
      rep(j, n)
      {
        char ccc;
        ifs >> ccc;
        v[i][j] = ccc - '0';
      }
    }
  }
}

void output_data(int case_num)
{
  if (exec_mode == 0) {
    // 標準出力
    rep(i, min((int)ans.size(), MAX_LENGTH))
    {
      cout << cc[ans[i]];
    }
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    rep(i, min((int)ans.size(), MAX_LENGTH))
    {
      ofs << cc[ans[i]];
    }
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

// スコア計算
ll CalcScore(vector<int>& vec)
{
  int times = 100;
  double ret = 0;
  rep(_, times)
  {
    int le = vec.size();
    int x = sx, y = sy;
    rep(i, le)
    {
      if (Rand01() < pro) {
        ;
      }
      else {
        int z = vec[i];
        if (z == 0) {
          if (x == 0 || v[x - 1][y]) continue;
        }
        if (z == 1) {
          if (y == 0 || h[x][y - 1]) continue;
        }
        if (z == 2) {
          if (x == 19 || v[x][y]) continue;
        }
        if (z == 3) {
          if (y == 19 || h[x][y]) continue;
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

  ans.clear();
  rep(i, MAX_LENGTH)
  {
    ans.push_back(Rand() % 4);
  }
  maxScore = CalcScore(ans);


  int dp[n][n];
  int dir[21][21];
  int dist[n][n];
  rep(i, n)
  {
    rep(j, n)
    {
      dp[i][j] = INF;
      dir[20][20] = -1;
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
        while (nx != 0 && v[nx - 1][ny] == 0) {
          nx--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 1) {
        while (ny != 0 && h[nx][ny - 1] == 0) {
          ny--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 2) {
        while (nx != 19 && v[nx][ny] == 0) {
          nx++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 3) {
        while (ny != 19 && h[nx][ny] == 0) {
          ny++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (nx == x && ny == y) continue;
      if (dp[nx][ny] > dp[x][y] + cnt + pro * 15) {
        dp[nx][ny] = dp[x][y] + cnt + pro * 15;
        dist[nx][ny] = cnt;
        dir[nx][ny] = i;
        tmp.first = dp[nx][ny];
        tmp.second.first = nx;
        tmp.second.second = ny;
        que.push(tmp);
      }
    }
  }

  rep(i, MAX_LENGTH + 10)
  {
    ans.push_back(rand() % 4);
  }


  if (dp[tx][ty] <= MAX_LENGTH) {
    ans.clear();
    int x = tx, y = ty;
    while (x != sx || y != sy) {
      int nx = x - dx[dir[x][y]] * dist[x][y];
      int ny = y - dy[dir[x][y]] * dist[x][y];
      rep(i, (dist[x][y] * (1.0 / (1.0 - pro))) + pro * 15)
      {
        ans.push_back(dir[x][y]);
      }
      x = nx;
      y = ny;
    }

    reverse(ans.begin(), ans.end());

    while (ans.size() < MAX_LENGTH) {
      ans.push_back(rand() % 4);
    }

    maxScore = CalcScore(ans);

  }
  int aaa = 0;
  real_maxScore = -1;
  if (ans.size() > MAX_LENGTH) {
    int ok = 0;
    rep(i, 10)
    {
      rep(j, 10)
      {
        ans.clear();
        int x = tx - i, y = ty - j;
        if (dp[x][y] >= MAX_LENGTH) continue;
        while (x != sx || y != sy) {
          int nx = x - dx[dir[x][y]] * dist[x][y];
          int ny = y - dy[dir[x][y]] * dist[x][y];
          rep(i, (dist[x][y] * (1.0 / (1.0 - pro))) + pro * 15)
          {
            ans.push_back(dir[x][y]);
          }
          x = nx;
          y = ny;
        }

        reverse(ans.begin(), ans.end());

        while (ans.size() < MAX_LENGTH) {
          ans.push_back(rand() % 4);
        }
        if (ans.size() == MAX_LENGTH) {
          ok = 1;
          aaa = 1;
          maxScore = CalcScore(ans);
          if (maxScore > real_maxScore) {
            real_maxScore = maxScore;
            real_ans = ans;
          }
        }
      }
    }

    if (ok == 0) {
      ans.clear();
      rep(i, MAX_LENGTH + 10)
      {
        ans.push_back(rand() % 4);
      }
    }
    else {
      ans = real_ans;
      maxScore = real_maxScore;
      if (maxScore == -1) {
        rep(i, MAX_LENGTH + 10)
        {
          ans.push_back(rand() % 4);
        }
      }
    }
  }

  int loop = 0;
  if (aaa || ans.size() > MAX_LENGTH) {
    // 愚直解
    if (aaa == 0) {
      ans.clear();
      rep(i, MAX_LENGTH)
      {
        ans.push_back(Rand() % 4);
      }
    }

    maxScore = CalcScore(ans);

    real_ans = ans;
    real_maxScore = maxScore;

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

      int x = Rand() % MAX_LENGTH;
      int y = Rand() % 4;
      int keepy = ans[x];
      ans[x] = y;

      int tmpScore = CalcScore(ans);

      int diffScore = tmpScore - maxScore;

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      double prob = exp((double)diffScore / temp);
      if (prob > Rand01()) {
        maxScore += diffScore;
        if (maxScore > real_maxScore) {
          real_maxScore = maxScore;
          real_ans = ans;
        }
      }
      else {
        // 元に戻す
        ans[x] = keepy;
      }
    }

    // 最高スコアを戻す
    ans = real_ans;
    maxScore = real_maxScore;
  }

  output_data(num);

  // デバッグ用
  if (exec_mode != 0) {
    cout << "ans.size() = " << ans.size() << endl;
    cout << "loop = " << loop << endl;
    cout << maxScore << endl;
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