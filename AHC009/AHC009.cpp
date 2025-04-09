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

  // 0以上1未満の小数をとる乱数
  static double rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* 変数 */
{
  // 入力用変数
  int n = 20;
  int sx, sy, tx, ty;
  double pro;
  int h[20][20];
  int v[20][20];

  // 解答用変数
  ll maxScore;
  vector<int> ans;

  // 焼きなまし用変数
  ll real_maxScore;
  vector<int> real_ans;

}  // namespace

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
      if (rand01() < pro) {
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

///////////////////////////////////////////////////////////////////////////////
// 始まったらやること
// 1. 入力部
// 2. maxScoreとans
// 3. 出力部
// 4. 愚直解
// 5. スコア計算関数
// 6. 貪欲解
// 7. 山登り
// 8. 焼きなまし
///////////////////////////////////////////////////////////////////////////////
int Solve(int mode, int num)
{
  srand((unsigned)time(NULL));
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();
  while (rand() % 100) {
    Rand();
  }

  // 入力部
  string fileNameIfs = "./in/000" + to_string(num) + ".txt";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    cin >> sx >> sy >> tx >> ty >> pro;
    rep(i, 20)
    {
      rep(j, 19)
      {
        char ccc;
        cin >> ccc;
        h[i][j] = ccc - '0';
      }
    }
    rep(i, 19)
    {
      rep(j, 20)
      {
        char ccc;
        cin >> ccc;
        v[i][j] = ccc - '0';
      }
    }
  }
  else {  // ファイル入力する
    mode = 1;
    ifs >> sx >> sy >> tx >> ty >> pro;
    rep(i, 20)
    {
      rep(j, 19)
      {
        char ccc;
        ifs >> ccc;
        h[i][j] = ccc - '0';
      }
    }
    rep(i, 19)
    {
      rep(j, 20)
      {
        char ccc;
        ifs >> ccc;
        v[i][j] = ccc - '0';
      }
    }
  }

  ans.clear();
  rep(i, 200)
  {
    ans.push_back(Rand() % 4);
  }
  maxScore = CalcScore(ans);


  int dp[20][20];
  int dir[21][21];
  int dist[20][20];
  rep(i, 20)
  {
    rep(j, 20)
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
    //cout << "     " << que.size() << endl;
    tmp = que.top();
    que.pop();
    int x = tmp.second.first;
    int y = tmp.second.second;
    //cout << x << y << endl;
    if (tmp.first > dp[x][y]) continue;
    if (dp[x][y] > 250) continue;
    rep(i, 4)
    {
      //cout << i;
      int nx = x, ny = y;
      int cnt = 0;
      if (i == 0) {
        while (nx != 0 && v[nx - 1][ny] == 0) {
          //cout << i;
          //cout << nx << ' ' << ny << endl;
          nx--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 1) {
        while (ny != 0 && h[nx][ny - 1] == 0) {
          //cout << i;
          ny--;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 2) {
        while (nx != 19 && v[nx][ny] == 0) {
          //cout << i;
          nx++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (i == 3) {
        while (ny != 19 && h[nx][ny] == 0) {
          //cout << i;
          ny++;
          cnt++;
          if (nx == tx && ny == ty) {
            break;
          }
        }
      }
      if (nx == x && ny == y) continue;
      //cout << "OK";
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

  rep(i, 210)
  {
    ans.push_back(rand() % 4);
  }


  if (dp[tx][ty] <= 200) {
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

    while (ans.size() < 200) {
      ans.push_back(rand() % 4);
    }

    maxScore = CalcScore(ans);

  }
  int aaa = 0;
  real_maxScore = -1;
  if (ans.size() > 200) {
    int ok = 0;
    rep(i, 10)
    {
      rep(j, 10)
      {
        ans.clear();
        int x = tx - i, y = ty - j;
        if (dp[x][y] >= 200) continue;
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

        while (ans.size() < 200) {
          ans.push_back(rand() % 4);
        }
        if (ans.size() == 200) {
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
      rep(i, 210)
      {
        ans.push_back(rand() % 4);
      }
    }
    else {
      ans = real_ans;
      maxScore = real_maxScore;
      if (maxScore == -1) {
        rep(i, 210)
        {
          ans.push_back(rand() % 4);
        }
      }
    }
  }

  int loop = 0;
  if (aaa || ans.size() > 200) {
    // 愚直解
    if (aaa == 0) {
      ans.clear();
      rep(i, 200)
      {
        ans.push_back(Rand() % 4);
      }
    }

    maxScore = CalcScore(ans);

    real_ans = ans;
    real_maxScore = maxScore;

    // 山登り解、焼きなまし解
    double now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double TL = 1.9;
    double start_temp = 2048;
    double end_temp = 0.0001;
    // int loop          = 0;
    int keep[1000][3];
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        end_time = clock();
        now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
      }
      if (now_time > TL) break;

      int x = Rand() % 200;
      int y = Rand() % 4;
      int keepy = ans[x];
      ans[x] = y;

      int tmpScore = CalcScore(ans);

      int diffScore = tmpScore - maxScore;

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      double prob = exp((double)diffScore / temp);
      // if (tmp > 0) {
      if (prob > rand01()) {
        // cout << tmpScore << endl;
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

  // 解の出力
  if (mode == 0) {
    rep(i, min((int)ans.size(), 200))
    {
      cout << cc[ans[i]];
    }
    cout << endl;
  }

  // デバッグ用
  if (mode != 0) {
    cout << "ans.size() = " << ans.size() << endl;
    cout << "loop = " << loop << endl;
    cout << maxScore << endl;
    end_time = clock();
    cout << (double)(end_time - start_time) / CLOCKS_PER_SEC << "sec." << endl;
  }

  // ファイル出力
  if (mode != 0) {
    string fileNameOfs = "sample_out.txt";
    ofstream ofs(fileNameOfs);
    rep(i, min((int)ans.size(), 200))
    {
      ofs << cc[ans[i]];
    }
    ofs << endl;
    ofs.close();
  }

  return 0;
}

int main()
{
  int mode = 0;

  if (mode == 0) {
    Solve(mode, 8);
  }
  else if (mode == 1) {

    srep(i, 1, 10)
    {
      cout << i << endl;
      Solve(mode, i);
    }
  }

  return 0;
}