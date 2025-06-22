#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<double, P> PDP;

#define MAX_N 200005
const int INF = 1001001001;

namespace /* 乱数 */
{
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
}

namespace /* グリッド用 */
{
  inline int ManhattanDistance(int x1, int y1, int x2, int y2)
  {
    return std::abs(x1 - x2) + std::abs(y1 - y2);
  }

  const int dx[4] = { -1, 0, 1, 0 };
  const int dy[4] = { 0, -1, 0, 1 };
  const string ULDR = "ULDR";
}

namespace /* 問題設定 */
{
  const int N = 30;
  const int Q = 1000;
}

namespace /* ハイパーパラメータ */
{
  double initialD = 3585;
  double start_temp = 103.48026;
  double end_temp = 6.74495e-07;
  int loopTimes = 3007;
  int addTimes = 5007;

  double SabunCostMultiple = 1.5;
}

namespace /* ローカル用 */
{
  double dReal[2][N + 1][N + 1]; // 正解の長さ
  double scoreSumGlobal; // スコア管理
}

namespace /* 変数 */
{
  PDP p;
}

namespace /* 焼きなまし用変数 */
{
  double dist_res[Q] = {};				// 受け取った長さ
}

namespace /* 行、列、切れ目の構造 */
{
  double up[N] = {}, down[N] = {}, l[N] = {}, r[N] = {};
  int cut_v[N] = {}, cut_h[N] = {}; // 0‾30をとる半開区間
  int vsum[Q][N + 1][N + 1];
  int hsum[Q][N + 1][N + 1];
  vector<int> turn_v[N], turn_h[N];

  double best_up[N], best_down[N], best_l[N], best_r[N];
  int best_cut_v[N], best_cut_h[N];

  double diff_sum = 0;
  double dist_est[Q] = {};

  vector<pair<double, double>> vec(1100);

  double dUD[N + 1][N + 1], dLR[N + 1][N + 1];
  vector<int> PathIDVectorUD[N][N], PathIDVectorLR[N][N];
}

void ClearGlobalVariables()
{
  scoreSumGlobal = 0;
  for (int i = 0; i < 2; ++i)for (int j = 0; j < N + 1; ++j)for (int k = 0; k < N + 1; ++k) dReal[i][j][k] = 0;
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < N + 1; ++j) {
      for (int k = 0; k < N + 1; ++k) {
        vsum[i][j][k] = 0;
        hsum[i][j][k] = 0;
      }
    }
  }
  for (int i = 0; i < N; ++i) {
    turn_v[i].clear();
    turn_h[i].clear();
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      PathIDVectorLR[i][j].clear();
      PathIDVectorUD[i][j].clear();
      dUD[i][j] = 0;
      dLR[i][j] = 0;
    }
  }
}

bool IsOutOfBounds(int nx, int ny)
{
  if (nx < 0 || N <= nx || ny < 0 || N <= ny) return 1;
  return 0;
}

int CalcScore(const int sx, const int sy, const string& ans, const double aValue, const double eValue, int k)
{
  int x = sx, y = sy;
  int res = 0;
  for (int i = 0; i < ans.size(); ++i) {
    if (ans[i] == 'U') {
      res += dReal[0][x][y];
      x--;
    }
    if (ans[i] == 'L') {
      res += dReal[1][x][y];
      y--;
    }
    if (ans[i] == 'D') {
      res += dReal[0][x + 1][y];
      x++;
    }
    if (ans[i] == 'R') {
      res += dReal[1][x][y + 1];
      y++;
    }
  }
  scoreSumGlobal += pow(0.998, 1000 - k) * aValue / res;
  return res * eValue;
}

double dp[N][N];
int nxt[N][N];

void AnnealingMode0(double temp);
void AnnealingMode1(double temp);
void AnnealingMode2(double temp);
void AnnealingMode3();
void AnnealingMode3Vertical(int rN);
void AnnealingMode3Horizontal(int rN);
void FinalAdjustment();

void Dijkstra(int sx, int sy, int gx, int gy)
{
  for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) dp[i][j] = INF;
  dp[sx][sy] = 0;
  priority_queue<PDP, vector<PDP>, greater<PDP>> que;
  p.first = 0;
  p.second.first = sx;
  p.second.second = sy;
  que.push(p);
  while (que.size()) {
    p = que.top();
    que.pop();
    int x = p.second.first;
    int y = p.second.second;
    double val = p.first;
    if (x == gx && y == gy) { break; }
    if (val > dp[x][y]) { continue; }
    for (int i = 0; i < 4; ++i) {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (IsOutOfBounds(nx, ny)) { continue; }
      double dd = 0;
      if (i == 0) {
        if (x <= cut_v[y]) dd = up[y];
        else dd = down[y];
      }
      if (i == 1) {
        if (y <= cut_h[x]) dd = l[x];
        else dd = r[x];
      }
      if (i == 2) {
        if (x + 1 <= cut_v[y]) dd = up[y];
        else dd = down[y];
      }
      if (i == 3) {
        if (y + 1 <= cut_h[x]) dd = l[x];
        else dd = r[x];
      }
      if (dp[nx][ny] > dp[x][y] + dd) {
        dp[nx][ny] = dp[x][y] + dd;
        nxt[nx][ny] = i;
        p.first = dp[nx][ny];
        p.second.first = nx;
        p.second.second = ny;
        que.push(p);
      }
    }
  }
}

void FinalAdjustment()
{
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      dUD[i][j] = 0;
      dLR[i][j] = 0;
    }
  }

  for (int iter = 0; iter < 20000; ++iter) {
    int dir = Rand() % 2;
    int x = 0, y = 0;
    if (dir == 0) {
      x = Rand() % 29 + 1;
      y = Rand() % 30;
    }
    else {
      x = Rand() % 30;
      y = Rand() % 29 + 1;
    }
    double delta = Rand01() * 200.0 - 100.0;
    if (dir == 0) {
      if (x <= cut_v[y]) {
        if (up[y] + delta < 1000 || 9000 < up[y] + delta) { continue; }
      }
      else {
        if (down[y] + delta < 1000 || 9000 < down[y] + delta) { continue; }
      }
    }
    else {
      if (y <= cut_h[x]) {
        if (l[x] + delta < 1000 || 9000 < l[x] + delta) { continue; }
      }
      else {
        if (r[x] + delta < 1000 || 9000 < r[x] + delta) { continue; }
      }
    }

    double diff = 0;

    if (dir == 0) {
      for (int i = 0; i < PathIDVectorUD[x][y].size(); ++i) {
        int t = PathIDVectorUD[x][y][i];
        diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + delta))) * (40000.0 / dist_res[t]);
      }
    }
    else {
      for (int i = 0; i < PathIDVectorLR[x][y].size(); ++i) {
        int t = PathIDVectorLR[x][y][i];
        diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + delta))) * (40000.0 / dist_res[t]);
      }
    }

    if (diff > 0) {
      if (dir == 0) {
        for (int i = 0; i < PathIDVectorUD[x][y].size(); ++i) {
          int t = PathIDVectorUD[x][y][i];
          dist_est[t] += delta;
        }
      }
      else {
        for (int i = 0; i < PathIDVectorLR[x][y].size(); ++i) {
          int t = PathIDVectorLR[x][y][i];
          dist_est[t] += delta;
        }
      }
      if (dir == 0) dUD[x][y] += delta;
      else dLR[x][y] += delta;
      diff_sum -= diff;
    }
  }
}

void Dijkstra2(int sx, int sy, int gx, int gy)
{
  for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) dp[i][j] = INF;
  dp[sx][sy] = 0;
  priority_queue<PDP, vector<PDP>, greater<PDP>> que;
  p.first = 0;
  p.second.first = sx;
  p.second.second = sy;
  que.push(p);
  while (que.size()) {
    p = que.top();
    que.pop();
    int x = p.second.first;
    int y = p.second.second;
    double val = p.first;
    if (x == gx && y == gy) { break; }
    if (val > dp[x][y]) { continue; }
    for (int i = 0; i < 4; ++i) {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (IsOutOfBounds(nx, ny)) { continue; }
      double dd = 0;
      if (i == 0) {
        if (x <= cut_v[y]) dd = up[y];
        else dd = down[y];
        dd += dUD[x][y];
      }
      if (i == 1) {
        if (y <= cut_h[x]) dd = l[x];
        else dd = r[x];
        dd += dLR[x][y];
      }
      if (i == 2) {
        if (x + 1 <= cut_v[y]) dd = up[y];
        else dd = down[y];
        dd += dUD[x + 1][y];
      }
      if (i == 3) {
        if (y + 1 <= cut_h[x]) dd = l[x];
        else dd = r[x];
        dd += dLR[x][y + 1];
      }

      dd = max(dd, 1.0);

      if (dp[nx][ny] > dp[x][y] + dd) {
        dp[nx][ny] = dp[x][y] + dd;
        nxt[nx][ny] = i;
        p.first = dp[nx][ny];
        p.second.first = nx;
        p.second.second = ny;
        que.push(p);
      }
    }
  }
}


/* メモ
  辺の番号は1‾29(1-index)
  縦辺(i-1,j)->(i,j)はd[0][i][j];
  横辺(i,j-1)->(i,j)はd[1][i][j];
*/

string ConstructPath(int sx, int sy, int gx, int gy, vector<int>& v)
{
  string ans;
  while (gx != sx || gy != sy) {
    v.push_back(nxt[gx][gy]);
    ans += ULDR[nxt[gx][gy]];
    int nx = gx - dx[nxt[gx][gy]];
    int ny = gy - dy[nxt[gx][gy]];
    gx = nx;
    gy = ny;
  }
  reverse(ans.begin(), ans.end());
  reverse(v.begin(), v.end());
  return ans;
}

void UpdatePathInfo(int turn, int sx, int sy, const vector<int>& v)
{
  int vpath[N + 1][N + 1];
  int hpath[N + 1][N + 1];
  for (int i = 0; i < N + 1; ++i) {
    for (int j = 0; j < N + 1; ++j) {
      vpath[i][j] = 0;
      hpath[i][j] = 0;
    }
  }

  for (int i = 0; i < v.size(); ++i) {
    int x = sx;
    int y = sy;
    if (v[i] == 2) x += 1;
    if (v[i] == 3) y += 1;

    if (v[i] % 2 == 0) {
      vpath[x][y] = 1;
      PathIDVectorUD[x][y].push_back(turn);
    }
    else {
      hpath[x][y] = 1;
      PathIDVectorLR[x][y].push_back(turn);
    }

    sx += dx[v[i]];
    sy += dy[v[i]];
  }

  for (int j = 0; j < N; ++j) {
    for (int i = 1; i < N; ++i) {
      vsum[turn][i][j] = vsum[turn][i - 1][j] + vpath[i][j];
    }
    if (vsum[turn][N - 1][j]) turn_v[j].push_back(turn);
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 1; j < N; ++j) {
      hsum[turn][i][j] = hsum[turn][i][j - 1] + hpath[i][j];
    }
    if (hsum[turn][i][N - 1]) turn_h[i].push_back(turn);
  }
}

void UpdateDiffSum(int turn)
{
  diff_sum = 0;
  for (int i = 0; i < Q; ++i) dist_est[i] = 0;

  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < turn_v[j].size(); ++i) {
      int t = turn_v[j][i];
      dist_est[t] += up[j] * vsum[t][cut_v[j]][j] + down[j] * (vsum[t][N - 1][j] - vsum[t][cut_v[j]][j]);
    }
  }
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < turn_h[j].size(); ++i) {
      int t = turn_h[j][i];
      dist_est[t] += l[j] * hsum[t][j][cut_h[j]] + r[j] * (hsum[t][j][N - 1] - hsum[t][j][cut_h[j]]);
    }
  }
  for (int i = 0; i < turn + 1; ++i) {
    diff_sum += std::abs(dist_res[i] - dist_est[i]) * (40000.0 / dist_res[i]);
  }

  for (int i = 0; i < N; ++i) {
    diff_sum += std::abs(up[i] - down[i]) * SabunCostMultiple;
    diff_sum += std::abs(l[i] - r[i]) * SabunCostMultiple;
  }
}

void SaveBestParams()
{
  for (int i = 0; i < N; ++i) {
    best_up[i] = up[i];
    best_down[i] = down[i];
    best_l[i] = l[i];
    best_r[i] = r[i];
    best_cut_v[i] = cut_v[i];
    best_cut_h[i] = cut_h[i];
  }
}

void RestoreBestParams()
{
  for (int i = 0; i < N; ++i) {
    up[i] = best_up[i];
    down[i] = best_down[i];
    l[i] = best_l[i];
    r[i] = best_r[i];
    cut_v[i] = best_cut_v[i];
    cut_h[i] = best_cut_h[i];
  }
}

int Solve(string inputFileNum)
{
  // 時間計測
  clock_t start, end;
  start = clock();

  int inputMode = 1;
  string fileNameIfs = (string)"./in/" + inputFileNum + ".txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs) { // 標準入力する
    inputMode = 0;
  }
  else {
    for (int i = 0; i < N; ++i) for (int j = 1; j < N; ++j) ifs >> dReal[1][i][j];
    for (int i = 1; i < N; ++i) for (int j = 0; j < N; ++j) ifs >> dReal[0][i][j];
  }

  string fileName = (string)"./out/" + inputFileNum + "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);

  for (int i = 0; i < N; ++i) {
    up[i] = initialD;
    down[i] = initialD;
    l[i] = initialD;
    r[i] = initialD;
    cut_v[i] = 15;
    cut_h[i] = 15;
  }
  int time_over = 0;
  int all_loop = 0;
  for (int turn = 0; turn < Q; ++turn) {
    int sx, sy, gx, gy;

    // 入力受け取り
    if (inputMode == 0) {
      cin >> sx >> sy >> gx >> gy;
    }
    else {
      ifs >> sx >> sy >> gx >> gy;
    }

    // キープ
    int ksx = sx, ksy = sy, kgx = gx, kgy = gy;

    // プライオリティーキューで最短路を求める（ダイクストラ）
    Dijkstra2(sx, sy, gx, gy);

    // 回答文字列生成
    vector<int> v;
    string ans = ConstructPath(sx, sy, gx, gy, v);
    gx = kgx; gy = kgy;
    if (inputMode == 0) {
      cout << ans << endl;
      fflush(stdout);
    }
    else {
      ofs << ans << endl;
    }

    // パスの長さの受け取り
    double dist = 0;
    if (inputMode == 0) {
      cin >> dist;
    }
    else {
      double aValue, eValue;
      ifs >> aValue >> eValue;
      dist = CalcScore(sx, sy, ans, aValue, eValue, turn + 1);
    }

    // 山登り用
    dist_res[turn] = dist;

    int length = ans.size();

    if (turn < 0) {
      continue;
    }
    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 1.8) {
      if (time_over == 0 && ifs) {
        cout << "TIME OVER" << ' ' << turn << endl;
        time_over = 1;
      }
      continue;
    }

    // 新しいデータ構造更新
    UpdatePathInfo(turn, sx, sy, v);
    sx = ksx; sy = ksy;

    UpdateDiffSum(turn);

    double best_diff_sum = diff_sum;
    SaveBestParams();

    // 焼きなましで予測値をさらに調整
    end = clock();
    double now_time = (double)(end - start) / CLOCKS_PER_SEC;
    double TL = 1.8;

    for (int loop = 0; loop < loopTimes; ++loop) {
      all_loop++;
      if (loop % 100 == 1) {
        end = clock();
        now_time = (double)(end - start) / CLOCKS_PER_SEC;
      }
      if (now_time > TL) {
        break;
      }

      double temp = start_temp + (end_temp - start_temp) * ((double)loop / (loopTimes + addTimes));


      int MODE = loop % 3;
      if (turn >= 200 && all_loop % 5247 == 1853) MODE = 3;

      // up, down, left, rightのいずれか一つを変更
      if (MODE == 0) {
        AnnealingMode0(temp);
      }

      // upとdown、もしくはleftとrightを変更
      if (MODE == 1) {
        AnnealingMode1(temp);
      }

      // 区切りを一か所変更
      if (MODE == 2) {
        AnnealingMode2(temp);
      }

      // 1ライン調整
      if (MODE == 3) {
        AnnealingMode3();
      }
    }


    // スコアが悪化したらロールバック
    if (best_diff_sum < diff_sum) {
      diff_sum = best_diff_sum;
      RestoreBestParams();
    }

    FinalAdjustment();
  }

  scoreSumGlobal *= 2312311;

  if (inputMode == 1) {
    string fileNameParam = (string)"./parameter/" + inputFileNum + "_param.txt";
    const char* cstrParam = fileNameParam.c_str();
    ofstream ofsParam(cstrParam);
    ofsParam << "up / down" << endl;
    for (int i = 0; i < N; ++i) ofsParam << i << ' ' << cut_v[i] << ' ' << up[i] << ' ' << down[i] << endl;
    ofsParam << "left / right" << endl;
    for (int i = 0; i < N; ++i) ofsParam << i << ' ' << cut_h[i] << ' ' << l[i] << ' ' << r[i] << endl;

    ofsParam << "dUD" << endl;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        ofsParam << dUD[i][j] << ' ';
      }
      ofsParam << endl;
    }
    ofsParam << "dLR" << endl;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        ofsParam << dLR[i][j] << ' ';
      }
      ofsParam << endl;
    }

    ofsParam.close();
  }

  ofs.close();

  return scoreSumGlobal;
}

void AnnealingMode0(double temp)
{
  int dir = Rand() % 4;
  int idx = Rand() % N;
  double delta = Rand01() * 20 - 10;

  if (dir == 0) if (up[idx] + delta < 1000 || 9000 < up[idx] + delta) { return; }
  if (dir == 1) if (l[idx] + delta < 1000 || 9000 < l[idx] + delta) { return; }
  if (dir == 2) if (down[idx] + delta < 1000 || 9000 < down[idx] + delta) { return; }
  if (dir == 3) if (r[idx] + delta < 1000 || 9000 < r[idx] + delta) { return; }

  double diff = 0;

  if (dir == 0) {
    for (int i = 0; i < turn_v[idx].size(); ++i) {
      int t = turn_v[idx][i];
      double d = delta * vsum[t][cut_v[idx]][idx];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
    diff += std::abs(up[idx] - down[idx]) * SabunCostMultiple - std::abs((up[idx] + delta) - down[idx]) * SabunCostMultiple;
  }
  if (dir == 1) {
    for (int i = 0; i < turn_h[idx].size(); ++i) {
      int t = turn_h[idx][i];
      double d = delta * hsum[t][idx][cut_h[idx]];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
    diff += std::abs(l[idx] - r[idx]) * SabunCostMultiple - std::abs((l[idx] + delta) - r[idx]) * SabunCostMultiple;
  }
  if (dir == 2) {
    for (int i = 0; i < turn_v[idx].size(); ++i) {
      int t = turn_v[idx][i];
      double d = delta * (vsum[t][N - 1][idx] - vsum[t][cut_v[idx]][idx]);
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
    diff += std::abs(up[idx] - down[idx]) * SabunCostMultiple - std::abs(up[idx] - (down[idx] + delta)) * SabunCostMultiple;
  }
  if (dir == 3) {
    for (int i = 0; i < turn_h[idx].size(); ++i) {
      int t = turn_h[idx][i];
      double d = delta * (hsum[t][idx][N - 1] - hsum[t][idx][cut_h[idx]]);
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
    diff += std::abs(l[idx] - r[idx]) * SabunCostMultiple - std::abs(l[idx] - (r[idx] + delta)) * SabunCostMultiple;
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (dir == 0) {
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        double d = delta * vsum[t][cut_v[idx]][idx];
        dist_est[t] += d;
      }
      up[idx] += delta;
    }
    if (dir == 1) {
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        double d = delta * hsum[t][idx][cut_h[idx]];
        dist_est[t] += d;
      }
      l[idx] += delta;
    }
    if (dir == 2) {
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        double d = delta * (vsum[t][N - 1][idx] - vsum[t][cut_v[idx]][idx]);
        dist_est[t] += d;
      }
      down[idx] += delta;
    }
    if (dir == 3) {
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        double d = delta * (hsum[t][idx][N - 1] - hsum[t][idx][cut_h[idx]]);
        dist_est[t] += d;
      }
      r[idx] += delta;
    }

    diff_sum -= diff;
  }
}

void AnnealingMode1(double temp)
{
  int dir = Rand() % 2;
  int idx = Rand() % N;
  double delta = Rand01() * 20 - 10;

  if (dir == 0) {
    if (up[idx] + delta < 1000 || 9000 < up[idx] + delta) { return; }
    if (down[idx] + delta < 1000 || 9000 < down[idx] + delta) { return; }
  }
  if (dir == 1) {
    if (l[idx] + delta < 1000 || 9000 < l[idx] + delta) { return; }
    if (r[idx] + delta < 1000 || 9000 < r[idx] + delta) { return; }
  }

  double diff = 0;

  if (dir == 0) {
    for (int i = 0; i < turn_v[idx].size(); ++i) {
      int t = turn_v[idx][i];
      double d = delta * vsum[t][N - 1][idx];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }
  if (dir == 1) {
    for (int i = 0; i < turn_h[idx].size(); ++i) {
      int t = turn_h[idx][i];
      double d = delta * hsum[t][idx][N - 1];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (dir == 0) {
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        double d = delta * vsum[t][N - 1][idx];
        dist_est[t] += d;
      }
      up[idx] += delta;
      down[idx] += delta;
    }
    if (dir == 1) {
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        double d = delta * hsum[t][idx][N - 1];
        dist_est[t] += d;
      }
      l[idx] += delta;
      r[idx] += delta;
    }

    diff_sum -= diff;
  }
}

void AnnealingMode2(double temp)
{
  int dir = Rand() % 2;
  int idx = Rand() % 30;
  int delta = Rand() % 30 + 1;
  if (Rand() % 2 == 1) delta *= -1;
  if (dir == 0) if (cut_v[idx] + delta < 3 || 27 <= cut_v[idx] + delta) { return; }
  if (dir == 1) if (cut_h[idx] + delta < 3 || 27 <= cut_h[idx] + delta) { return; }

  double diff = 0;

  if (dir == 0) {
    for (int i = 0; i < turn_v[idx].size(); ++i) {
      int t = turn_v[idx][i];
      double d = 0;
      if (delta > 0) {
        d = (up[idx] - down[idx]) * (vsum[t][cut_v[idx] + delta][idx] - vsum[t][cut_v[idx]][idx]);
      }
      else {
        d = (down[idx] - up[idx]) * (vsum[t][cut_v[idx]][idx] - vsum[t][cut_v[idx] + delta][idx]);
      }
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }
  if (dir == 1) {
    for (int i = 0; i < turn_h[idx].size(); ++i) {
      int t = turn_h[idx][i];
      double d = 0;
      if (delta > 0) {
        d = (l[idx] - r[idx]) * (hsum[t][idx][cut_h[idx] + delta] - hsum[t][idx][cut_h[idx]]);
      }
      else {
        d = (r[idx] - l[idx]) * (hsum[t][idx][cut_h[idx]] - hsum[t][idx][cut_h[idx] + delta]);
      }
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (dir == 0) {
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        double d = 0;
        if (delta > 0) {
          d = (up[idx] - down[idx]) * (vsum[t][cut_v[idx] + delta][idx] - vsum[t][cut_v[idx]][idx]);
        }
        else {
          d = (down[idx] - up[idx]) * (vsum[t][cut_v[idx]][idx] - vsum[t][cut_v[idx] + delta][idx]);
        }
        dist_est[t] += d;
      }
      cut_v[idx] += delta;
    }
    if (dir == 1) {
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        double d = 0;
        if (delta > 0) {
          d = (l[idx] - r[idx]) * (hsum[t][idx][cut_h[idx] + delta] - hsum[t][idx][cut_h[idx]]);
        }
        else {
          d = (r[idx] - l[idx]) * (hsum[t][idx][cut_h[idx]] - hsum[t][idx][cut_h[idx] + delta]);
        }
        dist_est[t] += d;
      }
      cut_h[idx] += delta;
    }

    diff_sum -= diff;
  }
}

void AnnealingMode3Vertical(int idx)
{
  vector<double> rem;
  for (int i = 0; i < turn_v[idx].size(); ++i) {
    int t = turn_v[idx][i];
    rem.push_back(dist_res[t] - (dist_est[t] - (up[idx] * vsum[t][cut_v[idx]][idx] + down[idx] * (vsum[t][N - 1][idx] - vsum[t][cut_v[idx]][idx]))));
  }
  double best_up = up[idx], best_down = down[idx];
  int best_cut = cut_v[idx];
  double maxDiff = 0;
  int changed = 0;
  for (int i = 0; i < turn_v[idx].size(); ++i) {
    int t = turn_v[idx][i];
    maxDiff += (abs(rem[i]) - std::abs(dist_res[t] - dist_est[t])) * (40000.0 / dist_res[t]);
  }
  double keep_diff = maxDiff;
  for (int cut = 3; cut < 27; ++cut) {
    for (int trial = 0; trial < 20; ++trial) {
      double u = Rand() % 8001 + 1000;

      double sum = 0;
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        if (vsum[t][N - 1][idx] - vsum[t][cut][idx] == 0) vec[i].first = 1001001;
        else vec[i].first = (rem[i] - u * vsum[t][cut][idx]) / (vsum[t][N - 1][idx] - vsum[t][cut][idx]);
        vec[i].second = (vsum[t][N - 1][idx] - vsum[t][cut][idx]) * (40000.0 / dist_res[t]);
        sum += vec[i].second;
      }
      if (sum == 0) { continue; }
      vec[turn_v[idx].size()].first = u;
      vec[turn_v[idx].size()].second = SabunCostMultiple;
      sum += SabunCostMultiple;
      sort(vec.begin(), vec.begin() + turn_v[idx].size() + 1);
      int ite = 0;
      double cnt = vec[ite].second;
      while (ite < turn_v[idx].size() + 1) {
        if (cnt < sum - cnt) {
          ite++;
          cnt += vec[ite].second;
        }
        else {
          break;
        }
      }
      double d = 0;
      double dval = vec[ite].first;
      dval = max(dval, 1000.0);
      dval = min(dval, 9000.0);
      for (int i = 0; i < turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        d += (abs(rem[i]) - std::abs(rem[i] - (vsum[t][N - 1][idx] - vsum[t][cut][idx]) * dval)) * (40000.0 / dist_res[t]);
      }
      d += std::abs(up[idx] - down[idx]) * SabunCostMultiple - std::abs(u - dval) * SabunCostMultiple;

      if (d > maxDiff) {
        maxDiff = d;
        best_up = u;
        best_down = dval;
        best_cut = cut;
        changed = 1;
      }
    }
  }
  if (changed == 0) { return; }
  int cdiff = best_cut - cut_v[idx];
  for (int i = 0; i < turn_v[idx].size(); ++i) {
    int t = turn_v[idx][i];
    double d = 0;
    if (cdiff > 0) {
      d = (up[idx] - down[idx]) * (vsum[t][cut_v[idx] + cdiff][idx] - vsum[t][cut_v[idx]][idx]);
    }
    else {
      d = (down[idx] - up[idx]) * (vsum[t][cut_v[idx]][idx] - vsum[t][cut_v[idx] + cdiff][idx]);
    }
    dist_est[t] += d;
  }
  cut_v[idx] = best_cut;
  double udiff = best_up - up[idx];
  for (int i = 0; i < turn_v[idx].size(); ++i) {
    int t = turn_v[idx][i];
    double d = udiff * vsum[t][cut_v[idx]][idx];
    dist_est[t] += d;
  }
  up[idx] = best_up;
  double ddiff = best_down - down[idx];
  for (int i = 0; i < turn_v[idx].size(); ++i) {
    int t = turn_v[idx][i];
    double d = ddiff * (vsum[t][N - 1][idx] - vsum[t][cut_v[idx]][idx]);
    dist_est[t] += d;
  }
  down[idx] = best_down;
  diff_sum -= (maxDiff - keep_diff);
}

void AnnealingMode3Horizontal(int idx)
{
  vector<double> rem;
  for (int i = 0; i < turn_h[idx].size(); ++i) {
    int t = turn_h[idx][i];
    rem.push_back(dist_res[t] - (dist_est[t] - (l[idx] * hsum[t][idx][cut_h[idx]] + r[idx] * (hsum[t][idx][N - 1] - hsum[t][idx][cut_h[idx]]))));
  }
  double best_l = l[idx], best_r = r[idx];
  int best_cut = cut_h[idx];
  double maxDiff = 0;
  int changed = 0;
  for (int i = 0; i < turn_h[idx].size(); ++i) {
    int t = turn_h[idx][i];
    maxDiff += (abs(rem[i]) - std::abs(dist_res[t] - dist_est[t])) * (40000.0 / dist_res[t]);
  }
  double keep_diff = maxDiff;
  for (int cut = 3; cut < 27; ++cut) {
    for (int trial = 0; trial < 20; ++trial) {
      double lval = Rand() % 8001 + 1000;
      double sum = 0;
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        if (hsum[t][idx][N - 1] - hsum[t][idx][cut] == 0) vec[i].first = 1001001;
        else vec[i].first = (rem[i] - lval * hsum[t][idx][cut]) / (hsum[t][idx][N - 1] - hsum[t][idx][cut]);
        vec[i].second = (hsum[t][idx][N - 1] - hsum[t][idx][cut]) * (40000.0 / dist_res[t]);
        sum += vec[i].second;
      }
      if (sum == 0) { continue; }
      vec[turn_h[idx].size()].first = lval;
      vec[turn_h[idx].size()].second = SabunCostMultiple;
      sum += SabunCostMultiple;
      sort(vec.begin(), vec.begin() + turn_h[idx].size() + 1);
      int ite = 0;
      double cnt = vec[ite].second;
      while (ite < turn_h[idx].size() + 1) {
        if (cnt < sum - cnt) {
          ite++;
          cnt += vec[ite].second;
        }
        else {
          break;
        }
      }
      double d = 0;
      double rval = vec[ite].first;
      rval = max(rval, 1000.0);
      rval = min(rval, 9000.0);
      for (int i = 0; i < turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        d += (abs(rem[i]) - std::abs(rem[i] - (hsum[t][idx][N - 1] - hsum[t][idx][cut]) * rval)) * (40000.0 / dist_res[t]);
      }
      d += std::abs(l[idx] - r[idx]) * SabunCostMultiple - std::abs(lval - rval) * SabunCostMultiple;

      if (d > maxDiff) {
        maxDiff = d;
        best_l = lval;
        best_r = rval;
        best_cut = cut;
        changed = 1;
      }
    }
  }
  if (changed == 0) { return; }
  int cdiff = best_cut - cut_h[idx];
  for (int i = 0; i < turn_h[idx].size(); ++i) {
    int t = turn_h[idx][i];
    double d = 0;
    if (cdiff > 0) {
      d = (l[idx] - r[idx]) * (hsum[t][idx][cut_h[idx] + cdiff] - hsum[t][idx][cut_h[idx]]);
    }
    else {
      d = (r[idx] - l[idx]) * (hsum[t][idx][cut_h[idx]] - hsum[t][idx][cut_h[idx] + cdiff]);
    }
    dist_est[t] += d;
  }
  cut_h[idx] = best_cut;
  double ldiff = best_l - l[idx];
  for (int i = 0; i < turn_h[idx].size(); ++i) {
    int t = turn_h[idx][i];
    double d = ldiff * hsum[t][idx][cut_h[idx]];
    dist_est[t] += d;
  }
  l[idx] = best_l;
  double rdiff = best_r - r[idx];
  for (int i = 0; i < turn_h[idx].size(); ++i) {
    int t = turn_h[idx][i];
    double d = rdiff * (hsum[t][idx][N - 1] - hsum[t][idx][cut_h[idx]]);
    dist_est[t] += d;
  }
  r[idx] = best_r;
  diff_sum -= (maxDiff - keep_diff);
}

void AnnealingMode3()
{
  int dir = Rand() % 2;
  int idx = Rand() % N;
  if (dir == 0) {
    AnnealingMode3Vertical(idx);
  }
  else {
    AnnealingMode3Horizontal(idx);
  }
}

int main()
{
  clock_t start, end;
  start = clock();

  int mode = 1;

  if (mode == 0) { 
    Solve("noinput");
  }
  else if (mode == 1) {
    vector<P> ranking;
    ll allScore = 0;
    for (int i = 0; i < 10; ++i) {
      string inputFileNum;
      inputFileNum += (char)((i % 10000) / 1000 + '0');
      inputFileNum += (char)((i % 1000) / 100 + '0');
      inputFileNum += (char)((i % 100) / 10 + '0');
      inputFileNum += (char)((i % 10) / 1 + '0');
      Solve(inputFileNum);
      cout << (int)scoreSumGlobal << " " << i << endl;
      ranking.push_back(P(scoreSumGlobal, i));
      allScore += scoreSumGlobal;
      ClearGlobalVariables();
    }
    sort(ranking.begin(), ranking.end());
    for (int i = 0; i < 10; ++i) cout << ranking[i].first << ' ' << ranking[i].second << endl;
    cout << allScore << endl;
  }

  end = clock();
  if (mode != 0) cout << (double)(end - start) / CLOCKS_PER_SEC << endl;
}
