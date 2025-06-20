#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
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
  const int n = 30;
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
  double dReal[2][n + 1][n + 1]; // 正解の長さ
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
  double up[n] = {}, down[n] = {}, l[n] = {}, r[n] = {};
  int cut_v[n] = {}, cut_h[n] = {}; // 0‾30をとる半開区間
  int vsum[Q][n + 1][n + 1];
  int hsum[Q][n + 1][n + 1];
  vector<int> turn_v[n], turn_h[n];

  double best_up[n], best_down[n], best_l[n], best_r[n];
  int best_cut_v[n], best_cut_h[n];

  double diff_sum = 0;
  double dist_est[Q] = {};

  vector<pair<double, double>> vec(1100);

  double dUD[n + 1][n + 1], dLR[n + 1][n + 1];
  vector<int> PathIDVectorUD[n][n], PathIDVectorLR[n][n];
}

void ClearGlobalVariables()
{
  scoreSumGlobal = 0;
  for (int i = 0; i < 2; ++i)for (int j = 0; j < n + 1; ++j)for (int k = 0; k < n + 1; ++k) dReal[i][j][k] = 0;
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < n + 1; ++j) {
      for (int k = 0; k < n + 1; ++k) {
        vsum[i][j][k] = 0;
        hsum[i][j][k] = 0;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    turn_v[i].clear();
    turn_h[i].clear();
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      PathIDVectorLR[i][j].clear();
      PathIDVectorUD[i][j].clear();
      dUD[i][j] = 0;
      dLR[i][j] = 0;
    }
  }
}

bool IsOutOfBounds(int nx, int ny)
{
  if (nx < 0 || n <= nx || ny < 0 || n <= ny) return 1;
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

double dp[n][n];
int nxt[n][n];

void AnnealingMode0(double temp);
void AnnealingMode1(double temp);
void AnnealingMode2(double temp);
void AnnealingMode3();
void AnnealingMode3Vertical(int rN);
void AnnealingMode3Horizontal(int rN);
void FinalAdjustment();

void Dijkstra(int sx, int sy, int gx, int gy)
{
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) dp[i][j] = INF;
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
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dUD[i][j] = 0;
      dLR[i][j] = 0;
    }
  }

  for (int _ = 0; _ < 20000; ++_) {
    int rv = Rand() % 2;
    int rx = 0, ry = 0;
    if (rv == 0) {
      rx = Rand() % 29 + 1;
      ry = Rand() % 30;
    }
    else {
      rx = Rand() % 30;
      ry = Rand() % 29 + 1;
    }
    double ra = Rand01() * 200.0 - 100.0;
    if (rv == 0) {
      if (rx <= cut_v[ry]) {
        if (up[ry] + ra < 1000 || 9000 < up[ry] + ra) { continue; }
      }
      else {
        if (down[ry] + ra < 1000 || 9000 < down[ry] + ra) { continue; }
      }
    }
    else {
      if (ry <= cut_h[rx]) {
        if (l[rx] + ra < 1000 || 9000 < l[rx] + ra) { continue; }
      }
      else {
        if (r[rx] + ra < 1000 || 9000 < r[rx] + ra) { continue; }
      }
    }

    double diff = 0;

    if (rv == 0) {
      for (int i = 0; i < PathIDVectorUD[rx][ry].size(); ++i) {
        int pathID = PathIDVectorUD[rx][ry][i];
        diff += (std::abs(dist_res[pathID] - dist_est[pathID]) - std::abs(dist_res[pathID] - (dist_est[pathID] + ra))) * (40000.0 / dist_res[pathID]);
      }
    }
    else {
      for (int i = 0; i < PathIDVectorLR[rx][ry].size(); ++i) {
        int pathID = PathIDVectorLR[rx][ry][i];
        diff += (std::abs(dist_res[pathID] - dist_est[pathID]) - std::abs(dist_res[pathID] - (dist_est[pathID] + ra))) * (40000.0 / dist_res[pathID]);
      }
    }

    if (diff > 0) {
      if (rv == 0) {
        for (int i = 0; i < PathIDVectorUD[rx][ry].size(); ++i) {
          int pathID = PathIDVectorUD[rx][ry][i];
          dist_est[pathID] += ra;
        }
      }
      else {
        for (int i = 0; i < PathIDVectorLR[rx][ry].size(); ++i) {
          int pathID = PathIDVectorLR[rx][ry][i];
          dist_est[pathID] += ra;
        }
      }
      if (rv == 0) dUD[rx][ry] += ra;
      else dLR[rx][ry] += ra;
      diff_sum -= diff;
    }
  }
}

void Dijkstra2(int sx, int sy, int gx, int gy)
{
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) dp[i][j] = INF;
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
    int ggx = gx - dx[nxt[gx][gy]];
    int ggy = gy - dy[nxt[gx][gy]];
    gx = ggx;
    gy = ggy;
  }
  reverse(ans.begin(), ans.end());
  reverse(v.begin(), v.end());
  return ans;
}

void UpdatePathInfo(int turn, int sx, int sy, const vector<int>& v)
{
  int vpath[n + 1][n + 1];
  int hpath[n + 1][n + 1];
  for (int i = 0; i < n + 1; ++i) {
    for (int j = 0; j < n + 1; ++j) {
      vpath[i][j] = 0;
      hpath[i][j] = 0;
    }
  }

  int ksx = sx, ksy = sy;
  for (int i = 0; i < v.size(); ++i) {
    int nsx = sx;
    int nsy = sy;
    if (v[i] == 2) nsx += 1;
    if (v[i] == 3) nsy += 1;

    if (v[i] % 2 == 0) {
      vpath[nsx][nsy] = 1;
      PathIDVectorUD[nsx][nsy].push_back(turn);
    }
    else {
      hpath[nsx][nsy] = 1;
      PathIDVectorLR[nsx][nsy].push_back(turn);
    }

    sx += dx[v[i]];
    sy += dy[v[i]];
  }

  for (int j = 0; j < n; ++j) {
    for (int i = 1; i < n; ++i) {
      vsum[turn][i][j] = vsum[turn][i - 1][j] + vpath[i][j];
    }
    if (vsum[turn][n - 1][j]) turn_v[j].push_back(turn);
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 1; j < n; ++j) {
      hsum[turn][i][j] = hsum[turn][i][j - 1] + hpath[i][j];
    }
    if (hsum[turn][i][n - 1]) turn_h[i].push_back(turn);
  }
}

void UpdateDiffSum(int turn)
{
  diff_sum = 0;
  for (int i = 0; i < Q; ++i) dist_est[i] = 0;

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < turn_v[j].size(); ++i) {
      int turnID = turn_v[j][i];
      dist_est[turnID] += up[j] * vsum[turnID][cut_v[j]][j] + down[j] * (vsum[turnID][n - 1][j] - vsum[turnID][cut_v[j]][j]);
    }
  }
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < turn_h[j].size(); ++i) {
      int turnID = turn_h[j][i];
      dist_est[turnID] += l[j] * hsum[turnID][j][cut_h[j]] + r[j] * (hsum[turnID][j][n - 1] - hsum[turnID][j][cut_h[j]]);
    }
  }
  for (int i = 0; i < turn + 1; ++i) {
    diff_sum += std::abs(dist_res[i] - dist_est[i]) * (40000.0 / dist_res[i]);
  }

  for (int i = 0; i < n; ++i) {
    diff_sum += std::abs(up[i] - down[i]) * SabunCostMultiple;
    diff_sum += std::abs(l[i] - r[i]) * SabunCostMultiple;
  }
}

void SaveBestParams()
{
  for (int i = 0; i < n; ++i) {
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
  for (int i = 0; i < n; ++i) {
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
    for (int i = 0; i < n; ++i) for (int j = 1; j < n; ++j) ifs >> dReal[1][i][j];
    for (int i = 1; i < n; ++i) for (int j = 0; j < n; ++j) ifs >> dReal[0][i][j];
  }

  string fileName = (string)"./out/" + inputFileNum + "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);

  for (int i = 0; i < n; ++i) {
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
    for (int i = 0; i < n; ++i) ofsParam << i << ' ' << cut_v[i] << ' ' << up[i] << ' ' << down[i] << endl;
    ofsParam << "left / right" << endl;
    for (int i = 0; i < n; ++i) ofsParam << i << ' ' << cut_h[i] << ' ' << l[i] << ' ' << r[i] << endl;

    ofsParam << "dUD" << endl;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        ofsParam << dUD[i][j] << ' ';
      }
      ofsParam << endl;
    }
    ofsParam << "dLR" << endl;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
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
  int rULDR = Rand() % 4;
  int rN = Rand() % n;
  double rA = Rand01() * 20 - 10;

  if (rULDR == 0) if (up[rN] + rA < 1000 || 9000 < up[rN] + rA) { return; }
  if (rULDR == 1) if (l[rN] + rA < 1000 || 9000 < l[rN] + rA) { return; }
  if (rULDR == 2) if (down[rN] + rA < 1000 || 9000 < down[rN] + rA) { return; }
  if (rULDR == 3) if (r[rN] + rA < 1000 || 9000 < r[rN] + rA) { return; }

  double diff = 0;

  if (rULDR == 0) {
    for (int i = 0; i < turn_v[rN].size(); ++i) {
      int turnID = turn_v[rN][i];
      double d = rA * vsum[turnID][cut_v[rN]][rN];
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
    diff += std::abs(up[rN] - down[rN]) * SabunCostMultiple - std::abs((up[rN] + rA) - down[rN]) * SabunCostMultiple;
  }
  if (rULDR == 1) {
    for (int i = 0; i < turn_h[rN].size(); ++i) {
      int turnID = turn_h[rN][i];
      double d = rA * hsum[turnID][rN][cut_h[rN]];
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
    diff += std::abs(l[rN] - r[rN]) * SabunCostMultiple - std::abs((l[rN] + rA) - r[rN]) * SabunCostMultiple;
  }
  if (rULDR == 2) {
    for (int i = 0; i < turn_v[rN].size(); ++i) {
      int turnID = turn_v[rN][i];
      double d = rA * (vsum[turnID][n - 1][rN] - vsum[turnID][cut_v[rN]][rN]);
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
    diff += std::abs(up[rN] - down[rN]) * SabunCostMultiple - std::abs(up[rN] - (down[rN] + rA)) * SabunCostMultiple;
  }
  if (rULDR == 3) {
    for (int i = 0; i < turn_h[rN].size(); ++i) {
      int turnID = turn_h[rN][i];
      double d = rA * (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut_h[rN]]);
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
    diff += std::abs(l[rN] - r[rN]) * SabunCostMultiple - std::abs(l[rN] - (r[rN] + rA)) * SabunCostMultiple;
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (rULDR == 0) {
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        double d = rA * vsum[turnID][cut_v[rN]][rN];
        dist_est[turnID] += d;
      }
    }
    if (rULDR == 1) {
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        double d = rA * hsum[turnID][rN][cut_h[rN]];
        dist_est[turnID] += d;
      }
    }
    if (rULDR == 2) {
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        double d = rA * (vsum[turnID][n - 1][rN] - vsum[turnID][cut_v[rN]][rN]);
        dist_est[turnID] += d;
      }
    }
    if (rULDR == 3) {
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        double d = rA * (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut_h[rN]]);
        dist_est[turnID] += d;
      }
    }
    if (rULDR == 0) up[rN] += rA;
    if (rULDR == 1) l[rN] += rA;
    if (rULDR == 2) down[rN] += rA;
    if (rULDR == 3) r[rN] += rA;
    diff_sum -= diff;
  }
}

void AnnealingMode1(double temp)
{
  int rUL = Rand() % 2;
  int rN = Rand() % n;
  double rA = Rand01() * 20 - 10;

  if (rUL == 0) if (up[rN] + rA < 1000 || 9000 < up[rN] + rA) { return; }
  if (rUL == 1) if (l[rN] + rA < 1000 || 9000 < l[rN] + rA) { return; }
  if (rUL == 0) if (down[rN] + rA < 1000 || 9000 < down[rN] + rA) { return; }
  if (rUL == 1) if (r[rN] + rA < 1000 || 9000 < r[rN] + rA) { return; }

  double diff = 0;

  if (rUL == 0) {
    for (int i = 0; i < turn_v[rN].size(); ++i) {
      int turnID = turn_v[rN][i];
      double d = rA * vsum[turnID][n - 1][rN];
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
  }
  if (rUL == 1) {
    for (int i = 0; i < turn_h[rN].size(); ++i) {
      int turnID = turn_h[rN][i];
      double d = rA * hsum[turnID][rN][n - 1];
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (rUL == 0) {
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        double d = rA * vsum[turnID][n - 1][rN];
        dist_est[turnID] += d;
      }
    }
    if (rUL == 1) {
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        double d = rA * hsum[turnID][rN][n - 1];
        dist_est[turnID] += d;
      }
    }
    if (rUL == 0) up[rN] += rA;
    if (rUL == 1) l[rN] += rA;
    if (rUL == 0) down[rN] += rA;
    if (rUL == 1) r[rN] += rA;
    diff_sum -= diff;
  }
}

void AnnealingMode2(double temp)
{
  int rUL = Rand() % 2;
  int rN = Rand() % 30;
  int rA = Rand() % 30 + 1;
  if (Rand() % 2 == 1) rA *= -1;
  if (rUL == 0) if (cut_v[rN] + rA < 3 || 27 <= cut_v[rN] + rA) { return; }
  if (rUL == 1) if (cut_h[rN] + rA < 3 || 27 <= cut_h[rN] + rA) { return; }

  double diff = 0;

  if (rUL == 0) {
    for (int i = 0; i < turn_v[rN].size(); ++i) {
      int turnID = turn_v[rN][i];
      double d = 0;
      if (rA > 0) {
        d = (up[rN] - down[rN]) * (vsum[turnID][cut_v[rN] + rA][rN] - vsum[turnID][cut_v[rN]][rN]);
      }
      else {
        d = (down[rN] - up[rN]) * (vsum[turnID][cut_v[rN]][rN] - vsum[turnID][cut_v[rN] + rA][rN]);
      }
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
  }
  if (rUL == 1) {
    for (int i = 0; i < turn_h[rN].size(); ++i) {
      int turnID = turn_h[rN][i];
      double d = 0;
      if (rA > 0) {
        d = (l[rN] - r[rN]) * (hsum[turnID][rN][cut_h[rN] + rA] - hsum[turnID][rN][cut_h[rN]]);
      }
      else {
        d = (r[rN] - l[rN]) * (hsum[turnID][rN][cut_h[rN]] - hsum[turnID][rN][cut_h[rN] + rA]);
      }
      diff += (std::abs(dist_res[turnID] - dist_est[turnID]) - std::abs(dist_res[turnID] - (dist_est[turnID] + d))) * (40000.0 / dist_res[turnID]);
    }
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (rUL == 0) {
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        double d = 0;
        if (rA > 0) {
          d = (up[rN] - down[rN]) * (vsum[turnID][cut_v[rN] + rA][rN] - vsum[turnID][cut_v[rN]][rN]);
        }
        else {
          d = (down[rN] - up[rN]) * (vsum[turnID][cut_v[rN]][rN] - vsum[turnID][cut_v[rN] + rA][rN]);
        }
        dist_est[turnID] += d;
      }
    }
    if (rUL == 1) {
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        double d = 0;
        if (rA > 0) {
          d = (l[rN] - r[rN]) * (hsum[turnID][rN][cut_h[rN] + rA] - hsum[turnID][rN][cut_h[rN]]);
        }
        else {
          d = (r[rN] - l[rN]) * (hsum[turnID][rN][cut_h[rN]] - hsum[turnID][rN][cut_h[rN] + rA]);
        }
        dist_est[turnID] += d;
      }
    }

    if (rUL == 0) cut_v[rN] += rA;
    if (rUL == 1) cut_h[rN] += rA;

    diff_sum -= diff;
  }
}

void AnnealingMode3Vertical(int rN)
{
  vector<double> amari;
  for (int i = 0; i < turn_v[rN].size(); ++i) {
    int turnID = turn_v[rN][i];
    amari.push_back(dist_res[turnID] - (dist_est[turnID] - (up[rN] * vsum[turnID][cut_v[rN]][rN] + down[rN] * (vsum[turnID][n - 1][rN] - vsum[turnID][cut_v[rN]][rN]))));
  }
  double arg_up = up[rN], arg_down = down[rN];
  int arg_cut_v = cut_v[rN];
  double maxDiff = 0;
  int changed = 0;
  for (int i = 0; i < turn_v[rN].size(); ++i) {
    int turnID = turn_v[rN][i];
    maxDiff += (abs(amari[i]) - std::abs(dist_res[turnID] - dist_est[turnID])) * (40000.0 / dist_res[turnID]);
  }
  double keep_diff = maxDiff;
  for (int cut = 3; cut < 27; ++cut) {
    for (int randomChallenge = 0; randomChallenge < 20; ++randomChallenge) {
      double rA = Rand() % 8001 + 1000;

      double countSum = 0;
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        if (vsum[turnID][n - 1][rN] - vsum[turnID][cut][rN] == 0) vec[i].first = 1001001;
        else vec[i].first = (amari[i] - rA * vsum[turnID][cut][rN]) / (vsum[turnID][n - 1][rN] - vsum[turnID][cut][rN]);
        vec[i].second = (vsum[turnID][n - 1][rN] - vsum[turnID][cut][rN]) * (40000.0 / dist_res[turnID]);
        countSum += vec[i].second;
      }
      if (countSum == 0) { continue; }
      vec[turn_v[rN].size()].first = rA;
      vec[turn_v[rN].size()].second = SabunCostMultiple;
      countSum += SabunCostMultiple;
      sort(vec.begin(), vec.begin() + turn_v[rN].size() + 1);
      int ite = 0;
      double countNow = vec[ite].second;
      while (ite < turn_v[rN].size() + 1) {
        if (countNow < countSum - countNow) {
          ite++;
          countNow += vec[ite].second;
        }
        else {
          break;
        }
      }
      double d = 0;
      double DownValue = vec[ite].first;
      DownValue = max(DownValue, 1000.0);
      DownValue = min(DownValue, 9000.0);
      for (int i = 0; i < turn_v[rN].size(); ++i) {
        int turnID = turn_v[rN][i];
        d += (abs(amari[i]) - std::abs(amari[i] - (vsum[turnID][n - 1][rN] - vsum[turnID][cut][rN]) * DownValue)) * (40000.0 / dist_res[turnID]);
      }
      d += std::abs(up[rN] - down[rN]) * SabunCostMultiple - std::abs(rA - DownValue) * SabunCostMultiple;

      if (d > maxDiff) {
        maxDiff = d;
        arg_up = rA;
        arg_down = DownValue;
        arg_cut_v = cut;
        changed = 1;
      }
    }
  }
  if (changed == 0) { return; }
  int cutDiff = arg_cut_v - cut_v[rN];
  for (int i = 0; i < turn_v[rN].size(); ++i) {
    int turnID = turn_v[rN][i];
    double d = 0;
    if (cutDiff > 0) {
      d = (up[rN] - down[rN]) * (vsum[turnID][cut_v[rN] + cutDiff][rN] - vsum[turnID][cut_v[rN]][rN]);
    }
    else {
      d = (down[rN] - up[rN]) * (vsum[turnID][cut_v[rN]][rN] - vsum[turnID][cut_v[rN] + cutDiff][rN]);
    }
    dist_est[turnID] += d;
  }
  cut_v[rN] = arg_cut_v;
  double upDiff = arg_up - up[rN];
  for (int i = 0; i < turn_v[rN].size(); ++i) {
    int turnID = turn_v[rN][i];
    double d = upDiff * vsum[turnID][cut_v[rN]][rN];
    dist_est[turnID] += d;
  }
  up[rN] = arg_up;
  double downDiff = arg_down - down[rN];
  for (int i = 0; i < turn_v[rN].size(); ++i) {
    int turnID = turn_v[rN][i];
    double d = downDiff * (vsum[turnID][n - 1][rN] - vsum[turnID][cut_v[rN]][rN]);
    dist_est[turnID] += d;
  }
  down[rN] = arg_down;
  diff_sum -= (maxDiff - keep_diff);
}

void AnnealingMode3Horizontal(int rN)
{
  vector<double> amari;
  for (int i = 0; i < turn_h[rN].size(); ++i) {
    int turnID = turn_h[rN][i];
    amari.push_back(dist_res[turnID] - (dist_est[turnID] - (l[rN] * hsum[turnID][rN][cut_h[rN]] + r[rN] * (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut_h[rN]]))));
  }
  double argmax_l = l[rN], argmax_r = r[rN];
  int arg_cut_h = cut_h[rN];
  double maxDiff = 0;
  int changed = 0;
  for (int i = 0; i < turn_h[rN].size(); ++i) {
    int turnID = turn_h[rN][i];
    maxDiff += (abs(amari[i]) - std::abs(dist_res[turnID] - dist_est[turnID])) * (40000.0 / dist_res[turnID]);
  }
  double keep_diff = maxDiff;
  for (int cut = 3; cut < 27; ++cut) {
    for (int randomChallenge = 0; randomChallenge < 20; ++randomChallenge) {
      double rA = Rand() % 8001 + 1000;
      double countSum = 0;
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        if (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut] == 0) vec[i].first = 1001001;
        else vec[i].first = (amari[i] - rA * hsum[turnID][rN][cut]) / (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut]);
        vec[i].second = (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut]) * (40000.0 / dist_res[turnID]);
        countSum += vec[i].second;
      }
      if (countSum == 0) { continue; }
      vec[turn_h[rN].size()].first = rA;
      vec[turn_h[rN].size()].second = SabunCostMultiple;
      countSum += SabunCostMultiple;
      sort(vec.begin(), vec.begin() + turn_h[rN].size() + 1);
      int ite = 0;
      double countNow = vec[ite].second;
      while (ite < turn_h[rN].size() + 1) {
        if (countNow < countSum - countNow) {
          ite++;
          countNow += vec[ite].second;
        }
        else {
          break;
        }
      }
      double d = 0;
      double RightValue = vec[ite].first;
      RightValue = max(RightValue, 1000.0);
      RightValue = min(RightValue, 9000.0);
      for (int i = 0; i < turn_h[rN].size(); ++i) {
        int turnID = turn_h[rN][i];
        d += (abs(amari[i]) - std::abs(amari[i] - (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut]) * RightValue)) * (40000.0 / dist_res[turnID]);
      }
      d += std::abs(l[rN] - r[rN]) * SabunCostMultiple - std::abs(rA - RightValue) * SabunCostMultiple;

      if (d > maxDiff) {
        maxDiff = d;
        argmax_l = rA;
        argmax_r = RightValue;
        arg_cut_h = cut;
        changed = 1;
      }
    }
  }
  if (changed == 0) { return; }
  int cutDiff = arg_cut_h - cut_h[rN];
  for (int i = 0; i < turn_h[rN].size(); ++i) {
    int turnID = turn_h[rN][i];
    double d = 0;
    if (cutDiff > 0) {
      d = (l[rN] - r[rN]) * (hsum[turnID][rN][cut_h[rN] + cutDiff] - hsum[turnID][rN][cut_h[rN]]);
    }
    else {
      d = (r[rN] - l[rN]) * (hsum[turnID][rN][cut_h[rN]] - hsum[turnID][rN][cut_h[rN] + cutDiff]);
    }
    dist_est[turnID] += d;
  }
  cut_h[rN] = arg_cut_h;
  double leftDiff = argmax_l - l[rN];
  for (int i = 0; i < turn_h[rN].size(); ++i) {
    int turnID = turn_h[rN][i];
    double d = leftDiff * hsum[turnID][rN][cut_h[rN]];
    dist_est[turnID] += d;
  }
  l[rN] = argmax_l;
  double rightDiff = argmax_r - r[rN];
  for (int i = 0; i < turn_h[rN].size(); ++i) {
    int turnID = turn_h[rN][i];
    double d = rightDiff * (hsum[turnID][rN][n - 1] - hsum[turnID][rN][cut_h[rN]]);
    dist_est[turnID] += d;
  }
  r[rN] = argmax_r;
  diff_sum -= (maxDiff - keep_diff);
}

void AnnealingMode3()
{
  int rUL = Rand() % 2;
  int rN = Rand() % n;
  if (rUL == 0) {
    AnnealingMode3Vertical(rN);
  }
  else {
    AnnealingMode3Horizontal(rN);
  }
}

int main()
{
  clock_t start, end;
  start = clock();

  int mode = 1;

  if (mode == 0) { // 提出用
    Solve("noinput");
  }
  else if (mode == 1) { // サンプル0‾99でチェック
    vector<P> ranking;
    ll allScore = 0;
    for (int i = 0; i < 100; ++i) {
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
