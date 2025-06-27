#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <queue>
#include <string>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<double, P> PDP;

#define MAX_N 200005
const int INF = 1001001001;

class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

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
  array<double, Q> dist_res = {};				// 受け取った長さ
}

namespace /* 行、列、切れ目の構造 */
{
  class Edge
  {
  public:
    array<double, N> up = {};
    array<double, N> down = {};
    array<double, N> left = {};
    array<double, N> right = {};
    array<int, N> cut_v = {};
    array<int, N> cut_h = {}; // 0‾30をとる半開区間

    void copyFrom(const Edge& other)
    {
      for (int i = 0; i < N; ++i) {
        up[i] = other.up[i];
        down[i] = other.down[i];
        left[i] = other.left[i];
        right[i] = other.right[i];
        cut_v[i] = other.cut_v[i];
        cut_h[i] = other.cut_h[i];
      }
    }
  };
  Edge edge;
  array<array<array<int, N + 1>, N + 1>, Q> vsum;
  array<array<array<int, N + 1>, N + 1>, Q> hsum;
  vector<int> turn_v[N], turn_h[N];

  Edge best_edge;

  double diff_sum = 0;
  array<double, Q> dist_est = {};

  vector<pair<double, double>> vec(1100);

  array<array<double, N + 1>, N + 1> dUD, dLR;
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
  for (int i = 0; i < (int)ans.size(); ++i) {
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
      if (x <= edge.cut_v[y]) {
        if (edge.up[y] + delta < 1000 || 9000 < edge.up[y] + delta) { continue; }
      }
      else {
        if (edge.down[y] + delta < 1000 || 9000 < edge.down[y] + delta) { continue; }
      }
    }
    else {
      if (y <= edge.cut_h[x]) {
        if (edge.left[x] + delta < 1000 || 9000 < edge.left[x] + delta) { continue; }
      }
      else {
        if (edge.right[x] + delta < 1000 || 9000 < edge.right[x] + delta) { continue; }
      }
    }

    double diff = 0;

    if (dir == 0) {
      for (int i = 0; i < (int)PathIDVectorUD[x][y].size(); ++i) {
        int t = PathIDVectorUD[x][y][i];
        diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + delta))) * (40000.0 / dist_res[t]);
      }
    }
    else {
      for (int i = 0; i < (int)PathIDVectorLR[x][y].size(); ++i) {
        int t = PathIDVectorLR[x][y][i];
        diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + delta))) * (40000.0 / dist_res[t]);
      }
    }

    if (diff > 0) {
      if (dir == 0) {
        for (int i = 0; i < (int)PathIDVectorUD[x][y].size(); ++i) {
          int t = PathIDVectorUD[x][y][i];
          dist_est[t] += delta;
        }
      }
      else {
        for (int i = 0; i < (int)PathIDVectorLR[x][y].size(); ++i) {
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
        if (x <= edge.cut_v[y]) dd = edge.up[y];
        else dd = edge.down[y];
        dd += dUD[x][y];
      }
      if (i == 1) {
        if (y <= edge.cut_h[x]) dd = edge.left[x];
        else dd = edge.right[x];
        dd += dLR[x][y];
      }
      if (i == 2) {
        if (x + 1 <= edge.cut_v[y]) dd = edge.up[y];
        else dd = edge.down[y];
        dd += dUD[x + 1][y];
      }
      if (i == 3) {
        if (y + 1 <= edge.cut_h[x]) dd = edge.left[x];
        else dd = edge.right[x];
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

  for (int i = 0; i < (int)v.size(); ++i) {
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
    for (int i = 0; i < (int)turn_v[j].size(); ++i) {
      int t = turn_v[j][i];
      dist_est[t] += edge.up[j] * vsum[t][edge.cut_v[j]][j] + edge.down[j] * (vsum[t][N - 1][j] - vsum[t][edge.cut_v[j]][j]);
    }
  }
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < (int)turn_h[j].size(); ++i) {
      int t = turn_h[j][i];
      dist_est[t] += edge.left[j] * hsum[t][j][edge.cut_h[j]] + edge.right[j] * (hsum[t][j][N - 1] - hsum[t][j][edge.cut_h[j]]);
    }
  }
  for (int i = 0; i < turn + 1; ++i) {
    diff_sum += std::abs(dist_res[i] - dist_est[i]) * (40000.0 / dist_res[i]);
  }

  for (int i = 0; i < N; ++i) {
    diff_sum += std::abs(edge.up[i] - edge.down[i]) * SabunCostMultiple;
    diff_sum += std::abs(edge.left[i] - edge.right[i]) * SabunCostMultiple;
  }
}

void SaveBestParams()
{
  best_edge.copyFrom(edge);
}

void RestoreBestParams()
{
  edge.copyFrom(best_edge);
}

int Solve(string case_num)
{
  // 時間計測
  Timer timer;
  timer.start();

  int inputMode = 1;
  string fileNameIfs = "./in/" + case_num + ".txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs) { // 標準入力する
    inputMode = 0;
  }
  else {
    for (int i = 0; i < N; ++i) for (int j = 1; j < N; ++j) ifs >> dReal[1][i][j];
    for (int i = 1; i < N; ++i) for (int j = 0; j < N; ++j) ifs >> dReal[0][i][j];
  }

  string fileName = "./out/" + case_num + "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);

  for (int i = 0; i < N; ++i) {
    edge.up[i] = initialD;
    edge.down[i] = initialD;
    edge.left[i] = initialD;
    edge.right[i] = initialD;
    edge.cut_v[i] = 15;
    edge.cut_h[i] = 15;
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

    // int length = ans.size();

    if (turn < 0) {
      continue;
    }
    if (timer.get_elapsed_time() > 1.8) {
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
    double now_time = timer.get_elapsed_time();
    double TL = 1.8;

    for (int loop = 0; loop < loopTimes; ++loop) {
      all_loop++;
      if (loop % 100 == 1) {
        now_time = timer.get_elapsed_time();
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
    string fileNameParam = (string)"./parameter/" + case_num + "_param.txt";
    const char* cstrParam = fileNameParam.c_str();
    ofstream ofsParam(cstrParam);
    ofsParam << "up / down" << endl;
    for (int i = 0; i < N; ++i) ofsParam << i << ' ' << edge.cut_v[i] << ' ' << edge.up[i] << ' ' << edge.down[i] << endl;
    ofsParam << "left / right" << endl;
    for (int i = 0; i < N; ++i) ofsParam << i << ' ' << edge.cut_h[i] << ' ' << edge.left[i] << ' ' << edge.right[i] << endl;

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

  // 方向に応じた参照を取得
  auto& targetEdge = (dir == 0) ? edge.up[idx] :
    (dir == 1) ? edge.left[idx] :
    (dir == 2) ? edge.down[idx] : edge.right[idx];

  // 範囲チェック
  if (targetEdge + delta < 1000 || 9000 < targetEdge + delta) { return; }

  double diff = 0;

  // 方向に応じた処理を共通化
  bool isVertical = (dir == 0 || dir == 2);
  bool isFirst = (dir == 0 || dir == 1);

  auto& turns = isVertical ? turn_v[idx] : turn_h[idx];
  auto& edge1 = isVertical ? edge.up[idx] : edge.left[idx];
  auto& edge2 = isVertical ? edge.down[idx] : edge.right[idx];
  auto& cut = isVertical ? edge.cut_v[idx] : edge.cut_h[idx];

  auto getMultiplier = [&](int t) {
    if (isVertical) {
      return isFirst ? vsum[t][cut][idx] : (vsum[t][N - 1][idx] - vsum[t][cut][idx]);
    }
    else {
      return isFirst ? hsum[t][idx][cut] : (hsum[t][idx][N - 1] - hsum[t][idx][cut]);
    }
    };

  // diff計算
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    double d = delta * getMultiplier(t);
    diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
  }

  // 隣接エッジとの差分コスト計算
  if (isFirst) {
    diff += std::abs(edge1 - edge2) * SabunCostMultiple - std::abs((edge1 + delta) - edge2) * SabunCostMultiple;
  }
  else {
    diff += std::abs(edge1 - edge2) * SabunCostMultiple - std::abs(edge1 - (edge2 + delta)) * SabunCostMultiple;
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    // 更新処理も共通化
    for (int i = 0; i < (int)turns.size(); ++i) {
      int t = turns[i];
      double d = delta * getMultiplier(t);
      dist_est[t] += d;
    }
    targetEdge += delta;
    diff_sum -= diff;
  }
}

void AnnealingMode1(double temp)
{
  int dir = Rand() % 2;
  int idx = Rand() % N;
  double delta = Rand01() * 20 - 10;

  if (dir == 0) {
    if (edge.up[idx] + delta < 1000 || 9000 < edge.up[idx] + delta) { return; }
    if (edge.down[idx] + delta < 1000 || 9000 < edge.down[idx] + delta) { return; }
  }
  if (dir == 1) {
    if (edge.left[idx] + delta < 1000 || 9000 < edge.left[idx] + delta) { return; }
    if (edge.right[idx] + delta < 1000 || 9000 < edge.right[idx] + delta) { return; }
  }

  double diff = 0;

  if (dir == 0) {
    for (int i = 0; i < (int)turn_v[idx].size(); ++i) {
      int t = turn_v[idx][i];
      double d = delta * vsum[t][N - 1][idx];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }
  if (dir == 1) {
    for (int i = 0; i < (int)turn_h[idx].size(); ++i) {
      int t = turn_h[idx][i];
      double d = delta * hsum[t][idx][N - 1];
      diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
    }
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    if (dir == 0) {
      for (int i = 0; i < (int)turn_v[idx].size(); ++i) {
        int t = turn_v[idx][i];
        double d = delta * vsum[t][N - 1][idx];
        dist_est[t] += d;
      }
      edge.up[idx] += delta;
      edge.down[idx] += delta;
    }
    if (dir == 1) {
      for (int i = 0; i < (int)turn_h[idx].size(); ++i) {
        int t = turn_h[idx][i];
        double d = delta * hsum[t][idx][N - 1];
        dist_est[t] += d;
      }
      edge.left[idx] += delta;
      edge.right[idx] += delta;
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

  // 方向に応じた参照を取得
  bool isVertical = (dir == 0);
  auto& targetCut = isVertical ? edge.cut_v[idx] : edge.cut_h[idx];

  // 範囲チェック
  if (targetCut + delta < 3 || 27 <= targetCut + delta) { return; }

  // 方向に応じた処理を共通化
  auto& turns = isVertical ? turn_v[idx] : turn_h[idx];
  auto& edge1 = isVertical ? edge.up[idx] : edge.left[idx];
  auto& edge2 = isVertical ? edge.down[idx] : edge.right[idx];

  auto getSum = [&](int t, int cutPos) {
    return isVertical ? vsum[t][cutPos][idx] : hsum[t][idx][cutPos];
    };

  // 差分計算を共通化
  auto calcDelta = [&](int t) {
    if (delta > 0) {
      return (edge1 - edge2) * (getSum(t, targetCut + delta) - getSum(t, targetCut));
    }
    else {
      return (edge2 - edge1) * (getSum(t, targetCut) - getSum(t, targetCut + delta));
    }
    };

  double diff = 0;

  // diff計算
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    double d = calcDelta(t);
    diff += (std::abs(dist_res[t] - dist_est[t]) - std::abs(dist_res[t] - (dist_est[t] + d))) * (40000.0 / dist_res[t]);
  }

  double prob = exp((double)diff / temp);

  if (prob > Rand01()) {
    // 更新処理も共通化
    for (int i = 0; i < (int)turns.size(); ++i) {
      int t = turns[i];
      double d = calcDelta(t);
      dist_est[t] += d;
    }
    targetCut += delta;
    diff_sum -= diff;
  }
}

template<bool IsVertical>
void AnnealingMode3Common(int idx)
{
  vector<double> rem;
  auto& turns = IsVertical ? turn_v[idx] : turn_h[idx];
  auto& edge1 = IsVertical ? edge.up : edge.left;
  auto& edge2 = IsVertical ? edge.down : edge.right;
  auto& cut = IsVertical ? edge.cut_v : edge.cut_h;

  auto getSum1 = [&](int t, int cutPos) {
    return IsVertical ? vsum[t][cutPos][idx] : hsum[t][idx][cutPos];
    };
  auto getSum2 = [&](int t) {
    return IsVertical ? vsum[t][N - 1][idx] : hsum[t][idx][N - 1];
    };

  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    rem.push_back(dist_res[t] - (dist_est[t] - (edge1[idx] * getSum1(t, cut[idx]) + edge2[idx] * (getSum2(t) - getSum1(t, cut[idx])))));
  }
  double best_edge1 = edge1[idx], best_edge2 = edge2[idx];
  int best_cut = cut[idx];
  double maxDiff = 0;
  int changed = 0;
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    maxDiff += (abs(rem[i]) - std::abs(dist_res[t] - dist_est[t])) * (40000.0 / dist_res[t]);
  }
  double keep_diff = maxDiff;
  for (int cutPos = 3; cutPos < 27; ++cutPos) {
    for (int trial = 0; trial < 20; ++trial) {
      double val1 = Rand() % 8001 + 1000;

      double sum = 0;
      for (int i = 0; i < (int)turns.size(); ++i) {
        int t = turns[i];
        if (getSum2(t) - getSum1(t, cutPos) == 0) vec[i].first = 1001001;
        else vec[i].first = (rem[i] - val1 * getSum1(t, cutPos)) / (getSum2(t) - getSum1(t, cutPos));
        vec[i].second = (getSum2(t) - getSum1(t, cutPos)) * (40000.0 / dist_res[t]);
        sum += vec[i].second;
      }
      if (sum == 0) { continue; }
      vec[turns.size()].first = val1;
      vec[turns.size()].second = SabunCostMultiple;
      sum += SabunCostMultiple;
      sort(vec.begin(), vec.begin() + turns.size() + 1);
      int ite = 0;
      double cnt = vec[ite].second;
      while (ite < (int)turns.size() + 1) {
        if (cnt < sum - cnt) {
          ite++;
          cnt += vec[ite].second;
        }
        else {
          break;
        }
      }
      double d = 0;
      double val2 = vec[ite].first;
      val2 = max(val2, 1000.0);
      val2 = min(val2, 9000.0);
      for (int i = 0; i < (int)turns.size(); ++i) {
        int t = turns[i];
        d += (abs(rem[i]) - std::abs(rem[i] - (getSum2(t) - getSum1(t, cutPos)) * val2)) * (40000.0 / dist_res[t]);
      }
      d += std::abs(edge1[idx] - edge2[idx]) * SabunCostMultiple - std::abs(val1 - val2) * SabunCostMultiple;

      if (d > maxDiff) {
        maxDiff = d;
        best_edge1 = val1;
        best_edge2 = val2;
        best_cut = cutPos;
        changed = 1;
      }
    }
  }
  if (changed == 0) { return; }
  int cdiff = best_cut - cut[idx];
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    double d = 0;
    if (cdiff > 0) {
      d = (edge1[idx] - edge2[idx]) * (getSum1(t, cut[idx] + cdiff) - getSum1(t, cut[idx]));
    }
    else {
      d = (edge2[idx] - edge1[idx]) * (getSum1(t, cut[idx]) - getSum1(t, cut[idx] + cdiff));
    }
    dist_est[t] += d;
  }
  cut[idx] = best_cut;
  double diff1 = best_edge1 - edge1[idx];
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    double d = diff1 * getSum1(t, cut[idx]);
    dist_est[t] += d;
  }
  edge1[idx] = best_edge1;
  double diff2 = best_edge2 - edge2[idx];
  for (int i = 0; i < (int)turns.size(); ++i) {
    int t = turns[i];
    double d = diff2 * (getSum2(t) - getSum1(t, cut[idx]));
    dist_est[t] += d;
  }
  edge2[idx] = best_edge2;
  diff_sum -= (maxDiff - keep_diff);
}

void AnnealingMode3Vertical(int idx)
{
  AnnealingMode3Common<true>(idx);
}

void AnnealingMode3Horizontal(int idx)
{
  AnnealingMode3Common<false>(idx);
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
  Timer timer;
  timer.start();

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

  if (mode != 0) cout << timer.get_elapsed_time() << endl;
}
