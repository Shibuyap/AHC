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

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

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
Timer globalTimer;

// U, L, D, R, UL, LD, DR, RU
const int dx[8] = { -1, 0, 1, 0, -1, 1, 1, -1 };
const int dy[8] = { 0, -1, 0, 1, -1, -1, 1, 1 };

namespace /* 乱数ライブラリ */
{
  static uint32_t rand32()
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


  static double rand_01()
  {
    return (rand32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* 変数 */
{
  // 入力用変数
  int N, M;
  int X[1000], Y[1000];
  int S;
  int W[64][64];

  // 解答用変数
  const int ANS_SIZE = 20000;

  class State
  {
  public:
    double max_score;
    int ans_count;
    int ans[ANS_SIZE][4][2];
    int deleted[ANS_SIZE];
    int del_count;
    int f[64][64];
    int edge[64][64][8];
    int deg[64][64];
    int rowCnt[64], colCnt[64];

    // クリア関数
    void clear()
    {
      max_score = 0;
      ans_count = 0;
      del_count = 0;
      for (int i = 0; i < ANS_SIZE; ++i) deleted[i] = 0;
      for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
        f[i][j] = false;
        for (int k = 0; k < 8; ++k) edge[i][j][k] = 0;
        deg[i][j] = 0;
      }
      for (int i = 0; i < 64; ++i) {
        rowCnt[i] = 0;
        colCnt[i] = 0;
      }
    }

    // 他のStateからコピー
    void copyFrom(const State& other)
    {
      max_score = other.max_score;
      ans_count = other.ans_count;
      del_count = other.del_count;

      for (int i = 0; i < ans_count; ++i) {
        for (int j = 0; j < 4; ++j) {
          ans[i][j][0] = other.ans[i][j][0];
          ans[i][j][1] = other.ans[i][j][1];
        }
        deleted[i] = other.deleted[i];
      }

      for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
        f[i][j] = other.f[i][j];
        for (int k = 0; k < 8; ++k) edge[i][j][k] = other.edge[i][j][k];
        deg[i][j] = other.deg[i][j];
      }

      for (int i = 0; i < 64; ++i) {
        rowCnt[i] = other.rowCnt[i];
        colCnt[i] = other.colCnt[i];
      }
    }

    // 長方形の辺を設定/解除するメソッド
    void setRectangleLines(int x[4], int y[4], int rectDir, bool setValue);

    // スコア計算メソッド
    double calcScore() const;

    // 初期化メソッド
    void init();

    // 長方形を追加するメソッド
    void addRectangle(int x[4], int y[4], int rectDir);

    // 長方形を削除するメソッド
    void removeRectangle(int x[4], int y[4], int rectDir);
  };

  State current_state;  // 現在の状態
  State seed_state;     // シード探索用の一時保存
  State best_state;     // 最良解の保存

  // その他
  int stats[20][2];

}  // namespace

void reset_stats()
{
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 2; ++j) {
      stats[i][j] = 0;
    }
  }
}

bool outOfBounds(int x, int y)
{
  if (x < 0 || N <= x || y < 0 || N <= y) {
    return true;
  }
  return false;
}

// State::calcScoreの実装
double State::calcScore() const
{
  double score = 1000000.0 * N * N / M;

  int sum = 0;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (f[i][j]) {
        sum += W[i][j];
      }
    }
  }

  score = score * sum / S;

  return score;
}

// 旧関数（互換性のため残す）
double calc_score()
{
  return current_state.calcScore();
}

void compact_answers()
{
  int tmpCount = 0;
  for (int i = 0; i < current_state.ans_count; ++i) {
    if (current_state.deleted[i]) {
      tmpCount++;
      current_state.deleted[i] = 0;
    }
    else {
      for (int j = 0; j < 4; ++j) for (int k = 0; k < 2; ++k) current_state.ans[i - tmpCount][j][k] = current_state.ans[i][j][k];
    }
  }
  if (tmpCount != current_state.del_count) {
    cerr << "error" << endl;
    cerr << tmpCount << ' ' << current_state.del_count << endl;
  }
  current_state.ans_count -= tmpCount;
  current_state.del_count = 0;
}

// ローカルで複数ケース試すための全て消す関数
void clear_all_multi_case()
{
  current_state.clear();
  best_state.clear();
  seed_state.clear();
  reset_stats();
}

// 初期状態作成（これを呼べばスタート位置に戻れることを想定、real_max_score等は戻さない）
// State::initの実装
void State::init()
{
  clear();
  for (int i = 0; i < M; ++i) {
    f[X[i]][Y[i]] = true;
    deg[X[i]][Y[i]] = 100;
    rowCnt[Y[i]]++;
    colCnt[X[i]]++;
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> N >> M;
    for (int i = 0; i < M; ++i) { cin >> X[i] >> Y[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N >> M;
    for (int i = 0; i < M; ++i) { ifs >> X[i] >> Y[i]; }
  }

  S = 0;
  int c = (N - 1) / 2;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      W[i][j] = (i - c) * (i - c) + (j - c) * (j - c) + 1;
      S += W[i][j];
    }
  }

  current_state.init();
}

// 解答出力
void output_data(int mode, int case_num)
{
  if (mode == 0) {
    cout << current_state.ans_count - current_state.del_count << endl;
    for (int i = 0; i < current_state.ans_count; ++i) {
      if (current_state.deleted[i]) { continue; }
      for (int j = 0; j < 4; ++j) for (int k = 0; k < 2; ++k) { cout << current_state.ans[i][j][k] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << current_state.ans_count - current_state.del_count << endl;
    for (int i = 0; i < current_state.ans_count; ++i) {
      if (current_state.deleted[i]) { continue; }
      for (int j = 0; j < 4; ++j) for (int k = 0; k < 2; ++k) { ofs << current_state.ans[i][j][k] << ' '; }
      ofs << endl;
    }

    ofs.close();
  }
}

/*
  8方向の順序のルール
  0 : 上 : U
  1 : 左 : L
  2 : 下 : D
  3 : 右 : R
  4 : 左上 : UL
  5 : 左下 : LD
  6 : 右下 : DR
  7 : 右上 : RU
*/

// 直線方向の隣接点を探索する共通関数
// 戻り値：距離（見つからない場合は-1、既に線がある場合は-2）
inline int findStr(int x, int y, int dir, int start, int end, int step, bool isVertical)
{
  for (int i = start; step > 0 ? i < end : i >= end; i += step) {
    int nx = isVertical ? i : x;
    int ny = isVertical ? y : i;

    if (current_state.f[nx][ny]) {
      int oppositeDir = (dir + 2) % 4; // 反対方向
      if (current_state.edge[nx][ny][oppositeDir]) return -2;
      return abs(isVertical ? nx - x : ny - y);
    }
  }
  return -1;
}

// 斜め方向の隣接点を探索する共通関数
// 戻り値：距離（見つからない場合は-1、既に線がある場合は-2）
inline int findDiag(int x, int y, int dir, int limit, int dx, int dy)
{
  // 斜め方向の反対方向の対応表
  // 4(UL) <-> 6(DR), 5(LD) <-> 7(RU)
  const int diagonalOpposite[8] = { 0, 0, 0, 0, 6, 7, 4, 5 };

  for (int i = 1; i < limit + 1; ++i) {
    int nx = x + i * dx;
    int ny = y + i * dy;

    if (current_state.f[nx][ny]) {
      int oppositeDir = diagonalOpposite[dir];
      if (current_state.edge[nx][ny][oppositeDir]) return -2;
      return i;
    }
  }
  return -1;
}

// 各方向の1番近い点が使えるかどうか
// 入力 z：方向
// 戻り値：距離
inline int find(int x, int y, int z)
{
  // 上
  if (z == 0) {
    return findStr(x, y, 0, x - 1, -1, -1, true);
  }

  // 左
  if (z == 1) {
    return findStr(x, y, 1, y - 1, -1, -1, false);
  }

  // 下
  if (z == 2) {
    return findStr(x, y, 2, x + 1, N, 1, true);
  }

  // 右
  if (z == 3) {
    return findStr(x, y, 3, y + 1, N, 1, false);
  }

  // 左上
  if (z == 4) {
    int ma = min(x, y);
    return findDiag(x, y, 4, ma, -1, -1);
  }

  // 左下
  if (z == 5) {
    int ma = min(N - 1 - x, y);
    return findDiag(x, y, 5, ma, 1, -1);
  }

  // 右下
  if (z == 6) {
    int ma = min(N - 1 - x, N - 1 - y);
    return findDiag(x, y, 6, ma, 1, 1);
  }

  // 右上
  if (z == 7) {
    int ma = min(x, N - 1 - y);
    return findDiag(x, y, 7, ma, -1, 1);
  }

  cerr << "ERROR find" << endl;
  return -2;
}

// 長方形の対角線上に邪魔な頂点がないかチェックする共通関数
inline bool CheckDiagonalPath(int xx, int yy, int start, int limit, int dx, int dy)
{
  for (int i = start; i < limit; ++i) {
    if (current_state.f[xx + i * dx][yy + i * dy]) return false;
  }
  return true;
}

// 長方形の辺上に邪魔な頂点がないかチェックする共通関数
inline bool CheckStraightPath(int x, int y, int start, int end, bool isVertical)
{
  if (isVertical) {
    for (int i = start; i != end; i += (start < end ? 1 : -1)) {
      if (current_state.f[i][y]) return false;
    }
  }
  else {
    for (int j = start; j != end; j += (start < end ? 1 : -1)) {
      if (current_state.f[x][j]) return false;
    }
  }
  return true;
}

// 直線方向の長方形チェック用構造体
struct RectCheckParamsStraight
{
  int dx1, dy1, dx2, dy2;  // 座標計算用
  int line1, line2;         // チェックする辺の方向
  int path1_start, path1_end, path1_vertical;
  int path2_start, path2_end, path2_vertical;
};

// 斜め方向の長方形チェック用構造体
struct RectCheckParamsDiagonal
{
  int dx1, dy1, dx2, dy2;  // 座標計算用
  int line1, line2;         // チェックする辺の方向
  int path1_dx, path1_dy, path1_limit;
  int path2_dx, path2_dy, path2_limit;
};

// 直線方向の長方形チェック共通関数
inline bool CheckRectangleStraight(int x, int y, int diff1, int diff2, const RectCheckParamsStraight& params)
{
  int xx = x + params.dx1 * diff1 + params.dx2 * diff2;
  int yy = y + params.dy1 * diff1 + params.dy2 * diff2;

  // その点がグリッド内か
  if (outOfBounds(xx, yy)) return false;

  // そこに点が存在しているか
  if (!current_state.f[xx][yy]) return false;

  // 辺が既に存在していないか
  if (current_state.edge[xx][yy][params.line1] || current_state.edge[xx][yy][params.line2]) return false;

  // 間に邪魔な頂点がないか
  int p1_x = params.path1_vertical ? 0 : xx;
  int p1_y = params.path1_vertical ? yy : 0;
  int p1_start = params.path1_vertical ? (xx + params.path1_start) : (yy + params.path1_start);
  int p1_end = params.path1_vertical ? (x + params.path1_end) : (y + params.path1_end);

  int p2_x = params.path2_vertical ? 0 : xx;
  int p2_y = params.path2_vertical ? yy : 0;
  int p2_start = params.path2_vertical ? (xx + params.path2_start) : (yy + params.path2_start);
  int p2_end = params.path2_vertical ? (x + params.path2_end) : (y + params.path2_end);

  if (!CheckStraightPath(p1_x, p1_y, p1_start, p1_end, params.path1_vertical)) return false;
  if (!CheckStraightPath(p2_x, p2_y, p2_start, p2_end, params.path2_vertical)) return false;

  return true;
}

// 斜め方向の長方形チェック共通関数
inline bool CheckRectangleDiagonal(int x, int y, int diff1, int diff2, const RectCheckParamsDiagonal& params)
{
  int xx = x + params.dx1 * diff1 + params.dx2 * diff2;
  int yy = y + params.dy1 * diff1 + params.dy2 * diff2;

  // その点がグリッド内か
  if (outOfBounds(xx, yy)) return false;

  // そこに点が存在しているか
  if (!current_state.f[xx][yy]) return false;

  // 辺が既に存在していないか
  if (current_state.edge[xx][yy][params.line1] || current_state.edge[xx][yy][params.line2]) return false;

  // 間に邪魔な頂点がないか (path1_limit, path2_limitは実際の値に置き換える必要がある)
  int limit1 = (params.path1_limit == -1) ? diff2 : diff1;
  int limit2 = (params.path2_limit == -1) ? diff2 : diff1;

  if (!CheckDiagonalPath(xx, yy, 1, limit1, params.path1_dx, params.path1_dy)) return false;
  if (!CheckDiagonalPath(xx, yy, 1, limit2, params.path2_dx, params.path2_dy)) return false;

  return true;
}

// 4点目の確認
inline bool canMake(int x, int y, int z, int diff1, int diff2)
{
  // 必ず反時計回り

  // 直線方向のパラメータ設定
  static const RectCheckParamsStraight straightParams[4] = {
    // 上左 (z=0)
    { -1, 0, 0, -1, 3, 2, 1, 0, false, 1, 0, true },
    // 左下 (z=1)
    { 0, -1, 1, 0, 0, 3, -1, 0, true, 1, 0, false },
    // 下右 (z=2)
    { 1, 0, 0, 1, 1, 0, -1, 0, false, -1, 0, true },
    // 右上 (z=3)
    { 0, 1, -1, 0, 2, 1, 1, 0, true, -1, 0, false }
  };

  if (z >= 0 && z <= 3) {
    return CheckRectangleStraight(x, y, diff1, diff2, straightParams[z]);
  }

  // 斜め方向のパラメータ設定
  static const RectCheckParamsDiagonal diagonalParams[4] = {
    // 左上・左下 (z=4) path1_limit=-1 は diff2を使う、path2_limit=1 は diff1を使う
    { -1, -1, 1, -1, 7, 6, -1, 1, -1, 1, 1, 1 },
    // 左下・右下 (z=5)
    { 1, -1, 1, 1, 4, 7, -1, -1, -1, -1, 1, 1 },
    // 右下・右上 (z=6)
    { 1, 1, -1, 1, 5, 4, 1, -1, -1, -1, -1, 1 },
    // 右上・左上 (z=7)
    { -1, 1, -1, -1, 6, 5, 1, 1, -1, 1, -1, 1 }
  };

  if (z >= 4 && z <= 7) {
    return CheckRectangleDiagonal(x, y, diff1, diff2, diagonalParams[z - 4]);
  }

  cerr << "ERROR canMake" << endl;
  return false;
}

// 長方形の辺の設定用構造体
struct RectLineSettings
{
  int line1_1, line1_2;  // 1番目の点の辺
  int line2_1, line2_2;  // 2番目の点の辺
  int line3_1, line3_2;  // 3番目の点の辺
  int line4_1, line4_2;  // 4番目の点の辺
};

// 各方向の辺の設定
static const RectLineSettings rectLineSettings[8] = {
  // 上左 (RectDir = 0)
  { 0, 1, 1, 2, 2, 3, 3, 0 },
  // 左下 (RectDir = 1)
  { 1, 2, 2, 3, 3, 0, 0, 1 },
  // 下右 (RectDir = 2)
  { 2, 3, 3, 0, 0, 1, 1, 2 },
  // 右上 (RectDir = 3)
  { 3, 0, 0, 1, 1, 2, 2, 3 },
  // 左上・左下 (RectDir = 4)
  { 4, 5, 5, 6, 6, 7, 7, 4 },
  // 左下・右下 (RectDir = 5)
  { 5, 6, 6, 7, 7, 4, 4, 5 },
  // 右下・右上 (RectDir = 6)
  { 6, 7, 7, 4, 4, 5, 5, 6 },
  // 右上・左上 (RectDir = 7)
  { 7, 4, 4, 5, 5, 6, 6, 7 }
};

// State::setRectangleLinesの実装
void State::setRectangleLines(int x[4], int y[4], int rectDir, bool setValue)
{
  const RectLineSettings& settings = rectLineSettings[rectDir];

  edge[x[0]][y[0]][settings.line1_1] = setValue;
  edge[x[0]][y[0]][settings.line1_2] = setValue;
  edge[x[1]][y[1]][settings.line2_1] = setValue;
  edge[x[1]][y[1]][settings.line2_2] = setValue;
  edge[x[2]][y[2]][settings.line3_1] = setValue;
  edge[x[2]][y[2]][settings.line3_2] = setValue;
  edge[x[3]][y[3]][settings.line4_1] = setValue;
  edge[x[3]][y[3]][settings.line4_2] = setValue;
}

// State::addRectangleの実装
void State::addRectangle(int x[4], int y[4], int rectDir)
{
  ans[ans_count][0][0] = x[0]; ans[ans_count][0][1] = y[0];
  ans[ans_count][1][0] = x[1]; ans[ans_count][1][1] = y[1];
  ans[ans_count][2][0] = x[2]; ans[ans_count][2][1] = y[2];
  ans[ans_count][3][0] = x[3]; ans[ans_count][3][1] = y[3];
  deleted[ans_count] = 0;
  ans_count++;

  f[x[0]][y[0]] = 1;
  colCnt[x[0]]++;
  rowCnt[y[0]]++;

  for (int i = 0; i < 4; ++i) {
    deg[x[i]][y[i]]++;
  }

  setRectangleLines(x, y, rectDir, true);
}

// State::removeRectangleの実装
void State::removeRectangle(int x[4], int y[4], int rectDir)
{
  f[x[0]][y[0]] = 0;
  colCnt[x[0]]--;
  rowCnt[y[0]]--;

  for (int i = 0; i < 4; ++i) {
    deg[x[i]][y[i]]--;
  }

  setRectangleLines(x, y, rectDir, false);
}

/*
  メモ
  - ある1点を用いて描ける四角形は8種類
  - 使用する可能性のある頂点も8個
*/

// ランダムに1点が足せるかどうか
void try_add(double temperature)
{
  int x = rand32() % N;
  int y = rand32() % N;

  if (current_state.f[x][y]) { return; }

  stats[1][1]++;

  int u = -1;
  if (current_state.rowCnt[y] != 0) {
    u = find(x, y, 0);
    if (u == -2) { return; }
  }
  int l = -1;
  if (current_state.colCnt[x] != 0) {
    l = find(x, y, 1);
    if (l == -2) { return; }
  }
  int d = -1;
  if (current_state.rowCnt[y] != 0) {
    d = find(x, y, 2);
    if (d == -2) { return; }
  }
  int r = -1;
  if (current_state.colCnt[x] != 0) {
    r = find(x, y, 3);
    if (r == -2) { return; }
  }
  int ul = find(x, y, 4);
  if (ul == -2) { return; }
  int ld = find(x, y, 5);
  if (ld == -2) { return; }
  int dr = find(x, y, 6);
  if (dr == -2) { return; }
  int ru = find(x, y, 7);
  if (ru == -2) { return; }

  // 8種類の長方形
  int dir = -1;
  int xx = -1, yy = -1;
  int x1 = -1, y1 = -1;
  int x3 = -1, y3 = -1;

  // 上左
  if (dir == -1 && u != -1 && l != -1 && canMake(x, y, 0, u, l)) {
    dir = 0;
    xx = x - u;
    yy = y - l;
    x1 = x - u;
    y1 = y;
    x3 = x;
    y3 = y - l;
  }
  // 左下
  if (dir == -1 && l != -1 && d != -1 && canMake(x, y, 1, l, d)) {
    dir = 1;
    xx = x + d;
    yy = y - l;
    x1 = x;
    y1 = y - l;
    x3 = x + d;
    y3 = y;
  }
  // 下右
  if (dir == -1 && d != -1 && r != -1 && canMake(x, y, 2, d, r)) {
    dir = 2;
    xx = x + d;
    yy = y + r;
    x1 = x + d;
    y1 = y;
    x3 = x;
    y3 = y + r;
  }
  // 右上
  if (dir == -1 && r != -1 && u != -1 && canMake(x, y, 3, r, u)) {
    dir = 3;
    xx = x - u;
    yy = y + r;
    x1 = x;
    y1 = y + r;
    x3 = x - u;
    y3 = y;
  }

  // 左上・左下
  if (dir == -1 && ul != -1 && ld != -1 &&
    canMake(x, y, 4, ul, ld)) {
    dir = 4;
    xx = x - ul + ld;
    yy = y - ul - ld;
    x1 = x - ul;
    y1 = y - ul;
    x3 = x + ld;
    y3 = y - ld;
  }
  // 左下・右下
  if (dir == -1 && ld != -1 && dr != -1 &&
    canMake(x, y, 5, ld, dr)) {
    dir = 5;
    xx = x + ld + dr;
    yy = y - ld + dr;
    x1 = x + ld;
    y1 = y - ld;
    x3 = x + dr;
    y3 = y + dr;
  }
  // 右下・右上
  if (dir == -1 && dr != -1 && ru != -1 &&
    canMake(x, y, 6, dr, ru)) {
    dir = 6;
    xx = x + dr - ru;
    yy = y + dr + ru;
    x1 = x + dr;
    y1 = y + dr;
    x3 = x - ru;
    y3 = y + ru;
  }
  // 右上・左上
  if (dir == -1 && ru != -1 && ul != -1 &&
    canMake(x, y, 7, ru, ul)) {
    dir = 7;
    xx = x - ru - ul;
    yy = y + ru - ul;
    x1 = x - ru;
    y1 = y + ru;
    x3 = x - ul;
    y3 = y - ul;
  }

  if (dir == -1) { return; }

  double diffScore = 1000000.0 * N * N / M * W[x][y] / S;

  double prob = exp(diffScore / temperature);
  if (prob > rand_01()) {
    stats[1][0]++;

    current_state.max_score += diffScore;

    // 長方形を追加
    int rectX[4] = { x, x1, xx, x3 };
    int rectY[4] = { y, y1, yy, y3 };
    current_state.addRectangle(rectX, rectY, dir);

    if (current_state.max_score > best_state.max_score) {
      best_state.copyFrom(current_state);
    }
  }
}

inline int GetDir(int x1, int y1, int x2, int y2)
{
  if (x2 < x1 && y2 == y1) return 0;
  if (x2 == x1 && y2 < y1) return 1;
  if (x2 > x1 && y2 == y1) return 2;
  if (x2 == x1 && y2 > y1) return 3;
  if (x2 < x1 && y2 < y1) return 4;
  if (x2 > x1 && y2 < y1) return 5;
  if (x2 > x1 && y2 > y1) return 6;
  if (x2 < x1 && y2 > y1) return 7;
  return -1;
}

// ランダムに1点選びほかに影響ないなら削除
void tryDel(double temperature)
{
  if (current_state.ans_count == 0) { return; }
  int ite = rand32() % current_state.ans_count;
  if (current_state.deleted[ite]) { return; }
  if (current_state.deg[current_state.ans[ite][0][0]][current_state.ans[ite][0][1]] > 1) { return; }

  stats[2][1]++;

  int x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = current_state.ans[ite][i][0];
    y[i] = current_state.ans[ite][i][1];
  }

  double diffScore = -1000000.0 * N * N / M * W[x[0]][y[0]] / S;

  double prob = exp(diffScore / temperature);
  if (prob > rand_01()) {
    stats[2][0]++;

    current_state.max_score += diffScore;

    // 長方形を削除
    int dir = GetDir(x[0], y[0], x[1], y[1]);
    current_state.removeRectangle(x, y, dir);

    current_state.deleted[ite] = 1;
    current_state.del_count++;
  }
}

// 焼きなまし処理の共通関数
// 戻り値: pair<loop, rollbackCount>
pair<int, int> SimulatedAnnealing(double timeLimit, double startTemp, double endTemp, const char* debugLabel = "")
{
  double nowTime = globalTimer.get_elapsed_time();
  double nowProgress = nowTime / timeLimit;
  int loop = 0;
  int rollbackCount = 0;

  while (true) {
    loop++;
    if (loop % 100 == 1) {
      nowTime = globalTimer.get_elapsed_time();
      nowProgress = nowTime / timeLimit;
    }
    if (nowProgress > 1.0) {
      break;
    }

    // 現在のスコアが悪いときは元に戻す
    if (current_state.max_score * 1.2 < best_state.max_score) {
      current_state.copyFrom(best_state);
      rollbackCount++;
    }

    if (current_state.del_count >= 10000) {
      compact_answers();
    }

    double temperature = startTemp + (endTemp - startTemp) * nowProgress;

    if (rand32() % 2 == 0) {
      try_add(temperature);
    }
    else {
      tryDel(temperature);
    }
  }

  if (debugLabel[0] != '\0') {
    cerr << debugLabel << " loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
  }

  return make_pair(loop, rollbackCount);
}

int Solve(int mode, int problemNum = 0)
{
  globalTimer.start();

  // 初期状態作成
  current_state.init();

  // 愚直解作成
  current_state.max_score = calc_score();
  best_state.copyFrom(current_state);
  seed_state.copyFrom(current_state);

  // シード作り
  int seedCount = 20;  // 0にするとシード作成を行わない
  for (int seed = 0; seed < seedCount; ++seed) {
    globalTimer.start();

    // 初期状態に戻す
    current_state.init();
    current_state.max_score = calc_score();

    // 焼きなまし
    SimulatedAnnealing(4.2 / seedCount, 200048, 0, "seed");

    // スコアが良ければシードを更新
    if (best_state.max_score > seed_state.max_score) {
      seed_state.copyFrom(best_state);
    }
  }

  // シードから戻す
  current_state.copyFrom(seed_state);
  best_state.copyFrom(current_state);

  // 焼きなまし
  globalTimer.start();
  pair<int, int> mainResult = SimulatedAnnealing(0.5, 20048, 0, "main");
  int loop = mainResult.first;
  int rollbackCount = mainResult.second;

  // 一番スコアの良い解
  current_state.copyFrom(best_state);

  compact_answers();

  calc_score();

  // デバッグログ
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", N = " << N << endl;
    cout << "ans_count = " << current_state.ans_count << ", del_count = " << current_state.del_count << endl;
    cout << "max_score = " << current_state.max_score << endl;
    cout << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    for (int i = 1; i < 5; ++i) {
      cout << "Method" << i << " = " << stats[i][0] << " / " << stats[i][1] << endl;
    }
    cout << endl;
  }

  cerr << loop << endl;
  return current_state.max_score;
}

int SolveOuter(int mode, int problemNum = 0)
{
  // 入力受け取り
  input_data(problemNum);

  int score = Solve(mode, problemNum);

  // 解答の出力
  output_data(mode, problemNum);

  return score;
}

int main()
{
  Timer mainTimer;
  mainTimer.start();

  int mode = 2;

  // 提出用
  if (mode == 0) {
    for (int i = 0; i < 1; ++i) {
      SolveOuter(mode, 3);
      clear_all_multi_case();
    }
  }
  // 1ケース試す
  else if (mode == 1) {
    SolveOuter(mode, 3);
  }
  // 複数ケース試す
  else if (mode == 2) {
    int scoreSum = 0;
    for (int i = 0; i < 10; ++i) {
      scoreSum += SolveOuter(mode, i);
      clear_all_multi_case();
    }
    cout << "scoreSum = " << scoreSum << endl;
  }

  cerr << mainTimer.get_elapsed_time() / CLOCKS_PER_SEC;
  return 0;
}
