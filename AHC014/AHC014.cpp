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

// U, L, D, R, UL, LD, DR, RU
const int dx[8] = { -1, 0, 1, 0, -1, 1, 1, -1 };
const int dy[8] = { 0, -1, 0, 1, -1, -1, 1, 1 };

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
    double maxScore;
    int ansSize;
    int ans[ANS_SIZE][4][2];
    int ansDelete[ANS_SIZE];
    int ansDeleteCount;
    int f[64][64];
    int line[64][64][8];
    int use[64][64];
    int cntH[64], cntW[64];

    // クリア関数
    void clear()
    {
      maxScore = 0;
      ansSize = 0;
      ansDeleteCount = 0;
      for (int i = 0; i < ANS_SIZE; ++i) ansDelete[i] = 0;
      for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
        f[i][j] = false;
        for (int k = 0; k < 8; ++k) line[i][j][k] = 0;
        use[i][j] = 0;
      }
      for (int i = 0; i < 64; ++i) {
        cntH[i] = 0;
        cntW[i] = 0;
      }
    }

    // 他のStateからコピー
    void copyFrom(const State& other)
    {
      maxScore = other.maxScore;
      ansSize = other.ansSize;
      ansDeleteCount = other.ansDeleteCount;

      for (int i = 0; i < ansSize; ++i) {
        for (int j = 0; j < 4; ++j) {
          ans[i][j][0] = other.ans[i][j][0];
          ans[i][j][1] = other.ans[i][j][1];
        }
        ansDelete[i] = other.ansDelete[i];
      }

      for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
        f[i][j] = other.f[i][j];
        for (int k = 0; k < 8; ++k) line[i][j][k] = other.line[i][j][k];
        use[i][j] = other.use[i][j];
      }

      for (int i = 0; i < 64; ++i) {
        cntH[i] = other.cntH[i];
        cntW[i] = other.cntW[i];
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
  int methodCount[20][2];

}  // namespace

void MethodCountReset()
{
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 2; ++j) { methodCount[i][j] = 0; }
  }
}

bool IsNGXY(int x, int y)
{
  if (x < 0 || N <= x || y < 0 || N <= y) return true;
  return false;
}

// State::calcScoreの実装
double State::calcScore() const
{
  double resd = 1000000.0 * N * N / M;

  int sum = 0;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (f[i][j]) {
        sum += W[i][j];
      }
    }
  }

  resd = resd * sum / S;

  return resd;
}

// 旧関数（互換性のため残す）
double CalcScore()
{
  return current_state.calcScore();
}

void RefleshAns()
{
  int tmpCount = 0;
  for (int i = 0; i < current_state.ansSize; ++i) {
    if (current_state.ansDelete[i]) {
      tmpCount++;
      current_state.ansDelete[i] = 0;
    }
    else {
      for (int j = 0; j < 4; ++j) for (int k = 0; k < 2; ++k) current_state.ans[i - tmpCount][j][k] = current_state.ans[i][j][k];
    }
  }
  if (tmpCount != current_state.ansDeleteCount) {
    cerr << "error" << endl;
    cerr << tmpCount << ' ' << current_state.ansDeleteCount << endl;
  }
  current_state.ansSize -= tmpCount;
  current_state.ansDeleteCount = 0;
}

// ローカルで複数ケース試すための全て消す関数
void AllClear_MultiCase()
{
  current_state.clear();
  best_state.clear();
  seed_state.clear();
  MethodCountReset();
}

// 初期状態作成（これを呼べばスタート位置に戻れることを想定、real_maxScore等は戻さない）
// State::initの実装
void State::init()
{
  clear();
  for (int i = 0; i < M; ++i) {
    f[X[i]][Y[i]] = true;
    use[X[i]][Y[i]] = 100;
    cntH[Y[i]]++;
    cntW[X[i]]++;
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
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
void Output(int mode, int problemNum)
{
  if (mode == 0) {
    cout << current_state.ansSize - current_state.ansDeleteCount << endl;
    for (int i = 0; i < current_state.ansSize; ++i) {
      if (current_state.ansDelete[i]) { continue; }
      for (int j = 0; j < 4; ++j) for (int k = 0; k < 2; ++k) { cout << current_state.ans[i][j][k] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    ofs << current_state.ansSize - current_state.ansDeleteCount << endl;
    for (int i = 0; i < current_state.ansSize; ++i) {
      if (current_state.ansDelete[i]) { continue; }
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
inline int FindNeighborStraight(int x, int y, int dir, int start, int end, int step, bool isVertical)
{
  for (int i = start; step > 0 ? i < end : i >= end; i += step) {
    int nx = isVertical ? i : x;
    int ny = isVertical ? y : i;

    if (current_state.f[nx][ny]) {
      int oppositeDir = (dir + 2) % 4; // 反対方向
      if (current_state.line[nx][ny][oppositeDir]) return -2;
      return abs(isVertical ? nx - x : ny - y);
    }
  }
  return -1;
}

// 斜め方向の隣接点を探索する共通関数
// 戻り値：距離（見つからない場合は-1、既に線がある場合は-2）
inline int FindNeighborDiagonal(int x, int y, int dir, int limit, int dx, int dy)
{
  // 斜め方向の反対方向の対応表
  // 4(UL) <-> 6(DR), 5(LD) <-> 7(RU)
  const int diagonalOpposite[8] = { 0, 0, 0, 0, 6, 7, 4, 5 };

  for (int i = 1; i < limit + 1; ++i) {
    int nx = x + i * dx;
    int ny = y + i * dy;

    if (current_state.f[nx][ny]) {
      int oppositeDir = diagonalOpposite[dir];
      if (current_state.line[nx][ny][oppositeDir]) return -2;
      return i;
    }
  }
  return -1;
}

// 各方向の1番近い点が使えるかどうか
// 入力 z：方向
// 戻り値：距離
inline int FindNeighborPoint(int x, int y, int z)
{
  // 上
  if (z == 0) {
    return FindNeighborStraight(x, y, 0, x - 1, -1, -1, true);
  }

  // 左
  if (z == 1) {
    return FindNeighborStraight(x, y, 1, y - 1, -1, -1, false);
  }

  // 下
  if (z == 2) {
    return FindNeighborStraight(x, y, 2, x + 1, N, 1, true);
  }

  // 右
  if (z == 3) {
    return FindNeighborStraight(x, y, 3, y + 1, N, 1, false);
  }

  // 左上
  if (z == 4) {
    int ma = min(x, y);
    return FindNeighborDiagonal(x, y, 4, ma, -1, -1);
  }

  // 左下
  if (z == 5) {
    int ma = min(N - 1 - x, y);
    return FindNeighborDiagonal(x, y, 5, ma, 1, -1);
  }

  // 右下
  if (z == 6) {
    int ma = min(N - 1 - x, N - 1 - y);
    return FindNeighborDiagonal(x, y, 6, ma, 1, 1);
  }

  // 右上
  if (z == 7) {
    int ma = min(x, N - 1 - y);
    return FindNeighborDiagonal(x, y, 7, ma, -1, 1);
  }

  cerr << "ERROR FindNeighborPoint" << endl;
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
  if (IsNGXY(xx, yy)) return false;

  // そこに点が存在しているか
  if (!current_state.f[xx][yy]) return false;

  // 辺が既に存在していないか
  if (current_state.line[xx][yy][params.line1] || current_state.line[xx][yy][params.line2]) return false;

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
  if (IsNGXY(xx, yy)) return false;

  // そこに点が存在しているか
  if (!current_state.f[xx][yy]) return false;

  // 辺が既に存在していないか
  if (current_state.line[xx][yy][params.line1] || current_state.line[xx][yy][params.line2]) return false;

  // 間に邪魔な頂点がないか (path1_limit, path2_limitは実際の値に置き換える必要がある)
  int limit1 = (params.path1_limit == -1) ? diff2 : diff1;
  int limit2 = (params.path2_limit == -1) ? diff2 : diff1;

  if (!CheckDiagonalPath(xx, yy, 1, limit1, params.path1_dx, params.path1_dy)) return false;
  if (!CheckDiagonalPath(xx, yy, 1, limit2, params.path2_dx, params.path2_dy)) return false;

  return true;
}

// 4点目の確認
inline bool CanMakeRectangle(int x, int y, int z, int diff1, int diff2)
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

  cerr << "ERROR CanMakeRectangle" << endl;
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

  line[x[0]][y[0]][settings.line1_1] = setValue;
  line[x[0]][y[0]][settings.line1_2] = setValue;
  line[x[1]][y[1]][settings.line2_1] = setValue;
  line[x[1]][y[1]][settings.line2_2] = setValue;
  line[x[2]][y[2]][settings.line3_1] = setValue;
  line[x[2]][y[2]][settings.line3_2] = setValue;
  line[x[3]][y[3]][settings.line4_1] = setValue;
  line[x[3]][y[3]][settings.line4_2] = setValue;
}

// State::addRectangleの実装
void State::addRectangle(int x[4], int y[4], int rectDir)
{
  ans[ansSize][0][0] = x[0]; ans[ansSize][0][1] = y[0];
  ans[ansSize][1][0] = x[1]; ans[ansSize][1][1] = y[1];
  ans[ansSize][2][0] = x[2]; ans[ansSize][2][1] = y[2];
  ans[ansSize][3][0] = x[3]; ans[ansSize][3][1] = y[3];
  ansDelete[ansSize] = 0;
  ansSize++;

  f[x[0]][y[0]] = 1;
  cntW[x[0]]++;
  cntH[y[0]]++;

  for (int i = 0; i < 4; ++i) {
    use[x[i]][y[i]]++;
  }

  setRectangleLines(x, y, rectDir, true);
}

// State::removeRectangleの実装
void State::removeRectangle(int x[4], int y[4], int rectDir)
{
  f[x[0]][y[0]] = 0;
  cntW[x[0]]--;
  cntH[y[0]]--;

  for (int i = 0; i < 4; ++i) {
    use[x[i]][y[i]]--;
  }

  setRectangleLines(x, y, rectDir, false);
}

/*
  メモ
  - ある1点を用いて描ける四角形は8種類
  - 使用する可能性のある頂点も8個
*/

// ランダムに1点が足せるかどうか
void Method1(double temperature)
{
  int x = Rand() % N;
  int y = Rand() % N;

  if (current_state.f[x][y]) { return; }

  methodCount[1][1]++;

  int u = -1;
  if (current_state.cntH[y] != 0) {
    u = FindNeighborPoint(x, y, 0);
    if (u == -2) { return; }
  }
  int l = -1;
  if (current_state.cntW[x] != 0) {
    l = FindNeighborPoint(x, y, 1);
    if (l == -2) { return; }
  }
  int d = -1;
  if (current_state.cntH[y] != 0) {
    d = FindNeighborPoint(x, y, 2);
    if (d == -2) { return; }
  }
  int r = -1;
  if (current_state.cntW[x] != 0) {
    r = FindNeighborPoint(x, y, 3);
    if (r == -2) { return; }
  }
  int ul = FindNeighborPoint(x, y, 4);
  if (ul == -2) { return; }
  int ld = FindNeighborPoint(x, y, 5);
  if (ld == -2) { return; }
  int dr = FindNeighborPoint(x, y, 6);
  if (dr == -2) { return; }
  int ru = FindNeighborPoint(x, y, 7);
  if (ru == -2) { return; }

  // 8種類の長方形
  int RectDir = -1;
  int xx = -1, yy = -1;
  int x1 = -1, y1 = -1;
  int x3 = -1, y3 = -1;

  // 上左
  if (RectDir == -1 && u != -1 && l != -1 && CanMakeRectangle(x, y, 0, u, l)) {
    RectDir = 0;
    xx = x - u;
    yy = y - l;
    x1 = x - u;
    y1 = y;
    x3 = x;
    y3 = y - l;
  }
  // 左下
  if (RectDir == -1 && l != -1 && d != -1 && CanMakeRectangle(x, y, 1, l, d)) {
    RectDir = 1;
    xx = x + d;
    yy = y - l;
    x1 = x;
    y1 = y - l;
    x3 = x + d;
    y3 = y;
  }
  // 下右
  if (RectDir == -1 && d != -1 && r != -1 && CanMakeRectangle(x, y, 2, d, r)) {
    RectDir = 2;
    xx = x + d;
    yy = y + r;
    x1 = x + d;
    y1 = y;
    x3 = x;
    y3 = y + r;
  }
  // 右上
  if (RectDir == -1 && r != -1 && u != -1 && CanMakeRectangle(x, y, 3, r, u)) {
    RectDir = 3;
    xx = x - u;
    yy = y + r;
    x1 = x;
    y1 = y + r;
    x3 = x - u;
    y3 = y;
  }

  // 左上・左下
  if (RectDir == -1 && ul != -1 && ld != -1 &&
    CanMakeRectangle(x, y, 4, ul, ld)) {
    RectDir = 4;
    xx = x - ul + ld;
    yy = y - ul - ld;
    x1 = x - ul;
    y1 = y - ul;
    x3 = x + ld;
    y3 = y - ld;
  }
  // 左下・右下
  if (RectDir == -1 && ld != -1 && dr != -1 &&
    CanMakeRectangle(x, y, 5, ld, dr)) {
    RectDir = 5;
    xx = x + ld + dr;
    yy = y - ld + dr;
    x1 = x + ld;
    y1 = y - ld;
    x3 = x + dr;
    y3 = y + dr;
  }
  // 右下・右上
  if (RectDir == -1 && dr != -1 && ru != -1 &&
    CanMakeRectangle(x, y, 6, dr, ru)) {
    RectDir = 6;
    xx = x + dr - ru;
    yy = y + dr + ru;
    x1 = x + dr;
    y1 = y + dr;
    x3 = x - ru;
    y3 = y + ru;
  }
  // 右上・左上
  if (RectDir == -1 && ru != -1 && ul != -1 &&
    CanMakeRectangle(x, y, 7, ru, ul)) {
    RectDir = 7;
    xx = x - ru - ul;
    yy = y + ru - ul;
    x1 = x - ru;
    y1 = y + ru;
    x3 = x - ul;
    y3 = y - ul;
  }

  if (RectDir == -1) { return; }

  double diffScore = 1000000.0 * N * N / M * W[x][y] / S;

  double prob = exp(diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[1][0]++;

    current_state.maxScore += diffScore;

    // 長方形を追加
    int rectX[4] = { x, x1, xx, x3 };
    int rectY[4] = { y, y1, yy, y3 };
    current_state.addRectangle(rectX, rectY, RectDir);

    if (current_state.maxScore > best_state.maxScore) {
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
void Method2(double temperature)
{
  if (current_state.ansSize == 0) { return; }
  int ite = Rand() % current_state.ansSize;
  if (current_state.ansDelete[ite]) { return; }
  if (current_state.use[current_state.ans[ite][0][0]][current_state.ans[ite][0][1]] > 1) { return; }

  methodCount[2][1]++;

  int x[4], y[4];
  for (int i = 0; i < 4; ++i) {
    x[i] = current_state.ans[ite][i][0];
    y[i] = current_state.ans[ite][i][1];
  }

  double diffScore = -1000000.0 * N * N / M * W[x[0]][y[0]] / S;

  double prob = exp(diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[2][0]++;

    current_state.maxScore += diffScore;

    // 長方形を削除
    int RectDir = GetDir(x[0], y[0], x[1], y[1]);
    current_state.removeRectangle(x, y, RectDir);

    current_state.ansDelete[ite] = 1;
    current_state.ansDeleteCount++;
  }
}

int Solve(int mode, int problemNum = 0)
{
  start_timer();

  // 初期状態作成
  current_state.init();

  // 愚直解作成
  current_state.maxScore = CalcScore();
  best_state.copyFrom(current_state);
  seed_state.copyFrom(current_state);

  // シード作り
  int seedCount = 20;  // 0にするとシード作成を行わない
  for (int tei = 0; tei < seedCount; ++tei) {
    start_timer();

    // 初期状態に戻す
    current_state.init();
    current_state.maxScore = CalcScore();

    // 焼きなまし
    double nowTime = get_elapsed_time();

    double TL = 4.2 / seedCount;
    double nowProgress = nowTime / TL;
    double startTemperature = 200048;
    double endTemperature = 0;
    int loop = 0;
    int rollbackCount = 0;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        nowTime = get_elapsed_time();
        nowProgress = nowTime / TL;
      }
      if (nowProgress > 1.0) {
        break;
      }

      // 現在のスコアが悪いときは元に戻す
      if (current_state.maxScore * 1.2 < best_state.maxScore) {
        current_state.copyFrom(best_state);
        rollbackCount++;
      }

      if (current_state.ansDeleteCount >= 10000) {
        RefleshAns();
      }

      double temperature = startTemperature + (endTemperature - startTemperature) * nowProgress;

      if (Rand() % 2 == 0) {
        Method1(temperature);
      }
      else {
        Method2(temperature);
      }
    }

    cerr << "seed loop = " << loop << ", rollbackCount = " << rollbackCount << endl;

    // スコアが良ければシードを更新
    current_state.copyFrom(best_state);
    if (current_state.maxScore > seed_state.maxScore) {
      seed_state.copyFrom(current_state);
    }
  }

  // シードから戻す
  current_state.copyFrom(seed_state);
  best_state.copyFrom(current_state);

  // 焼きなまし
  start_timer();
  double nowTime = get_elapsed_time();
  double TL = 0.5;
  double nowProgress = nowTime / TL;
  double startTemperature = 20048;
  double endTemperature = 0;
  int loop = 0;
  int rollbackCount = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      nowTime = get_elapsed_time();
      nowProgress = nowTime / TL;
    }
    if (nowProgress > 1.0) {
      break;
    }

    // 現在のスコアが悪いときは元に戻す
    if (current_state.maxScore * 1.2 < best_state.maxScore) {
      current_state.copyFrom(best_state);
      rollbackCount++;
    }

    if (current_state.ansDeleteCount >= 10000) {
      RefleshAns();
    }

    double temperature = startTemperature + (endTemperature - startTemperature) * nowProgress;

    if (Rand() % 2 == 0) {
      Method1(temperature);
    }
    else {
      Method2(temperature);
    }
  }

  cerr << "main loop = " << loop << ", rollbackCount = " << rollbackCount << endl;

  // 一番スコアの良い解
  current_state.copyFrom(best_state);

  RefleshAns();

  CalcScore();

  // デバッグログ
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", N = " << N << endl;
    cout << "ansSize = " << current_state.ansSize << ", ansDeleteCount = " << current_state.ansDeleteCount << endl;
    cout << "maxScore = " << current_state.maxScore << endl;
    cout << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    for (int i = 1; i < 5; ++i) {
      cout << "Method" << i << " = " << methodCount[i][0] << " / " << methodCount[i][1] << endl;
    }
    cout << endl;
  }

  cerr << loop << endl;
  return current_state.maxScore;
}

int SolveOuter(int mode, int problemNum = 0)
{
  // 入力受け取り
  Input(problemNum);

  int score = Solve(mode, problemNum);

  // 解答の出力
  Output(mode, problemNum);

  return score;
}

int main()
{
  clock_t mainStart, mainEnd;
  mainStart = clock();
  mainEnd = clock();

  int mode = 2;

  // 提出用
  if (mode == 0) {
    for (int i = 0; i < 1; ++i) {
      SolveOuter(mode, 3);
      AllClear_MultiCase();
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
      AllClear_MultiCase();
    }
    cout << "scoreSum = " << scoreSum << endl;
  }

  mainEnd = clock();
  cerr << (double)(mainEnd - mainStart) / CLOCKS_PER_SEC;
  return 0;
}
