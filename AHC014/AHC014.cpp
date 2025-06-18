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

  struct State
  {
    double maxScore;
    int ansSize;
    int ans[ANS_SIZE][4][2];
    int ansDelete[ANS_SIZE];
    int ansDeleteCount;
    int f[64][64];
    int line[64][64][8];
    int use[64][64];
    int cntH[64], cntW[64];
  };

  State current_state;
  State seed_state;
  State best_state;

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

// スコア計算
double CalcScore()
{
  double resd = 1000000.0 * N * N / M;

  int sum = 0;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (current_state.f[i][j]) {
        sum += W[i][j];
      }
    }
  }

  resd = resd * sum / S;

  return resd;
}

void NormalClear()
{
  current_state.maxScore = 0;
  current_state.ansSize = 0;
  current_state.ansDeleteCount = 0;
  for (int i = 0; i < ANS_SIZE; ++i) current_state.ansDelete[i] = 0;
  for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
    current_state.f[i][j] = false;
    for (int k = 0; k < 8; ++k) current_state.line[i][j][k] = 0;
    current_state.use[i][j] = 0;
  }
  for (int i = 0; i < 64; ++i) {
    current_state.cntH[i] = 0;
    current_state.cntW[i] = 0;
  }
}

void RealClear()
{
  best_state.maxScore = 0;
  best_state.ansSize = 0;
  best_state.ansDeleteCount = 0;
  for (int i = 0; i < ANS_SIZE; ++i) best_state.ansDelete[i] = 0;
  for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
    best_state.f[i][j] = false;
    for (int k = 0; k < 8; ++k) best_state.line[i][j][k] = 0;
    best_state.use[i][j] = 0;
  }
  for (int i = 0; i < 64; ++i) {
    best_state.cntH[i] = 0;
    best_state.cntW[i] = 0;
  }
}

void SeedClear()
{
  seed_state.maxScore = 0;
  seed_state.ansSize = 0;
  seed_state.ansDeleteCount = 0;
  for (int i = 0; i < ANS_SIZE; ++i) seed_state.ansDelete[i] = 0;
  for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) {
    seed_state.f[i][j] = false;
    for (int k = 0; k < 8; ++k) seed_state.line[i][j][k] = 0;
    seed_state.use[i][j] = 0;
  }
  for (int i = 0; i < 64; ++i) {
    seed_state.cntH[i] = 0;
    seed_state.cntW[i] = 0;
  }
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
  NormalClear();
  RealClear();
  SeedClear();
  MethodCountReset();
}

// 初期状態作成（これを呼べばスタート位置に戻れることを想定、real_maxScore等は戻さない）
void Init()
{
  NormalClear();
  for (int i = 0; i < M; ++i) {
    current_state.f[X[i]][Y[i]] = true;
    current_state.use[X[i]][Y[i]] = 100;
    current_state.cntH[Y[i]]++;
    current_state.cntW[X[i]]++;
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

  Init();
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

void CopyToReal()
{
  best_state = current_state;
}

void CopyToSeed()
{
  seed_state = current_state;
}

void RollBackFromReal()
{
  current_state = best_state;
}

void RollBackFromSeed()
{
  current_state = seed_state;
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

// 4点目の確認
inline bool CanMakeRectangle(int x, int y, int z, int diff1, int diff2)
{
  // 必ず反時計回り

  // 上左
  if (z == 0) {
    int xx = x - diff1;
    int yy = y - diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][3] || current_state.line[xx][yy][2]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckStraightPath(xx, 0, yy + 1, y, false)) return false;
    if (!CheckStraightPath(0, yy, xx + 1, x, true)) return false;

    return true;
  }

  // 左下
  if (z == 1) {
    int xx = x + diff2;
    int yy = y - diff1;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][0] || current_state.line[xx][yy][3]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckStraightPath(0, yy, xx - 1, x, true)) return false;
    if (!CheckStraightPath(xx, 0, yy + 1, y, false)) return false;

    return true;
  }

  // 下右
  if (z == 2) {
    int xx = x + diff1;
    int yy = y + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][1] || current_state.line[xx][yy][0]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckStraightPath(xx, 0, yy - 1, y, false)) return false;
    if (!CheckStraightPath(0, yy, xx - 1, x, true)) return false;

    return true;
  }

  // 右上
  if (z == 3) {
    int xx = x - diff2;
    int yy = y + diff1;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][2] || current_state.line[xx][yy][1]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckStraightPath(0, yy, xx + 1, x, true)) return false;
    if (!CheckStraightPath(xx, 0, yy - 1, y, false)) return false;

    return true;
  }

  // 左上・左下
  if (z == 4) {
    int xx = x - diff1 + diff2;
    int yy = y - diff1 - diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][7] || current_state.line[xx][yy][6]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckDiagonalPath(xx, yy, 1, diff2, -1, 1)) return false;
    if (!CheckDiagonalPath(xx, yy, 1, diff1, 1, 1)) return false;

    return true;
  }

  // 左下・右下
  if (z == 5) {
    int xx = x + diff1 + diff2;
    int yy = y - diff1 + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][4] || current_state.line[xx][yy][7]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckDiagonalPath(xx, yy, 1, diff2, -1, -1)) return false;
    if (!CheckDiagonalPath(xx, yy, 1, diff1, -1, 1)) return false;

    return true;
  }

  // 右下・右上
  if (z == 6) {
    int xx = x + diff1 - diff2;
    int yy = y + diff1 + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][5] || current_state.line[xx][yy][4]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckDiagonalPath(xx, yy, 1, diff2, 1, -1)) return false;
    if (!CheckDiagonalPath(xx, yy, 1, diff1, -1, -1)) return false;

    return true;
  }

  // 右上・左上
  if (z == 7) {
    int xx = x - diff1 - diff2;
    int yy = y + diff1 - diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!current_state.f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (current_state.line[xx][yy][6] || current_state.line[xx][yy][5]) return false;

    // 間に邪魔な頂点がないか
    if (!CheckDiagonalPath(xx, yy, 1, diff2, 1, 1)) return false;
    if (!CheckDiagonalPath(xx, yy, 1, diff1, 1, -1)) return false;

    return true;
  }

  cerr << "ERROR CanMakeRectangle" << endl;
  return false;
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

    current_state.ans[current_state.ansSize][0][0] = x;
    current_state.ans[current_state.ansSize][0][1] = y;
    current_state.ans[current_state.ansSize][1][0] = x1;
    current_state.ans[current_state.ansSize][1][1] = y1;
    current_state.ans[current_state.ansSize][2][0] = xx;
    current_state.ans[current_state.ansSize][2][1] = yy;
    current_state.ans[current_state.ansSize][3][0] = x3;
    current_state.ans[current_state.ansSize][3][1] = y3;
    current_state.ansDelete[current_state.ansSize] = 0;
    current_state.ansSize++;
    current_state.f[x][y] = 1;
    current_state.cntW[x]++;
    current_state.cntH[y]++;
    current_state.use[x][y]++;
    current_state.use[x1][y1]++;
    current_state.use[xx][yy]++;
    current_state.use[x3][y3]++;

    if (RectDir == 0) {
      current_state.line[x][y][0] = true;
      current_state.line[x][y][1] = true;
      current_state.line[x1][y1][1] = true;
      current_state.line[x1][y1][2] = true;
      current_state.line[xx][yy][2] = true;
      current_state.line[xx][yy][3] = true;
      current_state.line[x3][y3][3] = true;
      current_state.line[x3][y3][0] = true;
    }
    else if (RectDir == 1) {
      current_state.line[x][y][1] = true;
      current_state.line[x][y][2] = true;
      current_state.line[x1][y1][2] = true;
      current_state.line[x1][y1][3] = true;
      current_state.line[xx][yy][3] = true;
      current_state.line[xx][yy][0] = true;
      current_state.line[x3][y3][0] = true;
      current_state.line[x3][y3][1] = true;
    }
    else if (RectDir == 2) {
      current_state.line[x][y][2] = true;
      current_state.line[x][y][3] = true;
      current_state.line[x1][y1][3] = true;
      current_state.line[x1][y1][0] = true;
      current_state.line[xx][yy][0] = true;
      current_state.line[xx][yy][1] = true;
      current_state.line[x3][y3][1] = true;
      current_state.line[x3][y3][2] = true;
    }
    else if (RectDir == 3) {
      current_state.line[x][y][3] = true;
      current_state.line[x][y][0] = true;
      current_state.line[x1][y1][0] = true;
      current_state.line[x1][y1][1] = true;
      current_state.line[xx][yy][1] = true;
      current_state.line[xx][yy][2] = true;
      current_state.line[x3][y3][2] = true;
      current_state.line[x3][y3][3] = true;
    }
    else if (RectDir == 4) {
      current_state.line[x][y][4] = true;
      current_state.line[x][y][5] = true;
      current_state.line[x1][y1][5] = true;
      current_state.line[x1][y1][6] = true;
      current_state.line[xx][yy][6] = true;
      current_state.line[xx][yy][7] = true;
      current_state.line[x3][y3][7] = true;
      current_state.line[x3][y3][4] = true;
    }
    else if (RectDir == 5) {
      current_state.line[x][y][5] = true;
      current_state.line[x][y][6] = true;
      current_state.line[x1][y1][6] = true;
      current_state.line[x1][y1][7] = true;
      current_state.line[xx][yy][7] = true;
      current_state.line[xx][yy][4] = true;
      current_state.line[x3][y3][4] = true;
      current_state.line[x3][y3][5] = true;
    }
    else if (RectDir == 6) {
      current_state.line[x][y][6] = true;
      current_state.line[x][y][7] = true;
      current_state.line[x1][y1][7] = true;
      current_state.line[x1][y1][4] = true;
      current_state.line[xx][yy][4] = true;
      current_state.line[xx][yy][5] = true;
      current_state.line[x3][y3][5] = true;
      current_state.line[x3][y3][6] = true;
    }
    else if (RectDir == 7) {
      current_state.line[x][y][7] = true;
      current_state.line[x][y][4] = true;
      current_state.line[x1][y1][4] = true;
      current_state.line[x1][y1][5] = true;
      current_state.line[xx][yy][5] = true;
      current_state.line[xx][yy][6] = true;
      current_state.line[x3][y3][6] = true;
      current_state.line[x3][y3][7] = true;
    }

    if (current_state.maxScore > best_state.maxScore) {
      // RefleshAns();
      CopyToReal();
    }
  }
  else {
    // 元に戻す
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

    current_state.f[x[0]][y[0]] = 0;
    current_state.cntW[x[0]]--;
    current_state.cntH[y[0]]--;
    for (int i = 0; i < 4; ++i) current_state.use[x[i]][y[i]]--;

    int RectDir = GetDir(x[0], y[0], x[1], y[1]);
    if (RectDir == 0) {
      current_state.line[x[0]][y[0]][0] = false;
      current_state.line[x[0]][y[0]][1] = false;
      current_state.line[x[1]][y[1]][1] = false;
      current_state.line[x[1]][y[1]][2] = false;
      current_state.line[x[2]][y[2]][2] = false;
      current_state.line[x[2]][y[2]][3] = false;
      current_state.line[x[3]][y[3]][3] = false;
      current_state.line[x[3]][y[3]][0] = false;
    }
    else if (RectDir == 1) {
      current_state.line[x[0]][y[0]][1] = false;
      current_state.line[x[0]][y[0]][2] = false;
      current_state.line[x[1]][y[1]][2] = false;
      current_state.line[x[1]][y[1]][3] = false;
      current_state.line[x[2]][y[2]][3] = false;
      current_state.line[x[2]][y[2]][0] = false;
      current_state.line[x[3]][y[3]][0] = false;
      current_state.line[x[3]][y[3]][1] = false;
    }
    else if (RectDir == 2) {
      current_state.line[x[0]][y[0]][2] = false;
      current_state.line[x[0]][y[0]][3] = false;
      current_state.line[x[1]][y[1]][3] = false;
      current_state.line[x[1]][y[1]][0] = false;
      current_state.line[x[2]][y[2]][0] = false;
      current_state.line[x[2]][y[2]][1] = false;
      current_state.line[x[3]][y[3]][1] = false;
      current_state.line[x[3]][y[3]][2] = false;
    }
    else if (RectDir == 3) {
      current_state.line[x[0]][y[0]][3] = false;
      current_state.line[x[0]][y[0]][0] = false;
      current_state.line[x[1]][y[1]][0] = false;
      current_state.line[x[1]][y[1]][1] = false;
      current_state.line[x[2]][y[2]][1] = false;
      current_state.line[x[2]][y[2]][2] = false;
      current_state.line[x[3]][y[3]][2] = false;
      current_state.line[x[3]][y[3]][3] = false;
    }
    else if (RectDir == 4) {
      current_state.line[x[0]][y[0]][4] = false;
      current_state.line[x[0]][y[0]][5] = false;
      current_state.line[x[1]][y[1]][5] = false;
      current_state.line[x[1]][y[1]][6] = false;
      current_state.line[x[2]][y[2]][6] = false;
      current_state.line[x[2]][y[2]][7] = false;
      current_state.line[x[3]][y[3]][7] = false;
      current_state.line[x[3]][y[3]][4] = false;
    }
    else if (RectDir == 5) {
      current_state.line[x[0]][y[0]][5] = false;
      current_state.line[x[0]][y[0]][6] = false;
      current_state.line[x[1]][y[1]][6] = false;
      current_state.line[x[1]][y[1]][7] = false;
      current_state.line[x[2]][y[2]][7] = false;
      current_state.line[x[2]][y[2]][4] = false;
      current_state.line[x[3]][y[3]][4] = false;
      current_state.line[x[3]][y[3]][5] = false;
    }
    else if (RectDir == 6) {
      current_state.line[x[0]][y[0]][6] = false;
      current_state.line[x[0]][y[0]][7] = false;
      current_state.line[x[1]][y[1]][7] = false;
      current_state.line[x[1]][y[1]][4] = false;
      current_state.line[x[2]][y[2]][4] = false;
      current_state.line[x[2]][y[2]][5] = false;
      current_state.line[x[3]][y[3]][5] = false;
      current_state.line[x[3]][y[3]][6] = false;
    }
    else if (RectDir == 7) {
      current_state.line[x[0]][y[0]][7] = false;
      current_state.line[x[0]][y[0]][4] = false;
      current_state.line[x[1]][y[1]][4] = false;
      current_state.line[x[1]][y[1]][5] = false;
      current_state.line[x[2]][y[2]][5] = false;
      current_state.line[x[2]][y[2]][6] = false;
      current_state.line[x[3]][y[3]][6] = false;
      current_state.line[x[3]][y[3]][7] = false;
    }

    current_state.ansDelete[ite] = 1;
    current_state.ansDeleteCount++;
  }
  else {
    // 元に戻す
  }
}

int Solve(int mode, int problemNum = 0)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // 初期状態作成
  Init();

  // 愚直解作成
  current_state.maxScore = CalcScore();
  CopyToReal();
  CopyToSeed();

  // シード作り
  int seedCount = 20;  // 0にするとシード作成を行わない
  for (int tei = 0; tei < seedCount; ++tei) {
    startTime = clock();

    // 初期状態に戻す
    Init();
    current_state.maxScore = CalcScore();

    // 焼きなまし
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;

    double TL = 4.2 / seedCount;
    double nowProgress = nowTime / TL;
    double startTemperature = 200048;
    double endTemperature = 0;
    int loop = 0;
    int rollbackCount = 0;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
      }
      if (nowProgress > 1.0) { break; }

      // 現在のスコアが悪いときは元に戻す
      if (current_state.maxScore * 1.2 < best_state.maxScore) {
        RollBackFromReal();
        rollbackCount++;
      }

      if (current_state.ansDeleteCount >= 10000) {
        RefleshAns();
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      // メソッド選択
      int me = 1;
      if (Rand() % 2 == 0) {
        me = 2;
      }

      // 各メソッド処理

      if (me == 1) {
        Method1(temperature);
      }

      if (me == 2) {
        Method2(temperature);
      }
    }  // while文ここまで（シード作成）

    cerr << "seed loop = " << loop << ", rollbackCount = " << rollbackCount << endl;

    // スコアが良ければシードを更新
    RollBackFromReal();
    if (current_state.maxScore > seed_state.maxScore) {
      CopyToSeed();
    }

    // ここで消去するものがあれば消去する
  }

  // シードから戻す
  RollBackFromSeed();
  CopyToReal();

  // 焼きなまし
  startTime = clock();
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  double TL = 0.5;
  double nowProgress = nowTime / TL;
  double startTemperature = 20048;
  double endTemperature = 0;
  int loop = 0;
  int rollbackCount = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      endTime = clock();
      nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
    }
    if (nowProgress > 1.0) { break; }

    // 現在のスコアが悪いときは元に戻す
    if (current_state.maxScore * 1.2 < best_state.maxScore) {
      RollBackFromReal();
      rollbackCount++;
    }

    if (current_state.ansDeleteCount >= 10000) {
      RefleshAns();
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    // メソッド選択
    int me = 1;
    if (Rand() % 2 == 0) {
      me = 2;
    }

    // 各メソッド処理

    if (me == 1) {
      Method1(temperature);
    }

    if (me == 2) {
      Method2(temperature);
    }
  }  // while文ここまで（メインループ）

  cerr << "main loop = " << loop << ", rollbackCount = " << rollbackCount << endl;

  // 一番スコアの良い解
  RollBackFromReal();

  RefleshAns();

  CalcScore();

  // デバッグログ
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", N = " << N << endl;
    cout << "ansSize = " << current_state.ansSize << ", ansDeleteCount = " << current_state.ansDeleteCount
      << endl;
    cout << "maxScore = " << current_state.maxScore << endl;
    cout << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    for (int i = 1; i < 5; ++i) {
      cout << "Method" << i << " = " << methodCount[i][0] << " / "
        << methodCount[i][1] << endl;
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
