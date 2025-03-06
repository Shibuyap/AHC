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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = s; i > t; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

/*
いろいろ
const int INF = 1001001001;
const char cc[4] = {'U', 'L', 'D', 'R'};

*/

// U, L, D, R, UL, LD, DR, RU
const int dx[8] = { -1, 0, 1, 0, -1, 1, 1, -1 };
const int dy[8] = { 0, -1, 0, 1, -1, -1, 1, 1 };

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
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
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }
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
  double maxScore;
  int ansSize;
  int ans[ANS_SIZE][4][2];
  int ansDelete[ANS_SIZE];
  int ansDeleteCount;
  int f[64][64];
  int line[64][64][8];
  int use[64][64];
  int cntH[64], cntW[64];

  double real_maxScore;
  int real_ansSize;
  int real_ans[ANS_SIZE][4][2];
  int real_ansDelete[ANS_SIZE];
  int real_ansDeleteCount;
  bool real_f[64][64];
  bool real_line[64][64][8];
  int real_use[64][64];
  int real_cntH[64], real_cntW[64];

  double seed_maxScore;
  int seed_ansSize;
  int seed_ans[ANS_SIZE][4][2];
  int seed_ansDelete[ANS_SIZE];
  int seed_ansDeleteCount;
  bool seed_f[64][64];
  bool seed_line[64][64][8];
  int seed_use[64][64];
  int seed_cntH[64], seed_cntW[64];

  // その他
  int methodCount[20][2];

}  // namespace

void MethodCountReset()
{
  rep(i, 20)
  {
    rep(j, 2) { methodCount[i][j] = 0; }
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

  rep(i, N)
  {
    rep(j, N)
    {
      if (f[i][j]) {
        sum += W[i][j];
      }
    }
  }

  resd = resd * sum / S;

  return resd;
}

void NormalClear()
{
  maxScore = 0;
  ansSize = 0;
  ansDeleteCount = 0;
  rep(i, ANS_SIZE) ansDelete[i] = 0;
  rep(i, 64) rep(j, 64)
  {
    f[i][j] = false;
    rep(k, 8) line[i][j][k] = 0;
    use[i][j] = 0;
  }
  rep(i, 64)
  {
    cntH[i] = 0;
    cntW[i] = 0;
  }
}

void RealClear()
{
  real_maxScore = 0;
  real_ansSize = 0;
  real_ansDeleteCount = 0;
  rep(i, ANS_SIZE) real_ansDelete[i] = 0;
  rep(i, 64) rep(j, 64)
  {
    real_f[i][j] = false;
    rep(k, 8) real_line[i][j][k] = 0;
    real_use[i][j] = 0;
  }
  rep(i, 64)
  {
    real_cntH[i] = 0;
    real_cntW[i] = 0;
  }
}

void SeedClear()
{
  seed_maxScore = 0;
  seed_ansSize = 0;
  seed_ansDeleteCount = 0;
  rep(i, ANS_SIZE) seed_ansDelete[i] = 0;
  rep(i, 64) rep(j, 64)
  {
    seed_f[i][j] = false;
    rep(k, 8) seed_line[i][j][k] = 0;
    seed_use[i][j] = 0;
  }
  rep(i, 64)
  {
    seed_cntH[i] = 0;
    seed_cntW[i] = 0;
  }
}

void RefleshAns()
{
  int tmpCount = 0;
  rep(i, ansSize)
  {
    if (ansDelete[i]) {
      tmpCount++;
      ansDelete[i] = 0;
    }
    else {
      rep(j, 4) rep(k, 2) ans[i - tmpCount][j][k] = ans[i][j][k];
    }
  }
  if (tmpCount != ansDeleteCount) {
    cerr << "error" << endl;
    cerr << tmpCount << ' ' << ansDeleteCount << endl;
  }
  ansSize -= tmpCount;
  ansDeleteCount = 0;
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
  rep(i, M)
  {
    f[X[i]][Y[i]] = true;
    use[X[i]][Y[i]] = 100;
    cntH[Y[i]]++;
    cntW[X[i]]++;
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> N >> M;
    rep(i, M) { cin >> X[i] >> Y[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N >> M;
    rep(i, M) { ifs >> X[i] >> Y[i]; }
  }

  S = 0;
  int c = (N - 1) / 2;
  rep(i, N)
  {
    rep(j, N)
    {
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
    cout << ansSize - ansDeleteCount << endl;
    rep(i, ansSize)
    {
      if (ansDelete[i]) continue;
      rep(j, 4) rep(k, 2) { cout << ans[i][j][k] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    ofs << ansSize - ansDeleteCount << endl;
    rep(i, ansSize)
    {
      if (ansDelete[i]) continue;
      rep(j, 4) rep(k, 2) { ofs << ans[i][j][k] << ' '; }
      ofs << endl;
    }

    ofs.close();
  }
}

void CopyToReal()
{
  real_maxScore = maxScore;
  real_ansSize = ansSize;
  real_ansDeleteCount = ansDeleteCount;
  rep(i, ansSize)
  {
    rep(j, 4) rep(k, 2) { real_ans[i][j][k] = ans[i][j][k]; }
    real_ansDelete[i] = ansDelete[i];
  }
  rep(i, N) rep(j, N)
  {
    real_f[i][j] = f[i][j];
    rep(k, 8) real_line[i][j][k] = line[i][j][k];
    real_use[i][j] = use[i][j];
  }
  rep(i, 64)
  {
    real_cntH[i] = cntH[i];
    real_cntW[i] = cntW[i];
  }
}

void CopyToSeed()
{
  seed_maxScore = maxScore;
  seed_ansSize = ansSize;
  seed_ansDeleteCount = ansDeleteCount;
  rep(i, ansSize)
  {
    rep(j, 4) rep(k, 2) { seed_ans[i][j][k] = ans[i][j][k]; }
    seed_ansDelete[i] = ansDelete[i];
  }
  rep(i, N) rep(j, N)
  {
    seed_f[i][j] = f[i][j];
    rep(k, 8) seed_line[i][j][k] = line[i][j][k];
    seed_use[i][j] = use[i][j];
  }
  rep(i, 64)
  {
    seed_cntH[i] = cntH[i];
    seed_cntW[i] = cntW[i];
  }
}

void RollBackFromReal()
{
  maxScore = real_maxScore;
  ansSize = real_ansSize;
  ansDeleteCount = real_ansDeleteCount;
  rep(i, ansSize)
  {
    rep(j, 4) rep(k, 2) { ans[i][j][k] = real_ans[i][j][k]; }
    ansDelete[i] = real_ansDelete[i];
  }
  rep(i, N) rep(j, N)
  {
    f[i][j] = real_f[i][j];
    rep(k, 8) line[i][j][k] = real_line[i][j][k];
    use[i][j] = real_use[i][j];
  }
  rep(i, 64)
  {
    cntH[i] = real_cntH[i];
    cntW[i] = real_cntW[i];
  }
}

void RollBackFromSeed()
{
  maxScore = seed_maxScore;
  ansSize = seed_ansSize;
  ansDeleteCount = seed_ansDeleteCount;
  rep(i, ansSize)
  {
    rep(j, 4) rep(k, 2) { ans[i][j][k] = seed_ans[i][j][k]; }
    ansDelete[i] = seed_ansDelete[i];
  }
  rep(i, N) rep(j, N)
  {
    f[i][j] = seed_f[i][j];
    rep(k, 8) line[i][j][k] = seed_line[i][j][k];
    use[i][j] = seed_use[i][j];
  }
  rep(i, 64)
  {
    cntH[i] = seed_cntH[i];
    cntW[i] = seed_cntW[i];
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

// 各方向の1番近い点が使えるかどうか
// 入力 z：方向
// 戻り値：距離
inline int FindNeighborPoint(int x, int y, int z)
{
  // 上
  if (z == 0) {
    drep(i, x)
    {
      if (f[i][y]) {
        if (line[i][y][2]) return -2;
        return x - i;
      }
    }
    return -1;
  }

  // 左
  if (z == 1) {
    drep(j, y)
    {
      if (f[x][j]) {
        if (line[x][j][3]) return -2;
        return y - j;
      }
    }
    return -1;
  }

  // 下
  if (z == 2) {
    srep(i, x + 1, N)
    {
      if (f[i][y]) {
        if (line[i][y][0]) return -2;
        return i - x;
      }
    }
    return -1;
  }

  // 右
  if (z == 3) {
    srep(j, y + 1, N)
    {
      if (f[x][j]) {
        if (line[x][j][1]) return -2;
        return j - y;
      }
    }
    return -1;
  }

  // 左上
  if (z == 4) {
    int ma = min(x, y);
    srep(i, 1, ma + 1)
    {
      if (f[x - i][y - i]) {
        if (line[x - i][y - i][6]) return -2;
        return i;
      }
    }
    return -1;
  }

  // 左下
  if (z == 5) {
    int ma = min(N - 1 - x, y);
    srep(i, 1, ma + 1)
    {
      if (f[x + i][y - i]) {
        if (line[x + i][y - i][7]) return -2;
        return i;
      }
    }
    return -1;
  }

  // 右下
  if (z == 6) {
    int ma = min(N - 1 - x, N - 1 - y);
    srep(i, 1, ma + 1)
    {
      if (f[x + i][y + i]) {
        if (line[x + i][y + i][4]) return -2;
        return i;
      }
    }
    return -1;
  }

  // 右上
  if (z == 7) {
    int ma = min(x, N - 1 - y);
    srep(i, 1, ma + 1)
    {
      if (f[x - i][y + i]) {
        if (line[x - i][y + i][5]) return -2;
        return i;
      }
    }
    return -1;
  }

  cerr << "ERROR FindNeighborPoint" << endl;
  return -2;
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
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][3] || line[xx][yy][2]) return false;

    // 間に邪魔な頂点がないか
    srep(j, yy + 1, y) if (f[xx][j]) return false;
    srep(i, xx + 1, x) if (f[i][yy]) return false;

    return true;
  }

  // 左下
  if (z == 1) {
    int xx = x + diff2;
    int yy = y - diff1;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][0] || line[xx][yy][3]) return false;

    // 間に邪魔な頂点がないか
    dsrep(i, xx - 1, x) if (f[i][yy]) return false;
    srep(j, yy + 1, y) if (f[xx][j]) return false;

    return true;
  }

  // 下右
  if (z == 2) {
    int xx = x + diff1;
    int yy = y + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][1] || line[xx][yy][0]) return false;

    // 間に邪魔な頂点がないか
    dsrep(j, yy - 1, y) if (f[xx][j]) return false;
    dsrep(i, xx - 1, x) if (f[i][yy]) return false;

    return true;
  }

  // 右上
  if (z == 3) {
    int xx = x - diff2;
    int yy = y + diff1;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][2] || line[xx][yy][1]) return false;

    // 間に邪魔な頂点がないか
    srep(i, xx + 1, x) if (f[i][yy]) return false;
    dsrep(j, yy - 1, y) if (f[xx][j]) return false;

    return true;
  }

  // 左上・左下
  if (z == 4) {
    int xx = x - diff1 + diff2;
    int yy = y - diff1 - diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][7] || line[xx][yy][6]) return false;

    // 間に邪魔な頂点がないか
    srep(i, 1, diff2) if (f[xx - i][yy + i]) return false;
    srep(i, 1, diff1) if (f[xx + i][yy + i]) return false;

    return true;
  }

  // 左下・右下
  if (z == 5) {
    int xx = x + diff1 + diff2;
    int yy = y - diff1 + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][4] || line[xx][yy][7]) return false;

    // 間に邪魔な頂点がないか
    srep(i, 1, diff2) if (f[xx - i][yy - i]) return false;
    srep(i, 1, diff1) if (f[xx - i][yy + i]) return false;

    return true;
  }

  // 右下・右上
  if (z == 6) {
    int xx = x + diff1 - diff2;
    int yy = y + diff1 + diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][5] || line[xx][yy][4]) return false;

    // 間に邪魔な頂点がないか
    srep(i, 1, diff2) if (f[xx + i][yy - i]) return false;
    srep(i, 1, diff1) if (f[xx - i][yy - i]) return false;

    return true;
  }

  // 右上・左上
  if (z == 7) {
    int xx = x - diff1 - diff2;
    int yy = y + diff1 - diff2;

    // その点がグリッド内か
    if (IsNGXY(xx, yy)) return false;

    // そこに点が存在しているか
    if (!f[xx][yy]) return false;

    // 辺が既に存在していないか
    if (line[xx][yy][6] || line[xx][yy][5]) return false;

    // 間に邪魔な頂点がないか
    srep(i, 1, diff2) if (f[xx + i][yy + i]) return false;
    srep(i, 1, diff1) if (f[xx + i][yy - i]) return false;

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
  int x = randxor() % N;
  int y = randxor() % N;

  if (f[x][y]) return;

  methodCount[1][1]++;

  int u = -1;
  if (cntH[y] != 0) {
    u = FindNeighborPoint(x, y, 0);
    if (u == -2) return;
  }
  int l = -1;
  if (cntW[x] != 0) {
    l = FindNeighborPoint(x, y, 1);
    if (l == -2) return;
  }
  int d = -1;
  if (cntH[y] != 0) {
    d = FindNeighborPoint(x, y, 2);
    if (d == -2) return;
  }
  int r = -1;
  if (cntW[x] != 0) {
    r = FindNeighborPoint(x, y, 3);
    if (r == -2) return;
  }
  int ul = FindNeighborPoint(x, y, 4);
  if (ul == -2) return;
  int ld = FindNeighborPoint(x, y, 5);
  if (ld == -2) return;
  int dr = FindNeighborPoint(x, y, 6);
  if (dr == -2) return;
  int ru = FindNeighborPoint(x, y, 7);
  if (ru == -2) return;

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

  if (RectDir == -1) return;

  double diffScore = 1000000.0 * N * N / M * W[x][y] / S;

  double prob = exp(diffScore / temperature);
  if (prob > rand01()) {
    methodCount[1][0]++;

    maxScore += diffScore;

    ans[ansSize][0][0] = x;
    ans[ansSize][0][1] = y;
    ans[ansSize][1][0] = x1;
    ans[ansSize][1][1] = y1;
    ans[ansSize][2][0] = xx;
    ans[ansSize][2][1] = yy;
    ans[ansSize][3][0] = x3;
    ans[ansSize][3][1] = y3;
    ansDelete[ansSize] = 0;
    ansSize++;
    f[x][y] = 1;
    cntW[x]++;
    cntH[y]++;
    use[x][y]++;
    use[x1][y1]++;
    use[xx][yy]++;
    use[x3][y3]++;

    if (RectDir == 0) {
      line[x][y][0] = true;
      line[x][y][1] = true;
      line[x1][y1][1] = true;
      line[x1][y1][2] = true;
      line[xx][yy][2] = true;
      line[xx][yy][3] = true;
      line[x3][y3][3] = true;
      line[x3][y3][0] = true;
    }
    else if (RectDir == 1) {
      line[x][y][1] = true;
      line[x][y][2] = true;
      line[x1][y1][2] = true;
      line[x1][y1][3] = true;
      line[xx][yy][3] = true;
      line[xx][yy][0] = true;
      line[x3][y3][0] = true;
      line[x3][y3][1] = true;
    }
    else if (RectDir == 2) {
      line[x][y][2] = true;
      line[x][y][3] = true;
      line[x1][y1][3] = true;
      line[x1][y1][0] = true;
      line[xx][yy][0] = true;
      line[xx][yy][1] = true;
      line[x3][y3][1] = true;
      line[x3][y3][2] = true;
    }
    else if (RectDir == 3) {
      line[x][y][3] = true;
      line[x][y][0] = true;
      line[x1][y1][0] = true;
      line[x1][y1][1] = true;
      line[xx][yy][1] = true;
      line[xx][yy][2] = true;
      line[x3][y3][2] = true;
      line[x3][y3][3] = true;
    }
    else if (RectDir == 4) {
      line[x][y][4] = true;
      line[x][y][5] = true;
      line[x1][y1][5] = true;
      line[x1][y1][6] = true;
      line[xx][yy][6] = true;
      line[xx][yy][7] = true;
      line[x3][y3][7] = true;
      line[x3][y3][4] = true;
    }
    else if (RectDir == 5) {
      line[x][y][5] = true;
      line[x][y][6] = true;
      line[x1][y1][6] = true;
      line[x1][y1][7] = true;
      line[xx][yy][7] = true;
      line[xx][yy][4] = true;
      line[x3][y3][4] = true;
      line[x3][y3][5] = true;
    }
    else if (RectDir == 6) {
      line[x][y][6] = true;
      line[x][y][7] = true;
      line[x1][y1][7] = true;
      line[x1][y1][4] = true;
      line[xx][yy][4] = true;
      line[xx][yy][5] = true;
      line[x3][y3][5] = true;
      line[x3][y3][6] = true;
    }
    else if (RectDir == 7) {
      line[x][y][7] = true;
      line[x][y][4] = true;
      line[x1][y1][4] = true;
      line[x1][y1][5] = true;
      line[xx][yy][5] = true;
      line[xx][yy][6] = true;
      line[x3][y3][6] = true;
      line[x3][y3][7] = true;
    }

    if (maxScore > real_maxScore) {
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
  if (ansSize == 0) return;
  int ite = randxor() % ansSize;
  if (ansDelete[ite]) return;
  if (use[ans[ite][0][0]][ans[ite][0][1]] > 1) return;

  methodCount[2][1]++;

  int x[4], y[4];
  rep(i, 4)
  {
    x[i] = ans[ite][i][0];
    y[i] = ans[ite][i][1];
  }

  double diffScore = -1000000.0 * N * N / M * W[x[0]][y[0]] / S;

  double prob = exp(diffScore / temperature);
  if (prob > rand01()) {
    methodCount[2][0]++;

    maxScore += diffScore;

    f[x[0]][y[0]] = 0;
    cntW[x[0]]--;
    cntH[y[0]]--;
    rep(i, 4) use[x[i]][y[i]]--;

    int RectDir = GetDir(x[0], y[0], x[1], y[1]);
    if (RectDir == 0) {
      line[x[0]][y[0]][0] = false;
      line[x[0]][y[0]][1] = false;
      line[x[1]][y[1]][1] = false;
      line[x[1]][y[1]][2] = false;
      line[x[2]][y[2]][2] = false;
      line[x[2]][y[2]][3] = false;
      line[x[3]][y[3]][3] = false;
      line[x[3]][y[3]][0] = false;
    }
    else if (RectDir == 1) {
      line[x[0]][y[0]][1] = false;
      line[x[0]][y[0]][2] = false;
      line[x[1]][y[1]][2] = false;
      line[x[1]][y[1]][3] = false;
      line[x[2]][y[2]][3] = false;
      line[x[2]][y[2]][0] = false;
      line[x[3]][y[3]][0] = false;
      line[x[3]][y[3]][1] = false;
    }
    else if (RectDir == 2) {
      line[x[0]][y[0]][2] = false;
      line[x[0]][y[0]][3] = false;
      line[x[1]][y[1]][3] = false;
      line[x[1]][y[1]][0] = false;
      line[x[2]][y[2]][0] = false;
      line[x[2]][y[2]][1] = false;
      line[x[3]][y[3]][1] = false;
      line[x[3]][y[3]][2] = false;
    }
    else if (RectDir == 3) {
      line[x[0]][y[0]][3] = false;
      line[x[0]][y[0]][0] = false;
      line[x[1]][y[1]][0] = false;
      line[x[1]][y[1]][1] = false;
      line[x[2]][y[2]][1] = false;
      line[x[2]][y[2]][2] = false;
      line[x[3]][y[3]][2] = false;
      line[x[3]][y[3]][3] = false;
    }
    else if (RectDir == 4) {
      line[x[0]][y[0]][4] = false;
      line[x[0]][y[0]][5] = false;
      line[x[1]][y[1]][5] = false;
      line[x[1]][y[1]][6] = false;
      line[x[2]][y[2]][6] = false;
      line[x[2]][y[2]][7] = false;
      line[x[3]][y[3]][7] = false;
      line[x[3]][y[3]][4] = false;
    }
    else if (RectDir == 5) {
      line[x[0]][y[0]][5] = false;
      line[x[0]][y[0]][6] = false;
      line[x[1]][y[1]][6] = false;
      line[x[1]][y[1]][7] = false;
      line[x[2]][y[2]][7] = false;
      line[x[2]][y[2]][4] = false;
      line[x[3]][y[3]][4] = false;
      line[x[3]][y[3]][5] = false;
    }
    else if (RectDir == 6) {
      line[x[0]][y[0]][6] = false;
      line[x[0]][y[0]][7] = false;
      line[x[1]][y[1]][7] = false;
      line[x[1]][y[1]][4] = false;
      line[x[2]][y[2]][4] = false;
      line[x[2]][y[2]][5] = false;
      line[x[3]][y[3]][5] = false;
      line[x[3]][y[3]][6] = false;
    }
    else if (RectDir == 7) {
      line[x[0]][y[0]][7] = false;
      line[x[0]][y[0]][4] = false;
      line[x[1]][y[1]][4] = false;
      line[x[1]][y[1]][5] = false;
      line[x[2]][y[2]][5] = false;
      line[x[2]][y[2]][6] = false;
      line[x[3]][y[3]][6] = false;
      line[x[3]][y[3]][7] = false;
    }

    ansDelete[ite] = 1;
    ansDeleteCount++;
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
  maxScore = CalcScore();
  CopyToReal();
  CopyToSeed();

  // シード作り
  int seedCount = 20;  // 0にするとシード作成を行わない
  rep(tei, seedCount)
  {
    startTime = clock();

    // 初期状態に戻す
    Init();
    maxScore = CalcScore();

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
      if (nowProgress > 1.0) break;

      // 現在のスコアが悪いときは元に戻す
      if (maxScore * 1.2 < real_maxScore) {
        RollBackFromReal();
        rollbackCount++;
      }

      if (ansDeleteCount >= 10000) {
        RefleshAns();
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      // メソッド選択
      int me = 1;
      if (randxor() % 2 == 0) {
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

    // スコアが良ければシードを更新
    RollBackFromReal();
    if (maxScore > seed_maxScore) {
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
    if (nowProgress > 1.0) break;

    // 現在のスコアが悪いときは元に戻す
    if (maxScore * 1.2 < real_maxScore) {
      RollBackFromReal();
      rollbackCount++;
    }

    if (ansDeleteCount >= 10000) {
      RefleshAns();
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    // メソッド選択
    int me = 1;
    if (randxor() % 2 == 0) {
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

  // 一番スコアの良い解
  RollBackFromReal();

  RefleshAns();

  CalcScore();

  // デバッグログ
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", N = " << N << endl;
    cout << "ansSize = " << ansSize << ", ansDeleteCount = " << ansDeleteCount
      << endl;
    cout << "maxScore = " << maxScore << endl;
    cout << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    srep(i, 1, 5)
    {
      cout << "Method" << i << " = " << methodCount[i][0] << " / "
        << methodCount[i][1] << endl;
    }
    cout << endl;
  }

  cerr << loop << endl;
  return maxScore;
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

  // 乱数調整
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  int mode = 0;

  // 提出用
  if (mode == 0) {
    rep(i, 1)
    {
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
    rep(i, 100)
    {
      scoreSum += SolveOuter(mode, i);
      AllClear_MultiCase();
    }
    cout << "scoreSum = " << scoreSum << endl;
  }

  mainEnd = clock();
  cerr << (double)(mainEnd - mainStart) / CLOCKS_PER_SEC;
  return 0;
}
