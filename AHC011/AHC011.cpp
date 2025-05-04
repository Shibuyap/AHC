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
#define MAX_N 100

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

///////////////////////////////////////////
// 上：2
// 下：8
// 左：1
// 右：4
//
///////////////////////////////////////////

const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char cc[4] = { 'U', 'L', 'D', 'R' };

namespace /* 乱数ライブラリ */
{
  static uint32_t Rand() {
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


  static double Rand01() {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* 変数 */
{
  // 入力用変数
  int n, t;
  int board[10][10];
  int startX, startY;
  int cnt[16];

  // 解答用変数
  ll maxScore;
  vector<int> ans;
  int keepAns[2100];

  // 焼きなまし用変数
  double TL = 2.9;
  double start_temp = 2048;
  double end_temp = 0.0001;
  ll best_maxScore;
  vector<int> best_ans;

}  // namespace

namespace /* Union Find*/
{
  int parUF[MAX_N];   // 親
  int rankUF[MAX_N];  // 木の深さ
  int cntUF[MAX_N];   // 属する頂点の個数(親のみ正しい)

  // n要素で初期化
  void UFinit(const int nn) {
    for (int i = 0; i < nn; i++) {
      parUF[i] = i;
      rankUF[i] = 0;
      cntUF[i] = 1;
      int val = board[i / n][i % n];
      if (val == 1 || val == 2 || val == 4 || val == 8) {
        cntUF[i] = 1;
      }
    }
  }

  // 木の根を求める
  int findUF(int x) {
    if (parUF[x] == x) {
      return x;
    }
    else {
      return parUF[x] = findUF(parUF[x]);
    }
  }

  // xとyの属する集合を併合
  void uniteUF(int x, int y) {
    x = findUF(x);
    y = findUF(y);
    if (x == y) return;

    if (rankUF[x] < rankUF[y]) {
      parUF[x] = y;
      cntUF[y] += cntUF[x];
    }
    else {
      parUF[y] = x;
      cntUF[x] += cntUF[y];
      if (rankUF[x] == rankUF[y]) rankUF[x]++;
    }
  } /* Union Find*/

  // xとyが同じ集合に属するか否か
  bool sameUF(int x, int y) {
    return findUF(x) == findUF(y);
  }
}  // namespace

// スコア計算
int boardForCalc[10][10];
int CalcScore(const vector<int>& ope) {
  UFinit(n * n);

  int x = startX, y = startY;
  rep(i, n) {
    rep(j, n) {
      boardForCalc[i][j] = board[i][j];
    }
  }

  rep(i, ope.size()) {
    swap(boardForCalc[x][y], boardForCalc[x + dx[ope[i]]][y + dy[ope[i]]]);
    x += dx[ope[i]];
    y += dy[ope[i]];
  }

  // 横の繋がり
  rep(i, n) {
    rep(j, n - 1) {
      if ((boardForCalc[i][j] & (1 << 2)) && (boardForCalc[i][j + 1] & (1 << 0))) {
        uniteUF(i * n + j, i * n + j + 1);
      }
    }
  }

  // 縦の繋がり
  rep(i, n - 1) {
    rep(j, n) {
      if ((boardForCalc[i][j] & (1 << 3)) && (boardForCalc[i + 1][j] & (1 << 1))) {
        uniteUF(i * n + j, (i + 1) * n + j);
      }
    }
  }

  int res = 0;

  // スコア計算
  rep(i, n) {
    rep(j, n) {
      int ij = i * n + j;
      if (findUF(ij) == ij) {
        if (cntUF[ij] == n * n - 1) {
          res = 500000.0 * (2.0 - (double)ope.size() / t);
          return res;
        }
        else {
          res = max(res, (int)round(500000.0 * cntUF[ij] / (n * n - 1)));
        }
      }
    }
  }

  return res;
}

// 入力受け取り
void Input(int);
void ResetAll();

bool IsOK(int x, int y) {
  if (x < 0 || n <= x || y < 0 || n <= y) {
    return false;
  }
  return true;
}

bool IsOKRoute(const vector<int>& ope) {
  int x = startX;
  int y = startY;
  rep(i, ope.size()) {
    x += dx[ope[i]];
    y += dy[ope[i]];
    if (!IsOK(x, y)) {
      return false;
    }
  }
  return true;
}

// kから後ろを全リセット
void op_shuffle_suffix(double temp) {
  int ite = Rand() % t;

  int x = startX;
  int y = startY;
  rep(i, t) {
    keepAns[i] = ans[i];
    if (i < ite) {
      x += dx[ans[i]];
      y += dy[ans[i]];
    }
    else {
      int val = Rand() % 4;
      while (!IsOK(x + dx[val], y + dy[val])) {
        val = Rand() % 4;
      }

      ans[i] = val;
      x += dx[val];
      y += dy[val];
    }
  }

  int tmpScore = CalcScore(ans);

  int diffScore = tmpScore - maxScore;

  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    maxScore += diffScore;
    if (maxScore > best_maxScore) {
      best_maxScore = maxScore;
      best_ans = ans;
    }
  }
  else {
    // 元に戻す
    rep(i, t) {
      ans[i] = keepAns[i];
    }
  }
}

// kとkの直後をスワップ
void op_swap_adjacent(double temp) {
  int ite = Rand() % (t - 1);

  swap(ans[ite], ans[ite + 1]);

  if (!IsOKRoute(ans)) {
    swap(ans[ite], ans[ite + 1]);
    return;
  }

  int tmpScore = CalcScore(ans);

  int diffScore = tmpScore - maxScore;

  double prob = exp((double)diffScore / temp);
  if (prob > Rand01()) {
    maxScore += diffScore;
    if (maxScore > best_maxScore) {
      best_maxScore = maxScore;
      best_ans = ans;
    }
  }
  else {
    // 元に戻す
    swap(ans[ite], ans[ite + 1]);
  }
}

// 木を全探索
namespace
{
  int dfsBoard[10][10];
  int dfsCnt[16];
  void DfsInit() {
    rep(i, n) {
      rep(j, n) {
        dfsBoard[i][j] = -1;
      }
    }
    rep(i, 16) {
      dfsCnt[i] = cnt[i];
    }
  }

  bool CheckDfs(int x, int y) {
    if (x == 0) {
      if (dfsBoard[x][y] & (1 << 1)) {
        return false;
      }
    }
    else {
      if (dfsBoard[x][y] & (1 << 1)) {
        if (!(dfsBoard[x - 1][y] & (1 << 3))) {
          return false;
        }
      }
      else {
        if ((dfsBoard[x - 1][y] & (1 << 3))) {
          return false;
        }
      }
    }

    if (y == 0) {
      if (dfsBoard[x][y] & (1 << 0)) {
        return false;
      }
    }
    else {
      if (dfsBoard[x][y] & (1 << 0)) {
        if (!(dfsBoard[x][y - 1] & (1 << 2))) {
          return false;
        }
      }
      else {
        if ((dfsBoard[x][y - 1] & (1 << 2))) {
          return false;
        }
      }
    }

    int val = dfsBoard[x][y];
    if (val == 2 && x > 0 && dfsBoard[x - 1][y] == 8) {
      return false;
    }
    if (val == 1 && y > 0 && dfsBoard[x][y - 1] == 4) {
      return false;
    }

    return true;
  }

  bool CheckAllDfs() {
    UFinit(n * n);
    // 横の繋がり
    rep(i, n) {
      rep(j, n - 1) {
        if ((dfsBoard[i][j] & (1 << 2)) && (dfsBoard[i][j + 1] & (1 << 0))) {
          uniteUF(i * n + j, i * n + j + 1);
        }
      }
    }

    // 縦の繋がり
    rep(i, n - 1) {
      rep(j, n) {
        if ((dfsBoard[i][j] & (1 << 3)) && (dfsBoard[i + 1][j] & (1 << 1))) {
          uniteUF(i * n + j, (i + 1) * n + j);
        }
      }
    }

    if (cntUF[findUF(0)] == n * n - 1) {
      return true;
    }
    return false;
  }

  bool CheckAllDfs(const vector<vector<int>>& vec) {
    UFinit(n * n);
    // 横の繋がり
    rep(i, n) {
      rep(j, n - 1) {
        if ((vec[i][j] & (1 << 2)) && (vec[i][j + 1] & (1 << 0))) {
          uniteUF(i * n + j, i * n + j + 1);
        }
      }
    }

    // 縦の繋がり
    rep(i, n - 1) {
      rep(j, n) {
        if ((vec[i][j] & (1 << 3)) && (vec[i + 1][j] & (1 << 1))) {
          uniteUF(i * n + j, (i + 1) * n + j);
        }
      }
    }

    if (cntUF[findUF(0)] == n * n - 1) {
      return true;
    }
    return false;
  }

  int dfs_search(int ite) {
    int x = ite / n;
    int y = ite % n;
    if (ite == n * n - 1) {
      dfsBoard[x][y] = 0;
      if (!CheckDfs(x, y)) {
        return 0;
      }
      if (!CheckAllDfs()) {
        return 0;
      }
      return 1;
    }

    srep(i, 1, 16) {
      if (dfsCnt[i] == 0) {
        continue;
      }
      dfsBoard[x][y] = i;
      dfsCnt[i]--;
      if (!CheckDfs(x, y)) {
        dfsCnt[i]++;
        continue;
      }
      if (dfs_search(ite + 1)) {
        return 1;
      }
      dfsCnt[i]++;
    }
    return 0;
  }

  void dfs_print_board() {
    rep(i, n) {
      rep(j, n) {
        cout << dfsBoard[i][j] << ' ';
      }
      cout << endl;
    }
  }
}  // namespace

// 木を焼きなましで見つける
int aniBoard[10][10];
int bestMaxAniBoard[10][10];
int CalcAniScore() {
  UFinit(n * n);
  // 横の繋がり
  rep(i, n) {
    rep(j, n - 1) {
      if ((aniBoard[i][j] & (1 << 2)) && (aniBoard[i][j + 1] & (1 << 0))) {
        uniteUF(i * n + j, i * n + j + 1);
      }
    }
  }

  // 縦の繋がり
  rep(i, n - 1) {
    rep(j, n) {
      if ((aniBoard[i][j] & (1 << 3)) && (aniBoard[i + 1][j] & (1 << 1))) {
        uniteUF(i * n + j, (i + 1) * n + j);
      }
    }
  }

  int MaxSize = 0;
  rep(i, n) {
    rep(j, n) {
      MaxSize = max(MaxSize, cntUF[findUF(i * n + j)]);
    }
    if (MaxSize >= n * n / 2) {
      break;
    }
  }

  return MaxSize;
}

bool FindTreeAni(bool isReset = false) {
  clock_t startAniTime, endAniTime;
  const double AniTL = 0.1;
  startAniTime = clock();
  endAniTime = clock();
  double nowAniTime = (double)(endAniTime - startAniTime) / CLOCKS_PER_SEC;
  int loopAni = 0;
  double startAniTemp = 0.1;
  double endAniTemp = 0.0;

  if (isReset) {
    rep(i, n) {
      rep(j, n) {
        aniBoard[i][j] = board[i][j];
      }
    }
  }

  // カーソルは右下固定
  swap(aniBoard[startX][startY], aniBoard[n - 1][n - 1]);

  int maxAniScore = CalcAniScore();
  int bestMaxAniScore = maxAniScore;
  rep(i, n) {
    rep(j, n) {
      bestMaxAniBoard[i][j] = aniBoard[i][j];
    }
  }

  while (true) {
    if (loopAni % 100 == 1) {
      endAniTime = clock();
      nowAniTime = (double)(endAniTime - startAniTime) / CLOCKS_PER_SEC;
      if (nowAniTime > AniTL) break;
    }

    int x1 = Rand() % n;
    int y1 = Rand() % n;
    int x2 = Rand() % n;
    while (x1 == x2) {
      x2 = Rand() % n;
    }
    int y2 = Rand() % n;
    while (y1 == y2) {
      y2 = Rand() % n;
    }

    if (x1 == n - 1 && y1 == n - 1) {
      continue;
    }
    if (x2 == n - 1 && y2 == n - 1) {
      continue;
    }

    swap(aniBoard[x1][y1], aniBoard[x2][y2]);
    int newPoint = CalcAniScore();

    int diffScore = newPoint - maxAniScore;

    double temp = startAniTemp + (endAniTemp - startAniTemp) * nowAniTime / AniTL;
    double prob = exp((double)diffScore / temp);
    if (prob > Rand01()) {
      maxAniScore += diffScore;
      if (maxAniScore > bestMaxAniScore) {
        bestMaxAniScore = maxAniScore;
        rep(i, n) {
          rep(j, n) {
            bestMaxAniBoard[i][j] = aniBoard[i][j];
          }
        }
      }
    }
    else {
      // 元に戻す
      swap(aniBoard[x1][y1], aniBoard[x2][y2]);
    }

    loopAni++;
    if (bestMaxAniScore == n * n - 1) {
      break;
    }
  }

  maxAniScore = bestMaxAniScore;
  rep(i, n) {
    rep(j, n) {
      aniBoard[i][j] = bestMaxAniBoard[i][j];
    }
  }

  if (maxAniScore == n * n - 1) {
    return true;
  }
  return false;
}

// 作成した木からピースの種類ごとの番号を決定する
vector<int> kindNumbers[16];
vector<P> originNum[16];
void InitKindNumbers() {
  rep(i, n) {
    rep(j, n) {
      kindNumbers[dfsBoard[i][j]].push_back(i * n + j);
      originNum[board[i][j]].push_back(P(i, j));
    }
  }
}

int peaceNum[10][10];
// 作成した盤面の転倒数をチェックする
bool CheckInversion() {
  int tmpBoard[10][10];
  rep(i, n) {
    rep(j, n) {
      tmpBoard[i][j] = peaceNum[i][j];
    }
  }
  int cnt = 0;

  rep(i, n) {
    rep(j, n) {
      int num = i * n + j;
      int x = -1;
      int y = -1;
      rep(k, n) {
        rep(l, n) {
          if (tmpBoard[k][l] == num) {
            x = k;
            y = j;
          }
        }
      }
      while (y < i) {
        cnt++;
        swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
        y++;
      }
      while (y > i) {
        cnt++;
        swap(tmpBoard[x][y], tmpBoard[x][y - 1]);
        y--;
      }
      while (x < i) {
        cnt++;
        swap(tmpBoard[x][y], tmpBoard[x + 1][y]);
        x++;
      }
      while (x > i) {
        cnt++;
        swap(tmpBoard[x][y], tmpBoard[x - 1][y]);
        x--;
      }
    }
  }

  if (cnt % 2 == 0) {
    return true;
  }
  return false;
}

// ピースに番号を振る
void init_piece_numbers() {
  int ite[16] = {};
  rep(i, n) {
    rep(j, n) {
      peaceNum[i][j] = kindNumbers[board[i][j]][ite[board[i][j]]];
      ite[board[i][j]]++;
    }
  }

  // もし転倒数が奇数なら1箇所スワップする
  if (!CheckInversion()) {
    rep(i, 16) {
      if (kindNumbers[i].size() >= 2) {
        int num1 = kindNumbers[i][0];
        int num2 = kindNumbers[i][1];
        int x1, y1, x2, y2;
        rep(j, n) {
          rep(k, n) {
            if (peaceNum[j][k] == num1) {
              x1 = j;
              y1 = k;
            }
            if (peaceNum[j][k] == num2) {
              x2 = j;
              y2 = k;
            }
          }
        }
        swap(peaceNum[x1][y1], peaceNum[x2][y2]);
        break;
      }
    }
  }
}

pair<P, P> shuffle_same_kind_piece() {
  int ite = Rand() % 16;
  while (originNum[ite].size() <= 1) {
    ite = Rand() % 16;
  }

  int a = Rand() % originNum[ite].size();
  int b = Rand() % originNum[ite].size();
  while (a == b) {
    b = Rand() % originNum[ite].size();
  }

  int x1 = originNum[ite][a].first;
  int y1 = originNum[ite][a].second;
  int x2 = originNum[ite][b].first;
  int y2 = originNum[ite][b].second;
  swap(peaceNum[x1][y1], peaceNum[x2][y2]);

  return pair<P, P>({ {x1, y1}, {x2, y2} });
}

// 木を作成する手順を1つ作成する
void apply_move(int& x, int& y, int nd, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  ansDfs.push_back(nd);
  swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
  x += dx[nd];
  y += dy[nd];
}

void route_move_cursor(int& x, int& y, int& xx, int& yy, int tx, int ty, int ii, int jj, vector<vector<int>>& tmpBoard, vector<int>& ansDfs, int mode = 0) {
  int dp[10][10];
  int dir[10][10];
  if (mode == 0) {
    rep(i, n) {
      rep(j, n) {
        if (i < ii) {
          dp[i][j] = -1;
        }
        else if (i == ii && j < jj) {
          dp[i][j] = -1;
        }
        else {
          dp[i][j] = INF;
        }
      }
    }
  }
  else if (mode == 2) {
    rep(i, n) {
      rep(j, n) {
        if (i < ii) {
          dp[i][j] = -1;
        }
        else if (j < jj) {
          dp[i][j] = -1;
        }
        else {
          dp[i][j] = INF;
        }
      }
    }
  }

  queue<P> que;
  que.push(P(x, y));
  dp[x][y] = 0;
  dp[xx][yy] = -1;
  dir[x][y] = -1;
  while (que.size()) {
    int a = que.front().first;
    int b = que.front().second;
    que.pop();
    bool isFinish = false;
    rep(i, 4) {
      int na = a + dx[i];
      int nb = b + dy[i];
      if (IsOK(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == tx && nb == ty) {
          isFinish = true;
          break;
        }
      }
    }
    if (isFinish) {
      break;
    }
  }

  vector<int> rev;
  int revX = tx, revY = ty;
  while (revX != x || revY != y) {
    int revD = dir[revX][revY];
    rev.push_back(revD);
    revX -= dx[revD];
    revY -= dy[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }
  rep(i, 4) {
    if (x + dx[i] == xx && y + dy[i] == yy) {
      swap(tmpBoard[x][y], tmpBoard[xx][yy]);
      swap(x, xx);
      swap(y, yy);
      ansDfs.push_back(i);
      break;
    }
  }
}

// 3×2マスを使って上にマスを入れ替える
// (xx,yy) = (左上,右上)
void swap_vertical_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  while (y != yy + 1) {
    ansDfs.push_back(3);
    swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
    y++;
  }
  // 上左下右
  vector<int> order = { 0, 1, 2, 3, 2, 1, 0, 0, 3, 2, 1, 2, 3, 0, 0, 1, 2 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }
}

// 2×3マスを使って上にマスを入れ替える
// (xx,yy) = (左上,右上)
void swap_horizontal_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  while (x != n - 1) {
    apply_move(x, y, 2, tmpBoard, ansDfs);
  }
  // 上左下右
  vector<int> order = { 1, 0, 3, 2, 3, 0, 1, 1, 2, 3, 0, 3, 2, 1, 1, 0, 3 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }
}

void fix_last_two_in_row(int& x, int& y, int ii, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  // カーソル移動
  int dp[10][10];
  int dir[10][10];
  rep(i, n) {
    rep(j, n) {
      if (i <= ii) {
        dp[i][j] = -1;
      }
      else {
        dp[i][j] = INF;
      }
    }
  }

  queue<P> que;
  que.push(P(x, y));
  dp[x][y] = 0;
  dp[ii + 1][n - 2] = -1;
  dp[ii][n - 1] = INF;
  dir[x][y] = -1;
  while (que.size()) {
    int a = que.front().first;
    int b = que.front().second;
    que.pop();
    bool isFinish = false;
    rep(i, 4) {
      int na = a + dx[i];
      int nb = b + dy[i];
      if (IsOK(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == ii && nb == n - 1) {
          isFinish = true;
          break;
        }
      }
    }
    if (isFinish) {
      break;
    }
  }

  vector<int> rev;
  int revX = ii, revY = n - 1;
  while (revX != x || revY != y) {
    int revD = dir[revX][revY];
    rev.push_back(revD);
    revX -= dx[revD];
    revY -= dy[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }

  // 上左下右
  vector<int> order = { 1, 2 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }
}

void fix_last_two_in_col(int& x, int& y, int jj, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  // カーソル移動
  int dp[10][10];
  int dir[10][10];
  rep(i, n) {
    rep(j, n) {
      if (i < n - 2) {
        dp[i][j] = -1;
      }
      else if (j <= jj) {
        dp[i][j] = -1;
      }
      else {
        dp[i][j] = INF;
      }
    }
  }

  queue<P> que;
  que.push(P(x, y));
  dp[x][y] = 0;
  dp[n - 2][jj + 1] = -1;
  dp[n - 1][jj] = INF;
  dir[x][y] = -1;
  while (que.size()) {
    int a = que.front().first;
    int b = que.front().second;
    que.pop();
    bool isFinish = false;
    rep(i, 4) {
      int na = a + dx[i];
      int nb = b + dy[i];
      if (IsOK(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == n - 1 && nb == jj) {
          isFinish = true;
          break;
        }
      }
    }
    if (isFinish) {
      break;
    }
  }

  vector<int> rev;
  int revX = n - 1, revY = jj;
  while (revX != x || revY != y) {
    int revD = dir[revX][revY];
    rev.push_back(revD);
    revX -= dx[revD];
    revY -= dy[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }

  // 上左下右
  vector<int> order = { 0, 3 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + dx[nd]][y + dy[nd]]);
    x += dx[nd];
    y += dy[nd];
  }
}

vector<int> FindAnsDfs() {
  vector<vector<int>> tmpBoard(10, vector<int>(10));
  rep(i, n) {
    rep(j, n) {
      tmpBoard[i][j] = peaceNum[i][j];
    }
  }

  vector<int> ansDfs;

  int x = startX;
  int y = startY;

  // 1列目から下から3列目まで完成させる
  rep(i, n - 2) {
    rep(j, n - 1) {
      int num = i * n + j;
      if (j == n - 2) {
        num = i * n + j + 1;
      }
      int xx = -1, yy = -1;
      rep(k, n) {
        rep(l, n) {
          if (tmpBoard[k][l] == num) {
            xx = k;
            yy = l;
          }
        }
        if (xx != -1) {
          break;
        }
      }

      // 適切な位置にピースを移動させる
      while (xx != i || yy != j) {
        if (yy != j) {
          if (yy < j) {
            route_move_cursor(x, y, xx, yy, xx, yy + 1, i, j, tmpBoard, ansDfs);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx, yy - 1, i, j, tmpBoard, ansDfs);
          }
        }
        else if (xx != i) {
          route_move_cursor(x, y, xx, yy, xx - 1, yy, i, j, tmpBoard, ansDfs);
        }
      }
    }
    if (x == i && y == n - 1) {
      ansDfs.push_back(2);
      swap(tmpBoard[x][y], tmpBoard[x + 1][y]);
      x++;
    }

    // 最後の2個を揃える
    if (tmpBoard[i][n - 1] == i * n + n - 2) {
      swap_vertical_pair(x, y, i, n - 2, tmpBoard, ansDfs);
    }
    else {
      int num = i * n + n - 2;
      int xx = -1, yy = -1;
      rep(k, n) {
        rep(l, n) {
          if (tmpBoard[k][l] == num) {
            xx = k;
            yy = l;
          }
        }
        if (xx != -1) {
          break;
        }
      }
      // 適切な位置にピースを移動させる
      while (xx != i + 1 || yy != n - 2) {
        if (yy != n - 2) {
          if (yy < n - 2) {
            route_move_cursor(x, y, xx, yy, xx, yy + 1, i, n - 1, tmpBoard, ansDfs);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx, yy - 1, i, n - 1, tmpBoard, ansDfs);
          }
        }
        else if (xx != i + 1) {
          route_move_cursor(x, y, xx, yy, xx - 1, yy, i, n - 1, tmpBoard, ansDfs);
        }
      }
      fix_last_two_in_row(x, y, i, tmpBoard, ansDfs);
    }
  }

  // 下2列をそろえる
  rep(j, n - 2) {
    // 下のマスを上のマスの位置に持ってくる
    int num = (n - 1) * n + j;
    int xx = -1, yy = -1;
    rep(k, n) {
      rep(l, n) {
        if (tmpBoard[k][l] == num) {
          xx = k;
          yy = l;
        }
      }
      if (xx != -1) {
        break;
      }
    }

    // 適切な位置にピースを移動させる
    while (xx != n - 2 || yy != j) {
      if (yy != j) {
        if (yy < j) {
          route_move_cursor(x, y, xx, yy, xx, yy + 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
        else {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      else if (xx != n - 2) {
        route_move_cursor(x, y, xx, yy, xx - 1, yy, n - 2, j, tmpBoard, ansDfs);
      }
    }

    if (x == n - 1 && y == j) {
      ansDfs.push_back(3);
      swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
      y++;
    }

    // 最後の2個を揃える
    if (tmpBoard[n - 1][j] == (n - 2) * n + j) {
      swap_horizontal_pair(x, y, n - 2, j, tmpBoard, ansDfs);
    }
    else {
      int num = (n - 2) * n + j;
      int xx = -1, yy = -1;
      rep(k, n) {
        rep(l, n) {
          if (tmpBoard[k][l] == num) {
            xx = k;
            yy = l;
          }
        }
        if (xx != -1) {
          break;
        }
      }
      // 適切な位置にピースを移動させる
      while (xx != n - 2 || yy != j + 1) {
        if (xx != n - 2) {
          if (xx < n - 2) {
            route_move_cursor(x, y, xx, yy, xx + 1, yy, n - 2, j, tmpBoard, ansDfs, 2);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, n - 2, j, tmpBoard, ansDfs, 2);
          }
        }
        else if (yy != j + 1) {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      fix_last_two_in_col(x, y, j, tmpBoard, ansDfs);
    }
  }

  // カーソルを右下に持ってくる
  while (x != n - 1 || y != n - 1) {
    if (y != n - 1) {
      apply_move(x, y, 3, tmpBoard, ansDfs);
    }
    else {
      apply_move(x, y, 2, tmpBoard, ansDfs);
    }
  }

  return ansDfs;
}

vector<int> FindAnsDfs2() {
  vector<vector<int>> tmpBoard(10, vector<int>(10));
  rep(i, n) {
    rep(j, n) {
      tmpBoard[i][j] = board[i][j];
    }
  }

  vector<int> ansDfs;

  int x = startX;
  int y = startY;

  // 1列目から下から3列目まで完成させる
  rep(i, n - 2) {
    rep(j, n - 1) {
      int num = aniBoard[i][j];
      if (j == n - 2) {
        num = aniBoard[i][n - 1];
      }

      vector<int> vxx, vyy;
      rep(k, n) {
        rep(l, n) {
          if (k < i) {
            continue;
          }
          if (k == i && l < j) {
            continue;
          }
          if (tmpBoard[k][l] == num) {
            vxx.push_back(k);
            vyy.push_back(l);
          }
        }
      }

      int miniScore = INF;
      vector<int> miniVec;
      int keepX = x;
      int keepY = y;
      int miniX = -1;
      int miniY = -1;
      vector<vector<int>> keeptmpBoard;
      rep(k, vxx.size()) {
        x = keepX;
        y = keepY;
        int xx = vxx[k];
        int yy = vyy[k];
        vector<int> tmpVec;
        vector<vector<int>> tmptmpBoard = tmpBoard;
        // 適切な位置にピースを移動させる

        while (xx != i || yy != j) {
          if (yy != j) {
            if (yy < j) {
              route_move_cursor(x, y, xx, yy, xx, yy + 1, i, j, tmptmpBoard, tmpVec);
            }
            else {
              route_move_cursor(x, y, xx, yy, xx, yy - 1, i, j, tmptmpBoard, tmpVec);
            }
          }
          else if (xx != i) {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, i, j, tmptmpBoard, tmpVec);
          }
        }

        if (tmpVec.size() < miniScore || Rand() % 100 == 0) {
          miniScore = tmpVec.size();
          miniVec = tmpVec;
          keeptmpBoard = tmptmpBoard;
          miniX = x;
          miniY = y;
        }
      }

      for (auto& nd : miniVec) {
        ansDfs.push_back(nd);
      }
      x = miniX;
      y = miniY;
      tmpBoard = keeptmpBoard;
    }

    if (x == i && y == n - 1) {
      ansDfs.push_back(2);
      swap(tmpBoard[x][y], tmpBoard[x + 1][y]);
      x++;
    }

    // 最後の2個を揃える
    if (tmpBoard[i][n - 1] == aniBoard[i][n - 2]) {
      swap_vertical_pair(x, y, i, n - 2, tmpBoard, ansDfs);
    }
    else {
      int num = aniBoard[i][n - 2];
      vector<int> vxx, vyy;
      rep(k, n) {
        rep(l, n) {
          if (k <= i) {
            continue;
          }
          if (tmpBoard[k][l] == num) {
            vxx.push_back(k);
            vyy.push_back(l);
          }
        }
      }

      int miniScore = INF;
      vector<int> miniVec;
      int keepX = x;
      int keepY = y;
      int miniX = -1;
      int miniY = -1;
      vector<vector<int>> keeptmpBoard;
      rep(k, vxx.size()) {
        x = keepX;
        y = keepY;
        int xx = vxx[k];
        int yy = vyy[k];
        vector<int> tmpVec;
        vector<vector<int>> tmptmpBoard = tmpBoard;

        // 適切な位置にピースを移動させる
        while (xx != i + 1 || yy != n - 2) {
          if (yy != n - 2) {
            if (yy < n - 2) {
              route_move_cursor(x, y, xx, yy, xx, yy + 1, i, n - 1, tmptmpBoard, tmpVec);
            }
            else {
              route_move_cursor(x, y, xx, yy, xx, yy - 1, i, n - 1, tmptmpBoard, tmpVec);
            }
          }
          else if (xx != i + 1) {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, i, n - 1, tmptmpBoard, tmpVec);
          }
        }

        if (tmpVec.size() < miniScore || Rand() % 100 == 0) {
          miniScore = tmpVec.size();
          miniVec = tmpVec;
          keeptmpBoard = tmptmpBoard;
          miniX = x;
          miniY = y;
        }
      }

      for (auto& nd : miniVec) {
        ansDfs.push_back(nd);
      }
      x = miniX;
      y = miniY;
      tmpBoard = keeptmpBoard;

      fix_last_two_in_row(x, y, i, tmpBoard, ansDfs);
    }
  }

  // 下2列をそろえる
  rep(j, n - 2) {
    // 下のマスを上のマスの位置に持ってくる
    int num = aniBoard[n - 1][j];
    int xx = -1, yy = -1;
    rep(k, n) {
      rep(l, n) {
        if (k < n - 2) {
          continue;
        }
        if (l < j) {
          continue;
        }
        if (tmpBoard[k][l] == num) {
          xx = k;
          yy = l;
        }
      }
      if (xx != -1) {
        break;
      }
    }

    // 適切な位置にピースを移動させる
    while (xx != n - 2 || yy != j) {
      if (yy != j) {
        if (yy < j) {
          route_move_cursor(x, y, xx, yy, xx, yy + 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
        else {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      else if (xx != n - 2) {
        route_move_cursor(x, y, xx, yy, xx - 1, yy, n - 2, j, tmpBoard, ansDfs);
      }
    }

    if (x == n - 1 && y == j) {
      ansDfs.push_back(3);
      swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
      y++;
    }

    // 最後の2個を揃える
    if (tmpBoard[n - 1][j] == aniBoard[n - 2][j]) {
      swap_horizontal_pair(x, y, n - 2, j, tmpBoard, ansDfs);
    }
    else {
      int num = aniBoard[n - 2][j];
      int xx = -1, yy = -1;
      rep(k, n) {
        rep(l, n) {
          if (k < n - 2) {
            continue;
          }
          if (l <= j) {
            continue;
          }
          if (tmpBoard[k][l] == num) {
            xx = k;
            yy = l;
          }
        }
        if (xx != -1) {
          break;
        }
      }
      // 適切な位置にピースを移動させる
      while (xx != n - 2 || yy != j + 1) {
        if (xx != n - 2) {
          if (xx < n - 2) {
            route_move_cursor(x, y, xx, yy, xx + 1, yy, n - 2, j, tmpBoard, ansDfs, 2);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, n - 2, j, tmpBoard, ansDfs, 2);
          }
        }
        else if (yy != j + 1) {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, n - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      fix_last_two_in_col(x, y, j, tmpBoard, ansDfs);
    }
  }

  // カーソルを右下に持ってくる
  while (x != n - 1 || y != n - 1) {
    if (y != n - 1) {
      apply_move(x, y, 3, tmpBoard, ansDfs);
    }
    else {
      apply_move(x, y, 2, tmpBoard, ansDfs);
    }
  }

  if (!CheckAllDfs(tmpBoard)) {
    ansDfs.clear();
  }

  return ansDfs;
}

int exec_mode = 0; // 0: 標準出力, 1: ファイル出力
void output_data(int case_num) {
  if (exec_mode == 0) {
    // 標準出力
    rep(i, ans.size()) {
      cout << cc[ans[i]];
    }
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    rep(i, ans.size()) {
      ofs << cc[ans[i]];
    }
    ofs << endl;
  }
}

int Solve(int mode, int problemNum = 0) {
  start_timer();

  // 入力部
  Input(problemNum);

  // 木を1つ見つける
  bool isFind = false;
  {
    // 焼きなまし
    rep(_, 25) {
      bool isReset = true;
      if (FindTreeAni(isReset)) {
        isFind = true;
        break;
      }
    }
    rep(i, n) {
      rep(j, n) {
        dfsBoard[i][j] = aniBoard[i][j];
      }
    }
  }

  if (isFind) {
    // それぞれのピースに番号を付ける
    InitKindNumbers();
    init_piece_numbers();

    vector<int> ansDfs = FindAnsDfs2();
    ans = ansDfs;
    maxScore = CalcScore(ans);

    best_ans = ans;
    best_maxScore = maxScore;

    int loop = 0;
    double now_time = get_elapsed_time();
    while (true) {
      if (loop % 10 == 1) {
        now_time = get_elapsed_time();
        if (now_time > TL) break;
      }
      pair<P, P> pp[2];
      rep(i, 2) {
        pp[i] = shuffle_same_kind_piece();
      }

      vector<int> tmpAns = FindAnsDfs2();
      int tmpScore = CalcScore(tmpAns);

      int diffScore = tmpScore - maxScore;

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      double prob = exp((double)diffScore / temp);
      if (prob > Rand01()) {
        ans = tmpAns;
        maxScore = tmpScore;

        if (maxScore > best_maxScore) {
          best_maxScore = maxScore;
          best_ans = ans;
        }
      }
      else {
        // 元に戻す
        rep(i, 2) {
          swap(peaceNum[pp[i].first.first][pp[i].first.second], peaceNum[pp[i].second.first][pp[i].second.second]);
        }
      }
      loop++;
    }
  }
  else {
    // 愚直解
    {
      int x = startX;
      int y = startY;
      rep(i, t) {
        int val = Rand() % 4;
        while (!IsOK(x + dx[val], y + dy[val])) {
          val = Rand() % 4;
        }
        ans.push_back(val);
        x += dx[val];
        y += dy[val];
      }
    }

    maxScore = CalcScore(ans);

    best_ans = ans;
    best_maxScore = maxScore;

    // 山登り解、焼きなまし解
    double now_time = get_elapsed_time();
    int loop = 0;
    while (true) {
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
        if (now_time > TL) break;
      }

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      if (loop % 10 == 0) {
        op_shuffle_suffix(temp);
      }
      else {
        op_swap_adjacent(temp);
      }

      loop++;
    }

    // 最高スコアを戻す
    ans = best_ans;
    maxScore = best_maxScore;
  }

  // デバッグ用
  if (mode != 0) {
    cout << maxScore << endl;
    cout << get_elapsed_time() << "sec." << endl;
  }

  output_data(problemNum);

  return 0;
}

int main() {
  exec_mode = 10;

  if (exec_mode == 0) {
    Solve(exec_mode);
  }
  else if (exec_mode == 1) {
    Solve(exec_mode, 0);
  }
  else if (exec_mode == 10) {
    rep(_, 100) {
      Solve(exec_mode, _);
      ResetAll();
    }
  }

  return 0;
}

void Input(int problemNum) {
  std::ostringstream sout;
  sout << std::setfill('0') << std::setw(4) << problemNum;
  std::string numStr = sout.str();
  string fileNameIfs = "in/" + numStr + ".txt ";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    cin >> n >> t;
    rep(i, n) {
      string str;
      cin >> str;
      rep(j, n) {
        if ('0' <= str[j] && str[j] <= '9') {
          board[i][j] = str[j] - '0';
        }
        else {
          board[i][j] = str[j] - 'a' + 10;
        }
        if (board[i][j] == 0) {
          startX = i;
          startY = j;
        }
      }
    }
  }
  else {  // ファイル入力する
    ifs >> n >> t;
    rep(i, n) {
      string str;
      ifs >> str;
      rep(j, n) {
        if ('0' <= str[j] && str[j] <= '9') {
          board[i][j] = str[j] - '0';
        }
        else {
          board[i][j] = str[j] - 'a' + 10;
        }
        if (board[i][j] == 0) {
          startX = i;
          startY = j;
        }
      }
    }
  }

  rep(i, n) {
    rep(j, n) {
      cnt[board[i][j]]++;
    }
  }
}

void ResetAll() {
  ans.clear();
  best_ans.clear();
  rep(i, 16) {
    kindNumbers[i].clear();
    originNum[i].clear();
  }
}
