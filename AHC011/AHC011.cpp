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

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 100

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

///////////////////////////////////////////
// 上：2
// 下：8
// 左：1
// 右：4
//
///////////////////////////////////////////

const int INF = 1001001001;
const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U', 'L', 'D', 'R' };

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


  static double rand_unit()
  {
    return (rand32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* 変数 */
{
  // 入力用変数
  int board_size;
  int turn_limit;
  int board[10][10];
  int startX, startY;
  int kind_count[16];

  // 解答用変数
  ll cur_score;
  vector<int> route;
  int route_backup[2100];

  // 焼きなまし用変数
  double TIME_LIMIT = 2.9;
  double start_temp = 2048;
  double end_temp = 0.0001;
  ll best_score;
  vector<int> best_route;

}  // namespace

namespace /* Union Find*/
{
  int parent[MAX_N];   // 親
  int rank[MAX_N];  // 木の深さ
  int cnt[MAX_N];   // 属する頂点の個数(親のみ正しい)

  // n要素で初期化
  void init(const int n)
  {
    for (int i = 0; i < n; i++) {
      parent[i] = i;
      rank[i] = 0;
      cnt[i] = 1;
      int val = board[i / board_size][i % board_size];
      if (val == 1 || val == 2 || val == 4 || val == 8) {
        cnt[i] = 1;
      }
    }
  }

  // 木の根を求める
  int find(int x)
  {
    if (parent[x] == x) {
      return x;
    }
    else {
      return parent[x] = find(parent[x]);
    }
  }

  // xとyの属する集合を併合
  void unite(int x, int y)
  {
    x = find(x);
    y = find(y);
    if (x == y) { return; }

    if (rank[x] < rank[y]) {
      parent[x] = y;
      cnt[y] += cnt[x];
    }
    else {
      parent[y] = x;
      cnt[x] += cnt[y];
      if (rank[x] == rank[y]) rank[x]++;
    }
  } /* Union Find*/

  // xとyが同じ集合に属するか否か
  bool same(int x, int y)
  {
    return find(x) == find(y);
  }
}  // namespace

// スコア計算
int calc_board[10][10];
int calc_score(const vector<int>& ope)
{
  init(board_size * board_size);

  int x = startX, y = startY;
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      calc_board[i][j] = board[i][j];
    }
  }

  for (int i = 0; i < ope.size(); ++i) {
    swap(calc_board[x][y], calc_board[x + DX[ope[i]]][y + DY[ope[i]]]);
    x += DX[ope[i]];
    y += DY[ope[i]];
  }

  // 横の繋がり
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size - 1; ++j) {
      if ((calc_board[i][j] & (1 << 2)) && (calc_board[i][j + 1] & (1 << 0))) {
        unite(i * board_size + j, i * board_size + j + 1);
      }
    }
  }

  // 縦の繋がり
  for (int i = 0; i < board_size - 1; ++i) {
    for (int j = 0; j < board_size; ++j) {
      if ((calc_board[i][j] & (1 << 3)) && (calc_board[i + 1][j] & (1 << 1))) {
        unite(i * board_size + j, (i + 1) * board_size + j);
      }
    }
  }

  int res = 0;

  // スコア計算
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      int ij = i * board_size + j;
      if (find(ij) == ij) {
        if (cnt[ij] == board_size * board_size - 1) {
          res = 500000.0 * (2.0 - (double)ope.size() / turn_limit);
          return res;
        }
        else {
          res = max(res, (int)round(500000.0 * cnt[ij] / (board_size * board_size - 1)));
        }
      }
    }
  }

  return res;
}

void read_input(int case_num)
{
  std::ostringstream sout;
  sout << std::setfill('0') << std::setw(4) << case_num;
  std::string numStr = sout.str();
  string fileNameIfs = "in/" + numStr + ".txt ";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    cin >> board_size >> turn_limit;
    for (int i = 0; i < board_size; ++i) {
      string str;
      cin >> str;
      for (int j = 0; j < board_size; ++j) {
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
    ifs >> board_size >> turn_limit;
    for (int i = 0; i < board_size; ++i) {
      string str;
      ifs >> str;
      for (int j = 0; j < board_size; ++j) {
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

  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      kind_count[board[i][j]]++;
    }
  }
}

bool in_bounds(int x, int y)
{
  if (x < 0 || board_size <= x || y < 0 || board_size <= y) {
    return false;
  }
  return true;
}

bool route_in_bounds(const vector<int>& ope)
{
  int x = startX;
  int y = startY;
  for (int i = 0; i < ope.size(); ++i) {
    x += DX[ope[i]];
    y += DY[ope[i]];
    if (!in_bounds(x, y)) {
      return false;
    }
  }
  return true;
}

// kから後ろを全リセット
void shuffle_suffix(double temp)
{
  int ite = rand32() % turn_limit;

  int x = startX;
  int y = startY;
  for (int i = 0; i < turn_limit; ++i) {
    route_backup[i] = route[i];
    if (i < ite) {
      x += DX[route[i]];
      y += DY[route[i]];
    }
    else {
      int val = rand32() % 4;
      while (!in_bounds(x + DX[val], y + DY[val])) {
        val = rand32() % 4;
      }

      route[i] = val;
      x += DX[val];
      y += DY[val];
    }
  }

  int tmpScore = calc_score(route);

  int diffScore = tmpScore - cur_score;

  double prob = exp((double)diffScore / temp);
  if (prob > rand_unit()) {
    cur_score += diffScore;
    if (cur_score > best_score) {
      best_score = cur_score;
      best_route = route;
    }
  }
  else {
    // 元に戻す
    for (int i = 0; i < turn_limit; ++i) {
      route[i] = route_backup[i];
    }
  }
}

// kとkの直後をスワップ
void swap_adjacent(double temp)
{
  int ite = rand32() % (turn_limit - 1);

  swap(route[ite], route[ite + 1]);

  if (!route_in_bounds(route)) {
    swap(route[ite], route[ite + 1]);
    return;
  }

  int tmpScore = calc_score(route);

  int diffScore = tmpScore - cur_score;

  double prob = exp((double)diffScore / temp);
  if (prob > rand_unit()) {
    cur_score += diffScore;
    if (cur_score > best_score) {
      best_score = cur_score;
      best_route = route;
    }
  }
  else {
    // 元に戻す
    swap(route[ite], route[ite + 1]);
  }
}

// 木を全探索
namespace
{
  int tree_board[10][10];
  bool CheckAllDfs(const vector<vector<int>>& vec)
  {
    init(board_size * board_size);
    // 横の繋がり
    for (int i = 0; i < board_size; ++i) {
      for (int j = 0; j < board_size - 1; ++j) {
        if ((vec[i][j] & (1 << 2)) && (vec[i][j + 1] & (1 << 0))) {
          unite(i * board_size + j, i * board_size + j + 1);
        }
      }
    }

    // 縦の繋がり
    for (int i = 0; i < board_size - 1; ++i) {
      for (int j = 0; j < board_size; ++j) {
        if ((vec[i][j] & (1 << 3)) && (vec[i + 1][j] & (1 << 1))) {
          unite(i * board_size + j, (i + 1) * board_size + j);
        }
      }
    }

    if (cnt[find(0)] == board_size * board_size - 1) {
      return true;
    }
    return false;
  }
}  // namespace

// 木を焼きなましで見つける
int anneal_board[10][10];
int best_anneal_board[10][10];
int calc_tree_score()
{
  init(board_size * board_size);
  // 横の繋がり
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size - 1; ++j) {
      if ((anneal_board[i][j] & (1 << 2)) && (anneal_board[i][j + 1] & (1 << 0))) {
        unite(i * board_size + j, i * board_size + j + 1);
      }
    }
  }

  // 縦の繋がり
  for (int i = 0; i < board_size - 1; ++i) {
    for (int j = 0; j < board_size; ++j) {
      if ((anneal_board[i][j] & (1 << 3)) && (anneal_board[i + 1][j] & (1 << 1))) {
        unite(i * board_size + j, (i + 1) * board_size + j);
      }
    }
  }

  int MaxSize = 0;
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      MaxSize = max(MaxSize, cnt[find(i * board_size + j)]);
    }
    if (MaxSize >= board_size * board_size / 2) {
      break;
    }
  }

  return MaxSize;
}

bool find_tree_anneal(bool isReset = false)
{
  Timer timer;
  timer.start();
  const double ANNEAL_TIME_LIMIT = 0.1;
  double elapsed = timer.get_elapsed_time();
  int loop = 0;
  double start_temp_anneal = 0.1;
  double end_temp_anneal = 0.0;

  if (isReset) {
    for (int i = 0; i < board_size; ++i) {
      for (int j = 0; j < board_size; ++j) {
        anneal_board[i][j] = board[i][j];
      }
    }
  }

  // カーソルは右下固定
  swap(anneal_board[startX][startY], anneal_board[board_size - 1][board_size - 1]);

  int max_score = calc_tree_score();
  int best_max_score = max_score;
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      best_anneal_board[i][j] = anneal_board[i][j];
    }
  }

  while (true) {
    if (loop % 100 == 1) {
      elapsed = timer.get_elapsed_time();
      if (elapsed > ANNEAL_TIME_LIMIT) {
        break;
      }
    }

    int x1 = rand32() % board_size;
    int y1 = rand32() % board_size;
    int x2 = rand32() % board_size;
    while (x1 == x2) {
      x2 = rand32() % board_size;
    }
    int y2 = rand32() % board_size;
    while (y1 == y2) {
      y2 = rand32() % board_size;
    }

    if (x1 == board_size - 1 && y1 == board_size - 1) {
      continue;
    }
    if (x2 == board_size - 1 && y2 == board_size - 1) {
      continue;
    }

    swap(anneal_board[x1][y1], anneal_board[x2][y2]);
    int newPoint = calc_tree_score();

    int diffScore = newPoint - max_score;

    double temp = start_temp_anneal + (end_temp_anneal - start_temp_anneal) * elapsed / ANNEAL_TIME_LIMIT;
    double prob = exp((double)diffScore / temp);
    if (prob > rand_unit()) {
      max_score += diffScore;
      if (max_score > best_max_score) {
        best_max_score = max_score;
        for (int i = 0; i < board_size; ++i) {
          for (int j = 0; j < board_size; ++j) {
            best_anneal_board[i][j] = anneal_board[i][j];
          }
        }
      }
    }
    else {
      // 元に戻す
      swap(anneal_board[x1][y1], anneal_board[x2][y2]);
    }

    loop++;
    if (best_max_score == board_size * board_size - 1) {
      break;
    }
  }

  max_score = best_max_score;
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      anneal_board[i][j] = best_anneal_board[i][j];
    }
  }

  if (max_score == board_size * board_size - 1) {
    return true;
  }
  return false;
}

int piece_number[10][10];
// 作成した盤面の転倒数をチェックする
bool is_even_inversion()
{
  int tmp_board[10][10];
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      tmp_board[i][j] = piece_number[i][j];
    }
  }
  int cnt = 0;

  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      int num = i * board_size + j;
      int x = -1;
      int y = -1;
      for (int k = 0; k < board_size; ++k) {
        for (int l = 0; l < board_size; ++l) {
          if (tmp_board[k][l] == num) {
            x = k;
            y = j;
          }
        }
      }
      while (y < i) {
        cnt++;
        swap(tmp_board[x][y], tmp_board[x][y + 1]);
        y++;
      }
      while (y > i) {
        cnt++;
        swap(tmp_board[x][y], tmp_board[x][y - 1]);
        y--;
      }
      while (x < i) {
        cnt++;
        swap(tmp_board[x][y], tmp_board[x + 1][y]);
        x++;
      }
      while (x > i) {
        cnt++;
        swap(tmp_board[x][y], tmp_board[x - 1][y]);
        x--;
      }
    }
  }

  if (cnt % 2 == 0) {
    return true;
  }
  return false;
}

// 作成した木からピースの種類ごとの番号を決定する
vector<int> kindNumbers[16];
vector<P> origin_positions[16];
void init_kind_indices()
{
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      kindNumbers[tree_board[i][j]].push_back(i * board_size + j);
      origin_positions[board[i][j]].push_back(P(i, j));
    }
  }
}

void reset_state()
{
  route.clear();
  best_route.clear();
  for (int i = 0; i < 16; ++i) {
    kindNumbers[i].clear();
    origin_positions[i].clear();
  }
}

// ピースに番号を振る
void init_piece_numbers()
{
  int ite[16] = {};
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      piece_number[i][j] = kindNumbers[board[i][j]][ite[board[i][j]]];
      ite[board[i][j]]++;
    }
  }

  // もし転倒数が奇数なら1箇所スワップする
  if (!is_even_inversion()) {
    for (int i = 0; i < 16; ++i) {
      if (kindNumbers[i].size() >= 2) {
        int num1 = kindNumbers[i][0];
        int num2 = kindNumbers[i][1];
        int x1, y1, x2, y2;
        for (int j = 0; j < board_size; ++j) {
          for (int k = 0; k < board_size; ++k) {
            if (piece_number[j][k] == num1) {
              x1 = j;
              y1 = k;
            }
            if (piece_number[j][k] == num2) {
              x2 = j;
              y2 = k;
            }
          }
        }
        swap(piece_number[x1][y1], piece_number[x2][y2]);
        break;
      }
    }
  }
}

pair<P, P> shuffle_same_kind_piece()
{
  int ite = rand32() % 16;
  while (origin_positions[ite].size() <= 1) {
    ite = rand32() % 16;
  }

  int a = rand32() % origin_positions[ite].size();
  int b = rand32() % origin_positions[ite].size();
  while (a == b) {
    b = rand32() % origin_positions[ite].size();
  }

  int x1 = origin_positions[ite][a].first;
  int y1 = origin_positions[ite][a].second;
  int x2 = origin_positions[ite][b].first;
  int y2 = origin_positions[ite][b].second;
  swap(piece_number[x1][y1], piece_number[x2][y2]);

  return pair<P, P>({ {x1, y1}, {x2, y2} });
}

// 木を作成する手順を1つ作成する
void apply_move(int& x, int& y, int nd, vector<vector<int>>& tmpBoard, vector<int>& route_tmp)
{
  route_tmp.push_back(nd);
  swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
  x += DX[nd];
  y += DY[nd];
}

void route_move_cursor(int& x, int& y, int& xx, int& yy, int tx, int ty, int ii, int jj, vector<vector<int>>& tmpBoard, vector<int>& answer, int mode = 0)
{
  int dp[10][10];
  int dir[10][10];
  if (mode == 0) {
    for (int i = 0; i < board_size; ++i) {
      for (int j = 0; j < board_size; ++j) {
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
    for (int i = 0; i < board_size; ++i) {
      for (int j = 0; j < board_size; ++j) {
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
    bool finished = false;
    for (int i = 0; i < 4; ++i) {
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == tx && nb == ty) {
          finished = true;
          break;
        }
      }
    }
    if (finished) {
      break;
    }
  }

  vector<int> reverse_path;
  int revX = tx, revY = ty;
  while (revX != x || revY != y) {
    int rev_dir = dir[revX][revY];
    reverse_path.push_back(rev_dir);
    revX -= DX[rev_dir];
    revY -= DY[rev_dir];
  }

  reverse(reverse_path.begin(), reverse_path.end());

  for (auto nd : reverse_path) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
  for (int i = 0; i < 4; ++i) {
    if (x + DX[i] == xx && y + DY[i] == yy) {
      swap(tmpBoard[x][y], tmpBoard[xx][yy]);
      swap(x, xx);
      swap(y, yy);
      answer.push_back(i);
      break;
    }
  }
}

// 3×2マスを使って上にマスを入れ替える
// (xx,yy) = (左上,右上)
void swap_vertical_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& answer)
{
  while (y != yy + 1) {
    answer.push_back(3);
    swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
    y++;
  }
  // 上左下右
  vector<int> order = { 0, 1, 2, 3, 2, 1, 0, 0, 3, 2, 1, 2, 3, 0, 0, 1, 2 };
  for (auto nd : order) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

// 2×3マスを使って上にマスを入れ替える
// (xx,yy) = (左上,右上)
void swap_horizontal_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& answer)
{
  while (x != board_size - 1) {
    apply_move(x, y, 2, tmpBoard, answer);
  }
  // 上左下右
  vector<int> order = { 1, 0, 3, 2, 3, 0, 1, 1, 2, 3, 0, 3, 2, 1, 1, 0, 3 };
  for (auto nd : order) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

void fix_last_two_in_row(int& x, int& y, int ii, vector<vector<int>>& tmpBoard, vector<int>& answer)
{
  // カーソル移動
  int dp[10][10];
  int dir[10][10];
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
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
  dp[ii + 1][board_size - 2] = -1;
  dp[ii][board_size - 1] = INF;
  dir[x][y] = -1;
  while (que.size()) {
    int a = que.front().first;
    int b = que.front().second;
    que.pop();
    bool finished = false;
    for (int i = 0; i < 4; ++i) {
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == ii && nb == board_size - 1) {
          finished = true;
          break;
        }
      }
    }
    if (finished) {
      break;
    }
  }

  vector<int> reverse_path;
  int revX = ii, revY = board_size - 1;
  while (revX != x || revY != y) {
    int rev_dir = dir[revX][revY];
    reverse_path.push_back(rev_dir);
    revX -= DX[rev_dir];
    revY -= DY[rev_dir];
  }

  reverse(reverse_path.begin(), reverse_path.end());

  for (auto nd : reverse_path) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }

  // 上左下右
  vector<int> order = { 1, 2 };
  for (auto nd : order) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

void fix_last_two_in_col(int& x, int& y, int jj, vector<vector<int>>& tmpBoard, vector<int>& answer)
{
  // カーソル移動
  int dp[10][10];
  int dir[10][10];
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      if (i < board_size - 2) {
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
  dp[board_size - 2][jj + 1] = -1;
  dp[board_size - 1][jj] = INF;
  dir[x][y] = -1;
  while (que.size()) {
    int a = que.front().first;
    int b = que.front().second;
    que.pop();
    bool finished = false;
    for (int i = 0; i < 4; ++i) {
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == board_size - 1 && nb == jj) {
          finished = true;
          break;
        }
      }
    }
    if (finished) {
      break;
    }
  }

  vector<int> reverse_path;
  int revX = board_size - 1, revY = jj;
  while (revX != x || revY != y) {
    int rev_dir = dir[revX][revY];
    reverse_path.push_back(rev_dir);
    revX -= DX[rev_dir];
    revY -= DY[rev_dir];
  }

  reverse(reverse_path.begin(), reverse_path.end());

  for (auto nd : reverse_path) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }

  // 上左下右
  vector<int> order = { 0, 3 };
  for (auto nd : order) {
    answer.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

vector<int> build_route()
{
  vector<vector<int>> tmpBoard(10, vector<int>(10));
  for (int i = 0; i < board_size; ++i) {
    for (int j = 0; j < board_size; ++j) {
      tmpBoard[i][j] = board[i][j];
    }
  }

  vector<int> answer;

  int x = startX;
  int y = startY;

  // 1列目から下から3列目まで完成させる
  for (int i = 0; i < board_size - 2; ++i) {
    for (int j = 0; j < board_size - 1; ++j) {
      int num = anneal_board[i][j];
      if (j == board_size - 2) {
        num = anneal_board[i][board_size - 1];
      }

      vector<int> xs, ys;
      for (int k = 0; k < board_size; ++k) {
        for (int l = 0; l < board_size; ++l) {
          if (k < i) {
            continue;
          }
          if (k == i && l < j) {
            continue;
          }
          if (tmpBoard[k][l] == num) {
            xs.push_back(k);
            ys.push_back(l);
          }
        }
      }

      int min_score = INF;
      vector<int> min_moves;
      int keepX = x;
      int keepY = y;
      int minX = -1;
      int minY = -1;
      vector<vector<int>> keeptmpBoard;
      for (int k = 0; k < xs.size(); ++k) {
        x = keepX;
        y = keepY;
        int xx = xs[k];
        int yy = ys[k];
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

        if (tmpVec.size() < min_score || rand32() % 100 == 0) {
          min_score = tmpVec.size();
          min_moves = tmpVec;
          keeptmpBoard = tmptmpBoard;
          minX = x;
          minY = y;
        }
      }

      for (auto& nd : min_moves) {
        answer.push_back(nd);
      }
      x = minX;
      y = minY;
      tmpBoard = keeptmpBoard;
    }

    if (x == i && y == board_size - 1) {
      answer.push_back(2);
      swap(tmpBoard[x][y], tmpBoard[x + 1][y]);
      x++;
    }

    // 最後の2個を揃える
    if (tmpBoard[i][board_size - 1] == anneal_board[i][board_size - 2]) {
      swap_vertical_pair(x, y, i, board_size - 2, tmpBoard, answer);
    }
    else {
      int num = anneal_board[i][board_size - 2];
      vector<int> xs, ys;
      for (int k = 0; k < board_size; ++k) {
        for (int l = 0; l < board_size; ++l) {
          if (k <= i) {
            continue;
          }
          if (tmpBoard[k][l] == num) {
            xs.push_back(k);
            ys.push_back(l);
          }
        }
      }

      int min_score = INF;
      vector<int> min_moves;
      int keepX = x;
      int keepY = y;
      int minX = -1;
      int minY = -1;
      vector<vector<int>> keeptmpBoard;
      for (int k = 0; k < xs.size(); ++k) {
        x = keepX;
        y = keepY;
        int xx = xs[k];
        int yy = ys[k];
        vector<int> tmpVec;
        vector<vector<int>> tmptmpBoard = tmpBoard;

        // 適切な位置にピースを移動させる
        while (xx != i + 1 || yy != board_size - 2) {
          if (yy != board_size - 2) {
            if (yy < board_size - 2) {
              route_move_cursor(x, y, xx, yy, xx, yy + 1, i, board_size - 1, tmptmpBoard, tmpVec);
            }
            else {
              route_move_cursor(x, y, xx, yy, xx, yy - 1, i, board_size - 1, tmptmpBoard, tmpVec);
            }
          }
          else if (xx != i + 1) {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, i, board_size - 1, tmptmpBoard, tmpVec);
          }
        }

        if (tmpVec.size() < min_score || rand32() % 100 == 0) {
          min_score = tmpVec.size();
          min_moves = tmpVec;
          keeptmpBoard = tmptmpBoard;
          minX = x;
          minY = y;
        }
      }

      for (auto& nd : min_moves) {
        answer.push_back(nd);
      }
      x = minX;
      y = minY;
      tmpBoard = keeptmpBoard;

      fix_last_two_in_row(x, y, i, tmpBoard, answer);
    }
  }

  // 下2列をそろえる
  for (int j = 0; j < board_size - 2; ++j) {
    // 下のマスを上のマスの位置に持ってくる
    int num = anneal_board[board_size - 1][j];
    int xx = -1, yy = -1;
    for (int k = 0; k < board_size; ++k) {
      for (int l = 0; l < board_size; ++l) {
        if (k < board_size - 2) {
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
    while (xx != board_size - 2 || yy != j) {
      if (yy != j) {
        if (yy < j) {
          route_move_cursor(x, y, xx, yy, xx, yy + 1, board_size - 2, j, tmpBoard, answer, 2);
        }
        else {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, board_size - 2, j, tmpBoard, answer, 2);
        }
      }
      else if (xx != board_size - 2) {
        route_move_cursor(x, y, xx, yy, xx - 1, yy, board_size - 2, j, tmpBoard, answer);
      }
    }

    if (x == board_size - 1 && y == j) {
      answer.push_back(3);
      swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
      y++;
    }

    // 最後の2個を揃える
    if (tmpBoard[board_size - 1][j] == anneal_board[board_size - 2][j]) {
      swap_horizontal_pair(x, y, board_size - 2, j, tmpBoard, answer);
    }
    else {
      int num = anneal_board[board_size - 2][j];
      int xx = -1, yy = -1;
      for (int k = 0; k < board_size; ++k) {
        for (int l = 0; l < board_size; ++l) {
          if (k < board_size - 2) {
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
      while (xx != board_size - 2 || yy != j + 1) {
        if (xx != board_size - 2) {
          if (xx < board_size - 2) {
            route_move_cursor(x, y, xx, yy, xx + 1, yy, board_size - 2, j, tmpBoard, answer, 2);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, board_size - 2, j, tmpBoard, answer, 2);
          }
        }
        else if (yy != j + 1) {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, board_size - 2, j, tmpBoard, answer, 2);
        }
      }
      fix_last_two_in_col(x, y, j, tmpBoard, answer);
    }
  }

  // カーソルを右下に持ってくる
  while (x != board_size - 1 || y != board_size - 1) {
    if (y != board_size - 1) {
      apply_move(x, y, 3, tmpBoard, answer);
    }
    else {
      apply_move(x, y, 2, tmpBoard, answer);
    }
  }

  if (!CheckAllDfs(tmpBoard)) {
    answer.clear();
  }

  return answer;
}

int mode = 0; // 0: 標準出力, 1: ファイル出力
void output_data(int case_num)
{
  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < route.size(); ++i) {
      cout << DIR_CHAR[route[i]];
    }
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < route.size(); ++i) {
      ofs << DIR_CHAR[route[i]];
    }
    ofs << endl;
  }
}

int solve_case(int mode, int case_num = 0)
{
  Timer timer;
  timer.start();

  // 入力部
  read_input(case_num);

  // 木を1つ見つける
  bool found = false;
  {
    // 焼きなまし
    for (int _ = 0; _ < 25; ++_) {
      bool isReset = true;
      if (find_tree_anneal(isReset)) {
        found = true;
        break;
      }
    }
    for (int i = 0; i < board_size; ++i) {
      for (int j = 0; j < board_size; ++j) {
        tree_board[i][j] = anneal_board[i][j];
      }
    }
  }

  if (found) {
    // それぞれのピースに番号を付ける
    init_kind_indices();
    init_piece_numbers();

    vector<int> answer = build_route();
    route = answer;
    cur_score = calc_score(route);

    best_route = route;
    best_score = cur_score;

    int loop = 0;
    double now_time = timer.get_elapsed_time();
    while (true) {
      if (loop % 10 == 1) {
        now_time = timer.get_elapsed_time();
        if (now_time > TIME_LIMIT) {
          break;
        }
      }
      pair<P, P> swap_pairs[2];
      for (int i = 0; i < 2; ++i) {
        swap_pairs[i] = shuffle_same_kind_piece();
      }

      vector<int> tmpAns = build_route();
      int tmpScore = calc_score(tmpAns);

      int diffScore = tmpScore - cur_score;

      double temp = start_temp + (end_temp - start_temp) * now_time / TIME_LIMIT;
      double prob = exp((double)diffScore / temp);
      if (prob > rand_unit()) {
        route = tmpAns;
        cur_score = tmpScore;

        if (cur_score > best_score) {
          best_score = cur_score;
          best_route = route;
        }
      }
      else {
        // 元に戻す
        for (int i = 0; i < 2; ++i) {
          swap(piece_number[swap_pairs[i].first.first][swap_pairs[i].first.second], piece_number[swap_pairs[i].second.first][swap_pairs[i].second.second]);
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
      for (int i = 0; i < turn_limit; ++i) {
        int val = rand32() % 4;
        while (!in_bounds(x + DX[val], y + DY[val])) {
          val = rand32() % 4;
        }
        route.push_back(val);
        x += DX[val];
        y += DY[val];
      }
    }

    cur_score = calc_score(route);

    best_route = route;
    best_score = cur_score;

    // 山登り解、焼きなまし解
    double now_time = timer.get_elapsed_time();
    int loop = 0;
    while (true) {
      if (loop % 100 == 1) {
        now_time = timer.get_elapsed_time();
        if (now_time > TIME_LIMIT) {
          break;
        }
      }

      double temp = start_temp + (end_temp - start_temp) * now_time / TIME_LIMIT;
      if (loop % 10 == 0) {
        shuffle_suffix(temp);
      }
      else {
        swap_adjacent(temp);
      }

      loop++;
    }

    // 最高スコアを戻す
    route = best_route;
    cur_score = best_score;
  }

  // デバッグ用
  if (mode != 0) {
    cout << cur_score << endl;
    cout << timer.get_elapsed_time() << "sec." << endl;
  }

  output_data(case_num);

  return 0;
}

int main()
{
  mode = 10;

  if (mode == 0) {
    solve_case(mode);
  }
  else if (mode == 1) {
    solve_case(mode, 0);
  }
  else if (mode == 10) {
    for (int _ = 4; _ < 5; ++_) {
      solve_case(mode, _);
      reset_state();
    }
  }

  return 0;
}
