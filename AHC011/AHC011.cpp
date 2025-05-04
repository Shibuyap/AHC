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

// �^�C�}�[
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
// ��F2
// ���F8
// ���F1
// �E�F4
//
///////////////////////////////////////////

const int INF = 1001001001;
const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U', 'L', 'D', 'R' };

namespace /* �������C�u���� */
{
  static uint32_t rand32() {
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


  static double rand_unit() {
    return (rand32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* �ϐ� */
{
  // ���͗p�ϐ�
  int board_size;
  int turn_limit;
  int board[10][10];
  int startX, startY;
  int kind_count[16];

  // �𓚗p�ϐ�
  ll cur_score;
  vector<int> route;
  int route_backup[2100];

  // �Ă��Ȃ܂��p�ϐ�
  double TL = 2.9;
  double start_temp = 2048;
  double end_temp = 0.0001;
  ll best_score;
  vector<int> best_route;

}  // namespace

namespace /* Union Find*/
{
  int parUF[MAX_N];   // �e
  int rankUF[MAX_N];  // �؂̐[��
  int cntUF[MAX_N];   // �����钸�_�̌�(�e�̂ݐ�����)

  // n�v�f�ŏ�����
  void UFinit(const int nn) {
    for (int i = 0; i < nn; i++) {
      parUF[i] = i;
      rankUF[i] = 0;
      cntUF[i] = 1;
      int val = board[i / board_size][i % board_size];
      if (val == 1 || val == 2 || val == 4 || val == 8) {
        cntUF[i] = 1;
      }
    }
  }

  // �؂̍������߂�
  int findUF(int x) {
    if (parUF[x] == x) {
      return x;
    }
    else {
      return parUF[x] = findUF(parUF[x]);
    }
  }

  // x��y�̑�����W���𕹍�
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

  // x��y�������W���ɑ����邩�ۂ�
  bool sameUF(int x, int y) {
    return findUF(x) == findUF(y);
  }
}  // namespace

// �X�R�A�v�Z
int boardForCalc[10][10];
int calc_score(const vector<int>& ope) {
  UFinit(board_size * board_size);

  int x = startX, y = startY;
  rep(i, board_size) {
    rep(j, board_size) {
      boardForCalc[i][j] = board[i][j];
    }
  }

  rep(i, ope.size()) {
    swap(boardForCalc[x][y], boardForCalc[x + DX[ope[i]]][y + DY[ope[i]]]);
    x += DX[ope[i]];
    y += DY[ope[i]];
  }

  // ���̌q����
  rep(i, board_size) {
    rep(j, board_size - 1) {
      if ((boardForCalc[i][j] & (1 << 2)) && (boardForCalc[i][j + 1] & (1 << 0))) {
        uniteUF(i * board_size + j, i * board_size + j + 1);
      }
    }
  }

  // �c�̌q����
  rep(i, board_size - 1) {
    rep(j, board_size) {
      if ((boardForCalc[i][j] & (1 << 3)) && (boardForCalc[i + 1][j] & (1 << 1))) {
        uniteUF(i * board_size + j, (i + 1) * board_size + j);
      }
    }
  }

  int res = 0;

  // �X�R�A�v�Z
  rep(i, board_size) {
    rep(j, board_size) {
      int ij = i * board_size + j;
      if (findUF(ij) == ij) {
        if (cntUF[ij] == board_size * board_size - 1) {
          res = 500000.0 * (2.0 - (double)ope.size() / turn_limit);
          return res;
        }
        else {
          res = max(res, (int)round(500000.0 * cntUF[ij] / (board_size * board_size - 1)));
        }
      }
    }
  }

  return res;
}

void read_input(int problemNum) {
  std::ostringstream sout;
  sout << std::setfill('0') << std::setw(4) << problemNum;
  std::string numStr = sout.str();
  string fileNameIfs = "in/" + numStr + ".txt ";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // �W�����͂���
    cin >> board_size >> turn_limit;
    rep(i, board_size) {
      string str;
      cin >> str;
      rep(j, board_size) {
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
  else {  // �t�@�C�����͂���
    ifs >> board_size >> turn_limit;
    rep(i, board_size) {
      string str;
      ifs >> str;
      rep(j, board_size) {
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

  rep(i, board_size) {
    rep(j, board_size) {
      kind_count[board[i][j]]++;
    }
  }
}

bool in_bounds(int x, int y) {
  if (x < 0 || board_size <= x || y < 0 || board_size <= y) {
    return false;
  }
  return true;
}

bool route_in_bounds(const vector<int>& ope) {
  int x = startX;
  int y = startY;
  rep(i, ope.size()) {
    x += DX[ope[i]];
    y += DY[ope[i]];
    if (!in_bounds(x, y)) {
      return false;
    }
  }
  return true;
}

// k�������S���Z�b�g
void shuffle_suffix(double temp) {
  int ite = rand32() % turn_limit;

  int x = startX;
  int y = startY;
  rep(i, turn_limit) {
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
    // ���ɖ߂�
    rep(i, turn_limit) {
      route[i] = route_backup[i];
    }
  }
}

// k��k�̒�����X���b�v
void swap_adjacent(double temp) {
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
    // ���ɖ߂�
    swap(route[ite], route[ite + 1]);
  }
}

// �؂�S�T��
namespace
{
  int dfsBoard[10][10];
  bool CheckAllDfs(const vector<vector<int>>& vec) {
    UFinit(board_size * board_size);
    // ���̌q����
    rep(i, board_size) {
      rep(j, board_size - 1) {
        if ((vec[i][j] & (1 << 2)) && (vec[i][j + 1] & (1 << 0))) {
          uniteUF(i * board_size + j, i * board_size + j + 1);
        }
      }
    }

    // �c�̌q����
    rep(i, board_size - 1) {
      rep(j, board_size) {
        if ((vec[i][j] & (1 << 3)) && (vec[i + 1][j] & (1 << 1))) {
          uniteUF(i * board_size + j, (i + 1) * board_size + j);
        }
      }
    }

    if (cntUF[findUF(0)] == board_size * board_size - 1) {
      return true;
    }
    return false;
  }
}  // namespace

// �؂��Ă��Ȃ܂��Ō�����
int aniBoard[10][10];
int bestMaxAniBoard[10][10];
int calc_anneal_score() {
  UFinit(board_size * board_size);
  // ���̌q����
  rep(i, board_size) {
    rep(j, board_size - 1) {
      if ((aniBoard[i][j] & (1 << 2)) && (aniBoard[i][j + 1] & (1 << 0))) {
        uniteUF(i * board_size + j, i * board_size + j + 1);
      }
    }
  }

  // �c�̌q����
  rep(i, board_size - 1) {
    rep(j, board_size) {
      if ((aniBoard[i][j] & (1 << 3)) && (aniBoard[i + 1][j] & (1 << 1))) {
        uniteUF(i * board_size + j, (i + 1) * board_size + j);
      }
    }
  }

  int MaxSize = 0;
  rep(i, board_size) {
    rep(j, board_size) {
      MaxSize = max(MaxSize, cntUF[findUF(i * board_size + j)]);
    }
    if (MaxSize >= board_size * board_size / 2) {
      break;
    }
  }

  return MaxSize;
}

bool anneal_find_tree(bool isReset = false) {
  clock_t startAniTime, endAniTime;
  const double AniTL = 0.1;
  startAniTime = clock();
  endAniTime = clock();
  double nowAniTime = (double)(endAniTime - startAniTime) / CLOCKS_PER_SEC;
  int loopAni = 0;
  double startAniTemp = 0.1;
  double endAniTemp = 0.0;

  if (isReset) {
    rep(i, board_size) {
      rep(j, board_size) {
        aniBoard[i][j] = board[i][j];
      }
    }
  }

  // �J�[�\���͉E���Œ�
  swap(aniBoard[startX][startY], aniBoard[board_size - 1][board_size - 1]);

  int maxAniScore = calc_anneal_score();
  int bestMaxAniScore = maxAniScore;
  rep(i, board_size) {
    rep(j, board_size) {
      bestMaxAniBoard[i][j] = aniBoard[i][j];
    }
  }

  while (true) {
    if (loopAni % 100 == 1) {
      endAniTime = clock();
      nowAniTime = (double)(endAniTime - startAniTime) / CLOCKS_PER_SEC;
      if (nowAniTime > AniTL) break;
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

    swap(aniBoard[x1][y1], aniBoard[x2][y2]);
    int newPoint = calc_anneal_score();

    int diffScore = newPoint - maxAniScore;

    double temp = startAniTemp + (endAniTemp - startAniTemp) * nowAniTime / AniTL;
    double prob = exp((double)diffScore / temp);
    if (prob > rand_unit()) {
      maxAniScore += diffScore;
      if (maxAniScore > bestMaxAniScore) {
        bestMaxAniScore = maxAniScore;
        rep(i, board_size) {
          rep(j, board_size) {
            bestMaxAniBoard[i][j] = aniBoard[i][j];
          }
        }
      }
    }
    else {
      // ���ɖ߂�
      swap(aniBoard[x1][y1], aniBoard[x2][y2]);
    }

    loopAni++;
    if (bestMaxAniScore == board_size * board_size - 1) {
      break;
    }
  }

  maxAniScore = bestMaxAniScore;
  rep(i, board_size) {
    rep(j, board_size) {
      aniBoard[i][j] = bestMaxAniBoard[i][j];
    }
  }

  if (maxAniScore == board_size * board_size - 1) {
    return true;
  }
  return false;
}

int piece_number[10][10];
// �쐬�����Ֆʂ̓]�|�����`�F�b�N����
bool is_even_inversion() {
  int tmp_board[10][10];
  rep(i, board_size) {
    rep(j, board_size) {
      tmp_board[i][j] = piece_number[i][j];
    }
  }
  int cnt = 0;

  rep(i, board_size) {
    rep(j, board_size) {
      int num = i * board_size + j;
      int x = -1;
      int y = -1;
      rep(k, board_size) {
        rep(l, board_size) {
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

// �쐬�����؂���s�[�X�̎�ނ��Ƃ̔ԍ������肷��
vector<int> kindNumbers[16];
vector<P> origin_positions[16];
void init_kind_indices() {
  rep(i, board_size) {
    rep(j, board_size) {
      kindNumbers[dfsBoard[i][j]].push_back(i * board_size + j);
      origin_positions[board[i][j]].push_back(P(i, j));
    }
  }
}

void reset_state() {
  route.clear();
  best_route.clear();
  rep(i, 16) {
    kindNumbers[i].clear();
    origin_positions[i].clear();
  }
}

// �s�[�X�ɔԍ���U��
void init_piece_numbers() {
  int ite[16] = {};
  rep(i, board_size) {
    rep(j, board_size) {
      piece_number[i][j] = kindNumbers[board[i][j]][ite[board[i][j]]];
      ite[board[i][j]]++;
    }
  }

  // �����]�|������Ȃ�1�ӏ��X���b�v����
  if (!is_even_inversion()) {
    rep(i, 16) {
      if (kindNumbers[i].size() >= 2) {
        int num1 = kindNumbers[i][0];
        int num2 = kindNumbers[i][1];
        int x1, y1, x2, y2;
        rep(j, board_size) {
          rep(k, board_size) {
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

pair<P, P> shuffle_same_kind_piece() {
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

// �؂��쐬����菇��1�쐬����
void apply_move(int& x, int& y, int nd, vector<vector<int>>& tmpBoard, vector<int>& route_tmp) {
  route_tmp.push_back(nd);
  swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
  x += DX[nd];
  y += DY[nd];
}

void route_move_cursor(int& x, int& y, int& xx, int& yy, int tx, int ty, int ii, int jj, vector<vector<int>>& tmpBoard, vector<int>& ansDfs, int mode = 0) {
  int dp[10][10];
  int dir[10][10];
  if (mode == 0) {
    rep(i, board_size) {
      rep(j, board_size) {
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
    rep(i, board_size) {
      rep(j, board_size) {
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
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
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
    revX -= DX[revD];
    revY -= DY[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
  rep(i, 4) {
    if (x + DX[i] == xx && y + DY[i] == yy) {
      swap(tmpBoard[x][y], tmpBoard[xx][yy]);
      swap(x, xx);
      swap(y, yy);
      ansDfs.push_back(i);
      break;
    }
  }
}

// 3�~2�}�X���g���ď�Ƀ}�X�����ւ���
// (xx,yy) = (����,�E��)
void swap_vertical_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  while (y != yy + 1) {
    ansDfs.push_back(3);
    swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
    y++;
  }
  // �㍶���E
  vector<int> order = { 0, 1, 2, 3, 2, 1, 0, 0, 3, 2, 1, 2, 3, 0, 0, 1, 2 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

// 2�~3�}�X���g���ď�Ƀ}�X�����ւ���
// (xx,yy) = (����,�E��)
void swap_horizontal_pair(int& x, int& y, int xx, int yy, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  while (x != board_size - 1) {
    apply_move(x, y, 2, tmpBoard, ansDfs);
  }
  // �㍶���E
  vector<int> order = { 1, 0, 3, 2, 3, 0, 1, 1, 2, 3, 0, 3, 2, 1, 1, 0, 3 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

void fix_last_two_in_row(int& x, int& y, int ii, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  // �J�[�\���ړ�
  int dp[10][10];
  int dir[10][10];
  rep(i, board_size) {
    rep(j, board_size) {
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
    bool isFinish = false;
    rep(i, 4) {
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == ii && nb == board_size - 1) {
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
  int revX = ii, revY = board_size - 1;
  while (revX != x || revY != y) {
    int revD = dir[revX][revY];
    rev.push_back(revD);
    revX -= DX[revD];
    revY -= DY[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }

  // �㍶���E
  vector<int> order = { 1, 2 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

void fix_last_two_in_col(int& x, int& y, int jj, vector<vector<int>>& tmpBoard, vector<int>& ansDfs) {
  // �J�[�\���ړ�
  int dp[10][10];
  int dir[10][10];
  rep(i, board_size) {
    rep(j, board_size) {
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
    bool isFinish = false;
    rep(i, 4) {
      int na = a + DX[i];
      int nb = b + DY[i];
      if (in_bounds(na, nb) && dp[na][nb] > dp[a][b] + 1) {
        dp[na][nb] = dp[a][b] + 1;
        que.push(P(na, nb));
        dir[na][nb] = i;
        if (na == board_size - 1 && nb == jj) {
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
  int revX = board_size - 1, revY = jj;
  while (revX != x || revY != y) {
    int revD = dir[revX][revY];
    rev.push_back(revD);
    revX -= DX[revD];
    revY -= DY[revD];
  }

  reverse(rev.begin(), rev.end());

  for (auto nd : rev) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }

  // �㍶���E
  vector<int> order = { 0, 3 };
  for (auto nd : order) {
    ansDfs.push_back(nd);
    swap(tmpBoard[x][y], tmpBoard[x + DX[nd]][y + DY[nd]]);
    x += DX[nd];
    y += DY[nd];
  }
}

vector<int> build_route() {
  vector<vector<int>> tmpBoard(10, vector<int>(10));
  rep(i, board_size) {
    rep(j, board_size) {
      tmpBoard[i][j] = board[i][j];
    }
  }

  vector<int> ansDfs;

  int x = startX;
  int y = startY;

  // 1��ڂ��牺����3��ڂ܂Ŋ���������
  rep(i, board_size - 2) {
    rep(j, board_size - 1) {
      int num = aniBoard[i][j];
      if (j == board_size - 2) {
        num = aniBoard[i][board_size - 1];
      }

      vector<int> vxx, vyy;
      rep(k, board_size) {
        rep(l, board_size) {
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
        // �K�؂Ȉʒu�Ƀs�[�X���ړ�������

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

        if (tmpVec.size() < miniScore || rand32() % 100 == 0) {
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

    if (x == i && y == board_size - 1) {
      ansDfs.push_back(2);
      swap(tmpBoard[x][y], tmpBoard[x + 1][y]);
      x++;
    }

    // �Ō��2�𑵂���
    if (tmpBoard[i][board_size - 1] == aniBoard[i][board_size - 2]) {
      swap_vertical_pair(x, y, i, board_size - 2, tmpBoard, ansDfs);
    }
    else {
      int num = aniBoard[i][board_size - 2];
      vector<int> vxx, vyy;
      rep(k, board_size) {
        rep(l, board_size) {
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

        // �K�؂Ȉʒu�Ƀs�[�X���ړ�������
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

        if (tmpVec.size() < miniScore || rand32() % 100 == 0) {
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

  // ��2������낦��
  rep(j, board_size - 2) {
    // ���̃}�X����̃}�X�̈ʒu�Ɏ����Ă���
    int num = aniBoard[board_size - 1][j];
    int xx = -1, yy = -1;
    rep(k, board_size) {
      rep(l, board_size) {
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

    // �K�؂Ȉʒu�Ƀs�[�X���ړ�������
    while (xx != board_size - 2 || yy != j) {
      if (yy != j) {
        if (yy < j) {
          route_move_cursor(x, y, xx, yy, xx, yy + 1, board_size - 2, j, tmpBoard, ansDfs, 2);
        }
        else {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, board_size - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      else if (xx != board_size - 2) {
        route_move_cursor(x, y, xx, yy, xx - 1, yy, board_size - 2, j, tmpBoard, ansDfs);
      }
    }

    if (x == board_size - 1 && y == j) {
      ansDfs.push_back(3);
      swap(tmpBoard[x][y], tmpBoard[x][y + 1]);
      y++;
    }

    // �Ō��2�𑵂���
    if (tmpBoard[board_size - 1][j] == aniBoard[board_size - 2][j]) {
      swap_horizontal_pair(x, y, board_size - 2, j, tmpBoard, ansDfs);
    }
    else {
      int num = aniBoard[board_size - 2][j];
      int xx = -1, yy = -1;
      rep(k, board_size) {
        rep(l, board_size) {
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
      // �K�؂Ȉʒu�Ƀs�[�X���ړ�������
      while (xx != board_size - 2 || yy != j + 1) {
        if (xx != board_size - 2) {
          if (xx < board_size - 2) {
            route_move_cursor(x, y, xx, yy, xx + 1, yy, board_size - 2, j, tmpBoard, ansDfs, 2);
          }
          else {
            route_move_cursor(x, y, xx, yy, xx - 1, yy, board_size - 2, j, tmpBoard, ansDfs, 2);
          }
        }
        else if (yy != j + 1) {
          route_move_cursor(x, y, xx, yy, xx, yy - 1, board_size - 2, j, tmpBoard, ansDfs, 2);
        }
      }
      fix_last_two_in_col(x, y, j, tmpBoard, ansDfs);
    }
  }

  // �J�[�\�����E���Ɏ����Ă���
  while (x != board_size - 1 || y != board_size - 1) {
    if (y != board_size - 1) {
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

int exec_mode = 0; // 0: �W���o��, 1: �t�@�C���o��
void output_data(int case_num) {
  if (exec_mode == 0) {
    // �W���o��
    rep(i, route.size()) {
      cout << DIR_CHAR[route[i]];
    }
    cout << endl;
  }
  else {
    // �t�@�C���o��
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    rep(i, route.size()) {
      ofs << DIR_CHAR[route[i]];
    }
    ofs << endl;
  }
}

int solve_case(int mode, int problemNum = 0) {
  start_timer();

  // ���͕�
  read_input(problemNum);

  // �؂�1������
  bool isFind = false;
  {
    // �Ă��Ȃ܂�
    rep(_, 25) {
      bool isReset = true;
      if (anneal_find_tree(isReset)) {
        isFind = true;
        break;
      }
    }
    rep(i, board_size) {
      rep(j, board_size) {
        dfsBoard[i][j] = aniBoard[i][j];
      }
    }
  }

  if (isFind) {
    // ���ꂼ��̃s�[�X�ɔԍ���t����
    init_kind_indices();
    init_piece_numbers();

    vector<int> ansDfs = build_route();
    route = ansDfs;
    cur_score = calc_score(route);

    best_route = route;
    best_score = cur_score;

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

      vector<int> tmpAns = build_route();
      int tmpScore = calc_score(tmpAns);

      int diffScore = tmpScore - cur_score;

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
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
        // ���ɖ߂�
        rep(i, 2) {
          swap(piece_number[pp[i].first.first][pp[i].first.second], piece_number[pp[i].second.first][pp[i].second.second]);
        }
      }
      loop++;
    }
  }
  else {
    // �𒼉�
    {
      int x = startX;
      int y = startY;
      rep(i, turn_limit) {
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

    // �R�o����A�Ă��Ȃ܂���
    double now_time = get_elapsed_time();
    int loop = 0;
    while (true) {
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
        if (now_time > TL) break;
      }

      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      if (loop % 10 == 0) {
        shuffle_suffix(temp);
      }
      else {
        swap_adjacent(temp);
      }

      loop++;
    }

    // �ō��X�R�A��߂�
    route = best_route;
    cur_score = best_score;
  }

  // �f�o�b�O�p
  if (mode != 0) {
    cout << cur_score << endl;
    cout << get_elapsed_time() << "sec." << endl;
  }

  output_data(problemNum);

  return 0;
}

int main() {
  exec_mode = 10;

  if (exec_mode == 0) {
    solve_case(exec_mode);
  }
  else if (exec_mode == 1) {
    solve_case(exec_mode, 0);
  }
  else if (exec_mode == 10) {
    srep(_, 4, 5) {
      solve_case(exec_mode, _);
      reset_state();
    }
  }

  return 0;
}
