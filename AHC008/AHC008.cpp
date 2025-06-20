#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <queue>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

// いろいろ
const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char cc[4] = { 'U', 'L', 'D', 'R' };
const char c_pet[4] = { 'u', 'l', 'd', 'r' };

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
} // namespace

namespace /* 変数 */
{
  // 問題設定
  const int MAX_HUMAN = 10;
  const int MAX_PET = 20;
  const int N = 30;
  const int T = 300;

  // 入力用変数
  int n, m;
  vector<int> px, py, pt, hx, hy;
  int dog_cnt;

  // マスの状態
  /*
  0 : 通れる
  1 : 柵
  */
  int grid[N][N];

} // namespace

namespace /* ユーティリティ */
{
  bool ng_xy(int x, int y);
  bool can_fence(int x, int y);
  bool can_catch(int sx, int sy);
  P calc_route(int sx, int sy, int gx, int gy);
  int pet_catch(int id);

  // 人間と柵の干渉をチェックして座標を更新
  void check_collision(string& out)
  {
    for (int i = 0; i < m; ++i) {
      if (out[i] == 'U') {
        hx[i]--;
        if (grid[hx[i]][hy[i]] == 1) {
          hx[i]++;
          out[i] = '.';
        }
      }
      if (out[i] == 'D') {
        hx[i]++;
        if (grid[hx[i]][hy[i]] == 1) {
          hx[i]--;
          out[i] = '.';
        }
      }
      if (out[i] == 'L') {
        hy[i]--;
        if (grid[hx[i]][hy[i]] == 1) {
          hy[i]++;
          out[i] = '.';
        }
      }
      if (out[i] == 'R') {
        hy[i]++;
        if (grid[hx[i]][hy[i]] == 1) {
          hy[i]--;
          out[i] = '.';
        }
      }
    }
  }

  // ペットの移動を読み込んで座標を更新
  void read_pet_pos()
  {
    for (int i = 0; i < n; ++i) {
      string move;
      cin >> move;
      for (int j = 0; j < move.size(); ++j) {
        if (move[j] == 'U') {
          px[i]--;
        }
        if (move[j] == 'D') {
          px[i]++;
        }
        if (move[j] == 'L') {
          py[i]--;
        }
        if (move[j] == 'R') {
          py[i]++;
        }
      }
    }
  }

  // グリッドをデバッグ出力
  void print_grid()
  {
    for (int i = 0; i < N; ++i) {
      cout << "# ";
      for (int j = 0; j < N; ++j) {
        cout << grid[i][j];
      }
      cout << endl;
    }
  }

  // 左右の柵を作る処理
  void make_fence_both(int x, int y, string& out)
  {
    if (can_fence(x, y - 1)) {
      out += "l";
      grid[x][y - 1] = 1;
    }
    else if (can_fence(x, y + 1)) {
      out += "r";
      grid[x][y + 1] = 1;
    }
    else {
      out += ".";
    }
  }

  // 左の柵を作る処理
  void make_fence_left(int x, int y, string& out)
  {
    if (can_fence(x, y - 1)) {
      out += "l";
      grid[x][y - 1] = 1;
    }
    else {
      out += ".";
    }
  }

  // 右の柵を作る処理
  void make_fence_right(int x, int y, string& out)
  {
    if (can_fence(x, y + 1)) {
      out += "r";
      grid[x][y + 1] = 1;
    }
    else {
      out += ".";
    }
  }

  // ペットを取得できるかチェックして柵を作る
  int try_get_pet(int id, string& out)
  {
    int get_flag = -1;
    grid[hx[id]][hy[id]] = 1;
    for (int j = 0; j < 4; ++j) {
      int nx = hx[id] + dx[j];
      int ny = hy[id] + dy[j];
      if (ng_xy(nx, ny))
        continue;
      if (grid[nx][ny])
        continue;
      if (can_catch(nx, ny)) {
        get_flag = j;
        break;
      }
    }
    grid[hx[id]][hy[id]] = 0;

    if (get_flag != -1) {
      cout << "# " << "getpet" << endl;
      out += c_pet[get_flag];
      grid[hx[id] + dx[get_flag]][hy[id] + dy[get_flag]] = 1;
    }
    return get_flag;
  }

  // 中央の列に向かう動き
  void move_middle(int id, string& out)
  {
    if (hx[id] < N / 2) {
      out += "D";
    }
    else if (hx[id] > N / 2) {
      out += "U";
    }
    else {
      out += ".";
    }
  }

  // 左右に徘徊する動き
  void patrol(int id, int& dir, string& out)
  {
    if (dir == 0 && (hy[id] == 0 || grid[hx[id]][hy[id] - 1])) {
      dir = 1;
    }
    else if (dir == 1 && (hy[id] == N - 1 || grid[hx[id]][hy[id] + 1])) {
      dir = 0;
    }

    if (dir == 0) {
      out += "L";
    }
    else {
      out += "R";
    }
  }

  // 全員が目標位置に到達しているかチェック
  bool all_at_corners()
  {
    for (int i = 0; i < m; ++i) {
      if (i % 2 == 0) {
        if (hx[i] != 0 || hy[i] != 0) {
          return false;
        }
      }
      else {
        if (hx[i] != 0 || hy[i] != 29) {
          return false;
        }
      }
    }
    return true;
  }

  // 犬を捕獲可能な状態かチェック
  bool can_capture_dogs()
  {
    for (int i = 0; i < n; ++i) {
      if (pt[i] != 4)
        continue;
      if (px[i] != 0) {
        return false;
      }
      if (py[i] <= 1 || 28 <= py[i]) {
        return false;
      }
    }
    if (!can_fence(0, 1)) {
      return false;
    }
    if (!can_fence(0, 28)) {
      return false;
    }
    return true;
  }
  bool ng_xy(int x, int y)
  {
    return x < 0 || x >= N || y < 0 || y >= N;
  }

  // 最短距離計算
  // 戻り値 : (距離, 1歩目の方向)
  P calc_route(int sx, int sy, int gx, int gy)
  {
    int dp[N][N];
    int dir[N][N];
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        dp[i][j] = INF;
        dir[i][j] = -1;
      }
    }

    queue<P> que;
    que.push(P(sx, sy));
    dp[sx][sy] = 0;

    while (que.size()) {
      int x = que.front().first;
      int y = que.front().second;
      que.pop();
      if (x == gx && y == gy)
        break;
      for (int i = 0; i < 4; ++i) {
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (ng_xy(nx, ny))
          continue;
        if (grid[nx][ny] == 1)
          continue;
        if (dp[nx][ny] > dp[x][y] + 1) {
          dp[nx][ny] = dp[x][y] + 1;
          if (x == sx && y == sy) {
            dir[nx][ny] = i;
          }
          else {
            dir[nx][ny] = dir[x][y];
          }
          que.push(P(nx, ny));
        }
      }
    }

    return P(dp[gx][gy], dir[gx][gy]);
  }

  // 柵を作れるか
  bool can_fence(int x, int y)
  {
    if (ng_xy(x, y)) {
      return false;
    }
    if (grid[x][y] == 1) {
      return false;
    }
    // ペットとのマンハッタン距離が0or1ならNG
    for (int i = 0; i < n; ++i) {
      int dist = abs(x - px[i]) + abs(y - py[i]);
      if (dist <= 1) {
        return false;
      }
    }
    // 人間とのマンハッタン距離が0ならNG
    for (int i = 0; i < m; ++i) {
      if (hx[i] == x && hy[i] == y) {
        return false;
      }
    }
    return true;
  }

  // 人間idがこのターンペットを捕獲できるか
  // 1:上側で捕獲可能、2:下側で捕獲可能
  int pet_catch(int id)
  {
    if (can_fence(hx[id] - 1, hy[id])) {
      int ok = 1;
      for (int i = 0; i < hx[id]; ++i) {
        if (hy[id] != 0 && grid[i][hy[id] - 1] != 1) {
          ok = 0;
          break;
        }
        if (hy[id] != N - 1 && grid[i][hy[id] + 1] != 1) {
          ok = 0;
          break;
        }
        for (int j = 0; j < m; ++j) {
          if (hx[j] == i && hy[j] == hy[id]) {
            ok = 0;
            break;
          }
        }
        if (ok == 0)
          break;
      }

      if (ok) {
        for (int i = 0; i < hx[id]; ++i) {
          for (int j = 0; j < n; ++j) {
            if (px[j] == i && py[j] == hy[id]) {
              return 1;
            }
          }
        }
      }
    }

    if (can_fence(hx[id] + 1, hy[id])) {
      int ok = 1;
      for (int i = N - 1; i > hx[id]; i--) {
        if (hy[id] != 0 && grid[i][hy[id] - 1] != 1) {
          ok = 0;
          break;
        }
        if (hy[id] != N - 1 && grid[i][hy[id] + 1] != 1) {
          ok = 0;
          break;
        }
        for (int j = 0; j < m; ++j) {
          if (hx[j] == i && hy[j] == hy[id]) {
            ok = 0;
            break;
          }
        }
        if (ok == 0)
          break;
      }
      if (ok) {
        for (int i = N - 1; i > hx[id]; i--) {
          for (int j = 0; j < n; ++j) {
            if (px[j] == i && py[j] == hy[id]) {
              return 2;
            }
          }
        }
      }
    }

    return 0;
  }

  // (sx,sy)に柵を作ることで20マス以下の区画でペットを孤立させられるか
  bool can_catch(int sx, int sy)
  {
    if (!can_fence(sx, sy))
      return false;
    int personFlag = 0;
    int petFlag = 0;
    queue<P> que;
    que.push(P(sx, sy));
    vector<P> keep;
    keep.push_back(P(sx, sy));
    grid[sx][sy] = 1;
    while (que.size()) {
      int x = que.front().first;
      int y = que.front().second;
      que.pop();
      for (int i = 0; i < 4; ++i) {
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (ng_xy(nx, ny))
          continue;
        if (grid[nx][ny])
          continue;
        for (int j = 0; j < m; ++j) {
          if (hx[j] == nx && hy[j] == ny) {
            personFlag = 1;
          }
        }
        if (personFlag)
          break;
        for (int j = 0; j < n; ++j) {
          if (px[j] == nx && py[j] == ny) {
            petFlag = 1;
          }
        }
        que.push(P(nx, ny));
        keep.push_back(P(nx, ny));
        grid[nx][ny] = 1;
      }
      if (keep.size() > 20 || personFlag) {
        break;
      }
    }
    for (auto p : keep) {
      grid[p.first][p.second] = 0;
    }

    if (keep.size() > 20 || personFlag || !petFlag) {
      return false;
    }
    return true;
  }
}

// 櫛状に柵を作る解法
class SolveVer1
{
public:
  void Solve()
  {
    /*
    まず柵を作る
    柵が作り終わったら15列で待機
    */

    // 人間の状態
    /*
    0 : 何もしていない
    1 : 柵のスタート地点に向かっている
    2 : 柵を作っている
    */
    int humanMode[MAX_HUMAN] = {};
    // 目的の柵の列番号
    int myFence[MAX_HUMAN] = {};
    for (int i = 0; i < MAX_HUMAN; ++i) {
      myFence[i] = -1;
    }
    // 徘徊する方向
    // 0:LEFT, 1:RIGHT
    int dir[MAX_HUMAN] = {};

    // 柵の状況
    /*
    0 : 未着手
    1 : 着手中
    10 : 完成済み
    */
    int finishFenceCount = 0;
    int fence[N] = {};
    for (int i = 0; i < N; ++i) {
      if (i % 4 != 1) {
        fence[i] = 10;
        finishFenceCount++;
      }
    }

    for (int turn = 0; turn < T; ++turn) {
      string out;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        if (finishFenceCount < N) {
          if (humanMode[i] == 0) {
            for (int j = 0; j < N; ++j) {
              if (fence[j] == 0) {
                myFence[i] = j;
                humanMode[i] = 1;
                fence[j] = 1;
                break;
              }
            }
          }

          if (humanMode[i] == 1 && hx[i] == 0 && hy[i] == myFence[i]) {
            humanMode[i] = 2;
          }

          if (humanMode[i] == 1) {
            P p = calc_route(hx[i], hy[i], 0, myFence[i]);
            out += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hx[i] == N / 2) {
              out += "D";
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              make_fence_both(hx[i], hy[i], out);
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              make_fence_left(hx[i], hy[i], out);
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              make_fence_right(hx[i], hy[i], out);
            }
            else {
              // 両サイド柵あり
              if (hx[i] == N - 1) {
                humanMode[i] = 0;
                fence[myFence[i]] = 10;
                myFence[i] = -1;
                finishFenceCount++;
                out += ".";
              }
              else {
                out += "D";
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              out += "D";
            }
            else if (hx[i] > N / 2) {
              out += "U";
            }
            else {
              if (dir[i] == 0 && hy[i] == 0) {
                dir[i] = 1;
              }
              else if (dir[i] == 1 && hy[i] == N - 1) {
                dir[i] = 0;
              }

              if (pet_catch(i) == 1) {
                out += "u";
                grid[hx[i] - 1][hy[i]] = 1;
              }
              else if (pet_catch(i) == 2) {
                out += "d";
                grid[hx[i] + 1][hy[i]] = 1;
              }
              else if (dir[i] == 0) {
                out += "L";
              }
              else {
                out += "R";
              }
            }
          }
        }
        else {
          if (hx[i] != N / 2) {
            if (hx[i] < N / 2) {
              out += "D";
            }
            else {
              out += "U";
            }
          }
          else {
            if (dir[i] == 0 && hy[i] == 0) {
              dir[i] = 1;
            }
            else if (dir[i] == 1 && hy[i] == N - 1) {
              dir[i] = 0;
            }

            if (pet_catch(i) == 1) {
              out += "u";
              grid[hx[i] - 1][hy[i]] = 1;
            }
            else if (pet_catch(i) == 2) {
              out += "d";
              grid[hx[i] + 1][hy[i]] = 1;
            }
            else if (dir[i] == 0) {
              out += "L";
            }
            else {
              out += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      // このターンの人間と柵の干渉チェック
      // 柵作成を優先するためコメントアウト
      // for (int i = 0; i < m; ++i) {
      //   int x = -1, y = -1;
      //   if (out[i] == 'u') {
      //     x = hx[i] - 1;
      //     y = hy[i];
      //   }
      //   if (out[i] == 'd') {
      //     x = hx[i] + 1;
      //     y = hy[i];
      //   }
      //   if (out[i] == 'l') {
      //     x = hx[i];
      //     y = hy[i] - 1;
      //   }
      //   if (out[i] == 'r') {
      //     x = hx[i];
      //     y = hy[i] + 1;
      //   }
      //   if (x == -1) { continue; }
      //   for (int j = 0; j < m; ++j) {
      //     if (hx[j] == x && hy[j] == y) {
      //       out[i]  = '.';
      //       grid[x][y] = 0;
      //     }
      //   }
      // }

      cout << out << endl;
      fflush(stdout);

      read_pet_pos();
    }
  }
};

class SolveVer2
{
public:
  vector<vector<int>> sample;
  vector<P> route[11];
  int touch[11][1000][4];

  SolveVer2()
  {
    Init();
  }

  void Init()
  {
    sample = {
        {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3},
        {0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0},
        {10, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3},
        {0, 10, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0},
        {10, 0, 0, 10, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0},
        {0, 0, 10, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4},
        {0, 10, 0, 0, 9, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0},
        {10, 0, 0, 9, 0, 0, 9, 0, 0, 2, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0},
        {0, 0, 9, 0, 0, 9, 0, 0, 0, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4},
        {0, 9, 0, 0, 9, 0, 0, 9, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0},
        {9, 0, 0, 9, 0, 0, 9, 0, 0, 9, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0, 0},
        {0, 0, 9, 0, 0, 9, 0, 0, 9, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0, 0, 4},
        {0, 9, 0, 0, 9, 0, 0, 9, 0, 0, 8, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0},
        {9, 0, 0, 9, 0, 0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 3, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0},
        {0, 0, 9, 0, 0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 3, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 4},
        {0, 9, 0, 0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0},
        {9, 0, 0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0},
        {0, 0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 5},
        {0, 9, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 5, 0},
        {9, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 4, 0, 0, 4, 0, 0, 5, 0, 0},
        {0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 0, 0, 0, 4, 0, 0, 5, 0, 0, 5},
        {0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 0, 0, 0, 5, 0, 0, 5, 0},
        {8, 0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 5, 0, 0, 5, 0, 0},
        {0, 0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 0, 0, 0, 5, 0, 0, 5},
        {0, 8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 6, 0, 0, 0, 0, 0, 5, 0},
        {8, 0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 6, 0, 0, 6, 0, 0, 5, 0, 0},
        {0, 0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 6, 0, 0, 6, 0, 0, 0, 0, 0, 5},
        {0, 8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 6, 0, 0, 6, 0, 0, 6, 0, 0, 0, 0},
        {8, 0, 0, 8, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 7, 0, 0, 6, 0, 0, 6, 0, 0, 6, 0, 0, 6, 0, 0} };

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        for (int k = 0; k < 4; ++k) {
          touch[i][j][k] = 0;
        }
      }
    }

    // route[1]
    {
      int x = 3, y = 4;
    }
  }
};

// 犬を捕まえる
class SolveVer3
{
public:
  void Solve2(int startTurn)
  {
    /*
    まず柵を作る
    柵が作り終わったら15列で待機
    */

    // 人間の状態
    /*
    0 : 何もしていない
    1 : 柵のスタート地点に向かっている
    2 : 柵を作っている
    */
    int humanMode[MAX_HUMAN] = {};
    // 目的の柵の列番号
    int myFence[MAX_HUMAN] = {};
    for (int i = 0; i < MAX_HUMAN; ++i) {
      myFence[i] = -1;
    }
    // 徘徊する方向
    // 0:LEFT, 1:RIGHT
    int dir[MAX_HUMAN] = {};

    // 柵の状況
    /*
    0 : 未着手
    1 : 着手中
    10 : 完成済み
    */
    int finishFenceCount = 0;
    int fence[N] = {};
    for (int i = 0; i < N; ++i) {
      if (i % 4 != 1) {
        fence[i] = 10;
        finishFenceCount++;
      }
    }

    for (int turn = startTurn; turn < T; ++turn) {
      // cout << "# Solve2" << endl;
      string out;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int get_flag = try_get_pet(i, out);
        if (get_flag == -1 && finishFenceCount < N) {
          if (humanMode[i] == 0) {
            if (i % 2 == 0) {
              for (int j = 0; j < N; ++j) {
                if (fence[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fence[j] = 1;
                  break;
                }
              }
            }
            else {
              for (int j = N - 1; j >= 0; --j) {
                if (fence[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fence[j] = 1;
                  break;
                }
              }
            }
          }

          int startX = 2;
          if (myFence[i] == 1 || myFence[i] == 29) {
            startX = 3;
          }

          if (humanMode[i] == 1 && hx[i] == startX && hy[i] == myFence[i]) {
            humanMode[i] = 2;
          }

          cout << "#humanMode[" << i << "] = " << humanMode[i] << endl;
          if (humanMode[i] == 1) {
            P p = calc_route(hx[i], hy[i], startX, myFence[i]);
            out += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hx[i] == N / 2) {
              out += "D";
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              make_fence_both(hx[i], hy[i], out);
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              make_fence_left(hx[i], hy[i], out);
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              make_fence_right(hx[i], hy[i], out);
            }
            else {
              // 両サイド柵あり
              if (hx[i] == N - 1) {
                humanMode[i] = 0;
                fence[myFence[i]] = 10;
                myFence[i] = -1;
                finishFenceCount++;
                out += ".";
              }
              else {
                out += "D";
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              out += "D";
            }
            else if (hx[i] > N / 2) {
              out += "U";
            }
            else {
              patrol(i, dir[i], out);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            out += "D";
          }
          else if (hx[i] > N / 2) {
            out += "U";
          }
          else {
            if (dir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              dir[i] = 1;
            }
            else if (dir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              dir[i] = 0;
            }

            if (dir[i] == 0) {
              out += "L";
            }
            else {
              out += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      print_grid();

      cout << out << endl;
      fflush(stdout);

      read_pet_pos();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string out;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (can_capture_dogs()) {
          out += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          out += "..";
        }
        for (int i = 2; i < m; ++i) {
          out += ".";
        }
      }
      else {
        for (int i = 0; i < m; ++i) {
          if (finish == 0) {
            if (i == 0) {
              int gx = -1, gy = -1;
              for (int j = 14; j > 0; j--) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j - 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][1] == 0) {
                  gx = 2;
                  gy = 1 - 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy + 1)) {
                    out += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              for (int j = 15; j < 29; ++j) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j + 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][28] == 0) {
                  gx = 2;
                  gy = 28 + 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy - 1)) {
                    out += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = calc_route(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                out += ".";
              }
              else {
                out += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = calc_route(hx[i], hy[i], 0, 0);
            }
            else {
              p = calc_route(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              out += ".";
            }
            else {
              out += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      cout << out << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (all_at_corners()) {
          finish = 2;
        }
      }

      read_pet_pos();
    }

    Solve2(turn);
  }
};

// SolveVer3を改良
class SolveVer4
{
public:
  void Solve2(int startTurn)
  {
    /*
    まず柵を作る
    柵が作り終わったら15列で待機
    */

    // 人間の状態
    /*
    0 : 何もしていない
    1 : 柵のスタート地点に向かっている
    2 : 柵を作っている
    */
    int humanMode[MAX_HUMAN] = {};
    // 目的の柵の列番号
    int myFence[MAX_HUMAN] = {};
    for (int i = 0; i < MAX_HUMAN; ++i) {
      myFence[i] = -1;
    }
    // 徘徊する方向
    // 0:LEFT, 1:RIGHT
    int dir[MAX_HUMAN] = {};

    // 柵の状況
    /*
    0 : 未着手
    1 : 着手中
    10 : 完成済み
    */
    int finishFenceCount = 0;
    int fenceU[N] = {};
    int fenceD[N] = {};
    for (int i = 0; i < N; ++i) {
      if (i % 4 != 1) {
        fenceU[i] = 10;
        finishFenceCount++;
        fenceD[i] = 10;
        finishFenceCount++;
      }
    }

    for (int turn = startTurn; turn < T; ++turn) {
      // cout << "# Solve2" << endl;
      string out;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int get_flag = try_get_pet(i, out);
        if (get_flag == -1 && finishFenceCount < N * 2) {
          if (humanMode[i] == 0) {
            if (i % 2 == 0) {
              for (int j = 0; j < N; ++j) {
                if (fenceU[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fenceU[j] = 1;
                  break;
                }
                if (fenceD[j] == 0) {
                  myFence[i] = j + 100;
                  humanMode[i] = 1;
                  fenceD[j] = 1;
                  break;
                }
              }
            }
            else {
              for (int j = N - 1; j >= 0; --j) {
                if (fenceU[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fenceU[j] = 1;
                  break;
                }
                if (fenceD[j] == 0) {
                  myFence[i] = j + 100;
                  humanMode[i] = 1;
                  fenceD[j] = 1;
                  break;
                }
              }
            }
          }

          int startX = 14;

          if (myFence[i] >= 100) {
            startX = 16;
          }

          if (humanMode[i] == 1 && hx[i] == startX && hy[i] == myFence[i] % 100) {
            humanMode[i] = 2;
          }

          cout << "#humanMode[" << i << "] = " << humanMode[i] << endl;
          if (humanMode[i] == 1) {
            P p = calc_route(hx[i], hy[i], startX, myFence[i] % 100);
            out += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              if (can_fence(hx[i], hy[i] - 1)) {
                out += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else if (can_fence(hx[i], hy[i] + 1)) {
                out += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              if (can_fence(hx[i], hy[i] - 1)) {
                out += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              if (can_fence(hx[i], hy[i] + 1)) {
                out += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else {
              // 両サイド柵あり
              if (myFence[i] < 100) {
                int goalX = 2;
                if (myFence[i] == 1 || myFence[i] == 29) {
                  goalX = 3;
                }
                if (hx[i] == goalX) {
                  humanMode[i] = 0;
                  fenceU[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  out += "D";
                }
                else {
                  out += "U";
                }
              }
              else {
                if (hx[i] == N - 1) {
                  humanMode[i] = 0;
                  fenceD[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  out += "U";
                }
                else {
                  out += "D";
                }
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              out += "D";
            }
            else if (hx[i] > N / 2) {
              out += "U";
            }
            else {
              patrol(i, dir[i], out);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            out += "D";
          }
          else if (hx[i] > N / 2) {
            out += "U";
          }
          else {
            if (dir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              dir[i] = 1;
            }
            else if (dir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              dir[i] = 0;
            }

            if (dir[i] == 0) {
              out += "L";
            }
            else {
              out += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      print_grid();

      cout << out << endl;
      fflush(stdout);

      read_pet_pos();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string out;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (can_capture_dogs()) {
          out += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          out += "..";
        }
        for (int i = 2; i < m; ++i) {
          out += ".";
        }
      }
      else {
        for (int i = 0; i < m; ++i) {
          if (finish == 0) {
            if (i == 0) {
              int gx = -1, gy = -1;
              for (int j = 14; j > 0; j--) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j - 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][1] == 0) {
                  gx = 2;
                  gy = 1 - 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy + 1)) {
                    out += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              for (int j = 15; j < 29; ++j) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j + 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][28] == 0) {
                  gx = 2;
                  gy = 28 + 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy - 1)) {
                    out += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = calc_route(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                out += ".";
              }
              else {
                out += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = calc_route(hx[i], hy[i], 0, 0);
            }
            else {
              p = calc_route(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              out += ".";
            }
            else {
              out += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      cout << out << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (all_at_corners()) {
          finish = 2;
        }
      }

      read_pet_pos();
    }

    Solve2(turn);
  }
};

// SolveVer4と向きを逆に
class SolveVer5
{
public:
  void Solve2(int startTurn)
  {
    /*
    まず柵を作る
    柵が作り終わったら15列で待機
    */

    // 人間の状態
    /*
    0 : 何もしていない
    1 : 柵のスタート地点に向かっている
    2 : 柵を作っている
    */
    int humanMode[MAX_HUMAN] = {};
    // 目的の柵の列番号
    int myFence[MAX_HUMAN] = {};
    for (int i = 0; i < MAX_HUMAN; ++i) {
      myFence[i] = -1;
    }
    // 徘徊する方向
    // 0:LEFT, 1:RIGHT
    int dir[MAX_HUMAN] = {};

    // 柵の状況
    /*
    0 : 未着手
    1 : 着手中
    10 : 完成済み
    */
    int finishFenceCount = 0;
    int fenceU[N] = {};
    int fenceD[N] = {};
    for (int i = 0; i < N; ++i) {
      if (i % 4 != 1) {
        fenceU[i] = 10;
        finishFenceCount++;
        fenceD[i] = 10;
        finishFenceCount++;
      }
    }

    for (int turn = startTurn; turn < T; ++turn) {
      // cout << "# Solve2" << endl;
      string out;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int get_flag = try_get_pet(i, out);
        if (get_flag == -1 && finishFenceCount < N * 2) {
          if (humanMode[i] == 0) {
            if (i % 2 == 0) {
              for (int j = 0; j < N; ++j) {
                if (fenceU[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fenceU[j] = 1;
                  break;
                }
                if (fenceD[j] == 0) {
                  myFence[i] = j + 100;
                  humanMode[i] = 1;
                  fenceD[j] = 1;
                  break;
                }
              }
            }
            else {
              for (int j = N - 1; j >= 0; --j) {
                if (fenceU[j] == 0) {
                  myFence[i] = j;
                  humanMode[i] = 1;
                  fenceU[j] = 1;
                  break;
                }
                if (fenceD[j] == 0) {
                  myFence[i] = j + 100;
                  humanMode[i] = 1;
                  fenceD[j] = 1;
                  break;
                }
              }
            }
          }

          int startX = 14;

          if (myFence[i] >= 100) {
            startX = 16;
          }

          if (humanMode[i] == 1 && hx[i] == startX && hy[i] == myFence[i] % 100) {
            humanMode[i] = 2;
          }

          cout << "#humanMode[" << i << "] = " << humanMode[i] << endl;
          if (humanMode[i] == 1) {
            P p = calc_route(hx[i], hy[i], startX, myFence[i] % 100);
            out += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              if (can_fence(hx[i], hy[i] - 1)) {
                out += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else if (can_fence(hx[i], hy[i] + 1)) {
                out += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              if (can_fence(hx[i], hy[i] - 1)) {
                out += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              if (can_fence(hx[i], hy[i] + 1)) {
                out += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                out += ".";
              }
            }
            else {
              // 両サイド柵あり
              if (myFence[i] < 100) {
                int goalX = 2;
                if (myFence[i] == 1 || myFence[i] == 29) {
                  goalX = 3;
                }
                if (hx[i] == goalX) {
                  humanMode[i] = 0;
                  fenceU[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  out += "D";
                }
                else {
                  out += "U";
                }
              }
              else {
                if (hx[i] == N - 1) {
                  humanMode[i] = 0;
                  fenceD[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  out += "U";
                }
                else {
                  out += "D";
                }
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              out += "D";
            }
            else if (hx[i] > N / 2) {
              out += "U";
            }
            else {
              patrol(i, dir[i], out);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            out += "D";
          }
          else if (hx[i] > N / 2) {
            out += "U";
          }
          else {
            if (dir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              dir[i] = 1;
            }
            else if (dir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              dir[i] = 0;
            }

            if (dir[i] == 0) {
              out += "L";
            }
            else {
              out += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      print_grid();

      cout << out << endl;
      fflush(stdout);

      read_pet_pos();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string out;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (can_capture_dogs()) {
          out += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          out += "..";
        }
        for (int i = 2; i < m; ++i) {
          out += ".";
        }
      }
      else {
        for (int i = 0; i < m; ++i) {
          if (finish == 0) {
            if (i == 0) {
              int gx = -1, gy = -1;
              for (int j = 14; j > 0; j--) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j - 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][1] == 0) {
                  gx = 2;
                  gy = 1 - 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy + 1)) {
                    out += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              for (int j = 15; j < 29; ++j) {
                if (grid[1][j] == 0) {
                  gx = 1;
                  gy = j + 1;
                  break;
                }
              }
              if (gx == -1) {
                if (grid[2][28] == 0) {
                  gx = 2;
                  gy = 28 + 1;
                }
              }
              if (gx == -1) {
                out += ".";
              }
              else {
                P p = calc_route(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (can_fence(gx, gy - 1)) {
                    out += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    out += ".";
                  }
                }
                else {
                  out += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = calc_route(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                out += ".";
              }
              else {
                out += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = calc_route(hx[i], hy[i], 0, 0);
            }
            else {
              p = calc_route(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              out += ".";
            }
            else {
              out += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      check_collision(out);

      cout << out << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (all_at_corners()) {
          finish = 2;
        }
      }

      read_pet_pos();
    }

    Solve2(turn);
  }
};

// 入力受け取り
void Input()
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << 0 << ".txt";
  ifstream ifs(oss.str());
  if (!ifs.is_open()) { // 標準入力する
    cin >> n;
    for (int i = 0; i < n; ++i) {
      int x, y, t;
      cin >> x >> y >> t;
      x--;
      y--;
      px.push_back(x);
      py.push_back(y);
      pt.push_back(t);
    }
    cin >> m;
    for (int i = 0; i < m; ++i) {
      int x, y;
      cin >> x >> y;
      x--;
      y--;
      hx.push_back(x);
      hy.push_back(y);
    }
  }
  else { // ファイル入力する
    ifs >> n;
    for (int i = 0; i < n; ++i) {
      int x, y, t;
      ifs >> x >> y >> t;
      px.push_back(x);
      py.push_back(y);
      pt.push_back(t);
    }
    ifs >> m;
    for (int i = 0; i < m; ++i) {
      int x, y;
      ifs >> x >> y;
      hx.push_back(x);
      hy.push_back(y);
    }
  }

  for (int i = 0; i < n; ++i) {
    if (pt[i] == 4) {
      dog_cnt++;
    }
  }
}

int Solve(int mode)
{
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();

  // 入力受け取り
  Input();

  if (dog_cnt == 0) {
    SolveVer1 solve1;
    solve1.Solve();
  }
  else {
    SolveVer4 solve4;
    solve4.Solve();
  }

  return 0;
}

int main()
{
  int mode = 0;

  if (mode == 0) {
    Solve(mode);
  }
  else if (mode == 1) {
  }

  return 0;
}
