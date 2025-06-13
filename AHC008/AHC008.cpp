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

#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

// いろいろ
const int INF = 1001001001;
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char cc[4] = { 'U', 'L', 'D', 'R' };
const char cGetPet[4] = { 'u', 'l', 'd', 'r' };

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
  int dogCount;

  // マスの状態
  /*
  0 : 通れる
  1 : 柵
  */
  int grid[N][N];

} // namespace

namespace /* ユーティリティ */
{
  bool NgXY(int x, int y);
  bool CanMakeFence(int x, int y);
  bool CanCatchPet(int sx, int sy);
  P CalcRouteBfs(int sx, int sy, int gx, int gy);
  int PetCatch(int id);

  // 人間と柵の干渉をチェックして座標を更新
  void CheckHumanFenceCollision(string& outStr)
  {
    for (int i = 0; i < m; ++i) {
      if (outStr[i] == 'U') {
        hx[i]--;
        if (grid[hx[i]][hy[i]] == 1) {
          hx[i]++;
          outStr[i] = '.';
        }
      }
      if (outStr[i] == 'D') {
        hx[i]++;
        if (grid[hx[i]][hy[i]] == 1) {
          hx[i]--;
          outStr[i] = '.';
        }
      }
      if (outStr[i] == 'L') {
        hy[i]--;
        if (grid[hx[i]][hy[i]] == 1) {
          hy[i]++;
          outStr[i] = '.';
        }
      }
      if (outStr[i] == 'R') {
        hy[i]++;
        if (grid[hx[i]][hy[i]] == 1) {
          hy[i]--;
          outStr[i] = '.';
        }
      }
    }
  }

  // ペットの移動を読み込んで座標を更新
  void ReadAndUpdatePetPositions()
  {
    for (int i = 0; i < n; ++i) {
      string petMove;
      cin >> petMove;
      for (int j = 0; j < petMove.size(); ++j) {
        if (petMove[j] == 'U') {
          px[i]--;
        }
        if (petMove[j] == 'D') {
          px[i]++;
        }
        if (petMove[j] == 'L') {
          py[i]--;
        }
        if (petMove[j] == 'R') {
          py[i]++;
        }
      }
    }
  }

  // グリッドをデバッグ出力
  void PrintGrid()
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
  void TryMakeFenceBothSides(int x, int y, string& outStr)
  {
    if (CanMakeFence(x, y - 1)) {
      outStr += "l";
      grid[x][y - 1] = 1;
    }
    else if (CanMakeFence(x, y + 1)) {
      outStr += "r";
      grid[x][y + 1] = 1;
    }
    else {
      outStr += ".";
    }
  }

  // 左の柵を作る処理
  void TryMakeFenceLeft(int x, int y, string& outStr)
  {
    if (CanMakeFence(x, y - 1)) {
      outStr += "l";
      grid[x][y - 1] = 1;
    }
    else {
      outStr += ".";
    }
  }

  // 右の柵を作る処理
  void TryMakeFenceRight(int x, int y, string& outStr)
  {
    if (CanMakeFence(x, y + 1)) {
      outStr += "r";
      grid[x][y + 1] = 1;
    }
    else {
      outStr += ".";
    }
  }

  // ペットを取得できるかチェックして柵を作る
  int CheckAndTryGetPet(int id, string& outStr)
  {
    int getPetFlag = -1;
    grid[hx[id]][hy[id]] = 1;
    for (int j = 0; j < 4; ++j) {
      int nx = hx[id] + dx[j];
      int ny = hy[id] + dy[j];
      if (NgXY(nx, ny))
        continue;
      if (grid[nx][ny])
        continue;
      if (CanCatchPet(nx, ny)) {
        getPetFlag = j;
        break;
      }
    }
    grid[hx[id]][hy[id]] = 0;

    if (getPetFlag != -1) {
      cout << "# " << "getpet" << endl;
      outStr += cGetPet[getPetFlag];
      grid[hx[id] + dx[getPetFlag]][hy[id] + dy[getPetFlag]] = 1;
    }
    return getPetFlag;
  }

  // 中央の列に向かう動き
  void MoveToMiddleRow(int id, string& outStr)
  {
    if (hx[id] < N / 2) {
      outStr += "D";
    }
    else if (hx[id] > N / 2) {
      outStr += "U";
    }
    else {
      outStr += ".";
    }
  }

  // 左右に徘徊する動き
  void PatrolLeftRight(int id, int& humanDir, string& outStr)
  {
    if (humanDir == 0 && (hy[id] == 0 || grid[hx[id]][hy[id] - 1])) {
      humanDir = 1;
    }
    else if (humanDir == 1 && (hy[id] == N - 1 || grid[hx[id]][hy[id] + 1])) {
      humanDir = 0;
    }

    if (humanDir == 0) {
      outStr += "L";
    }
    else {
      outStr += "R";
    }
  }

  // 全員が目標位置に到達しているかチェック
  bool CheckAllHumansAtCorners()
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
  bool CanCaptureDogs()
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
    if (!CanMakeFence(0, 1)) {
      return false;
    }
    if (!CanMakeFence(0, 28)) {
      return false;
    }
    return true;
  }
  bool NgXY(int x, int y)
  {
    return x < 0 || x >= N || y < 0 || y >= N;
  }

  // 最短距離計算
  // 戻り値 : (距離, 1歩目の方向)
  P CalcRouteBfs(int sx, int sy, int gx, int gy)
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
        if (NgXY(nx, ny))
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
  bool CanMakeFence(int x, int y)
  {
    if (NgXY(x, y)) {
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
  int PetCatch(int id)
  {
    if (CanMakeFence(hx[id] - 1, hy[id])) {
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

    if (CanMakeFence(hx[id] + 1, hy[id])) {
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
  bool CanCatchPet(int sx, int sy)
  {
    if (!CanMakeFence(sx, sy))
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
        if (NgXY(nx, ny))
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
    int humanDir[MAX_HUMAN] = {};

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
      string outStr;
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
            P p = CalcRouteBfs(hx[i], hy[i], 0, myFence[i]);
            outStr += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hx[i] == N / 2) {
              outStr += "D";
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              TryMakeFenceBothSides(hx[i], hy[i], outStr);
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              TryMakeFenceLeft(hx[i], hy[i], outStr);
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              TryMakeFenceRight(hx[i], hy[i], outStr);
            }
            else {
              // 両サイド柵あり
              if (hx[i] == N - 1) {
                humanMode[i] = 0;
                fence[myFence[i]] = 10;
                myFence[i] = -1;
                finishFenceCount++;
                outStr += ".";
              }
              else {
                outStr += "D";
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              outStr += "D";
            }
            else if (hx[i] > N / 2) {
              outStr += "U";
            }
            else {
              if (humanDir[i] == 0 && hy[i] == 0) {
                humanDir[i] = 1;
              }
              else if (humanDir[i] == 1 && hy[i] == N - 1) {
                humanDir[i] = 0;
              }

              if (PetCatch(i) == 1) {
                outStr += "u";
                grid[hx[i] - 1][hy[i]] = 1;
              }
              else if (PetCatch(i) == 2) {
                outStr += "d";
                grid[hx[i] + 1][hy[i]] = 1;
              }
              else if (humanDir[i] == 0) {
                outStr += "L";
              }
              else {
                outStr += "R";
              }
            }
          }
        }
        else {
          if (hx[i] != N / 2) {
            if (hx[i] < N / 2) {
              outStr += "D";
            }
            else {
              outStr += "U";
            }
          }
          else {
            if (humanDir[i] == 0 && hy[i] == 0) {
              humanDir[i] = 1;
            }
            else if (humanDir[i] == 1 && hy[i] == N - 1) {
              humanDir[i] = 0;
            }

            if (PetCatch(i) == 1) {
              outStr += "u";
              grid[hx[i] - 1][hy[i]] = 1;
            }
            else if (PetCatch(i) == 2) {
              outStr += "d";
              grid[hx[i] + 1][hy[i]] = 1;
            }
            else if (humanDir[i] == 0) {
              outStr += "L";
            }
            else {
              outStr += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      // このターンの人間と柵の干渉チェック
      // 柵作成を優先するためコメントアウト
      // for (int i = 0; i < m; ++i) {
      //   int x = -1, y = -1;
      //   if (outStr[i] == 'u') {
      //     x = hx[i] - 1;
      //     y = hy[i];
      //   }
      //   if (outStr[i] == 'd') {
      //     x = hx[i] + 1;
      //     y = hy[i];
      //   }
      //   if (outStr[i] == 'l') {
      //     x = hx[i];
      //     y = hy[i] - 1;
      //   }
      //   if (outStr[i] == 'r') {
      //     x = hx[i];
      //     y = hy[i] + 1;
      //   }
      //   if (x == -1) continue;
      //   for (int j = 0; j < m; ++j) {
      //     if (hx[j] == x && hy[j] == y) {
      //       outStr[i]  = '.';
      //       grid[x][y] = 0;
      //     }
      //   }
      // }

      cout << outStr << endl;
      fflush(stdout);

      ReadAndUpdatePetPositions();
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
    int humanDir[MAX_HUMAN] = {};

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

    srep(turn, startTurn, T)
    {
      // cout << "# Solve2" << endl;
      string outStr;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int getPetFlag = CheckAndTryGetPet(i, outStr);
        if (getPetFlag == -1 && finishFenceCount < N) {
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
              drep(j, N)
              {
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
            P p = CalcRouteBfs(hx[i], hy[i], startX, myFence[i]);
            outStr += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hx[i] == N / 2) {
              outStr += "D";
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              TryMakeFenceBothSides(hx[i], hy[i], outStr);
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              TryMakeFenceLeft(hx[i], hy[i], outStr);
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              TryMakeFenceRight(hx[i], hy[i], outStr);
            }
            else {
              // 両サイド柵あり
              if (hx[i] == N - 1) {
                humanMode[i] = 0;
                fence[myFence[i]] = 10;
                myFence[i] = -1;
                finishFenceCount++;
                outStr += ".";
              }
              else {
                outStr += "D";
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              outStr += "D";
            }
            else if (hx[i] > N / 2) {
              outStr += "U";
            }
            else {
              PatrolLeftRight(i, humanDir[i], outStr);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            outStr += "D";
          }
          else if (hx[i] > N / 2) {
            outStr += "U";
          }
          else {
            if (humanDir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              humanDir[i] = 1;
            }
            else if (humanDir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              humanDir[i] = 0;
            }

            if (humanDir[i] == 0) {
              outStr += "L";
            }
            else {
              outStr += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      PrintGrid();

      cout << outStr << endl;
      fflush(stdout);

      ReadAndUpdatePetPositions();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string outStr;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (CanCaptureDogs()) {
          outStr += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          outStr += "..";
        }
        srep(i, 2, m)
        {
          outStr += ".";
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy + 1)) {
                    outStr += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              srep(j, 15, 29)
              {
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy - 1)) {
                    outStr += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = CalcRouteBfs(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                outStr += ".";
              }
              else {
                outStr += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = CalcRouteBfs(hx[i], hy[i], 0, 0);
            }
            else {
              p = CalcRouteBfs(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              outStr += ".";
            }
            else {
              outStr += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      cout << outStr << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (CheckAllHumansAtCorners()) {
          finish = 2;
        }
      }

      ReadAndUpdatePetPositions();
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
    int humanDir[MAX_HUMAN] = {};

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

    srep(turn, startTurn, T)
    {
      // cout << "# Solve2" << endl;
      string outStr;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int getPetFlag = CheckAndTryGetPet(i, outStr);
        if (getPetFlag == -1 && finishFenceCount < N * 2) {
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
              drep(j, N)
              {
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
            P p = CalcRouteBfs(hx[i], hy[i], startX, myFence[i] % 100);
            outStr += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              if (CanMakeFence(hx[i], hy[i] - 1)) {
                outStr += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else if (CanMakeFence(hx[i], hy[i] + 1)) {
                outStr += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                outStr += ".";
              }
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              if (CanMakeFence(hx[i], hy[i] - 1)) {
                outStr += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else {
                outStr += ".";
              }
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              if (CanMakeFence(hx[i], hy[i] + 1)) {
                outStr += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                outStr += ".";
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
                  outStr += "D";
                }
                else {
                  outStr += "U";
                }
              }
              else {
                if (hx[i] == N - 1) {
                  humanMode[i] = 0;
                  fenceD[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  outStr += "U";
                }
                else {
                  outStr += "D";
                }
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              outStr += "D";
            }
            else if (hx[i] > N / 2) {
              outStr += "U";
            }
            else {
              PatrolLeftRight(i, humanDir[i], outStr);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            outStr += "D";
          }
          else if (hx[i] > N / 2) {
            outStr += "U";
          }
          else {
            if (humanDir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              humanDir[i] = 1;
            }
            else if (humanDir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              humanDir[i] = 0;
            }

            if (humanDir[i] == 0) {
              outStr += "L";
            }
            else {
              outStr += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      PrintGrid();

      cout << outStr << endl;
      fflush(stdout);

      ReadAndUpdatePetPositions();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string outStr;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (CanCaptureDogs()) {
          outStr += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          outStr += "..";
        }
        srep(i, 2, m)
        {
          outStr += ".";
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy + 1)) {
                    outStr += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              srep(j, 15, 29)
              {
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy - 1)) {
                    outStr += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = CalcRouteBfs(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                outStr += ".";
              }
              else {
                outStr += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = CalcRouteBfs(hx[i], hy[i], 0, 0);
            }
            else {
              p = CalcRouteBfs(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              outStr += ".";
            }
            else {
              outStr += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      cout << outStr << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (CheckAllHumansAtCorners()) {
          finish = 2;
        }
      }

      ReadAndUpdatePetPositions();
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
    int humanDir[MAX_HUMAN] = {};

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

    srep(turn, startTurn, T)
    {
      // cout << "# Solve2" << endl;
      string outStr;
      // 人間の行動
      for (int i = 0; i < m; ++i) {
        int getPetFlag = CheckAndTryGetPet(i, outStr);
        if (getPetFlag == -1 && finishFenceCount < N * 2) {
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
              drep(j, N)
              {
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
            P p = CalcRouteBfs(hx[i], hy[i], startX, myFence[i] % 100);
            outStr += cc[p.second];
          }
          else if (humanMode[i] == 2) {
            if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0 && hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 左にも右にも柵がない
              if (CanMakeFence(hx[i], hy[i] - 1)) {
                outStr += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else if (CanMakeFence(hx[i], hy[i] + 1)) {
                outStr += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                outStr += ".";
              }
            }
            else if (hy[i] != 0 && grid[hx[i]][hy[i] - 1] == 0) {
              // 左に柵がない
              if (CanMakeFence(hx[i], hy[i] - 1)) {
                outStr += "l";
                grid[hx[i]][hy[i] - 1] = 1;
              }
              else {
                outStr += ".";
              }
            }
            else if (hy[i] != N - 1 && grid[hx[i]][hy[i] + 1] == 0) {
              // 右に柵がない
              if (CanMakeFence(hx[i], hy[i] + 1)) {
                outStr += "r";
                grid[hx[i]][hy[i] + 1] = 1;
              }
              else {
                outStr += ".";
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
                  outStr += "D";
                }
                else {
                  outStr += "U";
                }
              }
              else {
                if (hx[i] == N - 1) {
                  humanMode[i] = 0;
                  fenceD[myFence[i] % 100] = 10;
                  myFence[i] = -1;
                  finishFenceCount++;
                  outStr += "U";
                }
                else {
                  outStr += "D";
                }
              }
            }
          }
          else {
            if (hx[i] < N / 2) {
              outStr += "D";
            }
            else if (hx[i] > N / 2) {
              outStr += "U";
            }
            else {
              PatrolLeftRight(i, humanDir[i], outStr);
            }
          }
        }
        else {
          if (hx[i] < N / 2) {
            outStr += "D";
          }
          else if (hx[i] > N / 2) {
            outStr += "U";
          }
          else {
            if (humanDir[i] == 0 && (hy[i] == 0 || grid[hx[i]][hy[i] - 1])) {
              humanDir[i] = 1;
            }
            else if (humanDir[i] == 1 && (hy[i] == N - 1 || grid[hx[i]][hy[i] + 1])) {
              humanDir[i] = 0;
            }

            if (humanDir[i] == 0) {
              outStr += "L";
            }
            else {
              outStr += "R";
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      PrintGrid();

      cout << outStr << endl;
      fflush(stdout);

      ReadAndUpdatePetPositions();
    }
  }

  void Solve()
  {
    int finish = 0;
    int blockCount = 0;
    int turn = 0;
    for (turn = 0; turn < T; turn++) {
      string outStr;
      if (finish == 3) {
        break;
      }
      else if (finish == 2) {
        if (CanCaptureDogs()) {
          outStr += "rl";
          grid[0][1] = 1;
          grid[0][28] = 1;
          finish = 3;
        }
        else {
          outStr += "..";
        }
        srep(i, 2, m)
        {
          outStr += ".";
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy + 1)) {
                    outStr += "r";
                    grid[gx][gy + 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else if (i == 1) {
              int gx = -1, gy = -1;
              srep(j, 15, 29)
              {
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
                outStr += ".";
              }
              else {
                P p = CalcRouteBfs(hx[i], hy[i], gx, gy);
                if (p.first == 0) {
                  if (CanMakeFence(gx, gy - 1)) {
                    outStr += "l";
                    grid[gx][gy - 1] = 1;
                    blockCount++;
                  }
                  else {
                    outStr += ".";
                  }
                }
                else {
                  outStr += cc[p.second];
                }
              }
            }
            else {
              int targetY = (i % 2 == 0) ? 0 : 29;
              P p = CalcRouteBfs(hx[i], hy[i], 0, targetY);
              if (p.first == 0) {
                outStr += ".";
              }
              else {
                outStr += cc[p.second];
              }
            }
          }
          if (finish == 1) {
            P p;
            if (i % 2 == 0) {
              p = CalcRouteBfs(hx[i], hy[i], 0, 0);
            }
            else {
              p = CalcRouteBfs(hx[i], hy[i], 0, 29);
            }
            if (p.first == 0) {
              outStr += ".";
            }
            else {
              outStr += cc[p.second];
            }
          }
        }
      }

      // このターンの人間と柵の干渉チェック
      CheckHumanFenceCollision(outStr);

      cout << outStr << endl;
      cout << "# finish = " << finish << endl;
      fflush(stdout);

      if (finish == 0 && blockCount == 30) {
        finish = 1;
      }
      else if (finish == 1) {
        if (CheckAllHumansAtCorners()) {
          finish = 2;
        }
      }

      ReadAndUpdatePetPositions();
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
      dogCount++;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/*
    メモ
    cat ./in/0001.txt | ./tester.exe ./a.exe > 0001_out.txt
*/
///////////////////////////////////////////////////////////////////////////////
int Solve(int mode)
{
  srand((unsigned)time(NULL));
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();
  while (rand() % 100) {
    Rand();
  }

  // 入力受け取り
  Input();

  if (dogCount == 0) {
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
