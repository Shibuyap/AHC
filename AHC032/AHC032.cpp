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
}  // namespace

double TL = 1.9;
int mode;
clock_t startTime, endTime;

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int n = 9;
const int m = 20;
const int T = 81;
const int MOD = 998244353;
int board[n][n];
int initial_board[n][n];
int magic_pattern[m][3][3];
int current_solution[T][3];
ll current_score;

int best_board[n][n];
ll best_solution[T][3];
ll best_score;

void CopyToBest()
{
  best_score = current_score;
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < 3; ++j) {
      best_solution[i][j] = current_solution[i][j];
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      best_board[i][j] = board[i][j];
    }
  }
}

void CopyFromBest()
{
  current_score = best_score;
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < 3; ++j) {
      current_solution[i][j] = best_solution[i][j];
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      board[i][j] = best_board[i][j];
    }
  }
}

// スコア計算
ll CalcScore()
{
  ll res = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      res += board[i][j] % MOD;
    }
  }

  return res;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < 3; ++j) {
      current_solution[i][j] = -1;
    }
  }
}

void InitializeAns()
{
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < 3; ++j) {
      current_solution[i][j] = -1;
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      board[i][j] = initial_board[i][j];
    }
  }
  current_score = CalcScore();
}

// 入力受け取り
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    int _n, _m, _k;
    cin >> _n >> _m >> _k;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        cin >> board[i][j];
      }
    }
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          cin >> magic_pattern[i][j][k];
        }
      }
    }
  }
  // ファイル入力する
  else {
    int _n, _m, _k;
    ifs >> _n >> _m >> _k;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        ifs >> board[i][j];
      }
    }
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          ifs >> magic_pattern[i][j][k];
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      initial_board[i][j] = board[i][j];
    }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// 初期解生成
void Initialize()
{
  for (int i = 0; i < T; ++i) {
    current_solution[i][0] = -1;
  }
  current_score = CalcScore();
  CopyToBest();
}

// 解答出力
void Output(ofstream& ofs)
{
  int L = 0;
  for (int i = 0; i < T; ++i) {
    if (current_solution[i][0] == -1) continue;
    L++;
  }
  if (mode == 0) {
    cout << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current_solution[i][0] == -1) continue;
      for (int j = 0; j < 3; ++j) cout << current_solution[i][j] << ' ';
      cout << endl;
    }
  }
  else {
    ofs << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current_solution[i][0] == -1) continue;
      for (int j = 0; j < 3; ++j) ofs << current_solution[i][j] << ' ';
      ofs << endl;
    }
  }
}

double nowTime = 0;
double startTemp = 1001001001;
double endTemp = 0.1;

int use[3][3];
void Rule1(int x, int y, int dir1, int dir2)
{
  for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) use[k][l] = 0;
  if (dir1 == 0 && dir2 == 0) {
    if (x == n - 3) {
      if (y == n - 3) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            use[i][j] = 1;
          }
        }
      }
      else {
        for (int i = 0; i < 3; ++i) use[i][0] = 1;
      }
    }
    else {
      if (y == n - 3) {
        for (int l = 0; l < 3; ++l) use[0][l] = 1;
      }
      else {
        use[0][0] = 1;
      }
    }
  }
  if (dir1 == 0 && dir2 == 1) {
    if (x == n - 3) {
      if (y == 0) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            use[i][j] = 1;
          }
        }
      }
      else {
        for (int i = 0; i < 3; ++i) use[i][2] = 1;
      }
    }
    else {
      if (y == 0) {
        for (int l = 0; l < 3; ++l) use[0][l] = 1;
      }
      else {
        use[0][2] = 1;
      }
    }
  }
  if (dir1 == 1 && dir2 == 0) {
    if (x == 0) {
      if (y == n - 3) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            use[i][j] = 1;
          }
        }
      }
      else {
        for (int i = 0; i < 3; ++i) use[i][0] = 1;
      }
    }
    else {
      if (y == n - 3) {
        for (int l = 0; l < 3; ++l) use[2][l] = 1;
      }
      else {
        use[2][0] = 1;
      }
    }
  }
  if (dir1 == 1 && dir2 == 1) {
    if (x == 0) {
      if (y == 0) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            use[i][j] = 1;
          }
        }
      }
      else {
        for (int i = 0; i < 3; ++i) use[i][2] = 1;
      }
    }
    else {
      if (y == 0) {
        for (int l = 0; l < 3; ++l) use[2][l] = 1;
      }
      else {
        use[2][2] = 1;
      }
    }
  }
}

void Rule2(int i, int j, int dir1)
{
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {
      use[p][q] = 0;
      if (dir1 == 0) {
        if (i == n - 3) use[p][q] = 1;
        if (p == 0) use[p][q] = 1;
      }
      else {
        if (i == 0) use[p][q] = 1;
        if (p == 2) use[p][q] = 1;
      }
    }
  }
}

ll maxSum;
ll ma[3][3];
ll maAnssArr[10];
int maAnsCount = 0;

ll now[3][3];
ll anssArr[10];
void Method2DFS(int mm, int cnt, int lim)
{
  if (cnt == lim) return;
  ll keep[3][3];
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {
      keep[p][q] = now[p][q];
    }
  }
  for (int i = mm; i < m; ++i) {
    anssArr[cnt] = i;
    cnt++;
    ll tmpSum = 0;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        now[p][q] = (now[p][q] + magic_pattern[i][p][q]) % MOD;
        if (use[p][q]) tmpSum += now[p][q];
      }
    }
    if (tmpSum > maxSum) {
      maxSum = tmpSum;
      for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) {
          ma[p][q] = now[p][q];
        }
      }
      for (int j = 0; j < cnt; ++j) maAnssArr[j] = anssArr[j];
      maAnsCount = cnt;
    }

    Method2DFS(i, cnt, lim);

    cnt--;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        now[p][q] = keep[p][q];
      }
    }
  }
}

int keepA[n][n];
int keepAns[110][3];
int baseA[n][n];

void Method4(double timeLimit)
{
  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) break;

    ll hosyou = Rand() % 200000000 + MOD - 200000000;

    int ng = 0;
    InitializeAns();
    int cnt = 0;
    int dir1 = Rand() % 2;
    int dir2 = Rand() % 2;
    // dir1     = 0;
    dir2 = 0;
    for (int ii = 0; ii < n - 2; ++ii) {
      int i = ii;
      if (dir1) i = n - 3 - ii;
      if (ii == n - 3 && cnt + 3 * 6 + 4 > T) {
        ng = 1;
        break;
      }

      int nowCnt = cnt;
      ll maPosSum = 0;
      int maCntTail = 0;
      for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
          baseA[p][q] = board[p][q];
        }
      }

      for (int jjj = 0; jjj < n - 2; ++jjj) {
        int ng2 = 0;
        for (int jj = 0; jj < jjj; ++jj) {
          int j = jj;
          dir2 = 0;
          if (dir2) j = n - 3 - jj;
          maAnsCount = 0;
          maxSum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              ma[k][l] = board[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          Rule1(i, j, dir1, dir2);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use[k][l]) {
                maxSum += ma[k][l];
              }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  now[p][q] = board[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            Method2DFS(0, 0, 4);
          }

          for (int k = 0; k < maAnsCount; ++k) {
            int ansM = maAnssArr[k];
            current_solution[cnt][0] = ansM;
            current_solution[cnt][1] = i;
            current_solution[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              board[i + k][j + l] = ma[k][l];
            }
          }
        }
        if (ng) break;

        for (int jj = n - 3; jj > jjj; jj--) {
          int j = jj;
          dir2 = 1;
          maAnsCount = 0;
          maxSum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              ma[k][l] = board[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          Rule1(i, j, dir1, dir2);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use[k][l]) {
                maxSum += ma[k][l];
              }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  now[p][q] = board[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            Method2DFS(0, 0, 4);
          }

          for (int k = 0; k < maAnsCount; ++k) {
            int ansM = maAnssArr[k];
            current_solution[cnt][0] = ansM;
            current_solution[cnt][1] = i;
            current_solution[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              board[i + k][j + l] = ma[k][l];
            }
          }
        }

        {
          int j = jjj;
          maAnsCount = 0;
          maxSum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              ma[k][l] = board[i + k][j + l];
              now[k][l] = ma[k][l];
            }
          }
          // Rule1(i, j, dir1, dir2);
          Rule2(i, j, dir1);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use[k][l]) {
                maxSum += ma[k][l];
              }
            }
          }

          if (useCount == 1) {
            Method2DFS(0, 0, 1);
            if (maxSum < hosyou && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  now[p][q] = board[i + p][j + q];
                }
              }
              Method2DFS(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            Method2DFS(0, 0, 3);
          }
          else {
            if (cnt + 4 > T) {
              ng = 1;
              break;
            }
            int num = min(6, T - cnt);
            Method2DFS(0, 0, num);
          }

          for (int k = 0; k < maAnsCount; ++k) {
            int ansM = maAnssArr[k];
            current_solution[cnt][0] = ansM;
            current_solution[cnt][1] = i;
            current_solution[cnt][2] = j;
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              board[i + k][j + l] = ma[k][l];
            }
          }
        }

        if (ng) break;

        ll tmpPosSum = 0;

        if (dir1 == 0) {
          for (int j = 0; j < n; ++j) {
            tmpPosSum += board[i][j];
          }
          if (i == n - 3) {
            for (int j = 0; j < n; ++j) {
              tmpPosSum += board[i + 1][j];
              tmpPosSum += board[i + 2][j];
            }
          }
        }
        else {
          for (int j = 0; j < n; ++j) {
            tmpPosSum += board[i + 2][j];
          }
          if (i == 0) {
            for (int j = 0; j < n; ++j) {
              tmpPosSum += board[i][j];
              tmpPosSum += board[i + 1][j];
            }
          }
        }
        if (tmpPosSum > maPosSum) {
          maPosSum = tmpPosSum;
          for (int t = nowCnt; t < cnt; ++t) {
            for (int k = 0; k < 3; ++k) keepAns[t][k] = current_solution[t][k];
          }
          for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
              keepA[p][q] = board[p][q];
            }
          }
          maCntTail = cnt;
        }

        cnt = nowCnt;
        for (int p = 0; p < n; ++p) {
          for (int q = 0; q < n; ++q) {
            board[p][q] = baseA[p][q];
          }
        }
      }

      if (ng) break;

      for (int t = nowCnt; t < maCntTail; ++t) {
        for (int k = 0; k < 3; ++k) {
          current_solution[t][k] = keepAns[t][k];
        }
      }
      cnt = maCntTail;
      for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
          board[p][q] = keepA[p][q];
        }
      }
    }
    if (ng) continue;

    for (int t = cnt; t < T; ++t) current_solution[t][0] = -1;
    current_score = CalcScore();
    if (current_score > best_score) {
      CopyToBest();
    }
  }

  if (mode != 0) cout << loopCount << endl;
}

ll Solve(int probNum)
{
  startTime = clock();
  endTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();
  Method4(TL);
  CopyFromBest();

  CopyFromBest();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
