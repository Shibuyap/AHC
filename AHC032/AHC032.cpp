#include <algorithm>
#include <array>
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

int initial_board[n][n];
int magic_pattern[m][3][3];

class State
{
public:
  array<array<int, n>, n> board;
  array<array<int, 3>, T> solution; // solution[i][0]: magic index, solution[i][1]: x, solution[i][2]: y
  ll score;
  State()
  {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        board[i][j] = 0;
      }
    }
    for (int i = 0; i < T; ++i) {
      solution[i][0] = -1;
      solution[i][1] = -1;
      solution[i][2] = -1;
    }
    score = 0;
  }
};

// スコア計算
ll CalcScore(const State& current)
{
  ll res = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      res += current.board[i][j] % MOD;
    }
  }

  return res;
}

void InitializeAns(State& current)
{
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < 3; ++j) {
      current.solution[i][j] = -1;
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      current.board[i][j] = initial_board[i][j];
    }
  }
  current.score = CalcScore(current);
}

// 入力受け取り
State Input(int problemNum)
{
  State current;

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
        cin >> current.board[i][j];
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
        ifs >> current.board[i][j];
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
      initial_board[i][j] = current.board[i][j];
    }
  }

  return current;
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
void Initialize(State& current)
{
  for (int i = 0; i < T; ++i) {
    current.solution[i][0] = -1;
  }
  current.score = CalcScore(current);
}

// 解答出力
void Output(int probNum, const State& current)
{
  ofstream ofs;
  OpenOfs(probNum, ofs);

  int L = 0;
  for (int i = 0; i < T; ++i) {
    if (current.solution[i][0] == -1) continue;
    L++;
  }

  if (mode == 0) {
    cout << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current.solution[i][0] == -1) continue;
      for (int j = 0; j < 3; ++j) cout << current.solution[i][j] << ' ';
      cout << endl;
    }
  }
  else {
    ofs << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current.solution[i][0] == -1) continue;
      for (int j = 0; j < 3; ++j) ofs << current.solution[i][j] << ' ';
      ofs << endl;
    }
  }

  if (ofs.is_open()) {
    ofs.close();
  }
}

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

void Method4(double timeLimit, State& current)
{
  State best = current;

  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) break;

    ll hosyou = Rand() % 200000000 + MOD - 200000000;

    int ng = 0;
    InitializeAns(current);
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
          baseA[p][q] = current.board[p][q];
        }
      }

      for (int jjj = 0; jjj < n - 2; ++jjj) {
        for (int jj = 0; jj < jjj; ++jj) {
          int j = jj;
          dir2 = 0;
          if (dir2) j = n - 3 - jj;
          maAnsCount = 0;
          maxSum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              ma[k][l] = current.board[i + k][j + l];
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
                  now[p][q] = current.board[i + p][j + q];
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
            if (cnt < T) {
              current.solution[cnt][0] = ansM;
              current.solution[cnt][1] = i;
              current.solution[cnt][2] = j;
            }
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = ma[k][l];
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
              ma[k][l] = current.board[i + k][j + l];
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
                  now[p][q] = current.board[i + p][j + q];
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
            if (cnt < T) {
              current.solution[cnt][0] = ansM;
              current.solution[cnt][1] = i;
              current.solution[cnt][2] = j;
            }
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = ma[k][l];
            }
          }
        }

        {
          int j = jjj;
          maAnsCount = 0;
          maxSum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              ma[k][l] = current.board[i + k][j + l];
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
                  now[p][q] = current.board[i + p][j + q];
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
            if (cnt < T) {
              current.solution[cnt][0] = ansM;
              current.solution[cnt][1] = i;
              current.solution[cnt][2] = j;
            }
            cnt++;
          }
          if (cnt > T) {
            ng = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = ma[k][l];
            }
          }
        }

        if (ng) break;

        ll tmpPosSum = 0;

        if (dir1 == 0) {
          for (int j = 0; j < n; ++j) {
            tmpPosSum += current.board[i][j];
          }
          if (i == n - 3) {
            for (int j = 0; j < n; ++j) {
              tmpPosSum += current.board[i + 1][j];
              tmpPosSum += current.board[i + 2][j];
            }
          }
        }
        else {
          for (int j = 0; j < n; ++j) {
            tmpPosSum += current.board[i + 2][j];
          }
          if (i == 0) {
            for (int j = 0; j < n; ++j) {
              tmpPosSum += current.board[i][j];
              tmpPosSum += current.board[i + 1][j];
            }
          }
        }
        if (tmpPosSum > maPosSum) {
          maPosSum = tmpPosSum;
          for (int t = nowCnt; t < cnt; ++t) {
            for (int k = 0; k < 3; ++k) keepAns[t][k] = current.solution[t][k];
          }
          for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
              keepA[p][q] = current.board[p][q];
            }
          }
          maCntTail = cnt;
        }

        cnt = nowCnt;
        for (int p = 0; p < n; ++p) {
          for (int q = 0; q < n; ++q) {
            current.board[p][q] = baseA[p][q];
          }
        }
      }

      if (ng) break;

      for (int t = nowCnt; t < maCntTail; ++t) {
        for (int k = 0; k < 3; ++k) {
          current.solution[t][k] = keepAns[t][k];
        }
      }
      cnt = maCntTail;
      for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
          current.board[p][q] = keepA[p][q];
        }
      }
    }
    if (ng) continue;

    for (int t = cnt; t < T; ++t) current.solution[t][0] = -1;
    current.score = CalcScore(current);
    if (current.score > best.score) {
      best = current;
    }
  }

  current = best;

  if (mode != 0) cout << "Method 4 : " << loopCount << endl;
}

ll SolveCase(int probNum)
{
  startTime = clock();
  endTime = clock();

  State current = Input(probNum);

  Initialize(current);
  Method4(TL, current);

  Output(probNum, current);

  cerr << "Problem " << probNum << " solved." << endl;
  ll score = 0;
  if (mode != 0) {
    score = CalcScore(current);
    cerr << "Score for problem " << probNum << ": " << score << endl;
  }

  cerr << "Time taken for problem " << probNum << ": " << GetNowTime() << " seconds." << endl;
  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    SolveCase(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = SolveCase(i);
      cerr << "Score for problem " << i << ": " << score << endl;
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
