#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <utility>

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
  array<array<int, 3>, T> turns; // turns[i][0]: magic index, turns[i][1]: x, turns[i][2]: y
  ll score;
  State()
  {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        board[i][j] = 0;
      }
    }
    for (int i = 0; i < T; ++i) {
      turns[i][0] = -1;
      turns[i][1] = -1;
      turns[i][2] = -1;
    }
    score = 0;
  }

  void reset()
  {
    for (int i = 0; i < T; ++i) {
      for (int j = 0; j < 3; ++j) {
        turns[i][j] = -1;
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        board[i][j] = initial_board[i][j];
      }
    }
    score = calc_score();
  }

  // スコア計算
  ll calc_score()
  {
    score = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        score += board[i][j] % MOD;
      }
    }
    return score;
  }
};

// 入力受け取り
State Input(int problemNum)
{
  State current;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

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

// 解答出力
void Output(int probNum, const State& current)
{
  ofstream ofs;
  OpenOfs(probNum, ofs);

  int L = 0;
  for (int i = 0; i < T; ++i) {
    if (current.turns[i][0] == -1) { continue; }
    L++;
  }

  if (mode == 0) {
    cout << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current.turns[i][0] == -1) { continue; }
      for (int j = 0; j < 3; ++j) cout << current.turns[i][j] << ' ';
      cout << endl;
    }
  }
  else {
    ofs << L << endl;
    for (int i = 0; i < T; ++i) {
      if (current.turns[i][0] == -1) { continue; }
      for (int j = 0; j < 3; ++j) ofs << current.turns[i][j] << ' ';
      ofs << endl;
    }
  }

  if (ofs.is_open()) {
    ofs.close();
  }
}

// ── 3×3 ブロック内で「スコア対象になるマス」をマークする
// dirV = 0 … 上→下に走査、1 … 下→上
// dirH = 0 … 左→右に走査、1 … 右→左
// 盤端に達した方向に応じて列全体／行全体／全面を有効化
void mark_edge_mask_2(int x /*0-based*/, int y /*0-based*/,
  int dirV, int dirH,
  int(&mask)[3][3])           // 出力: 0 or 1
{
  // 1) 全消し
  for (auto& row : mask) std::fill(std::begin(row), std::end(row), 0);

  // 2) アンカー（通常 1 マスだけ）
  const int anchorRow = (dirV == 0 ? 0 : 2);   // 上端 or 下端
  const int anchorCol = (dirH == 0 ? 0 : 2);   // 左端 or 右端
  mask[anchorRow][anchorCol] = 1;

  // 3) 盤の端に当たったら行／列ごと拡張
  const bool hitVEdge = (dirV == 0 ? x == n - 3 : x == 0);
  const bool hitHEdge = (dirH == 0 ? y == n - 3 : y == 0);

  // -- 縦方向の端
  if (hitVEdge) {
    for (int r = 0; r < 3; ++r) {
      mask[r][anchorCol] = 1;
    }
  }

  // -- 横方向の端
  if (hitHEdge) {
    for (int c = 0; c < 3; ++c) {
      mask[anchorRow][c] = 1;
    }
  }

  // 4) 角 (両端) の場合は 3×3 全面
  if (hitVEdge && hitHEdge)
    for (auto& row : mask) std::fill(std::begin(row), std::end(row), 1);
}

int use_mask[3][3];

void mark_line_mask(int i, int j, int dir1)
{
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {
      use_mask[p][q] = 0;
      if (dir1 == 0) {
        if (i == n - 3) use_mask[p][q] = 1;
        if (p == 0) use_mask[p][q] = 1;
      }
      else {
        if (i == 0) use_mask[p][q] = 1;
        if (p == 2) use_mask[p][q] = 1;
      }
    }
  }
}

ll best_sum;
ll best_sum_grid[3][3];
ll best_pattern_ids[10];
int best_pattern_count = 0;

ll work_grid[3][3];
ll pattern_ids[10];
void search_best_patterns_dfs(int startIdx, int depth, int maxDepth)
{
  if (depth == maxDepth) { return; }
  ll keep[3][3];
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {
      keep[p][q] = work_grid[p][q];
    }
  }
  for (int i = startIdx; i < m; ++i) {
    pattern_ids[depth] = i;
    depth++;
    ll tmp_sum = 0;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        work_grid[p][q] = (work_grid[p][q] + magic_pattern[i][p][q]) % MOD;
        if (use_mask[p][q]) tmp_sum += work_grid[p][q];
      }
    }
    if (tmp_sum > best_sum) {
      best_sum = tmp_sum;
      for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) {
          best_sum_grid[p][q] = work_grid[p][q];
        }
      }
      for (int j = 0; j < depth; ++j) best_pattern_ids[j] = pattern_ids[j];
      best_pattern_count = depth;
    }

    search_best_patterns_dfs(i, depth, maxDepth);

    depth--;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        work_grid[p][q] = keep[p][q];
      }
    }
  }
}

int backupBoard[n][n];
int keepAns[110][3];
int originalBoard[n][n];

void heuristicSweep(double timeLimit, State& current)
{
  State best = current;

  int loopCount = 0;
  while (true) {
    loopCount++;
    double nowTime = GetNowTime();
    if (nowTime > timeLimit) {
      break;
    }

    ll scoreThreshold = Rand() % 200000000 + MOD - 200000000;

    int abortFlag = 0;
    current.reset();
    int turnCount = 0;
    int dir1 = Rand() % 2;
    int dir2 = Rand() % 2;
    // dir1     = 0;
    dir2 = 0;
    for (int rowIdx = 0; rowIdx < n - 2; ++rowIdx) {
      int i = rowIdx;
      if (dir1) i = n - 3 - rowIdx;
      if (rowIdx == n - 3 && turnCount + 3 * 6 + 4 > T) {
        abortFlag = 1;
        break;
      }

      int nowCnt = turnCount;
      ll maxPosSum = 0;
      int maxCntTail = 0;
      for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
          originalBoard[p][q] = current.board[p][q];
        }
      }

      for (int colBoundary = 0; colBoundary < n - 2; ++colBoundary) {
        for (int colIdx = 0; colIdx < colBoundary; ++colIdx) {
          int j = colIdx;
          dir2 = 0;
          if (dir2) j = n - 3 - colIdx;
          best_pattern_count = 0;
          best_sum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              best_sum_grid[k][l] = current.board[i + k][j + l];
              work_grid[k][l] = best_sum_grid[k][l];
            }
          }
          mark_edge_mask_2(i, j, dir1, dir2, use_mask);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use_mask[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use_mask[k][l]) {
                best_sum += best_sum_grid[k][l];
              }
            }
          }

          if (useCount == 1) {
            search_best_patterns_dfs(0, 0, 1);
            if (best_sum < scoreThreshold && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  work_grid[p][q] = current.board[i + p][j + q];
                }
              }
              search_best_patterns_dfs(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            search_best_patterns_dfs(0, 0, 3);
          }
          else {
            if (turnCount + 4 > T) {
              abortFlag = 1;
              break;
            }
            search_best_patterns_dfs(0, 0, 4);
          }

          for (int k = 0; k < best_pattern_count; ++k) {
            int ansM = best_pattern_ids[k];
            if (turnCount < T) {
              current.turns[turnCount][0] = ansM;
              current.turns[turnCount][1] = i;
              current.turns[turnCount][2] = j;
            }
            turnCount++;
          }
          if (turnCount > T) {
            abortFlag = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = best_sum_grid[k][l];
            }
          }
        }
        if (abortFlag) { break; }

        for (int colIdx = n - 3; colIdx > colBoundary; colIdx--) {
          int j = colIdx;
          dir2 = 1;
          best_pattern_count = 0;
          best_sum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              best_sum_grid[k][l] = current.board[i + k][j + l];
              work_grid[k][l] = best_sum_grid[k][l];
            }
          }
          mark_edge_mask_2(i, j, dir1, dir2, use_mask);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use_mask[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use_mask[k][l]) {
                best_sum += best_sum_grid[k][l];
              }
            }
          }

          if (useCount == 1) {
            search_best_patterns_dfs(0, 0, 1);
            if (best_sum < scoreThreshold && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  work_grid[p][q] = current.board[i + p][j + q];
                }
              }
              search_best_patterns_dfs(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            search_best_patterns_dfs(0, 0, 3);
          }
          else {
            if (turnCount + 4 > T) {
              abortFlag = 1;
              break;
            }
            search_best_patterns_dfs(0, 0, 4);
          }

          for (int k = 0; k < best_pattern_count; ++k) {
            int ansM = best_pattern_ids[k];
            if (turnCount < T) {
              current.turns[turnCount][0] = ansM;
              current.turns[turnCount][1] = i;
              current.turns[turnCount][2] = j;
            }
            turnCount++;
          }
          if (turnCount > T) {
            abortFlag = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = best_sum_grid[k][l];
            }
          }
        }

        {
          int j = colBoundary;
          best_pattern_count = 0;
          best_sum = 0;
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              best_sum_grid[k][l] = current.board[i + k][j + l];
              work_grid[k][l] = best_sum_grid[k][l];
            }
          }

          mark_line_mask(i, j, dir1);
          int useCount = 0;
          for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) useCount += use_mask[k][l];
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              if (use_mask[k][l]) {
                best_sum += best_sum_grid[k][l];
              }
            }
          }

          if (useCount == 1) {
            search_best_patterns_dfs(0, 0, 1);
            if (best_sum < scoreThreshold && Rand() % 3 != 0) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  work_grid[p][q] = current.board[i + p][j + q];
                }
              }
              search_best_patterns_dfs(0, 0, 2);
            }
          }
          else if (useCount <= 3) {
            search_best_patterns_dfs(0, 0, 3);
          }
          else {
            if (turnCount + 4 > T) {
              abortFlag = 1;
              break;
            }
            int num = min(6, T - turnCount);
            search_best_patterns_dfs(0, 0, num);
          }

          for (int k = 0; k < best_pattern_count; ++k) {
            int ansM = best_pattern_ids[k];
            if (turnCount < T) {
              current.turns[turnCount][0] = ansM;
              current.turns[turnCount][1] = i;
              current.turns[turnCount][2] = j;
            }
            turnCount++;
          }
          if (turnCount > T) {
            abortFlag = 1;
            break;
          }
          for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
              current.board[i + k][j + l] = best_sum_grid[k][l];
            }
          }
        }

        if (abortFlag) { break; }

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
        if (tmpPosSum > maxPosSum) {
          maxPosSum = tmpPosSum;
          for (int t = nowCnt; t < turnCount; ++t) {
            for (int k = 0; k < 3; ++k) keepAns[t][k] = current.turns[t][k];
          }
          for (int p = 0; p < n; ++p) {
            for (int q = 0; q < n; ++q) {
              backupBoard[p][q] = current.board[p][q];
            }
          }
          maxCntTail = turnCount;
        }

        turnCount = nowCnt;
        for (int p = 0; p < n; ++p) {
          for (int q = 0; q < n; ++q) {
            current.board[p][q] = originalBoard[p][q];
          }
        }
      }

      if (abortFlag) { break; }

      for (int t = nowCnt; t < maxCntTail; ++t) {
        for (int k = 0; k < 3; ++k) {
          current.turns[t][k] = keepAns[t][k];
        }
      }
      turnCount = maxCntTail;
      for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
          current.board[p][q] = backupBoard[p][q];
        }
      }
    }
    if (abortFlag) { continue; }

    for (int t = turnCount; t < T; ++t) current.turns[t][0] = -1;
    current.calc_score();
    if (current.score > best.score) {
      best = current;
    }
  }

  current = best;

  if (mode != 0) cout << "Method 4 : " << loopCount << endl;
}

ll solveOneProblem(int probNum)
{
  startTime = clock();
  endTime = clock();

  State current = Input(probNum);

  current.reset();
  heuristicSweep(TL, current);

  Output(probNum, current);

  ll score = 0;
  if (mode != 0) {
    score = current.calc_score();
  }

  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    solveOneProblem(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = solveOneProblem(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
