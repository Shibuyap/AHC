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
Timer timer;

const int MAX_D = 55;
const int MAX_N = 56;
const int MAX_LINECOUNT = 60;

// スコア計算関連の定数
const int SCORE_MULTIPLIER = 100;
const int GRID_WIDTH = 1000;
const int INVALID_SIZE = 1001001;

// 最適化パラメータ
const int MAX_CANDIDATES = 30;
const int RANDOM_TRIAL_COUNT = 32;
const int SMALL_DATA_THRESHOLD = 5;
const int DEFAULT_TRIAL_LIMIT = 10;

// 焼きなまし法のパラメータ
const double ANNEALING_START_TEMP = 10000.1;
const double ANNEALING_END_TEMP = 0.1;
const int ANNEALING_TEMP_COEFFICIENT = 100;

class Solution
{
public:
  int ans[MAX_D][MAX_N][4];
  ll ansScore;
  int ansLinePos[MAX_D][MAX_LINECOUNT];
  int ansLineCount[MAX_D];
  int ansBaseLineCount;

  Solution()
  {
    clear();
  }

  void clear()
  {
    memset(ans, 0, sizeof(ans));
    ansScore = 0;
    memset(ansLinePos, 0, sizeof(ansLinePos));
    memset(ansLineCount, 0, sizeof(ansLineCount));
    ansBaseLineCount = 0;
  }

  void copyFrom(const Solution& other)
  {
    memcpy(ans, other.ans, sizeof(ans));
    ansScore = other.ansScore;
    memcpy(ansLinePos, other.ansLinePos, sizeof(ansLinePos));
    memcpy(ansLineCount, other.ansLineCount, sizeof(ansLineCount));
    ansBaseLineCount = other.ansBaseLineCount;
  }

  void copyTo(Solution& other) const
  {
    other.copyFrom(*this);
  }

  // 座標を取得するヘルパー関数
  struct Rectangle
  {
    int startX, startY, endX, endY;
  };

  Rectangle getRect(int day, int element) const
  {
    return {
      ans[day][element][0],
      ans[day][element][1],
      ans[day][element][2],
      ans[day][element][3]
    };
  }

  void setRect(int day, int element, int sx, int sy, int ex, int ey)
  {
    ans[day][element][0] = sx;
    ans[day][element][1] = sy;
    ans[day][element][2] = ex;
    ans[day][element][3] = ey;
  }
};

// 列スケジュール管理クラス
class ColumnSchedule
{
public:
  int columnNum[MAX_D][MAX_N];  // 各要素がどの列に属するか
  int schedules[MAX_D][MAX_LINECOUNT][MAX_N];  // 各列に属する要素のリスト
  int schedulesCount[MAX_D][MAX_LINECOUNT];  // 各列の要素数
  int schedulesPosition[MAX_D][MAX_LINECOUNT][MAX_N];  // 各要素の位置

  ColumnSchedule()
  {
    clear();
  }

  void clear()
  {
    memset(columnNum, 0, sizeof(columnNum));
    memset(schedules, 0, sizeof(schedules));
    memset(schedulesCount, 0, sizeof(schedulesCount));
    memset(schedulesPosition, 0, sizeof(schedulesPosition));
  }

  void copyFrom(const ColumnSchedule& other)
  {
    memcpy(columnNum, other.columnNum, sizeof(columnNum));
    memcpy(schedules, other.schedules, sizeof(schedules));
    memcpy(schedulesCount, other.schedulesCount, sizeof(schedulesCount));
    memcpy(schedulesPosition, other.schedulesPosition, sizeof(schedulesPosition));
  }

  // 要素を特定の列に追加
  void addElementToColumn(int day, int column, int element, int position)
  {
    columnNum[day][element] = column;
    int idx = schedulesCount[day][column];
    schedules[day][column][idx] = element;
    schedulesPosition[day][column][idx] = position;
    schedulesCount[day][column]++;
  }

  // 要素のインデックスを検索
  int findElementIndex(int day, int column, int element) const
  {
    for (int k = 0; k < schedulesCount[day][column]; k++) {
      if (schedules[day][column][k] == element) {
        return k;
      }
    }
    return -1;
  }
};

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

  // 配列シャッフル
  void FisherYates(int* data, int elementCount)
  {
    for (int i = elementCount - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      int temp = data[i];
      data[i] = data[j];
      data[j] = temp;
    }
  }

  // std::array用のオーバーロード
  template<size_t N>
  void FisherYates(std::array<int, N>& data, int elementCount)
  {
    for (int i = elementCount - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      std::swap(data[i], data[j]);
    }
  }
}  // namespace

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

// 2つのソート済み配列をマージしながらスコアを計算する関数
ll calculateMergedScore(const vector<P>& before, const vector<P>& now)
{
  ll score = 0;
  int itr1 = 0, itr2 = 0;
  int cnt1 = 0, cnt2 = 0;
  int pos = 0;

  while (itr1 < before.size() || itr2 < now.size()) {
    int val;
    if (itr1 == before.size()) {
      val = now[itr2].first;
      cnt2 += now[itr2].second;
      itr2++;
    }
    else if (itr2 == now.size()) {
      val = before[itr1].first;
      cnt1 += before[itr1].second;
      itr1++;
    }
    else if (before[itr1].first <= now[itr2].first) {
      val = before[itr1].first;
      cnt1 += before[itr1].second;
      itr1++;
    }
    else {
      val = now[itr2].first;
      cnt2 += now[itr2].second;
      itr2++;
    }

    if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) {
      score += (ll)val - pos;
    }
    pos = val;
  }

  return score;
}

// 要素が指定された幅に収まるために必要な高さを計算
inline int calculateRequiredSize(int elementSize, int width)
{
  return (elementSize - 1) / width + 1;
}

// 焼きなまし法の受理確率を計算
inline bool acceptByAnnealing(double diffScore, double startTemp, double endTemp, double currentTime, double timeLimit)
{
  double temp = startTemp + (endTemp - startTemp) * currentTime / timeLimit;
  double prob = exp(diffScore * ANNEALING_TEMP_COEFFICIENT / temp);
  return prob > Rand01();
}

// 隣接日の差分スコアを計算（CalcDiffScore2の後に定義）
int calculateAdjacentDaysDiffScore(int day, int lineNum);

// シャッフル試行回数の制限をチェック
inline bool shouldBreakShuffle(int size, int iteration)
{
  if (size == 1 && iteration >= 1) return true;
  if (size == 2 && iteration >= 2) return true;
  if (size == 3 && iteration >= 6) return true;
  return false;
}

// スコアを更新して最適解を保存（グローバル変数の後に定義）
void updateScoreAndSaveBest(ll diffScore);

// 列内の要素位置を再計算（グローバル変数の後に定義）
void recalculateColumnPositions(int day, int column);

// 最適な配置位置を見つける構造体
struct PlacementResult
{
  int position;
  int need;
  int amari;
  int over;
  bool found;
};

// 要素に対して最適な配置位置を見つける（グローバル変数の後に定義）
PlacementResult findBestPlacement(int day, int element, const std::array<int, MAX_LINECOUNT>& currentPos, int startLine, int endLine);

// 要素が指定された幅に収まるかチェック（グローバル変数の後に定義）
bool canFitInColumn(int day, int element, int width);

// 要素を配置するのに必要な余白を計算（グローバル変数の後に定義）
int calculateMargin(int day, int element, int width, int height);

// 横線のカウントを計算（グローバル変数の後に定義）
int calculateHorizontalLineCount(int day, int column);

// 要素の交換による差分スコアを計算（グローバル変数の後に定義）
int calculateSwapDiffScore(int day, int elem1, int elem2);

const double TIME_LIMIT = 2.8;
double TL = TIME_LIMIT;
int mode;

// Constants are already defined in Solution.hpp
int lineMaxLimit = 100;

const int w = GRID_WIDTH;
int dayCount, elementCount;  // dayCount -> dayCount, elementCount -> elementCount
double emptySpaceRatio;  // emptySpaceRatio -> emptySpaceRatio
std::array<std::array<int, MAX_N>, MAX_D> elementSizes;  // a -> elementSizes
std::array<int, MAX_N> maxElementSize;  // maxA -> maxElementSize
std::array<int, MAX_D> dailyTotalSize;  // sumA -> dailyTotalSize
std::array<int, MAX_D> daysDifficultySorted;

int bestPairCount[MAX_D][w + 10];
int bestPairValue[MAX_D][w + 10][10];
int bestPairs[MAX_D][w + 10][10][2];

Solution ans;

Solution temp_ans;

Solution best_ans;

Solution backup_ans;

int keep31Count = 0;
int keep31KeepSize = 2;
std::array<Solution, 10> keep31_ans;

void CopyToBackupAns()
{
  backup_ans.copyFrom(best_ans);
}

void CopyToBestAns()
{
  best_ans.copyFrom(ans);
}

void CopyFromBackupAns()
{
  best_ans.copyFrom(backup_ans);
}

void CopyFromBestAns()
{
  ans.copyFrom(best_ans);
}

void CopyToKeep31(int idx)
{
  keep31_ans[idx].copyFrom(ans);
}

void CopyFromKeep31(int idx)
{
  ans.copyFrom(keep31_ans[idx]);
}

void CopyToTemp()
{
  temp_ans.copyFrom(ans);
}

void CopyFromTemp()
{
  ans.copyFrom(temp_ans);
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ans.ansScore = INF;
  best_ans.ansScore = INF;
}

void InitBestPairs()
{
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 1; j < w + 1; ++j) {
      bestPairCount[i][j] = 0;
      for (int k = 0; k < 10; ++k) {
        bestPairValue[i][j][k] = 0;
      }
      int n2 = 0;
      int n2Need = calculateRequiredSize(elementSizes[i][n2], j);
      if (n2Need > w) { continue; }
      for (int n1 = elementCount - 1; n1 >= 0; --n1) {
        if (n2 == elementCount) { break; }
        int val = elementSizes[i][n1];
        int n1Need = (val - 1) / j + 1;
        if (n1Need + n2Need > w) { continue; }
        while (n2 < elementCount - 1) {
          int nxtNeed = calculateRequiredSize(elementSizes[i][n2 + 1], j);
          if (n1Need + nxtNeed <= w) {
            n2++;
            n2Need = nxtNeed;
          }
          else {
            break;
          }
        }
        for (int k = 0; k < 10; ++k) {
          if (n2 - k < 0) { break; }
          if (n2 - k == n1) { continue; }
          int valSum = val + elementSizes[i][n2 - k];
          int junni = bestPairCount[i][j];
          if (junni == 10) {
            if (valSum > bestPairValue[i][j][9]) {
              bestPairValue[i][j][9] = valSum;
              bestPairs[i][j][9][0] = n1;
              bestPairs[i][j][9][1] = n2 - k;
              junni = 9;
            }
          }
          if (junni == 10) { break; }
          while (junni > 0) {
            if (bestPairValue[i][j][junni] > bestPairValue[i][j][junni - 1]) {
              int swa = bestPairValue[i][j][junni - 1];
              bestPairValue[i][j][junni - 1] = bestPairValue[i][j][junni];
              bestPairValue[i][j][junni] = swa;
              swa = bestPairs[i][j][junni - 1][0];
              bestPairs[i][j][junni - 1][0] = bestPairs[i][j][junni][0];
              bestPairs[i][j][junni][0] = swa;
              swa = bestPairs[i][j][junni - 1][1];
              bestPairs[i][j][junni - 1][1] = bestPairs[i][j][junni][1];
              bestPairs[i][j][junni][1] = swa;
              junni--;
            }
            else {
              break;
            }
          }
          if (bestPairCount[i][j] < 10) bestPairCount[i][j]++;
        }
      }
    }
  }
}

// 入力受け取り
void Input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    int www;
    cin >> www >> dayCount >> elementCount;
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < elementCount; ++j) {
        cin >> elementSizes[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    int www;
    ifs >> www >> dayCount >> elementCount;
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < elementCount; ++j) {
        ifs >> elementSizes[i][j];
      }
    }
  }

  for (int j = 0; j < elementCount; ++j) {
    maxElementSize[j] = 0;
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      maxElementSize[j] = max(maxElementSize[j], elementSizes[i][j]);
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    dailyTotalSize[i] = 0;
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      dailyTotalSize[i] += elementSizes[i][j];
    }
  }
  vector<P> sumsPairs;
  for (int i = 0; i < dayCount; ++i) {
    sumsPairs.emplace_back(dailyTotalSize[i], i);
  }
  sort(sumsPairs.begin(), sumsPairs.end());
  for (int i = 0; i < dayCount; ++i) {
    daysDifficultySorted[i] = sumsPairs[i].second;
  }

  int allSum = 0;
  for (int i = 0; i < dayCount; ++i) {
    allSum += dailyTotalSize[i];
  }
  emptySpaceRatio = (double)(w * w * dayCount - allSum) / (w * w * dayCount);
}

// 出力ファイルストリームオープン
void OpenOfs(int case_num, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < elementCount; ++j) {
        for (int k = 0; k < 4; ++k) {
          cout << ans.ans[i][j][k] << ' ';
        }
        cout << endl;
      }
    }
  }
  else {
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < elementCount; ++j) {
        for (int k = 0; k < 4; ++k) {
          ofs << ans.ans[i][j][k] << ' ';
        }
        ofs << endl;
      }
    }
  }
}

bool IsInvalidSolution()
{
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      if (ans.ans[i][j][0] < 0) return true;
      if (ans.ans[i][j][1] < 0) return true;
      if (ans.ans[i][j][0] >= ans.ans[i][j][2]) return true;
      if (ans.ans[i][j][1] >= ans.ans[i][j][3]) return true;
      if (ans.ans[i][j][2] > w) return true;
      if (ans.ans[i][j][3] > w) return true;
      if ((ans.ans[i][j][2] - ans.ans[i][j][0]) * (ans.ans[i][j][3] - ans.ans[i][j][1]) < elementSizes[i][j]) return true;
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      for (int k = j + 1; k < elementCount; ++k) {
        if (max(ans.ans[i][j][0], ans.ans[i][k][0]) < min(ans.ans[i][j][2], ans.ans[i][k][2])) {
          if (max(ans.ans[i][j][1], ans.ans[i][k][1]) < min(ans.ans[i][j][3], ans.ans[i][k][3])) { return true; }
        }
      }
    }
  }
  return false;
}

// スコア計算
ll CalcScore()
{
  ll totalScore = 1;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      auto rect = ans.getRect(i, j);
      totalScore += (ll)max(0, elementSizes[i][j] - (rect.endX - rect.startX) * (rect.endY - rect.startY)) * SCORE_MULTIPLIER;
    }
  }
  {
    vector<P> before;
    for (int i = 0; i < dayCount; ++i) {
      vector<P> now;
      for (int j = 0; j < elementCount; ++j) {
        auto rect = ans.getRect(i, j);
        if (rect.startY != 0 && rect.startY != w) {
          now.emplace_back(rect.startY * w + rect.startX, 1);
          now.emplace_back(rect.startY * w + rect.endX, -1);
        }
        if (rect.endY != 0 && rect.endY != w) {
          now.emplace_back(rect.endY * w + rect.startX, 1);
          now.emplace_back(rect.endY * w + rect.endX, -1);
        }
      }
      sort(now.begin(), now.end());
      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < before.size() || itr2 < now.size()) {
          if (itr1 == before.size()) {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
          else if (itr2 == now.size()) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else if (before[itr1].first <= now[itr2].first) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
        }
      }
      before = now;
    }
  }
  {
    vector<P> before;
    for (int i = 0; i < dayCount; ++i) {
      vector<P> now;
      for (int j = 0; j < elementCount; ++j) {
        auto rect = ans.getRect(i, j);
        if (rect.startX != 0 && rect.startX != w) {
          now.emplace_back(rect.startX * w + rect.startY, 1);
          now.emplace_back(rect.startX * w + rect.endY, -1);
        }
        if (rect.endX != 0 && rect.endX != w) {
          now.emplace_back(rect.endX * w + rect.startY, 1);
          now.emplace_back(rect.endX * w + rect.endY, -1);
        }
      }
      sort(now.begin(), now.end());
      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < before.size() || itr2 < now.size()) {
          if (itr1 == before.size()) {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
          else if (itr2 == now.size()) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else if (before[itr1].first <= now[itr2].first) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
        }
      }
      before = now;
    }
  }
  return totalScore;
}

std::array<int, 1000> CalcScoreForMethod3BeforeArr;
std::array<int, 1000> CalcScoreForMethod3NowArr;
std::array<int, 1000> CalcScoreForMethod3Used;
ll CalcScoreForMethod3()
{
  ll totalScore = 1;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      auto rect = ans.getRect(i, j);
      totalScore += (ll)max(0, elementSizes[i][j] - (rect.endX - rect.startX) * (rect.endY - rect.startY)) * SCORE_MULTIPLIER;
    }
  }

  // 縦線
  {
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < ans.ansLineCount[i] + 1; ++j) {
        if (j == 0 || j == ans.ansLineCount[i]) {
          CalcScoreForMethod3Used[j] = 1;
        }
        else {
          CalcScoreForMethod3Used[j] = 0;
        }
      }
      for (int j = 0; j < elementCount; ++j) {
        for (int k = 0; k < ans.ansLineCount[i]; ++k) {
          if (ans.ans[i][j][1] == ans.ansLinePos[i][k]) {
            CalcScoreForMethod3Used[k] = 1;
            CalcScoreForMethod3Used[k + 1] = 1;
          }
        }
      }
      for (int j = 0; j < ans.ansLineCount[i] + 1; ++j) {
        if (CalcScoreForMethod3Used[j] == 0) { totalScore += w * 2; }
      }
    }
  }

  // 横線
  {
    int beforeTail = 0;
    for (int i = 0; i < dayCount; ++i) {
      int nowTail = 0;
      for (int j = 0; j < elementCount; ++j) {
        auto rect = ans.getRect(i, j);
        if (rect.startX != 0 && rect.startX != w) {
          CalcScoreForMethod3NowArr[nowTail] = (rect.startX * w + rect.startY) * 10 + 1;
          nowTail++;
          CalcScoreForMethod3NowArr[nowTail] = (rect.startX * w + rect.endY) * 10 - 1;
          nowTail++;
        }
        if (rect.endX != 0 && rect.endX != w) {
          CalcScoreForMethod3NowArr[nowTail] = (rect.endX * w + rect.startY) * 10 + 1;
          nowTail++;
          CalcScoreForMethod3NowArr[nowTail] = (rect.endX * w + rect.endY) * 10 - 1;
          nowTail++;
        }
      }

      sort(CalcScoreForMethod3NowArr.begin(), CalcScoreForMethod3NowArr.begin() + nowTail);

      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < beforeTail || itr2 < nowTail) {
          if (itr1 == beforeTail) {
            int val = (CalcScoreForMethod3NowArr[itr2] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += CalcScoreForMethod3NowArr[itr2] - val * 10;
            itr2++;
          }
          else if (itr2 == nowTail) {
            int val = (CalcScoreForMethod3BeforeArr[itr1] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += CalcScoreForMethod3BeforeArr[itr1] - val * 10;
            itr1++;
          }
          else if ((CalcScoreForMethod3BeforeArr[itr1] + 1) / 10 <= (CalcScoreForMethod3NowArr[itr2] + 1) / 10) {
            int val = (CalcScoreForMethod3BeforeArr[itr1] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt1 += CalcScoreForMethod3BeforeArr[itr1] - val * 10;
            itr1++;
          }
          else {
            int val = (CalcScoreForMethod3NowArr[itr2] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { totalScore += (ll)val - pos; }
            pos = val;
            cnt2 += CalcScoreForMethod3NowArr[itr2] - val * 10;
            itr2++;
          }
        }
      }
      beforeTail = nowTail;
      for (int j = 0; j < nowTail; ++j) {
        CalcScoreForMethod3BeforeArr[j] = CalcScoreForMethod3NowArr[j];
      }
    }
  }
  return totalScore;
}

void Initialize()
{
  for (int i = 0; i < dayCount; ++i) {
    int now = 0;
    for (int j = 0; j < elementCount; ++j) {
      int need = calculateRequiredSize(elementSizes[i][j], w);
      int newY = min(now + need, w - (elementCount - 1 - j));
      ans.ans[i][j][0] = 0;
      ans.ans[i][j][1] = now;
      ans.ans[i][j][2] = w;
      ans.ans[i][j][3] = newY;
      now = newY;
    }
  }

  ans.ansScore = CalcScore();
  CopyToBestAns();
}

void Method1()
{
  int ok = 1;
  for (int i = 0; i < dayCount; ++i) {
    int x = 0;
    int y = 0;
    std::array<int, MAX_N> used = {};
    for (int j = 0; j < elementCount; ++j) {
      int minAmari = INT_INF;
      int bestJ = -1;
      int bestX2 = -1;
      int bestY2 = -1;
      int newx = -1;
      int newy = -1;
      for (int jj = 0; jj < elementCount; ++jj) {
        if (used[jj]) { continue; }
        int lenx = w - x;
        int leny = w - y;
        if (elementSizes[i][jj] > lenx * leny) { continue; }
        int need1 = calculateRequiredSize(elementSizes[i][jj], leny);
        int need2 = calculateRequiredSize(elementSizes[i][jj], lenx);
        int amari1 = need1 * leny - elementSizes[i][jj];
        int amari2 = need2 * lenx - elementSizes[i][jj];
        if (amari1 <= minAmari) {
          minAmari = amari1;
          bestJ = jj;
          bestX2 = x + need1;
          bestY2 = w;
          newx = x + need1;
          newy = y;
        }
        if (amari2 <= minAmari) {
          minAmari = amari2;
          bestJ = jj;
          bestX2 = w;
          bestY2 = y + need2;
          newx = x;
          newy = y + need2;
        }
      }
      if (bestJ == -1) {
        ok = 0;
        break;
      }
      ans.ans[i][bestJ][0] = x;
      ans.ans[i][bestJ][1] = y;
      ans.ans[i][bestJ][2] = bestX2;
      ans.ans[i][bestJ][3] = bestY2;
      used[bestJ] = 1;
      x = newx;
      y = newy;
    }
    if (ok == 0) { break; }
  }
  if (ok) {
    ans.ansScore = CalcScore();
    if (ans.ansScore < best_ans.ansScore) { CopyToBestAns(); }
  }
}

void MethodPerfect()
{
  ll maxASum = 0;
  for (int i = 0; i < elementCount; ++i) {
    maxASum += maxElementSize[i];
  }
  if (maxASum > w * w) { return; }

  int ok = 1;
  for (int i = 0; i < dayCount; ++i) {
    int x = 0;
    int y = 0;
    std::array<int, MAX_N> used = {};
    for (int j = 0; j < elementCount; ++j) {
      int minAmari = INT_INF;
      int bestJ = -1;
      int bestX2 = -1;
      int bestY2 = -1;
      int newx = -1;
      int newy = -1;
      for (int jj = 0; jj < elementCount; ++jj) {
        if (used[jj]) { continue; }
        int lenx = w - x;
        int leny = w - y;
        if (maxElementSize[jj] > lenx * leny) {
          ok = 0;
          break;
        }
        int need1 = (maxElementSize[jj] - 1) / leny + 1;
        int need2 = (maxElementSize[jj] - 1) / lenx + 1;
        int amari1 = need1 * leny - maxElementSize[jj];
        int amari2 = need2 * lenx - maxElementSize[jj];
        if (amari1 <= minAmari) {
          minAmari = amari1;
          bestJ = jj;
          bestX2 = x + need1;
          bestY2 = w;
          newx = x + need1;
          newy = y;
        }
        if (amari2 <= minAmari) {
          minAmari = amari2;
          bestJ = jj;
          bestX2 = w;
          bestY2 = y + need2;
          newx = x;
          newy = y + need2;
        }
      }
      if (ok == 0) { break; }
      if (bestJ >= 0) {
        ans.ans[i][bestJ][0] = x;
        ans.ans[i][bestJ][1] = y;
        ans.ans[i][bestJ][2] = bestX2;
        ans.ans[i][bestJ][3] = bestY2;
        used[bestJ] = 1;
        x = newx;
        y = newy;
      }
    }
    if (ok == 0) { break; }
  }

  if (ok) {
    ans.ansScore = CalcScore();
    CopyToBestAns();
  }
}

int preCalcScheduleSizes[MAX_D][MAX_N][MAX_LINECOUNT];
std::array<int, MAX_LINECOUNT> widths = {};

ColumnSchedule columnSchedule;
ColumnSchedule best_columnSchedule;

// realに格納
void CopyToBest_M42()
{
  best_columnSchedule.copyFrom(columnSchedule);
}

// realから戻す
void CoptToCurrent_M42()
{
  columnSchedule.copyFrom(best_columnSchedule);
}

// スコアを更新して最適解を保存の実装
inline void updateScoreAndSaveBest(ll diffScore)
{
  ans.ansScore -= diffScore;
  if (ans.ansScore < best_ans.ansScore) {
    CopyToBestAns();
    CopyToBest_M42();
  }
}

// 列内の要素位置を再計算の実装
inline void recalculateColumnPositions(int day, int column)
{
  columnSchedule.schedulesPosition[day][column][0] = 0;
  for (int k = 0; k < columnSchedule.schedulesCount[day][column]; k++) {
    int element = columnSchedule.schedules[day][column][k];
    columnSchedule.schedulesPosition[day][column][k + 1] =
      columnSchedule.schedulesPosition[day][column][k] + preCalcScheduleSizes[day][element][column];
  }
  columnSchedule.schedulesPosition[day][column][columnSchedule.schedulesCount[day][column]] = w;
}

int CalcAdjacentDayScore(int day1, int day2, int lineNum)
{
  int diffScore2 = 0;
  int ite1 = 1, ite2 = 1;
  while (ite1 < columnSchedule.schedulesCount[day1][lineNum] && ite2 < columnSchedule.schedulesCount[day2][lineNum]) {
    int num1 = columnSchedule.schedulesPosition[day1][lineNum][ite1];
    int num2 = columnSchedule.schedulesPosition[day2][lineNum][ite2];
    if (num1 == num2) {
      diffScore2 -= widths[lineNum] * 2;
      ite1++;
      ite2++;
    }
    else if (num1 < num2) {
      ite1++;
    }
    else {
      ite2++;
    }
  }
  return diffScore2;
}

// 隣接日の差分スコアを計算の実装
inline int calculateAdjacentDaysDiffScore(int day, int lineNum)
{
  int score = 0;
  if (day > 0) {
    score += CalcAdjacentDayScore(day, day - 1, lineNum);
  }
  if (day < dayCount - 1) {
    score += CalcAdjacentDayScore(day, day + 1, lineNum);
  }
  return score;
}

// 要素が指定された幅に収まるかチェックの実装
inline bool canFitInColumn(int day, int element, int width)
{
  return width * w >= elementSizes[day][element];
}

// 要素を配置するのに必要な余白を計算の実装
inline int calculateMargin(int day, int element, int width, int height)
{
  return width - ((elementSizes[day][element] - 1) / height + 1);
}

// 横線のカウントを計算の実装
inline int calculateHorizontalLineCount(int day, int column)
{
  if (day == 0 || day == dayCount - 1) {
    return max(0, columnSchedule.schedulesCount[day][column] - 1);
  }

  int count = max(0, columnSchedule.schedulesCount[day - 1][column] - 1) +
    max(0, columnSchedule.schedulesCount[day][column] - 1);

  // 隣接日で同じ位置にある要素を除外
  int ite1 = 1, ite2 = 1;
  while (ite1 < columnSchedule.schedulesCount[day - 1][column] &&
    ite2 < columnSchedule.schedulesCount[day][column]) {
    if (columnSchedule.schedulesPosition[day - 1][column][ite1] ==
      columnSchedule.schedulesPosition[day][column][ite2]) {
      count -= 2;
      ite1++;
      ite2++;
    }
    else if (columnSchedule.schedulesPosition[day - 1][column][ite1] <
      columnSchedule.schedulesPosition[day][column][ite2]) {
      ite1++;
    }
    else {
      ite2++;
    }
  }

  return count;
}

// 要素の交換による差分スコアを計算の実装
inline int calculateSwapDiffScore(int day, int elem1, int elem2)
{
  int line1 = columnSchedule.columnNum[day][elem1];
  int line2 = columnSchedule.columnNum[day][elem2];

  if (line1 == line2) return 0;

  int idx1 = columnSchedule.findElementIndex(day, line1, elem1);
  int idx2 = columnSchedule.findElementIndex(day, line2, elem2);

  int height1 = columnSchedule.schedulesPosition[day][line1][idx1 + 1] -
    columnSchedule.schedulesPosition[day][line1][idx1];
  int height2 = columnSchedule.schedulesPosition[day][line2][idx2 + 1] -
    columnSchedule.schedulesPosition[day][line2][idx2];

  return (preCalcScheduleSizes[day][elem1][line1] - preCalcScheduleSizes[day][elem1][line2]) * widths[line1] +
    (preCalcScheduleSizes[day][elem2][line2] - preCalcScheduleSizes[day][elem2][line1]) * widths[line2];
}

// 要素に対して最適な配置位置を見つけるの実装
inline PlacementResult findBestPlacement(int day, int element, const std::array<int, MAX_LINECOUNT>& currentPos, int startLine, int endLine)
{
  PlacementResult result = { -1, 0, INT_INF, INT_INF, false };

  for (int k = endLine - 1; k >= startLine; k--) {
    int width = ans.ansLinePos[day][k + 1] - ans.ansLinePos[day][k];
    if (!canFitInColumn(day, element, width)) { break; }

    int need = calculateRequiredSize(elementSizes[day][element], width);
    int over = 0;
    int amari = need * width - elementSizes[day][element];

    if (currentPos[k] + need > w) {
      over = (currentPos[k] + need - w) * width;
    }

    if (over < result.over || (over == result.over && amari <= result.amari)) {
      result.position = k;
      result.need = need;
      result.amari = amari;
      result.over = over;
      result.found = true;
    }
  }

  return result;
}

int shuffleMembers[MAX_N];
int shuffleTmpPosition[MAX_N];
int shuffleNeighborPos[MAX_N * 2];
constexpr std::array<std::array<int, 3>, 6> shuffleThree = { { {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0} } };
int CalcDiffScore3(int iter, int memberCount, int neighborPosCount, int day, int lineNum, int margin)
{
  int diffScore3 = 0;
  if (memberCount == 0) { return diffScore3; }
  if (iter >= 1) {
    if (memberCount == 1) {
      ;
    }
    else if (memberCount == 2) {
      int swa = shuffleMembers[0];
      shuffleMembers[0] = shuffleMembers[1];
      shuffleMembers[1] = swa;
    }
    else if (memberCount == 3) {
      int tmp[3] = {};
      for (int j = 0; j < 3; ++j) {
        tmp[shuffleThree[iter - 1][j]] = shuffleMembers[j];
      }
      for (int j = 0; j < 3; ++j) {
        shuffleMembers[j] = tmp[shuffleThree[iter][j]];
      }
    }
    else {
      FisherYates(shuffleMembers, memberCount);
    }
  }
  int tmpMargin = margin;
  int ite1 = 0;
  int sum = preCalcScheduleSizes[day][shuffleMembers[0]][lineNum];
  shuffleTmpPosition[0] = 0;
  int ite2 = 1;
  shuffleTmpPosition[ite2] = sum;
  while (ite1 < neighborPosCount && ite2 < memberCount) {
    int pos1 = shuffleNeighborPos[ite1];
    if (pos1 < sum) {
      ite1++;
    }
    else if (sum + tmpMargin < pos1) {
      sum += preCalcScheduleSizes[day][shuffleMembers[ite2]][lineNum];
      ite2++;
      shuffleTmpPosition[ite2] = sum;
    }
    else {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      if (ite1 < neighborPosCount) {
        int nxtPos1 = shuffleNeighborPos[ite1];
        if (nxtPos1 == pos1) {
          diffScore3 += widths[lineNum] * 2;
          ite1++;
        }
      }
      tmpMargin -= pos1 - sum;
      sum = pos1;
      shuffleTmpPosition[ite2] = sum;
      sum += preCalcScheduleSizes[day][shuffleMembers[ite2]][lineNum];
      ite2++;
      shuffleTmpPosition[ite2] = sum;
    }
  }
  while (ite2 < memberCount) {
    sum += preCalcScheduleSizes[day][shuffleMembers[ite2]][lineNum];
    ite2++;
    shuffleTmpPosition[ite2] = sum;
  }
  shuffleTmpPosition[ite2] = w;
  return diffScore3;
}

int shuffle2Members[MAX_N];
int shuffle2TmpPosition[MAX_N];
int shuffle2NeighborNewPos1[MAX_N];
int shuffle2NeighborNewPos2[MAX_N];

inline void advanceSum(int& sum, int& ite1, int day, int lineNum) {
  sum += preCalcScheduleSizes[day][shuffle2Members[ite1]][lineNum];
  ite1++;
  shuffle2TmpPosition[ite1] = sum;
}

inline void updateSumAndMargin(int& tmpMargin, int& sum, int& ite1, int pos, int day, int lineNum) {
  tmpMargin -= pos - sum;
  sum = pos;
  shuffle2TmpPosition[ite1] = sum;
  advanceSum(sum, ite1, day, lineNum);
}

inline bool checkDayMinus2Connection(int day, int lineNum, int pos, int& diffScore3) {
  if (1 < day) {
    for (int k = 0; k < columnSchedule.schedulesCount[day - 2][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[day - 2][lineNum][k] == pos) {
        diffScore3 -= widths[lineNum] * 2;
        return true;
      }
    }
  }
  return false;
}

inline void setPos2IfDayMinus2Connected(int day, int lineNum, int& pos2, int pos22) {
  if (1 < day) {
    for (int k = 0; k < columnSchedule.schedulesCount[day - 2][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[day - 2][lineNum][k] == pos22) {
        pos2 = pos22;
        break;
      }
    }
  }
}

inline void addScoreIfDayMinus2Connected(int day, int lineNum, int pos, int& diffScore3) {
  if (1 < day) {
    for (int k = 0; k < columnSchedule.schedulesCount[day - 2][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[day - 2][lineNum][k] == pos) {
        diffScore3 += widths[lineNum] * 2;
        break;
      }
    }
  }
}

inline void updatePos3Common(int& diffScore3, int& tmpMargin, int& sum, int& ite1, int& ite3, int pos3, int day, int lineNum)
{
  shuffle2NeighborNewPos2[ite3] = pos3;
  if (day < elementCount - 2) {
    for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos3) {
        diffScore3 += widths[lineNum] * 2;
        break;
      }
    }
  }
  tmpMargin -= pos3 - sum;
  sum = pos3;
  shuffle2TmpPosition[ite1] = sum;
  advanceSum(sum, ite1, day, lineNum);
}

inline void updatePos2Common(int& diffScore3, int& tmpMargin, int& sum, int& ite1, int& ite2, int pos2, int day, int lineNum) {
  shuffle2NeighborNewPos1[ite2] = pos2;
  addScoreIfDayMinus2Connected(day, lineNum, pos2, diffScore3);
  tmpMargin -= pos2 - sum;
  sum = pos2;
  shuffle2TmpPosition[ite1] = sum;
  advanceSum(sum, ite1, day, lineNum);
}

int CalcDiffScore3_2(int iter, int memberCount, int day, int lineNum, int margin)
{
  int beforeCount = 0;
  if (day > 0) {
    beforeCount = columnSchedule.schedulesCount[day - 1][lineNum];
    for (int k = 0; k < beforeCount + 1; ++k) {
      shuffle2NeighborNewPos1[k] = columnSchedule.schedulesPosition[day - 1][lineNum][k];
    }
  }
  int afterCount = 0;
  if (day < dayCount - 1) {
    afterCount = columnSchedule.schedulesCount[day + 1][lineNum];
    for (int k = 0; k < afterCount + 1; ++k) {
      shuffle2NeighborNewPos2[k] = columnSchedule.schedulesPosition[day + 1][lineNum][k];
    }
  }

  int diffScore3 = 0;

  if (memberCount == 0) { return diffScore3; }
  if (iter >= 1) {
    if (memberCount == 1) {
      ;
    }
    else if (memberCount == 2) {
      int swa = shuffle2Members[0];
      shuffle2Members[0] = shuffle2Members[1];
      shuffle2Members[1] = swa;
    }
    else if (memberCount == 3) {
      int tmp[3] = {};
      for (int j = 0; j < 3; ++j) {
        tmp[shuffleThree[iter - 1][j]] = shuffle2Members[j];
      }
      for (int j = 0; j < 3; ++j) {
        shuffle2Members[j] = tmp[shuffleThree[iter][j]];
      }
    }
    else {
      FisherYates(shuffle2Members, memberCount);
    }
  }
  int tmpMargin = margin;
  int ite1 = 1;
  int sum = preCalcScheduleSizes[day][shuffle2Members[0]][lineNum];
  shuffle2TmpPosition[0] = 0;
  shuffle2TmpPosition[ite1] = sum;

  int ite2 = 1;
  int ite3 = 1;

  while (ite1 < memberCount) {
    int num2 = 0, pos2 = 0, pos22 = 0;
    if (ite2 < beforeCount) {
      num2 = columnSchedule.schedules[day - 1][lineNum][ite2 - 1];
      pos2 = columnSchedule.schedulesPosition[day - 1][lineNum][ite2 - 1] + preCalcScheduleSizes[day - 1][num2][lineNum];
      pos22 = columnSchedule.schedulesPosition[day - 1][lineNum][ite2];
      pos2 = min(pos22, max(pos2, sum + 1));
      pos2 = pos22;
      setPos2IfDayMinus2Connected(day, lineNum, pos2, pos22);
    }
    int num3 = 0, pos3 = 0, pos32 = 0;
    if (ite3 < afterCount) {
      num3 = columnSchedule.schedules[day + 1][lineNum][ite3 - 1];
      pos3 = columnSchedule.schedulesPosition[day + 1][lineNum][ite3 - 1] + preCalcScheduleSizes[day + 1][num3][lineNum];
      pos32 = columnSchedule.schedulesPosition[day + 1][lineNum][ite3];
      pos3 = min(pos32, max(pos3, sum + 1));
      pos3 = pos32;
      if (day < elementCount - 2) {
        for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
          if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos32) {
            pos3 = pos32;
            break;
          }
        }
      }
    }
    if (ite2 < beforeCount && ite3 < afterCount) {
      if (pos22 < sum) {
        ite2++;
        continue;
      }
      if (pos32 < sum) {
        ite3++;
        continue;
      }
      if (sum + tmpMargin < pos2 && sum + tmpMargin < pos3) {
        advanceSum(sum, ite1, day, lineNum);
        continue;
      }
      if (pos2 <= sum + tmpMargin && pos3 <= sum + tmpMargin) {
        if (pos22 == pos32) {
          bool isConnect1 = checkDayMinus2Connection(day, lineNum, pos22, diffScore3);
          bool isConnect2 = false;
          if (day < elementCount - 2) {
            for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
              if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos32) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect2 = true;
                break;
              }
            }
          }
          if (pos22 <= sum + tmpMargin && pos32 <= sum + tmpMargin && (isConnect1 || isConnect2)) {
            // nPos = pos22 = pos32
            if (isConnect1) diffScore3 += widths[lineNum] * 2;
            if (isConnect2) diffScore3 += widths[lineNum] * 2;
            updateSumAndMargin(tmpMargin, sum, ite1, pos22, day, lineNum);
            ite2++;
            ite3++;
          }
          else {
            // nPos = min(pos2,pos3)
            if (pos2 <= pos3) {
              // nPos = pos2
              updatePos2Common(diffScore3, tmpMargin, sum, ite1, ite2, pos2, day, lineNum);
              ite2++;
            }
            else {
              // nPos = pos3
              updatePos3Common(diffScore3, tmpMargin, sum, ite1, ite3, pos3, day, lineNum);
              ite3++;
            }
          }
        }
        else if (pos22 < pos32) {
          bool isConnect = checkDayMinus2Connection(day, lineNum, pos22, diffScore3);
          if (pos22 <= sum + tmpMargin && isConnect) {
            // nPos = pos22
            if (isConnect) diffScore3 += widths[lineNum] * 2;
            updateSumAndMargin(tmpMargin, sum, ite1, pos22, day, lineNum);
          }
          else {
            // nPos = pos2
            updatePos2Common(diffScore3, tmpMargin, sum, ite1, ite2, pos2, day, lineNum);
          }
          ite2++;
        }
        else {
          bool isConnect = false;
          if (day < elementCount - 2) {
            for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
              if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos32) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect = true;
                break;
              }
            }
          }
          if (pos32 <= sum + tmpMargin && isConnect) {
            // nPos = pos32
            if (isConnect) diffScore3 += widths[lineNum] * 2;
            updateSumAndMargin(tmpMargin, sum, ite1, pos32, day, lineNum);
          }
          else {
            // nPos = pos3
            updatePos3Common(diffScore3, tmpMargin, sum, ite1, ite3, pos3, day, lineNum);
          }
          ite3++;
        }
      }
      else if (pos2 <= sum + tmpMargin) {
        if (pos22 < sum) {
          ite2++;
          continue;
        }
        if (sum + tmpMargin < pos2) {
          advanceSum(sum, ite1, day, lineNum);
          continue;
        }
        bool isConnect = checkDayMinus2Connection(day, lineNum, pos22, diffScore3);
        if (isConnect) {
          if (pos22 <= sum + tmpMargin) {
            diffScore3 += widths[lineNum] * 2;
            updateSumAndMargin(tmpMargin, sum, ite1, pos22, day, lineNum);
          }
          else {
            updatePos2Common(diffScore3, tmpMargin, sum, ite1, ite2, pos2, day, lineNum);
          }
        }
        else {
          updatePos2Common(diffScore3, tmpMargin, sum, ite1, ite2, pos2, day, lineNum);
        }
        ite2++;
      }
      else {
        if (pos32 < sum) {
          ite3++;
          continue;
        }
        if (sum + tmpMargin < pos3) {
          advanceSum(sum, ite1, day, lineNum);
          continue;
        }
        bool isConnect = false;
        if (day < elementCount - 2) {
          for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
            if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos32) {
              diffScore3 -= widths[lineNum] * 2;
              isConnect = true;
              break;
            }
          }
        }
        if (isConnect) {
          if (pos32 <= sum + tmpMargin) {
            diffScore3 += widths[lineNum] * 2;
            updateSumAndMargin(tmpMargin, sum, ite1, pos32, day, lineNum);
          }
          else {
            updatePos3Common(diffScore3, tmpMargin, sum, ite1, ite3, pos3, day, lineNum);
          }
        }
        else {
          updatePos3Common(diffScore3, tmpMargin, sum, ite1, ite3, pos3, day, lineNum);
        }
        ite3++;
      }
    }
    else if (ite2 < beforeCount) {
      if (pos22 < sum) {
        ite2++;
        continue;
      }
      if (sum + tmpMargin < pos2) {
        advanceSum(sum, ite1, day, lineNum);
        continue;
      }
      bool isConnect = checkDayMinus2Connection(day, lineNum, pos22, diffScore3);
      if (isConnect && pos22 <= sum + tmpMargin) {
        diffScore3 += widths[lineNum] * 2;
        updateSumAndMargin(tmpMargin, sum, ite1, pos22, day, lineNum);
      }
      else {
        updatePos2Common(diffScore3, tmpMargin, sum, ite1, ite2, pos2, day, lineNum);
      }
      ite2++;
    }
    else if (ite3 < afterCount) {
      if (pos32 < sum) {
        ite3++;
        continue;
      }
      if (sum + tmpMargin < pos3) {
        advanceSum(sum, ite1, day, lineNum);
        continue;
      }

      bool isConnect = false;
      if (day < elementCount - 2) {
        for (int k = 0; k < columnSchedule.schedulesCount[day + 2][lineNum]; ++k) {
          if (columnSchedule.schedulesPosition[day + 2][lineNum][k] == pos32) {
            diffScore3 -= widths[lineNum] * 2;
            isConnect = true;
            break;
          }
        }
      }

      if (sum + tmpMargin < pos32 || !isConnect) {
        updatePos3Common(diffScore3, tmpMargin, sum, ite1, ite3, pos3, day, lineNum);
      }
      else {
        if (isConnect) diffScore3 += widths[lineNum] * 2;
        updateSumAndMargin(tmpMargin, sum, ite1, pos32, day, lineNum);
      }

      ite3++;
    }
    else {
      advanceSum(sum, ite1, day, lineNum);
    }
  }

  shuffle2TmpPosition[memberCount] = w;

  ite1 = 1;
  ite2 = 1;
  ite3 = 1;
  while (ite1 < memberCount && ite2 < beforeCount) {
    if (shuffle2TmpPosition[ite1] == shuffle2NeighborNewPos1[ite2]) {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      ite2++;
    }
    else if (shuffle2TmpPosition[ite1] < shuffle2NeighborNewPos1[ite2]) {
      ite1++;
    }
    else {
      ite2++;
    }
  }

  ite1 = 1;
  ite2 = 1;
  ite3 = 1;
  while (ite1 < memberCount && ite3 < afterCount) {
    if (shuffle2TmpPosition[ite1] == shuffle2NeighborNewPos2[ite3]) {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      ite3++;
    }
    else if (shuffle2TmpPosition[ite1] < shuffle2NeighborNewPos2[ite3]) {
      ite1++;
    }
    else {
      ite3++;
    }
  }

  return diffScore3;
}

double swap_start_temp = ANNEALING_START_TEMP;
double swap_end_temp = ANNEALING_END_TEMP;
double swap_startTime = 0;
double swap_nowTime = 0;
double swap_timeLimit = 0;
int swapMasks[32] = {};

std::vector<int> swapCandidates;
std::vector<int> swapNonCandidates;
int swapNeighborPosCount = 0;

int swap_tmpAnsNextLine[MAX_N];
int swap_tmpAnsNextLinePosition[MAX_N];
int swap_tmpAnsNextCount = 0;
int swap_tmpAnsCurrentLine[MAX_N];
int swap_tmpAnsCurrentLinePosition[MAX_N];
int swap_tmpAnsCurrentCount = 0;

int swap_tmpAnsNextLineBeforePosition[MAX_N];
int swap_tmpAnsNextLineAfterPosition[MAX_N];
int swap_tmpAnsCurrentLineBeforePosition[MAX_N];
int swap_tmpAnsCurrentLineAfterPosition[MAX_N];

// 1対多スワップ
int swapCount;
void Method4_3_1()
{
  int d = Rand() % dayCount;
  int n = Rand() % elementCount;
  int lineNum = columnSchedule.columnNum[d][n];
  int lineCapacity = w;
  for (int k = 0; k < columnSchedule.schedulesCount[d][lineNum]; ++k) {
    int num = columnSchedule.schedules[d][lineNum][k];
    if (num != n) { lineCapacity -= preCalcScheduleSizes[d][num][lineNum]; }
  }

  int nextLine = Rand() % ans.ansLineCount[d];
  while (nextLine == lineNum) {
    nextLine = Rand() % ans.ansLineCount[d];
  }
  int beforeCount = columnSchedule.schedulesCount[d][lineNum];
  int afterCount = columnSchedule.schedulesCount[d][nextLine];
  swapCandidates.clear();
  swapCandidates.reserve(MAX_CANDIDATES);
  swapNonCandidates.clear();
  swapNonCandidates.reserve(MAX_CANDIDATES);
  int nextLineCapacity = w;
  int maxNextLineCapacity = w;
  for (int k = 0; k < columnSchedule.schedulesCount[d][nextLine]; ++k) {
    int num = columnSchedule.schedules[d][nextLine][k];
    nextLineCapacity -= preCalcScheduleSizes[d][num][nextLine];
    if (preCalcScheduleSizes[d][num][lineNum] > lineCapacity) {
      maxNextLineCapacity -= preCalcScheduleSizes[d][num][nextLine];
      swapNonCandidates.push_back(num);
    }
    else {
      swapCandidates.push_back(num);
    }
  }

  if (swapCandidates.size() > MAX_CANDIDATES) { return; }
  if (preCalcScheduleSizes[d][n][nextLine] > maxNextLineCapacity) { return; }

  int maskCount = 0;
  if (swapCandidates.size() == 0) {
    swapMasks[maskCount] = 0;
    maskCount++;
  }
  else if (swapCandidates.size() <= SMALL_DATA_THRESHOLD) {
    for (int mask = 0; mask < (1 << swapCandidates.size()); ++mask) {
      // if (mask == 0) { continue; }
      swapMasks[maskCount] = mask;
      maskCount++;
    }
  }
  else {
    for (int trial = 0; trial < RANDOM_TRIAL_COUNT; ++trial) {
      int mask = Rand() % ((1 << swapCandidates.size()));
      swapMasks[maskCount] = mask;
      maskCount++;
    }
  }
  FisherYates(swapMasks, maskCount);
  for (int idx = 0; idx < maskCount; ++idx) {
    int mask = swapMasks[idx];
    int nextSpace = 0;
    int needSpace = 0;
    for (int jj = 0; jj < swapCandidates.size(); ++jj) {
      if (mask & (1 << jj)) {
        int j = swapCandidates[jj];
        nextSpace += preCalcScheduleSizes[d][j][nextLine];
        needSpace += preCalcScheduleSizes[d][j][lineNum];
        if (needSpace > lineCapacity) { break; }
      }
    }
    if (needSpace <= lineCapacity && preCalcScheduleSizes[d][n][nextLine] <= nextLineCapacity + nextSpace) {
      // スワップ可能
      int moveCount = 0;
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if (mask & (1 << jj)) { moveCount++; }
      }
      int diffScore1 = 0;
      if (swapCandidates.size() == 0 || mask == 0) {
        if (beforeCount > 1) {
          if (d == 0 || d == dayCount - 1) {
            diffScore1 += widths[lineNum];
          }
          else {
            diffScore1 += widths[lineNum] * 2;
          }
        }
        if (afterCount > 0) {
          if (d == 0 || d == dayCount - 1) {
            diffScore1 -= widths[nextLine];
          }
          else {
            diffScore1 -= widths[nextLine] * 2;
          }
        }
      }
      else {
        diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]) * 2;
        if (d == 0 || d == dayCount - 1) { diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]); }
      }

      int diffScore2 = calculateAdjacentDaysDiffScore(d, lineNum) + calculateAdjacentDaysDiffScore(d, nextLine);

      // 10回シャッフルnextLine
      swapNeighborPosCount = 0;
      if (d > 0) {
        for (int k = 1; k < columnSchedule.schedulesCount[d - 1][nextLine]; ++k) {
          shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d - 1][nextLine][k];
          swapNeighborPosCount++;
        }
      }
      if (d < dayCount - 1) {
        for (int k = 1; k < columnSchedule.schedulesCount[d + 1][nextLine]; ++k) {
          shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d + 1][nextLine][k];
          swapNeighborPosCount++;
        }
      }
      sort(shuffleNeighborPos, shuffleNeighborPos + swapNeighborPosCount);
      swapNonCandidates.push_back(n);
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if ((mask & (1 << jj)) == 0) {
          int j = swapCandidates[jj];
          swapNonCandidates.push_back(j);
        }
      }
      for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
        shuffleMembers[i] = swapNonCandidates[i];
      }
      int margin = w;
      for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
        int num = swapNonCandidates[i];
        margin -= preCalcScheduleSizes[d][num][nextLine];
      }
      int diffScore3 = -INT_INF;
      for (int iter = 0; iter < 10; ++iter) {
        if (swapNonCandidates.size() == 1 && iter >= 1) { break; }
        if (swapNonCandidates.size() == 2 && iter >= 2) { break; }
        if (swapNonCandidates.size() == 3 && iter >= 6) { break; }
        int tmpDiffScore3 = CalcDiffScore3(iter, swapNonCandidates.size(), swapNeighborPosCount, d, nextLine, margin);
        if (tmpDiffScore3 > diffScore3) {
          diffScore3 = tmpDiffScore3;
          for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
            swap_tmpAnsNextLine[i] = shuffleMembers[i];
            swap_tmpAnsNextLinePosition[i] = shuffleTmpPosition[i];
          }
          swap_tmpAnsNextLinePosition[swapNonCandidates.size()] = shuffleTmpPosition[swapNonCandidates.size()];
          swap_tmpAnsNextCount = swapNonCandidates.size();
        }
      }

      // 10回シャッフルlineNum
      swapNeighborPosCount = 0;
      if (d > 0) {
        for (int k = 1; k < columnSchedule.schedulesCount[d - 1][lineNum]; ++k) {
          shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d - 1][lineNum][k];
          swapNeighborPosCount++;
        }
      }
      if (d < dayCount - 1) {
        for (int k = 1; k < columnSchedule.schedulesCount[d + 1][lineNum]; ++k) {
          shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d + 1][lineNum][k];
          swapNeighborPosCount++;
        }
      }
      sort(shuffleNeighborPos, shuffleNeighborPos + swapNeighborPosCount);
      int tmpKouhos[MAX_N];
      int tmpKouhoCount = 0;
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if (mask & (1 << jj)) {
          int j = swapCandidates[jj];
          tmpKouhos[tmpKouhoCount] = j;
          tmpKouhoCount++;
        }
      }
      for (int j = 0; j < tmpKouhoCount; ++j) {
        swapCandidates[j] = tmpKouhos[j];
      }
      swapCandidates.resize(tmpKouhoCount);
      for (int k = 0; k < columnSchedule.schedulesCount[d][lineNum]; ++k) {
        int num = columnSchedule.schedules[d][lineNum][k];
        if (num != n) {
          swapCandidates.push_back(num);
        }
      }
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        shuffleMembers[i] = swapCandidates[i];
      }
      margin = w;
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        int num = swapCandidates[i];
        margin -= preCalcScheduleSizes[d][num][lineNum];
      }
      int diffScore4 = -INT_INF;
      for (int iter = 0; iter < 10; ++iter) {
        if (shouldBreakShuffle(swapCandidates.size(), iter)) { break; }
        int tmpDiffScore4 = CalcDiffScore3(iter, swapCandidates.size(), swapNeighborPosCount, d, lineNum, margin);
        if (tmpDiffScore4 > diffScore4) {
          diffScore4 = tmpDiffScore4;
          for (int i = 0; i < (int)swapCandidates.size(); ++i) {
            swap_tmpAnsCurrentLine[i] = shuffleMembers[i];
            swap_tmpAnsCurrentLinePosition[i] = shuffleTmpPosition[i];
          }
          swap_tmpAnsCurrentLinePosition[swapCandidates.size()] = shuffleTmpPosition[swapCandidates.size()];
          swap_tmpAnsCurrentCount = swapCandidates.size();
        }
      }

      int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

      double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * swap_nowTime / swap_timeLimit);
      const double prob = exp((double)totalDiffScore * ANNEALING_TEMP_COEFFICIENT / temp);

      if (prob > Rand01()) {
        swapCount++;
        columnSchedule.schedulesCount[d][nextLine] = swap_tmpAnsNextCount;
        for (int i = 0; i < swap_tmpAnsNextCount; ++i) {
          int num = swap_tmpAnsNextLine[i];
          columnSchedule.columnNum[d][num] = nextLine;
          columnSchedule.schedules[d][nextLine][i] = num;
          columnSchedule.schedulesPosition[d][nextLine][i] = swap_tmpAnsNextLinePosition[i];
        }
        columnSchedule.schedulesPosition[d][nextLine][swap_tmpAnsNextCount] = swap_tmpAnsNextLinePosition[swap_tmpAnsNextCount];

        columnSchedule.schedulesCount[d][lineNum] = swap_tmpAnsCurrentCount;
        for (int i = 0; i < swap_tmpAnsCurrentCount; ++i) {
          int num = swap_tmpAnsCurrentLine[i];
          columnSchedule.columnNum[d][num] = lineNum;
          columnSchedule.schedules[d][lineNum][i] = num;
          columnSchedule.schedulesPosition[d][lineNum][i] = swap_tmpAnsCurrentLinePosition[i];
        }
        columnSchedule.schedulesPosition[d][lineNum][swap_tmpAnsCurrentCount] = swap_tmpAnsCurrentLinePosition[swap_tmpAnsCurrentCount];

        updateScoreAndSaveBest(totalDiffScore);
      }
      else {
        ;
      }

      break;
    }
  }
}

// 1対多スワップ（柔軟）
void Method4_3_2()
{
  int d = Rand() % dayCount;
  int n = Rand() % elementCount;
  int lineNum = columnSchedule.columnNum[d][n];
  int lineCapacity = w;
  for (int k = 0; k < columnSchedule.schedulesCount[d][lineNum]; ++k) {
    int num = columnSchedule.schedules[d][lineNum][k];
    if (num != n) { lineCapacity -= preCalcScheduleSizes[d][num][lineNum]; }
  }

  int nextLine = Rand() % ans.ansLineCount[d];
  while (nextLine == lineNum) {
    nextLine = Rand() % ans.ansLineCount[d];
  }
  int beforeCount = columnSchedule.schedulesCount[d][lineNum];
  int afterCount = columnSchedule.schedulesCount[d][nextLine];
  swapCandidates.clear();
  swapCandidates.reserve(MAX_CANDIDATES);
  swapNonCandidates.clear();
  swapNonCandidates.reserve(MAX_CANDIDATES);
  int nextLineCapacity = w;
  int maxNextLineCapacity = w;
  for (int k = 0; k < columnSchedule.schedulesCount[d][nextLine]; ++k) {
    int num = columnSchedule.schedules[d][nextLine][k];
    nextLineCapacity -= preCalcScheduleSizes[d][num][nextLine];
    if (preCalcScheduleSizes[d][num][lineNum] > lineCapacity) {
      maxNextLineCapacity -= preCalcScheduleSizes[d][num][nextLine];
      swapNonCandidates.push_back(num);
    }
    else {
      swapCandidates.push_back(num);
    }
  }

  if (swapCandidates.size() > MAX_CANDIDATES) { return; }
  if (swapCandidates.size() == 0) { return; }
  if (preCalcScheduleSizes[d][n][nextLine] > maxNextLineCapacity) { return; }

  int maskCount = 0;
  if (swapCandidates.size() == 0) {
    swapMasks[maskCount] = 0;
    maskCount++;
  }
  else if (swapCandidates.size() <= SMALL_DATA_THRESHOLD) {
    for (int mask = 0; mask < (1 << swapCandidates.size()); ++mask) {
      if (mask == 0) { continue; }
      swapMasks[maskCount] = mask;
      maskCount++;
    }
  }
  else {
    for (int trial = 0; trial < RANDOM_TRIAL_COUNT; ++trial) {
      int mask = Rand() % ((1 << swapCandidates.size()) - 1) + 1;
      swapMasks[maskCount] = mask;
      maskCount++;
    }
  }
  FisherYates(swapMasks, maskCount);
  for (int idx = 0; idx < maskCount; ++idx) {
    int mask = swapMasks[idx];
    int nextSpace = 0;
    int needSpace = 0;
    for (int jj = 0; jj < swapCandidates.size(); ++jj) {
      if (mask & (1 << jj)) {
        int j = swapCandidates[jj];
        nextSpace += preCalcScheduleSizes[d][j][nextLine];
        needSpace += preCalcScheduleSizes[d][j][lineNum];
        if (needSpace > lineCapacity) { break; }
      }
    }
    if (needSpace <= lineCapacity && preCalcScheduleSizes[d][n][nextLine] <= nextLineCapacity + nextSpace) {
      // スワップ可能
      int moveCount = 0;
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if (mask & (1 << jj)) { moveCount++; }
      }
      int diffScore1 = 0;
      if (swapCandidates.size() == 0) {
        if (beforeCount > 1) {
          diffScore1 = widths[lineNum] * 2;
          if (d == 0 || d == dayCount - 1) { diffScore1 = widths[lineNum]; }
        }
      }
      else {
        diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]) * 2;
        if (d == 0 || d == dayCount - 1) { diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]); }
      }

      int diffScore2 = calculateAdjacentDaysDiffScore(d, lineNum) + calculateAdjacentDaysDiffScore(d, nextLine);

      // 10回シャッフルnextLine
      swapNonCandidates.push_back(n);
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if ((mask & (1 << jj)) == 0) {
          int j = swapCandidates[jj];
          swapNonCandidates.push_back(j);
        }
      }
      for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
        shuffle2Members[i] = swapNonCandidates[i];
      }
      int margin = w;
      for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
        int num = swapNonCandidates[i];
        margin -= preCalcScheduleSizes[d][num][nextLine];
      }
      int diffScore3 = -INT_INF;
      for (int iter = 0; iter < 10; ++iter) {
        if (swapNonCandidates.size() == 1 && iter >= 1) { break; }
        if (swapNonCandidates.size() == 2 && iter >= 2) { break; }
        if (swapNonCandidates.size() == 3 && iter >= 6) { break; }
        int tmpDiffScore3 = CalcDiffScore3_2(iter, swapNonCandidates.size(), d, nextLine, 0);
        if (tmpDiffScore3 > diffScore3) {
          diffScore3 = tmpDiffScore3;
          for (int i = 0; i < (int)swapNonCandidates.size(); ++i) {
            swap_tmpAnsNextLine[i] = shuffle2Members[i];
            swap_tmpAnsNextLinePosition[i] = shuffle2TmpPosition[i];
          }
          swap_tmpAnsNextLinePosition[swapNonCandidates.size()] = shuffle2TmpPosition[swapNonCandidates.size()];
          swap_tmpAnsNextCount = swapNonCandidates.size();
        }

        if (d > 0) {
          beforeCount = columnSchedule.schedulesCount[d - 1][nextLine];
          for (int k = 0; k < beforeCount + 1; ++k) {
            swap_tmpAnsNextLineBeforePosition[k] = shuffle2NeighborNewPos1[k];
          }
        }
        if (d < dayCount - 1) {
          int afterCount = columnSchedule.schedulesCount[d + 1][nextLine];
          for (int k = 0; k < afterCount + 1; ++k) {
            swap_tmpAnsNextLineAfterPosition[k] = shuffle2NeighborNewPos2[k];
          }
        }
      }

      // 10回シャッフルlineNum
      int tmpKouhos[MAX_N];
      int tmpKouhoCount = 0;
      for (int jj = 0; jj < swapCandidates.size(); ++jj) {
        if (mask & (1 << jj)) {
          int j = swapCandidates[jj];
          tmpKouhos[tmpKouhoCount] = j;
          tmpKouhoCount++;
        }
      }
      for (int j = 0; j < tmpKouhoCount; ++j) {
        swapCandidates[j] = tmpKouhos[j];
      }
      swapCandidates.resize(tmpKouhoCount);
      for (int k = 0; k < columnSchedule.schedulesCount[d][lineNum]; ++k) {
        int num = columnSchedule.schedules[d][lineNum][k];
        if (num != n) {
          swapCandidates.push_back(num);
        }
      }
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        shuffle2Members[i] = swapCandidates[i];
      }
      margin = w;
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        int num = swapCandidates[i];
        margin -= preCalcScheduleSizes[d][num][lineNum];
      }
      int diffScore4 = -INT_INF;
      for (int iter = 0; iter < 10; ++iter) {
        if (shouldBreakShuffle(swapCandidates.size(), iter)) { break; }
        int tmpDiffScore4 = CalcDiffScore3_2(iter, swapCandidates.size(), d, lineNum, 0);
        if (tmpDiffScore4 > diffScore4) {
          diffScore4 = tmpDiffScore4;
          for (int i = 0; i < (int)swapCandidates.size(); ++i) {
            swap_tmpAnsCurrentLine[i] = shuffle2Members[i];
            swap_tmpAnsCurrentLinePosition[i] = shuffle2TmpPosition[i];
          }
          swap_tmpAnsCurrentLinePosition[swapCandidates.size()] = shuffle2TmpPosition[swapCandidates.size()];
          swap_tmpAnsCurrentCount = swapCandidates.size();

          if (d > 0) {
            beforeCount = columnSchedule.schedulesCount[d - 1][lineNum];
            for (int k = 0; k < beforeCount + 1; ++k) {
              swap_tmpAnsCurrentLineBeforePosition[k] = shuffle2NeighborNewPos1[k];
            }
          }
          if (d < dayCount - 1) {
            int afterCount = columnSchedule.schedulesCount[d + 1][lineNum];
            for (int k = 0; k < afterCount + 1; ++k) {
              swap_tmpAnsCurrentLineAfterPosition[k] = shuffle2NeighborNewPos2[k];
            }
          }
        }
      }

      int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

      double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * swap_nowTime / swap_timeLimit);
      const double prob = exp((double)totalDiffScore * ANNEALING_TEMP_COEFFICIENT / temp);

      if (prob > Rand01()) {
        columnSchedule.schedulesCount[d][nextLine] = swap_tmpAnsNextCount;
        for (int i = 0; i < swap_tmpAnsNextCount; ++i) {
          int num = swap_tmpAnsNextLine[i];
          columnSchedule.columnNum[d][num] = nextLine;
          columnSchedule.schedules[d][nextLine][i] = num;
          columnSchedule.schedulesPosition[d][nextLine][i] = swap_tmpAnsNextLinePosition[i];
        }
        columnSchedule.schedulesPosition[d][nextLine][swap_tmpAnsNextCount] = swap_tmpAnsNextLinePosition[swap_tmpAnsNextCount];

        columnSchedule.schedulesCount[d][lineNum] = swap_tmpAnsCurrentCount;
        for (int i = 0; i < swap_tmpAnsCurrentCount; ++i) {
          int num = swap_tmpAnsCurrentLine[i];
          columnSchedule.columnNum[d][num] = lineNum;
          columnSchedule.schedules[d][lineNum][i] = num;
          columnSchedule.schedulesPosition[d][lineNum][i] = swap_tmpAnsCurrentLinePosition[i];
        }
        columnSchedule.schedulesPosition[d][lineNum][swap_tmpAnsCurrentCount] = swap_tmpAnsCurrentLinePosition[swap_tmpAnsCurrentCount];

        if (d > 0) {
          beforeCount = columnSchedule.schedulesCount[d - 1][lineNum];
          for (int k = 0; k < beforeCount + 1; ++k) {
            columnSchedule.schedulesPosition[d - 1][lineNum][k] = swap_tmpAnsCurrentLineBeforePosition[k];
          }
        }
        if (d < dayCount - 1) {
          int afterCount = columnSchedule.schedulesCount[d + 1][lineNum];
          for (int k = 0; k < afterCount + 1; ++k) {
            columnSchedule.schedulesPosition[d + 1][lineNum][k] = swap_tmpAnsCurrentLineAfterPosition[k];
          }
        }
        if (d > 0) {
          beforeCount = columnSchedule.schedulesCount[d - 1][nextLine];
          for (int k = 0; k < beforeCount + 1; ++k) {
            columnSchedule.schedulesPosition[d - 1][nextLine][k] = swap_tmpAnsNextLineBeforePosition[k];
          }
        }
        if (d < dayCount - 1) {
          int afterCount = columnSchedule.schedulesCount[d + 1][nextLine];
          for (int k = 0; k < afterCount + 1; ++k) {
            columnSchedule.schedulesPosition[d + 1][nextLine][k] = swap_tmpAnsNextLineAfterPosition[k];
          }
        }

        updateScoreAndSaveBest(totalDiffScore);
      }
      else {
        ;
      }

      break;
    }
  }
}

// 列スワップ
void Method4_3_3()
{
  int d = Rand() % dayCount;
  int line1 = Rand() % ans.ansLineCount[d];
  int line2 = Rand() % ans.ansLineCount[d];
  while (line2 == line1) {
    line2 = Rand() % ans.ansLineCount[d];
  }
  if (columnSchedule.schedulesCount[d][line1] == 0 || columnSchedule.schedulesCount[d][line2] == 0) { return; }

  {
    int line1NextSum = 0;
    for (int k = 0; k < columnSchedule.schedulesCount[d][line1]; ++k) {
      int num = columnSchedule.schedules[d][line1][k];
      line1NextSum += preCalcScheduleSizes[d][num][line2];
    }
    if (line1NextSum > w) { return; }
  }
  {
    int line2NextSum = 0;
    for (int k = 0; k < columnSchedule.schedulesCount[d][line2]; ++k) {
      int num = columnSchedule.schedules[d][line2][k];
      line2NextSum += preCalcScheduleSizes[d][num][line1];
    }
    if (line2NextSum > w) { return; }
  }

  int diffScore1 = (columnSchedule.schedulesCount[d][line1] - columnSchedule.schedulesCount[d][line2]) * (widths[line1] - widths[line2]) * 2;
  if (d == 0 || d == dayCount - 1) diffScore1 = (columnSchedule.schedulesCount[d][line1] - columnSchedule.schedulesCount[d][line2]) * (widths[line1] - widths[line2]);

  int diffScore2 = 0;
  if (d > 0) {
    diffScore2 += CalcAdjacentDayScore(d, d - 1, line1);
    diffScore2 += CalcAdjacentDayScore(d, d - 1, line2);
  }
  if (d < dayCount - 1) {
    diffScore2 += CalcAdjacentDayScore(d, d + 1, line1);
    diffScore2 += CalcAdjacentDayScore(d, d + 1, line2);
  }

  // 10回シャッフルline2
  swapNeighborPosCount = 0;
  if (d > 0) {
    for (int k = 1; k < columnSchedule.schedulesCount[d - 1][line2]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d - 1][line2][k];
      swapNeighborPosCount++;
    }
  }
  if (d < dayCount - 1) {
    for (int k = 1; k < columnSchedule.schedulesCount[d + 1][line2]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d + 1][line2][k];
      swapNeighborPosCount++;
    }
  }
  sort(shuffleNeighborPos, shuffleNeighborPos + swapNeighborPosCount);

  swapCandidates.clear();
  swapCandidates.reserve(columnSchedule.schedulesCount[d][line1]);
  for (int k = 0; k < columnSchedule.schedulesCount[d][line1]; ++k) {
    swapCandidates.push_back(columnSchedule.schedules[d][line1][k]);
  }
  int margin = w;
  for (int i = 0; i < (int)swapCandidates.size(); ++i) {
    shuffleMembers[i] = columnSchedule.schedules[d][line1][i];
    margin -= preCalcScheduleSizes[d][shuffleMembers[i]][line2];
  }

  int diffScore3 = -1;
  for (int iter = 0; iter < 10; ++iter) {
    if (shouldBreakShuffle(swapCandidates.size(), iter)) { break; }
    int tmpDiffScore3 = CalcDiffScore3(iter, swapCandidates.size(), swapNeighborPosCount, d, line2, margin);
    if (tmpDiffScore3 > diffScore3) {
      diffScore3 = tmpDiffScore3;
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        swap_tmpAnsNextLine[i] = shuffleMembers[i];
        swap_tmpAnsNextLinePosition[i] = shuffleTmpPosition[i];
      }
      swap_tmpAnsNextLinePosition[swapCandidates.size()] = shuffleTmpPosition[swapCandidates.size()];
      swap_tmpAnsNextCount = swapCandidates.size();
    }
  }

  // 10回シャッフルline1
  swapNeighborPosCount = 0;
  if (d > 0) {
    for (int k = 1; k < columnSchedule.schedulesCount[d - 1][line1]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d - 1][line1][k];
      swapNeighborPosCount++;
    }
  }
  if (d < dayCount - 1) {
    for (int k = 1; k < columnSchedule.schedulesCount[d + 1][line1]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d + 1][line1][k];
      swapNeighborPosCount++;
    }
  }
  sort(shuffleNeighborPos, shuffleNeighborPos + swapNeighborPosCount);

  swapCandidates.clear();
  swapCandidates.reserve(columnSchedule.schedulesCount[d][line2]);
  for (int k = 0; k < columnSchedule.schedulesCount[d][line2]; ++k) {
    swapCandidates.push_back(columnSchedule.schedules[d][line2][k]);
  }
  margin = w;
  for (int i = 0; i < (int)swapCandidates.size(); ++i) {
    shuffleMembers[i] = columnSchedule.schedules[d][line2][i];
    margin -= preCalcScheduleSizes[d][shuffleMembers[i]][line1];
  }

  int diffScore4 = -1;
  for (int iter = 0; iter < 10; ++iter) {
    if (shouldBreakShuffle(swapCandidates.size(), iter)) { break; }
    int tmpDiffScore4 = CalcDiffScore3(iter, swapCandidates.size(), swapNeighborPosCount, d, line1, margin);
    if (tmpDiffScore4 > diffScore4) {
      diffScore4 = tmpDiffScore4;
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        swap_tmpAnsCurrentLine[i] = shuffleMembers[i];
        swap_tmpAnsCurrentLinePosition[i] = shuffleTmpPosition[i];
      }
      swap_tmpAnsCurrentLinePosition[swapCandidates.size()] = shuffleTmpPosition[swapCandidates.size()];
      swap_tmpAnsCurrentCount = swapCandidates.size();
    }
  }

  int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

  double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * swap_nowTime / swap_timeLimit);
  const double prob = exp((double)totalDiffScore * 10000 / temp);

  if (totalDiffScore >= 0) {
    // if (prob > Rand01()) {
    // cout << "OK";
    columnSchedule.schedulesCount[d][line2] = swap_tmpAnsNextCount;
    for (int i = 0; i < swap_tmpAnsNextCount; ++i) {
      int num = swap_tmpAnsNextLine[i];
      columnSchedule.columnNum[d][num] = line2;
      columnSchedule.schedules[d][line2][i] = num;
      columnSchedule.schedulesPosition[d][line2][i] = swap_tmpAnsNextLinePosition[i];
    }
    columnSchedule.schedulesPosition[d][line2][swap_tmpAnsNextCount] = swap_tmpAnsNextLinePosition[swap_tmpAnsNextCount];

    columnSchedule.schedulesCount[d][line1] = swap_tmpAnsCurrentCount;
    for (int i = 0; i < swap_tmpAnsCurrentCount; ++i) {
      int num = swap_tmpAnsCurrentLine[i];
      columnSchedule.columnNum[d][num] = line1;
      columnSchedule.schedules[d][line1][i] = num;
      columnSchedule.schedulesPosition[d][line1][i] = swap_tmpAnsCurrentLinePosition[i];
    }
    columnSchedule.schedulesPosition[d][line1][swap_tmpAnsCurrentCount] = swap_tmpAnsCurrentLinePosition[swap_tmpAnsCurrentCount];

    ans.ansScore -= totalDiffScore;
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
      CopyToBest_M42();
    }
  }
}

// 横線を移動
void Method4_3_4_2()
{
  int d = Rand() % dayCount;
  int lineNum = Rand() % ans.ansBaseLineCount;
  if (columnSchedule.schedulesCount[d][lineNum] <= 1) { return; }
  int raIndex = Rand() % columnSchedule.schedulesCount[d][lineNum];
  int dir = Rand() % 2;
  if (raIndex == 0) { dir = 1; }
  if (raIndex == columnSchedule.schedulesCount[d][lineNum] - 1) { dir = 0; }
  int beforeLinePos = columnSchedule.schedulesPosition[d][lineNum][raIndex];
  if (dir == 1) { beforeLinePos = columnSchedule.schedulesPosition[d][lineNum][raIndex + 1]; }

  int margin = w;
  int startDay = d;
  int endDay = d;
  for (int i = d; i >= 0; --i) {
    int lineIndex = -1;
    for (int j = 1; j < columnSchedule.schedulesCount[i][lineNum]; ++j) {
      if (columnSchedule.schedulesPosition[i][lineNum][j] == beforeLinePos) {
        lineIndex = j;
        break;
      }
    }
    if (lineIndex == -1) { break; }
    startDay = i;
    if (dir == 0) {
      int n = columnSchedule.schedules[i][lineNum][lineIndex];
      margin = min(margin, (columnSchedule.schedulesPosition[i][lineNum][lineIndex + 1] - columnSchedule.schedulesPosition[i][lineNum][lineIndex]) - preCalcScheduleSizes[i][n][lineNum]);
    }
    else {
      int n = columnSchedule.schedules[i][lineNum][lineIndex - 1];
      margin = min(margin, (columnSchedule.schedulesPosition[i][lineNum][lineIndex] - columnSchedule.schedulesPosition[i][lineNum][lineIndex - 1]) - preCalcScheduleSizes[i][n][lineNum]);
    }
    if (margin == 0) { return; }
  }
  for (int i = d + 1; i < dayCount; ++i) {
    int lineIndex = -1;
    for (int j = 1; j < columnSchedule.schedulesCount[i][lineNum]; ++j) {
      if (columnSchedule.schedulesPosition[i][lineNum][j] == beforeLinePos) {
        lineIndex = j;
        break;
      }
    }
    if (lineIndex == -1) { break; }
    endDay = i;
    if (dir == 0) {
      int n = columnSchedule.schedules[i][lineNum][lineIndex];
      margin = min(margin, (columnSchedule.schedulesPosition[i][lineNum][lineIndex + 1] - columnSchedule.schedulesPosition[i][lineNum][lineIndex]) - preCalcScheduleSizes[i][n][lineNum]);
    }
    else {
      int n = columnSchedule.schedules[i][lineNum][lineIndex - 1];
      margin = min(margin, (columnSchedule.schedulesPosition[i][lineNum][lineIndex] - columnSchedule.schedulesPosition[i][lineNum][lineIndex - 1]) - preCalcScheduleSizes[i][n][lineNum]);
    }
    if (margin == 0) { return; }
  }

  int moveAmount = Rand() % margin + 1;
  int afterLinePos = beforeLinePos + moveAmount;
  if (dir == 1) { afterLinePos = beforeLinePos - moveAmount; }

  int diffScore = 0;
  if (startDay > 0) {
    for (int k = 0; k < columnSchedule.schedulesCount[startDay - 1][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[startDay - 1][lineNum][k] == beforeLinePos) diffScore -= widths[lineNum] * 2;
      if (columnSchedule.schedulesPosition[startDay - 1][lineNum][k] == afterLinePos) diffScore += widths[lineNum] * 2;
    }
  }
  if (endDay < dayCount - 1) {
    for (int k = 0; k < columnSchedule.schedulesCount[endDay + 1][lineNum]; ++k) {
      if (columnSchedule.schedulesPosition[endDay + 1][lineNum][k] == beforeLinePos) diffScore -= widths[lineNum] * 2;
      if (columnSchedule.schedulesPosition[endDay + 1][lineNum][k] == afterLinePos) diffScore += widths[lineNum] * 2;
    }
  }

  if (diffScore >= 0) {
    for (int i = startDay; i < endDay + 1; ++i) {
      int ok = 0;
      for (int j = 1; j < columnSchedule.schedulesCount[i][lineNum]; ++j) {
        if (columnSchedule.schedulesPosition[i][lineNum][j] == beforeLinePos) {
          columnSchedule.schedulesPosition[i][lineNum][j] = afterLinePos;
          ok = 1;
          break;
        }
      }
    }
    ans.ansScore -= diffScore;
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
      CopyToBest_M42();
    }
  }
}

void Method4_3_5()
{
  int lineNum = Rand() % ans.ansBaseLineCount;
  int dir = Rand() % 2;
  if (lineNum == 0) { dir = 1; }
  if (lineNum == ans.ansBaseLineCount - 1) { dir = 0; }
  int margin = widths[lineNum] - 1;
  for (int i = 0; i < dayCount; ++i) {
    for (int k = 0; k < columnSchedule.schedulesCount[i][lineNum]; ++k) {
      int num = columnSchedule.schedules[i][lineNum][k];
      int tmpMargin = widths[lineNum] - ((elementSizes[i][num] - 1) / (columnSchedule.schedulesPosition[i][lineNum][k + 1] - columnSchedule.schedulesPosition[i][lineNum][k]) + 1);
      if (tmpMargin < margin) margin = tmpMargin;
    }
    if (margin <= 0) { break; }
  }
  if (margin == 0) { return; }
  int moveAmount = Rand() % margin + 1;
  int diffScore = 0;
  for (int i = 1; i < dayCount; ++i) {
    int ite1 = 1;
    int ite2 = 1;
    while (ite1 < columnSchedule.schedulesCount[i - 1][lineNum] || ite2 < columnSchedule.schedulesCount[i][lineNum]) {
      if (ite1 == columnSchedule.schedulesCount[i - 1][lineNum]) {
        diffScore += moveAmount;
        ite2++;
      }
      else if (ite2 == columnSchedule.schedulesCount[i][lineNum]) {
        diffScore += moveAmount;
        ite1++;
      }
      else {
        if (columnSchedule.schedulesPosition[i - 1][lineNum][ite1] == columnSchedule.schedulesPosition[i][lineNum][ite2]) {
          ite1++;
          ite2++;
        }
        else if (columnSchedule.schedulesPosition[i - 1][lineNum][ite1] < columnSchedule.schedulesPosition[i][lineNum][ite2]) {
          diffScore += moveAmount;
          ite1++;
        }
        else {
          diffScore += moveAmount;
          ite2++;
        }
      }
    }
    if (dir == 0) {
      ite1 = 1;
      ite2 = 1;
      while (ite1 < columnSchedule.schedulesCount[i - 1][lineNum - 1] || ite2 < columnSchedule.schedulesCount[i][lineNum - 1]) {
        if (ite1 == columnSchedule.schedulesCount[i - 1][lineNum - 1]) {
          diffScore -= moveAmount;
          ite2++;
        }
        else if (ite2 == columnSchedule.schedulesCount[i][lineNum - 1]) {
          diffScore -= moveAmount;
          ite1++;
        }
        else {
          if (columnSchedule.schedulesPosition[i - 1][lineNum - 1][ite1] == columnSchedule.schedulesPosition[i][lineNum - 1][ite2]) {
            ite1++;
            ite2++;
          }
          else if (columnSchedule.schedulesPosition[i - 1][lineNum - 1][ite1] < columnSchedule.schedulesPosition[i][lineNum - 1][ite2]) {
            diffScore -= moveAmount;
            ite1++;
          }
          else {
            diffScore -= moveAmount;
            ite2++;
          }
        }
      }
    }
    else {
      ite1 = 1;
      ite2 = 1;
      while (ite1 < columnSchedule.schedulesCount[i - 1][lineNum + 1] || ite2 < columnSchedule.schedulesCount[i][lineNum + 1]) {
        if (ite1 == columnSchedule.schedulesCount[i - 1][lineNum + 1]) {
          diffScore -= moveAmount;
          ite2++;
        }
        else if (ite2 == columnSchedule.schedulesCount[i][lineNum + 1]) {
          diffScore -= moveAmount;
          ite1++;
        }
        else {
          if (columnSchedule.schedulesPosition[i - 1][lineNum + 1][ite1] == columnSchedule.schedulesPosition[i][lineNum + 1][ite2]) {
            ite1++;
            ite2++;
          }
          else if (columnSchedule.schedulesPosition[i - 1][lineNum + 1][ite1] < columnSchedule.schedulesPosition[i][lineNum + 1][ite2]) {
            diffScore -= moveAmount;
            ite1++;
          }
          else {
            diffScore -= moveAmount;
            ite2++;
          }
        }
      }
    }
  }

  double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * swap_nowTime / swap_timeLimit);
  const double prob = exp((double)diffScore * 1 / temp);

  if (prob > Rand01()) {
    if (dir == 0) {
      widths[lineNum] -= moveAmount;
      widths[lineNum - 1] += moveAmount;
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < elementCount; ++j) {
          preCalcScheduleSizes[i][j][lineNum] = calculateRequiredSize(elementSizes[i][j], widths[lineNum]);
          preCalcScheduleSizes[i][j][lineNum - 1] = calculateRequiredSize(elementSizes[i][j], widths[lineNum - 1]);
        }
      }
    }
    else {
      widths[lineNum] -= moveAmount;
      widths[lineNum + 1] += moveAmount;
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < elementCount; ++j) {
          preCalcScheduleSizes[i][j][lineNum] = calculateRequiredSize(elementSizes[i][j], widths[lineNum]);
          preCalcScheduleSizes[i][j][lineNum + 1] = calculateRequiredSize(elementSizes[i][j], widths[lineNum + 1]);
        }
      }
    }
    ans.ansScore -= diffScore;
    for (int i = 0; i < dayCount; ++i) {
      ans.ansLinePos[i][0] = 0;
      for (int j = 0; j < ans.ansBaseLineCount; ++j) {
        ans.ansLinePos[i][j + 1] = ans.ansLinePos[i][j] + widths[j];
      }
    }
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
      CopyToBest_M42();
    }
  }
}

// 交換できるやつを交換
void Method4_3_6()
{
  int d = Rand() % dayCount;
  int n = Rand() % elementCount;
  int lineNum = columnSchedule.columnNum[d][n];
  int lineIndex = 0;
  int nSpace = 0;
  lineIndex = columnSchedule.findElementIndex(d, lineNum, n);
  if (lineIndex != -1) {
    nSpace = (columnSchedule.schedulesPosition[d][lineNum][lineIndex + 1] - columnSchedule.schedulesPosition[d][lineNum][lineIndex]) * widths[lineNum];
  }
  int now = Rand() % elementCount;
  for (int jisoo = 0; jisoo < elementCount; ++jisoo) {
    now = (now + 97) % elementCount;
    if (now == n) { continue; }
    int nextLine = columnSchedule.columnNum[d][now];
    int nextLineIndex = 0;
    int nowSpace = 0;
    for (int i = 0; i < columnSchedule.schedulesCount[d][nextLine]; ++i) {
      if (columnSchedule.schedules[d][nextLine][i] == now) {
        nextLineIndex = i;
        nowSpace = (columnSchedule.schedulesPosition[d][nextLine][i + 1] - columnSchedule.schedulesPosition[d][nextLine][i]) * widths[nextLine];
        break;
      }
    }
    if (elementSizes[d][n] <= nowSpace && elementSizes[d][now] <= nSpace) {
      columnSchedule.columnNum[d][n] = nextLine;
      columnSchedule.columnNum[d][now] = lineNum;
      columnSchedule.schedules[d][lineNum][lineIndex] = now;
      columnSchedule.schedules[d][nextLine][nextLineIndex] = n;
      break;
    }
  }
}

// 列シャッフル
void Method4_3_7()
{
  int d = Rand() % dayCount;
  int n = Rand() % elementCount;
  int lineNum = columnSchedule.columnNum[d][n];
  int lineCapacity = w;
  if (columnSchedule.schedulesCount[d][lineNum] == 1) { return; }

  int diffScore2 = 0;
  if (d > 0) { diffScore2 += CalcAdjacentDayScore(d, d - 1, lineNum); }
  if (d < dayCount - 1) { diffScore2 += CalcAdjacentDayScore(d, d + 1, lineNum); }
  // 10回シャッフルlineNum
  swapNeighborPosCount = 0;
  if (d > 0) {
    for (int k = 1; k < columnSchedule.schedulesCount[d - 1][lineNum]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d - 1][lineNum][k];
      swapNeighborPosCount++;
    }
  }
  if (d < dayCount - 1) {
    for (int k = 1; k < columnSchedule.schedulesCount[d + 1][lineNum]; ++k) {
      shuffleNeighborPos[swapNeighborPosCount] = columnSchedule.schedulesPosition[d + 1][lineNum][k];
      swapNeighborPosCount++;
    }
  }

  swapCandidates.clear();
  swapCandidates.reserve(columnSchedule.schedulesCount[d][lineNum]);
  for (int k = 0; k < columnSchedule.schedulesCount[d][lineNum]; ++k) {
    swapCandidates.push_back(columnSchedule.schedules[d][lineNum][k]);
  }
  int margin = w;
  for (int i = 0; i < (int)swapCandidates.size(); ++i) {
    shuffleMembers[i] = columnSchedule.schedules[d][lineNum][i];
    margin -= preCalcScheduleSizes[d][shuffleMembers[i]][lineNum];
  }

  int diffScore4 = -1;
  for (int iter = 0; iter < 10; ++iter) {
    if (shouldBreakShuffle(swapCandidates.size(), iter)) { break; }
    int tmpDiffScore4 = CalcDiffScore3(iter, swapCandidates.size(), swapNeighborPosCount, d, lineNum, margin);
    if (tmpDiffScore4 > diffScore4) {
      diffScore4 = tmpDiffScore4;
      for (int i = 0; i < (int)swapCandidates.size(); ++i) {
        swap_tmpAnsCurrentLine[i] = shuffleMembers[i];
        swap_tmpAnsCurrentLinePosition[i] = shuffleTmpPosition[i];
      }
      swap_tmpAnsCurrentLinePosition[swapCandidates.size()] = shuffleTmpPosition[swapCandidates.size()];
      swap_tmpAnsCurrentCount = swapCandidates.size();
    }
  }

  int totalDiffScore = diffScore2 + diffScore4;

  double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * swap_nowTime / swap_timeLimit);
  const double prob = exp((double)totalDiffScore * ANNEALING_TEMP_COEFFICIENT / temp);

  if (prob > Rand01()) {
    columnSchedule.schedulesCount[d][lineNum] = swap_tmpAnsCurrentCount;
    for (int i = 0; i < swap_tmpAnsCurrentCount; ++i) {
      int num = swap_tmpAnsCurrentLine[i];
      columnSchedule.columnNum[d][num] = lineNum;
      columnSchedule.schedules[d][lineNum][i] = num;
      columnSchedule.schedulesPosition[d][lineNum][i] = swap_tmpAnsCurrentLinePosition[i];
    }
    columnSchedule.schedulesPosition[d][lineNum][swap_tmpAnsCurrentCount] = swap_tmpAnsCurrentLinePosition[swap_tmpAnsCurrentCount];

    ans.ansScore -= totalDiffScore;
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
      CopyToBest_M42();
    }
  }
  else {
    ;
  }
}

// 幅を不要な列に移動
int horizontalLineCount[MAX_LINECOUNT];
int lineNumbers[MAX_LINECOUNT];
void Method4_3_8()
{
  for (int j = 0; j < ans.ansBaseLineCount; ++j) lineNumbers[j] = j;
  FisherYates(lineNumbers, ans.ansBaseLineCount);
  int minLineNum = -1;
  int minLineCount = INT_INF;
  int taisyouLineNum = lineNumbers[0];
  for (int jj = 0; jj < 2; ++jj) {
    int j = lineNumbers[jj];
    horizontalLineCount[j] = 0;
    for (int i = 0; i < dayCount; ++i) {
      horizontalLineCount[j] += calculateHorizontalLineCount(i, j);
    }
    if (horizontalLineCount[j] < minLineCount) {
      minLineCount = horizontalLineCount[j];
      minLineNum = j;
      if (minLineCount == 0) { break; }
    }
  }
  if (minLineNum == taisyouLineNum) { return; }

  // 対象列の中身をなるべく小さくする
  // marginを計算する
  int margin = widths[taisyouLineNum] - 1;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < columnSchedule.schedulesCount[i][taisyouLineNum]; ++j) {
      int num = columnSchedule.schedules[i][taisyouLineNum][j];
      int height = columnSchedule.schedulesPosition[i][taisyouLineNum][j + 1] - columnSchedule.schedulesPosition[i][taisyouLineNum][j];
      for (int k = 0; k < elementCount; ++k) {
        if (k == num) { break; }
        if (columnSchedule.columnNum[i][k] == taisyouLineNum) { continue; }
        if (preCalcScheduleSizes[i][k][taisyouLineNum] <= height) {
          int nextLine = columnSchedule.columnNum[i][k];
          int nextLineIndex = -1;
          for (int l = 0; l < columnSchedule.schedulesCount[i][nextLine]; ++l) {
            if (columnSchedule.schedules[i][nextLine][l] == k) {
              nextLineIndex = l;
              break;
            }
          }
          if (preCalcScheduleSizes[i][num][nextLine] > columnSchedule.schedulesPosition[i][nextLine][nextLineIndex + 1] - columnSchedule.schedulesPosition[i][nextLine][nextLineIndex]) { continue; }
          swap(columnSchedule.schedules[i][taisyouLineNum][j], columnSchedule.schedules[i][nextLine][nextLineIndex]);
          swap(columnSchedule.columnNum[i][num], columnSchedule.columnNum[i][k]);
          break;
        }
      }
      int newNum = columnSchedule.schedules[i][taisyouLineNum][j];
      // margin計算
      margin = min(margin, calculateMargin(i, newNum, widths[taisyouLineNum], height));
    }
  }

  if (margin == 0) { return; }
  int diffScore = margin * (horizontalLineCount[taisyouLineNum] - minLineCount);

  // ansLinePosとpreCalcScheduleSizesとwidthsを更新
  widths[taisyouLineNum] -= margin;
  widths[minLineNum] += margin;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < ans.ansBaseLineCount; ++j) {
      ans.ansLinePos[i][j + 1] = ans.ansLinePos[i][j] + widths[j];
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      preCalcScheduleSizes[i][j][taisyouLineNum] = calculateRequiredSize(elementSizes[i][j], widths[taisyouLineNum]);
      preCalcScheduleSizes[i][j][minLineNum] = calculateRequiredSize(elementSizes[i][j], widths[minLineNum]);
    }
  }

  ans.ansScore -= diffScore;
  if (ans.ansScore < best_ans.ansScore) {
    CopyToBestAns();
    CopyToBest_M42();
  }
}

// ソート
int sortArray[MAX_N];
int sortLineNum[MAX_N];
int sortLineIndex[MAX_N];
void Method4_3_9()
{
  int d = Rand() % dayCount;
  int cnt = 0;
  for (int j = 0; j < ans.ansBaseLineCount; ++j) {
    for (int k = 0; k < columnSchedule.schedulesCount[d][j]; ++k) {
      sortLineNum[cnt] = j;
      sortLineIndex[cnt] = k;
      sortArray[cnt] = widths[j] * (columnSchedule.schedulesPosition[d][j][k + 1] - columnSchedule.schedulesPosition[d][j][k]) * 100 + cnt;
      cnt++;
    }
  }
  sort(sortArray, sortArray + cnt);

  for (int j = 0; j < ans.ansBaseLineCount; ++j) {
    for (int k = 0; k < columnSchedule.schedulesCount[d][j]; ++k) {
      columnSchedule.schedules[d][j][k] = -1;
    }
  }

  int ite = 0;
  for (int j = 0; j < elementCount; ++j) {
    while (true) {
      if (elementSizes[d][j] <= sortArray[ite] / 100) {
        int posNum = sortArray[ite] % 100;
        int lineNum = sortLineNum[posNum];
        int lineIndex = sortLineIndex[posNum];
        columnSchedule.columnNum[d][j] = lineNum;
        columnSchedule.schedules[d][lineNum][lineIndex] = j;
        ite++;
        break;
      }
      else {
        ite++;
      }
    }
  }
}

void Method4_3(double timeLimit)
{
  swap_startTime = timer.get_elapsed_time();
  swap_timeLimit = timeLimit;
  swapCount = 0;

  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      int ng = 1;
      for (int k = 0; k < ans.ansLineCount[i]; ++k) {
        if (ans.ans[i][j][1] == ans.ansLinePos[i][k] && ans.ans[i][j][3] == ans.ansLinePos[i][k + 1]) {
          columnSchedule.columnNum[i][j] = k;
          break;
        }
      }
    }
  }

  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < ans.ansLineCount[i]; ++k) {
        int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
        if (width == 0) {
          preCalcScheduleSizes[i][j][k] = INVALID_SIZE;
        }
        else {
          preCalcScheduleSizes[i][j][k] = calculateRequiredSize(elementSizes[i][j], width);
        }
      }
    }
  }

  // 初期解作成
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < ans.ansLineCount[i]; ++j) {
      columnSchedule.schedulesCount[i][j] = 0;
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = elementCount - 1; j >= 0; --j) {
      int lineNum = columnSchedule.columnNum[i][j];
      columnSchedule.schedules[i][lineNum][columnSchedule.schedulesCount[i][lineNum]] = j;
      columnSchedule.schedulesCount[i][lineNum]++;
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < ans.ansLineCount[i]; ++j) {
      columnSchedule.schedulesPosition[i][j][0] = 0;
      for (int k = 0; k < columnSchedule.schedulesCount[i][j]; ++k) {
        int num = columnSchedule.schedules[i][j][k];
        columnSchedule.schedulesPosition[i][j][k + 1] = columnSchedule.schedulesPosition[i][j][k] + preCalcScheduleSizes[i][num][j];
        if (k == columnSchedule.schedulesCount[i][j] - 1) { columnSchedule.schedulesPosition[i][j][k + 1] = w; }
      }
    }
  }

  for (int j = 0; j < ans.ansBaseLineCount; ++j) {
    widths[j] = ans.ansLinePos[0][j + 1] - ans.ansLinePos[0][j];
  }

  // realに格納
  CopyToBest_M42();

  int loopCount = 0;
  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      swap_nowTime = timer.get_elapsed_time();
      if (swap_nowTime > timeLimit) {
        break;
      }
    }

    int ra = Rand() % 101;
    if (ra < 0) {
      if (swap_nowTime < (swap_startTime + timeLimit) / 2) {
        continue;
      }
      Method4_3_2();
    }
    else if (ra < 69) {
      Method4_3_1();
    }
    else if (ra < 70) {
      Method4_3_3();
    }
    else if (ra < 90) {
      Method4_3_4_2();
    }
    else if (ra < 91) {
      Method4_3_5();
    }
    else if (ra < 100) {
      Method4_3_6();
    }
    else if (ra < 100) {
      Method4_3_8();
    }
    else if (ra < 100) {
      Method4_3_7();
    }
    else if (ra < 101) {
      if (Rand() % 100 != 0) { continue; }
      Method4_3_9();
    }
  }

  CopyFromBestAns();
  CoptToCurrent_M42();

  // ansScheduleLineNumからansを作成
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < ans.ansBaseLineCount; ++j) {
      for (int k = 0; k < columnSchedule.schedulesCount[i][j]; ++k) {
        int num = columnSchedule.schedules[i][j][k];
        ans.ans[i][num][0] = columnSchedule.schedulesPosition[i][j][k];
        ans.ans[i][num][2] = columnSchedule.schedulesPosition[i][j][k + 1];
        int posNum = columnSchedule.columnNum[i][num];
        ans.ans[i][num][1] = ans.ansLinePos[i][posNum];
        ans.ans[i][num][3] = ans.ansLinePos[i][posNum + 1];
      }
    }
  }

  ans.ansScore = CalcScore();
  CopyToBestAns();

  CopyFromBestAns();
  CoptToCurrent_M42();
}

int oshiiLineCount[MAX_D];
int oshiiLinePos[MAX_D][MAX_LINECOUNT];
int M3_alreadyUsed[MAX_D][MAX_N];
int M3_alreadyCount = 0;
int Method3_Oshii()
{
  for (int i = 0; i < dayCount; ++i) {
    if (oshiiLineCount[i] == -1) { return 0; }
  }
  // 作り直し
  {
    ans.ansBaseLineCount = oshiiLineCount[0];
    for (int i = 0; i < dayCount; ++i) {
      ans.ansLineCount[i] = oshiiLineCount[i];
      for (int j = 0; j < oshiiLineCount[i] + 1; ++j) {
        ans.ansLinePos[i][j] = oshiiLinePos[i][j];
      }
    }

    std::array<int, MAX_LINECOUNT> now = {};
    int ng = 0;
    for (int ii = dayCount - 1; ii >= 0; --ii) {
      int i = daysDifficultySorted[ii];
      for (int j = 0; j < ans.ansLineCount[i]; ++j) {
        now[j] = 0;
      }

      for (int j = elementCount - 1; j >= 0; --j) {
        if (M3_alreadyUsed[i][j]) { continue; }
        int minAmari = INT_INF;
        int minOver = INT_INF;
        int posNum = -1;
        int tmpNeed = 0;
        for (int k = ans.ansLineCount[i] - 1; k >= M3_alreadyCount; --k) {
          int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
          if (!canFitInColumn(i, j, width)) { break; }
          int need = calculateRequiredSize(elementSizes[i][j], width);
          int over = 0;
          int amari = need * width - elementSizes[i][j];
          if (now[k] + need > w) {
            over = (now[k] + need - w) * width;
            amari = 0;
          }

          if (over < minOver || (over == minOver && amari <= minAmari)) {
            minAmari = amari;
            minOver = over;
            posNum = k;
            tmpNeed = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans.ans[i][j][0] = now[posNum];
        ans.ans[i][j][2] = now[posNum] + tmpNeed;
        ans.ans[i][j][1] = ans.ansLinePos[i][posNum];
        ans.ans[i][j][3] = ans.ansLinePos[i][posNum + 1];
        now[posNum] += tmpNeed;
      }
      if (ng) { break; }
    }

    if (ng) {
      CopyFromBestAns();
      return 0;
    }
  }

  // 調整
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < ans.ansLineCount[i]; ++k) {
        if (ans.ans[i][j][1] == ans.ansLinePos[i][k] && ans.ans[i][j][3] == ans.ansLinePos[i][k + 1]) {
          columnSchedule.columnNum[i][j] = k;
          break;
        }
      }
    }
  }

  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < ans.ansLineCount[i]; ++k) {
        int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
        if (width == 0) {
          preCalcScheduleSizes[i][j][k] = INVALID_SIZE;
        }
        else {
          preCalcScheduleSizes[i][j][k] = calculateRequiredSize(elementSizes[i][j], width);
        }
      }
    }
  }

  int ng = 0;
  int kouhos[MAX_N];
  int kouhoCount = 0;
  int loopCount = 0;
  for (int ii = 0; ii < dayCount; ++ii) {
    int i = daysDifficultySorted[ii];
    int ok = 1;
    for (int j = 0; j < elementCount; ++j) {
      if (ans.ans[i][j][2] > w) {
        ok = 0;
        break;
      }
    }
    if (ok) { continue; }
    int nowSum[MAX_LINECOUNT] = {};
    for (int j = 0; j < elementCount; ++j) {
      nowSum[columnSchedule.columnNum[i][j]] = max(nowSum[columnSchedule.columnNum[i][j]], ans.ans[i][j][2]);
    }
    int ngCount = 0;
    for (int j = 0; j < elementCount; ++j) {
      if (nowSum[j] > w) ngCount++;
    }

    // 逆引き作成
    for (int j = 0; j < ans.ansLineCount[i]; ++j) {
      columnSchedule.schedulesCount[i][j] = 0;
    }
    for (int j = 0; j < elementCount; ++j) {
      int lineNum = columnSchedule.columnNum[i][j];
      columnSchedule.schedules[i][lineNum][columnSchedule.schedulesCount[i][lineNum]] = j;
      columnSchedule.schedulesCount[i][lineNum]++;
    }

    for (int leSserafim = 0; leSserafim < 10000; ++leSserafim) {
      int n = Rand() % elementCount;
      int lineNum = columnSchedule.columnNum[i][n];
      int nextLine = Rand() % ans.ansBaseLineCount;
      while (nextLine == lineNum) {
        nextLine = Rand() % ans.ansBaseLineCount;
      }
      kouhoCount = 0;
      if (nowSum[lineNum] <= w && nowSum[nextLine] <= w && preCalcScheduleSizes[i][n][nextLine] > w) { continue; }
      for (int j = 0; j < columnSchedule.schedulesCount[i][nextLine]; ++j) {
        if (nowSum[lineNum] <= w && nowSum[nextLine] <= w && preCalcScheduleSizes[i][columnSchedule.schedules[i][nextLine][j]][lineNum] > w) { continue; }
        kouhos[kouhoCount] = columnSchedule.schedules[i][nextLine][j];
        kouhoCount++;
      }

      if (kouhoCount == 0 || kouhoCount > MAX_CANDIDATES) { continue; }

      if (kouhoCount <= SMALL_DATA_THRESHOLD) {
        for (int chaewon = (1 << kouhoCount) - 1; chaewon >= 0; --chaewon) {
          if (chaewon == 0) { continue; }
          int nextSpace = 0;
          int needSpace = 0;
          for (int jj = 0; jj < kouhoCount; ++jj) {
            if (chaewon & (1 << jj)) {
              int j = kouhos[jj];
              nextSpace += preCalcScheduleSizes[i][j][nextLine];
              needSpace += preCalcScheduleSizes[i][j][lineNum];
            }
          }
          int diffOver = max(0, nowSum[lineNum] - w) + max(0, nowSum[nextLine] - w);
          diffOver -= max(0, nowSum[lineNum] - preCalcScheduleSizes[i][n][lineNum] + needSpace - w) + max(0, nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][n][nextLine] - w);
          if (diffOver >= 0) {
            // スワップ
            if (nowSum[lineNum] > w) ngCount--;
            if (nowSum[nextLine] > w) ngCount--;
            nowSum[lineNum] = nowSum[lineNum] - preCalcScheduleSizes[i][n][lineNum] + needSpace;
            nowSum[nextLine] = nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][n][nextLine];
            if (nowSum[lineNum] > w) ngCount++;
            if (nowSum[nextLine] > w) ngCount++;
            columnSchedule.columnNum[i][n] = nextLine;
            for (int jj = 0; jj < kouhoCount; ++jj) {
              if (chaewon & (1 << jj)) {
                int j = kouhos[jj];
                columnSchedule.columnNum[i][j] = lineNum;
              }
            }

            // 逆引き更新
            for (int j = 0; j < ans.ansLineCount[i]; ++j) {
              columnSchedule.schedulesCount[i][j] = 0;
            }
            for (int j = 0; j < elementCount; ++j) {
              int lineNum = columnSchedule.columnNum[i][j];
              columnSchedule.schedules[i][lineNum][columnSchedule.schedulesCount[i][lineNum]] = j;
              columnSchedule.schedulesCount[i][lineNum]++;
            }

            break;
          }
        }
      }
      else {
        for (int eunchae = 0; eunchae < 32; ++eunchae) {
          int chaewon = Rand() % ((1 << kouhoCount) - 1) + 1;
          if (chaewon == 0) { continue; }
          int nextSpace = 0;
          int needSpace = 0;
          for (int jj = 0; jj < kouhoCount; ++jj) {
            if (chaewon & (1 << jj)) {
              int j = kouhos[jj];
              nextSpace += preCalcScheduleSizes[i][j][nextLine];
              needSpace += preCalcScheduleSizes[i][j][lineNum];
            }
          }
          int diffOver = max(0, nowSum[lineNum] - w) + max(0, nowSum[nextLine] - w);
          diffOver -= max(0, nowSum[lineNum] - preCalcScheduleSizes[i][n][lineNum] + needSpace - w) + max(0, nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][n][nextLine] - w);
          if (diffOver >= 0) {
            // スワップ
            if (nowSum[lineNum] > w) ngCount--;
            if (nowSum[nextLine] > w) ngCount--;
            nowSum[lineNum] = nowSum[lineNum] - preCalcScheduleSizes[i][n][lineNum] + needSpace;
            nowSum[nextLine] = nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][n][nextLine];
            if (nowSum[lineNum] > w) ngCount++;
            if (nowSum[nextLine] > w) ngCount++;
            columnSchedule.columnNum[i][n] = nextLine;
            for (int jj = 0; jj < kouhoCount; ++jj) {
              if (chaewon & (1 << jj)) {
                int j = kouhos[jj];
                columnSchedule.columnNum[i][j] = lineNum;
              }
            }

            // 逆引き更新
            for (int j = 0; j < ans.ansLineCount[i]; ++j) {
              columnSchedule.schedulesCount[i][j] = 0;
            }
            for (int j = 0; j < elementCount; ++j) {
              int lineNum = columnSchedule.columnNum[i][j];
              columnSchedule.schedules[i][lineNum][columnSchedule.schedulesCount[i][lineNum]] = j;
              columnSchedule.schedulesCount[i][lineNum]++;
            }

            break;
          }
        }
      }

      if (ngCount == 0) { break; }
    }

    if (ngCount > 0) {
      ng = 1;
      break;
    }
  }

  if (ng) {
    CopyFromBestAns();
    return 0;
  }

  // ansScheduleLineNumからansを作成
  {
    std::array<int, MAX_LINECOUNT> now = {};
    std::array<int, MAX_LINECOUNT> lastNum = {};
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < ans.ansLineCount[i]; ++j) {
        now[j] = 0;
        lastNum[j] = -1;
      }
      for (int j = elementCount - 1; j >= 0; --j) {
        int posNum = columnSchedule.columnNum[i][j];
        int need = preCalcScheduleSizes[i][j][posNum];
        ans.ans[i][j][0] = now[posNum];
        ans.ans[i][j][2] = now[posNum] + need;
        ans.ans[i][j][1] = ans.ansLinePos[i][posNum];
        ans.ans[i][j][3] = ans.ansLinePos[i][posNum + 1];
        now[posNum] += need;
        lastNum[posNum] = j;
      }
      for (int j = 0; j < ans.ansLineCount[i]; ++j) {
        if (lastNum[j] == -1) { continue; }
        ans.ans[i][lastNum[j]][2] = w;
      }
    }
    ans.ansScore = CalcScore();
    int ret = 1;
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
      ret = 2;
    }
    CopyFromBestAns();
    return ret;
  }
}

void Method3_Oshii2()
{
  for (int i = 0; i < dayCount; ++i) {
    if (oshiiLineCount[i] == -1) { return; }
  }

  // 作り直し
  {
    ans.ansBaseLineCount = oshiiLineCount[0];
    for (int i = 0; i < dayCount; ++i) {
      ans.ansLineCount[i] = oshiiLineCount[i];
      for (int j = 0; j < oshiiLineCount[i] + 1; ++j) {
        ans.ansLinePos[i][j] = oshiiLinePos[i][j];
      }
    }

    std::array<int, MAX_LINECOUNT> now = {};
    int ng = 0;
    for (int ii = dayCount - 1; ii >= 0; --ii) {
      int i = daysDifficultySorted[ii];
      for (int j = 0; j < ans.ansLineCount[i]; ++j) {
        now[j] = 0;
      }

      for (int j = elementCount - 1; j >= 0; --j) {
        int minAmari = INT_INF;
        int minOver = INT_INF;
        int posNum = -1;
        int tmpNeed = 0;
        for (int k = ans.ansLineCount[i] - 1; k >= 0; --k) {
          int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
          if (!canFitInColumn(i, j, width)) { break; }
          int need = calculateRequiredSize(elementSizes[i][j], width);
          int over = 0;
          int amari = need * width - elementSizes[i][j];
          if (now[k] + need > w) {
            over = (now[k] + need - w) * width;
            amari = 0;
          }

          if (over < minOver || (over == minOver && amari <= minAmari)) {
            minAmari = amari;
            minOver = over;
            posNum = k;
            tmpNeed = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans.ans[i][j][0] = now[posNum];
        ans.ans[i][j][2] = now[posNum] + tmpNeed;
        ans.ans[i][j][1] = ans.ansLinePos[i][posNum];
        ans.ans[i][j][3] = ans.ansLinePos[i][posNum + 1];
        now[posNum] += tmpNeed;
      }
      if (ng) { break; }
    }

    if (ng) {
      CopyFromBestAns();
      return;
    }
  }

  // 調整
  {
    int ng = 0;
    std::array<int, MAX_LINECOUNT> now = {};
    std::array<int, MAX_LINECOUNT> lastNum = {};
    for (int ii = 0; ii < dayCount; ++ii) {
      int i = daysDifficultySorted[ii];
      int ok = 1;
      for (int j = 0; j < elementCount; ++j) {
        if (ans.ans[i][j][2] > w) {
          ok = 0;
          break;
        }
      }
      if (ok) { continue; }

      ok = 0;
      for (int newjeans = 1; newjeans < ans.ansLineCount[i]; ++newjeans) {
        int keepPos = ans.ansLinePos[i][newjeans];
        ans.ansLinePos[i][newjeans] = ans.ansLinePos[i][newjeans + 1];

        for (int j = 0; j < ans.ansLineCount[i]; ++j) {
          now[j] = 0;
          lastNum[j] = -1;
        }

        int hanni = 1;
        for (int j = elementCount - 1; j >= 0; --j) {
          int minAmari = INT_INF;
          int posNum = -1;
          int tmpNeed = 0;
          for (int k = ans.ansLineCount[i] - 1; k >= 0; --k) {
            int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
            if (!canFitInColumn(i, j, width)) { continue; }
            int need = calculateRequiredSize(elementSizes[i][j], width);
            if (now[k] + need > w) { continue; }
            int amari = need * width - elementSizes[i][j];
            if (amari <= minAmari) {
              minAmari = amari;
              posNum = k;
              tmpNeed = need;
            }
          }

          if (posNum == -1) {
            hanni = 0;
            break;
          }

          ans.ans[i][j][0] = now[posNum];
          ans.ans[i][j][2] = now[posNum] + tmpNeed;
          ans.ans[i][j][1] = ans.ansLinePos[i][posNum];
          ans.ans[i][j][3] = ans.ansLinePos[i][posNum + 1];
          now[posNum] += tmpNeed;
          lastNum[posNum] = j;
        }
        if (hanni == 1) {
          for (int j = 0; j < ans.ansLineCount[i]; ++j) {
            if (lastNum[j] == -1) { continue; }
            ans.ans[i][lastNum[j]][2] = w;
          }
          ok = 1;
          break;
        }
        else {
          ans.ansLinePos[i][newjeans] = keepPos;
        }
      }
      if (ok == 0) {
        ng = 1;
        break;
      }
    }

    if (ng) {
      CopyFromBestAns();
      return;
    }
  }

  // ansScheduleLineNumからansを作成
  {
    ans.ansScore = CalcScore();
    if (ans.ansScore < best_ans.ansScore) {
      CopyToBestAns();
    }
    CopyFromBestAns();
  }
}

int oshiiDecideMethod = 0;
int oshiiMinMax = INT_INF;
int oshiiMinNGCount = INT_INF;

int Method3_Normal(int loopCount)
{
  const int MIN_LINECOUNT = 2;
  ans.ansBaseLineCount = Rand() % elementCount + MIN_LINECOUNT;
  if (best_ans.ansBaseLineCount >= MIN_LINECOUNT && Rand() % 2 == 0) {
    ans.ansBaseLineCount = best_ans.ansBaseLineCount + 1;
  }
  else if (loopCount >= 100000 && best_ans.ansBaseLineCount != -1) {
    ans.ansBaseLineCount = min(ans.ansBaseLineCount, best_ans.ansBaseLineCount + 1);
    if (Rand() % 10 != 0) {
      ans.ansBaseLineCount = max(ans.ansBaseLineCount, best_ans.ansBaseLineCount - 1);
    }
  }
  if (loopCount >= 100000 && best_ans.ansBaseLineCount == -1) {
    ans.ansBaseLineCount = Rand() % 3 + MIN_LINECOUNT;
  }

  if (ans.ansBaseLineCount <= M3_alreadyCount) {
    return 0;
  }

  if (ans.ansBaseLineCount > lineMaxLimit) {
    return 0;
  }

  int startW = ans.ansLinePos[0][M3_alreadyCount];
  int nokoriW = w - startW;
  int startLine = M3_alreadyCount;
  int nokoriLine = ans.ansBaseLineCount - M3_alreadyCount;

  {
    int ra = Rand() % 100;
    if (ra < 25) {
      int hosyouWidth = Rand() % 21;
      if (startW + nokoriLine * hosyouWidth > 800) {
        return 0;
      }
      int ok2 = 0;
      for (int _ = 0; _ < 10; ++_) {
        ans.ansLinePos[0][ans.ansBaseLineCount] = (startW + nokoriW - nokoriLine * hosyouWidth);
        for (int i = startLine + 1; i < ans.ansBaseLineCount; ++i) {
          ans.ansLinePos[0][i] = startW + Rand() % (nokoriW - nokoriLine * hosyouWidth);
        }
        sort(ans.ansLinePos[0] + startLine, ans.ansLinePos[0] + ans.ansBaseLineCount);
        int ok = 1;
        for (int i = 0; i < ans.ansBaseLineCount; ++i) {
          if (ans.ansLinePos[0][i] >= ans.ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) {
        return 0;
      }
      for (int i = startLine; i < ans.ansBaseLineCount; ++i) {
        ans.ansLinePos[0][i + 1] += hosyouWidth * (i + 1 - startLine);
      }
    }
    else if (ra < 35) {
      int ok2 = 0;
      for (int _ = 0; _ < 10; ++_) {
        ans.ansLinePos[0][ans.ansBaseLineCount] = w;
        int maxNeed = (maxElementSize[elementCount - 1] - 1) / w + 1;
        ans.ansLinePos[0][ans.ansBaseLineCount - 1] = w - maxNeed;
        if (ans.ansLinePos[0][ans.ansBaseLineCount - 1] <= ans.ansLinePos[0][startLine]) {
          break;
        }
        int ng = 0;
        for (int i = startLine + 1; i < ans.ansBaseLineCount - 1; ++i) {
          ans.ansLinePos[0][i] = startW + (w - maxNeed - startW) * i / (ans.ansBaseLineCount - 1 - startLine);
          ans.ansLinePos[0][i] += Rand() % 31 - 15;
          if (ans.ansLinePos[0][i] <= ans.ansLinePos[0][i - 1]) {
            ng = 1;
            break;
          }
        }
        if (ng) {
          continue;
        }
        sort(ans.ansLinePos[0] + startLine, ans.ansLinePos[0] + ans.ansBaseLineCount);
        int ok = 1;
        for (int i = 0; i < ans.ansBaseLineCount; ++i) {
          if (ans.ansLinePos[0][i] == ans.ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) {
        return 0;
      }
    }
    else if (ra < 50) {
      int ok2 = 0;
      for (int _ = 0; _ < 10; ++_) {
        ans.ansLinePos[0][ans.ansBaseLineCount] = w;
        int maxNeed = (maxElementSize[elementCount - 1] - 1) / w + 1;
        ans.ansLinePos[0][ans.ansBaseLineCount - 1] = w - maxNeed;
        if (ans.ansLinePos[0][ans.ansBaseLineCount - 1] <= ans.ansLinePos[0][startLine]) { break; }
        int ng = 0;
        for (int i = startLine + 1; i < ans.ansBaseLineCount - 1; ++i) {
          ans.ansLinePos[0][i] = startW + Rand() % (nokoriW - maxNeed);
        }
        sort(ans.ansLinePos[0] + startLine, ans.ansLinePos[0] + ans.ansBaseLineCount);
        int ok = 1;
        for (int i = 0; i < ans.ansBaseLineCount; ++i) {
          if (ans.ansLinePos[0][i] == ans.ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) {
        return 0;
      }
    }
    else {
      int ok2 = 0;
      for (int _ = 0; _ < 10; ++_) {
        ans.ansLinePos[0][ans.ansBaseLineCount] = w;
        for (int i = 1; i < ans.ansBaseLineCount; ++i) {
          ans.ansLinePos[0][i] = startW + Rand() % nokoriW;
        }
        sort(ans.ansLinePos[0] + startLine, ans.ansLinePos[0] + ans.ansBaseLineCount);
        int ok = 1;
        for (int i = 0; i < ans.ansBaseLineCount; ++i) {
          if (ans.ansLinePos[0][i] >= ans.ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) return 0;
    }
  }

  for (int i = 0; i < ans.ansBaseLineCount; ++i) {
    if (ans.ansLinePos[0][i] >= ans.ansLinePos[0][i + 1]) {
      return 0;
    }
  }
  if (ans.ansBaseLineCount < best_ans.ansBaseLineCount - 1) {
    return 0;
  }

  {
    std::array<int, MAX_LINECOUNT> widths = {};
    for (int i = 0; i < nokoriLine; ++i) {
      widths[i] = ans.ansLinePos[0][startLine + i + 1] - ans.ansLinePos[0][startLine + i];
    }
    sort(widths.begin(), widths.begin() + nokoriLine);
    if (widths[nokoriLine - 1] * w < maxElementSize[elementCount - 1]) {
      return 0;
    }
    for (int i = 0; i < nokoriLine; ++i) {
      ans.ansLinePos[0][startLine + i + 1] = ans.ansLinePos[0][startLine + i] + widths[i];
    }
  }
  for (int i = 0; i < dayCount; ++i) {
    ans.ansLineCount[i] = ans.ansBaseLineCount;
  }
  for (int i = 1; i < dayCount; ++i) {
    for (int j = 0; j < ans.ansBaseLineCount + 1; ++j) {
      ans.ansLinePos[i][j] = ans.ansLinePos[0][j];
    }
  }

  int ng = 0;
  std::array<int, MAX_LINECOUNT> now = {};
  std::array<int, MAX_LINECOUNT> lastNum = {};
  int tmpOshiiMax = 0;
  int tmpOshiiNGMax = 0;
  int numsOrder[MAX_N];
  for (int ii = dayCount - 1; ii >= 0; --ii) {
    int i = daysDifficultySorted[ii];
    for (int j = 0; j < ans.ansLineCount[i]; ++j) {
      now[j] = 0;
      lastNum[j] = -1;
    }

    int tmpOshiiSum = 0;
    int tmpOshiiNGCount = 0;
    for (int j = 0; j < elementCount; ++j) numsOrder[j] = j;
    if (Rand() % 5 == 0) {
      int cnt = Rand() % 30;
      for (int _ = 0; _ < cnt; ++_) {
        int n = Rand() % (elementCount - 1);
        swap(numsOrder[n], numsOrder[n + 1]);
      }
    }
    for (int jj = elementCount - 1; jj >= 0; --jj) {
      int j = numsOrder[jj];
      if (M3_alreadyUsed[i][j]) {
        continue;
      }
      int minAmari = INT_INF;
      int minOver = INT_INF;
      int posNum = -1;
      int tmpNeed = 0;
      for (int k = ans.ansLineCount[i] - 1; k >= startLine; --k) {
        int width = ans.ansLinePos[i][k + 1] - ans.ansLinePos[i][k];
        if (!canFitInColumn(i, j, width)) {
          break;
        }
        int need = (elementSizes[i][j] - 1) / width + 1;
        int over = 0;
        int amari = need * width - elementSizes[i][j];
        if (now[k] + need > w) {
          over = (now[k] + need - w) * width;
          amari = 0;
        }
        if (over < minOver || (over == minOver && amari <= minAmari)) {
          minAmari = amari;
          minOver = over;
          posNum = k;
          tmpNeed = need;
        }
      }

      if (posNum == -1) {
        ng = 1;
        break;
      }
      if (minOver > 0) {
        if (best_ans.ansBaseLineCount != -1 && ans.ansBaseLineCount <= best_ans.ansBaseLineCount) {
          ng = 1;
          break;
        }
        if (best_ans.ansBaseLineCount == -1 && ans.ansBaseLineCount < MIN_LINECOUNT) {
          ng = 1;
          break;
        }
      }

      ans.ans[i][j][0] = now[posNum];
      ans.ans[i][j][2] = now[posNum] + tmpNeed;
      ans.ans[i][j][1] = ans.ansLinePos[i][posNum];
      ans.ans[i][j][3] = ans.ansLinePos[i][posNum + 1];
      now[posNum] += tmpNeed;
      lastNum[posNum] = j;

      if (oshiiDecideMethod == 0) {
        tmpOshiiSum += minOver;
        if (tmpOshiiSum > oshiiMinMax) {
          for (int iii = ii; iii < dayCount - 1; ++iii) {
            swap(daysDifficultySorted[iii], daysDifficultySorted[iii + 1]);
          }
          ng = 1;
          break;
        }
      }
      else if (oshiiDecideMethod == 1) {
        if (minOver > 0) { tmpOshiiNGCount++; }
        if (tmpOshiiNGCount > oshiiMinNGCount) {
          ng = 1;
          break;
        }
      }
    }
    if (ng == 1) {
      break;
    }

    for (int j = startLine; j < ans.ansLineCount[i]; ++j) {
      if (lastNum[j] == -1) {
        continue;
      }
      ans.ans[i][lastNum[j]][2] = w;
    }

    tmpOshiiMax = max(tmpOshiiMax, tmpOshiiSum);
    tmpOshiiNGMax = max(tmpOshiiNGMax, tmpOshiiNGCount);
  }

  if (oshiiDecideMethod == 0) {
    if (tmpOshiiMax > 0) { ng = 2; }
  }
  else if (oshiiDecideMethod == 1) {
    if (tmpOshiiNGMax > 0) { ng = 2; }
  }

  if (ng == 0) {
    int ret = 1;
    ans.ansScore = CalcScoreForMethod3();
    if (ans.ansScore < best_ans.ansScore) {
      ret = 2;
      CopyToBestAns();
      for (int i = 0; i < dayCount; ++i) {
        oshiiLineCount[i] = -1;
      }
      oshiiMinMax = INT_INF;
      oshiiMinNGCount = INT_INF;
    }
    return ret;
  }
  else if (ng == 2) {
    if (oshiiDecideMethod == 0) {
      if (tmpOshiiMax < oshiiMinMax) {
        for (int i = 0; i < dayCount; ++i) {
          oshiiLineCount[i] = ans.ansLineCount[i];
          for (int j = 0; j < ans.ansLineCount[i] + 1; ++j) {
            oshiiLinePos[i][j] = ans.ansLinePos[i][j];
          }
        }
        oshiiMinMax = tmpOshiiMax;
      }
    }
    else if (oshiiDecideMethod == 1) {
      if (tmpOshiiNGMax < oshiiMinNGCount) {
        for (int i = 0; i < dayCount; ++i) {
          oshiiLineCount[i] = ans.ansLineCount[i];
          for (int j = 0; j < ans.ansLineCount[i] + 1; ++j) {
            oshiiLinePos[i][j] = ans.ansLinePos[i][j];
          }
        }
        oshiiMinNGCount = tmpOshiiNGMax;
      }
    }
    return 0;
  }
  return 0;
}

void Method3_1(double timeLimit)
{
  keep31Count = 0;
  int loopCount = 0;
  ans.ansBaseLineCount = -1;
  best_ans.ansBaseLineCount = -1;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < MAX_LINECOUNT; ++j) {
      ans.ansLinePos[i][j] = 0;
    }
    ans.ansLineCount[i] = -1;
    best_ans.ansLineCount[i] = -1;
    oshiiLineCount[i] = -1;
  }
  oshiiMinMax = INT_INF;
  oshiiMinNGCount = INT_INF;

  M3_alreadyCount = 0;
  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      M3_alreadyUsed[i][j] = 0;
    }
  }

  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      if (timer.get_elapsed_time() > timeLimit) {
        break;
      }
    }

    // already更新
    if (false) {
      M3_alreadyCount = 0;
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < elementCount; ++j) {
          M3_alreadyUsed[i][j] = 0;
        }
      }

      if (best_ans.ansBaseLineCount == -1) {
        M3_alreadyCount = Rand() % 3;
      }
      else {
        if (best_ans.ansBaseLineCount <= SMALL_DATA_THRESHOLD) {
          M3_alreadyCount = Rand() % 3;
        }
        else {
          M3_alreadyCount = Rand() % (elementCount / 2 + 1);
        }
      }

      if (M3_alreadyCount == 0) {
        continue;
      }
      for (int blackpink = 0; blackpink < M3_alreadyCount; ++blackpink) {
        int baseWidth = Rand() % 100 + 50;
        if (M3_alreadyCount >= 5) { baseWidth = Rand() % 100 + 20; }
        if (blackpink > 0 && ans.ansLinePos[0][blackpink - 1] + baseWidth + 20 > w) {
          M3_alreadyCount = 0;
          for (int i = 0; i < dayCount; ++i) {
            for (int j = 0; j < elementCount; ++j) {
              M3_alreadyUsed[i][j] = 0;
            }
          }
          break;
        }

        int maxWidth = -1;
        double maxValue = 0;
        for (int lisa = 0; lisa < 20; ++lisa) {
          int width = baseWidth + lisa;
          int minValue = INT_INF;
          int ng = 0;
          for (int i = 0; i < dayCount; ++i) {
            int ok = 0;
            for (int k = 0; k < bestPairCount[i][width]; ++k) {
              if (M3_alreadyUsed[i][bestPairs[i][width][k][0]] == 0 && M3_alreadyUsed[i][bestPairs[i][width][k][1]] == 0) {
                ok = 1;
                if (bestPairValue[i][width][k] < minValue) { minValue = bestPairValue[i][width][k]; }
                break;
              }
            }
            if (ok == 0) {
              ng = 1;
              break;
            }
          }
          if (ng) {
            continue;
          }
          double tmpValue = (double)minValue / (width * w);
          if (tmpValue > maxValue) {
            maxValue = tmpValue;
            maxWidth = width;
          }
        }
        if (maxWidth == -1) {
          M3_alreadyCount = 0;
          for (int i = 0; i < dayCount; ++i) {
            for (int j = 0; j < elementCount; ++j) {
              M3_alreadyUsed[i][j] = 0;
            }
          }
          break;
        }

        ans.ansLinePos[0][blackpink + 1] = ans.ansLinePos[0][blackpink] + maxWidth;
        for (int i = 0; i < dayCount; ++i) {
          for (int k = 0; k < bestPairCount[i][maxWidth]; ++k) {
            if (M3_alreadyUsed[i][bestPairs[i][maxWidth][k][0]] == 0 && M3_alreadyUsed[i][bestPairs[i][maxWidth][k][1]] == 0) {
              int j1 = bestPairs[i][maxWidth][k][0];
              int j2 = bestPairs[i][maxWidth][k][1];
              M3_alreadyUsed[i][j1] = 1;
              M3_alreadyUsed[i][j2] = 1;

              int need1 = (elementSizes[i][j1] - 1) / maxWidth + 1;
              ans.ans[i][j1][0] = 0;
              ans.ans[i][j1][2] = need1;
              ans.ans[i][j1][1] = ans.ansLinePos[0][blackpink];
              ans.ans[i][j1][3] = ans.ansLinePos[0][blackpink + 1];
              ans.ans[i][j2][0] = need1;
              ans.ans[i][j2][2] = w;
              ans.ans[i][j2][1] = ans.ansLinePos[0][blackpink];
              ans.ans[i][j2][3] = ans.ansLinePos[0][blackpink + 1];
              break;
            }
          }
        }
      }

      continue;
    }

    int ra = Rand() % 10000;
    int isOK = 0;
    if (ra < 9990) {
      isOK = Method3_Normal(loopCount);
    }
    else if (ra < 10000) {
      isOK = Method3_Oshii();
      for (int i = 0; i < dayCount; ++i) {
        oshiiLineCount[i] = -1;
      }
      oshiiMinMax = INT_INF;
      oshiiMinNGCount = INT_INF;
    }

    if (isOK == 2) {
      CopyToKeep31(keep31Count % keep31KeepSize);
      keep31Count++;
    }
  }

  CopyFromBestAns();
}

// 縦線をずらす
void Method3_2(double timeLimit)
{
  // 適用可能かチェック
  {
    for (int i = 0; i < dayCount - 1; ++i) {
      for (int j = 0; j < ans.ansLineCount[i]; ++j) {
        if (ans.ansLinePos[i][j] != ans.ansLinePos[i + 1][j]) {
          return;
        }
      }
    }
  }

  int loopCount = 0;
  double start_temp = 100.1;
  double end_temp = 0.0;

  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      if (timer.get_elapsed_time() > timeLimit) {
        break;
      }
    }

    int ra = Rand() % (ans.ansBaseLineCount - 1) + 1;
    int ra2 = 0;
    while (ra2 == 0) {
      ra2 = Rand() % 11 - 5;
    }
    ans.ansLinePos[0][ra] += ra2;
    if (ans.ansLinePos[0][ra - 1] >= ans.ansLinePos[0][ra] || ans.ansLinePos[0][ra] >= ans.ansLinePos[0][ra + 1]) {
      ans.ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ra >= 2 && ans.ansLinePos[0][ra - 1] - ans.ansLinePos[0][ra - 2] > ans.ansLinePos[0][ra] - ans.ansLinePos[0][ra - 1]) {
      ans.ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ans.ansLinePos[0][ra] - ans.ansLinePos[0][ra - 1] > ans.ansLinePos[0][ra + 1] - ans.ansLinePos[0][ra]) {
      ans.ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ra <= ans.ansBaseLineCount - 2 && ans.ansLinePos[0][ra + 1] - ans.ansLinePos[0][ra] > ans.ansLinePos[0][ra + 2] - ans.ansLinePos[0][ra + 1]) {
      ans.ansLinePos[0][ra] -= ra2;
      continue;
    }

    int ng = 0;
    std::array<int, MAX_LINECOUNT> now = {};
    std::array<int, MAX_LINECOUNT> lastNum = {};
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < ans.ansBaseLineCount; ++j) {
        now[j] = 0;
        lastNum[j] = -1;
      }

      for (int j = elementCount - 1; j >= 0; --j) {
        int minAmari = INT_INF;
        int posNum = -1;
        int tmpNeed = 0;
        for (int k = ans.ansBaseLineCount - 1; k >= 0; --k) {
          int width = ans.ansLinePos[0][k + 1] - ans.ansLinePos[0][k];
          if (!canFitInColumn(i, j, width)) {
            break;
          }
          int need = calculateRequiredSize(elementSizes[i][j], width);
          if (now[k] + need > w) {
            continue;
          }
          int amari = need * width - elementSizes[i][j];
          if (amari <= minAmari) {
            minAmari = amari;
            posNum = k;
            tmpNeed = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans.ans[i][j][0] = now[posNum];
        ans.ans[i][j][2] = now[posNum] + tmpNeed;
        ans.ans[i][j][1] = ans.ansLinePos[0][posNum];
        ans.ans[i][j][3] = ans.ansLinePos[0][posNum + 1];
        now[posNum] += tmpNeed;
        lastNum[posNum] = j;
      }
      if (ng) {
        break;
      }

      for (int j = 0; j < ans.ansBaseLineCount; ++j) {
        if (lastNum[j] == -1) {
          continue;
        }
        ans.ans[i][lastNum[j]][2] = w;
      }
    }

    if (ng) {
      ans.ansLinePos[0][ra] -= ra2;
      continue;
    }

    int beforeScore = ans.ansScore;
    ans.ansScore = CalcScoreForMethod3();
    int diffScore = beforeScore - ans.ansScore;

    double temp = (start_temp + (end_temp - start_temp) * timer.get_elapsed_time() / timeLimit);
    const double prob = exp((double)diffScore / temp);

    if (prob > Rand01()) {
      for (int i = 1; i < dayCount; ++i) {
        ans.ansLinePos[i][ra] = ans.ansLinePos[0][ra];
      }
      if (ans.ansScore <= best_ans.ansScore) { CopyToBestAns(); }
    }
    else {
      ans.ansScore = beforeScore;
      ans.ansLinePos[0][ra] -= ra2;
    }
  }

  CopyFromBestAns();
}

void Method6_ColumnShuffle(double timeLimit)
{
  std::array<int, MAX_LINECOUNT> widths = {};
  for (int j = 0; j < ans.ansBaseLineCount; ++j) {
    widths[j] = ans.ansLinePos[0][j + 1] - ans.ansLinePos[0][j];
  }

  for (int i = 0; i < dayCount; ++i) {
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < ans.ansBaseLineCount; ++k) {
        if (ans.ans[i][j][1] == ans.ansLinePos[i][k] && ans.ans[i][j][3] == ans.ansLinePos[i][k + 1]) {
          columnSchedule.columnNum[i][j] = k;
          break;
        }
      }
    }
  }

  int v[MAX_LINECOUNT] = {};
  for (int i = 0; i < ans.ansBaseLineCount; ++i) {
    v[i] = i;
  }
  for (int ningning = 0; ningning < 100; ++ningning) {
    FisherYates(v, ans.ansBaseLineCount);
    int nextLinePos[MAX_LINECOUNT] = {};
    int nextArgPos[MAX_LINECOUNT] = {};
    for (int i = 0; i < ans.ansBaseLineCount; ++i) {
      nextLinePos[i + 1] = nextLinePos[i] + widths[v[i]];
      nextArgPos[v[i]] = i;
    }
    for (int i = 0; i < dayCount; ++i) {
      for (int j = 0; j < elementCount; ++j) {
        int nextLine = nextArgPos[columnSchedule.columnNum[i][j]];
        ans.ans[i][j][1] = nextLinePos[nextLine];
        ans.ans[i][j][3] = nextLinePos[nextLine + 1];
      }
    }
    ans.ansScore = CalcScore();
    if (ans.ansScore < best_ans.ansScore) { CopyToBestAns(); }
    if (ningning % 10 == 9 && timer.get_elapsed_time() > timeLimit) {
      break;
    }
  }
  CopyFromBestAns();
}

void Method7()
{
  CopyFromBestAns();
  CoptToCurrent_M42();

  for (int i = 0; i < ans.ansBaseLineCount; ++i) {
    widths[i] = ans.ansLinePos[0][i + 1] - ans.ansLinePos[0][i];
  }

  int yokoLineCount[MAX_LINECOUNT] = {};
  for (int j = 0; j < ans.ansBaseLineCount; ++j) {
    for (int i = 1; i < dayCount; ++i) {
      yokoLineCount[j] += max(0, columnSchedule.schedulesCount[i - 1][j] - 1) + max(0, columnSchedule.schedulesCount[i][j] - 1);
      int ite1 = 1;
      int ite2 = 1;
      while (ite1 < columnSchedule.schedulesCount[i - 1][j] && ite2 < columnSchedule.schedulesCount[i][j]) {
        if (columnSchedule.schedulesPosition[i - 1][j][ite1] == columnSchedule.schedulesPosition[i][j][ite2]) {
          yokoLineCount[j] -= 2;
          ite1++;
          ite2++;
        }
        else if (columnSchedule.schedulesPosition[i - 1][j][ite1] < columnSchedule.schedulesPosition[i][j][ite2]) {
          ite1++;
        }
        else {
          ite2++;
        }
      }
    }
  }

  for (int line1 = 0; line1 < ans.ansBaseLineCount; ++line1) {
    for (int line2 = 0; line2 < ans.ansBaseLineCount; ++line2) {
      if (line1 == line2) {
        continue;
      }
      if (yokoLineCount[line1] <= yokoLineCount[line2]) {
        continue;
      }
      int margin = widths[line1] - 1;
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < columnSchedule.schedulesCount[i][line1]; ++j) {
          int num = columnSchedule.schedules[i][line1][j];
          int height = columnSchedule.schedulesPosition[i][line1][j + 1] - columnSchedule.schedulesPosition[i][line1][j];
          for (int k = 0; k < elementCount; ++k) {
            if (k == num) {
              break;
            }
            if (columnSchedule.columnNum[i][k] == line1) {
              continue;
            }
            if (preCalcScheduleSizes[i][k][line1] <= height) {
              int nextLine = columnSchedule.columnNum[i][k];
              int nextLineIndex = -1;
              for (int l = 0; l < columnSchedule.schedulesCount[i][nextLine]; ++l) {
                if (columnSchedule.schedules[i][nextLine][l] == k) {
                  nextLineIndex = l;
                  break;
                }
              }
              if (preCalcScheduleSizes[i][num][nextLine] > columnSchedule.schedulesPosition[i][nextLine][nextLineIndex + 1] - columnSchedule.schedulesPosition[i][nextLine][nextLineIndex]) {
                continue;
              }
              swap(columnSchedule.schedules[i][line1][j], columnSchedule.schedules[i][nextLine][nextLineIndex]);
              swap(columnSchedule.columnNum[i][num], columnSchedule.columnNum[i][k]);
              break;
            }
          }
          int newNum = columnSchedule.schedules[i][line1][j];
          // margin計算
          margin = min(margin, calculateMargin(i, newNum, widths[line1], height));
        }
      }

      if (margin == 0) {
        continue;
      }
      int diffScore = margin * (yokoLineCount[line1] - yokoLineCount[line2]);

      // ansLinePosとpreCalcScheduleSizesとwidthsを更新
      widths[line1] -= margin;
      widths[line2] += margin;
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < ans.ansBaseLineCount; ++j) {
          ans.ansLinePos[i][j + 1] = ans.ansLinePos[i][j] + widths[j];
        }
      }
      for (int i = 0; i < dayCount; ++i) {
        for (int j = 0; j < elementCount; ++j) {
          preCalcScheduleSizes[i][j][line1] = (elementSizes[i][j] - 1) / widths[line1] + 1;
          preCalcScheduleSizes[i][j][line2] = (elementSizes[i][j] - 1) / widths[line2] + 1;
        }
      }

      ans.ansScore -= diffScore;
    }
  }

  if (ans.ansScore < best_ans.ansScore) {
    CopyToBestAns();
    CopyToBest_M42();
  }
}

int M8ans[1010][MAX_D][MAX_N][4];
int M8ansOK[1010][MAX_D];
int M8ansScore[1010][MAX_D];
int M8ansLineCount[1010];
int M8ansLinePos[1010][MAX_LINECOUNT];
int M8ansansCount;
int M8ansNum[MAX_D];
int M8ansansCountEachDay[MAX_D];
int M8ansansNumEachDay[MAX_D][1010];

void CopyM8ToAns()
{
  for (int i = 0; i < dayCount; ++i) {
    int num = M8ansNum[i];
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < 4; ++k) {
        ans.ans[i][j][k] = M8ans[num][i][j][k];
      }
    }
  }
}

void Method8(double timeLimit)
{
  const int TMP_NUM = 1000;
  const int NG_NUM = 1005;
  M8ansansCount = 0;
  for (int i = 0; i < dayCount; ++i) {
    M8ansNum[i] = -1;
  }

  for (int i = 0; i < dayCount; ++i) {
    M8ansNum[i] = NG_NUM;
    M8ansScore[NG_NUM][i] = INVALID_SIZE;
    for (int j = 0; j < elementCount; ++j) {
      for (int k = 0; k < 4; ++k) {
        M8ans[NG_NUM][i][j][k] = ans.ans[i][j][k];
      }
    }
  }

  const double M8_startTime = timer.get_elapsed_time();
  double M8_NowTime = timer.get_elapsed_time();
  int loopCount = 0;
  while (M8ansansCount < 1000) {
    loopCount++;
    if (loopCount % 25 == 0) {
      M8_NowTime = timer.get_elapsed_time();
      if (M8_NowTime > M8_startTime + (timeLimit - M8_startTime) / 2) {
        break;
      }
    }

    // 縦線作成
    ans.ansBaseLineCount = Rand() % 5 + 2;
    int ok2 = 0;
    for (int _2 = 0; _2 < 10; ++_2) {
      ans.ansLinePos[0][ans.ansBaseLineCount] = w;
      for (int i = 1; i < ans.ansBaseLineCount; ++i) {
        ans.ansLinePos[0][i] = Rand() % w;
      }
      sort(ans.ansLinePos[0], ans.ansLinePos[0] + ans.ansBaseLineCount);
      int ok = 1;
      for (int i = 0; i < ans.ansBaseLineCount; ++i) {
        if (ans.ansLinePos[0][i] >= ans.ansLinePos[0][i + 1]) {
          ok = 0;
          break;
        }
      }
      if (ok) {
        ok2 = 1;
        break;
      }
    }
    if (ok2 == 0) {
      continue;
    }

    {
      std::array<int, MAX_LINECOUNT> widths = {};
      for (int i = 0; i < ans.ansBaseLineCount; ++i) {
        widths[i] = ans.ansLinePos[0][i + 1] - ans.ansLinePos[0][i];
      }
      sort(widths.begin(), widths.begin() + ans.ansBaseLineCount);
      for (int i = 0; i < ans.ansBaseLineCount; ++i) {
        ans.ansLinePos[0][i + 1] = ans.ansLinePos[0][i] + widths[i];
      }
    }

    M8ansLineCount[TMP_NUM] = ans.ansBaseLineCount;
    for (int i = 0; i < ans.ansBaseLineCount + 1; ++i) {
      M8ansLinePos[TMP_NUM][i] = ans.ansLinePos[0][i];
    }

    int okCount = 0;
    std::array<int, MAX_LINECOUNT> now = {};
    std::array<int, MAX_LINECOUNT> lastNum = {};
    int numsOrder[MAX_N];
    for (int ii = dayCount - 1; ii >= 0; --ii) {
      int i = daysDifficultySorted[ii];
      for (int j = 0; j < ans.ansBaseLineCount; ++j) {
        now[j] = 0;
        lastNum[j] = -1;
      }

      for (int j = 0; j < elementCount; ++j) numsOrder[j] = j;
      if (Rand() % 5 == 0) {
        int cnt = Rand() % 30;
        for (int _ = 0; _ < cnt; ++_) {
          int n = Rand() % (elementCount - 1);
          swap(numsOrder[n], numsOrder[n + 1]);
        }
      }
      int ok = 1;
      for (int jj = elementCount - 1; jj >= 0; --jj) {
        int j = numsOrder[jj];
        int minAmari = INT_INF;
        int posNum = -1;
        int tmpNeed = 0;
        for (int k = ans.ansBaseLineCount - 1; k >= 0; --k) {
          int width = ans.ansLinePos[0][k + 1] - ans.ansLinePos[0][k];
          if (!canFitInColumn(i, j, width)) { break; }
          int need = calculateRequiredSize(elementSizes[i][j], width);
          int amari = need * width - elementSizes[i][j];
          if (now[k] + need > w) {
            continue;
          }
          if (amari < minAmari) {
            minAmari = amari;
            posNum = k;
            tmpNeed = need;
          }
        }

        if (posNum == -1) {
          ok = 0;
          break;
        }

        M8ans[TMP_NUM][i][j][0] = now[posNum];
        M8ans[TMP_NUM][i][j][2] = now[posNum] + tmpNeed;
        M8ans[TMP_NUM][i][j][1] = ans.ansLinePos[0][posNum];
        M8ans[TMP_NUM][i][j][3] = ans.ansLinePos[0][posNum + 1];
        now[posNum] += tmpNeed;
        lastNum[posNum] = j;
      }
      if (ok == 1) {
        okCount++;
        M8ansOK[TMP_NUM][i] = 1;
        for (int j = 0; j < ans.ansBaseLineCount; ++j) {
          if (lastNum[j] == -1) { continue; }
          M8ans[TMP_NUM][i][lastNum[j]][2] = w;
        }
      }
      else {
        M8ansOK[TMP_NUM][i] = 0;
      }
    }

    if (okCount >= 2) {
      M8ansLineCount[M8ansansCount] = ans.ansBaseLineCount;
      for (int i = 0; i < ans.ansBaseLineCount + 1; ++i) {
        M8ansLinePos[M8ansansCount][i] = M8ansLinePos[TMP_NUM][i];
      }
      for (int i = 0; i < dayCount; ++i) {
        M8ansOK[M8ansansCount][i] = M8ansOK[TMP_NUM][i];
        if (M8ansOK[M8ansansCount][i]) {
          int score = 0;
          for (int j = 0; j < elementCount; ++j) {
            for (int k = 0; k < 4; ++k) {
              M8ans[M8ansansCount][i][j][k] = M8ans[TMP_NUM][i][j][k];
            }
            if (M8ans[M8ansansCount][i][j][0] != 0 && M8ans[M8ansansCount][i][j][0] != w) { score += M8ans[M8ansansCount][i][j][3] - M8ans[M8ansansCount][i][j][1]; }
            if (M8ans[M8ansansCount][i][j][2] != 0 && M8ans[M8ansansCount][i][j][2] != w) { score += M8ans[M8ansansCount][i][j][3] - M8ans[M8ansansCount][i][j][1]; }
          }
          score /= 2;
          score += (M8ansLineCount[M8ansansCount] - 1) * w;
          M8ansScore[M8ansansCount][i] = score;
          if (score < M8ansScore[M8ansNum[i]][i]) { M8ansNum[i] = M8ansansCount; }
        }
      }
      M8ansansCount++;
    }
  }

  CopyM8ToAns();
  ans.ansScore = CalcScore();
  if (ans.ansScore < best_ans.ansScore) { CopyToBestAns(); }
  for (int i = 0; i < dayCount; ++i) {
    M8ansansCountEachDay[i] = 0;
    for (int j = 0; j < M8ansansCount; ++j) {
      if (M8ansOK[j][i]) {
        M8ansansNumEachDay[i][M8ansansCountEachDay[i]] = j;
        M8ansansCountEachDay[i]++;
      }
    }
  }

  double M8_startTime2 = timer.get_elapsed_time();
  loopCount = 0;
  M8_NowTime = timer.get_elapsed_time();
  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      M8_NowTime = timer.get_elapsed_time();
      if (M8_NowTime > timeLimit) {
        break;
      }
    }

    if (loopCount % 12345 == 0) {
      CopyM8ToAns();
      ans.ansScore = CalcScore();
      if (ans.ansScore < best_ans.ansScore) { CopyToBestAns(); }
    }

    int d = Rand() % dayCount;
    if (M8ansansCountEachDay[d] == 0) { continue; }
    int num = M8ansansNumEachDay[d][Rand() % M8ansansCountEachDay[d]];
    if (num == M8ansNum[d]) { continue; }

    int beforeNum = M8ansNum[d];
    int diffScore = (M8ansScore[beforeNum][d] - M8ansScore[num][d]) * 2;
    if (d == 0 || d == dayCount - 1) { diffScore = (M8ansScore[beforeNum][d] - M8ansScore[num][d]); }

    if (d > 0) {
      int beforeDayNum = M8ansNum[d - 1];
      if (beforeNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[beforeNum] && ite2 < M8ansLineCount[beforeDayNum]) {
          if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            diffScore -= w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
      if (num < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[num] && ite2 < M8ansLineCount[beforeDayNum]) {
          if (M8ansLinePos[num][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            diffScore += w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[num][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
    }
    if (d < dayCount - 1) {
      int afterDayNum = M8ansNum[d + 1];
      if (beforeNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[beforeNum] && ite2 < M8ansLineCount[afterDayNum]) {
          if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            diffScore -= w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
      if (num < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[num] && ite2 < M8ansLineCount[afterDayNum]) {
          if (M8ansLinePos[num][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            diffScore += w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[num][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
    }

    double timeRatio = (M8_NowTime - M8_startTime2) / (timeLimit - M8_startTime2);
    double temp = (swap_start_temp + (swap_end_temp - swap_start_temp) * timeRatio);
    const double prob = exp((double)diffScore * ANNEALING_TEMP_COEFFICIENT / temp);

    if (prob > Rand01()) { M8ansNum[d] = num; }
  }

  CopyM8ToAns();
  ans.ansScore = CalcScore();
  if (ans.ansScore < best_ans.ansScore) {
    CopyToBestAns();
  }
}

int isFind = 0;
void Method4(int setCount)
{
  CopyToBackupAns();

  CopyFromBestAns();
  CopyToTemp();

  double outerTL = TL;
  double m4StartTime = timer.get_elapsed_time();
  for (int twice = 0; twice < setCount; ++twice) {
    CopyFromTemp();
    CopyToBestAns();
    double nowTime = timer.get_elapsed_time();
    double innerTL = (outerTL - m4StartTime) * (twice + 1) / setCount + m4StartTime;
    Method3_1(nowTime + (innerTL - nowTime) * 0.5);

    CopyFromBestAns();

    if (ans.ansBaseLineCount == -1) {
      isFind = 0;
    }
    else {
      isFind = 1;
    }

    if (ans.ansScore == 1) {
      return;
    }

    if (ans.ansBaseLineCount == -1) {
      Method8(nowTime + (innerTL - nowTime) * 0.98);
      if (best_ans.ansScore < backup_ans.ansScore) {
        CopyToBackupAns();
      }
    }
    else {
      Method3_2(nowTime + (innerTL - nowTime) * 0.6);

      CopyFromBestAns();
      Method4_3(nowTime + (innerTL - nowTime) * 0.98);

      // 列位置シャッフル
      CopyFromBestAns();
      Method6_ColumnShuffle(nowTime + (innerTL - nowTime) * 1.0);

      if (best_ans.ansScore < backup_ans.ansScore) {
        CopyToBackupAns();
      }
    }
  }

  CopyFromBackupAns();
  CopyFromBestAns();
}

void Method5()
{
  double TL31 = TL * 0.5;
  Method3_1(TL31);

  if (ans.ansBaseLineCount == -1) {
    return;
  }

  CopyToBackupAns();

  keep31Count = min(keep31KeepSize, keep31Count);
  double TL2 = TL * 0.7;
  for (int tzuyu = 0; tzuyu < keep31Count; ++tzuyu) {
    double nowTime = timer.get_elapsed_time();
    double innerTL = (TL2 - nowTime) * (tzuyu + 1) / keep31KeepSize + nowTime;
    double innerTL32 = (innerTL - nowTime) / 5.0 + nowTime;
    CopyFromKeep31(tzuyu);
    CopyToBestAns();
    Method3_2(innerTL32);
    CopyFromBestAns();
    Method4_3(innerTL);
    if (best_ans.ansScore < backup_ans.ansScore) {
      CopyToBackupAns();
    }
  }

  CopyFromBackupAns();

  CopyFromBestAns();
  Method4_3(TL);
  if (best_ans.ansScore < backup_ans.ansScore) {
    CopyToBackupAns();
  }

  CopyFromBackupAns();
  CopyFromBestAns();
}

ll Solve(int case_num)
{
  timer.start();
  timer.get_elapsed_time();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(case_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(case_num, ofs);

  // 初期解生成
  Initialize();

  Method1();

  MethodPerfect();

  Method4(1);

  // 解答を出力
  CopyFromBestAns();
  Output(ofs);

  if (ofs.is_open()) { ofs.close(); }

  if (mode != 0) {
    bool isInvalidSolution = IsInvalidSolution();
  }

  ll score = 0;
  if (mode != 0) { score = CalcScore(); }
  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    vector<ll> scores;
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(i);
      sum += score;
      int maxASum = 0;
      for (int j = 0; j < elementCount; ++j) {
        maxASum += maxElementSize[j];
      }
      if (true || (isFind == 0 && ans.ansScore > 1)) {
        cout << i << ", "
          << dayCount << ", "
          << elementCount << ", "
          << 1 - emptySpaceRatio << ", "
          << maxASum << ", "
          << isFind << ", "
          << ans.ansBaseLineCount << ", "
          << score << ", "
          << sum << ' '
          << timer.get_elapsed_time() << endl;
      }
      scores.push_back(score);
    }
    for (auto score : scores) {
      cout << score << endl;
    }
  }

  return 0;
}
