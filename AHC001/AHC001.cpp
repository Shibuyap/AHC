#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

using namespace std;
using namespace chrono;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 205

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

static uint32_t xorshift()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

// 0以上1未満の小数をとる乱数
static double rand01()
{
  return (xorshift() + 0.5) * (1.0 / UINT_MAX);
}

// ハイパラはここにおく
int hyperParam1;
int beamWidth = 1;
int innerLoops = 1;
array<int, 20> modeCount;

int shuffles[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

double timeLimit = 4.8;

struct Point
{
  int x;
  int y;
};

struct Rect
{
  Point topLeft;
  Point bottomRight;
};

class State
{
public:
  array<Rect, MAX_N> rects;
  int score = -1;
};

enum Direction { HORIZONTAL = 0, VERTICAL = 1 };

int allLoopTimes = 1;
int numRects;
array<Point, MAX_N> points;
array<int, MAX_N> targetSizes;
State currentState;
array<int, MAX_N> rectAreas;
State bestState;

inline void calcArea(int idx)
{
  rectAreas[idx] = (currentState.rects[idx].bottomRight.x - currentState.rects[idx].topLeft.x) * (currentState.rects[idx].bottomRight.y - currentState.rects[idx].topLeft.y);
}

inline void readInput(int case_num)
{
  // 入力
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());
  if (!ifs.is_open()) { // 標準入力する
    cin >> numRects;
    for (int i = 0; i < numRects; ++i) cin >> points[i].x >> points[i].y >> targetSizes[i];
  }
  else { // ファイル入力する
    ifs >> numRects;
    for (int i = 0; i < numRects; ++i) ifs >> points[i].x >> points[i].y >> targetSizes[i];
  }
}

// ファイル出力
void writeOutput(int case_num)
{
  std::ostringstream oss;
  oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ofstream ofs(oss.str());

  if (ofs.is_open()) {
    for (int i = 0; i < numRects; ++i) {
      ofs << currentState.rects[i].topLeft.x << ' ' << currentState.rects[i].topLeft.y << ' ' << currentState.rects[i].bottomRight.x << ' ' << currentState.rects[i].bottomRight.y << endl;
    }

    ofs.close();
  }
}

inline void writeErrorLog(int fileNum)
{
  string fileName = to_string(fileNum);
  fileName += "_out_ERROR.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);
  for (int i = 0; i < numRects; ++i) {
    ofs << currentState.rects[i].topLeft.x << ' ' << currentState.rects[i].topLeft.y << ' ' << currentState.rects[i].bottomRight.x << ' ' << currentState.rects[i].bottomRight.y << endl;
  }
  ofs.close();
}

array<double, MAX_N> rectScores;
double totalScore;
inline int calcScore(int ite)
{
  if (ite == -1) {
    double scoreSum = 0;
    for (int i = 0; i < numRects; ++i) {
      calcArea(i);
      rectScores[i] = 1.0 - (1.0 - (double)min(targetSizes[i], rectAreas[i]) / (double)max(targetSizes[i], rectAreas[i])) * (1.0 - (double)min(targetSizes[i], rectAreas[i]) / (double)max(targetSizes[i], rectAreas[i]));
      scoreSum += rectScores[i];
    }
    totalScore = scoreSum;
    scoreSum /= numRects;
    scoreSum *= 1e9;
    return round(scoreSum);
  }
  else {
    double scoreSum = totalScore;
    scoreSum -= rectScores[ite];
    calcArea(ite);
    rectScores[ite] = 1.0 - (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite])) * (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite]));
    scoreSum += rectScores[ite];
    totalScore = scoreSum;
    scoreSum /= numRects;
    scoreSum *= 1e9;
    return round(scoreSum);
  }
}

inline int checkOverlap(int i, int j)
{
  int overlapCount = 0;
  if (currentState.rects[i].topLeft.x <= currentState.rects[j].topLeft.x && currentState.rects[j].topLeft.x < currentState.rects[i].bottomRight.x) overlapCount++;
  else if (currentState.rects[j].topLeft.x <= currentState.rects[i].topLeft.x && currentState.rects[i].topLeft.x < currentState.rects[j].bottomRight.x) overlapCount++;
  if (currentState.rects[i].topLeft.y <= currentState.rects[j].topLeft.y && currentState.rects[j].topLeft.y < currentState.rects[i].bottomRight.y) overlapCount++;
  else if (currentState.rects[j].topLeft.y <= currentState.rects[i].topLeft.y && currentState.rects[i].topLeft.y < currentState.rects[j].bottomRight.y) overlapCount++;
  return overlapCount == 2;
}

// 矩形のY軸方向の重なりをチェック
inline int checkYOverlap(const Rect& rect1, const Rect& rect2)
{
  if (rect1.topLeft.y <= rect2.topLeft.y && rect2.topLeft.y < rect1.bottomRight.y) return 1;
  if (rect2.topLeft.y <= rect1.topLeft.y && rect1.topLeft.y < rect2.bottomRight.y) return 1;
  return 0;
}

// 矩形のX軸方向の重なりをチェック
inline int checkXOverlap(const Rect& rect1, const Rect& rect2)
{
  if (rect1.topLeft.x <= rect2.topLeft.x && rect2.topLeft.x < rect1.bottomRight.x) return 1;
  if (rect2.topLeft.x <= rect1.topLeft.x && rect1.topLeft.x < rect2.bottomRight.x) return 1;
  return 0;
}

// 点とY範囲の重なりをチェック
inline int checkPointYOverlap(int y, int topY, int bottomY)
{
  if (y <= topY && topY < y + 1) return 1;
  if (topY <= y && y < bottomY) return 1;
  return 0;
}

// 点とX範囲の重なりをチェック
inline int checkPointXOverlap(int x, int leftX, int rightX)
{
  if (x <= leftX && leftX < x + 1) return 1;
  if (leftX <= x && x < rightX) return 1;
  return 0;
}

array<int, MAX_N> sortedByX, sortedByY;
array<int, MAX_N> indexInSortedX, indexInSortedY;
inline void initSortArrays()
{
  vector<P> sortPairs;
  for (int i = 0; i < numRects; ++i) {
    sortPairs.emplace_back(P(points[i].x, i));
  }
  sort(sortPairs.begin(), sortPairs.end());
  for (int i = 0; i < numRects; ++i) {
    sortedByX[i] = sortPairs[i].second;
    indexInSortedX[sortPairs[i].second] = i;
  }

  sortPairs.clear();
  for (int i = 0; i < numRects; ++i) {
    sortPairs.emplace_back(P(points[i].y, i));
  }
  sort(sortPairs.begin(), sortPairs.end());
  for (int i = 0; i < numRects; ++i) {
    sortedByY[i] = sortPairs[i].second;
    indexInSortedY[sortPairs[i].second] = i;
  }
}

// 座標が範囲内かチェックする共通関数
inline int isInRange(int coord)
{
  return 0 <= coord && coord <= 10000;
}

// 矩形の座標が範囲内かチェック
inline int isRectInRange(const Rect& rect)
{
  return isInRange(rect.topLeft.x) && isInRange(rect.topLeft.y) &&
    isInRange(rect.bottomRight.x) && isInRange(rect.bottomRight.y);
}

// 0~10000を出ていないか
// 面積は1以上か
// (x[i]+0.5,y[i]+0.5)を含んでいるか
// 重なりがないか
inline int isValid2(int ite)
{
  if (ite == -1) {
    for (int i = 0; i < numRects; ++i) {
      if (!isRectInRange(currentState.rects[i])) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (currentState.rects[i].bottomRight.x <= currentState.rects[i].topLeft.x) return 0;
      if (currentState.rects[i].bottomRight.y <= currentState.rects[i].topLeft.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (points[i].x < currentState.rects[i].topLeft.x || currentState.rects[i].bottomRight.x <= points[i].x) return 0;
      if (points[i].y < currentState.rects[i].topLeft.y || currentState.rects[i].bottomRight.y <= points[i].y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      for (int j = i + 1; j < numRects; ++j) {
        if (checkOverlap(i, j)) return 0;
      }
    }
  }
  else {
    if (!isRectInRange(currentState.rects[ite])) return 0;
    if (currentState.rects[ite].bottomRight.x <= currentState.rects[ite].topLeft.x) return 0;
    if (currentState.rects[ite].bottomRight.y <= currentState.rects[ite].topLeft.y) return 0;
    if (points[ite].x < currentState.rects[ite].topLeft.x || currentState.rects[ite].bottomRight.x <= points[ite].x) return 0;
    if (points[ite].y < currentState.rects[ite].topLeft.y || currentState.rects[ite].bottomRight.y <= points[ite].y) return 0;
    for (int i = 0; i < numRects; ++i) {
      if (i == ite) { continue; }
      if (checkOverlap(i, ite)) return 0;
    }
  }
  return 1;
}

inline int isValid(int ite)
{
  if (ite == -1) {
    for (int i = 0; i < numRects; ++i) {
      if (!isRectInRange(currentState.rects[i])) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (currentState.rects[i].bottomRight.x <= currentState.rects[i].topLeft.x) return 0;
      if (currentState.rects[i].bottomRight.y <= currentState.rects[i].topLeft.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (points[i].x < currentState.rects[i].topLeft.x || currentState.rects[i].bottomRight.x <= points[i].x) return 0;
      if (points[i].y < currentState.rects[i].topLeft.y || currentState.rects[i].bottomRight.y <= points[i].y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      for (int j = i + 1; j < numRects; ++j) {
        if (checkOverlap(i, j)) return 0;
      }
    }
  }
  else {
    if (!isRectInRange(currentState.rects[ite])) return 0;
    if (currentState.rects[ite].bottomRight.x <= currentState.rects[ite].topLeft.x) return 0;
    if (currentState.rects[ite].bottomRight.y <= currentState.rects[ite].topLeft.y) return 0;
    if (points[ite].x < currentState.rects[ite].topLeft.x || currentState.rects[ite].bottomRight.x <= points[ite].x) return 0;
    if (points[ite].y < currentState.rects[ite].topLeft.y || currentState.rects[ite].bottomRight.y <= points[ite].y) return 0;
    int argX = indexInSortedX[ite];
    int nowLeft = currentState.rects[ite].topLeft.y;
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) return 0;
      if (currentState.rects[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, currentState.rects[i].bottomRight.y);
        if (nowLeft >= currentState.rects[ite].bottomRight.y) { break; }
      }
    }
    nowLeft = currentState.rects[ite].topLeft.y;
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) return 0;
      if (currentState.rects[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, currentState.rects[i].bottomRight.y);
        if (nowLeft >= currentState.rects[ite].bottomRight.y) { break; }
      }
    }
  }
  return 1;
}

inline void initRect(Rect& rect, Point& point)
{
  rect.topLeft.x = point.x;
  rect.topLeft.y = point.y;
  rect.bottomRight.x = point.x + 1;
  rect.bottomRight.y = point.y + 1;
}

// 方向に応じて座標を取得するヘルパー関数
inline int getCoord(const Point& p, Direction dir)
{
  return dir == HORIZONTAL ? p.x : p.y;
}

inline int& getCoordRef(Point& p, Direction dir)
{
  return dir == HORIZONTAL ? p.x : p.y;
}

// エッジタイプから矩形の座標への参照を取得
inline int& getRectCoordByEdge(Rect& rect, int edgeType)
{
  switch (edgeType) {
    case 0: return rect.topLeft.x;
    case 1: return rect.topLeft.y;
    case 2: return rect.bottomRight.x;
    case 3: return rect.bottomRight.y;
    default: return rect.topLeft.x; // エラー回避
  }
}

inline int getSortedIndex(int ite, Direction dir)
{
  return dir == HORIZONTAL ? indexInSortedX[ite] : indexInSortedY[ite];
}

inline int getSortedRect(int idx, Direction dir)
{
  return dir == HORIZONTAL ? sortedByX[idx] : sortedByY[idx];
}

// ポイント座標による更新ヘルパー関数
inline void updateRectBoundsFromPoint(Rect& rect, int i, int ite, Direction dir)
{
  if (getCoord(points[i], dir) <= getCoord(points[ite], dir)) {
    getCoordRef(rect.topLeft, dir) = max(
      getCoord(rect.topLeft, dir),
      getCoord(points[i], dir) + 1);
  }
  else {
    getCoordRef(rect.bottomRight, dir) = min(
      getCoord(rect.bottomRight, dir),
      getCoord(points[i], dir));
  }
}

// 一方向の拡張処理を共通化（ポイント版）
inline void expandLargeInDirection(Rect& rect, int ite, Direction primaryDir)
{
  int argIdx = getSortedIndex(ite, primaryDir);

  // 後方探索
  for (int ii = argIdx - 1; ii >= 0; --ii) {
    int i = getSortedRect(ii, primaryDir);
    if (getCoord(points[i], primaryDir) == getCoord(points[ite], primaryDir)) continue;

    int hasOverlap = (primaryDir == HORIZONTAL) ?
      checkPointYOverlap(points[i].y, rect.topLeft.y, rect.bottomRight.y) :
      checkPointXOverlap(points[i].x, rect.topLeft.x, rect.bottomRight.x);

    if (hasOverlap) {
      updateRectBoundsFromPoint(rect, i, ite, primaryDir);
    }
  }

  // 前方探索
  for (int ii = argIdx + 1; ii < numRects; ++ii) {
    int i = getSortedRect(ii, primaryDir);
    if (getCoord(points[i], primaryDir) == getCoord(points[ite], primaryDir)) continue;

    int hasOverlap = (primaryDir == HORIZONTAL) ?
      checkPointYOverlap(points[i].y, rect.topLeft.y, rect.bottomRight.y) :
      checkPointXOverlap(points[i].x, rect.topLeft.x, rect.bottomRight.x);

    if (hasOverlap) {
      updateRectBoundsFromPoint(rect, i, ite, primaryDir);
    }
  }
}

// 座標更新ヘルパー関数
inline void updateRectBounds(Rect& rect, int i, int ite, Direction dir)
{
  if (getCoord(points[i], dir) <= getCoord(points[ite], dir)) {
    getCoordRef(rect.topLeft, dir) = max(
      getCoord(rect.topLeft, dir),
      getCoord(currentState.rects[i].bottomRight, dir));
  }
  else {
    getCoordRef(rect.bottomRight, dir) = min(
      getCoord(rect.bottomRight, dir),
      getCoord(currentState.rects[i].topLeft, dir));
  }
}

// 一方向の拡張処理を共通化
inline void expandInDirection(Rect& rect, int ite, Direction primaryDir)
{
  Direction secondaryDir = (primaryDir == HORIZONTAL) ? VERTICAL : HORIZONTAL;

  // 主方向の初期化
  getCoordRef(rect.topLeft, primaryDir) = 0;
  getCoordRef(rect.bottomRight, primaryDir) = 10000;

  int argIdx = getSortedIndex(ite, primaryDir);

  // 後方探索
  for (int ii = argIdx - 1; ii >= 0; --ii) {
    int i = getSortedRect(ii, primaryDir);
    int hasOverlap = (primaryDir == HORIZONTAL) ?
      checkYOverlap(currentState.rects[i], rect) : checkXOverlap(currentState.rects[i], rect);

    if (hasOverlap) {
      updateRectBounds(rect, i, ite, primaryDir);
      break;
    }
  }

  // 前方探索
  for (int ii = argIdx + 1; ii < numRects; ++ii) {
    int i = getSortedRect(ii, primaryDir);
    int hasOverlap = (primaryDir == HORIZONTAL) ?
      checkYOverlap(currentState.rects[i], rect) : checkXOverlap(currentState.rects[i], rect);

    if (hasOverlap) {
      updateRectBounds(rect, i, ite, primaryDir);
      break;
    }
  }
}

// 第2方向の拡張処理を共通化（nowLeft追跡あり）
inline void expandSecondDirection(Rect& rect, int ite, Direction primaryDir)
{
  Direction secondaryDir = (primaryDir == HORIZONTAL) ? VERTICAL : HORIZONTAL;

  // 副方向の初期化
  getCoordRef(rect.topLeft, secondaryDir) = 0;
  getCoordRef(rect.bottomRight, secondaryDir) = 10000;

  int argIdx = getSortedIndex(ite, secondaryDir);
  int nowLeft = getCoord(rect.topLeft, primaryDir);

  // 後方探索
  for (int ii = argIdx - 1; ii >= 0; --ii) {
    int i = getSortedRect(ii, secondaryDir);
    int hasOverlap = (secondaryDir == VERTICAL) ?
      checkXOverlap(currentState.rects[i], rect) : checkYOverlap(currentState.rects[i], rect);

    if (hasOverlap) {
      updateRectBounds(rect, i, ite, secondaryDir);

      if (getCoord(currentState.rects[i].topLeft, primaryDir) <= nowLeft) {
        nowLeft = max(nowLeft, getCoord(currentState.rects[i].bottomRight, primaryDir));
        if (getCoord(rect.bottomRight, primaryDir) <= nowLeft) { break; }
      }
    }
  }

  // 前方探索
  nowLeft = getCoord(rect.topLeft, primaryDir);
  for (int ii = argIdx + 1; ii < numRects; ++ii) {
    int i = getSortedRect(ii, secondaryDir);
    int hasOverlap = (secondaryDir == VERTICAL) ?
      checkXOverlap(currentState.rects[i], rect) : checkYOverlap(currentState.rects[i], rect);

    if (hasOverlap) {
      updateRectBounds(rect, i, ite, secondaryDir);

      if (getCoord(currentState.rects[i].topLeft, primaryDir) <= nowLeft) {
        nowLeft = max(nowLeft, getCoord(currentState.rects[i].bottomRight, primaryDir));
        if (getCoord(rect.bottomRight, primaryDir) <= nowLeft) { break; }
      }
    }
  }
}

// エッジ順序を取得
inline void getShuffledEdgeOrder(int edgeOrder[4])
{
  int shuffleIndex = xorshift() % 24;
  for (int j = 0; j < 4; ++j) {
    edgeOrder[j] = shuffles[shuffleIndex][j];
  }
}

// 矩形の妥当性チェックとリセット
inline bool validateAndResetRect(Rect& rect, int ite)
{
  bool isInvalid = false;

  if (!isRectInRange(rect)) isInvalid = true;
  if (rect.bottomRight.x <= rect.topLeft.x) isInvalid = true;
  if (rect.bottomRight.y <= rect.topLeft.y) isInvalid = true;
  if (points[ite].x < rect.topLeft.x || rect.bottomRight.x <= points[ite].x) isInvalid = true;
  if (points[ite].y < rect.topLeft.y || rect.bottomRight.y <= points[ite].y) isInvalid = true;

  if (isInvalid) {
    initRect(rect, points[ite]);
  }

  return !isInvalid;
}

// 矩形サイズ調整ヘルパー関数
inline void adjustRectEdge(Rect& rect, int edgeType, int targetSize, int pointCoord, int adjustAmount, bool clampDiff = false)
{
  bool isHorizontal = (edgeType == 0 || edgeType == 2);
  bool isTopLeft = (edgeType == 0 || edgeType == 1);

  // 目標サイズの計算
  int currentOtherDim = isHorizontal ?
    (rect.bottomRight.y - rect.topLeft.y) :
    (rect.bottomRight.x - rect.topLeft.x);
  int maxDim = targetSize / currentOtherDim + adjustAmount;

  // 現在の寸法と差分
  int currentDim = isHorizontal ?
    (rect.bottomRight.x - rect.topLeft.x) :
    (rect.bottomRight.y - rect.topLeft.y);
  int diff = currentDim - maxDim;
  if (clampDiff && diff < 0) diff = 0;

  // 調整可能な容量
  int capacity;
  if (isTopLeft) {
    capacity = isHorizontal ?
      (pointCoord - rect.topLeft.x) :
      (pointCoord - rect.topLeft.y);
  }
  else {
    capacity = isHorizontal ?
      (rect.bottomRight.x - (pointCoord + 1)) :
      (rect.bottomRight.y - (pointCoord + 1));
  }

  // 実際の調整量
  int adjustment = (capacity >= diff) ? diff : capacity;

  // 座標の更新
  if (isTopLeft) {
    if (isHorizontal) {
      rect.topLeft.x += adjustment;
    }
    else {
      rect.topLeft.y += adjustment;
    }
  }
  else {
    if (isHorizontal) {
      rect.bottomRight.x -= adjustment;
    }
    else {
      rect.bottomRight.y -= adjustment;
    }
  }
}

// 矩形をターゲットサイズに調整
inline void adjustRectToTargetSize(Rect& rect, int ite, bool clampDiff = false)
{
  int edgeOrder[4];
  getShuffledEdgeOrder(edgeOrder);

  int adjustAmount = xorshift() % 2;

  for (int i = 0; i < 4; ++i) {
    int area = (rect.bottomRight.x - rect.topLeft.x) * (rect.bottomRight.y - rect.topLeft.y);
    if (area <= targetSizes[ite]) break;

    int pointCoord = (edgeOrder[i] == 0 || edgeOrder[i] == 2) ? points[ite].x : points[ite].y;
    adjustRectEdge(rect, edgeOrder[i], targetSizes[ite], pointCoord, adjustAmount, clampDiff);
  }
}

inline Rect expandRect(int ite)
{
  Rect expandedRect;
  initRect(expandedRect, points[ite]);

  Direction firstDir = (Direction)(xorshift() % 2);
  Direction secondDir = (firstDir == HORIZONTAL) ? VERTICAL : HORIZONTAL;

  // 第1方向の拡張
  expandInDirection(expandedRect, ite, firstDir);

  // 第2方向の拡張
  expandInDirection(expandedRect, ite, secondDir);
  expandSecondDirection(expandedRect, ite, firstDir);

  // サイズ調整
  adjustRectToTargetSize(expandedRect, ite, false);

  // 妥当性チェック
  validateAndResetRect(expandedRect, ite);
  return expandedRect;
}

inline void saveBest()
{
  bestState.score = currentState.score;
  for (int i = 0; i < numRects; ++i) {
    bestState.rects[i] = currentState.rects[i];
  }
}

inline void extendWithTemp(int ite, double temp)
{
  Rect prevRect = currentState.rects[ite];

  Rect expandedRect = expandRect(ite);
  currentState.rects[ite] = expandedRect;

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  const double prob = exp((double)scoreDiff / temp);

  if (prob > rand01()) {
    modeCount[4]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    currentState.rects[ite] = prevRect;
    calcScore(ite);
  }
}

array<Rect, MAX_N> secondBestRects;
int secondBestScore = -1;

inline void resetRects(int numToReset)
{
  // numToReset個つぶす
  for (int i = 0; i < (numToReset); ++i) {
    int ite = xorshift() % numRects;
    initRect(currentState.rects[ite], points[ite]);
  }

  currentState.score = calcScore(-1);
  if (currentState.score > bestState.score) {
    saveBest();
  }
}

inline void resetWorstRects(int numWorstToReset)
{
  vector<pair<double, int>> v;
  for (int i = 0; i < numRects; ++i) {
    v.emplace_back(pair<double, int>(rectScores[i], i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < (numWorstToReset); ++i) {
    int ite = v[i].second;
    initRect(currentState.rects[ite], points[ite]);
  }

  currentState.score = calcScore(-1);
  if (currentState.score > bestState.score) {
    saveBest();
  }
}

inline void createHole(int hole = 100)
{
  int ite = xorshift() % numRects;
  vector<int> affectedRects;
  affectedRects.emplace_back(ite);
  currentState.rects[ite].topLeft.x -= hole;
  currentState.rects[ite].topLeft.y -= hole;
  currentState.rects[ite].bottomRight.x += hole;
  currentState.rects[ite].bottomRight.y += hole;
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) affectedRects.emplace_back(i);
  }
  int numAffected = affectedRects.size();
  for (int i = 0; i < (numAffected); ++i) {
    initRect(currentState.rects[affectedRects[i]], points[affectedRects[i]]);
  }

  currentState.score = calcScore(-1);
  if (currentState.score > bestState.score) {
    saveBest();
  }
}

inline Rect expandRectLarge(int ite)
{
  Rect largeExpandedRect;
  largeExpandedRect.topLeft.x = max(0, (int)(points[ite].x - xorshift() % 1000));
  largeExpandedRect.topLeft.y = max(0, (int)(points[ite].y - xorshift() % 1000));
  largeExpandedRect.bottomRight.x = min(10000, (int)(points[ite].x + 1 + xorshift() % 1000));
  largeExpandedRect.bottomRight.y = min(10000, (int)(points[ite].y + 1 + xorshift() % 1000));

  Direction firstDir = (Direction)(xorshift() % 2);
  Direction secondDir = (firstDir == HORIZONTAL) ? VERTICAL : HORIZONTAL;

  // 両方向に拡張
  expandLargeInDirection(largeExpandedRect, ite, firstDir);
  expandLargeInDirection(largeExpandedRect, ite, secondDir);

  // 共通のエッジ順序取得
  int edgeOrder[4];
  getShuffledEdgeOrder(edgeOrder);

  // 共通のサイズ調整
  adjustRectToTargetSize(largeExpandedRect, ite, true);

  // 共通の検証とリセット
  if (validateAndResetRect(largeExpandedRect, ite)) {
    initRect(largeExpandedRect, points[ite]);
  }

  return largeExpandedRect;
}

inline void extendLarge(int ite)
{
  Rect largeExpandedRect = expandRectLarge(ite);
  currentState.rects[ite] = largeExpandedRect;

  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) {
      initRect(currentState.rects[i], points[i]);
    }
  }

  currentState.score = calcScore(-1);
  if (currentState.score > bestState.score) {
    saveBest();
  }
}

inline void changeSingleEdge(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 101 - 50;
  int edgeType = xorshift() % 4;

  // 座標を変更
  getRectCoordByEdge(currentState.rects[ite], edgeType) += diff;

  if (isValid(ite) == 0) {
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  const double prob = exp((double)scoreDiff / temp);
  if (prob > rand01()) {
    modeCount[0]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    calcScore(ite);
  }
}

inline void changeAllEdges(int ite, double temp)
{
  int deltas[4];
  for (int i = 0; i < 4; ++i) {
    deltas[i] = xorshift() % 101 - 50;
    getRectCoordByEdge(currentState.rects[ite], i) += deltas[i];
  }

  if (isValid(ite) == 0) {
    for (int i = 0; i < 4; ++i) {
      getRectCoordByEdge(currentState.rects[ite], i) -= deltas[i];
    }
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  const double prob = exp((double)scoreDiff / temp);
  if (prob > rand01()) {
    modeCount[3]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    for (int i = 0; i < 4; ++i) {
      getRectCoordByEdge(currentState.rects[ite], i) -= deltas[i];
    }
    calcScore(ite);
  }
}

inline void slideRect(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 101 - 50;
  Direction dir = (Direction)(xorshift() % 2);

  // 両端を同じ方向に移動
  getCoordRef(currentState.rects[ite].topLeft, dir) += diff;
  getCoordRef(currentState.rects[ite].bottomRight, dir) += diff;

  if (isValid(ite) == 0) {
    getCoordRef(currentState.rects[ite].topLeft, dir) -= diff;
    getCoordRef(currentState.rects[ite].bottomRight, dir) -= diff;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  const double prob = exp((double)scoreDiff / temp);
  if (newScore >= currentState.score) {
    modeCount[1]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    getCoordRef(currentState.rects[ite].topLeft, dir) -= diff;
    getCoordRef(currentState.rects[ite].bottomRight, dir) -= diff;
    calcScore(ite);
  }
}

inline void changeAspectRatio(int ite, double temp)
{
  int widthRatio = xorshift() % 9 + 1; // 1 ~ 9;
  int heightRatio = 10 - widthRatio;

  int totalRatio = widthRatio * heightRatio;
  int scaleFactor = sqrt(targetSizes[ite] / totalRatio);
  if (scaleFactor == 0) { return; }

  int width = widthRatio * scaleFactor;
  int height = heightRatio * scaleFactor;

  Rect prevRect = currentState.rects[ite];

  int minX = max(0, points[ite].x - (width - 1));
  int maxX = min(points[ite].x, 10000 - width);
  int rangeX = maxX - minX + 1;
  if (rangeX < 1) { return; }

  int minY = max(0, points[ite].y - (height - 1));
  int maxY = min(points[ite].y, 10000 - height);
  int rangeY = maxY - minY + 1;
  if (rangeY < 1) { return; }

  currentState.rects[ite].topLeft.x = xorshift() % rangeX + minX;
  currentState.rects[ite].bottomRight.x = currentState.rects[ite].topLeft.x + rangeX;
  currentState.rects[ite].topLeft.y = xorshift() % rangeY + minY;
  currentState.rects[ite].bottomRight.y = currentState.rects[ite].topLeft.y + rangeY;

  if (isValid(ite) == 0) {
    currentState.rects[ite] = prevRect;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  const double prob = exp((double)scoreDiff / temp);
  if (newScore >= currentState.score) {
    modeCount[2]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    currentState.rects[ite] = prevRect;
    calcScore(ite);
  }
}

inline int isSelfInvalid(int ite)
{
  if (!isRectInRange(currentState.rects[ite])) return 1;
  if (currentState.rects[ite].bottomRight.x <= currentState.rects[ite].topLeft.x) return 1;
  if (currentState.rects[ite].bottomRight.y <= currentState.rects[ite].topLeft.y) return 1;
  if (points[ite].x < currentState.rects[ite].topLeft.x || currentState.rects[ite].bottomRight.x <= points[ite].x) return 1;
  if (points[ite].y < currentState.rects[ite].topLeft.y || currentState.rects[ite].bottomRight.y <= points[ite].y) return 1;
  return 0;
}

inline int canShiftBoundary(int ite, int edgeType)
{
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) {
      if (edgeType == 0) currentState.rects[i].bottomRight.x = currentState.rects[ite].topLeft.x;
      if (edgeType == 1) currentState.rects[i].bottomRight.y = currentState.rects[ite].topLeft.y;
      if (edgeType == 2) currentState.rects[i].topLeft.x = currentState.rects[ite].bottomRight.x;
      if (edgeType == 3) currentState.rects[i].topLeft.y = currentState.rects[ite].bottomRight.y;

      if (isSelfInvalid(i)) return 0;
    }
  }
  return 1;
}

array<int, MAX_N> overlappingRects;
int overlapCount;

// エッジタイプから方向情報を取得
inline Direction getEdgePrimaryDir(int edgeType)
{
  return (edgeType == 0 || edgeType == 2) ? HORIZONTAL : VERTICAL;
}

inline Direction getEdgeSecondaryDir(int edgeType)
{
  return (edgeType == 0 || edgeType == 2) ? VERTICAL : HORIZONTAL;
}

inline bool isBackwardSearch(int edgeType)
{
  return edgeType < 2;
}

// 重なり探索の共通処理
inline void findOverlapsInDirection(int ite, int edgeType)
{
  Direction primaryDir = getEdgePrimaryDir(edgeType);
  Direction secondaryDir = getEdgeSecondaryDir(edgeType);
  bool backward = isBackwardSearch(edgeType);

  int argIdx = getSortedIndex(ite, primaryDir);
  int nowLeft = getCoord(currentState.rects[ite].topLeft, secondaryDir);
  int nowRight = getCoord(currentState.rects[ite].bottomRight, secondaryDir);

  int startIdx = backward ? argIdx - 1 : argIdx + 1;
  int endIdx = backward ? -1 : numRects;
  int step = backward ? -1 : 1;

  for (int ii = startIdx; ii != endIdx; ii += step) {
    int i = getSortedRect(ii, primaryDir);
    if (checkOverlap(i, ite)) {
      // 無効な重なりチェック
      bool invalidOverlap = false;
      if (edgeType == 0 && getCoord(currentState.rects[ite].topLeft, primaryDir) <= getCoord(points[i], primaryDir)) {
        invalidOverlap = true;
      }
      else if (edgeType == 1 && getCoord(currentState.rects[ite].topLeft, primaryDir) <= getCoord(points[i], primaryDir)) {
        invalidOverlap = true;
      }
      else if (edgeType == 2 && getCoord(points[i], primaryDir) < getCoord(currentState.rects[ite].bottomRight, primaryDir)) {
        invalidOverlap = true;
      }
      else if (edgeType == 3 && getCoord(points[i], primaryDir) < getCoord(currentState.rects[ite].bottomRight, primaryDir)) {
        invalidOverlap = true;
      }

      if (invalidOverlap) {
        overlappingRects[0] = -1;
        overlapCount = 1;
        return;
      }

      overlappingRects[overlapCount] = i;
      overlapCount++;
    }

    // nowLeft/nowRightの更新
    if (getCoord(currentState.rects[i].topLeft, secondaryDir) <= nowLeft) {
      nowLeft = max(nowLeft, getCoord(currentState.rects[i].bottomRight, secondaryDir));
      if (nowLeft >= nowRight) { break; }
    }
    if (nowRight <= getCoord(currentState.rects[i].bottomRight, secondaryDir)) {
      nowRight = min(nowRight, getCoord(currentState.rects[i].topLeft, secondaryDir));
      if (nowLeft >= nowRight) { break; }
    }
  }
}

inline void findOverlaps(int ite, int edgeType)
{
  overlapCount = 0;
  findOverlapsInDirection(ite, edgeType);
}

array<Rect, MAX_N> prevRects;
inline void shiftBoundary(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 50 + 1;
  int edgeType = xorshift() % 4;

  if (edgeType < 2) diff *= -1;

  getRectCoordByEdge(currentState.rects[ite], edgeType) += diff;

  if (isSelfInvalid(ite)) {
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    return;
  }

  findOverlaps(ite, edgeType);
  int numOverlaps = overlapCount;

  if (numOverlaps > 0 && overlappingRects[0] == -1) {
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    return;
  }


  for (int i = 0; i < (numOverlaps); ++i) {
    prevRects[i] = currentState.rects[overlappingRects[i]];
  }

  int isValidShift = 1;
  for (int i = 0; i < (numOverlaps); ++i) {
    // 隣接する矩形の境界を調整
    if (edgeType < 2) {
      getRectCoordByEdge(currentState.rects[overlappingRects[i]], edgeType + 2) = getRectCoordByEdge(currentState.rects[ite], edgeType);
    }
    else {
      getRectCoordByEdge(currentState.rects[overlappingRects[i]], edgeType - 2) = getRectCoordByEdge(currentState.rects[ite], edgeType);
    }
    if (isSelfInvalid(overlappingRects[i])) isValidShift = 0;
  }

  if (isValidShift == 0) {
    for (int i = 0; i < (numOverlaps); ++i) {
      currentState.rects[overlappingRects[i]] = prevRects[i];
    }
    // 元に戻す
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    return;
  }

  for (int i = 0; i < (numOverlaps); ++i) calcScore(overlappingRects[i]);
  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentState.score;
  double prob = exp((double)scoreDiff / temp);

  if (prob > rand01()) {
    modeCount[5]++;
    currentState.score = newScore;
    if (currentState.score > bestState.score) {
      saveBest();
    }
  }
  else {
    for (int i = 0; i < (numOverlaps); ++i) {
      currentState.rects[overlappingRects[i]] = prevRects[i];
      calcScore(overlappingRects[i]);
    }
    // 元に戻す
    getRectCoordByEdge(currentState.rects[ite], edgeType) -= diff;
    calcScore(ite);
  }
}

inline void initSolution()
{
  for (int i = 0; i < numRects; ++i) {
    initRect(currentState.rects[i], points[i]);
  }

  currentState.score = calcScore(-1);
  if (currentState.score > bestState.score) {
    saveBest();
  }
}

array<Rect, MAX_N> finalBestRects;
int finalBestScore = -1;

Rect rects2[100][MAX_N];
Rect rects4[100][MAX_N];
array<int, 100> beamScores = {};

array<Rect, MAX_N> multiStartBestRects;
int multiStartBestScore = -1;

inline void multiStartSearch()
{
  Timer timer;
  for (int multiStartIter = 0; multiStartIter < (5); ++multiStartIter) {

    // 初期解
    // 左上(x,y)、右下(x+1,y+1)
    for (int i = 0; i < numRects; ++i) {
      initRect(currentState.rects[i], points[i]);
    }

    int innerIterations = 5;
    for (int innerIter = 0; innerIter < (innerIterations); ++innerIter) {
      timer.start();

      // 初期スコア計算
      currentState.score = calcScore(-1);
      saveBest();

      // 焼きなまし
      timer.start();
      double elapsedTime = timer.get_elapsed_time();
      double localTimeLimit = (0.10 / (double)innerIterations) / allLoopTimes;
      double startTemp = 2048;
      double endTemp = 0.1;
      double temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
      int loop = 0;
      while (true) {
        loop++;
        if (loop % 100 == 1) {
          elapsedTime = timer.get_elapsed_time();
          if (elapsedTime > localTimeLimit)break;
          temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
        }

        int mode = loop % 4;
        if (mode == 0) {
          int ite = xorshift() % numRects;
          changeSingleEdge(ite, temp);
        }
        else if (mode == 1) {
          int ite = xorshift() % numRects;
          slideRect(ite, temp);
        }
        else if (elapsedTime > 2.0 / innerIterations && mode == 2) {
          int ite = xorshift() % numRects;
          changeAspectRatio(ite, temp);
        }
        else if (mode == -3) {
          int ite = xorshift() % numRects;
          changeAllEdges(ite, temp);
        }
      }

      // 焼きなまし戻す
      currentState.score = bestState.score;
      for (int i = 0; i < numRects; ++i) {
        currentState.rects[i] = bestState.rects[i];
      }
      calcScore(-1);

      if (currentState.score > secondBestScore) {
        secondBestScore = currentState.score;
        for (int i = 0; i < numRects; ++i) {
          secondBestRects[i] = currentState.rects[i];
        }
      }
    }

    // secondBestScore戻す
    currentState.score = secondBestScore;
    secondBestScore = 0;
    for (int i = 0; i < numRects; ++i) {
      currentState.rects[i] = secondBestRects[i];
    }
    calcScore(-1);

    // 初期スコア計算
    currentState.score = calcScore(-1);
    saveBest();

    // 焼きなまし(2回目)
    Timer timer2;
    timer2.start();
    double elapsedTime = timer2.get_elapsed_time();
    double localTimeLimit = 0.02 / allLoopTimes;
    double startTemp = 50048;
    double endTemp = 0.1;
    double temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
    int loop = 0;
    int useSecondPhase = multiStartIter % 2;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        elapsedTime = timer2.get_elapsed_time();
        if (elapsedTime > localTimeLimit)break;
        temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
      }

      int mode = loop % 5;
      if (mode == 0) { // a,b,c,dのうち1つ変更
        int ite = xorshift() % numRects;
        changeSingleEdge(ite, temp);
      }
      else if (useSecondPhase && mode == 1) { // 位置をスライド
        int ite = xorshift() % numRects;
        slideRect(ite, temp);
      }
      else if (mode == 2 && useSecondPhase && elapsedTime > 2.0 / innerIterations) { // ランダムにアスペクト比を変更
        int ite = xorshift() % numRects;
        changeAspectRatio(ite, temp);
      }
      else if (mode == -3) { // a,b,c,dを4つ同時に変更
        int ite = xorshift() % numRects;
        changeAllEdges(ite, temp);
      }
      else if (mode == 4) { // 膨らまして縮める
        int ite = xorshift() % numRects;
        extendWithTemp(ite, temp);
      }
    }

    // 焼きなまし戻す
    currentState.score = bestState.score;
    for (int i = 0; i < numRects; ++i) {
      currentState.rects[i] = bestState.rects[i];
    }
    calcScore(-1);

    if (currentState.score > multiStartBestScore) {
      multiStartBestScore = currentState.score;
      for (int i = 0; i < numRects; ++i) {
        multiStartBestRects[i] = currentState.rects[i];
      }
    }
  }

  // 元に戻しておく
  currentState.score = 0;
  bestState.score = 0;
  secondBestScore = 0;
  for (int i = 0; i < numRects; ++i) {
    initRect(currentState.rects[i], points[i]);
    bestState.rects[i] = currentState.rects[i];
    secondBestRects[i] = currentState.rects[i];
  }
}

int solve(int isSubmission, int fileNum)
{
  auto startClock = system_clock::now();
  Timer mainTimer;
  Timer loopTimer;
  mainTimer.start();

  readInput(fileNum);

  initSortArrays();


  for (int allLoop = 0; allLoop < (allLoopTimes); ++allLoop) {

    multiStartSearch();

    // multiStartBestScore戻す
    currentState.score = multiStartBestScore;
    for (int i = 0; i < numRects; ++i) {
      currentState.rects[i] = multiStartBestRects[i];
    }
    calcScore(-1);

    // 初期スコア計算
    currentState.score = calcScore(-1);
    saveBest();

    int parentCount = 1;

    for (int beamIdx = 0; beamIdx < (parentCount); ++beamIdx) {
      for (int j = 0; j < numRects; ++j) {
        rects2[beamIdx][j] = currentState.rects[j];
      }
    }

    // 焼きなまし2
    // 筋のいいやつを追う
    int mainIterations = 250 / allLoopTimes;
    for (int mainIter = 0; mainIter < (mainIterations); ++mainIter) {
      for (int i = 0; i < (6); ++i) modeCount[i] = 0;

      int innerLoopCount = 1;
      for (int beamIdx = 0; beamIdx < (innerLoopCount); ++beamIdx) {
        int idx = beamIdx % parentCount;
        for (int i = 0; i < numRects; ++i) {
          currentState.rects[i] = rects2[idx][i];
        }

        // 初期スコア計算
        currentState.score = calcScore(-1);
        saveBest();

        // 焼きなまし(2回目)
        loopTimer.start();
        startClock = system_clock::now();
        double elapsedTime = loopTimer.get_elapsed_time();
        double localTimeLimit = (((timeLimit - 0.7) / mainIterations) / innerLoopCount) / allLoopTimes;
        double startTemp = 20048.0;
        double endTemp = 0.1;
        double temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
        int loop = 0;
        int useSecondPhase = mainIter % 2;

        while (true) {
          loop++;
          if (loop % 100 == 1) {
            const double time = duration_cast<microseconds>(system_clock::now() - startClock).count() * 1e-6;
            if (time > localTimeLimit) { break; }
            const double progressRatio = time / localTimeLimit;   // 進捗。開始時が0.0、終了時が1.0
            temp = startTemp + (endTemp - startTemp) * progressRatio;
          }


          int mode = loop % 6;

          if (mode == -1) { // a,b,c,dのうち1つ変更
            int ite = xorshift() % numRects;
            changeSingleEdge(ite, temp);
          }
          else if (mode == 1) { // 位置をスライド
            int ite = xorshift() % numRects;
            slideRect(ite, temp);
          }
          else if (mode == -2 && useSecondPhase && elapsedTime > 2.0 / mainIterations) { // ランダムにアスペクト比を変更
            int ite = xorshift() % numRects;
            changeAspectRatio(ite, temp);
          }
          else if (mode == -3) { // a,b,c,dを4つ同時に変更
            int ite = xorshift() % numRects;
            changeAllEdges(ite, temp);
          }
          else if (mode == 4) { // 膨らまして縮める
            int ite = xorshift() % numRects;
            extendWithTemp(ite, temp);
          }
          else if (mode == 5) { // 境界をずらす
            int ite = xorshift() % numRects;
            shiftBoundary(ite, temp);
          }

          if (loop % 34567 == 1120) {
            int ite = xorshift() % numRects;
            extendLarge(ite);
          }

          // 計算誤差解消?
          if (loop % 10000 == 1) {
            currentState.score = calcScore(-1);
          }
        }

        // 焼きなまし戻す
        currentState.score = bestState.score;
        for (int i = 0; i < numRects; ++i) {
          currentState.rects[i] = bestState.rects[i];
        }
        calcScore(-1);
        if (currentState.score > secondBestScore) {
          secondBestScore = currentState.score;
          for (int i = 0; i < numRects; ++i) {
            secondBestRects[i] = currentState.rects[i];
          }
        }

        // ビームサーチの次の種にする
        beamScores[beamIdx] = currentState.score;
        for (int i = 0; i < numRects; ++i) {
          rects4[beamIdx][i] = currentState.rects[i];
        }
      }

      // 次の世代に継承
      vector<P> vBeam;
      for (int beamIdx = 0; beamIdx < (innerLoopCount); ++beamIdx) vBeam.emplace_back(P(beamScores[beamIdx], beamIdx));
      sort(vBeam.begin(), vBeam.end(), greater<P>());


      for (int ii = 0; ii < (parentCount); ++ii) {
        int i = vBeam[ii].second;
        for (int j = 0; j < numRects; ++j) {
          rects2[ii][j] = rects4[i][j];
        }
      }

      // 提出時以下は消す
      if (isSubmission == 0 && mainIter % 10 == 0) {
        cout << "mainIter = " << mainIter;
        cout << ", vBeam[0] = (" << vBeam[0].first << ", " << vBeam[0].second << ")" << endl;
      }

      // エスケープ
      if (mainTimer.get_elapsed_time() > timeLimit) { break; }
    }

    // secondBestScore戻す
    currentState.score = secondBestScore;
    for (int i = 0; i < numRects; ++i) {
      currentState.rects[i] = secondBestRects[i];
    }
    calcScore(-1);

    if (isSubmission == 0) {
      cout << "currentState.score = " << currentState.score << endl;
    }

    const int MOD = 1000000007;
    if (isSubmission == 0 && currentState.score > MOD) {
      cout << "ERROR" << endl;
      writeErrorLog(fileNum);
    }

    // real_real_real入れる
    if (currentState.score > finalBestScore && currentState.score < 1000000007) {
      finalBestScore = currentState.score;
      for (int i = 0; i < numRects; ++i) {
        finalBestRects[i] = currentState.rects[i];
      }
    }


    // すべて白紙にリセットする
    currentState.score = 0;
    bestState.score = 0;
    secondBestScore = 0;
    for (int i = 0; i < numRects; ++i) {
      initRect(currentState.rects[i], points[i]);
      bestState.rects[i] = currentState.rects[i];
      secondBestRects[i] = currentState.rects[i];
    }
  }

  // finalBestScore戻す
  currentState.score = finalBestScore;
  for (int i = 0; i < numRects; ++i) {
    currentState.rects[i] = finalBestRects[i];
  }
  calcScore(-1);

  // 最終出力
  if (isSubmission) {
    for (int i = 0; i < numRects; ++i) {
      cout << currentState.rects[i].topLeft.x << ' ' << currentState.rects[i].topLeft.y << ' ' << currentState.rects[i].bottomRight.x << ' ' << currentState.rects[i].bottomRight.y << endl;
    }
  }
  else {
    writeOutput(fileNum);
  }

  // 提出時以下は消す
  if (isSubmission == 0) {
    cout << "file No. = " << fileNum << ", currentState.score = " << currentState.score << endl;
  }

  if (isSubmission == 0 && currentState.score > 1000000007) {
    writeErrorLog(fileNum);
  }

  return currentState.score;
}

inline void clearRect(Rect& r)
{
  r.topLeft.x = 0;
  r.topLeft.y = 0;
  r.bottomRight.x = 0;
  r.bottomRight.y = 0;
}

inline void clearAll()
{
  numRects = 0;
  currentState.score = -1;
  bestState.score = -1;
  totalScore = 0;
  secondBestScore = -1;
  multiStartBestScore = -1;
  finalBestScore = -1;
  for (int i = 0; i < (MAX_N); ++i) {
    points[i].x = 0, points[i].y = 0, targetSizes[i] = 0;
    currentState.rects[i].topLeft.x = 0, currentState.rects[i].topLeft.y = 0, currentState.rects[i].bottomRight.x = 0, currentState.rects[i].bottomRight.y = 0;
    clearRect(currentState.rects[i]);
    rectAreas[i] = 0;
    clearRect(bestState.rects[i]);
    rectScores[i] = 0;
    sortedByX[i] = 0, sortedByY[i] = 0;
    indexInSortedX[i] = 0, indexInSortedY[i] = 0;
    clearRect(secondBestRects[i]);
    clearRect(multiStartBestRects[i]);
    clearRect(finalBestRects[i]);
  }
}

int main()
{
  int mode = 2;
  if (mode == 0) {
    solve(1, 0);
  }
  else if (mode == 1) {
    // コードテスト用
    solve(0, 0);
  }
  else if (mode == 2) {
    // スコア確認用
    for (int i = 0; i < 10; ++i) {
      clearAll();
      solve(0, i);
    }
  }

  return 0;
}
