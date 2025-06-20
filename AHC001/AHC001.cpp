#include <algorithm>
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
int modeCount[20];

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

int allLoopTimes = 1;
int numRects;
Point points[MAX_N];
int targetSizes[MAX_N];
Rect rectangles[MAX_N];
int rectAreas[MAX_N];
Rect bestRects[MAX_N];

inline void calcArea(int idx)
{
  rectAreas[idx] = (rectangles[idx].bottomRight.x - rectangles[idx].topLeft.x) * (rectangles[idx].bottomRight.y - rectangles[idx].topLeft.y);
}

inline void readInput(int fileNum)
{
  // 入力
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << fileNum << ".txt";
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
      ofs << rectangles[i].topLeft.x << ' ' << rectangles[i].topLeft.y << ' ' << rectangles[i].bottomRight.x << ' ' << rectangles[i].bottomRight.y << endl;
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
    ofs << rectangles[i].topLeft.x << ' ' << rectangles[i].topLeft.y << ' ' << rectangles[i].bottomRight.x << ' ' << rectangles[i].bottomRight.y << endl;
  }
  ofs.close();
}

int currentScore = -1;
int bestScore = -1;

double rectScores[MAX_N];
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
    scoreSum /= (double)numRects;
    scoreSum *= 1000000000.0;
    return round(scoreSum);
  }
  else {
    double scoreSum = totalScore;
    scoreSum -= rectScores[ite];
    calcArea(ite);
    rectScores[ite] = 1.0 - (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite])) * (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite]));
    scoreSum += rectScores[ite];
    totalScore = scoreSum;
    scoreSum /= (double)numRects;
    scoreSum *= 1000000000.0;
    return round(scoreSum);
  }
}

inline int checkOverlap(int i, int j)
{
  int overlapCount = 0;
  if (rectangles[i].topLeft.x <= rectangles[j].topLeft.x && rectangles[j].topLeft.x < rectangles[i].bottomRight.x) overlapCount++;
  else if (rectangles[j].topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < rectangles[j].bottomRight.x) overlapCount++;
  if (rectangles[i].topLeft.y <= rectangles[j].topLeft.y && rectangles[j].topLeft.y < rectangles[i].bottomRight.y) overlapCount++;
  else if (rectangles[j].topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < rectangles[j].bottomRight.y) overlapCount++;
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

int sortedByX[MAX_N], sortedByY[MAX_N];
int indexInSortedX[MAX_N], indexInSortedY[MAX_N];
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

// 0~10000を出ていないか
// 面積は1以上か
// (x[i]+0.5,y[i]+0.5)を含んでいるか
// 重なりがないか
inline int isValid2(int ite)
{
  if (ite == -1) {
    for (int i = 0; i < numRects; ++i) {
      if (rectangles[i].topLeft.x < 0 || 10000 < rectangles[i].topLeft.x) return 0;
      if (rectangles[i].topLeft.y < 0 || 10000 < rectangles[i].topLeft.y) return 0;
      if (rectangles[i].bottomRight.x < 0 || 10000 < rectangles[i].bottomRight.x) return 0;
      if (rectangles[i].bottomRight.y < 0 || 10000 < rectangles[i].bottomRight.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (rectangles[i].bottomRight.x <= rectangles[i].topLeft.x) return 0;
      if (rectangles[i].bottomRight.y <= rectangles[i].topLeft.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (points[i].x < rectangles[i].topLeft.x || rectangles[i].bottomRight.x <= points[i].x) return 0;
      if (points[i].y < rectangles[i].topLeft.y || rectangles[i].bottomRight.y <= points[i].y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      for (int j = i + 1; j < numRects; ++j) {
        if (checkOverlap(i, j)) return 0;
      }
    }
  }
  else {
    if (rectangles[ite].topLeft.x < 0 || 10000 < rectangles[ite].topLeft.x) return 0;
    if (rectangles[ite].topLeft.y < 0 || 10000 < rectangles[ite].topLeft.y) return 0;
    if (rectangles[ite].bottomRight.x < 0 || 10000 < rectangles[ite].bottomRight.x) return 0;
    if (rectangles[ite].bottomRight.y < 0 || 10000 < rectangles[ite].bottomRight.y) return 0;
    if (rectangles[ite].bottomRight.x <= rectangles[ite].topLeft.x) return 0;
    if (rectangles[ite].bottomRight.y <= rectangles[ite].topLeft.y) return 0;
    if (points[ite].x < rectangles[ite].topLeft.x || rectangles[ite].bottomRight.x <= points[ite].x) return 0;
    if (points[ite].y < rectangles[ite].topLeft.y || rectangles[ite].bottomRight.y <= points[ite].y) return 0;
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
      if (rectangles[i].topLeft.x < 0 || 10000 < rectangles[i].topLeft.x) return 0;
      if (rectangles[i].topLeft.y < 0 || 10000 < rectangles[i].topLeft.y) return 0;
      if (rectangles[i].bottomRight.x < 0 || 10000 < rectangles[i].bottomRight.x) return 0;
      if (rectangles[i].bottomRight.y < 0 || 10000 < rectangles[i].bottomRight.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (rectangles[i].bottomRight.x <= rectangles[i].topLeft.x) return 0;
      if (rectangles[i].bottomRight.y <= rectangles[i].topLeft.y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      if (points[i].x < rectangles[i].topLeft.x || rectangles[i].bottomRight.x <= points[i].x) return 0;
      if (points[i].y < rectangles[i].topLeft.y || rectangles[i].bottomRight.y <= points[i].y) return 0;
    }
    for (int i = 0; i < numRects; ++i) {
      for (int j = i + 1; j < numRects; ++j) {
        if (checkOverlap(i, j)) return 0;
      }
    }
  }
  else {
    if (rectangles[ite].topLeft.x < 0 || 10000 < rectangles[ite].topLeft.x) return 0;
    if (rectangles[ite].topLeft.y < 0 || 10000 < rectangles[ite].topLeft.y) return 0;
    if (rectangles[ite].bottomRight.x < 0 || 10000 < rectangles[ite].bottomRight.x) return 0;
    if (rectangles[ite].bottomRight.y < 0 || 10000 < rectangles[ite].bottomRight.y) return 0;
    if (rectangles[ite].bottomRight.x <= rectangles[ite].topLeft.x) return 0;
    if (rectangles[ite].bottomRight.y <= rectangles[ite].topLeft.y) return 0;
    if (points[ite].x < rectangles[ite].topLeft.x || rectangles[ite].bottomRight.x <= points[ite].x) return 0;
    if (points[ite].y < rectangles[ite].topLeft.y || rectangles[ite].bottomRight.y <= points[ite].y) return 0;
    int argX = indexInSortedX[ite];
    int nowLeft = rectangles[ite].topLeft.y;
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) return 0;
      if (rectangles[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
        if (nowLeft >= rectangles[ite].bottomRight.y) { break; }
      }
    }
    nowLeft = rectangles[ite].topLeft.y;
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) return 0;
      if (rectangles[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
        if (nowLeft >= rectangles[ite].bottomRight.y) { break; }
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

Rect expandedRect;
inline void expandRect(int ite)
{
  initRect(expandedRect, points[ite]);

  int orientationFlag = xorshift() % 2;
  if (orientationFlag == 0) {
    expandedRect.topLeft.x = 0;
    expandedRect.bottomRight.x = 10000;
    int argX = indexInSortedX[ite];
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      int hasOverlap = checkYOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          expandedRect.topLeft.x = max(expandedRect.topLeft.x, rectangles[i].bottomRight.x);
        }
        else {
          expandedRect.bottomRight.x = min(expandedRect.bottomRight.x, rectangles[i].topLeft.x);
        }
        break;
      }
    }
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      int hasOverlap = checkYOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          expandedRect.topLeft.x = max(expandedRect.topLeft.x, rectangles[i].bottomRight.x);
        }
        else {
          expandedRect.bottomRight.x = min(expandedRect.bottomRight.x, rectangles[i].topLeft.x);
        }
        break;
      }
    }

    expandedRect.topLeft.y = 0;
    expandedRect.bottomRight.y = 10000;
    int argY = indexInSortedY[ite];
    int nowLeft = expandedRect.topLeft.x;
    for (int ii = argY - 1; ii >= 0; --ii) {
      int i = sortedByY[ii];
      int hasOverlap = checkXOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          expandedRect.topLeft.y = max(expandedRect.topLeft.y, rectangles[i].bottomRight.y);
        }
        else {
          expandedRect.bottomRight.y = min(expandedRect.bottomRight.y, rectangles[i].topLeft.y);
        }
        if (rectangles[i].topLeft.x <= nowLeft) {
          nowLeft = max(nowLeft, rectangles[i].bottomRight.x);
          if (expandedRect.bottomRight.x <= nowLeft) { break; }
        }
      }
    }

    nowLeft = expandedRect.topLeft.x;
    for (int ii = argY + 1; ii < numRects; ++ii) {
      int i = sortedByY[ii];
      int hasOverlap = checkXOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          expandedRect.topLeft.y = max(expandedRect.topLeft.y, rectangles[i].bottomRight.y);
        }
        else {
          expandedRect.bottomRight.y = min(expandedRect.bottomRight.y, rectangles[i].topLeft.y);
        }
        if (rectangles[i].topLeft.x <= nowLeft) {
          nowLeft = max(nowLeft, rectangles[i].bottomRight.x);
          if (expandedRect.bottomRight.x <= nowLeft) { break; }
        }
      }
    }
  }
  else {
    expandedRect.topLeft.y = 0;
    expandedRect.bottomRight.y = 10000;
    int argY = indexInSortedY[ite];
    for (int ii = argY - 1; ii >= 0; --ii) {
      int i = sortedByY[ii];
      int hasOverlap = checkXOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          expandedRect.topLeft.y = max(expandedRect.topLeft.y, rectangles[i].bottomRight.y);
        }
        else {
          expandedRect.bottomRight.y = min(expandedRect.bottomRight.y, rectangles[i].topLeft.y);
        }
        break;
      }
    }
    for (int ii = argY + 1; ii < numRects; ++ii) {
      int i = sortedByY[ii];
      int hasOverlap = checkXOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          expandedRect.topLeft.y = max(expandedRect.topLeft.y, rectangles[i].bottomRight.y);
        }
        else {
          expandedRect.bottomRight.y = min(expandedRect.bottomRight.y, rectangles[i].topLeft.y);
        }
        break;
      }
    }

    expandedRect.topLeft.x = 0;
    expandedRect.bottomRight.x = 10000;
    int argX = indexInSortedX[ite];
    int nowLeft = expandedRect.topLeft.y;
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      int hasOverlap = checkYOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          expandedRect.topLeft.x = max(expandedRect.topLeft.x, rectangles[i].bottomRight.x);
        }
        else {
          expandedRect.bottomRight.x = min(expandedRect.bottomRight.x, rectangles[i].topLeft.x);
        }
        if (rectangles[i].topLeft.y <= nowLeft) {
          nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
          if (expandedRect.bottomRight.y <= nowLeft) { break; }
        }
      }
    }
    nowLeft = expandedRect.topLeft.y;
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      int hasOverlap = checkYOverlap(rectangles[i], expandedRect);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          expandedRect.topLeft.x = max(expandedRect.topLeft.x, rectangles[i].bottomRight.x);
        }
        else {
          expandedRect.bottomRight.x = min(expandedRect.bottomRight.x, rectangles[i].topLeft.x);
        }
        if (rectangles[i].topLeft.y <= nowLeft) {
          nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
          if (expandedRect.bottomRight.y <= nowLeft) { break; }
        }
      }
    }
  }

  int edgeOrder[4] = {};
  int shuffleIndex = xorshift() % 24;
  for (int j = 0; j < (4); ++j) edgeOrder[j] = shuffles[shuffleIndex][j];

  int adjustAmount = xorshift() % 2;

  for (int i = 0; i < (4); ++i) {
    int area = (expandedRect.bottomRight.x - expandedRect.topLeft.x) * (expandedRect.bottomRight.y - expandedRect.topLeft.y);
    if (area <= targetSizes[ite]) { break; }
    if (edgeOrder[i] == 0) {
      int maxWidth = targetSizes[ite] / (expandedRect.bottomRight.y - expandedRect.topLeft.y) + adjustAmount;
      int diff = (expandedRect.bottomRight.x - expandedRect.topLeft.x) - maxWidth;
      int capacity = points[ite].x - expandedRect.topLeft.x;
      if (capacity >= diff) {
        expandedRect.topLeft.x += diff;
      }
      else {
        expandedRect.topLeft.x += capacity;
      }
    }
    else if (edgeOrder[i] == 1) {
      int maxHeight = targetSizes[ite] / (expandedRect.bottomRight.x - expandedRect.topLeft.x) + adjustAmount;
      int diff = (expandedRect.bottomRight.y - expandedRect.topLeft.y) - maxHeight;
      int capacity = points[ite].y - expandedRect.topLeft.y;
      if (capacity >= diff) {
        expandedRect.topLeft.y += diff;
      }
      else {
        expandedRect.topLeft.y += capacity;
      }
    }
    else if (edgeOrder[i] == 2) {
      int maxWidth = targetSizes[ite] / (expandedRect.bottomRight.y - expandedRect.topLeft.y) + adjustAmount;
      int diff = (expandedRect.bottomRight.x - expandedRect.topLeft.x) - maxWidth;
      int capacity = expandedRect.bottomRight.x - (points[ite].x + 1);
      if (capacity >= diff) {
        expandedRect.bottomRight.x -= diff;
      }
      else {
        expandedRect.bottomRight.x -= capacity;
      }
    }
    else if (edgeOrder[i] == 3) {
      int maxHeight = targetSizes[ite] / (expandedRect.bottomRight.x - expandedRect.topLeft.x) + adjustAmount;
      int diff = (expandedRect.bottomRight.y - expandedRect.topLeft.y) - maxHeight;
      int capacity = expandedRect.bottomRight.y - (points[ite].y + 1);
      if (capacity >= diff) {
        expandedRect.bottomRight.y -= diff;
      }
      else {
        expandedRect.bottomRight.y -= capacity;
      }
    }
  }

  int isInvalid = 0;
  if (expandedRect.topLeft.x < 0 || 10000 < expandedRect.topLeft.x) isInvalid = 1;
  if (expandedRect.topLeft.y < 0 || 10000 < expandedRect.topLeft.y) isInvalid = 1;
  if (expandedRect.bottomRight.x < 0 || 10000 < expandedRect.bottomRight.x) isInvalid = 1;
  if (expandedRect.bottomRight.y < 0 || 10000 < expandedRect.bottomRight.y) isInvalid = 1;
  if (expandedRect.bottomRight.x <= expandedRect.topLeft.x) isInvalid = 1;
  if (expandedRect.bottomRight.y <= expandedRect.topLeft.y) isInvalid = 1;
  if (points[ite].x < expandedRect.topLeft.x || expandedRect.bottomRight.x <= points[ite].x) isInvalid = 1;
  if (points[ite].y < expandedRect.topLeft.y || expandedRect.bottomRight.y <= points[ite].y) isInvalid = 1;

  if (isInvalid) {
    initRect(expandedRect, points[ite]);
  }
}

inline void saveBest()
{
  bestScore = currentScore;
  for (int i = 0; i < numRects; ++i) {
    bestRects[i] = rectangles[i];
  }
}

inline void extendWithTemp(int ite, double temp)
{
  Rect prevRect = rectangles[ite];

  expandRect(ite);
  rectangles[ite] = expandedRect;

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  const double prob = exp((double)scoreDiff / temp);

  if (prob > rand01()) {
    modeCount[4]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite] = prevRect;
    calcScore(ite);
  }
}

Rect secondBestRects[MAX_N];
int secondBestScore = -1;

inline void resetRects(int numToReset)
{
  // numToReset個つぶす
  for (int i = 0; i < (numToReset); ++i) {
    int ite = xorshift() % numRects;
    initRect(rectangles[ite], points[ite]);
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
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
    initRect(rectangles[ite], points[ite]);
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
    saveBest();
  }
}

inline void createHole(int hole = 100)
{
  int ite = xorshift() % numRects;
  vector<int> affectedRects;
  affectedRects.emplace_back(ite);
  rectangles[ite].topLeft.x -= hole;
  rectangles[ite].topLeft.y -= hole;
  rectangles[ite].bottomRight.x += hole;
  rectangles[ite].bottomRight.y += hole;
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) affectedRects.emplace_back(i);
  }
  int numAffected = affectedRects.size();
  for (int i = 0; i < (numAffected); ++i) {
    initRect(rectangles[affectedRects[i]], points[affectedRects[i]]);
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
    saveBest();
  }
}

Rect largeExpandedRect;
inline void expandRectLarge(int ite)
{
  largeExpandedRect.topLeft.x = max(0, (int)(points[ite].x - xorshift() % 1000));
  largeExpandedRect.topLeft.y = max(0, (int)(points[ite].y - xorshift() % 1000));
  largeExpandedRect.bottomRight.x = min(10000, (int)(points[ite].x + 1 + xorshift() % 1000));
  largeExpandedRect.bottomRight.y = min(10000, (int)(points[ite].y + 1 + xorshift() % 1000));

  int isVerticalFirst = xorshift() % 2;

  if (isVerticalFirst == 0) {
    int argX = indexInSortedX[ite];
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (points[i].x == points[ite].x) { continue; }
      int hasOverlap = checkPointYOverlap(points[i].y, largeExpandedRect.topLeft.y, largeExpandedRect.bottomRight.y);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          largeExpandedRect.topLeft.x = max(largeExpandedRect.topLeft.x, points[i].x + 1);
        }
        else {
          largeExpandedRect.bottomRight.x = min(largeExpandedRect.bottomRight.x, points[i].x);
        }
      }
    }
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      if (points[i].x == points[ite].x) { continue; }
      int hasOverlap = checkPointYOverlap(points[i].y, largeExpandedRect.topLeft.y, largeExpandedRect.bottomRight.y);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          largeExpandedRect.topLeft.x = max(largeExpandedRect.topLeft.x, points[i].x + 1);
        }
        else {
          largeExpandedRect.bottomRight.x = min(largeExpandedRect.bottomRight.x, points[i].x);
        }
      }
    }

    int argY = indexInSortedY[ite];
    for (int ii = argY - 1; ii >= 0; --ii) {
      int i = sortedByY[ii];
      if (points[i].y == points[ite].y) { continue; }
      int hasOverlap = checkPointXOverlap(points[i].x, largeExpandedRect.topLeft.x, largeExpandedRect.bottomRight.x);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          largeExpandedRect.topLeft.y = max(largeExpandedRect.topLeft.y, points[i].y + 1);
        }
        else {
          largeExpandedRect.bottomRight.y = min(largeExpandedRect.bottomRight.y, points[i].y);
        }
      }
    }
    for (int ii = argY + 1; ii < numRects; ++ii) {
      int i = sortedByY[ii];
      if (points[i].y == points[ite].y) { continue; }
      int hasOverlap = checkPointXOverlap(points[i].x, largeExpandedRect.topLeft.x, largeExpandedRect.bottomRight.x);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          largeExpandedRect.topLeft.y = max(largeExpandedRect.topLeft.y, points[i].y + 1);
        }
        else {
          largeExpandedRect.bottomRight.y = min(largeExpandedRect.bottomRight.y, points[i].y);
        }
      }
    }
  }
  else {
    int argY = indexInSortedY[ite];
    for (int ii = argY - 1; ii >= 0; --ii) {
      int i = sortedByY[ii];
      if (points[i].y == points[ite].y) { continue; }
      int hasOverlap = checkPointXOverlap(points[i].x, largeExpandedRect.topLeft.x, largeExpandedRect.bottomRight.x);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          largeExpandedRect.topLeft.y = max(largeExpandedRect.topLeft.y, points[i].y + 1);
        }
        else {
          largeExpandedRect.bottomRight.y = min(largeExpandedRect.bottomRight.y, points[i].y);
        }
      }
    }
    for (int ii = argY + 1; ii < numRects; ++ii) {
      int i = sortedByY[ii];
      if (points[i].y == points[ite].y) { continue; }
      int hasOverlap = checkPointXOverlap(points[i].x, largeExpandedRect.topLeft.x, largeExpandedRect.bottomRight.x);
      if (hasOverlap) {
        if (points[i].y <= points[ite].y) {
          largeExpandedRect.topLeft.y = max(largeExpandedRect.topLeft.y, points[i].y + 1);
        }
        else {
          largeExpandedRect.bottomRight.y = min(largeExpandedRect.bottomRight.y, points[i].y);
        }
      }
    }

    int argX = indexInSortedX[ite];
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (points[i].x == points[ite].x) { continue; }
      int hasOverlap = checkPointYOverlap(points[i].y, largeExpandedRect.topLeft.y, largeExpandedRect.bottomRight.y);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          largeExpandedRect.topLeft.x = max(largeExpandedRect.topLeft.x, points[i].x + 1);
        }
        else {
          largeExpandedRect.bottomRight.x = min(largeExpandedRect.bottomRight.x, points[i].x);
        }
      }
    }
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      if (points[i].x == points[ite].x) { continue; }
      int hasOverlap = checkPointYOverlap(points[i].y, largeExpandedRect.topLeft.y, largeExpandedRect.bottomRight.y);
      if (hasOverlap) {
        if (points[i].x <= points[ite].x) {
          largeExpandedRect.topLeft.x = max(largeExpandedRect.topLeft.x, points[i].x + 1);
        }
        else {
          largeExpandedRect.bottomRight.x = min(largeExpandedRect.bottomRight.x, points[i].x);
        }
      }
    }
  }

  int edgeOrder[4] = {};
  int shuffleIndex = xorshift() % 24;
  for (int j = 0; j < (4); ++j) edgeOrder[j] = shuffles[shuffleIndex][j];

  int adjustAmount = xorshift() % 2;

  for (int i = 0; i < (4); ++i) {
    int area = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) * (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y);
    if (area <= targetSizes[ite]) { break; }
    if (edgeOrder[i] == 0) {
      int maxWidth = targetSizes[ite] / (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) + adjustAmount;
      int diff = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) - maxWidth;
      if (diff < 0) diff = 0;
      int capacity = points[ite].x - largeExpandedRect.topLeft.x;
      if (capacity >= diff) {
        largeExpandedRect.topLeft.x += diff;
      }
      else {
        largeExpandedRect.topLeft.x += capacity;
      }
    }
    else if (edgeOrder[i] == 1) {
      int maxHeight = targetSizes[ite] / (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) + adjustAmount;
      int diff = (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) - maxHeight;
      if (diff < 0) diff = 0;
      int capacity = points[ite].y - largeExpandedRect.topLeft.y;
      if (capacity >= diff) {
        largeExpandedRect.topLeft.y += diff;
      }
      else {
        largeExpandedRect.topLeft.y += capacity;
      }
    }
    else if (edgeOrder[i] == 2) {
      int maxWidth = targetSizes[ite] / (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) + adjustAmount;
      int diff = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) - maxWidth;
      if (diff < 0) diff = 0;
      int capacity = largeExpandedRect.bottomRight.x - (points[ite].x + 1);
      if (capacity >= diff) {
        largeExpandedRect.bottomRight.x -= diff;
      }
      else {
        largeExpandedRect.bottomRight.x -= capacity;
      }
    }
    else if (edgeOrder[i] == 3) {
      int maxHeight = targetSizes[ite] / (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) + adjustAmount;
      int diff = (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) - maxHeight;
      if (diff < 0) diff = 0;
      int capacity = largeExpandedRect.bottomRight.y - (points[ite].y + 1);
      if (capacity >= diff) {
        largeExpandedRect.bottomRight.y -= diff;
      }
      else {
        largeExpandedRect.bottomRight.y -= capacity;
      }
    }
  }

  int isInvalid = 0;
  if (largeExpandedRect.topLeft.x < 0 || 10000 < largeExpandedRect.topLeft.x) isInvalid = 1;
  if (largeExpandedRect.topLeft.y < 0 || 10000 < largeExpandedRect.topLeft.y) isInvalid = 1;
  if (largeExpandedRect.bottomRight.x < 0 || 10000 < largeExpandedRect.bottomRight.x) isInvalid = 1;
  if (largeExpandedRect.bottomRight.y < 0 || 10000 < largeExpandedRect.bottomRight.y) isInvalid = 1;
  if (largeExpandedRect.bottomRight.x <= largeExpandedRect.topLeft.x) isInvalid = 1;
  if (largeExpandedRect.bottomRight.y <= largeExpandedRect.topLeft.y) isInvalid = 1;
  if (points[ite].x < largeExpandedRect.topLeft.x || largeExpandedRect.bottomRight.x <= points[ite].x) isInvalid = 1;
  if (points[ite].y < largeExpandedRect.topLeft.y || largeExpandedRect.bottomRight.y <= points[ite].y) isInvalid = 1;

  if (isInvalid) {
    initRect(largeExpandedRect, points[ite]);
  }

}

inline void extendLarge(int ite)
{
  expandRectLarge(ite);
  rectangles[ite] = largeExpandedRect;

  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) {
      initRect(rectangles[i], points[i]);
    }
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
    saveBest();
  }
}

inline void changeSingleEdge(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 101 - 50;
  int edgeType = xorshift() % 4;

  if (edgeType == 0) rectangles[ite].topLeft.x += diff;
  if (edgeType == 1) rectangles[ite].topLeft.y += diff;
  if (edgeType == 2) rectangles[ite].bottomRight.x += diff;
  if (edgeType == 3) rectangles[ite].bottomRight.y += diff;

  if (isValid(ite) == 0) {
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  const double prob = exp((double)scoreDiff / temp);
  if (prob > rand01()) {
    modeCount[0]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    calcScore(ite);
  }
}

inline void changeAllEdges(int ite, double temp)
{
  int deltaTopX = xorshift() % 101 - 50;
  int deltaTopY = xorshift() % 101 - 50;
  int deltaBottomX = xorshift() % 101 - 50;
  int deltaBottomY = xorshift() % 101 - 50;

  rectangles[ite].topLeft.x += deltaTopX;
  rectangles[ite].topLeft.y += deltaTopY;
  rectangles[ite].bottomRight.x += deltaBottomX;
  rectangles[ite].bottomRight.y += deltaBottomY;

  if (isValid(ite) == 0) {
    rectangles[ite].topLeft.x -= deltaTopX;
    rectangles[ite].topLeft.y -= deltaTopY;
    rectangles[ite].bottomRight.x -= deltaBottomX;
    rectangles[ite].bottomRight.y -= deltaBottomY;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  const double prob = exp((double)scoreDiff / temp);
  if (prob > rand01()) {
    modeCount[3]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite].topLeft.x -= deltaTopX;
    rectangles[ite].topLeft.y -= deltaTopY;
    rectangles[ite].bottomRight.x -= deltaBottomX;
    rectangles[ite].bottomRight.y -= deltaBottomY;
    calcScore(ite);
  }
}

inline void slideRect(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 101 - 50;
  int direction = xorshift() % 2;

  if (direction == 0) {
    rectangles[ite].topLeft.x += diff;
    rectangles[ite].bottomRight.x += diff;
  }
  if (direction == 1) {
    rectangles[ite].topLeft.y += diff;
    rectangles[ite].bottomRight.y += diff;
  }

  if (isValid(ite) == 0) {
    if (direction == 0) {
      rectangles[ite].topLeft.x -= diff;
      rectangles[ite].bottomRight.x -= diff;
    }
    if (direction == 1) {
      rectangles[ite].topLeft.y -= diff;
      rectangles[ite].bottomRight.y -= diff;
    }
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  const double prob = exp((double)scoreDiff / temp);
  if (newScore >= currentScore) {
    modeCount[1]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    if (direction == 0) {
      rectangles[ite].topLeft.x -= diff;
      rectangles[ite].bottomRight.x -= diff;
    }
    if (direction == 1) {
      rectangles[ite].topLeft.y -= diff;
      rectangles[ite].bottomRight.y -= diff;
    }
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

  Rect prevRect = rectangles[ite];

  int minX = max(0, points[ite].x - (width - 1));
  int maxX = min(points[ite].x, 10000 - width);
  int rangeX = maxX - minX + 1;
  if (rangeX < 1) { return; }

  int minY = max(0, points[ite].y - (height - 1));
  int maxY = min(points[ite].y, 10000 - height);
  int rangeY = maxY - minY + 1;
  if (rangeY < 1) { return; }

  rectangles[ite].topLeft.x = xorshift() % rangeX + minX;
  rectangles[ite].bottomRight.x = rectangles[ite].topLeft.x + rangeX;
  rectangles[ite].topLeft.y = xorshift() % rangeY + minY;
  rectangles[ite].bottomRight.y = rectangles[ite].topLeft.y + rangeY;

  if (isValid(ite) == 0) {
    rectangles[ite] = prevRect;
    return;
  }

  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  const double prob = exp((double)scoreDiff / temp);
  if (newScore >= currentScore) {
    modeCount[2]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite] = prevRect;
    calcScore(ite);
  }
}

inline int isSelfInvalid(int ite)
{
  if (rectangles[ite].topLeft.x < 0 || 10000 < rectangles[ite].topLeft.x) return 1;
  if (rectangles[ite].topLeft.y < 0 || 10000 < rectangles[ite].topLeft.y) return 1;
  if (rectangles[ite].bottomRight.x < 0 || 10000 < rectangles[ite].bottomRight.x) return 1;
  if (rectangles[ite].bottomRight.y < 0 || 10000 < rectangles[ite].bottomRight.y) return 1;
  if (rectangles[ite].bottomRight.x <= rectangles[ite].topLeft.x) return 1;
  if (rectangles[ite].bottomRight.y <= rectangles[ite].topLeft.y) return 1;
  if (points[ite].x < rectangles[ite].topLeft.x || rectangles[ite].bottomRight.x <= points[ite].x) return 1;
  if (points[ite].y < rectangles[ite].topLeft.y || rectangles[ite].bottomRight.y <= points[ite].y) return 1;
  return 0;
}

inline int canShiftBoundary(int ite, int edgeType)
{
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) {
      if (edgeType == 0) rectangles[i].bottomRight.x = rectangles[ite].topLeft.x;
      if (edgeType == 1) rectangles[i].bottomRight.y = rectangles[ite].topLeft.y;
      if (edgeType == 2) rectangles[i].topLeft.x = rectangles[ite].bottomRight.x;
      if (edgeType == 3) rectangles[i].topLeft.y = rectangles[ite].bottomRight.y;

      if (isSelfInvalid(i)) return 0;
    }
  }
  return 1;
}

int overlappingRects[MAX_N];
int overlapCount;
inline void findOverlaps(int ite, int edgeType)
{
  overlapCount = 0;
  if (edgeType == 0) {
    int argX = indexInSortedX[ite];
    int nowLeft = rectangles[ite].topLeft.y;
    int nowRight = rectangles[ite].bottomRight.y;
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) {
        if (rectangles[ite].topLeft.x <= points[i].x) {
          overlappingRects[0] = -1;
          overlapCount = 1;
          return;
        }
        overlappingRects[overlapCount] = i;
        overlapCount++;
      }
      if (rectangles[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
        if (nowLeft >= nowRight) { break; }
      }
      if (nowRight <= rectangles[i].bottomRight.y) {
        nowRight = min(nowRight, rectangles[i].topLeft.y);
        if (nowLeft >= nowRight) { break; }
      }
    }
  }
  if (edgeType == 1) {
    int argY = indexInSortedY[ite];
    int nowLeft = rectangles[ite].topLeft.x;
    int nowRight = rectangles[ite].bottomRight.x;
    for (int ii = argY - 1; ii >= 0; --ii) {
      int i = sortedByY[ii];
      if (checkOverlap(i, ite)) {
        if (rectangles[ite].topLeft.y <= points[i].y) {
          overlappingRects[0] = -1;
          overlapCount = 1;
          return;
        }
        overlappingRects[overlapCount] = i;
        overlapCount++;
      }
      if (rectangles[i].topLeft.x <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.x);
        if (nowLeft >= nowRight) { break; }
      }
      if (nowRight <= rectangles[i].bottomRight.x) {
        nowRight = min(nowRight, rectangles[i].topLeft.x);
        if (nowLeft >= nowRight) { break; }
      }
    }
  }
  if (edgeType == 2) {
    int argX = indexInSortedX[ite];
    int nowLeft = rectangles[ite].topLeft.y;
    int nowRight = rectangles[ite].bottomRight.y;
    for (int ii = argX + 1; ii < numRects; ++ii) {
      int i = sortedByX[ii];
      if (checkOverlap(i, ite)) {
        if (points[i].x < rectangles[ite].bottomRight.x) {
          overlappingRects[0] = -1;
          overlapCount = 1;
          return;
        }
        overlappingRects[overlapCount] = i;
        overlapCount++;
      }
      if (rectangles[i].topLeft.y <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.y);
        if (nowLeft >= nowRight) { break; }
      }
      if (nowRight <= rectangles[i].bottomRight.y) {
        nowRight = min(nowRight, rectangles[i].topLeft.y);
        if (nowLeft >= nowRight) { break; }
      }
    }
  }
  if (edgeType == 3) {
    int argY = indexInSortedY[ite];
    int nowLeft = rectangles[ite].topLeft.x;
    int nowRight = rectangles[ite].bottomRight.x;
    for (int ii = argY + 1; ii < numRects; ++ii) {
      int i = sortedByY[ii];
      if (checkOverlap(i, ite)) {
        if (points[i].y < rectangles[ite].bottomRight.y) {
          overlappingRects[0] = -1;
          overlapCount = 1;
          return;
        }
        overlappingRects[overlapCount] = i;
        overlapCount++;
      }
      if (rectangles[i].topLeft.x <= nowLeft) {
        nowLeft = max(nowLeft, rectangles[i].bottomRight.x);
        if (nowLeft >= nowRight) { break; }
      }
      if (nowRight <= rectangles[i].bottomRight.x) {
        nowRight = min(nowRight, rectangles[i].topLeft.x);
        if (nowLeft >= nowRight) { break; }
      }
    }
  }
}

Rect prevRects[MAX_N];
inline void shiftBoundary(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 50 + 1;
  int edgeType = xorshift() % 4;

  if (edgeType < 2) diff *= -1;

  if (edgeType == 0) rectangles[ite].topLeft.x += diff;
  if (edgeType == 1) rectangles[ite].topLeft.y += diff;
  if (edgeType == 2) rectangles[ite].bottomRight.x += diff;
  if (edgeType == 3) rectangles[ite].bottomRight.y += diff;

  if (isSelfInvalid(ite)) {
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  findOverlaps(ite, edgeType);
  int numOverlaps = overlapCount;

  if (numOverlaps > 0 && overlappingRects[0] == -1) {
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }


  for (int i = 0; i < (numOverlaps); ++i) {
    prevRects[i] = rectangles[overlappingRects[i]];
  }

  int isValidShift = 1;
  for (int i = 0; i < (numOverlaps); ++i) {
    if (edgeType == 0) rectangles[overlappingRects[i]].bottomRight.x = rectangles[ite].topLeft.x;
    if (edgeType == 1) rectangles[overlappingRects[i]].bottomRight.y = rectangles[ite].topLeft.y;
    if (edgeType == 2) rectangles[overlappingRects[i]].topLeft.x = rectangles[ite].bottomRight.x;
    if (edgeType == 3) rectangles[overlappingRects[i]].topLeft.y = rectangles[ite].bottomRight.y;
    if (isSelfInvalid(overlappingRects[i])) isValidShift = 0;
  }

  if (isValidShift == 0) {
    for (int i = 0; i < (numOverlaps); ++i) {
      rectangles[overlappingRects[i]] = prevRects[i];
    }
    // 元に戻す
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  for (int i = 0; i < (numOverlaps); ++i) calcScore(overlappingRects[i]);
  int newScore = calcScore(ite);

  int scoreDiff = newScore - currentScore;
  double prob = exp((double)scoreDiff / temp);

  if (prob > rand01()) {
    modeCount[5]++;
    currentScore = newScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    for (int i = 0; i < (numOverlaps); ++i) {
      rectangles[overlappingRects[i]] = prevRects[i];
      calcScore(overlappingRects[i]);
    }
    // 元に戻す
    if (edgeType == 0) rectangles[ite].topLeft.x -= diff;
    if (edgeType == 1) rectangles[ite].topLeft.y -= diff;
    if (edgeType == 2) rectangles[ite].bottomRight.x -= diff;
    if (edgeType == 3) rectangles[ite].bottomRight.y -= diff;
    calcScore(ite);
  }
}

inline void initSolution()
{
  for (int i = 0; i < numRects; ++i) {
    initRect(rectangles[i], points[i]);
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
    saveBest();
  }
}

Rect finalBestRects[MAX_N];
int finalBestScore = -1;

Rect rects2[100][MAX_N];
Rect rects4[100][MAX_N];
int beamScores[100] = {};

Rect multiStartBestRects[MAX_N];
int multiStartBestScore = -1;

inline void multiStartSearch()
{
  clock_t start, end;
  for (int multiStartIter = 0; multiStartIter < (5); ++multiStartIter) {

    // 初期解
    // 左上(x,y)、右下(x+1,y+1)
    for (int i = 0; i < numRects; ++i) {
      initRect(rectangles[i], points[i]);
    }

    int innerIterations = 5;
    for (int innerIter = 0; innerIter < (innerIterations); ++innerIter) {
      start = clock();

      // 初期スコア計算
      currentScore = calcScore(-1);
      saveBest();

      // 焼きなまし
      start = clock();
      end = clock();
      double elapsedTime = ((double)end - start) / CLOCKS_PER_SEC;
      double localTimeLimit = (0.10 / (double)innerIterations) / allLoopTimes;
      double startTemp = 2048;
      double endTemp = 0.1;
      double temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
      int loop = 0;
      while (true) {
        loop++;
        if (loop % 100 == 1) {
          end = clock();
          elapsedTime = ((double)end - start) / CLOCKS_PER_SEC;
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
      currentScore = bestScore;
      for (int i = 0; i < numRects; ++i) {
        rectangles[i] = bestRects[i];
      }
      calcScore(-1);

      if (currentScore > secondBestScore) {
        secondBestScore = currentScore;
        for (int i = 0; i < numRects; ++i) {
          secondBestRects[i] = rectangles[i];
        }
      }
    }

    // secondBestScore戻す
    currentScore = secondBestScore;
    secondBestScore = 0;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = secondBestRects[i];
    }
    calcScore(-1);

    // 初期スコア計算
    currentScore = calcScore(-1);
    bestScore = currentScore;

    for (int i = 0; i < numRects; ++i) {
      bestRects[i] = rectangles[i];
    }


    // 焼きなまし(2回目)
    clock_t start, end;
    start = clock();
    end = clock();
    double elapsedTime = ((double)end - start) / CLOCKS_PER_SEC;
    double localTimeLimit = 0.02 / allLoopTimes;
    double startTemp = 50048;
    double endTemp = 0.1;
    double temp = startTemp + (endTemp - startTemp) * elapsedTime / localTimeLimit;
    int loop = 0;
    int useSecondPhase = multiStartIter % 2;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        end = clock();
        elapsedTime = ((double)end - start) / CLOCKS_PER_SEC;
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
    currentScore = bestScore;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = bestRects[i];
    }
    calcScore(-1);

    if (currentScore > multiStartBestScore) {
      multiStartBestScore = currentScore;
      for (int i = 0; i < numRects; ++i) {
        multiStartBestRects[i] = rectangles[i];
      }
    }
  }

  // 元に戻しておく
  currentScore = 0;
  bestScore = 0;
  secondBestScore = 0;
  for (int i = 0; i < numRects; ++i) {
    initRect(rectangles[i], points[i]);
    bestRects[i] = rectangles[i];
    secondBestRects[i] = rectangles[i];
  }
}

int solve(int isSubmission, int fileNum)
{
  auto startClock = system_clock::now();
  clock_t start, end;
  clock_t overallStart = clock();

  readInput(fileNum);

  initSortArrays();


  for (int allLoop = 0; allLoop < (allLoopTimes); ++allLoop) {

    multiStartSearch();

    // multiStartBestScore戻す
    currentScore = multiStartBestScore;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = multiStartBestRects[i];
    }
    calcScore(-1);

    // 初期スコア計算
    currentScore = calcScore(-1);
    bestScore = currentScore;

    for (int i = 0; i < numRects; ++i) {
      bestRects[i] = rectangles[i];
    }


    int parentCount = 1;

    for (int beamIdx = 0; beamIdx < (parentCount); ++beamIdx) {
      for (int j = 0; j < numRects; ++j) {
        rects2[beamIdx][j] = rectangles[j];
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
          rectangles[i] = rects2[idx][i];
        }

        // 初期スコア計算
        currentScore = calcScore(-1);
        bestScore = currentScore;

        for (int i = 0; i < numRects; ++i) {
          bestRects[i] = rectangles[i];
        }

        // 焼きなまし(2回目)
        start = clock();
        end = clock();
        startClock = system_clock::now();
        double elapsedTime = ((double)end - start) / CLOCKS_PER_SEC;
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
            currentScore = calcScore(-1);
          }
        }

        // 焼きなまし戻す
        currentScore = bestScore;
        for (int i = 0; i < numRects; ++i) {
          rectangles[i] = bestRects[i];
        }
        calcScore(-1);
        if (currentScore > secondBestScore) {
          secondBestScore = currentScore;
          for (int i = 0; i < numRects; ++i) {
            secondBestRects[i] = rectangles[i];
          }
        }

        // ビームサーチの次の種にする
        beamScores[beamIdx] = currentScore;
        for (int i = 0; i < numRects; ++i) {
          rects4[beamIdx][i] = rectangles[i];
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
      end = clock();
      if (((double)end - overallStart) / CLOCKS_PER_SEC > timeLimit) { break; }
    }

    // secondBestScore戻す
    currentScore = secondBestScore;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = secondBestRects[i];
    }
    calcScore(-1);

    if (isSubmission == 0) {
      cout << "currentScore = " << currentScore << endl;
    }

    const int MOD = 1000000007;
    if (isSubmission == 0 && currentScore > MOD) {
      cout << "ERROR" << endl;
      writeErrorLog(fileNum);
    }

    // real_real_real入れる
    if (currentScore > finalBestScore && currentScore < 1000000007) {
      finalBestScore = currentScore;
      for (int i = 0; i < numRects; ++i) {
        finalBestRects[i] = rectangles[i];
      }
    }


    // すべて白紙にリセットする
    currentScore = 0;
    bestScore = 0;
    secondBestScore = 0;
    for (int i = 0; i < numRects; ++i) {
      initRect(rectangles[i], points[i]);
      bestRects[i] = rectangles[i];
      secondBestRects[i] = rectangles[i];
    }
  }

  // finalBestScore戻す
  currentScore = finalBestScore;
  for (int i = 0; i < numRects; ++i) {
    rectangles[i] = finalBestRects[i];
  }
  calcScore(-1);

  // 最終出力
  if (isSubmission) {
    for (int i = 0; i < numRects; ++i) {
      cout << rectangles[i].topLeft.x << ' ' << rectangles[i].topLeft.y << ' ' << rectangles[i].bottomRight.x << ' ' << rectangles[i].bottomRight.y << endl;
    }
  }
  else {
    writeOutput(fileNum);
  }

  // 提出時以下は消す
  if (isSubmission == 0) {
    cout << "file No. = " << fileNum << ", currentScore = " << currentScore << endl;
  }

  if (isSubmission == 0 && currentScore > 1000000007) {
    writeErrorLog(fileNum);
  }

  return currentScore;
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
  currentScore = -1;
  bestScore = -1;
  totalScore = 0;
  secondBestScore = -1;
  multiStartBestScore = -1;
  finalBestScore = -1;
  for (int i = 0; i < (MAX_N); ++i) {
    points[i].x = 0, points[i].y = 0, targetSizes[i] = 0;
    rectangles[i].topLeft.x = 0, rectangles[i].topLeft.y = 0, rectangles[i].bottomRight.x = 0, rectangles[i].bottomRight.y = 0;
    clearRect(rectangles[i]);
    rectAreas[i] = 0;
    clearRect(bestRects[i]);
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
