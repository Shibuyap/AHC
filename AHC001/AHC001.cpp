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
    double sum = 0;
    for (int i = 0; i < numRects; ++i) {
      calcArea(i);
      rectScores[i] = 1.0 - (1.0 - (double)min(targetSizes[i], rectAreas[i]) / (double)max(targetSizes[i], rectAreas[i])) * (1.0 - (double)min(targetSizes[i], rectAreas[i]) / (double)max(targetSizes[i], rectAreas[i]));
      sum += rectScores[i];
    }
    totalScore = sum;
    sum /= (double)numRects;
    sum *= 1000000000.0;
    return round(sum);
  }
  else {
    double sum = totalScore;
    sum -= rectScores[ite];
    calcArea(ite);
    rectScores[ite] = 1.0 - (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite])) * (1.0 - (double)min(targetSizes[ite], rectAreas[ite]) / (double)max(targetSizes[ite], rectAreas[ite]));
    sum += rectScores[ite];
    totalScore = sum;
    sum /= (double)numRects;
    sum *= 1000000000.0;
    return round(sum);
  }
}

inline int checkOverlap(int i, int j)
{
  int cnt = 0;
  if (rectangles[i].topLeft.x <= rectangles[j].topLeft.x && rectangles[j].topLeft.x < rectangles[i].bottomRight.x) cnt++;
  else if (rectangles[j].topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < rectangles[j].bottomRight.x) cnt++;
  if (rectangles[i].topLeft.y <= rectangles[j].topLeft.y && rectangles[j].topLeft.y < rectangles[i].bottomRight.y) cnt++;
  else if (rectangles[j].topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < rectangles[j].bottomRight.y) cnt++;
  return cnt == 2;
}

int sortedByX[MAX_N], sortedByY[MAX_N];
int indexInSortedX[MAX_N], indexInSortedY[MAX_N];
inline void initSortArrays()
{
  vector<P> v;
  for (int i = 0; i < numRects; ++i) {
    v.emplace_back(P(points[i].x, i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < numRects; ++i) {
    sortedByX[i] = v[i].second;
    indexInSortedX[v[i].second] = i;
  }

  v.clear();
  for (int i = 0; i < numRects; ++i) {
    v.emplace_back(P(points[i].y, i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < numRects; ++i) {
    sortedByY[i] = v[i].second;
    indexInSortedY[v[i].second] = i;
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

  int flagTateYoko = xorshift() % 2;
  if (flagTateYoko == 0) {
    expandedRect.topLeft.x = 0;
    expandedRect.bottomRight.x = 10000;
    int argX = indexInSortedX[ite];
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      int flagKasanari = 0;
      if (rectangles[i].topLeft.y <= expandedRect.topLeft.y && expandedRect.topLeft.y < rectangles[i].bottomRight.y) flagKasanari = 1;
      if (expandedRect.topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < expandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.y <= expandedRect.topLeft.y && expandedRect.topLeft.y < rectangles[i].bottomRight.y) flagKasanari = 1;
      if (expandedRect.topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < expandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.x <= expandedRect.topLeft.x && expandedRect.topLeft.x < rectangles[i].bottomRight.x) flagKasanari = 1;
      if (expandedRect.topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < expandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.x <= expandedRect.topLeft.x && expandedRect.topLeft.x < rectangles[i].bottomRight.x) flagKasanari = 1;
      if (expandedRect.topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < expandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.x <= expandedRect.topLeft.x && expandedRect.topLeft.x < rectangles[i].bottomRight.x) flagKasanari = 1;
      if (expandedRect.topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < expandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.x <= expandedRect.topLeft.x && expandedRect.topLeft.x < rectangles[i].bottomRight.x) flagKasanari = 1;
      if (expandedRect.topLeft.x <= rectangles[i].topLeft.x && rectangles[i].topLeft.x < expandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.y <= expandedRect.topLeft.y && expandedRect.topLeft.y < rectangles[i].bottomRight.y) flagKasanari = 1;
      if (expandedRect.topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < expandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (rectangles[i].topLeft.y <= expandedRect.topLeft.y && expandedRect.topLeft.y < rectangles[i].bottomRight.y) flagKasanari = 1;
      if (expandedRect.topLeft.y <= rectangles[i].topLeft.y && rectangles[i].topLeft.y < expandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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

  int shuf[4] = {};
  int shuffleSeed = xorshift() % 24;
  for (int j = 0; j < (4); ++j) shuf[j] = shuffles[shuffleSeed][j];

  int yure = xorshift() % 2;

  for (int i = 0; i < (4); ++i) {
    int S = (expandedRect.bottomRight.x - expandedRect.topLeft.x) * (expandedRect.bottomRight.y - expandedRect.topLeft.y);
    if (S <= targetSizes[ite]) { break; }
    if (shuf[i] == 0) {
      int ma = targetSizes[ite] / (expandedRect.bottomRight.y - expandedRect.topLeft.y) + yure;
      int diff = (expandedRect.bottomRight.x - expandedRect.topLeft.x) - ma;
      int capa = points[ite].x - expandedRect.topLeft.x;
      if (capa >= diff) {
        expandedRect.topLeft.x += diff;
      }
      else {
        expandedRect.topLeft.x += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = targetSizes[ite] / (expandedRect.bottomRight.x - expandedRect.topLeft.x) + yure;
      int diff = (expandedRect.bottomRight.y - expandedRect.topLeft.y) - ma;
      int capa = points[ite].y - expandedRect.topLeft.y;
      if (capa >= diff) {
        expandedRect.topLeft.y += diff;
      }
      else {
        expandedRect.topLeft.y += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = targetSizes[ite] / (expandedRect.bottomRight.y - expandedRect.topLeft.y) + yure;
      int diff = (expandedRect.bottomRight.x - expandedRect.topLeft.x) - ma;
      int capa = expandedRect.bottomRight.x - (points[ite].x + 1);
      if (capa >= diff) {
        expandedRect.bottomRight.x -= diff;
      }
      else {
        expandedRect.bottomRight.x -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = targetSizes[ite] / (expandedRect.bottomRight.x - expandedRect.topLeft.x) + yure;
      int diff = (expandedRect.bottomRight.y - expandedRect.topLeft.y) - ma;
      int capa = expandedRect.bottomRight.y - (points[ite].y + 1);
      if (capa >= diff) {
        expandedRect.bottomRight.y -= diff;
      }
      else {
        expandedRect.bottomRight.y -= capa;
      }
    }
  }

  int ng = 0;
  if (expandedRect.topLeft.x < 0 || 10000 < expandedRect.topLeft.x) ng = 1;
  if (expandedRect.topLeft.y < 0 || 10000 < expandedRect.topLeft.y) ng = 1;
  if (expandedRect.bottomRight.x < 0 || 10000 < expandedRect.bottomRight.x) ng = 1;
  if (expandedRect.bottomRight.y < 0 || 10000 < expandedRect.bottomRight.y) ng = 1;
  if (expandedRect.bottomRight.x <= expandedRect.topLeft.x) ng = 1;
  if (expandedRect.bottomRight.y <= expandedRect.topLeft.y) ng = 1;
  if (points[ite].x < expandedRect.topLeft.x || expandedRect.bottomRight.x <= points[ite].x) ng = 1;
  if (points[ite].y < expandedRect.topLeft.y || expandedRect.bottomRight.y <= points[ite].y) ng = 1;

  if (ng) {
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
  Rect keep = rectangles[ite];

  expandRect(ite);
  rectangles[ite] = expandedRect;

  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  const double prob = exp((double)tmp / temp);

  if (prob > rand01()) {
    modeCount[4]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite] = keep;
    calcScore(ite);
  }
}

Rect best_best_rects[MAX_N];
int real_real_maxScore = -1;

inline void resetRects(int n_crush)
{
  // n_crush個つぶす
  for (int i = 0; i < (n_crush); ++i) {
    int ite = xorshift() % numRects;
    initRect(rectangles[ite], points[ite]);
  }

  currentScore = calcScore(-1);
  if (currentScore > bestScore) {
    saveBest();
  }
}

inline void resetWorstRects(int n_worst)
{
  vector<pair<double, int>> v;
  for (int i = 0; i < numRects; ++i) {
    v.emplace_back(pair<double, int>(rectScores[i], i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < (n_worst); ++i) {
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
  vector<int> keep;
  keep.emplace_back(ite);
  rectangles[ite].topLeft.x -= hole;
  rectangles[ite].topLeft.y -= hole;
  rectangles[ite].bottomRight.x += hole;
  rectangles[ite].bottomRight.y += hole;
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) keep.emplace_back(i);
  }
  int keepSize = keep.size();
  for (int i = 0; i < (keepSize); ++i) {
    initRect(rectangles[keep[i]], points[keep[i]]);
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

  int tateyoko = xorshift() % 2;

  if (tateyoko == 0) {
    int argX = indexInSortedX[ite];
    for (int ii = argX - 1; ii >= 0; --ii) {
      int i = sortedByX[ii];
      if (points[i].x == points[ite].x) { continue; }
      int flagKasanari = 0;
      if (points[i].y <= largeExpandedRect.topLeft.y && largeExpandedRect.topLeft.y < points[i].y + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.y <= points[i].y && points[i].y < largeExpandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].y <= largeExpandedRect.topLeft.y && largeExpandedRect.topLeft.y < points[i].y + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.y <= points[i].y && points[i].y < largeExpandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].x <= largeExpandedRect.topLeft.x && largeExpandedRect.topLeft.x < points[i].x + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.x <= points[i].x && points[i].x < largeExpandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].x <= largeExpandedRect.topLeft.x && largeExpandedRect.topLeft.x < points[i].x + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.x <= points[i].x && points[i].x < largeExpandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].x <= largeExpandedRect.topLeft.x && largeExpandedRect.topLeft.x < points[i].x + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.x <= points[i].x && points[i].x < largeExpandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].x <= largeExpandedRect.topLeft.x && largeExpandedRect.topLeft.x < points[i].x + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.x <= points[i].x && points[i].x < largeExpandedRect.bottomRight.x) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].y <= largeExpandedRect.topLeft.y && largeExpandedRect.topLeft.y < points[i].y + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.y <= points[i].y && points[i].y < largeExpandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
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
      int flagKasanari = 0;
      if (points[i].y <= largeExpandedRect.topLeft.y && largeExpandedRect.topLeft.y < points[i].y + 1) flagKasanari = 1;
      if (largeExpandedRect.topLeft.y <= points[i].y && points[i].y < largeExpandedRect.bottomRight.y) flagKasanari = 1;
      if (flagKasanari) {
        if (points[i].x <= points[ite].x) {
          largeExpandedRect.topLeft.x = max(largeExpandedRect.topLeft.x, points[i].x + 1);
        }
        else {
          largeExpandedRect.bottomRight.x = min(largeExpandedRect.bottomRight.x, points[i].x);
        }
      }
    }
  }

  int shuf[4] = {};
  int shuffleSeed = xorshift() % 24;
  for (int j = 0; j < (4); ++j) shuf[j] = shuffles[shuffleSeed][j];

  int yure = xorshift() % 2;

  for (int i = 0; i < (4); ++i) {
    int S = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) * (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y);
    if (S <= targetSizes[ite]) { break; }
    if (shuf[i] == 0) {
      int ma = targetSizes[ite] / (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) + yure;
      int diff = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) - ma;
      if (diff < 0) diff = 0;
      int capa = points[ite].x - largeExpandedRect.topLeft.x;
      if (capa >= diff) {
        largeExpandedRect.topLeft.x += diff;
      }
      else {
        largeExpandedRect.topLeft.x += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = targetSizes[ite] / (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) + yure;
      int diff = (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) - ma;
      if (diff < 0) diff = 0;
      int capa = points[ite].y - largeExpandedRect.topLeft.y;
      if (capa >= diff) {
        largeExpandedRect.topLeft.y += diff;
      }
      else {
        largeExpandedRect.topLeft.y += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = targetSizes[ite] / (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) + yure;
      int diff = (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) - ma;
      if (diff < 0) diff = 0;
      int capa = largeExpandedRect.bottomRight.x - (points[ite].x + 1);
      if (capa >= diff) {
        largeExpandedRect.bottomRight.x -= diff;
      }
      else {
        largeExpandedRect.bottomRight.x -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = targetSizes[ite] / (largeExpandedRect.bottomRight.x - largeExpandedRect.topLeft.x) + yure;
      int diff = (largeExpandedRect.bottomRight.y - largeExpandedRect.topLeft.y) - ma;
      if (diff < 0) diff = 0;
      int capa = largeExpandedRect.bottomRight.y - (points[ite].y + 1);
      if (capa >= diff) {
        largeExpandedRect.bottomRight.y -= diff;
      }
      else {
        largeExpandedRect.bottomRight.y -= capa;
      }
    }
  }

  int ng = 0;
  if (largeExpandedRect.topLeft.x < 0 || 10000 < largeExpandedRect.topLeft.x) ng = 1;
  if (largeExpandedRect.topLeft.y < 0 || 10000 < largeExpandedRect.topLeft.y) ng = 1;
  if (largeExpandedRect.bottomRight.x < 0 || 10000 < largeExpandedRect.bottomRight.x) ng = 1;
  if (largeExpandedRect.bottomRight.y < 0 || 10000 < largeExpandedRect.bottomRight.y) ng = 1;
  if (largeExpandedRect.bottomRight.x <= largeExpandedRect.topLeft.x) ng = 1;
  if (largeExpandedRect.bottomRight.y <= largeExpandedRect.topLeft.y) ng = 1;
  if (points[ite].x < largeExpandedRect.topLeft.x || largeExpandedRect.bottomRight.x <= points[ite].x) ng = 1;
  if (points[ite].y < largeExpandedRect.topLeft.y || largeExpandedRect.bottomRight.y <= points[ite].y) ng = 1;

  if (ng) {
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
  int abcd = xorshift() % 4;

  if (abcd == 0) rectangles[ite].topLeft.x += diff;
  if (abcd == 1) rectangles[ite].topLeft.y += diff;
  if (abcd == 2) rectangles[ite].bottomRight.x += diff;
  if (abcd == 3) rectangles[ite].bottomRight.y += diff;

  if (isValid(ite) == 0) {
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  const double prob = exp((double)tmp / temp);
  if (prob > rand01()) {
    modeCount[0]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
    calcScore(ite);
  }
}

inline void changeAllEdges(int ite, double temp)
{
  int diffA = xorshift() % 101 - 50;
  int diffB = xorshift() % 101 - 50;
  int diffC = xorshift() % 101 - 50;
  int diffD = xorshift() % 101 - 50;

  rectangles[ite].topLeft.x += diffA;
  rectangles[ite].topLeft.y += diffB;
  rectangles[ite].bottomRight.x += diffC;
  rectangles[ite].bottomRight.y += diffD;

  if (isValid(ite) == 0) {
    rectangles[ite].topLeft.x -= diffA;
    rectangles[ite].topLeft.y -= diffB;
    rectangles[ite].bottomRight.x -= diffC;
    rectangles[ite].bottomRight.y -= diffD;
    return;
  }

  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  const double prob = exp((double)tmp / temp);
  if (prob > rand01()) {
    modeCount[3]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite].topLeft.x -= diffA;
    rectangles[ite].topLeft.y -= diffB;
    rectangles[ite].bottomRight.x -= diffC;
    rectangles[ite].bottomRight.y -= diffD;
    calcScore(ite);
  }
}

inline void slideRect(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 101 - 50;
  int ab = xorshift() % 2;

  if (ab == 0) {
    rectangles[ite].topLeft.x += diff;
    rectangles[ite].bottomRight.x += diff;
  }
  if (ab == 1) {
    rectangles[ite].topLeft.y += diff;
    rectangles[ite].bottomRight.y += diff;
  }

  if (isValid(ite) == 0) {
    if (ab == 0) {
      rectangles[ite].topLeft.x -= diff;
      rectangles[ite].bottomRight.x -= diff;
    }
    if (ab == 1) {
      rectangles[ite].topLeft.y -= diff;
      rectangles[ite].bottomRight.y -= diff;
    }
    return;
  }

  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  const double prob = exp((double)tmp / temp);
  if (tmpScore >= currentScore) {
    modeCount[1]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    if (ab == 0) {
      rectangles[ite].topLeft.x -= diff;
      rectangles[ite].bottomRight.x -= diff;
    }
    if (ab == 1) {
      rectangles[ite].topLeft.y -= diff;
      rectangles[ite].bottomRight.y -= diff;
    }
    calcScore(ite);
  }
}

inline void changeAspectRatio(int ite, double temp)
{
  int yokoRatio = xorshift() % 9 + 1; // 1 ~ 9;
  int tateRatio = 10 - yokoRatio;

  int S = yokoRatio * tateRatio;
  int mul = sqrt(targetSizes[ite] / S);
  if (mul == 0) { return; }

  int yoko = yokoRatio * mul;
  int tate = tateRatio * mul;

  Rect keep = rectangles[ite];

  int leftA = max(0, points[ite].x - (yoko - 1));
  int rightA = min(points[ite].x, 10000 - yoko);
  int rangeA = rightA - leftA + 1;
  if (rangeA < 1) { return; }

  int leftB = max(0, points[ite].y - (tate - 1));
  int rightB = min(points[ite].y, 10000 - tate);
  int rangeB = rightB - leftB + 1;
  if (rangeB < 1) { return; }

  rectangles[ite].topLeft.x = xorshift() % rangeA + leftA;
  rectangles[ite].bottomRight.x = rectangles[ite].topLeft.x + rangeA;
  rectangles[ite].topLeft.y = xorshift() % rangeB + leftB;
  rectangles[ite].bottomRight.y = rectangles[ite].topLeft.y + rangeB;

  if (isValid(ite) == 0) {
    rectangles[ite] = keep;
    return;
  }

  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  const double prob = exp((double)tmp / temp);
  if (tmpScore >= currentScore) {
    modeCount[2]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    // 元に戻す
    rectangles[ite] = keep;
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

inline int canShiftBoundary(int ite, int abcd)
{
  for (int i = 0; i < numRects; ++i) {
    if (i == ite) { continue; }
    if (checkOverlap(i, ite)) {
      if (abcd == 0) rectangles[i].bottomRight.x = rectangles[ite].topLeft.x;
      if (abcd == 1) rectangles[i].bottomRight.y = rectangles[ite].topLeft.y;
      if (abcd == 2) rectangles[i].topLeft.x = rectangles[ite].bottomRight.x;
      if (abcd == 3) rectangles[i].topLeft.y = rectangles[ite].bottomRight.y;

      if (isSelfInvalid(i)) return 0;
    }
  }
  return 1;
}

int overlappingRects[MAX_N];
int overlapCount;
inline void findOverlaps(int ite, int abcd)
{
  overlapCount = 0;
  if (abcd == 0) {
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
  if (abcd == 1) {
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
  if (abcd == 2) {
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
  if (abcd == 3) {
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

int keepvA[MAX_N], keepvB[MAX_N], keepvC[MAX_N], keepvD[MAX_N];
inline void shiftBoundary(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = xorshift() % 50 + 1;
  int abcd = xorshift() % 4;

  if (abcd < 2) diff *= -1;

  if (abcd == 0) rectangles[ite].topLeft.x += diff;
  if (abcd == 1) rectangles[ite].topLeft.y += diff;
  if (abcd == 2) rectangles[ite].bottomRight.x += diff;
  if (abcd == 3) rectangles[ite].bottomRight.y += diff;

  if (isSelfInvalid(ite)) {
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  findOverlaps(ite, abcd);
  int vn = overlapCount;

  if (vn > 0 && overlappingRects[0] == -1) {
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }


  for (int i = 0; i < (vn); ++i) {
    keepvA[i] = rectangles[overlappingRects[i]].topLeft.x;
    keepvB[i] = rectangles[overlappingRects[i]].topLeft.y;
    keepvC[i] = rectangles[overlappingRects[i]].bottomRight.x;
    keepvD[i] = rectangles[overlappingRects[i]].bottomRight.y;
  }

  int ok = 1;
  for (int i = 0; i < (vn); ++i) {
    if (abcd == 0) rectangles[overlappingRects[i]].bottomRight.x = rectangles[ite].topLeft.x;
    if (abcd == 1) rectangles[overlappingRects[i]].bottomRight.y = rectangles[ite].topLeft.y;
    if (abcd == 2) rectangles[overlappingRects[i]].topLeft.x = rectangles[ite].bottomRight.x;
    if (abcd == 3) rectangles[overlappingRects[i]].topLeft.y = rectangles[ite].bottomRight.y;
    if (isSelfInvalid(overlappingRects[i])) ok = 0;
  }

  if (ok == 0) {
    for (int i = 0; i < (vn); ++i) {
      rectangles[overlappingRects[i]].topLeft.x = keepvA[i];
      rectangles[overlappingRects[i]].topLeft.y = keepvB[i];
      rectangles[overlappingRects[i]].bottomRight.x = keepvC[i];
      rectangles[overlappingRects[i]].bottomRight.y = keepvD[i];
    }
    // 元に戻す
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
    return;
  }

  for (int i = 0; i < (vn); ++i) calcScore(overlappingRects[i]);
  int tmpScore = calcScore(ite);

  int tmp = tmpScore - currentScore;
  double prob = exp((double)tmp / temp);

  if (prob > rand01()) {
    modeCount[5]++;
    currentScore = tmpScore;
    if (currentScore > bestScore) {
      saveBest();
    }
  }
  else {
    for (int i = 0; i < (vn); ++i) {
      rectangles[overlappingRects[i]].topLeft.x = keepvA[i];
      rectangles[overlappingRects[i]].topLeft.y = keepvB[i];
      rectangles[overlappingRects[i]].bottomRight.x = keepvC[i];
      rectangles[overlappingRects[i]].bottomRight.y = keepvD[i];
      calcScore(overlappingRects[i]);
    }
    // 元に戻す
    if (abcd == 0) rectangles[ite].topLeft.x -= diff;
    if (abcd == 1) rectangles[ite].topLeft.y -= diff;
    if (abcd == 2) rectangles[ite].bottomRight.x -= diff;
    if (abcd == 3) rectangles[ite].bottomRight.y -= diff;
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

Rect best_best_best_rects[MAX_N];
int real_real_real_maxScore = -1;

int a2[100][MAX_N], b2[100][MAX_N], c2[100][MAX_N], d2[100][MAX_N];
int a4[100][MAX_N], b4[100][MAX_N], c4[100][MAX_N], d4[100][MAX_N];
int maxScore4[100] = {};

int ui_tei_a[MAX_N], ui_tei_b[MAX_N], ui_tei_c[MAX_N], ui_tei_d[MAX_N];
int ui_tei_maxScore = -1;

inline void multiStartSearch()
{
  clock_t start, end;
  for (int ui_tei = 0; ui_tei < (5); ++ui_tei) {

    // 初期解
    // 左上(x,y)、右下(x+1,y+1)
    for (int i = 0; i < numRects; ++i) {
      initRect(rectangles[i], points[i]);
    }

    int T = 5;
    for (int _ = 0; _ < (T); ++_) {
      start = clock();

      // 初期スコア計算
      currentScore = calcScore(-1);
      saveBest();

      // 焼きなまし
      start = clock();
      end = clock();
      double now_time = ((double)end - start) / CLOCKS_PER_SEC;
      double TL = (0.10 / (double)T) / allLoopTimes;
      double start_temp = 2048;
      double end_temp = 0.1;
      double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      int loop = 0;
      while (true) {
        loop++;
        if (loop % 100 == 1) {
          end = clock();
          now_time = ((double)end - start) / CLOCKS_PER_SEC;
          if (now_time > TL)break;
          temp = start_temp + (end_temp - start_temp) * now_time / TL;
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
        else if (now_time > 2.0 / T && mode == 2) {
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

      if (currentScore > real_real_maxScore) {
        real_real_maxScore = currentScore;
        for (int i = 0; i < numRects; ++i) {
          best_best_rects[i] = rectangles[i];
        }
      }
    }

    // real_real_maxScore戻す
    currentScore = real_real_maxScore;
    real_real_maxScore = 0;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = best_best_rects[i];
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
    double now_time = ((double)end - start) / CLOCKS_PER_SEC;
    double TL = 0.02 / allLoopTimes;
    double start_temp = 50048;
    double end_temp = 0.1;
    double temp = start_temp + (end_temp - start_temp) * now_time / TL;
    int loop = 0;
    int kouhan = ui_tei % 2;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        end = clock();
        now_time = ((double)end - start) / CLOCKS_PER_SEC;
        if (now_time > TL)break;
        temp = start_temp + (end_temp - start_temp) * now_time / TL;
      }

      int mode = loop % 5;
      if (mode == 0) { // a,b,c,dのうち1つ変更
        int ite = xorshift() % numRects;
        changeSingleEdge(ite, temp);
      }
      else if (kouhan && mode == 1) { // 位置をスライド
        int ite = xorshift() % numRects;
        slideRect(ite, temp);
      }
      else if (mode == 2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
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

    if (currentScore > ui_tei_maxScore) {
      ui_tei_maxScore = currentScore;
      for (int i = 0; i < numRects; ++i) {
        ui_tei_a[i] = rectangles[i].topLeft.x;
        ui_tei_b[i] = rectangles[i].topLeft.y;
        ui_tei_c[i] = rectangles[i].bottomRight.x;
        ui_tei_d[i] = rectangles[i].bottomRight.y;
      }
    }
  }

  // 元に戻しておく
  currentScore = 0;
  bestScore = 0;
  real_real_maxScore = 0;
  for (int i = 0; i < numRects; ++i) {
    initRect(rectangles[i], points[i]);
    bestRects[i] = rectangles[i];
    best_best_rects[i] = rectangles[i];
  }
}

int solve(int teisyutu, int fileNum)
{
  auto startClock = system_clock::now();
  clock_t start, end;
  clock_t real_start = clock();

  readInput(fileNum);

  initSortArrays();


  for (int allLoop = 0; allLoop < (allLoopTimes); ++allLoop) {

    multiStartSearch();

    // ui_tei_maxScore戻す
    currentScore = ui_tei_maxScore;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i].topLeft.x = ui_tei_a[i];
      rectangles[i].topLeft.y = ui_tei_b[i];
      rectangles[i].bottomRight.x = ui_tei_c[i];
      rectangles[i].bottomRight.y = ui_tei_d[i];
    }
    calcScore(-1);

    // 初期スコア計算
    currentScore = calcScore(-1);
    bestScore = currentScore;

    for (int i = 0; i < numRects; ++i) {
      bestRects[i] = rectangles[i];
    }


    int oya = 1;

    for (int asai = 0; asai < (oya); ++asai) {
      for (int j = 0; j < numRects; ++j) {
        a2[asai][j] = rectangles[j].topLeft.x;
        b2[asai][j] = rectangles[j].topLeft.y;
        c2[asai][j] = rectangles[j].bottomRight.x;
        d2[asai][j] = rectangles[j].bottomRight.y;
      }
    }

    // 焼きなまし2
    // 筋のいいやつを追う
    int T = 250 / allLoopTimes;
    for (int _ = 0; _ < (T); ++_) {
      for (int i = 0; i < (6); ++i) modeCount[i] = 0;

      int TT = 1;
      for (int asai = 0; asai < (TT); ++asai) {
        int idx = asai % oya;
        for (int i = 0; i < numRects; ++i) {
          rectangles[i].topLeft.x = a2[idx][i];
          rectangles[i].topLeft.y = b2[idx][i];
          rectangles[i].bottomRight.x = c2[idx][i];
          rectangles[i].bottomRight.y = d2[idx][i];
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
        double now_time = ((double)end - start) / CLOCKS_PER_SEC;
        double TL = (((timeLimit - 0.7) / T) / TT) / allLoopTimes;
        double start_temp = 20048.0;
        double end_temp = 0.1;
        double temp = start_temp + (end_temp - start_temp) * now_time / TL;
        int loop = 0;
        int kouhan = _ % 2;

        while (true) {
          loop++;
          if (loop % 100 == 1) {
            const double time = duration_cast<microseconds>(system_clock::now() - startClock).count() * 1e-6;
            if (time > TL) { break; }
            const double progressRatio = time / TL;   // 進捗。開始時が0.0、終了時が1.0
            temp = start_temp + (end_temp - start_temp) * progressRatio;
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
          else if (mode == -2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
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
        if (currentScore > real_real_maxScore) {
          real_real_maxScore = currentScore;
          for (int i = 0; i < numRects; ++i) {
            best_best_rects[i] = rectangles[i];
          }
        }

        // ビームサーチの次の種にする
        maxScore4[asai] = currentScore;
        for (int i = 0; i < numRects; ++i) {
          a4[asai][i] = rectangles[i].topLeft.x;
          b4[asai][i] = rectangles[i].topLeft.y;
          c4[asai][i] = rectangles[i].bottomRight.x;
          d4[asai][i] = rectangles[i].bottomRight.y;
        }
      }

      // 次の世代に継承
      vector<P> vBeam;
      for (int asai = 0; asai < (TT); ++asai) vBeam.emplace_back(P(maxScore4[asai], asai));
      sort(vBeam.begin(), vBeam.end(), greater<P>());


      for (int ii = 0; ii < (oya); ++ii) {
        int i = vBeam[ii].second;
        for (int j = 0; j < numRects; ++j) {
          a2[ii][j] = a4[i][j];
          b2[ii][j] = b4[i][j];
          c2[ii][j] = c4[i][j];
          d2[ii][j] = d4[i][j];
        }
      }

      // 提出時以下は消す
      if (teisyutu == 0 && _ % 10 == 0) {
        cout << "_ = " << _;
        cout << ", vBeam[0] = (" << vBeam[0].first << ", " << vBeam[0].second << ")" << endl;
      }

      // エスケープ
      end = clock();
      if (((double)end - real_start) / CLOCKS_PER_SEC > timeLimit) { break; }
    }

    // real_real_maxScore戻す
    currentScore = real_real_maxScore;
    for (int i = 0; i < numRects; ++i) {
      rectangles[i] = best_best_rects[i];
    }
    calcScore(-1);

    if (teisyutu == 0) {
      cout << "currentScore = " << currentScore << endl;
    }

    const int MOD = 1000000007;
    if (teisyutu == 0 && currentScore > MOD) {
      cout << "ERROR" << endl;
      writeErrorLog(fileNum);
    }

    // real_real_real入れる
    if (currentScore > real_real_real_maxScore && currentScore < 1000000007) {
      real_real_real_maxScore = currentScore;
      for (int i = 0; i < numRects; ++i) {
        best_best_best_rects[i] = rectangles[i];
      }
    }


    // すべて白紙にリセットする
    currentScore = 0;
    bestScore = 0;
    real_real_maxScore = 0;
    for (int i = 0; i < numRects; ++i) {
      initRect(rectangles[i], points[i]);
      bestRects[i] = rectangles[i];
      best_best_rects[i] = rectangles[i];
    }
  }

  // real_real_real_maxScore戻す
  currentScore = real_real_real_maxScore;
  for (int i = 0; i < numRects; ++i) {
    rectangles[i] = best_best_best_rects[i];
  }
  calcScore(-1);

  // 最終出力
  if (teisyutu) {
    for (int i = 0; i < numRects; ++i) {
      cout << rectangles[i].topLeft.x << ' ' << rectangles[i].topLeft.y << ' ' << rectangles[i].bottomRight.x << ' ' << rectangles[i].bottomRight.y << endl;
    }
  }
  else {
    writeOutput(fileNum);
  }

  // 提出時以下は消す
  if (teisyutu == 0) {
    cout << "file No. = " << fileNum << ", currentScore = " << currentScore << endl;
  }

  if (teisyutu == 0 && currentScore > 1000000007) {
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
  real_real_maxScore = -1;
  ui_tei_maxScore = -1;
  real_real_real_maxScore = -1;
  for (int i = 0; i < (MAX_N); ++i) {
    points[i].x = 0, points[i].y = 0, targetSizes[i] = 0;
    rectangles[i].topLeft.x = 0, rectangles[i].topLeft.y = 0, rectangles[i].bottomRight.x = 0, rectangles[i].bottomRight.y = 0;
    clearRect(rectangles[i]);
    rectAreas[i] = 0;
    clearRect(bestRects[i]);
    rectScores[i] = 0;
    sortedByX[i] = 0, sortedByY[i] = 0;
    indexInSortedX[i] = 0, indexInSortedY[i] = 0;
    clearRect(best_best_rects[i]);
    ui_tei_a[i] = 0, ui_tei_b[i] = 0, ui_tei_c[i] = 0, ui_tei_d[i] = 0;
    clearRect(best_best_best_rects[i]);
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
