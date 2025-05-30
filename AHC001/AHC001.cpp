#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
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

#define srep(i,s,t) for (int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
using namespace chrono;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 205

static uint32_t Rand() {
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
static double Rand01() {
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}

// ハイパラはここにおく
int haipara_hoge;
int haipara_oya = 1;
int haipara_TT = 1;
int modeCount[20];

int shuffles[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

/*
いろいろ
int x = Rand()%100;
double x = Rand01();
mt19937 get_rand_mt;
shuffle(vec.begin(), vec.end(), get_rand_mt);
int dx[4] = {-1, 0, 1, 0};
int dy[4] = {0, -1, 0, 1};
char cc[4] = {'U','L','D','R'};
*/

double realTL = 4.8;

struct Point {
  int x;
  int y;
};

struct Rect {
  Point p1;
  Point p2;
};

int allLoopTimes = 1;
int n;
Point target_points[MAX_N];
int target_sizes[MAX_N];
Rect rects[MAX_N];
int area_sizes[MAX_N];
Rect best_rects[MAX_N];

inline void calc_area(int idx) {
  area_sizes[idx] = (rects[idx].p2.x - rects[idx].p1.x) * (rects[idx].p2.y - rects[idx].p1.y);
}

inline void nyuuryokuInit(int fileNum) {
  // 入力
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << fileNum << ".txt";
  ifstream ifs(oss.str());
  if (!ifs.is_open()) { // 標準入力する
    cin >> n;
    for (int i = 0; i < (n); ++i) cin >> target_points[i].x >> target_points[i].y >> target_sizes[i];
  }
  else { // ファイル入力する
    ifs >> n;
    for (int i = 0; i < (n); ++i) ifs >> target_points[i].x >> target_points[i].y >> target_sizes[i];
  }
}

inline void FileKakikomi(int fileNum) {
  string fileName = to_string(fileNum);
  fileName += "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);
  for (int i = 0; i < (n); ++i) {
    ofs << rects[i].p1.x << ' ' << rects[i].p1.y << ' ' << rects[i].p2.x << ' ' << rects[i].p2.y << endl;
  }
  ofs.close();
}

inline void FileKakikomiERROR(int fileNum) {
  string fileName = to_string(fileNum);
  fileName += "_out_ERROR.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);
  for (int i = 0; i < (n); ++i) {
    ofs << rects[i].p1.x << ' ' << rects[i].p1.y << ' ' << rects[i].p2.x << ' ' << rects[i].p2.y << endl;
  }
  ofs.close();
}

int maxScore = -1;
int real_maxScore = -1;

double p[MAX_N];
double pSum;
inline int calc(int ite) {
  if (ite == -1) {
    double sum = 0;
    for (int i = 0; i < (n); ++i) {
      calc_area(i);
      p[i] = 1.0 - (1.0 - (double)min(target_sizes[i], area_sizes[i]) / (double)max(target_sizes[i], area_sizes[i])) * (1.0 - (double)min(target_sizes[i], area_sizes[i]) / (double)max(target_sizes[i], area_sizes[i]));
      sum += p[i];
    }
    pSum = sum;
    sum /= (double)n;
    sum *= 1000000000.0;
    return round(sum);
  }
  else {
    double sum = pSum;
    sum -= p[ite];
    calc_area(ite);
    p[ite] = 1.0 - (1.0 - (double)min(target_sizes[ite], area_sizes[ite]) / (double)max(target_sizes[ite], area_sizes[ite])) * (1.0 - (double)min(target_sizes[ite], area_sizes[ite]) / (double)max(target_sizes[ite], area_sizes[ite]));
    sum += p[ite];
    pSum = sum;
    sum /= (double)n;
    sum *= 1000000000.0;
    return round(sum);
  }
}

inline int kasanarihantei(int i, int j) {
  int cnt = 0;
  if (rects[i].p1.x <= rects[j].p1.x && rects[j].p1.x < rects[i].p2.x) cnt++;
  else if (rects[j].p1.x <= rects[i].p1.x && rects[i].p1.x < rects[j].p2.x) cnt++;
  if (rects[i].p1.y <= rects[j].p1.y && rects[j].p1.y < rects[i].p2.y) cnt++;
  else if (rects[j].p1.y <= rects[i].p1.y && rects[i].p1.y < rects[j].p2.y) cnt++;
  return cnt == 2;
}

int sort_x[MAX_N], sort_y[MAX_N];
int arg_sort_x[MAX_N], arg_sort_y[MAX_N];
inline void sortInit() {
  vector<P> v;
  for (int i = 0; i < (n); ++i) {
    v.emplace_back(P(target_points[i].x, i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < (n); ++i) {
    sort_x[i] = v[i].second;
    arg_sort_x[v[i].second] = i;
  }

  v.clear();
  for (int i = 0; i < (n); ++i) {
    v.emplace_back(P(target_points[i].y, i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < (n); ++i) {
    sort_y[i] = v[i].second;
    arg_sort_y[v[i].second] = i;
  }
}

// 0~10000を出ていないか
// 面積は1以上か
// (x[i]+0.5,y[i]+0.5)を含んでいるか
// 重なりがないか
inline int isOK2(int ite) {
  if (ite == -1) {
    for (int i = 0; i < (n); ++i) {
      if (rects[i].p1.x < 0 || 10000 < rects[i].p1.x) return 0;
      if (rects[i].p1.y < 0 || 10000 < rects[i].p1.y) return 0;
      if (rects[i].p2.x < 0 || 10000 < rects[i].p2.x) return 0;
      if (rects[i].p2.y < 0 || 10000 < rects[i].p2.y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (rects[i].p2.x <= rects[i].p1.x) return 0;
      if (rects[i].p2.y <= rects[i].p1.y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (target_points[i].x < rects[i].p1.x || rects[i].p2.x <= target_points[i].x) return 0;
      if (target_points[i].y < rects[i].p1.y || rects[i].p2.y <= target_points[i].y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      srep(j, i + 1, n) {
        if (kasanarihantei(i, j)) return 0;
      }
    }
  }
  else {
    if (rects[ite].p1.x < 0 || 10000 < rects[ite].p1.x) return 0;
    if (rects[ite].p1.y < 0 || 10000 < rects[ite].p1.y) return 0;
    if (rects[ite].p2.x < 0 || 10000 < rects[ite].p2.x) return 0;
    if (rects[ite].p2.y < 0 || 10000 < rects[ite].p2.y) return 0;
    if (rects[ite].p2.x <= rects[ite].p1.x) return 0;
    if (rects[ite].p2.y <= rects[ite].p1.y) return 0;
    if (target_points[ite].x < rects[ite].p1.x || rects[ite].p2.x <= target_points[ite].x) return 0;
    if (target_points[ite].y < rects[ite].p1.y || rects[ite].p2.y <= target_points[ite].y) return 0;
    for (int i = 0; i < (n); ++i) {
      if (i == ite) continue;
      if (kasanarihantei(i, ite)) return 0;
    }
  }
  return 1;
}

inline int isOK(int ite) {
  if (ite == -1) {
    for (int i = 0; i < (n); ++i) {
      if (rects[i].p1.x < 0 || 10000 < rects[i].p1.x) return 0;
      if (rects[i].p1.y < 0 || 10000 < rects[i].p1.y) return 0;
      if (rects[i].p2.x < 0 || 10000 < rects[i].p2.x) return 0;
      if (rects[i].p2.y < 0 || 10000 < rects[i].p2.y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (rects[i].p2.x <= rects[i].p1.x) return 0;
      if (rects[i].p2.y <= rects[i].p1.y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (target_points[i].x < rects[i].p1.x || rects[i].p2.x <= target_points[i].x) return 0;
      if (target_points[i].y < rects[i].p1.y || rects[i].p2.y <= target_points[i].y) return 0;
    }
    for (int i = 0; i < (n); ++i) {
      srep(j, i + 1, n) {
        if (kasanarihantei(i, j)) return 0;
      }
    }
  }
  else {
    if (rects[ite].p1.x < 0 || 10000 < rects[ite].p1.x) return 0;
    if (rects[ite].p1.y < 0 || 10000 < rects[ite].p1.y) return 0;
    if (rects[ite].p2.x < 0 || 10000 < rects[ite].p2.x) return 0;
    if (rects[ite].p2.y < 0 || 10000 < rects[ite].p2.y) return 0;
    if (rects[ite].p2.x <= rects[ite].p1.x) return 0;
    if (rects[ite].p2.y <= rects[ite].p1.y) return 0;
    if (target_points[ite].x < rects[ite].p1.x || rects[ite].p2.x <= target_points[ite].x) return 0;
    if (target_points[ite].y < rects[ite].p1.y || rects[ite].p2.y <= target_points[ite].y) return 0;
    int argX = arg_sort_x[ite];
    int nowLeft = rects[ite].p1.y;
    drep(ii, argX) {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (rects[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.y);
        if (nowLeft >= rects[ite].p2.y) break;
      }
    }
    nowLeft = rects[ite].p1.y;
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (rects[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.y);
        if (nowLeft >= rects[ite].p2.y) break;
      }
    }
  }
  return 1;
}

inline void init_rect(Rect& rect, Point& point) {
  rect.p1.x = point.x;
  rect.p1.y = point.y;
  rect.p2.x = point.x + 1;
  rect.p2.y = point.y + 1;
}

Rect vExtend;
inline void hukuramashi(int ite) {
  init_rect(vExtend, target_points[ite]);

  int flagTateYoko = Rand() % 2;
  if (flagTateYoko == 0) {
    vExtend.p1.x = 0;
    vExtend.p2.x = 10000;
    int argX = arg_sort_x[ite];
    drep(ii, argX) {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rects[i].p1.y <= vExtend.p1.y && vExtend.p1.y < rects[i].p2.y) flagKasanari = 1;
      if (vExtend.p1.y <= rects[i].p1.y && rects[i].p1.y < vExtend.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtend.p1.x = max(vExtend.p1.x, rects[i].p2.x);
        }
        else {
          vExtend.p2.x = min(vExtend.p2.x, rects[i].p1.x);
        }
        break;
      }
    }
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rects[i].p1.y <= vExtend.p1.y && vExtend.p1.y < rects[i].p2.y) flagKasanari = 1;
      if (vExtend.p1.y <= rects[i].p1.y && rects[i].p1.y < vExtend.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtend.p1.x = max(vExtend.p1.x, rects[i].p2.x);
        }
        else {
          vExtend.p2.x = min(vExtend.p2.x, rects[i].p1.x);
        }
        break;
      }
    }

    vExtend.p1.y = 0;
    vExtend.p2.y = 10000;
    int argY = arg_sort_y[ite];
    int nowLeft = vExtend.p1.x;
    drep(ii, argY) {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rects[i].p1.x <= vExtend.p1.x && vExtend.p1.x < rects[i].p2.x) flagKasanari = 1;
      if (vExtend.p1.x <= rects[i].p1.x && rects[i].p1.x < vExtend.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtend.p1.y = max(vExtend.p1.y, rects[i].p2.y);
        }
        else {
          vExtend.p2.y = min(vExtend.p2.y, rects[i].p1.y);
        }
        if (rects[i].p1.x <= nowLeft) {
          nowLeft = max(nowLeft, rects[i].p2.x);
          if (vExtend.p2.x <= nowLeft) break;
        }
      }
    }

    nowLeft = vExtend.p1.x;
    srep(ii, argY + 1, n) {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rects[i].p1.x <= vExtend.p1.x && vExtend.p1.x < rects[i].p2.x) flagKasanari = 1;
      if (vExtend.p1.x <= rects[i].p1.x && rects[i].p1.x < vExtend.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtend.p1.y = max(vExtend.p1.y, rects[i].p2.y);
        }
        else {
          vExtend.p2.y = min(vExtend.p2.y, rects[i].p1.y);
        }
        if (rects[i].p1.x <= nowLeft) {
          nowLeft = max(nowLeft, rects[i].p2.x);
          if (vExtend.p2.x <= nowLeft) break;
        }
      }
    }
  }
  else {
    vExtend.p1.y = 0;
    vExtend.p2.y = 10000;
    int argY = arg_sort_y[ite];
    drep(ii, argY) {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rects[i].p1.x <= vExtend.p1.x && vExtend.p1.x < rects[i].p2.x) flagKasanari = 1;
      if (vExtend.p1.x <= rects[i].p1.x && rects[i].p1.x < vExtend.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtend.p1.y = max(vExtend.p1.y, rects[i].p2.y);
        }
        else {
          vExtend.p2.y = min(vExtend.p2.y, rects[i].p1.y);
        }
        break;
      }
    }
    srep(ii, argY + 1, n) {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rects[i].p1.x <= vExtend.p1.x && vExtend.p1.x < rects[i].p2.x) flagKasanari = 1;
      if (vExtend.p1.x <= rects[i].p1.x && rects[i].p1.x < vExtend.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtend.p1.y = max(vExtend.p1.y, rects[i].p2.y);
        }
        else {
          vExtend.p2.y = min(vExtend.p2.y, rects[i].p1.y);
        }
        break;
      }
    }

    vExtend.p1.x = 0;
    vExtend.p2.x = 10000;
    int argX = arg_sort_x[ite];
    int nowLeft = vExtend.p1.y;
    drep(ii, argX) {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rects[i].p1.y <= vExtend.p1.y && vExtend.p1.y < rects[i].p2.y) flagKasanari = 1;
      if (vExtend.p1.y <= rects[i].p1.y && rects[i].p1.y < vExtend.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtend.p1.x = max(vExtend.p1.x, rects[i].p2.x);
        }
        else {
          vExtend.p2.x = min(vExtend.p2.x, rects[i].p1.x);
        }
        if (rects[i].p1.y <= nowLeft) {
          nowLeft = max(nowLeft, rects[i].p2.y);
          if (vExtend.p2.y <= nowLeft) break;
        }
      }
    }
    nowLeft = vExtend.p1.y;
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rects[i].p1.y <= vExtend.p1.y && vExtend.p1.y < rects[i].p2.y) flagKasanari = 1;
      if (vExtend.p1.y <= rects[i].p1.y && rects[i].p1.y < vExtend.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtend.p1.x = max(vExtend.p1.x, rects[i].p2.x);
        }
        else {
          vExtend.p2.x = min(vExtend.p2.x, rects[i].p1.x);
        }
        if (rects[i].p1.y <= nowLeft) {
          nowLeft = max(nowLeft, rects[i].p2.y);
          if (vExtend.p2.y <= nowLeft) break;
        }
      }
    }
  }

  int shuf[4] = {};
  int shuffleSeed = Rand() % 24;
  for (int j = 0; j < (4); ++j) shuf[j] = shuffles[shuffleSeed][j];

  int yure = Rand() % 2;

  for (int i = 0; i < (4); ++i) {
    int S = (vExtend.p2.x - vExtend.p1.x) * (vExtend.p2.y - vExtend.p1.y);
    if (S <= target_sizes[ite]) break;
    if (shuf[i] == 0) {
      int ma = target_sizes[ite] / (vExtend.p2.y - vExtend.p1.y) + yure;
      int diff = (vExtend.p2.x - vExtend.p1.x) - ma;
      int capa = target_points[ite].x - vExtend.p1.x;
      if (capa >= diff) {
        vExtend.p1.x += diff;
      }
      else {
        vExtend.p1.x += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = target_sizes[ite] / (vExtend.p2.x - vExtend.p1.x) + yure;
      int diff = (vExtend.p2.y - vExtend.p1.y) - ma;
      int capa = target_points[ite].y - vExtend.p1.y;
      if (capa >= diff) {
        vExtend.p1.y += diff;
      }
      else {
        vExtend.p1.y += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = target_sizes[ite] / (vExtend.p2.y - vExtend.p1.y) + yure;
      int diff = (vExtend.p2.x - vExtend.p1.x) - ma;
      int capa = vExtend.p2.x - (target_points[ite].x + 1);
      if (capa >= diff) {
        vExtend.p2.x -= diff;
      }
      else {
        vExtend.p2.x -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = target_sizes[ite] / (vExtend.p2.x - vExtend.p1.x) + yure;
      int diff = (vExtend.p2.y - vExtend.p1.y) - ma;
      int capa = vExtend.p2.y - (target_points[ite].y + 1);
      if (capa >= diff) {
        vExtend.p2.y -= diff;
      }
      else {
        vExtend.p2.y -= capa;
      }
    }
  }

  int ng = 0;
  if (vExtend.p1.x < 0 || 10000 < vExtend.p1.x) ng = 1;
  if (vExtend.p1.y < 0 || 10000 < vExtend.p1.y) ng = 1;
  if (vExtend.p2.x < 0 || 10000 < vExtend.p2.x) ng = 1;
  if (vExtend.p2.y < 0 || 10000 < vExtend.p2.y) ng = 1;
  if (vExtend.p2.x <= vExtend.p1.x) ng = 1;
  if (vExtend.p2.y <= vExtend.p1.y) ng = 1;
  if (target_points[ite].x < vExtend.p1.x || vExtend.p2.x <= target_points[ite].x) ng = 1;
  if (target_points[ite].y < vExtend.p1.y || vExtend.p2.y <= target_points[ite].y) ng = 1;

  if (ng) {
    init_rect(vExtend, target_points[ite]);
  }
}

inline void store_best() {
  real_maxScore = maxScore;
  for (int i = 0; i < (n); ++i) {
    best_rects[i] = rects[i];
  }
}

inline void Extend(int ite, double temp) {
  Rect keep = rects[ite];

  hukuramashi(ite);
  rects[ite] = vExtend;

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);
  //if (tmpScore >= maxScore) {
  if (prob > Rand01()) {
    modeCount[4]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    // 元に戻す
    rects[ite] = keep;
    calc(ite);
  }
}

Rect best_best_rects[MAX_N];
int real_real_maxScore = -1;

int ui_tei_a[MAX_N], ui_tei_b[MAX_N], ui_tei_c[MAX_N], ui_tei_d[MAX_N];
int ui_tei_maxScore = -1;

inline void Tubusu(int tubusu) {
  // tubusu個つぶす
  for (int i = 0; i < (tubusu); ++i) {
    int ite = Rand() % n;
    init_rect(rects[ite], target_points[ite]);
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    store_best();
  }
}

inline void TubusuWorst(int tubusu_worst) {
  vector<pair<double, int>> v;
  for (int i = 0; i < (n); ++i) {
    v.emplace_back(pair<double, int>(p[i], i));
  }
  sort(v.begin(), v.end());
  for (int i = 0; i < (tubusu_worst); ++i) {
    int ite = v[i].second;
    init_rect(rects[ite], target_points[ite]);
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    store_best();
  }
}

inline void AnaWoAkeru(int hole = 100) {
  int ite = Rand() % n;
  vector<int> keep;
  keep.emplace_back(ite);
  rects[ite].p1.x -= hole;
  rects[ite].p1.y -= hole;
  rects[ite].p2.x += hole;
  rects[ite].p2.y += hole;
  for (int i = 0; i < (n); ++i) {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) keep.emplace_back(i);
  }
  int keepSize = keep.size();
  for (int i = 0; i < (keepSize); ++i) {
    init_rect(rects[keep[i]], target_points[keep[i]]);
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    store_best();
  }
}

Rect vExtendKing;
inline void hukuramashiKing(int ite) {
  vExtendKing.p1.x = max(0, (int)(target_points[ite].x - Rand() % 1000));
  vExtendKing.p1.y = max(0, (int)(target_points[ite].y - Rand() % 1000));
  vExtendKing.p2.x = min(10000, (int)(target_points[ite].x + 1 + Rand() % 1000));
  vExtendKing.p2.y = min(10000, (int)(target_points[ite].y + 1 + Rand() % 1000));

  int tateyoko = Rand() % 2;

  if (tateyoko == 0) {
    int argX = arg_sort_x[ite];
    drep(ii, argX) {
      int i = sort_x[ii];
      if (target_points[i].x == target_points[ite].x) continue;
      int flagKasanari = 0;
      if (target_points[i].y <= vExtendKing.p1.y && vExtendKing.p1.y < target_points[i].y + 1) flagKasanari = 1;
      if (vExtendKing.p1.y <= target_points[i].y && target_points[i].y < vExtendKing.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtendKing.p1.x = max(vExtendKing.p1.x, target_points[i].x + 1);
        }
        else {
          vExtendKing.p2.x = min(vExtendKing.p2.x, target_points[i].x);
        }
      }
    }
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      if (target_points[i].x == target_points[ite].x) continue;
      int flagKasanari = 0;
      if (target_points[i].y <= vExtendKing.p1.y && vExtendKing.p1.y < target_points[i].y + 1) flagKasanari = 1;
      if (vExtendKing.p1.y <= target_points[i].y && target_points[i].y < vExtendKing.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtendKing.p1.x = max(vExtendKing.p1.x, target_points[i].x + 1);
        }
        else {
          vExtendKing.p2.x = min(vExtendKing.p2.x, target_points[i].x);
        }
      }
    }

    int argY = arg_sort_y[ite];
    drep(ii, argY) {
      int i = sort_y[ii];
      if (target_points[i].y == target_points[ite].y) continue;
      int flagKasanari = 0;
      if (target_points[i].x <= vExtendKing.p1.x && vExtendKing.p1.x < target_points[i].x + 1) flagKasanari = 1;
      if (vExtendKing.p1.x <= target_points[i].x && target_points[i].x < vExtendKing.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtendKing.p1.y = max(vExtendKing.p1.y, target_points[i].y + 1);
        }
        else {
          vExtendKing.p2.y = min(vExtendKing.p2.y, target_points[i].y);
        }
      }
    }
    srep(ii, argY + 1, n) {
      int i = sort_y[ii];
      if (target_points[i].y == target_points[ite].y) continue;
      int flagKasanari = 0;
      if (target_points[i].x <= vExtendKing.p1.x && vExtendKing.p1.x < target_points[i].x + 1) flagKasanari = 1;
      if (vExtendKing.p1.x <= target_points[i].x && target_points[i].x < vExtendKing.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtendKing.p1.y = max(vExtendKing.p1.y, target_points[i].y + 1);
        }
        else {
          vExtendKing.p2.y = min(vExtendKing.p2.y, target_points[i].y);
        }
      }
    }
  }
  else {
    int argY = arg_sort_y[ite];
    drep(ii, argY) {
      int i = sort_y[ii];
      if (target_points[i].y == target_points[ite].y) continue;
      int flagKasanari = 0;
      if (target_points[i].x <= vExtendKing.p1.x && vExtendKing.p1.x < target_points[i].x + 1) flagKasanari = 1;
      if (vExtendKing.p1.x <= target_points[i].x && target_points[i].x < vExtendKing.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtendKing.p1.y = max(vExtendKing.p1.y, target_points[i].y + 1);
        }
        else {
          vExtendKing.p2.y = min(vExtendKing.p2.y, target_points[i].y);
        }
      }
    }
    srep(ii, argY + 1, n) {
      int i = sort_y[ii];
      if (target_points[i].y == target_points[ite].y) continue;
      int flagKasanari = 0;
      if (target_points[i].x <= vExtendKing.p1.x && vExtendKing.p1.x < target_points[i].x + 1) flagKasanari = 1;
      if (vExtendKing.p1.x <= target_points[i].x && target_points[i].x < vExtendKing.p2.x) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].y <= target_points[ite].y) {
          vExtendKing.p1.y = max(vExtendKing.p1.y, target_points[i].y + 1);
        }
        else {
          vExtendKing.p2.y = min(vExtendKing.p2.y, target_points[i].y);
        }
      }
    }

    int argX = arg_sort_x[ite];
    drep(ii, argX) {
      int i = sort_x[ii];
      if (target_points[i].x == target_points[ite].x) continue;
      int flagKasanari = 0;
      if (target_points[i].y <= vExtendKing.p1.y && vExtendKing.p1.y < target_points[i].y + 1) flagKasanari = 1;
      if (vExtendKing.p1.y <= target_points[i].y && target_points[i].y < vExtendKing.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtendKing.p1.x = max(vExtendKing.p1.x, target_points[i].x + 1);
        }
        else {
          vExtendKing.p2.x = min(vExtendKing.p2.x, target_points[i].x);
        }
      }
    }
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      if (target_points[i].x == target_points[ite].x) continue;
      int flagKasanari = 0;
      if (target_points[i].y <= vExtendKing.p1.y && vExtendKing.p1.y < target_points[i].y + 1) flagKasanari = 1;
      if (vExtendKing.p1.y <= target_points[i].y && target_points[i].y < vExtendKing.p2.y) flagKasanari = 1;
      if (flagKasanari) {
        if (target_points[i].x <= target_points[ite].x) {
          vExtendKing.p1.x = max(vExtendKing.p1.x, target_points[i].x + 1);
        }
        else {
          vExtendKing.p2.x = min(vExtendKing.p2.x, target_points[i].x);
        }
      }
    }
  }



  int shuf[4] = {};
  int shuffleSeed = Rand() % 24;
  for (int j = 0; j < (4); ++j) shuf[j] = shuffles[shuffleSeed][j];

  int yure = Rand() % 2;

  for (int i = 0; i < (4); ++i) {
    int S = (vExtendKing.p2.x - vExtendKing.p1.x) * (vExtendKing.p2.y - vExtendKing.p1.y);
    if (S <= target_sizes[ite]) break;
    if (shuf[i] == 0) {
      int ma = target_sizes[ite] / (vExtendKing.p2.y - vExtendKing.p1.y) + yure;
      int diff = (vExtendKing.p2.x - vExtendKing.p1.x) - ma;
      if (diff < 0) diff = 0;
      int capa = target_points[ite].x - vExtendKing.p1.x;
      if (capa >= diff) {
        vExtendKing.p1.x += diff;
      }
      else {
        vExtendKing.p1.x += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = target_sizes[ite] / (vExtendKing.p2.x - vExtendKing.p1.x) + yure;
      int diff = (vExtendKing.p2.y - vExtendKing.p1.y) - ma;
      if (diff < 0) diff = 0;
      int capa = target_points[ite].y - vExtendKing.p1.y;
      if (capa >= diff) {
        vExtendKing.p1.y += diff;
      }
      else {
        vExtendKing.p1.y += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = target_sizes[ite] / (vExtendKing.p2.y - vExtendKing.p1.y) + yure;
      int diff = (vExtendKing.p2.x - vExtendKing.p1.x) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing.p2.x - (target_points[ite].x + 1);
      if (capa >= diff) {
        vExtendKing.p2.x -= diff;
      }
      else {
        vExtendKing.p2.x -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = target_sizes[ite] / (vExtendKing.p2.x - vExtendKing.p1.x) + yure;
      int diff = (vExtendKing.p2.y - vExtendKing.p1.y) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing.p2.y - (target_points[ite].y + 1);
      if (capa >= diff) {
        vExtendKing.p2.y -= diff;
      }
      else {
        vExtendKing.p2.y -= capa;
      }
    }
  }

  int ng = 0;
  if (vExtendKing.p1.x < 0 || 10000 < vExtendKing.p1.x) ng = 1;
  if (vExtendKing.p1.y < 0 || 10000 < vExtendKing.p1.y) ng = 1;
  if (vExtendKing.p2.x < 0 || 10000 < vExtendKing.p2.x) ng = 1;
  if (vExtendKing.p2.y < 0 || 10000 < vExtendKing.p2.y) ng = 1;
  if (vExtendKing.p2.x <= vExtendKing.p1.x) ng = 1;
  if (vExtendKing.p2.y <= vExtendKing.p1.y) ng = 1;
  if (target_points[ite].x < vExtendKing.p1.x || vExtendKing.p2.x <= target_points[ite].x) ng = 1;
  if (target_points[ite].y < vExtendKing.p1.y || vExtendKing.p2.y <= target_points[ite].y) ng = 1;

  if (ng) {
    init_rect(vExtendKing, target_points[ite]);
  }

}

inline void ExtendKing(int ite) {
  hukuramashiKing(ite);
  rects[ite] = vExtendKing;

  for (int i = 0; i < (n); ++i) {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      init_rect(rects[i], target_points[i]);
    }
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    store_best();
  }
}

inline void oneChange(int ite, double temp) {
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int abcd = Rand() % 4;

  if (abcd == 0) rects[ite].p1.x += diff;
  if (abcd == 1) rects[ite].p1.y += diff;
  if (abcd == 2) rects[ite].p2.x += diff;
  if (abcd == 3) rects[ite].p2.y += diff;

  if (isOK(ite) == 0) {
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);
  //if (tmpScore >= maxScore) {
  if (prob > Rand01()) {
    modeCount[0]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    // 元に戻す
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    calc(ite);
  }
}

inline void fourChange(int ite, double temp) {
  int diffA = Rand() % 101 - 50;
  int diffB = Rand() % 101 - 50;
  int diffC = Rand() % 101 - 50;
  int diffD = Rand() % 101 - 50;

  rects[ite].p1.x += diffA;
  rects[ite].p1.y += diffB;
  rects[ite].p2.x += diffC;
  rects[ite].p2.y += diffD;

  if (isOK(ite) == 0) {
    rects[ite].p1.x -= diffA;
    rects[ite].p1.y -= diffB;
    rects[ite].p2.x -= diffC;
    rects[ite].p2.y -= diffD;
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);
  //if (tmpScore >= maxScore) {
  if (prob > Rand01()) {
    modeCount[3]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    // 元に戻す
    rects[ite].p1.x -= diffA;
    rects[ite].p1.y -= diffB;
    rects[ite].p2.x -= diffC;
    rects[ite].p2.y -= diffD;
    calc(ite);
  }
}

inline void Slide(int ite, double temp) {
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int ab = Rand() % 2;

  if (ab == 0) {
    rects[ite].p1.x += diff;
    rects[ite].p2.x += diff;
  }
  if (ab == 1) {
    rects[ite].p1.y += diff;
    rects[ite].p2.y += diff;
  }

  if (isOK(ite) == 0) {
    if (ab == 0) {
      rects[ite].p1.x -= diff;
      rects[ite].p2.x -= diff;
    }
    if (ab == 1) {
      rects[ite].p1.y -= diff;
      rects[ite].p2.y -= diff;
    }
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);
  if (tmpScore >= maxScore) {
    //if (prob > Rand01()) {
    modeCount[1]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    // 元に戻す
    if (ab == 0) {
      rects[ite].p1.x -= diff;
      rects[ite].p2.x -= diff;
    }
    if (ab == 1) {
      rects[ite].p1.y -= diff;
      rects[ite].p2.y -= diff;
    }
    calc(ite);
  }
}

inline void aspectChange(int ite, double temp) {
  int yokoRatio = Rand() % 9 + 1; // 1 ~ 9;
  int tateRatio = 10 - yokoRatio;

  int S = yokoRatio * tateRatio;
  int mul = sqrt(target_sizes[ite] / S);
  if (mul == 0) return;

  int yoko = yokoRatio * mul;
  int tate = tateRatio * mul;

  Rect keep = rects[ite];

  int leftA = max(0, target_points[ite].x - (yoko - 1));
  int rightA = min(target_points[ite].x, 10000 - yoko);
  int rangeA = rightA - leftA + 1;
  if (rangeA < 1) return;

  int leftB = max(0, target_points[ite].y - (tate - 1));
  int rightB = min(target_points[ite].y, 10000 - tate);
  int rangeB = rightB - leftB + 1;
  if (rangeB < 1) return;

  rects[ite].p1.x = Rand() % rangeA + leftA;
  rects[ite].p2.x = rects[ite].p1.x + rangeA;
  rects[ite].p1.y = Rand() % rangeB + leftB;
  rects[ite].p2.y = rects[ite].p1.y + rangeB;

  if (isOK(ite) == 0) {
    rects[ite] = keep;
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);
  if (tmpScore >= maxScore) {
    //if (prob > Rand01()) {
    modeCount[2]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    // 元に戻す
    rects[ite] = keep;
    calc(ite);
  }
}

inline int selfNg(int ite) {
  if (rects[ite].p1.x < 0 || 10000 < rects[ite].p1.x) return 1;
  if (rects[ite].p1.y < 0 || 10000 < rects[ite].p1.y) return 1;
  if (rects[ite].p2.x < 0 || 10000 < rects[ite].p2.x) return 1;
  if (rects[ite].p2.y < 0 || 10000 < rects[ite].p2.y) return 1;
  if (rects[ite].p2.x <= rects[ite].p1.x) return 1;
  if (rects[ite].p2.y <= rects[ite].p1.y) return 1;
  if (target_points[ite].x < rects[ite].p1.x || rects[ite].p2.x <= target_points[ite].x) return 1;
  if (target_points[ite].y < rects[ite].p1.y || rects[ite].p2.y <= target_points[ite].y) return 1;
  return 0;
}

inline int dokasuOK(int ite, int abcd) {
  for (int i = 0; i < (n); ++i) {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      if (abcd == 0) rects[i].p2.x = rects[ite].p1.x;
      if (abcd == 1) rects[i].p2.y = rects[ite].p1.y;
      if (abcd == 2) rects[i].p1.x = rects[ite].p2.x;
      if (abcd == 3) rects[i].p1.y = rects[ite].p2.y;

      if (selfNg(i)) return 0;
    }
  }
  return 1;
}

int arrKasanari[MAX_N];
int kasanariCount;
inline void kasanaritati(int ite, int abcd) {
  kasanariCount = 0;
  if (abcd == 0) {
    int argX = arg_sort_x[ite];
    int nowLeft = rects[ite].p1.y;
    int nowRight = rects[ite].p2.y;
    drep(ii, argX) {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (rects[ite].p1.x <= target_points[i].x) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rects[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.y);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rects[i].p2.y) {
        nowRight = min(nowRight, rects[i].p1.y);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 1) {
    int argY = arg_sort_y[ite];
    int nowLeft = rects[ite].p1.x;
    int nowRight = rects[ite].p2.x;
    drep(ii, argY) {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (rects[ite].p1.y <= target_points[i].y) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rects[i].p1.x <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.x);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rects[i].p2.x) {
        nowRight = min(nowRight, rects[i].p1.x);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 2) {
    int argX = arg_sort_x[ite];
    int nowLeft = rects[ite].p1.y;
    int nowRight = rects[ite].p2.y;
    srep(ii, argX + 1, n) {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (target_points[i].x < rects[ite].p2.x) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rects[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.y);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rects[i].p2.y) {
        nowRight = min(nowRight, rects[i].p1.y);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 3) {
    int argY = arg_sort_y[ite];
    int nowLeft = rects[ite].p1.x;
    int nowRight = rects[ite].p2.x;
    srep(ii, argY + 1, n) {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (target_points[i].y < rects[ite].p2.y) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rects[i].p1.x <= nowLeft) {
        nowLeft = max(nowLeft, rects[i].p2.x);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rects[i].p2.x) {
        nowRight = min(nowRight, rects[i].p1.x);
        if (nowLeft >= nowRight) break;
      }
    }
  }
}

int keepvA[MAX_N], keepvB[MAX_N], keepvC[MAX_N], keepvD[MAX_N];
inline void zurasi2(int ite, double temp) {
  int diff = 0;
  while (diff == 0) diff = Rand() % 50 + 1;
  int abcd = Rand() % 4;

  if (abcd < 2) diff *= -1;

  if (abcd == 0) rects[ite].p1.x += diff;
  if (abcd == 1) rects[ite].p1.y += diff;
  if (abcd == 2) rects[ite].p2.x += diff;
  if (abcd == 3) rects[ite].p2.y += diff;

  if (selfNg(ite)) {
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    return;
  }

  kasanaritati(ite, abcd);
  int vn = kasanariCount;

  if (vn > 0 && arrKasanari[0] == -1) {
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    return;
  }


  for (int i = 0; i < (vn); ++i) {
    keepvA[i] = rects[arrKasanari[i]].p1.x;
    keepvB[i] = rects[arrKasanari[i]].p1.y;
    keepvC[i] = rects[arrKasanari[i]].p2.x;
    keepvD[i] = rects[arrKasanari[i]].p2.y;
  }

  int ok = 1;
  for (int i = 0; i < (vn); ++i) {
    if (abcd == 0) rects[arrKasanari[i]].p2.x = rects[ite].p1.x;
    if (abcd == 1) rects[arrKasanari[i]].p2.y = rects[ite].p1.y;
    if (abcd == 2) rects[arrKasanari[i]].p1.x = rects[ite].p2.x;
    if (abcd == 3) rects[arrKasanari[i]].p1.y = rects[ite].p2.y;
    if (selfNg(arrKasanari[i])) ok = 0;
  }

  if (ok == 0) {
    for (int i = 0; i < (vn); ++i) {
      rects[arrKasanari[i]].p1.x = keepvA[i];
      rects[arrKasanari[i]].p1.y = keepvB[i];
      rects[arrKasanari[i]].p2.x = keepvC[i];
      rects[arrKasanari[i]].p2.y = keepvD[i];
    }
    // 元に戻す
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    return;
  }

  for (int i = 0; i < (vn); ++i) calc(arrKasanari[i]);
  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  double prob = exp((double)tmp / temp);

  if (prob > Rand01()) {
    modeCount[5]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      store_best();
    }
  }
  else {
    for (int i = 0; i < (vn); ++i) {
      rects[arrKasanari[i]].p1.x = keepvA[i];
      rects[arrKasanari[i]].p1.y = keepvB[i];
      rects[arrKasanari[i]].p2.x = keepvC[i];
      rects[arrKasanari[i]].p2.y = keepvD[i];
      calc(arrKasanari[i]);
    }
    // 元に戻す
    if (abcd == 0) rects[ite].p1.x -= diff;
    if (abcd == 1) rects[ite].p1.y -= diff;
    if (abcd == 2) rects[ite].p2.x -= diff;
    if (abcd == 3) rects[ite].p2.y -= diff;
    calc(ite);
  }
}

inline void shokiInit() {
  for (int i = 0; i < (n); ++i) {
    init_rect(rects[i], target_points[i]);
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    store_best();
  }
}

Rect best_best_best_rects[MAX_N];
int real_real_real_maxScore = -1;

int a2[100][MAX_N], b2[100][MAX_N], c2[100][MAX_N], d2[100][MAX_N];
int a4[100][MAX_N], b4[100][MAX_N], c4[100][MAX_N], d4[100][MAX_N];
int maxScore4[100] = {};

inline void Ui_Tei() {
  clock_t start, end;
  for (int ui_tei = 0; ui_tei < (5); ++ui_tei) {

    // 初期解
    // 左上(x,y)、右下(x+1,y+1)
    for (int i = 0; i < (n); ++i) {
      init_rect(rects[i], target_points[i]);
    }

    int T = 5;
    for (int _ = 0; _ < (T); ++_) {
      start = clock();

      // 初期スコア計算
      maxScore = calc(-1);
      store_best();

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
          int ite = Rand() % n;
          oneChange(ite, temp);
        }
        else if (mode == 1) {
          int ite = Rand() % n;
          Slide(ite, temp);
        }
        else if (now_time > 2.0 / T && mode == 2) {
          int ite = Rand() % n;
          aspectChange(ite, temp);
        }
        else if (mode == -3) {
          int ite = Rand() % n;
          fourChange(ite, temp);
        }
      }

      // 焼きなまし戻す
      maxScore = real_maxScore;
      for (int i = 0; i < (n); ++i) {
        rects[i] = best_rects[i];
      }
      calc(-1);

      if (maxScore > real_real_maxScore) {
        real_real_maxScore = maxScore;
        for (int i = 0; i < (n); ++i) {
          best_best_rects[i] = rects[i];
        }
      }
    }

    // real_real_maxScore戻す
    maxScore = real_real_maxScore;
    real_real_maxScore = 0;
    for (int i = 0; i < (n); ++i) {
      rects[i] = best_best_rects[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    for (int i = 0; i < (n); ++i) {
      best_rects[i] = rects[i];
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
        int ite = Rand() % n;
        oneChange(ite, temp);
      }
      else if (kouhan && mode == 1) { // 位置をスライド
        int ite = Rand() % n;
        Slide(ite, temp);
      }
      else if (mode == 2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
        int ite = Rand() % n;
        aspectChange(ite, temp);
      }
      else if (mode == -3) { // a,b,c,dを4つ同時に変更
        int ite = Rand() % n;
        fourChange(ite, temp);
      }
      else if (mode == 4) { // 膨らまして縮める
        int ite = Rand() % n;
        Extend(ite, temp);
      }
    }

    // 焼きなまし戻す
    maxScore = real_maxScore;
    for (int i = 0; i < (n); ++i) {
      rects[i] = best_rects[i];
    }
    calc(-1);

    if (maxScore > ui_tei_maxScore) {
      ui_tei_maxScore = maxScore;
      for (int i = 0; i < (n); ++i) {
        ui_tei_a[i] = rects[i].p1.x;
        ui_tei_b[i] = rects[i].p1.y;
        ui_tei_c[i] = rects[i].p2.x;
        ui_tei_d[i] = rects[i].p2.y;
      }
    }
  }

  // 元に戻しておく
  maxScore = 0;
  real_maxScore = 0;
  real_real_maxScore = 0;
  for (int i = 0; i < (n); ++i) {
    init_rect(rects[i], target_points[i]);
    best_rects[i] = rects[i];
    best_best_rects[i] = rects[i];
  }
}

int solve(int teisyutu, int fileNum) {
  auto startClock = system_clock::now();
  clock_t start, end;
  clock_t real_start = clock();

  nyuuryokuInit(fileNum);

  sortInit();


  for (int allLoop = 0; allLoop < (allLoopTimes); ++allLoop) {

    Ui_Tei();

    // ui_tei_maxScore戻す
    maxScore = ui_tei_maxScore;
    for (int i = 0; i < (n); ++i) {
      rects[i].p1.x = ui_tei_a[i];
      rects[i].p1.y = ui_tei_b[i];
      rects[i].p2.x = ui_tei_c[i];
      rects[i].p2.y = ui_tei_d[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    for (int i = 0; i < (n); ++i) {
      best_rects[i] = rects[i];
    }


    int oya = 1;

    for (int asai = 0; asai < (oya); ++asai) {
      for (int j = 0; j < (n); ++j) {
        a2[asai][j] = rects[j].p1.x;
        b2[asai][j] = rects[j].p1.y;
        c2[asai][j] = rects[j].p2.x;
        d2[asai][j] = rects[j].p2.y;
      }
    }

    // 焼きなまし2
    // 筋のいいやつを追う
    int T = 250 / allLoopTimes;
    for (int _ = 0; _ < (T); ++_) {
      for (int i = 0; i < (6); ++i) modeCount[i] = 0;

      int TT = 1;
      for (int asai = 0; asai < (TT); ++asai) {
        int kiyoshi = asai % oya;
        for (int i = 0; i < (n); ++i) {
          rects[i].p1.x = a2[kiyoshi][i];
          rects[i].p1.y = b2[kiyoshi][i];
          rects[i].p2.x = c2[kiyoshi][i];
          rects[i].p2.y = d2[kiyoshi][i];
        }

        // 初期スコア計算
        maxScore = calc(-1);
        real_maxScore = maxScore;

        for (int i = 0; i < (n); ++i) {
          best_rects[i] = rects[i];
        }

        // 焼きなまし(2回目)
        start = clock();
        end = clock();
        startClock = system_clock::now();
        double now_time = ((double)end - start) / CLOCKS_PER_SEC;
        double TL = (((realTL - 0.7) / T) / TT) / allLoopTimes;
        double start_temp = 20048.0;
        double end_temp = 0.1;
        double temp = start_temp + (end_temp - start_temp) * now_time / TL;
        int loop = 0;
        int kouhan = _ % 2;
        int tubusuFrequency = 30003;

        while (true) {
          loop++;
          if (loop % 100 == 1) {
            const double time = duration_cast<microseconds>(system_clock::now() - startClock).count() * 1e-6;
            if (time > TL) break;
            const double progressRatio = time / TL;   // 進捗。開始時が0.0、終了時が1.0
            temp = start_temp + (end_temp - start_temp) * progressRatio;
          }


          int mode = loop % 6;

          if (mode == -1) { // a,b,c,dのうち1つ変更
            int ite = Rand() % n;
            oneChange(ite, temp);
          }
          else if (mode == 1) { // 位置をスライド
            int ite = Rand() % n;
            Slide(ite, temp);
          }
          else if (mode == -2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
            int ite = Rand() % n;
            aspectChange(ite, temp);
          }
          else if (mode == -3) { // a,b,c,dを4つ同時に変更
            int ite = Rand() % n;
            fourChange(ite, temp);
          }
          else if (mode == 4) { // 膨らまして縮める
            int ite = Rand() % n;
            Extend(ite, temp);
          }
          else if (mode == 5) { // 境界をずらす
            int ite = Rand() % n;
            zurasi2(ite, temp);
          }

          if (loop % 34567 == 1120) {
            int ite = Rand() % n;
            ExtendKing(ite);
          }

          // 計算誤差解消?
          if (loop % 10000 == 1) {
            maxScore = calc(-1);
          }
        }

        // 焼きなまし戻す
        maxScore = real_maxScore;
        for (int i = 0; i < (n); ++i) {
          rects[i] = best_rects[i];
        }
        calc(-1);
        if (maxScore > real_real_maxScore) {
          real_real_maxScore = maxScore;
          for (int i = 0; i < (n); ++i) {
            best_best_rects[i] = rects[i];
          }
        }

        // ビームサーチの次の種にする
        maxScore4[asai] = maxScore;
        for (int i = 0; i < (n); ++i) {
          a4[asai][i] = rects[i].p1.x;
          b4[asai][i] = rects[i].p1.y;
          c4[asai][i] = rects[i].p2.x;
          d4[asai][i] = rects[i].p2.y;
        }
      }

      // 次の世代に継承
      vector<P> vBeam;
      for (int asai = 0; asai < (TT); ++asai) vBeam.emplace_back(P(maxScore4[asai], asai));
      sort(vBeam.begin(), vBeam.end(), greater<P>());


      for (int ii = 0; ii < (oya); ++ii) {
        int i = vBeam[ii].second;
        for (int j = 0; j < (n); ++j) {
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
      if (((double)end - real_start) / CLOCKS_PER_SEC > realTL) break;
    }

    // real_real_maxScore戻す
    maxScore = real_real_maxScore;
    for (int i = 0; i < (n); ++i) {
      rects[i] = best_best_rects[i];
    }
    calc(-1);

    if (teisyutu == 0) {
      cout << "maxScore = " << maxScore << endl;
    }

    const int MOD = 1000000007;
    if (teisyutu == 0 && maxScore > MOD) {
      cout << "ERROR" << endl;
      FileKakikomiERROR(fileNum);
    }

    // real_real_real入れる
    if (maxScore > real_real_real_maxScore && maxScore < 1000000007) {
      real_real_real_maxScore = maxScore;
      for (int i = 0; i < (n); ++i) {
        best_best_best_rects[i] = rects[i];
      }
    }


    // すべて白紙にリセットする
    maxScore = 0;
    real_maxScore = 0;
    real_real_maxScore = 0;
    for (int i = 0; i < (n); ++i) {
      init_rect(rects[i], target_points[i]);
      best_rects[i] = rects[i];
      best_best_rects[i] = rects[i];
    }
  }


  // real_real_real_maxScore戻す
  maxScore = real_real_real_maxScore;
  for (int i = 0; i < (n); ++i) {
    rects[i] = best_best_best_rects[i];
  }
  calc(-1);

  // 最終出力
  if (teisyutu) {
    for (int i = 0; i < (n); ++i) {
      cout << rects[i].p1.x << ' ' << rects[i].p1.y << ' ' << rects[i].p2.x << ' ' << rects[i].p2.y << endl;
    }
  }

  FileKakikomi(fileNum);

  // 提出時以下は消す
  if (teisyutu == 0) {
    cout << "file No. = " << fileNum << ", maxScore = " << maxScore << endl;
  }

  if (teisyutu == 0 && maxScore > 1000000007) {
    FileKakikomiERROR(fileNum);
  }

  return maxScore;
}

inline void clear_rect(Rect& r) {
  r.p1.x = 0;
  r.p1.y = 0;
  r.p2.x = 0;
  r.p2.y = 0;
}

inline void AllClear() {
  n = 0;
  maxScore = -1;
  real_maxScore = -1;
  pSum = 0;
  real_real_maxScore = -1;
  ui_tei_maxScore = -1;
  real_real_real_maxScore = -1;
  for (int i = 0; i < (MAX_N); ++i) {
    target_points[i].x = 0, target_points[i].y = 0, target_sizes[i] = 0;
    rects[i].p1.x = 0, rects[i].p1.y = 0, rects[i].p2.x = 0, rects[i].p2.y = 0;
    clear_rect(rects[i]);
    area_sizes[i] = 0;
    clear_rect(best_rects[i]);
    p[i] = 0;
    sort_x[i] = 0, sort_y[i] = 0;
    arg_sort_x[i] = 0, arg_sort_y[i] = 0;
    clear_rect(best_best_rects[i]);
    ui_tei_a[i] = 0, ui_tei_b[i] = 0, ui_tei_c[i] = 0, ui_tei_d[i] = 0;
    clear_rect(best_best_best_rects[i]);
  }
}

int main() {
  int teisyutu = 0;

  if (teisyutu) {
    solve(teisyutu, 0);
  }
  else {
    int mode = 0;
    if (mode == 0) { // コードテスト用
      solve(teisyutu, 0);
    }

    if (mode == 1) { // スコア確認用
      for (int i = 0; i < (1000); ++i) {
        srep(i, 0, 50) {
          for (int j = 0; j < (10); ++j) {
            AllClear();
            solve(teisyutu, i);
          }
        }
      }
    }

    if (mode == 2) { // ハイパラいじり用
      string fileName = "haipara_hoge.txt";
      const char* cstr = fileName.c_str();
      ofstream ofs(cstr);


      vector<int> arrHaiparaOya = { 1,2,4,8 };
      int m = arrHaiparaOya.size();

      for (int k = 0; k < (m); ++k) {
        for (int l = 0; l < (4); ++l) {
          haipara_oya = arrHaiparaOya[k];
          if (l == 0) haipara_TT = haipara_oya;
          if (l == 1) haipara_TT = haipara_oya + 1;
          if (l == 2) haipara_TT = haipara_oya * 2;
          if (l == 3) haipara_TT = haipara_oya * 4;

          ll sum = 0;
          srep(i, 0, 50) {
            for (int j = 0; j < (1); ++j) {
              AllClear();
              sum += solve(teisyutu, i);
            }
          }
          cout << "haipara_oya = " << haipara_oya;
          cout << ", haipara_TT = " << haipara_TT;
          cout << ", sum = " << sum << endl;
          ofs << "haipara_oya = " << haipara_oya;
          ofs << ", haipara_TT = " << haipara_TT;
          ofs << ", sum = " << sum << endl;
        }
      }

      ofs.close();
    }

  }

  return 0;
}


