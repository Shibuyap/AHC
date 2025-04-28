#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
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
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for (int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
using namespace chrono;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 205

static uint32_t Rand()
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
static double Rand01()
{
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

struct Point
{
  int x;
  int y;
};

struct Rect
{
  Point p1;
  Point p2;
};

int allLoopTimes = 1;
int n;
Point point[MAX_N];
int target_sizes[MAX_N];
Rect rect[MAX_N];
int area_sizes[MAX_N];
int real_a[MAX_N], real_b[MAX_N], real_c[MAX_N], real_d[MAX_N];

inline void nyuuryokuInit(int fileNum)
{
  // 入力
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << fileNum << ".txt";
  ifstream ifs(oss.str());
  if (!ifs.is_open()) { // 標準入力する
    cin >> n;
    rep(i, n) cin >> point[i].x >> point[i].y >> target_sizes[i];
  }
  else { // ファイル入力する
    ifs >> n;
    rep(i, n) ifs >> point[i].x >> point[i].y >> target_sizes[i];
  }
}

inline void FileKakikomi(int fileNum)
{
  string fileName = to_string(fileNum);
  fileName += "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);
  rep(i, n)
  {
    ofs << rect[i].p1.x << ' ' << rect[i].p1.y << ' ' << rect[i].p2.x << ' ' << rect[i].p2.y << endl;
  }
  ofs.close();
}

inline void FileKakikomiERROR(int fileNum)
{
  string fileName = to_string(fileNum);
  fileName += "_out_ERROR.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);
  rep(i, n)
  {
    ofs << rect[i].p1.x << ' ' << rect[i].p1.y << ' ' << rect[i].p2.x << ' ' << rect[i].p2.y << endl;
  }
  ofs.close();
}

int maxScore = -1;
int real_maxScore = -1;

double p[MAX_N];
double pSum;
inline int calc(int ite)
{
  if (ite == -1) {
    double sum = 0;
    rep(i, n)
    {
      area_sizes[i] = (rect[i].p2.x - rect[i].p1.x) * (rect[i].p2.y - rect[i].p1.y);
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
    area_sizes[ite] = (rect[ite].p2.x - rect[ite].p1.x) * (rect[ite].p2.y - rect[ite].p1.y);
    p[ite] = 1.0 - (1.0 - (double)min(target_sizes[ite], area_sizes[ite]) / (double)max(target_sizes[ite], area_sizes[ite])) * (1.0 - (double)min(target_sizes[ite], area_sizes[ite]) / (double)max(target_sizes[ite], area_sizes[ite]));
    sum += p[ite];
    pSum = sum;
    sum /= (double)n;
    sum *= 1000000000.0;
    return round(sum);
  }
}

inline int kasanarihantei(int i, int j)
{
  int cnt = 0;
  if (rect[i].p1.x <= rect[j].p1.x && rect[j].p1.x < rect[i].p2.x) cnt++;
  else if (rect[j].p1.x <= rect[i].p1.x && rect[i].p1.x < rect[j].p2.x) cnt++;
  if (rect[i].p1.y <= rect[j].p1.y && rect[j].p1.y < rect[i].p2.y) cnt++;
  else if (rect[j].p1.y <= rect[i].p1.y && rect[i].p1.y < rect[j].p2.y) cnt++;
  return cnt == 2;
}

int sort_x[MAX_N], sort_y[MAX_N];
int arg_sort_x[MAX_N], arg_sort_y[MAX_N];
inline void sortInit()
{
  vector<P> v;
  rep(i, n)
  {
    v.emplace_back(P(point[i].x, i));
  }
  sort(v.begin(), v.end());
  rep(i, n)
  {
    sort_x[i] = v[i].second;
    arg_sort_x[v[i].second] = i;
  }

  v.clear();
  rep(i, n)
  {
    v.emplace_back(P(point[i].y, i));
  }
  sort(v.begin(), v.end());
  rep(i, n)
  {
    sort_y[i] = v[i].second;
    arg_sort_y[v[i].second] = i;
  }
}

// 0~10000を出ていないか
// 面積は1以上か
// (x[i]+0.5,y[i]+0.5)を含んでいるか
// 重なりがないか
inline int isOK2(int ite)
{
  if (ite == -1) {
    rep(i, n)
    {
      if (rect[i].p1.x < 0 || 10000 < rect[i].p1.x) return 0;
      if (rect[i].p1.y < 0 || 10000 < rect[i].p1.y) return 0;
      if (rect[i].p2.x < 0 || 10000 < rect[i].p2.x) return 0;
      if (rect[i].p2.y < 0 || 10000 < rect[i].p2.y) return 0;
    }
    rep(i, n)
    {
      if (rect[i].p2.x <= rect[i].p1.x) return 0;
      if (rect[i].p2.y <= rect[i].p1.y) return 0;
    }
    rep(i, n)
    {
      if (point[i].x < rect[i].p1.x || rect[i].p2.x <= point[i].x) return 0;
      if (point[i].y < rect[i].p1.y || rect[i].p2.y <= point[i].y) return 0;
    }
    rep(i, n)
    {
      srep(j, i + 1, n)
      {
        if (kasanarihantei(i, j)) return 0;
      }
    }
  }
  else {
    if (rect[ite].p1.x < 0 || 10000 < rect[ite].p1.x) return 0;
    if (rect[ite].p1.y < 0 || 10000 < rect[ite].p1.y) return 0;
    if (rect[ite].p2.x < 0 || 10000 < rect[ite].p2.x) return 0;
    if (rect[ite].p2.y < 0 || 10000 < rect[ite].p2.y) return 0;
    if (rect[ite].p2.x <= rect[ite].p1.x) return 0;
    if (rect[ite].p2.y <= rect[ite].p1.y) return 0;
    if (point[ite].x < rect[ite].p1.x || rect[ite].p2.x <= point[ite].x) return 0;
    if (point[ite].y < rect[ite].p1.y || rect[ite].p2.y <= point[ite].y) return 0;
    rep(i, n)
    {
      if (i == ite) continue;
      if (kasanarihantei(i, ite)) return 0;
    }
  }
  return 1;
}

inline int isOK(int ite)
{
  if (ite == -1) {
    rep(i, n)
    {
      if (rect[i].p1.x < 0 || 10000 < rect[i].p1.x) return 0;
      if (rect[i].p1.y < 0 || 10000 < rect[i].p1.y) return 0;
      if (rect[i].p2.x < 0 || 10000 < rect[i].p2.x) return 0;
      if (rect[i].p2.y < 0 || 10000 < rect[i].p2.y) return 0;
    }
    rep(i, n)
    {
      if (rect[i].p2.x <= rect[i].p1.x) return 0;
      if (rect[i].p2.y <= rect[i].p1.y) return 0;
    }
    rep(i, n)
    {
      if (point[i].x < rect[i].p1.x || rect[i].p2.x <= point[i].x) return 0;
      if (point[i].y < rect[i].p1.y || rect[i].p2.y <= point[i].y) return 0;
    }
    rep(i, n)
    {
      srep(j, i + 1, n)
      {
        if (kasanarihantei(i, j)) return 0;
      }
    }
  }
  else {
    if (rect[ite].p1.x < 0 || 10000 < rect[ite].p1.x) return 0;
    if (rect[ite].p1.y < 0 || 10000 < rect[ite].p1.y) return 0;
    if (rect[ite].p2.x < 0 || 10000 < rect[ite].p2.x) return 0;
    if (rect[ite].p2.y < 0 || 10000 < rect[ite].p2.y) return 0;
    if (rect[ite].p2.x <= rect[ite].p1.x) return 0;
    if (rect[ite].p2.y <= rect[ite].p1.y) return 0;
    if (point[ite].x < rect[ite].p1.x || rect[ite].p2.x <= point[ite].x) return 0;
    if (point[ite].y < rect[ite].p1.y || rect[ite].p2.y <= point[ite].y) return 0;
    int argX = arg_sort_x[ite];
    int nowLeft = rect[ite].p1.y;
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (rect[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.y);
        if (nowLeft >= rect[ite].p2.y) break;
      }
    }
    nowLeft = rect[ite].p1.y;
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (rect[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.y);
        if (nowLeft >= rect[ite].p2.y) break;
      }
    }
  }
  return 1;
}

int vExtend[4];
inline void hukuramashi(int ite)
{
  vExtend[0] = point[ite].x;
  vExtend[1] = point[ite].y;
  vExtend[2] = point[ite].x + 1;
  vExtend[3] = point[ite].y + 1;

  int flagTateYoko = Rand() % 2;
  if (flagTateYoko == 0) {
    vExtend[0] = 0;
    vExtend[2] = 10000;
    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rect[i].p1.y <= vExtend[1] && vExtend[1] < rect[i].p2.y) flagKasanari = 1;
      if (vExtend[1] <= rect[i].p1.y && rect[i].p1.y < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtend[0] = max(vExtend[0], rect[i].p2.x);
        }
        else {
          vExtend[2] = min(vExtend[2], rect[i].p1.x);
        }
        break;
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rect[i].p1.y <= vExtend[1] && vExtend[1] < rect[i].p2.y) flagKasanari = 1;
      if (vExtend[1] <= rect[i].p1.y && rect[i].p1.y < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtend[0] = max(vExtend[0], rect[i].p2.x);
        }
        else {
          vExtend[2] = min(vExtend[2], rect[i].p1.x);
        }
        break;
      }
    }

    vExtend[1] = 0;
    vExtend[3] = 10000;
    int argY = arg_sort_y[ite];
    int nowLeft = vExtend[0];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rect[i].p1.x <= vExtend[0] && vExtend[0] < rect[i].p2.x) flagKasanari = 1;
      if (vExtend[0] <= rect[i].p1.x && rect[i].p1.x < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtend[1] = max(vExtend[1], rect[i].p2.y);
        }
        else {
          vExtend[3] = min(vExtend[3], rect[i].p1.y);
        }
        if (rect[i].p1.x <= nowLeft) {
          nowLeft = max(nowLeft, rect[i].p2.x);
          if (vExtend[2] <= nowLeft) break;
        }
      }
    }

    nowLeft = vExtend[0];
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rect[i].p1.x <= vExtend[0] && vExtend[0] < rect[i].p2.x) flagKasanari = 1;
      if (vExtend[0] <= rect[i].p1.x && rect[i].p1.x < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtend[1] = max(vExtend[1], rect[i].p2.y);
        }
        else {
          vExtend[3] = min(vExtend[3], rect[i].p1.y);
        }
        if (rect[i].p1.x <= nowLeft) {
          nowLeft = max(nowLeft, rect[i].p2.x);
          if (vExtend[2] <= nowLeft) break;
        }
      }
    }
  }
  else {
    vExtend[1] = 0;
    vExtend[3] = 10000;
    int argY = arg_sort_y[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rect[i].p1.x <= vExtend[0] && vExtend[0] < rect[i].p2.x) flagKasanari = 1;
      if (vExtend[0] <= rect[i].p1.x && rect[i].p1.x < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtend[1] = max(vExtend[1], rect[i].p2.y);
        }
        else {
          vExtend[3] = min(vExtend[3], rect[i].p1.y);
        }
        break;
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (rect[i].p1.x <= vExtend[0] && vExtend[0] < rect[i].p2.x) flagKasanari = 1;
      if (vExtend[0] <= rect[i].p1.x && rect[i].p1.x < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtend[1] = max(vExtend[1], rect[i].p2.y);
        }
        else {
          vExtend[3] = min(vExtend[3], rect[i].p1.y);
        }
        break;
      }
    }

    vExtend[0] = 0;
    vExtend[2] = 10000;
    int argX = arg_sort_x[ite];
    int nowLeft = vExtend[1];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rect[i].p1.y <= vExtend[1] && vExtend[1] < rect[i].p2.y) flagKasanari = 1;
      if (vExtend[1] <= rect[i].p1.y && rect[i].p1.y < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtend[0] = max(vExtend[0], rect[i].p2.x);
        }
        else {
          vExtend[2] = min(vExtend[2], rect[i].p1.x);
        }
        if (rect[i].p1.y <= nowLeft) {
          nowLeft = max(nowLeft, rect[i].p2.y);
          if (vExtend[3] <= nowLeft) break;
        }
      }
    }
    nowLeft = vExtend[1];
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (rect[i].p1.y <= vExtend[1] && vExtend[1] < rect[i].p2.y) flagKasanari = 1;
      if (vExtend[1] <= rect[i].p1.y && rect[i].p1.y < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtend[0] = max(vExtend[0], rect[i].p2.x);
        }
        else {
          vExtend[2] = min(vExtend[2], rect[i].p1.x);
        }
        if (rect[i].p1.y <= nowLeft) {
          nowLeft = max(nowLeft, rect[i].p2.y);
          if (vExtend[3] <= nowLeft) break;
        }
      }
    }
  }

  int shuf[4] = {};
  int shuffleSeed = Rand() % 24;
  rep(j, 4) shuf[j] = shuffles[shuffleSeed][j];

  int yure = Rand() % 2;

  rep(i, 4)
  {
    int S = (vExtend[2] - vExtend[0]) * (vExtend[3] - vExtend[1]);
    if (S <= target_sizes[ite]) break;
    if (shuf[i] == 0) {
      int ma = target_sizes[ite] / (vExtend[3] - vExtend[1]) + yure;
      int diff = (vExtend[2] - vExtend[0]) - ma;
      int capa = point[ite].x - vExtend[0];
      if (capa >= diff) {
        vExtend[0] += diff;
      }
      else {
        vExtend[0] += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = target_sizes[ite] / (vExtend[2] - vExtend[0]) + yure;
      int diff = (vExtend[3] - vExtend[1]) - ma;
      int capa = point[ite].y - vExtend[1];
      if (capa >= diff) {
        vExtend[1] += diff;
      }
      else {
        vExtend[1] += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = target_sizes[ite] / (vExtend[3] - vExtend[1]) + yure;
      int diff = (vExtend[2] - vExtend[0]) - ma;
      int capa = vExtend[2] - (point[ite].x + 1);
      if (capa >= diff) {
        vExtend[2] -= diff;
      }
      else {
        vExtend[2] -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = target_sizes[ite] / (vExtend[2] - vExtend[0]) + yure;
      int diff = (vExtend[3] - vExtend[1]) - ma;
      int capa = vExtend[3] - (point[ite].y + 1);
      if (capa >= diff) {
        vExtend[3] -= diff;
      }
      else {
        vExtend[3] -= capa;
      }
    }
  }

  int ng = 0;
  if (vExtend[0] < 0 || 10000 < vExtend[0]) ng = 1;
  if (vExtend[1] < 0 || 10000 < vExtend[1]) ng = 1;
  if (vExtend[2] < 0 || 10000 < vExtend[2]) ng = 1;
  if (vExtend[3] < 0 || 10000 < vExtend[3]) ng = 1;
  if (vExtend[2] <= vExtend[0]) ng = 1;
  if (vExtend[3] <= vExtend[1]) ng = 1;
  if (point[ite].x < vExtend[0] || vExtend[2] <= point[ite].x) ng = 1;
  if (point[ite].y < vExtend[1] || vExtend[3] <= point[ite].y) ng = 1;

  if (ng) {
    vExtend[0] = point[ite].x; vExtend[2] = point[ite].x + 1;
    vExtend[1] = point[ite].y; vExtend[3] = point[ite].y + 1;
  }

}

inline void Extend(int ite, double temp)
{
  int keepA = rect[ite].p1.x;
  int keepB = rect[ite].p1.y;
  int keepC = rect[ite].p2.x;
  int keepD = rect[ite].p2.y;

  hukuramashi(ite);
  rect[ite].p1.x = vExtend[0];
  rect[ite].p1.y = vExtend[1];
  rect[ite].p2.x = vExtend[2];
  rect[ite].p2.y = vExtend[3];

  int tmpScore = calc(ite);

  // int tmp = tmpScore - maxScore;
  // double prob = exp((double)tmp / temp);
  // if(prob > Rand01())

  if (tmpScore >= maxScore) {
    modeCount[4]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    // 元に戻す
    rect[ite].p1.x = keepA;
    rect[ite].p1.y = keepB;
    rect[ite].p2.x = keepC;
    rect[ite].p2.y = keepD;
    calc(ite);
  }
}

int real_real_a[MAX_N], real_real_b[MAX_N], real_real_c[MAX_N], real_real_d[MAX_N];
int real_real_maxScore = -1;

int ui_tei_a[MAX_N], ui_tei_b[MAX_N], ui_tei_c[MAX_N], ui_tei_d[MAX_N];
int ui_tei_maxScore = -1;

inline void Tubusu(int tubusu)
{
  rep(i, tubusu)
  { // tubusu個つぶす
    int ite = Rand() % n;
    rect[ite].p1.x = point[ite].x;
    rect[ite].p1.y = point[ite].y;
    rect[ite].p2.x = point[ite].x + 1;
    rect[ite].p2.y = point[ite].y + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }
  }
}

inline void TubusuWorst(int tubusu_worst)
{
  vector<pair<double, int>> v;
  rep(i, n)
  {
    v.emplace_back(pair<double, int>(p[i], i));
  }
  sort(v.begin(), v.end());
  rep(i, tubusu_worst)
  {
    int ite = v[i].second;
    rect[ite].p1.x = point[ite].x;
    rect[ite].p1.y = point[ite].y;
    rect[ite].p2.x = point[ite].x + 1;
    rect[ite].p2.y = point[ite].y + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }
  }
}

inline void AnaWoAkeru(int hole)
{
  int ite = Rand() % n;
  vector<int> keep;
  keep.emplace_back(ite);
  rect[ite].p1.x -= 100;
  rect[ite].p1.y -= 100;
  rect[ite].p2.x += 100;
  rect[ite].p2.y += 100;
  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) keep.emplace_back(i);
  }
  int keepSize = keep.size();
  rep(i, keepSize)
  {
    rect[keep[i]].p1.x = point[keep[i]].x;
    rect[keep[i]].p1.y = point[keep[i]].y;
    rect[keep[i]].p2.x = point[keep[i]].x + 1;
    rect[keep[i]].p2.y = point[keep[i]].y + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }
  }
}



int vExtendKing[4];
inline void hukuramashiKing(int ite)
{
  vExtendKing[0] = max(0, (int)(point[ite].x - Rand() % 1000));
  vExtendKing[1] = max(0, (int)(point[ite].y - Rand() % 1000));
  vExtendKing[2] = min(10000, (int)(point[ite].x + 1 + Rand() % 1000));
  vExtendKing[3] = min(10000, (int)(point[ite].y + 1 + Rand() % 1000));

  int tateyoko = Rand() % 2;

  if (tateyoko == 0) {
    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (point[i].x == point[ite].x) continue;
      int flagKasanari = 0;
      if (point[i].y <= vExtendKing[1] && vExtendKing[1] < point[i].y + 1) flagKasanari = 1;
      if (vExtendKing[1] <= point[i].y && point[i].y < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtendKing[0] = max(vExtendKing[0], point[i].x + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], point[i].x);
        }
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (point[i].x == point[ite].x) continue;
      int flagKasanari = 0;
      if (point[i].y <= vExtendKing[1] && vExtendKing[1] < point[i].y + 1) flagKasanari = 1;
      if (vExtendKing[1] <= point[i].y && point[i].y < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtendKing[0] = max(vExtendKing[0], point[i].x + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], point[i].x);
        }
      }
    }

    int argY = arg_sort_y[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (point[i].y == point[ite].y) continue;
      int flagKasanari = 0;
      if (point[i].x <= vExtendKing[0] && vExtendKing[0] < point[i].x + 1) flagKasanari = 1;
      if (vExtendKing[0] <= point[i].x && point[i].x < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtendKing[1] = max(vExtendKing[1], point[i].y + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], point[i].y);
        }
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (point[i].y == point[ite].y) continue;
      int flagKasanari = 0;
      if (point[i].x <= vExtendKing[0] && vExtendKing[0] < point[i].x + 1) flagKasanari = 1;
      if (vExtendKing[0] <= point[i].x && point[i].x < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtendKing[1] = max(vExtendKing[1], point[i].y + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], point[i].y);
        }
      }
    }
  }
  else {
    int argY = arg_sort_y[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (point[i].y == point[ite].y) continue;
      int flagKasanari = 0;
      if (point[i].x <= vExtendKing[0] && vExtendKing[0] < point[i].x + 1) flagKasanari = 1;
      if (vExtendKing[0] <= point[i].x && point[i].x < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtendKing[1] = max(vExtendKing[1], point[i].y + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], point[i].y);
        }
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (point[i].y == point[ite].y) continue;
      int flagKasanari = 0;
      if (point[i].x <= vExtendKing[0] && vExtendKing[0] < point[i].x + 1) flagKasanari = 1;
      if (vExtendKing[0] <= point[i].x && point[i].x < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].y <= point[ite].y) {
          vExtendKing[1] = max(vExtendKing[1], point[i].y + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], point[i].y);
        }
      }
    }

    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (point[i].x == point[ite].x) continue;
      int flagKasanari = 0;
      if (point[i].y <= vExtendKing[1] && vExtendKing[1] < point[i].y + 1) flagKasanari = 1;
      if (vExtendKing[1] <= point[i].y && point[i].y < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtendKing[0] = max(vExtendKing[0], point[i].x + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], point[i].x);
        }
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (point[i].x == point[ite].x) continue;
      int flagKasanari = 0;
      if (point[i].y <= vExtendKing[1] && vExtendKing[1] < point[i].y + 1) flagKasanari = 1;
      if (vExtendKing[1] <= point[i].y && point[i].y < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (point[i].x <= point[ite].x) {
          vExtendKing[0] = max(vExtendKing[0], point[i].x + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], point[i].x);
        }
      }
    }
  }



  int shuf[4] = {};
  int shuffleSeed = Rand() % 24;
  rep(j, 4) shuf[j] = shuffles[shuffleSeed][j];

  int yure = Rand() % 2;

  rep(i, 4)
  {
    int S = (vExtendKing[2] - vExtendKing[0]) * (vExtendKing[3] - vExtendKing[1]);
    if (S <= target_sizes[ite]) break;
    if (shuf[i] == 0) {
      int ma = target_sizes[ite] / (vExtendKing[3] - vExtendKing[1]) + yure;
      int diff = (vExtendKing[2] - vExtendKing[0]) - ma;
      if (diff < 0) diff = 0;
      int capa = point[ite].x - vExtendKing[0];
      if (capa >= diff) {
        vExtendKing[0] += diff;
      }
      else {
        vExtendKing[0] += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = target_sizes[ite] / (vExtendKing[2] - vExtendKing[0]) + yure;
      int diff = (vExtendKing[3] - vExtendKing[1]) - ma;
      if (diff < 0) diff = 0;
      int capa = point[ite].y - vExtendKing[1];
      if (capa >= diff) {
        vExtendKing[1] += diff;
      }
      else {
        vExtendKing[1] += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = target_sizes[ite] / (vExtendKing[3] - vExtendKing[1]) + yure;
      int diff = (vExtendKing[2] - vExtendKing[0]) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing[2] - (point[ite].x + 1);
      if (capa >= diff) {
        vExtendKing[2] -= diff;
      }
      else {
        vExtendKing[2] -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = target_sizes[ite] / (vExtendKing[2] - vExtendKing[0]) + yure;
      int diff = (vExtendKing[3] - vExtendKing[1]) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing[3] - (point[ite].y + 1);
      if (capa >= diff) {
        vExtendKing[3] -= diff;
      }
      else {
        vExtendKing[3] -= capa;
      }
    }
  }

  int ng = 0;
  if (vExtendKing[0] < 0 || 10000 < vExtendKing[0]) ng = 1;
  if (vExtendKing[1] < 0 || 10000 < vExtendKing[1]) ng = 1;
  if (vExtendKing[2] < 0 || 10000 < vExtendKing[2]) ng = 1;
  if (vExtendKing[3] < 0 || 10000 < vExtendKing[3]) ng = 1;
  if (vExtendKing[2] <= vExtendKing[0]) ng = 1;
  if (vExtendKing[3] <= vExtendKing[1]) ng = 1;
  if (point[ite].x < vExtendKing[0] || vExtendKing[2] <= point[ite].x) ng = 1;
  if (point[ite].y < vExtendKing[1] || vExtendKing[3] <= point[ite].y) ng = 1;

  if (ng) {
    vExtendKing[0] = point[ite].x; vExtendKing[2] = point[ite].x + 1;
    vExtendKing[1] = point[ite].y; vExtendKing[3] = point[ite].y + 1;
  }

}

inline void ExtendKing(int ite)
{
  int keepA = rect[ite].p1.x;
  int keepB = rect[ite].p1.y;
  int keepC = rect[ite].p2.x;
  int keepD = rect[ite].p2.y;

  hukuramashiKing(ite);
  rect[ite].p1.x = vExtendKing[0];
  rect[ite].p1.y = vExtendKing[1];
  rect[ite].p2.x = vExtendKing[2];
  rect[ite].p2.y = vExtendKing[3];


  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      rect[i].p1.x = point[i].x;
      rect[i].p1.y = point[i].y;
      rect[i].p2.x = point[i].x + 1;
      rect[i].p2.y = point[i].y + 1;
    }
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }
  }
}

inline void oneChange(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int abcd = Rand() % 4;

  if (abcd == 0) rect[ite].p1.x += diff;
  if (abcd == 1) rect[ite].p1.y += diff;
  if (abcd == 2) rect[ite].p2.x += diff;
  if (abcd == 3) rect[ite].p2.y += diff;

  if (isOK(ite) == 0) {
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  const double prob = exp((double)tmp / temp);

  if (prob > Rand01()) {
    modeCount[0]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    // 元に戻す
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    calc(ite);
  }
}

inline void fourChange(int ite, double temp)
{
  int diffA = Rand() % 101 - 50;
  int diffB = Rand() % 101 - 50;
  int diffC = Rand() % 101 - 50;
  int diffD = Rand() % 101 - 50;

  rect[ite].p1.x += diffA;
  rect[ite].p1.y += diffB;
  rect[ite].p2.x += diffC;
  rect[ite].p2.y += diffD;

  if (isOK(ite) == 0) {
    rect[ite].p1.x -= diffA;
    rect[ite].p1.y -= diffB;
    rect[ite].p2.x -= diffC;
    rect[ite].p2.y -= diffD;
    return;
  }

  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;

  const double prob = exp((double)tmp / temp);
  if (prob > Rand01()) {
    modeCount[3]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    // 元に戻す
    rect[ite].p1.x -= diffA;
    rect[ite].p1.y -= diffB;
    rect[ite].p2.x -= diffC;
    rect[ite].p2.y -= diffD;
    calc(ite);
  }
}

inline void Slide(int ite)
{
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int ab = Rand() % 2;

  if (ab == 0) {
    rect[ite].p1.x += diff;
    rect[ite].p2.x += diff;
  }
  if (ab == 1) {
    rect[ite].p1.y += diff;
    rect[ite].p2.y += diff;
  }

  if (isOK(ite) == 0) {
    if (ab == 0) {
      rect[ite].p1.x -= diff;
      rect[ite].p2.x -= diff;
    }
    if (ab == 1) {
      rect[ite].p1.y -= diff;
      rect[ite].p2.y -= diff;
    }
    return;
  }

  int tmpScore = calc(ite);

  if (tmpScore >= maxScore) {
    modeCount[1]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    // 元に戻す
    if (ab == 0) {
      rect[ite].p1.x -= diff;
      rect[ite].p2.x -= diff;
    }
    if (ab == 1) {
      rect[ite].p1.y -= diff;
      rect[ite].p2.y -= diff;
    }
    calc(ite);
  }
}

inline void aspectChange(int ite)
{
  int yokoRatio = Rand() % 9 + 1; // 1 ~ 9;
  int tateRatio = 10 - yokoRatio;

  int S = yokoRatio * tateRatio;
  int mul = sqrt(target_sizes[ite] / S);
  if (mul == 0) return;

  int yoko = yokoRatio * mul;
  int tate = tateRatio * mul;

  int keepA = rect[ite].p1.x;
  int keepB = rect[ite].p1.y;
  int keepC = rect[ite].p2.x;
  int keepD = rect[ite].p2.y;

  int leftA = max(0, point[ite].x - (yoko - 1));
  int rightA = min(point[ite].x, 10000 - yoko);
  int rangeA = rightA - leftA + 1;
  if (rangeA < 1) return;

  int leftB = max(0, point[ite].y - (tate - 1));
  int rightB = min(point[ite].y, 10000 - tate);
  int rangeB = rightB - leftB + 1;
  if (rangeB < 1) return;

  rect[ite].p1.x = Rand() % rangeA + leftA;
  rect[ite].p2.x = rect[ite].p1.x + rangeA;
  rect[ite].p1.y = Rand() % rangeB + leftB;
  rect[ite].p2.y = rect[ite].p1.y + rangeB;

  if (isOK(ite) == 0) {
    rect[ite].p1.x = keepA;
    rect[ite].p1.y = keepB;
    rect[ite].p2.x = keepC;
    rect[ite].p2.y = keepD;
    return;
  }

  int tmpScore = calc(ite);

  if (tmpScore >= maxScore) {
    modeCount[2]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    // 元に戻す
    rect[ite].p1.x = keepA;
    rect[ite].p1.y = keepB;
    rect[ite].p2.x = keepC;
    rect[ite].p2.y = keepD;
    calc(ite);
  }
}

inline int selfNg(int ite)
{
  if (rect[ite].p1.x < 0 || 10000 < rect[ite].p1.x) return 1;
  if (rect[ite].p1.y < 0 || 10000 < rect[ite].p1.y) return 1;
  if (rect[ite].p2.x < 0 || 10000 < rect[ite].p2.x) return 1;
  if (rect[ite].p2.y < 0 || 10000 < rect[ite].p2.y) return 1;
  if (rect[ite].p2.x <= rect[ite].p1.x) return 1;
  if (rect[ite].p2.y <= rect[ite].p1.y) return 1;
  if (point[ite].x < rect[ite].p1.x || rect[ite].p2.x <= point[ite].x) return 1;
  if (point[ite].y < rect[ite].p1.y || rect[ite].p2.y <= point[ite].y) return 1;
  return 0;
}

inline int dokasuOK(int ite, int abcd)
{
  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      if (abcd == 0) rect[i].p2.x = rect[ite].p1.x;
      if (abcd == 1) rect[i].p2.y = rect[ite].p1.y;
      if (abcd == 2) rect[i].p1.x = rect[ite].p2.x;
      if (abcd == 3) rect[i].p1.y = rect[ite].p2.y;

      if (selfNg(i)) return 0;
    }
  }
  return 1;
}

int arrKasanari[MAX_N];
int kasanariCount;
inline void kasanaritati(int ite, int abcd)
{
  kasanariCount = 0;
  if (abcd == 0) {
    int argX = arg_sort_x[ite];
    int nowLeft = rect[ite].p1.y;
    int nowRight = rect[ite].p2.y;
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (rect[ite].p1.x <= point[i].x) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rect[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.y);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rect[i].p2.y) {
        nowRight = min(nowRight, rect[i].p1.y);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 1) {
    int argY = arg_sort_y[ite];
    int nowLeft = rect[ite].p1.x;
    int nowRight = rect[ite].p2.x;
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (rect[ite].p1.y <= point[i].y) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rect[i].p1.x <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.x);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rect[i].p2.x) {
        nowRight = min(nowRight, rect[i].p1.x);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 2) {
    int argX = arg_sort_x[ite];
    int nowLeft = rect[ite].p1.y;
    int nowRight = rect[ite].p2.y;
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (point[i].x < rect[ite].p2.x) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rect[i].p1.y <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.y);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rect[i].p2.y) {
        nowRight = min(nowRight, rect[i].p1.y);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 3) {
    int argY = arg_sort_y[ite];
    int nowLeft = rect[ite].p1.x;
    int nowRight = rect[ite].p2.x;
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (point[i].y < rect[ite].p2.y) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (rect[i].p1.x <= nowLeft) {
        nowLeft = max(nowLeft, rect[i].p2.x);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= rect[i].p2.x) {
        nowRight = min(nowRight, rect[i].p1.x);
        if (nowLeft >= nowRight) break;
      }
    }
  }
}

int keepvA[MAX_N], keepvB[MAX_N], keepvC[MAX_N], keepvD[MAX_N];
inline void zurasi2(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = Rand() % 50 + 1;
  int abcd = Rand() % 4;

  if (abcd < 2) diff *= -1;

  if (abcd == 0) rect[ite].p1.x += diff;
  if (abcd == 1) rect[ite].p1.y += diff;
  if (abcd == 2) rect[ite].p2.x += diff;
  if (abcd == 3) rect[ite].p2.y += diff;

  if (selfNg(ite)) {
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    return;
  }

  kasanaritati(ite, abcd);
  int vn = kasanariCount;

  if (vn > 0 && arrKasanari[0] == -1) {
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    return;
  }


  rep(i, vn)
  {
    keepvA[i] = rect[arrKasanari[i]].p1.x;
    keepvB[i] = rect[arrKasanari[i]].p1.y;
    keepvC[i] = rect[arrKasanari[i]].p2.x;
    keepvD[i] = rect[arrKasanari[i]].p2.y;
  }

  int ok = 1;
  rep(i, vn)
  {
    if (abcd == 0) rect[arrKasanari[i]].p2.x = rect[ite].p1.x;
    if (abcd == 1) rect[arrKasanari[i]].p2.y = rect[ite].p1.y;
    if (abcd == 2) rect[arrKasanari[i]].p1.x = rect[ite].p2.x;
    if (abcd == 3) rect[arrKasanari[i]].p1.y = rect[ite].p2.y;
    if (selfNg(arrKasanari[i])) ok = 0;
  }

  if (ok == 0) {
    rep(i, vn)
    {
      rect[arrKasanari[i]].p1.x = keepvA[i];
      rect[arrKasanari[i]].p1.y = keepvB[i];
      rect[arrKasanari[i]].p2.x = keepvC[i];
      rect[arrKasanari[i]].p2.y = keepvD[i];
    }
    // 元に戻す
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    return;
  }

  rep(i, vn) calc(arrKasanari[i]);
  int tmpScore = calc(ite);

  int tmp = tmpScore - maxScore;
  double prob = exp((double)tmp / temp);

  if (prob > Rand01()) {
    modeCount[5]++;
    maxScore = tmpScore;
    if (maxScore > real_maxScore) {
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }
    }
  }
  else {
    rep(i, vn)
    {
      rect[arrKasanari[i]].p1.x = keepvA[i];
      rect[arrKasanari[i]].p1.y = keepvB[i];
      rect[arrKasanari[i]].p2.x = keepvC[i];
      rect[arrKasanari[i]].p2.y = keepvD[i];
      calc(arrKasanari[i]);
    }
    // 元に戻す
    if (abcd == 0) rect[ite].p1.x -= diff;
    if (abcd == 1) rect[ite].p1.y -= diff;
    if (abcd == 2) rect[ite].p2.x -= diff;
    if (abcd == 3) rect[ite].p2.y -= diff;
    calc(ite);
  }
}

inline void shokiInit()
{
  rep(i, n)
  {
    rect[i].p1.x = point[i].x;
    rect[i].p1.y = point[i].y;
    rect[i].p2.x = point[i].x + 1;
    rect[i].p2.y = point[i].y + 1;
  }
  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }
  }
}

int real_real_real_a[MAX_N], real_real_real_b[MAX_N], real_real_real_c[MAX_N], real_real_real_d[MAX_N];
int real_real_real_maxScore = -1;

int a2[100][MAX_N], b2[100][MAX_N], c2[100][MAX_N], d2[100][MAX_N];
int a4[100][MAX_N], b4[100][MAX_N], c4[100][MAX_N], d4[100][MAX_N];
int maxScore4[100] = {};

inline void Ui_Tei()
{
  clock_t start, end;
  rep(ui_tei, 5)
  {

    // 初期解
    // 左上(x,y)、右下(x+1,y+1)
    rep(i, n)
    {
      rect[i].p1.x = point[i].x; rect[i].p2.x = point[i].x + 1;
      rect[i].p1.y = point[i].y; rect[i].p2.y = point[i].y + 1;
    }

    int T = 5;
    rep(_, T)
    {
      start = clock();

      // 初期スコア計算
      maxScore = calc(-1);
      real_maxScore = maxScore;
      rep(i, n)
      {
        real_a[i] = rect[i].p1.x;
        real_b[i] = rect[i].p1.y;
        real_c[i] = rect[i].p2.x;
        real_d[i] = rect[i].p2.y;
      }

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
          Slide(ite);
        }
        else if (now_time > 2.0 / T && mode == 2) {
          int ite = Rand() % n;
          aspectChange(ite);
        }
        else if (mode == -3) {
          int ite = Rand() % n;
          fourChange(ite, temp);
        }
      }

      // 焼きなまし戻す
      maxScore = real_maxScore;
      rep(i, n)
      {
        rect[i].p1.x = real_a[i];
        rect[i].p1.y = real_b[i];
        rect[i].p2.x = real_c[i];
        rect[i].p2.y = real_d[i];
      }
      calc(-1);

      if (maxScore > real_real_maxScore) {
        real_real_maxScore = maxScore;
        rep(i, n)
        {
          real_real_a[i] = rect[i].p1.x;
          real_real_b[i] = rect[i].p1.y;
          real_real_c[i] = rect[i].p2.x;
          real_real_d[i] = rect[i].p2.y;
        }
      }
    }

    // real_real_maxScore戻す
    maxScore = real_real_maxScore;
    real_real_maxScore = 0;
    rep(i, n)
    {
      rect[i].p1.x = real_real_a[i];
      rect[i].p1.y = real_real_b[i];
      rect[i].p2.x = real_real_c[i];
      rect[i].p2.y = real_real_d[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
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
        Slide(ite);
      }
      else if (mode == 2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
        int ite = Rand() % n;
        aspectChange(ite);
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
    rep(i, n)
    {
      rect[i].p1.x = real_a[i];
      rect[i].p1.y = real_b[i];
      rect[i].p2.x = real_c[i];
      rect[i].p2.y = real_d[i];
    }
    calc(-1);

    if (maxScore > ui_tei_maxScore) {
      ui_tei_maxScore = maxScore;
      rep(i, n)
      {
        ui_tei_a[i] = rect[i].p1.x;
        ui_tei_b[i] = rect[i].p1.y;
        ui_tei_c[i] = rect[i].p2.x;
        ui_tei_d[i] = rect[i].p2.y;
      }
    }
  }

  // 元に戻しておく
  maxScore = 0;
  real_maxScore = 0;
  real_real_maxScore = 0;
  rep(i, n)
  {
    rect[i].p1.x = point[i].x; rect[i].p2.x = point[i].x + 1;
    rect[i].p1.y = point[i].y; rect[i].p2.y = point[i].y + 1;
    real_a[i] = rect[i].p1.x; real_b[i] = rect[i].p1.y; real_c[i] = rect[i].p2.x; real_d[i] = rect[i].p2.y;
    real_real_a[i] = rect[i].p1.x; real_real_b[i] = rect[i].p1.y; real_real_c[i] = rect[i].p2.x; real_real_d[i] = rect[i].p2.y;
  }
}

int solve(int teisyutu, int fileNum)
{
  auto startClock = system_clock::now();
  clock_t start, end;
  clock_t real_start = clock();

  nyuuryokuInit(fileNum);

  sortInit();


  rep(allLoop, allLoopTimes)
  {

    Ui_Tei();

    // ui_tei_maxScore戻す
    maxScore = ui_tei_maxScore;
    rep(i, n)
    {
      rect[i].p1.x = ui_tei_a[i];
      rect[i].p1.y = ui_tei_b[i];
      rect[i].p2.x = ui_tei_c[i];
      rect[i].p2.y = ui_tei_d[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    rep(i, n)
    {
      real_a[i] = rect[i].p1.x;
      real_b[i] = rect[i].p1.y;
      real_c[i] = rect[i].p2.x;
      real_d[i] = rect[i].p2.y;
    }


    int oya = 1;

    rep(asai, oya)
    {
      rep(j, n)
      {
        a2[asai][j] = rect[j].p1.x;
        b2[asai][j] = rect[j].p1.y;
        c2[asai][j] = rect[j].p2.x;
        d2[asai][j] = rect[j].p2.y;
      }
    }

    // 焼きなまし2
    // 筋のいいやつを追う
    int T = 250 / allLoopTimes;
    rep(_, T)
    {
      rep(i, 6) modeCount[i] = 0;

      int TT = 1;
      rep(asai, TT)
      {
        int kiyoshi = asai % oya;
        rep(i, n)
        {
          rect[i].p1.x = a2[kiyoshi][i];
          rect[i].p1.y = b2[kiyoshi][i];
          rect[i].p2.x = c2[kiyoshi][i];
          rect[i].p2.y = d2[kiyoshi][i];
        }

        // 初期スコア計算
        maxScore = calc(-1);
        real_maxScore = maxScore;

        rep(i, n)
        {
          real_a[i] = rect[i].p1.x;
          real_b[i] = rect[i].p1.y;
          real_c[i] = rect[i].p2.x;
          real_d[i] = rect[i].p2.y;
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
            Slide(ite);
          }
          else if (mode == -2 && kouhan && now_time > 2.0 / T) { // ランダムにアスペクト比を変更
            int ite = Rand() % n;
            aspectChange(ite);
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
            int tmpScore = calc(-1);
            // cout << fixed << setprecision(10) << tmpScore - maxScore << endl;
            maxScore = tmpScore;
          }

        }

        // 焼きなまし戻す
        maxScore = real_maxScore;
        rep(i, n)
        {
          rect[i].p1.x = real_a[i];
          rect[i].p1.y = real_b[i];
          rect[i].p2.x = real_c[i];
          rect[i].p2.y = real_d[i];
        }
        calc(-1);
        if (maxScore > real_real_maxScore) {
          real_real_maxScore = maxScore;
          rep(i, n)
          {
            real_real_a[i] = rect[i].p1.x;
            real_real_b[i] = rect[i].p1.y;
            real_real_c[i] = rect[i].p2.x;
            real_real_d[i] = rect[i].p2.y;
          }
        }

        // ビームサーチの次の種にする
        maxScore4[asai] = maxScore;
        rep(i, n)
        {
          a4[asai][i] = rect[i].p1.x;
          b4[asai][i] = rect[i].p1.y;
          c4[asai][i] = rect[i].p2.x;
          d4[asai][i] = rect[i].p2.y;
        }

        // cout << loop << ' ';
        // cout << temp << ' ' << temp2 << endl;
        // cout << real_real_maxScore << endl;
      }

      // 次の世代に継承
      vector<P> vBeam;
      rep(asai, TT) vBeam.emplace_back(P(maxScore4[asai], asai));
      sort(vBeam.begin(), vBeam.end(), greater<P>());


      rep(ii, oya)
      {
        int i = vBeam[ii].second;
        rep(j, n)
        {
          a2[ii][j] = a4[i][j];
          b2[ii][j] = b4[i][j];
          c2[ii][j] = c4[i][j];
          d2[ii][j] = d4[i][j];
        }
      }

      /*
      rep(i,6) cout << modeCount[i] << ' ';
      cout << endl;
      */


      // 提出時以下は消す
      if (teisyutu == 0 && _ % 10 == 0) {
        cout << "_ = " << _;
        cout << ", vBeam[0] = (" << vBeam[0].first << ", " << vBeam[0].second << ")" << endl;
        // cout << ", vBeam[1] = (" << vBeam[1].first << ", " << vBeam[1].second << ")" << endl;
      }


      // FileKakikomi(fileNum);

      // エスケープ
      end = clock();
      if (((double)end - real_start) / CLOCKS_PER_SEC > realTL) break;
    }

    // real_real_maxScore戻す
    maxScore = real_real_maxScore;
    rep(i, n)
    {
      rect[i].p1.x = real_real_a[i];
      rect[i].p1.y = real_real_b[i];
      rect[i].p2.x = real_real_c[i];
      rect[i].p2.y = real_real_d[i];
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
      rep(i, n)
      {
        real_real_real_a[i] = rect[i].p1.x;
        real_real_real_b[i] = rect[i].p1.y;
        real_real_real_c[i] = rect[i].p2.x;
        real_real_real_d[i] = rect[i].p2.y;
      }
    }


    // すべて白紙にリセットする
    maxScore = 0;
    real_maxScore = 0;
    real_real_maxScore = 0;
    rep(i, n)
    {
      rect[i].p1.x = point[i].x; rect[i].p2.x = point[i].x + 1;
      rect[i].p1.y = point[i].y; rect[i].p2.y = point[i].y + 1;
      real_a[i] = rect[i].p1.x; real_b[i] = rect[i].p1.y; real_c[i] = rect[i].p2.x; real_d[i] = rect[i].p2.y;
      real_real_a[i] = rect[i].p1.x; real_real_b[i] = rect[i].p1.y; real_real_c[i] = rect[i].p2.x; real_real_d[i] = rect[i].p2.y;
    }
  }


  // real_real_real_maxScore戻す
  maxScore = real_real_real_maxScore;
  rep(i, n)
  {
    rect[i].p1.x = real_real_real_a[i];
    rect[i].p1.y = real_real_real_b[i];
    rect[i].p2.x = real_real_real_c[i];
    rect[i].p2.y = real_real_real_d[i];
  }
  calc(-1);

  // 最終出力
  if (teisyutu) {
    rep(i, n)
    {
      cout << rect[i].p1.x << ' ' << rect[i].p1.y << ' ' << rect[i].p2.x << ' ' << rect[i].p2.y << endl;
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

inline void AllClear()
{
  n = 0;
  maxScore = -1;
  real_maxScore = -1;
  pSum = 0;
  real_real_maxScore = -1;
  ui_tei_maxScore = -1;
  real_real_real_maxScore = -1;
  rep(i, MAX_N)
  {
    point[i].x = 0, point[i].y = 0, target_sizes[i] = 0;
    rect[i].p1.x = 0, rect[i].p1.y = 0, rect[i].p2.x = 0, rect[i].p2.y = 0;
    area_sizes[i] = 0;
    real_a[i] = 0, real_b[i] = 0, real_c[i] = 0, real_d[i] = 0;
    p[i] = 0;
    sort_x[i] = 0, sort_y[i] = 0;
    arg_sort_x[i] = 0, arg_sort_y[i] = 0;
    real_real_a[i] = 0, real_real_b[i] = 0, real_real_c[i] = 0, real_real_d[i] = 0;
    ui_tei_a[i] = 0, ui_tei_b[i] = 0, ui_tei_c[i] = 0, ui_tei_d[i] = 0;
    real_real_real_a[i] = 0, real_real_real_b[i] = 0, real_real_real_c[i] = 0, real_real_real_d[i] = 0;
  }
}

int main()
{
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
      rep(i, 1000)
      {
        srep(i, 0, 50)
        {
          rep(j, 10)
          {
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

      rep(k, m)
      {
        rep(l, 4)
        {
          haipara_oya = arrHaiparaOya[k];
          if (l == 0) haipara_TT = haipara_oya;
          if (l == 1) haipara_TT = haipara_oya + 1;
          if (l == 2) haipara_TT = haipara_oya * 2;
          if (l == 3) haipara_TT = haipara_oya * 4;

          ll sum = 0;
          srep(i, 0, 50)
          {
            rep(j, 1)
            {
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


