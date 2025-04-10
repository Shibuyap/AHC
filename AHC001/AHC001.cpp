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
#include<queue>
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
int allLoopTimes = 1;
int n;
int x[MAX_N], y[MAX_N], r[MAX_N];
int a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
int s[MAX_N];
int real_a[MAX_N], real_b[MAX_N], real_c[MAX_N], real_d[MAX_N];

inline void nyuuryokuInit(int fileNum)
{
  // 入力
  string fileNameIfs = to_string(fileNum);
  fileNameIfs += ".txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs.is_open()) { // 標準入力する
    cin >> n;
    rep(i, n) cin >> x[i] >> y[i] >> r[i];
  }
  else { // ファイル入力する
    ifs >> n;
    rep(i, n) ifs >> x[i] >> y[i] >> r[i];
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
    ofs << a[i] << ' ' << b[i] << ' ' << c[i] << ' ' << d[i] << endl;
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
    ofs << a[i] << ' ' << b[i] << ' ' << c[i] << ' ' << d[i] << endl;
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
      s[i] = (c[i] - a[i]) * (d[i] - b[i]);
      p[i] = 1.0 - (1.0 - (double)min(r[i], s[i]) / (double)max(r[i], s[i])) * (1.0 - (double)min(r[i], s[i]) / (double)max(r[i], s[i]));
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
    s[ite] = (c[ite] - a[ite]) * (d[ite] - b[ite]);
    p[ite] = 1.0 - (1.0 - (double)min(r[ite], s[ite]) / (double)max(r[ite], s[ite])) * (1.0 - (double)min(r[ite], s[ite]) / (double)max(r[ite], s[ite]));
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
  if (a[i] <= a[j] && a[j] < c[i]) cnt++;
  else if (a[j] <= a[i] && a[i] < c[j]) cnt++;
  if (b[i] <= b[j] && b[j] < d[i]) cnt++;
  else if (b[j] <= b[i] && b[i] < d[j]) cnt++;
  return cnt == 2;
}

int sort_x[MAX_N], sort_y[MAX_N];
int arg_sort_x[MAX_N], arg_sort_y[MAX_N];
inline void sortInit()
{
  vector<P> v;
  rep(i, n)
  {
    v.emplace_back(P(x[i], i));
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
    v.emplace_back(P(y[i], i));
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
      if (a[i] < 0 || 10000 < a[i]) return 0;
      if (b[i] < 0 || 10000 < b[i]) return 0;
      if (c[i] < 0 || 10000 < c[i]) return 0;
      if (d[i] < 0 || 10000 < d[i]) return 0;
    }
    rep(i, n)
    {
      if (c[i] <= a[i]) return 0;
      if (d[i] <= b[i]) return 0;
    }
    rep(i, n)
    {
      if (x[i] < a[i] || c[i] <= x[i]) return 0;
      if (y[i] < b[i] || d[i] <= y[i]) return 0;
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
    if (a[ite] < 0 || 10000 < a[ite]) return 0;
    if (b[ite] < 0 || 10000 < b[ite]) return 0;
    if (c[ite] < 0 || 10000 < c[ite]) return 0;
    if (d[ite] < 0 || 10000 < d[ite]) return 0;
    if (c[ite] <= a[ite]) return 0;
    if (d[ite] <= b[ite]) return 0;
    if (x[ite] < a[ite] || c[ite] <= x[ite]) return 0;
    if (y[ite] < b[ite] || d[ite] <= y[ite]) return 0;
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
      if (a[i] < 0 || 10000 < a[i]) return 0;
      if (b[i] < 0 || 10000 < b[i]) return 0;
      if (c[i] < 0 || 10000 < c[i]) return 0;
      if (d[i] < 0 || 10000 < d[i]) return 0;
    }
    rep(i, n)
    {
      if (c[i] <= a[i]) return 0;
      if (d[i] <= b[i]) return 0;
    }
    rep(i, n)
    {
      if (x[i] < a[i] || c[i] <= x[i]) return 0;
      if (y[i] < b[i] || d[i] <= y[i]) return 0;
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
    if (a[ite] < 0 || 10000 < a[ite]) return 0;
    if (b[ite] < 0 || 10000 < b[ite]) return 0;
    if (c[ite] < 0 || 10000 < c[ite]) return 0;
    if (d[ite] < 0 || 10000 < d[ite]) return 0;
    if (c[ite] <= a[ite]) return 0;
    if (d[ite] <= b[ite]) return 0;
    if (x[ite] < a[ite] || c[ite] <= x[ite]) return 0;
    if (y[ite] < b[ite] || d[ite] <= y[ite]) return 0;
    int argX = arg_sort_x[ite];
    int nowLeft = b[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (b[i] <= nowLeft) {
        nowLeft = max(nowLeft, d[i]);
        if (nowLeft >= d[ite]) break;
      }
    }
    nowLeft = b[ite];
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) return 0;
      if (b[i] <= nowLeft) {
        nowLeft = max(nowLeft, d[i]);
        if (nowLeft >= d[ite]) break;
      }
    }
  }
  return 1;
}

int vExtend[4];
inline void hukuramashi(int ite)
{
  vExtend[0] = x[ite];
  vExtend[1] = y[ite];
  vExtend[2] = x[ite] + 1;
  vExtend[3] = y[ite] + 1;

  int flagTateYoko = Rand() % 2;
  if (flagTateYoko == 0) {
    vExtend[0] = 0;
    vExtend[2] = 10000;
    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (b[i] <= vExtend[1] && vExtend[1] < d[i]) flagKasanari = 1;
      if (vExtend[1] <= b[i] && b[i] < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtend[0] = max(vExtend[0], c[i]);
        }
        else {
          vExtend[2] = min(vExtend[2], a[i]);
        }
        break;
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (b[i] <= vExtend[1] && vExtend[1] < d[i]) flagKasanari = 1;
      if (vExtend[1] <= b[i] && b[i] < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtend[0] = max(vExtend[0], c[i]);
        }
        else {
          vExtend[2] = min(vExtend[2], a[i]);
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
      if (a[i] <= vExtend[0] && vExtend[0] < c[i]) flagKasanari = 1;
      if (vExtend[0] <= a[i] && a[i] < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtend[1] = max(vExtend[1], d[i]);
        }
        else {
          vExtend[3] = min(vExtend[3], b[i]);
        }
        if (a[i] <= nowLeft) {
          nowLeft = max(nowLeft, c[i]);
          if (vExtend[2] <= nowLeft) break;
        }
      }
    }

    nowLeft = vExtend[0];
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (a[i] <= vExtend[0] && vExtend[0] < c[i]) flagKasanari = 1;
      if (vExtend[0] <= a[i] && a[i] < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtend[1] = max(vExtend[1], d[i]);
        }
        else {
          vExtend[3] = min(vExtend[3], b[i]);
        }
        if (a[i] <= nowLeft) {
          nowLeft = max(nowLeft, c[i]);
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
      if (a[i] <= vExtend[0] && vExtend[0] < c[i]) flagKasanari = 1;
      if (vExtend[0] <= a[i] && a[i] < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtend[1] = max(vExtend[1], d[i]);
        }
        else {
          vExtend[3] = min(vExtend[3], b[i]);
        }
        break;
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      int flagKasanari = 0;
      if (a[i] <= vExtend[0] && vExtend[0] < c[i]) flagKasanari = 1;
      if (vExtend[0] <= a[i] && a[i] < vExtend[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtend[1] = max(vExtend[1], d[i]);
        }
        else {
          vExtend[3] = min(vExtend[3], b[i]);
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
      if (b[i] <= vExtend[1] && vExtend[1] < d[i]) flagKasanari = 1;
      if (vExtend[1] <= b[i] && b[i] < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtend[0] = max(vExtend[0], c[i]);
        }
        else {
          vExtend[2] = min(vExtend[2], a[i]);
        }
        if (b[i] <= nowLeft) {
          nowLeft = max(nowLeft, d[i]);
          if (vExtend[3] <= nowLeft) break;
        }
      }
    }
    nowLeft = vExtend[1];
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      int flagKasanari = 0;
      if (b[i] <= vExtend[1] && vExtend[1] < d[i]) flagKasanari = 1;
      if (vExtend[1] <= b[i] && b[i] < vExtend[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtend[0] = max(vExtend[0], c[i]);
        }
        else {
          vExtend[2] = min(vExtend[2], a[i]);
        }
        if (b[i] <= nowLeft) {
          nowLeft = max(nowLeft, d[i]);
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
    if (S <= r[ite]) break;
    if (shuf[i] == 0) {
      int ma = r[ite] / (vExtend[3] - vExtend[1]) + yure;
      int diff = (vExtend[2] - vExtend[0]) - ma;
      int capa = x[ite] - vExtend[0];
      if (capa >= diff) {
        vExtend[0] += diff;
      }
      else {
        vExtend[0] += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = r[ite] / (vExtend[2] - vExtend[0]) + yure;
      int diff = (vExtend[3] - vExtend[1]) - ma;
      int capa = y[ite] - vExtend[1];
      if (capa >= diff) {
        vExtend[1] += diff;
      }
      else {
        vExtend[1] += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = r[ite] / (vExtend[3] - vExtend[1]) + yure;
      int diff = (vExtend[2] - vExtend[0]) - ma;
      int capa = vExtend[2] - (x[ite] + 1);
      if (capa >= diff) {
        vExtend[2] -= diff;
      }
      else {
        vExtend[2] -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = r[ite] / (vExtend[2] - vExtend[0]) + yure;
      int diff = (vExtend[3] - vExtend[1]) - ma;
      int capa = vExtend[3] - (y[ite] + 1);
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
  if (x[ite] < vExtend[0] || vExtend[2] <= x[ite]) ng = 1;
  if (y[ite] < vExtend[1] || vExtend[3] <= y[ite]) ng = 1;

  if (ng) {
    vExtend[0] = x[ite]; vExtend[2] = x[ite] + 1;
    vExtend[1] = y[ite]; vExtend[3] = y[ite] + 1;
  }

}

inline void Extend(int ite, double temp)
{
  int keepA = a[ite];
  int keepB = b[ite];
  int keepC = c[ite];
  int keepD = d[ite];

  hukuramashi(ite);
  a[ite] = vExtend[0];
  b[ite] = vExtend[1];
  c[ite] = vExtend[2];
  d[ite] = vExtend[3];

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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    // 元に戻す
    a[ite] = keepA;
    b[ite] = keepB;
    c[ite] = keepC;
    d[ite] = keepD;
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
    a[ite] = x[ite];
    b[ite] = y[ite];
    c[ite] = x[ite] + 1;
    d[ite] = y[ite] + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
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
    a[ite] = x[ite];
    b[ite] = y[ite];
    c[ite] = x[ite] + 1;
    d[ite] = y[ite] + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
    }
  }
}

inline void AnaWoAkeru(int hole)
{
  int ite = Rand() % n;
  vector<int> keep;
  keep.emplace_back(ite);
  a[ite] -= 100;
  b[ite] -= 100;
  c[ite] += 100;
  d[ite] += 100;
  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) keep.emplace_back(i);
  }
  int keepSize = keep.size();
  rep(i, keepSize)
  {
    a[keep[i]] = x[keep[i]];
    b[keep[i]] = y[keep[i]];
    c[keep[i]] = x[keep[i]] + 1;
    d[keep[i]] = y[keep[i]] + 1;
  }

  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
    }
  }
}



int vExtendKing[4];
inline void hukuramashiKing(int ite)
{
  vExtendKing[0] = max(0, (int)(x[ite] - Rand() % 1000));
  vExtendKing[1] = max(0, (int)(y[ite] - Rand() % 1000));
  vExtendKing[2] = min(10000, (int)(x[ite] + 1 + Rand() % 1000));
  vExtendKing[3] = min(10000, (int)(y[ite] + 1 + Rand() % 1000));

  int tateyoko = Rand() % 2;

  if (tateyoko == 0) {
    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (x[i] == x[ite]) continue;
      int flagKasanari = 0;
      if (y[i] <= vExtendKing[1] && vExtendKing[1] < y[i] + 1) flagKasanari = 1;
      if (vExtendKing[1] <= y[i] && y[i] < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtendKing[0] = max(vExtendKing[0], x[i] + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], x[i]);
        }
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (x[i] == x[ite]) continue;
      int flagKasanari = 0;
      if (y[i] <= vExtendKing[1] && vExtendKing[1] < y[i] + 1) flagKasanari = 1;
      if (vExtendKing[1] <= y[i] && y[i] < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtendKing[0] = max(vExtendKing[0], x[i] + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], x[i]);
        }
      }
    }

    int argY = arg_sort_y[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (y[i] == y[ite]) continue;
      int flagKasanari = 0;
      if (x[i] <= vExtendKing[0] && vExtendKing[0] < x[i] + 1) flagKasanari = 1;
      if (vExtendKing[0] <= x[i] && x[i] < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtendKing[1] = max(vExtendKing[1], y[i] + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], y[i]);
        }
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (y[i] == y[ite]) continue;
      int flagKasanari = 0;
      if (x[i] <= vExtendKing[0] && vExtendKing[0] < x[i] + 1) flagKasanari = 1;
      if (vExtendKing[0] <= x[i] && x[i] < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtendKing[1] = max(vExtendKing[1], y[i] + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], y[i]);
        }
      }
    }
  }
  else {
    int argY = arg_sort_y[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (y[i] == y[ite]) continue;
      int flagKasanari = 0;
      if (x[i] <= vExtendKing[0] && vExtendKing[0] < x[i] + 1) flagKasanari = 1;
      if (vExtendKing[0] <= x[i] && x[i] < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtendKing[1] = max(vExtendKing[1], y[i] + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], y[i]);
        }
      }
    }
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (y[i] == y[ite]) continue;
      int flagKasanari = 0;
      if (x[i] <= vExtendKing[0] && vExtendKing[0] < x[i] + 1) flagKasanari = 1;
      if (vExtendKing[0] <= x[i] && x[i] < vExtendKing[2]) flagKasanari = 1;
      if (flagKasanari) {
        if (y[i] <= y[ite]) {
          vExtendKing[1] = max(vExtendKing[1], y[i] + 1);
        }
        else {
          vExtendKing[3] = min(vExtendKing[3], y[i]);
        }
      }
    }

    int argX = arg_sort_x[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (x[i] == x[ite]) continue;
      int flagKasanari = 0;
      if (y[i] <= vExtendKing[1] && vExtendKing[1] < y[i] + 1) flagKasanari = 1;
      if (vExtendKing[1] <= y[i] && y[i] < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtendKing[0] = max(vExtendKing[0], x[i] + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], x[i]);
        }
      }
    }
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (x[i] == x[ite]) continue;
      int flagKasanari = 0;
      if (y[i] <= vExtendKing[1] && vExtendKing[1] < y[i] + 1) flagKasanari = 1;
      if (vExtendKing[1] <= y[i] && y[i] < vExtendKing[3]) flagKasanari = 1;
      if (flagKasanari) {
        if (x[i] <= x[ite]) {
          vExtendKing[0] = max(vExtendKing[0], x[i] + 1);
        }
        else {
          vExtendKing[2] = min(vExtendKing[2], x[i]);
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
    if (S <= r[ite]) break;
    if (shuf[i] == 0) {
      int ma = r[ite] / (vExtendKing[3] - vExtendKing[1]) + yure;
      int diff = (vExtendKing[2] - vExtendKing[0]) - ma;
      if (diff < 0) diff = 0;
      int capa = x[ite] - vExtendKing[0];
      if (capa >= diff) {
        vExtendKing[0] += diff;
      }
      else {
        vExtendKing[0] += capa;
      }
    }
    else if (shuf[i] == 1) {
      int ma = r[ite] / (vExtendKing[2] - vExtendKing[0]) + yure;
      int diff = (vExtendKing[3] - vExtendKing[1]) - ma;
      if (diff < 0) diff = 0;
      int capa = y[ite] - vExtendKing[1];
      if (capa >= diff) {
        vExtendKing[1] += diff;
      }
      else {
        vExtendKing[1] += capa;
      }
    }
    else if (shuf[i] == 2) {
      int ma = r[ite] / (vExtendKing[3] - vExtendKing[1]) + yure;
      int diff = (vExtendKing[2] - vExtendKing[0]) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing[2] - (x[ite] + 1);
      if (capa >= diff) {
        vExtendKing[2] -= diff;
      }
      else {
        vExtendKing[2] -= capa;
      }
    }
    else if (shuf[i] == 3) {
      int ma = r[ite] / (vExtendKing[2] - vExtendKing[0]) + yure;
      int diff = (vExtendKing[3] - vExtendKing[1]) - ma;
      if (diff < 0) diff = 0;
      int capa = vExtendKing[3] - (y[ite] + 1);
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
  if (x[ite] < vExtendKing[0] || vExtendKing[2] <= x[ite]) ng = 1;
  if (y[ite] < vExtendKing[1] || vExtendKing[3] <= y[ite]) ng = 1;

  if (ng) {
    vExtendKing[0] = x[ite]; vExtendKing[2] = x[ite] + 1;
    vExtendKing[1] = y[ite]; vExtendKing[3] = y[ite] + 1;
  }

}

inline void ExtendKing(int ite)
{
  int keepA = a[ite];
  int keepB = b[ite];
  int keepC = c[ite];
  int keepD = d[ite];

  hukuramashiKing(ite);
  a[ite] = vExtendKing[0];
  b[ite] = vExtendKing[1];
  c[ite] = vExtendKing[2];
  d[ite] = vExtendKing[3];


  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      a[i] = x[i];
      b[i] = y[i];
      c[i] = x[i] + 1;
      d[i] = y[i] + 1;
    }
  }

  maxScore = calc(-1);
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
    }
  }
}

inline void oneChange(int ite, double temp)
{
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int abcd = Rand() % 4;

  if (abcd == 0) a[ite] += diff;
  if (abcd == 1) b[ite] += diff;
  if (abcd == 2) c[ite] += diff;
  if (abcd == 3) d[ite] += diff;

  if (isOK(ite) == 0) {
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    // 元に戻す
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
    calc(ite);
  }
}

inline void fourChange(int ite, double temp)
{
  int diffA = Rand() % 101 - 50;
  int diffB = Rand() % 101 - 50;
  int diffC = Rand() % 101 - 50;
  int diffD = Rand() % 101 - 50;

  a[ite] += diffA;
  b[ite] += diffB;
  c[ite] += diffC;
  d[ite] += diffD;

  if (isOK(ite) == 0) {
    a[ite] -= diffA;
    b[ite] -= diffB;
    c[ite] -= diffC;
    d[ite] -= diffD;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    // 元に戻す
    a[ite] -= diffA;
    b[ite] -= diffB;
    c[ite] -= diffC;
    d[ite] -= diffD;
    calc(ite);
  }
}

inline void Slide(int ite)
{
  int diff = 0;
  while (diff == 0) diff = Rand() % 101 - 50;
  int ab = Rand() % 2;

  if (ab == 0) {
    a[ite] += diff;
    c[ite] += diff;
  }
  if (ab == 1) {
    b[ite] += diff;
    d[ite] += diff;
  }

  if (isOK(ite) == 0) {
    if (ab == 0) {
      a[ite] -= diff;
      c[ite] -= diff;
    }
    if (ab == 1) {
      b[ite] -= diff;
      d[ite] -= diff;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    // 元に戻す
    if (ab == 0) {
      a[ite] -= diff;
      c[ite] -= diff;
    }
    if (ab == 1) {
      b[ite] -= diff;
      d[ite] -= diff;
    }
    calc(ite);
  }
}

inline void aspectChange(int ite)
{
  int yokoRatio = Rand() % 9 + 1; // 1 ~ 9;
  int tateRatio = 10 - yokoRatio;

  int S = yokoRatio * tateRatio;
  int mul = sqrt(r[ite] / S);
  if (mul == 0) return;

  int yoko = yokoRatio * mul;
  int tate = tateRatio * mul;

  int keepA = a[ite];
  int keepB = b[ite];
  int keepC = c[ite];
  int keepD = d[ite];

  int leftA = max(0, x[ite] - (yoko - 1));
  int rightA = min(x[ite], 10000 - yoko);
  int rangeA = rightA - leftA + 1;
  if (rangeA < 1) return;

  int leftB = max(0, y[ite] - (tate - 1));
  int rightB = min(y[ite], 10000 - tate);
  int rangeB = rightB - leftB + 1;
  if (rangeB < 1) return;

  a[ite] = Rand() % rangeA + leftA;
  c[ite] = a[ite] + rangeA;
  b[ite] = Rand() % rangeB + leftB;
  d[ite] = b[ite] + rangeB;

  if (isOK(ite) == 0) {
    a[ite] = keepA;
    b[ite] = keepB;
    c[ite] = keepC;
    d[ite] = keepD;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    // 元に戻す
    a[ite] = keepA;
    b[ite] = keepB;
    c[ite] = keepC;
    d[ite] = keepD;
    calc(ite);
  }
}

inline int selfNg(int ite)
{
  if (a[ite] < 0 || 10000 < a[ite]) return 1;
  if (b[ite] < 0 || 10000 < b[ite]) return 1;
  if (c[ite] < 0 || 10000 < c[ite]) return 1;
  if (d[ite] < 0 || 10000 < d[ite]) return 1;
  if (c[ite] <= a[ite]) return 1;
  if (d[ite] <= b[ite]) return 1;
  if (x[ite] < a[ite] || c[ite] <= x[ite]) return 1;
  if (y[ite] < b[ite] || d[ite] <= y[ite]) return 1;
  return 0;
}

inline int dokasuOK(int ite, int abcd)
{
  rep(i, n)
  {
    if (i == ite) continue;
    if (kasanarihantei(i, ite)) {
      if (abcd == 0) c[i] = a[ite];
      if (abcd == 1) d[i] = b[ite];
      if (abcd == 2) a[i] = c[ite];
      if (abcd == 3) b[i] = d[ite];

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
    int nowLeft = b[ite];
    int nowRight = d[ite];
    drep(ii, argX)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (a[ite] <= x[i]) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (b[i] <= nowLeft) {
        nowLeft = max(nowLeft, d[i]);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= d[i]) {
        nowRight = min(nowRight, b[i]);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 1) {
    int argY = arg_sort_y[ite];
    int nowLeft = a[ite];
    int nowRight = c[ite];
    drep(ii, argY)
    {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (b[ite] <= y[i]) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (a[i] <= nowLeft) {
        nowLeft = max(nowLeft, c[i]);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= c[i]) {
        nowRight = min(nowRight, a[i]);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 2) {
    int argX = arg_sort_x[ite];
    int nowLeft = b[ite];
    int nowRight = d[ite];
    srep(ii, argX + 1, n)
    {
      int i = sort_x[ii];
      if (kasanarihantei(i, ite)) {
        if (x[i] < c[ite]) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (b[i] <= nowLeft) {
        nowLeft = max(nowLeft, d[i]);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= d[i]) {
        nowRight = min(nowRight, b[i]);
        if (nowLeft >= nowRight) break;
      }
    }
  }
  if (abcd == 3) {
    int argY = arg_sort_y[ite];
    int nowLeft = a[ite];
    int nowRight = c[ite];
    srep(ii, argY + 1, n)
    {
      int i = sort_y[ii];
      if (kasanarihantei(i, ite)) {
        if (y[i] < d[ite]) {
          arrKasanari[0] = -1;
          kasanariCount = 1;
          return;
        }
        arrKasanari[kasanariCount] = i;
        kasanariCount++;
      }
      if (a[i] <= nowLeft) {
        nowLeft = max(nowLeft, c[i]);
        if (nowLeft >= nowRight) break;
      }
      if (nowRight <= c[i]) {
        nowRight = min(nowRight, a[i]);
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

  if (abcd == 0) a[ite] += diff;
  if (abcd == 1) b[ite] += diff;
  if (abcd == 2) c[ite] += diff;
  if (abcd == 3) d[ite] += diff;

  if (selfNg(ite)) {
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
    return;
  }

  kasanaritati(ite, abcd);
  int vn = kasanariCount;

  if (vn > 0 && arrKasanari[0] == -1) {
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
    return;
  }


  rep(i, vn)
  {
    keepvA[i] = a[arrKasanari[i]];
    keepvB[i] = b[arrKasanari[i]];
    keepvC[i] = c[arrKasanari[i]];
    keepvD[i] = d[arrKasanari[i]];
  }

  int ok = 1;
  rep(i, vn)
  {
    if (abcd == 0) c[arrKasanari[i]] = a[ite];
    if (abcd == 1) d[arrKasanari[i]] = b[ite];
    if (abcd == 2) a[arrKasanari[i]] = c[ite];
    if (abcd == 3) b[arrKasanari[i]] = d[ite];
    if (selfNg(arrKasanari[i])) ok = 0;
  }

  if (ok == 0) {
    rep(i, vn)
    {
      a[arrKasanari[i]] = keepvA[i];
      b[arrKasanari[i]] = keepvB[i];
      c[arrKasanari[i]] = keepvC[i];
      d[arrKasanari[i]] = keepvD[i];
    }
    // 元に戻す
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
      }
    }
  }
  else {
    rep(i, vn)
    {
      a[arrKasanari[i]] = keepvA[i];
      b[arrKasanari[i]] = keepvB[i];
      c[arrKasanari[i]] = keepvC[i];
      d[arrKasanari[i]] = keepvD[i];
      calc(arrKasanari[i]);
    }
    // 元に戻す
    if (abcd == 0) a[ite] -= diff;
    if (abcd == 1) b[ite] -= diff;
    if (abcd == 2) c[ite] -= diff;
    if (abcd == 3) d[ite] -= diff;
    calc(ite);
  }
}

inline void shokiInit()
{
  rep(i, n)
  {
    a[i] = x[i];
    b[i] = y[i];
    c[i] = x[i] + 1;
    d[i] = y[i] + 1;
  }
  int tmpScore = calc(-1);

  maxScore = tmpScore;
  if (maxScore > real_maxScore) {
    real_maxScore = maxScore;
    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
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
      a[i] = x[i]; c[i] = x[i] + 1;
      b[i] = y[i]; d[i] = y[i] + 1;
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
        real_a[i] = a[i];
        real_b[i] = b[i];
        real_c[i] = c[i];
        real_d[i] = d[i];
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
        a[i] = real_a[i];
        b[i] = real_b[i];
        c[i] = real_c[i];
        d[i] = real_d[i];
      }
      calc(-1);

      if (maxScore > real_real_maxScore) {
        real_real_maxScore = maxScore;
        rep(i, n)
        {
          real_real_a[i] = a[i];
          real_real_b[i] = b[i];
          real_real_c[i] = c[i];
          real_real_d[i] = d[i];
        }
      }
    }

    // real_real_maxScore戻す
    maxScore = real_real_maxScore;
    real_real_maxScore = 0;
    rep(i, n)
    {
      a[i] = real_real_a[i];
      b[i] = real_real_b[i];
      c[i] = real_real_c[i];
      d[i] = real_real_d[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
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
      a[i] = real_a[i];
      b[i] = real_b[i];
      c[i] = real_c[i];
      d[i] = real_d[i];
    }
    calc(-1);

    if (maxScore > ui_tei_maxScore) {
      ui_tei_maxScore = maxScore;
      rep(i, n)
      {
        ui_tei_a[i] = a[i];
        ui_tei_b[i] = b[i];
        ui_tei_c[i] = c[i];
        ui_tei_d[i] = d[i];
      }
    }
  }

  // 元に戻しておく
  maxScore = 0;
  real_maxScore = 0;
  real_real_maxScore = 0;
  rep(i, n)
  {
    a[i] = x[i]; c[i] = x[i] + 1;
    b[i] = y[i]; d[i] = y[i] + 1;
    real_a[i] = a[i]; real_b[i] = b[i]; real_c[i] = c[i]; real_d[i] = d[i];
    real_real_a[i] = a[i]; real_real_b[i] = b[i]; real_real_c[i] = c[i]; real_real_d[i] = d[i];
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
      a[i] = ui_tei_a[i];
      b[i] = ui_tei_b[i];
      c[i] = ui_tei_c[i];
      d[i] = ui_tei_d[i];
    }
    calc(-1);

    // 初期スコア計算
    maxScore = calc(-1);
    real_maxScore = maxScore;

    rep(i, n)
    {
      real_a[i] = a[i];
      real_b[i] = b[i];
      real_c[i] = c[i];
      real_d[i] = d[i];
    }


    int oya = 1;

    rep(asai, oya)
    {
      rep(j, n)
      {
        a2[asai][j] = a[j];
        b2[asai][j] = b[j];
        c2[asai][j] = c[j];
        d2[asai][j] = d[j];
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
          a[i] = a2[kiyoshi][i];
          b[i] = b2[kiyoshi][i];
          c[i] = c2[kiyoshi][i];
          d[i] = d2[kiyoshi][i];
        }

        // 初期スコア計算
        maxScore = calc(-1);
        real_maxScore = maxScore;

        rep(i, n)
        {
          real_a[i] = a[i];
          real_b[i] = b[i];
          real_c[i] = c[i];
          real_d[i] = d[i];
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
          a[i] = real_a[i];
          b[i] = real_b[i];
          c[i] = real_c[i];
          d[i] = real_d[i];
        }
        calc(-1);
        if (maxScore > real_real_maxScore) {
          real_real_maxScore = maxScore;
          rep(i, n)
          {
            real_real_a[i] = a[i];
            real_real_b[i] = b[i];
            real_real_c[i] = c[i];
            real_real_d[i] = d[i];
          }
        }

        // ビームサーチの次の種にする
        maxScore4[asai] = maxScore;
        rep(i, n)
        {
          a4[asai][i] = a[i];
          b4[asai][i] = b[i];
          c4[asai][i] = c[i];
          d4[asai][i] = d[i];
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
      a[i] = real_real_a[i];
      b[i] = real_real_b[i];
      c[i] = real_real_c[i];
      d[i] = real_real_d[i];
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
        real_real_real_a[i] = a[i];
        real_real_real_b[i] = b[i];
        real_real_real_c[i] = c[i];
        real_real_real_d[i] = d[i];
      }
    }


    // すべて白紙にリセットする
    maxScore = 0;
    real_maxScore = 0;
    real_real_maxScore = 0;
    rep(i, n)
    {
      a[i] = x[i]; c[i] = x[i] + 1;
      b[i] = y[i]; d[i] = y[i] + 1;
      real_a[i] = a[i]; real_b[i] = b[i]; real_c[i] = c[i]; real_d[i] = d[i];
      real_real_a[i] = a[i]; real_real_b[i] = b[i]; real_real_c[i] = c[i]; real_real_d[i] = d[i];
    }
  }


  // real_real_real_maxScore戻す
  maxScore = real_real_real_maxScore;
  rep(i, n)
  {
    a[i] = real_real_real_a[i];
    b[i] = real_real_real_b[i];
    c[i] = real_real_real_c[i];
    d[i] = real_real_real_d[i];
  }
  calc(-1);

  // 最終出力
  if (teisyutu) {
    rep(i, n)
    {
      cout << a[i] << ' ' << b[i] << ' ' << c[i] << ' ' << d[i] << endl;
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
    x[i] = 0, y[i] = 0, r[i] = 0;
    a[i] = 0, b[i] = 0, c[i] = 0, d[i] = 0;
    s[i] = 0;
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
  int teisyutu = 1;

  if (teisyutu) {
    solve(teisyutu, 1121);
  }
  else {
    int mode = 0;
    if (mode == 0) { // コードテスト用
      solve(teisyutu, 1121);
    }

    if (mode == 1) { // スコア確認用
      rep(i, 1000)
      {
        int start = 1120, end = 1126;
        srep(i, start, end + 1)
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
          int start = 1120, end = 1126;
          srep(i, start, end + 1)
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


