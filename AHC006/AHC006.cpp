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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

/*
いろいろ
const int INF = 1001001001;
const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, -1, 0, 1};
const char cc[4] = {'U','L','D','R'};
*/

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t        = x ^ (x << 11);
    x        = y;
    y        = z;
    z        = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0以上1未満の小数をとる乱数
  static double rand01()
  {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* 変数 */
{
  // 入力用変数
  int n = 1000;
  int m = 60;
  vector<int> a(1000), b(1000), c(1000), d(1000);
  vector<int> argA;

  // 解答用変数
  ll maxScore;
  ll minTime;
  vector<P> ans;
  vector<int> argAns(140);
  vector<int> argAns2(140);
  vector<int> use(1100);
  vector<int> pair_(2100);

  // 焼きなまし用変数
  ll real_maxScore;
  ll real_minTime;
  vector<P> real_ans;
  vector<int> real_argAns(140);
  vector<int> real_argAns2(140);
  vector<int> real_use(1100);
  vector<int> real_pair_(2100);

  void ResetParam()
  {
    argA.clear();
    maxScore = 0;
    minTime  = 0;
    ans.clear();
    argAns.resize(140);
    argAns2.resize(140);
    use.resize(1100);
    pair_.resize(2100);
    real_maxScore = 0;
    real_minTime  = 0;
    real_ans.clear();
  }

}  // namespace

// スコア計算
ll CalcScore(ll time)
{
  ll res = round(100000000.0 / (1000.0 + time));
  return res;
}

ll CalcTime()
{
  ll timeSum = 0;
  srep(i, 1, ans.size())
  {
    timeSum += abs(ans[i].first - ans[i - 1].first) + abs(ans[i].second - ans[i - 1].second);
  }
  return timeSum;
}

bool IsGood(int ite, int L, int R, int U, int D)
{
  if (a[ite] < L || R < a[ite]) return false;
  if (b[ite] < U || D < b[ite]) return false;
  if (c[ite] < L || R < c[ite]) return false;
  if (d[ite] < U || D < d[ite]) return false;
  return true;
}
bool IsGood2(int ite, int L1, int R1, int U1, int D1, int L2, int R2, int U2, int D2)
{
  if (a[ite] < L1 || R1 < a[ite]) return false;
  if (b[ite] < U1 || D1 < b[ite]) return false;
  if (c[ite] < L2 || R2 < c[ite]) return false;
  if (d[ite] < U2 || D2 < d[ite]) return false;
  return true;
}

void FilterInput(clock_t& start_time)
{
  vector<int> aa, bb, cc, dd;
  int mi   = 1001001;
  int loop = 0;
  clock_t end_time;
  double now_time;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      end_time = clock();
      now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    }
    if (now_time > 0.2) break;

    vector<int> aaa, bbb, ccc, ddd, argAA;
    int L1 = randxor() % 801;
    int R1 = randxor() % 801;
    int U1 = randxor() % 801;
    int D1 = randxor() % 801;
    int L2 = randxor() % 801;
    int R2 = randxor() % 801;
    int U2 = randxor() % 801;
    int D2 = randxor() % 801;
    if (L1 > R1) swap(L1, R1);
    if (U1 > D1) swap(U1, D1);
    if (L2 > R2) swap(L2, R2);
    if (U2 > D2) swap(U2, D2);

    // int sz = (R1 - L1) * (D1 - U1) + (R2 - L2) * (D2 - U2);
    int sz = (R1 - L1) + (D1 - U1);
    if (sz >= mi) continue;

    rep(i, n)
    {
      if (IsGood(i, L1, R1, U1, D1)) {
        // if (IsGood2(i, L1, R1, U1, D1, L2, R2, U2, D2)) {
        aaa.push_back(a[i]);
        bbb.push_back(b[i]);
        ccc.push_back(c[i]);
        ddd.push_back(d[i]);
        argAA.push_back(i);
      }
    }
    if (aaa.size() >= m + 1) {
      aa   = aaa;
      bb   = bbb;
      cc   = ccc;
      dd   = ddd;
      argA = argAA;
      mi   = sz;
    }
  }

  if (aa.size() >= m) {
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    n = a.size();
  }
}

///////////////////////////////////////////////////////////////////////////////
// 始まったらやること
// 1. 入力部
// 2. maxScoreとans
// 3. 出力部
// 4. 愚直解
// 5. スコア計算関数
// 6. 貪欲解
// 7. 山登り
// 8. 焼きなまし
///////////////////////////////////////////////////////////////////////////////
int Solve(int mode)
{
  srand((unsigned)time(NULL));
  clock_t start_time, end_time;
  start_time = clock();
  end_time   = clock();
  while (rand() % 100) {
    randxor();
  }

  // 入力部
  string fileNameIfs = "./in/0000.txt";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    rep(i, n)
    {
      cin >> a[i] >> b[i] >> c[i] >> d[i];
    }
  }
  else {  // ファイル入力する
    rep(i, n)
    {
      ifs >> a[i] >> b[i] >> c[i] >> d[i];
    }
  }
  ifs.close();

  FilterInput(start_time);

  // 愚直解
  ans.push_back(P(400, 400));
  rep(i, m)
  {
    ans.push_back(P(a[i], b[i]));
    argAns[i] = i;
    use[i]    = 1;
  }
  rep(i, m)
  {
    ans.push_back(P(c[i], d[i]));
  }
  ans.push_back(P(400, 400));
  minTime  = CalcTime();
  maxScore = CalcScore(minTime);

  real_ans      = ans;
  real_minTime  = minTime;
  real_maxScore = maxScore;

  // 山登り解、焼きなまし解
  double now_time   = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  double TL         = 1.8;
  double start_temp = 48;
  double end_temp   = 0.0001;
  int loop          = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      end_time = clock();
      now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    }
    if (now_time > 0.4) break;

    int x = randxor() % n;

    while (use[x]) {
      x = randxor() % n;
    }

    int ite = randxor() % m;

    int keepNum = argAns[ite];
    int diff    = 0;
    diff += abs(a[x] - ans[ite].first) + abs(b[x] - ans[ite].second);
    diff += abs(a[x] - ans[ite + 2].first) + abs(b[x] - ans[ite + 2].second);
    diff += abs(c[x] - ans[ite + m].first) + abs(d[x] - ans[ite + m].second);
    diff += abs(c[x] - ans[ite + m + 2].first) + abs(d[x] - ans[ite + m + 2].second);

    diff -= abs(ans[ite + 1].first - ans[ite].first) + abs(ans[ite + 1].second - ans[ite].second);
    diff -= abs(ans[ite + 1].first - ans[ite + 2].first) + abs(ans[ite + 1].second - ans[ite + 2].second);
    diff -= abs(ans[ite + m + 1].first - ans[ite + m].first) + abs(ans[ite + m + 1].second - ans[ite + m].second);
    diff -= abs(ans[ite + m + 1].first - ans[ite + m + 2].first) + abs(ans[ite + m + 1].second - ans[ite + m + 2].second);

    int tmpTime   = minTime + diff;
    int diffScore = -diff;

    double temp = start_temp + (end_temp - start_temp) * now_time / TL;
    double prob = exp((double)diffScore / temp);
    if (prob > rand01()) {
      minTime          = tmpTime;
      use[argAns[ite]] = 0;
      use[x]           = 1;
      argAns[ite]      = x;
      ans[ite + 1]     = P(a[x], b[x]);
      ans[ite + m + 1] = P(c[x], d[x]);
      maxScore         = CalcScore(minTime);
      if (maxScore > real_maxScore) {
        real_minTime  = minTime;
        real_maxScore = maxScore;
        real_ans      = ans;
        real_argAns   = argAns;
        real_use      = use;
      }
    }
    else {
      // 元に戻す
      ;
    }
  }

  // 最高スコアを戻す
  ans      = real_ans;
  minTime  = real_minTime;
  maxScore = real_maxScore;
  argAns   = real_argAns;
  use      = real_use;

  rep(i, m)
  {
    pair_[argAns[i]]        = i + m;
    pair_[argAns[i] + 1000] = i;
  }
  rep(i, m)
  {
    argAns[i + m] = argAns[i] + 1000;
  }
  real_argAns = argAns;
  real_pair_  = pair_;

  argAns2      = argAns;
  real_argAns2 = argAns2;

  // それぞれをTSP
  rep(ui_tei, 10)
  {
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        end_time = clock();
        now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
      }
      if (now_time > (TL / 15) * (6 + ui_tei)) break;

      int ite1 = randxor() % (m * 2);
      int ite2 = randxor() % (m * 2);
      while (ite1 == ite2) {
        ite2 = randxor() % (m * 2);
      }
      if (ite1 > ite2) swap(ite1, ite2);
      if (ite2 - ite1 == 1) continue;
      // if (ite1 == 0 && ite2 == m - 1) continue;

      if (argAns[ite1] + 1000 == argAns[ite2]) continue;
      if (argAns[ite1] < 1000) {
        int pa = pair_[argAns[ite1]];
        if (pa < ite2) { continue; }
      }
      if (argAns[ite2] >= 1000) {
        int pa = pair_[argAns[ite2]];
        if (ite1 < pa) { continue; }
      }

      int diff = 0;

      if (true || loop % 2 == 0) {
        diff += abs(ans[ite1 + 1].first - ans[ite2].first) + abs(ans[ite1 + 1].second - ans[ite2].second);
        diff += abs(ans[ite1 + 1].first - ans[ite2 + 2].first) + abs(ans[ite1 + 1].second - ans[ite2 + 2].second);
        diff -= abs(ans[ite1 + 1].first - ans[ite1].first) + abs(ans[ite1 + 1].second - ans[ite1].second);
        diff -= abs(ans[ite1 + 1].first - ans[ite1 + 2].first) + abs(ans[ite1 + 1].second - ans[ite1 + 2].second);

        diff += abs(ans[ite2 + 1].first - ans[ite1].first) + abs(ans[ite2 + 1].second - ans[ite1].second);
        diff += abs(ans[ite2 + 1].first - ans[ite1 + 2].first) + abs(ans[ite2 + 1].second - ans[ite1 + 2].second);
        diff -= abs(ans[ite2 + 1].first - ans[ite2].first) + abs(ans[ite2 + 1].second - ans[ite2].second);
        diff -= abs(ans[ite2 + 1].first - ans[ite2 + 2].first) + abs(ans[ite2 + 1].second - ans[ite2 + 2].second);

        int tmpTime   = minTime + diff;
        int diffScore = -diff;

        double temp = start_temp + (end_temp - start_temp) * now_time / TL;
        double prob = exp((double)diffScore / temp);
        if (prob > rand01()) {
          minTime = tmpTime;
          if (argAns[ite1] < 1000) {
            pair_[argAns[ite1] + 1000] = ite2;
          }
          else {
            pair_[argAns[ite1] - 1000] = ite2;
          }
          if (argAns[ite2] < 1000) {
            pair_[argAns[ite2] + 1000] = ite1;
          }
          else {
            pair_[argAns[ite2] - 1000] = ite1;
          }
          swap(argAns[ite1], argAns[ite2]);
          swap(ans[ite1 + 1], ans[ite2 + 1]);
          maxScore = CalcScore(minTime);
          if (maxScore > real_maxScore) {
            real_minTime  = minTime;
            real_maxScore = maxScore;
            real_ans      = ans;
            real_argAns   = argAns;
            real_use      = use;
            real_pair_    = pair_;
          }
        }
        else {
          // 元に戻す
          ;
        }
      }
      // else {
      //      diff += abs(ans[ite1 + m + 1].first - ans[ite2 + m].first) + abs(ans[ite1 + m + 1].second - ans[ite2 + m].second);
      //      diff += abs(ans[ite1 + m + 1].first - ans[ite2 + m + 2].first) + abs(ans[ite1 + m + 1].second - ans[ite2 + m + 2].second);
      //      diff -= abs(ans[ite1 + m + 1].first - ans[ite1 + m].first) + abs(ans[ite1 + m + 1].second - ans[ite1 + m].second);
      //      diff -= abs(ans[ite1 + m + 1].first - ans[ite1 + m + 2].first) + abs(ans[ite1 + m + 1].second - ans[ite1 + m + 2].second);

      //     diff += abs(ans[ite2 + m + 1].first - ans[ite1 + m].first) + abs(ans[ite2 + m + 1].second - ans[ite1 + m].second);
      //     diff += abs(ans[ite2 + m + 1].first - ans[ite1 + m + 2].first) + abs(ans[ite2 + m + 1].second - ans[ite1 + m + 2].second);
      //     diff -= abs(ans[ite2 + m + 1].first - ans[ite2 + m].first) + abs(ans[ite2 + m + 1].second - ans[ite2 + m].second);
      //     diff -= abs(ans[ite2 + m + 1].first - ans[ite2 + m + 2].first) + abs(ans[ite2 + m + 1].second - ans[ite2 + m + 2].second);

      //     int tmpTime   = minTime + diff;
      //     int diffScore = -diff;

      //     double temp = start_temp + (end_temp - start_temp) * now_time / TL;
      //     double prob = exp((double)diffScore / temp);
      //     if (prob > rand01()) {
      //         minTime = tmpTime;
      //         swap(argAns2[ite1], argAns2[ite2]);
      //         swap(ans[ite1 + m + 1], ans[ite2 + m + 1]);
      //         maxScore = CalcScore(minTime);
      //         if (maxScore > real_maxScore) {
      //             real_minTime  = minTime;
      //             real_maxScore = maxScore;
      //             real_ans      = ans;
      //             real_argAns2  = argAns2;
      //             real_use      = use;
      //         }
      //     } else {
      //         // 元に戻す
      //         ;
      //     }
      // }
    }

    // 最高スコアを戻す
    ans      = real_ans;
    minTime  = real_minTime;
    maxScore = real_maxScore;
    argAns   = real_argAns;
    argAns2  = real_argAns2;
    use      = real_use;
    pair_    = real_pair_;
  }

  TL = 1.95;

  const int INF  = 1001001001;
  int mi         = INF;
  argAns2        = argAns;
  vector<P> ans2 = ans;

  random_device seed_gen;
  mt19937 engine(seed_gen());

  vector<int> v;
  rep(i, m)
  {
    v.push_back(i);
  }

  while (true) {
    loop++;
    if (loop % 100 == 1) {
      end_time = clock();
      now_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    }
    if (now_time > TL) break;

    std::shuffle(v.begin(), v.end(), engine);
    vector<int> vv;
    rep(i, 50) vv.push_back(v[i]);
    sort(vv.begin(), vv.end());

    set<int> useIte;
    int cnt = 0;
    int now = 0;
    rep(i, m * 2)
    {
      if (argAns[i] < 1000) {
        if (cnt == vv[now]) {
          useIte.insert(argAns[i]);
          now++;
        }
        cnt++;
      }
    }

    int tmpTime = 0;
    int nowx    = 400;
    int nowy    = 400;
    vector<P> ans3;
    ans3.push_back(P(400, 400));
    rep(i, m * 2)
    {
      int ite = argAns[i];
      if (ite >= 1000) ite -= 1000;
      if (useIte.find(ite) != useIte.end()) {
        ans3.push_back(ans[i + 1]);
        tmpTime += abs(nowx - ans[i + 1].first) + abs(nowy - ans[i + 1].second);
        nowx = ans[i + 1].first;
        nowy = ans[i + 1].second;
      }
    }
    ans3.push_back(P(400, 400));

    if (tmpTime < mi) {
      mi   = tmpTime;
      ans2 = ans3;
      argAns2.clear();

      for (auto ite : useIte) {
        argAns2.push_back(ite);
      }
    }
  }

  ans    = ans2;
  argAns = argAns2;

  minTime  = CalcTime();
  maxScore = CalcScore(minTime);

  // 解の出力
  if (mode == 0) {
    cout << 50;
    rep(i, 50)
    {
      if (argAns[i] < 1000) { cout << " " << argA[argAns[i]] + 1; }
    }
    cout << endl;
    cout << ans.size();
    rep(i, ans.size())
    {
      cout << " " << ans[i].first << " " << ans[i].second;
    }
    cout << endl;
  }

  // デバッグ用
  if (mode != 0) {
    cout << "loop = " << loop << endl;
    cout << maxScore << endl;
    end_time = clock();
    cout << (double)(end_time - start_time) / CLOCKS_PER_SEC << "sec." << endl;
  }

  // ファイル出力
  if (true) {
    string fileNameOfs = "0000_out.txt";
    ofstream ofs(fileNameOfs);

    ofs << 50;
    rep(i, 50)
    {
      if (argAns[i] < 1000) { ofs << " " << argA[argAns[i]] + 1; }
    }
    ofs << endl;
    ofs << ans.size();
    rep(i, ans.size())
    {
      ofs << " " << ans[i].first << " " << ans[i].second;
    }
    ofs << endl;

    ofs.close();
  }

  // cout << maxScore << endl;

  return maxScore;
}

int main()
{
  int mode = 0;
  m        = 55;
  if (mode == 0) {
    Solve(mode);
  }
  else if (mode == 1) {
    Solve(mode);
  }

  return 0;
}
