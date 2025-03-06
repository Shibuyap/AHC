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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

int mode;

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
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
  // 0以上1未満の小数をとる乱数
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }
}  // namespace

const ll INF = 1001001001001001001;
struct edge
{
  int u, v, cost;
  int id;
};

// 定数
vector<edge> edges;
ll n, m, k;
ll x[110], y[110];
ll u[330], v[330], w[330];
ll a[5500], b[5500];
ll distn[110][5500];
ll distk[5500][110];
ll ndist[110][110];
ll ncost[110][110];
vector<P> nears[5500];
vector<P> nearsn[110];

// 更新する変数
ll maxScore;
ll real_maxScore;
ll p[110];
ll B[330];
ll real_p[110];
ll real_B[330];
int maxPowers[110] = {};
int real_maxPowers[110] = {};

// Union-Find
namespace
{
  void init(vector<int>& par, vector<int>& rank, int n)
  {
    for (int i = 0; i < n; i++) {
      par[i] = i;
      rank[i] = 0;
    }
  }

  // 木の根を求める
  int find(vector<int>& par, int x)
  {
    if (par[x] == x) {
      return x;
    }
    else {
      return par[x] = find(par, par[x]);
    }
  }

  // xとyの属する集合を併合
  void unite(vector<int>& par, vector<int>& rank, int x, int y)
  {
    x = find(par, x);
    y = find(par, y);
    if (x == y) return;

    if (rank[x] < rank[y]) {
      par[x] = y;
    }
    else {
      par[y] = x;
      if (rank[x] == rank[y]) rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool same(vector<int>& par, int x, int y)
  {
    return find(par, x) == find(par, y);
  }
}  // namespace

bool comp(const edge& e1, const edge& e2) { return e1.cost < e2.cost; }

int isNeed[110];
int real_isNeed[110];
int isUse[5500];
long long int kruscal(int V)
{
  int needCount = 0;
  rep(i, V) { needCount += isNeed[i]; }
  int E = edges.size();
  rep(i, E) isUse[i] = 0;

  vector<int> par(V + 10);   // 親
  vector<int> rank(V + 10);  // 木の深さ
  init(par, rank, V + 10);   // Union-Findの初期化
  long long int res = 0;
  int uniteCount = 0;
  for (int i = 0; i < E; i++) {
    edge e = edges[i];
    if (!isNeed[e.u] || !isNeed[e.v]) {
      continue;
    }
    if (!same(par, e.u, e.v)) {
      isUse[e.id] = 1;
      unite(par, rank, e.u, e.v);
      res += e.cost;
      uniteCount++;
    }
  }

  if (uniteCount != needCount - 1) res = INF;
  return res;
}

ll OuterKruscal() { return kruscal(n); }

ll OuterKruscal2()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  rep(i, n) { isNeed[i] = 1; }
  ll mi = OuterKruscal();
  int loop = 0;
  double TL = 0.1;
  while (true) {
    loop++;
    {
      endTime = clock();
      double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      double nowProgress = nowTime / TL;
      if (nowProgress > 1.0) break;
    }

    int num = randxor() % (n - 1) + 1;
    if (p[num] != 0) {
      continue;
    }

    isNeed[num] = 1 - isNeed[num];
    ll tmp = OuterKruscal();
    if (tmp <= mi) {
      // cout << mi << ' ' << tmp << endl;
      mi = tmp;
    }
    else {
      isNeed[num] = 1 - isNeed[num];
    }
  }

  // if (mode != 0) {
  //   cout << "kruscal loop = " << loop << endl;
  // }

  return mi;
}



bool Check()
{
  rep(i, k)
  {
    int ok = 0;
    int sz = nears[i].size();
    rep(j, sz)
    {
      int jj = nears[i][j].second;
      if (distk[i][jj] <= p[jj]) {
        ok = 1;
        break;
      }
    }
    if (!ok) {
      return false;
    }
  }
  return true;
}

bool Check2(int nn)
{
  for (auto ii : nearsn[nn]) {
    int i = ii.second;
    int ok = 0;
    int sz = nears[i].size();
    rep(j, sz)
    {
      int jj = nears[i][j].second;
      if (distk[i][jj] <= p[jj]) {
        ok = 1;
        break;
      }
    }
    if (!ok) {
      return false;
    }
  }
  return true;
}

ll keepW;
ll CalcScore(bool isALL = true)
{
  if (!Check()) {
    return -1;
  }

  double S = 0;
  rep(i, n) { S += p[i] * p[i]; }
  if (isALL) {
    keepW = 0;
    rep(i, m)
    {
      if (b[i]) {
        S += w[i];
        keepW += w[i];
      }
    }
  }
  else {
    S += keepW;
  }


  ll point = round(1e6 * (1.0 + 1e8 / (S + 1e7)));
  return point;
}

ll CalcScore2(int nn, bool isALL = true)
{
  if (!Check2(nn)) {
    return -1;
  }

  double S = 0;
  rep(i, n) { S += p[i] * p[i]; }
  if (isALL) {
    keepW = 0;
    rep(i, m)
    {
      if (b[i]) {
        S += w[i];
        keepW += w[i];
      }
    }
  }
  else {
    S += keepW;
  }


  ll point = round(1e6 * (1.0 + 1e8 / (S + 1e7)));
  return point;
}

ll Dist(ll x1, ll y1, ll x2, ll y2)
{
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

ll DistSquare(ll x1, ll y1, ll x2, ll y2)
{
  if (Dist(x1, y1, x2, y2) > 25000000) {
    return INF;
  }
  if (Dist(x1, y1, x2, y2) == 0) {
    return 0;
  }
  int ng = 0, ok = 5000;
  while (ng + 1 < ok) {
    int mid = (ok + ng) / 2;
    if (Dist(x1, y1, x2, y2) <= mid * mid) {
      ok = mid;
    }
    else {
      ng = mid;
    }
  }
  return ok;
}

void Init()
{
  rep(i, n) { p[i] = 0; }
  rep(i, m) { b[i] = 0; }
  rep(i, n) { isNeed[i] = 1; }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> n >> m >> k;
    rep(i, n) { cin >> x[i] >> y[i]; }
    rep(i, m) { cin >> u[i] >> v[i] >> w[i]; }
    rep(i, k) { cin >> a[i] >> b[i]; }

  }
  // ファイル入力する
  else {
    ifs >> n >> m >> k;
    rep(i, n) { ifs >> x[i] >> y[i]; }
    rep(i, m) { ifs >> u[i] >> v[i] >> w[i]; }
    rep(i, k) { ifs >> a[i] >> b[i]; }
  }

  rep(i, n)
  {
    rep(j, n)
    {
      ndist[i][j] = INF;
      ncost[i][j] = INF;
    }
  }

  rep(i, m)
  {
    u[i]--;
    v[i]--;
    ncost[u[i]][v[i]] = w[i];
    ncost[v[i]][u[i]] = w[i];
  }

  rep(i, k) { nears[i].clear(); }

  rep(i, n)
  {
    rep(j, k)
    {
      distn[i][j] = DistSquare(x[i], y[i], a[j], b[j]);
      distk[j][i] = distn[i][j];
      if (distk[j][i] <= 5000) {
        nears[j].push_back(P(distk[j][i], i));
      }
    }
  }

  rep(i, k) { sort(nears[i].begin(), nears[i].end()); }

  rep(i, m)
  {
    edge e;
    e.u = u[i];
    e.v = v[i];
    e.cost = w[i];
    e.id = i;
    edges.push_back(e);
  }

  sort(edges.begin(), edges.end(), comp);  // edge.costが小さい順にソートする

  rep(i, n)
  {
    rep(j, k)
    {
      if (distn[i][j] <= 5000) {
        nearsn[i].push_back(P(distn[i][j], j));
      }
    }
  }

  rep(i, n)
  {
    sort(nearsn[i].begin(), nearsn[i].end());
  }

  Init();
}

// 解答出力
void Output(int mode, int problemNum)
{
  if (mode == 0) {
    rep(i, n) { cout << p[i] << ' '; }
    cout << endl;
    rep(i, m) { cout << b[i] << ' '; }
    cout << endl;

  }
  else {
    // ファイル出力
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    rep(i, n) { ofs << p[i] << ' '; }
    ofs << endl;
    rep(i, m) { ofs << b[i] << ' '; }
    ofs << endl;

    ofs.close();
  }
}

void SetReal()
{
  real_maxScore = maxScore;
  rep(i, n) { real_p[i] = p[i]; }
  rep(i, m) { real_B[i] = B[i]; }
  rep(i, n) { real_maxPowers[i] = maxPowers[i]; }
  rep(i, n)
  {
    real_isNeed[i] = isNeed[i];
  }
}

// ランダムに1つ拡大縮小する
void Method1(double temperature)
{
  int num = randxor() % n;
  int pre = p[num];

  p[num] += randxor() % 51 - 25;
  p[num] = max(0LL, p[num]);
  p[num] = min(5000LL, p[num]);

  ll tmpScore = CalcScore();

  ll diffScore = tmpScore - maxScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > rand01()) {
    maxScore += diffScore;
    if (maxScore > real_maxScore) {
      SetReal();
    }
  }
  else {
    // 元に戻す
    p[num] = pre;
  }
}

void Solve1()
{
  // 解1
  rep(i, n) { p[i] = 5000; }
  rep(i, m) { b[i] = 1; }
}

void Solve2()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  OuterKruscal();
  rep(i, m) { b[i] = isUse[i]; }

  int flag[5500] = {};
  vector<P> farest;
  rep(i, k) { farest.push_back(nears[i][0]); }
  sort(farest.begin(), farest.end());

  drep(i, k)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, k)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  // 焼きなまし
  maxScore = CalcScore();
  real_maxScore = maxScore;
  rep(i, n) { real_p[i] = p[i]; }
  rep(i, m) { real_B[i] = B[i]; }

  double TL = 1.5;
  int loop = 0;
  double startTemperature = 200048;
  double endTemperature = 0;
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  double nowProgress = nowTime / TL;
  while (true) {
    loop++;
    if (loop % 100 == 0) {
      endTime = clock();
      nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
      if (nowProgress > 1.0) break;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    Method1(temperature);
  }

  if (mode != 0) {
    cout << loop << endl;
  }

  maxScore = real_maxScore;
  rep(i, n) { p[i] = real_p[i]; }
  rep(i, m) { B[i] = real_B[i]; }

  OuterKruscal2();
  rep(i, m) { b[i] = isUse[i]; }
}

void Solve3()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  OuterKruscal();
  rep(i, m) { b[i] = isUse[i]; }

  int flag[5500] = {};

  rep(i, k)
  {
    maxPowers[nears[i][0].second] =
      max(maxPowers[nears[i][0].second], nears[i][0].first);
  }
  vector<P> farest;
  rep(i, n) { farest.push_back(P(maxPowers[i], i)); }
  sort(farest.begin(), farest.end());

  drep(i, n)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, k)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  // 焼きなまし
  maxScore = CalcScore();
  real_maxScore = maxScore;
  rep(i, n) { real_p[i] = p[i]; }
  rep(i, m) { real_B[i] = B[i]; }

  rep(i, n) { real_maxPowers[i] = maxPowers[i]; }

  double TL_ALL = 1.5;
  int SET_COUNT = 2;
  rep(haibara, SET_COUNT)
  {
    startTime = clock();
    double TL = TL_ALL / SET_COUNT;

    int loop = 0;
    double startTemperature = 1000;
    double endTemperature = 0;
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
    double nowProgress = nowTime / TL;
    while (true) {
      loop++;
      if (loop % 10 == 0) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
        if (nowProgress > 1.0) break;
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      int pre_maxPowers[110];
      rep(i, n)
      {
        p[i] = 0;
        pre_maxPowers[i] = maxPowers[i];
      }

      int num = randxor() % n;


      maxPowers[num] += randxor() % 101 - 50;
      if (randxor() % 10 == 0) {
        maxPowers[num] += randxor() % 1001 - 500;
      }
      if (randxor() % 100 == 0) {
        maxPowers[num] += randxor() % 10001 - 5000;
      }
      maxPowers[num] = max(0, maxPowers[num]);
      maxPowers[num] = min(5000, maxPowers[num]);

      if (rand() % 10 == 0) {
        int num2 = randxor() % n;
        maxPowers[num2] += randxor() % 101 - 50;
        if (randxor() % 10 == 0) {
          maxPowers[num2] += randxor() % 1001 - 500;
        }
        if (randxor() % 100 == 0) {
          maxPowers[num2] += randxor() % 10001 - 5000;
        }
        maxPowers[num2] = max(0, maxPowers[num2]);
        maxPowers[num2] = min(5000, maxPowers[num2]);
      }


      rep(i, k) { flag[i] = 0; }

      farest.clear();
      rep(i, n) { farest.push_back(P(maxPowers[i], i)); }
      sort(farest.begin(), farest.end());

      drep(i, n)
      {
        int num = farest[i].second;
        int power = farest[i].first;
        if (p[num] != 0) {
          continue;
        }
        ll need = 0;
        for (const auto& pa : nearsn[num]) {
          int j = pa.second;
          ll kyori = pa.first;
          if (kyori > power) {
            break;
          }
          if (!flag[j]) {
            flag[j] = 1;
            need = max(need, kyori);
          }
        }
        //rep(j, k) {
        //  if (!flag[j] && distk[j][num] <= power) {
        //    flag[j] = 1;
        //    need = max(need, distk[j][num]);
        //  }
        //}
        p[num] = need;
      }

      ll tmpScore = CalcScore(false);

      ll diffScore = tmpScore - maxScore;

      double prob = exp((double)diffScore / temperature);
      if (prob > rand01()) {
        maxScore += diffScore;
        if (maxScore > real_maxScore) {
          if (mode != 0) {
            //cout << maxScore << ' ' << real_maxScore << endl;
          }
          SetReal();
        }
      }
      else {
        // 元に戻す
        rep(i, n) { maxPowers[i] = pre_maxPowers[i]; }
      }
    }

    if (mode != 0) {
      cout << loop << endl;
    }

    maxScore = real_maxScore;
    rep(i, n) { p[i] = real_p[i]; }
    rep(i, m) { B[i] = real_B[i]; }
    rep(i, n)
    {
      maxPowers[i] = real_maxPowers[i];
    }
  }

  OuterKruscal2();
  rep(i, m) { b[i] = isUse[i]; }
}


ll OuterKruscal3()
{

  ll tmp = OuterKruscal();

  ll res = 0;
  if (tmp != INF) {

    res = CalcScore();
    cout << res << endl;
  }


  return res;
}

void Solve4()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  OuterKruscal();
  rep(i, m) { b[i] = isUse[i]; }

  int flag[5500] = {};

  rep(i, k)
  {
    maxPowers[nears[i][0].second] =
      max(maxPowers[nears[i][0].second], nears[i][0].first);
  }
  vector<P> farest;
  rep(i, n) { farest.push_back(P(maxPowers[i], i)); }
  sort(farest.begin(), farest.end());

  drep(i, n)
  {
    int num = farest[i].second;
    int power = farest[i].first;
    if (p[num] != 0) {
      continue;
    }
    ll need = 0;
    rep(j, k)
    {
      if (!flag[j] && distk[j][num] <= power) {
        flag[j] = 1;
        need = max(need, distk[j][num]);
      }
    }
    p[num] = need;
  }

  OuterKruscal2();
  rep(i, m) { b[i] = isUse[i]; }

  // 焼きなまし
  maxScore = CalcScore();
  real_maxScore = maxScore;
  rep(i, n) { real_p[i] = p[i]; }
  rep(i, m) { real_B[i] = B[i]; }
  rep(i, n) { real_maxPowers[i] = maxPowers[i]; }
  rep(i, n)
  {
    real_isNeed[i] = isNeed[i];
  }

  double TL_ALL = 1.5;
  int SET_COUNT = 2;
  rep(haibara, SET_COUNT)
  {
    startTime = clock();
    double TL = TL_ALL / SET_COUNT;

    int loop = 0;
    double startTemperature = 1000;
    double endTemperature = 0;
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
    double nowProgress = nowTime / TL;
    while (true) {
      loop++;
      if (loop % 10 == 0) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
        if (nowProgress > 1.0) break;
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      if (randxor() % 10 != 0) {
        int pre_maxPowers[110];
        int pre_p[110];
        rep(i, n)
        {
          pre_p[i] = p[i];
          p[i] = 0;
          pre_maxPowers[i] = maxPowers[i];
        }

        int num = randxor() % n;


        maxPowers[num] += randxor() % 101 - 50;
        maxPowers[num] = max(0, maxPowers[num]);
        maxPowers[num] = min(5000, maxPowers[num]);

        rep(i, k) { flag[i] = 0; }

        farest.clear();
        rep(i, n) { farest.push_back(P(maxPowers[i], i)); }
        sort(farest.begin(), farest.end());

        drep(i, n)
        {
          int num = farest[i].second;
          int power = farest[i].first;
          if (p[num] != 0) {
            continue;
          }
          ll need = 0;
          for (const auto& pa : nearsn[num]) {
            int j = pa.second;
            ll kyori = pa.first;
            if (kyori > power) {
              break;
            }
            if (!flag[j]) {
              flag[j] = 1;
              need = max(need, kyori);
            }
          }
          //rep(j, k) {
          //  if (!flag[j] && distk[j][num] <= power) {
          //    flag[j] = 1;
          //    need = max(need, distk[j][num]);
          //  }
          //}
          p[num] = need;
        }

        int isUpdate = 0;
        vector<int> isNeedKeep;
        rep(i, n)
        {
          if (pre_p[i] == 0 && p[i] > 0 && isNeed[i] == 0) {
            isNeed[i] = 1;
            isNeedKeep.push_back(i);
            isUpdate = 1;
            break;
          }
        }
        ll tmpScore = 0;
        if (isUpdate) {
          tmpScore= OuterKruscal3();
        }
        else {
          tmpScore = CalcScore(false);
        }


        ll diffScore = tmpScore - maxScore;

        double prob = exp((double)diffScore / temperature);
        if (prob > rand01()) {
          maxScore += diffScore;
          if (maxScore > real_maxScore) {
            if (mode != 0) {
              //cout << maxScore << ' ' << real_maxScore << endl;
            }
            SetReal();
          }
        }
        else {
          // 元に戻す
          rep(i, n) { maxPowers[i] = pre_maxPowers[i]; }
          for (auto ite : isNeedKeep) {
            isNeed[ite] = 0;
          }
        }
      }
      else {
        int num = randxor() % (n - 1) + 1;
        if (p[num] != 0) {
          continue;
        }
        isNeed[num] = 1 - isNeed[num];
        ll tmpScore = OuterKruscal3();

        ll diffScore = tmpScore - maxScore;

        double prob = exp((double)diffScore / temperature);
        if (prob > rand01()) {
          maxScore += diffScore;
          if (maxScore > real_maxScore) {
            if (mode != 0) {
              cout << maxScore << ' ' << real_maxScore << endl;
            }
            SetReal();
          }
        }
        else {
          // 元に戻す
          isNeed[num] = 1 - isNeed[num];
        }
      }
    }

    if (mode != 0) {
      cout << loop << endl;
    }

    maxScore = real_maxScore;
    rep(i, n) { p[i] = real_p[i]; }
    rep(i, m) { B[i] = real_B[i]; }
    rep(i, n)
    {
      maxPowers[i] = real_maxPowers[i];
    }
    rep(i, n)
    {
      isNeed[i] = real_isNeed[i];
    }
  }

  OuterKruscal2();
  rep(i, m) { b[i] = isUse[i]; }
}


double Solve(int mode, int problemNum = 0)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // Solve1();
  // Solve2();
  Solve3();
  //Solve4();

  Check();

  if (mode != 0) {
    cout << CalcScore() << endl;
  }
  return CalcScore();
}

double SolveOuter(int mode, int problemNum = 0)
{
  // 入力受け取り
  Input(problemNum);

  double score = Solve(mode, problemNum);

  // 解答の出力
  Output(mode, problemNum);

  return score;
}

int main()
{
  // 乱数調整
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 0;

  // 提出用
  if (mode == 0) {
    SolveOuter(mode);
  }
  // 1ケース試す
  else if (mode == 1) {
    SolveOuter(mode, 9);
  }
  // 複数ケース試す
  else if (mode == 2) {
    rep(i, 10) { SolveOuter(mode, i); }
  }

  return 0;
}
