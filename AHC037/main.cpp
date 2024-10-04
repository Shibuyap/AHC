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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

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

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j   = randxor() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF      = 1001001001001001001;
const int INT_INF = 1001001001;
const int MA      = 1000000000;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
clock_t startTime, endTime;

double GetNowTime()
{
  endTime        = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int n     = 1000;
const int MAX_N = 5000;
ll L;

ll a[6000], b[6000];
ll ans[5500][4];
int ansN;
int ansCount;
ll ansCost;

ll real_a[6000], real_b[6000];
ll real_ans[5500][4];
int real_ansN;
int real_ansCount;
ll real_ansCost;

void CopyToReal()
{
  real_ansN     = ansN;
  real_ansCount = ansCount;
  real_ansCost  = ansCost;

  rep(i, ansN)
  {
    real_a[i] = a[i];
    real_b[i] = b[i];
  }

  rep(i, ansCount)
  {
    rep(j, 4)
    {
      real_ans[i][j] = ans[i][j];
    }
  }
}

void CopyToAns()
{
  ansN     = real_ansN;
  ansCount = real_ansCount;
  ansCost  = real_ansCost;

  rep(i, ansN)
  {
    a[i] = real_a[i];
    b[i] = real_b[i];
  }

  rep(i, ansCount)
  {
    rep(j, 4)
    {
      ans[i][j] = real_ans[i][j];
    }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ansCount = 0;
  ansCost  = 0;
}

// 入力受け取り
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
    int nnn;
    cin >> nnn;
    rep(i, n)
    {
      cin >> a[i] >> b[i];
    }
  }
  // ファイル入力する
  else {
    int nnn;
    ifs >> nnn;
    rep(i, n)
    {
      ifs >> a[i] >> b[i];
    }
  }

  L = 0;
  rep(i, n)
  {
    L = max(L, a[i]);
    L = max(L, b[i]);
  }

  vector<P> vp;
  rep(i, n)
  {
    vp.push_back(P(a[i], b[i]));
  }
  sort(vp.begin(), vp.end());
  rep(i, n)
  {
    a[i] = vp[i].first;
    b[i] = vp[i].second;
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
ll CalcScore()
{
  ll res = round(1.0 * n * L / (ansCost + 1) * 1000000.0);
  return res;
}

void Ans(int x1, int y1, int x2, int y2)
{
  ans[ansCount][0] = x1;
  ans[ansCount][1] = y1;
  ans[ansCount][2] = x2;
  ans[ansCount][3] = y2;
  ansCount++;
  ansCost += (x2 - x1) + (y2 - y1);
}

void EliminateDup()
{
  set<P> se;
  int f[5500] = {};
  rep(i, ansN)
  {
    P p = P(a[i], b[i]);
    if (se.find(p) == se.end()) {
      se.insert(p);
    }
    else {
      f[i] = 1;
    }
  }

  int nn = ansN;
  ansN   = 0;
  rep(i, nn)
  {
    if (f[i]) continue;
    a[ansN] = a[i];
    b[ansN] = b[i];
    ansN++;
  }
}

int pa[6000] = {};
void Initialize3()
{
  ansCost  = 0;
  ansCount = 0;

  rep(i, ansN)
  {
    ll mi  = a[i] + b[i];
    ll idx = -1;
    rep(j, ansN)
    {
      if (i == j) continue;
      if (a[j] <= a[i] && b[j] <= b[i] && (a[i] + b[i]) - (a[j] + b[j]) < mi) {
        mi  = (a[i] + b[i]) - (a[j] + b[j]);
        idx = j;
      }
    }
    pa[i] = idx;
  }

  queue<int> que;
  rep(i, ansN)
  {
    if (pa[i] == -1) {
      Ans(0, 0, a[i], b[i]);
      que.push(i);
    }
  }
  while (que.size()) {
    int i = que.front();
    que.pop();
    rep(j, ansN)
    {
      if (pa[j] == i) {
        Ans(a[i], b[i], a[j], b[j]);
        que.push(j);
      }
    }
  }
}

int sonCount[6000];
int del[6000];
void Initialize32(int dell = 1)
{
  ansCost  = 0;
  ansCount = 0;

  rep(i, ansN)
  {
    sonCount[i] = 0;
    del[i]      = 0;
  }

  rep(i, ansN)
  {
    ll mi  = a[i] + b[i];
    ll idx = -1;
    rep(j, ansN)
    {
      if (i == j) continue;
      if (a[j] <= a[i] && b[j] <= b[i] && (a[i] + b[i]) - (a[j] + b[j]) < mi) {
        mi  = (a[i] + b[i]) - (a[j] + b[j]);
        idx = j;
      }
    }
    pa[i] = idx;

    if (idx != -1) {
      sonCount[idx]++;
    }
  }

  queue<int> que;
  srep(i, n, ansN)
  {
    if (sonCount[i] == 0) {
      del[i] = 1;
      que.push(i);
    }
  }
  while (que.size()) {
    int i = que.front();
    que.pop();
    if (pa[i] != -1) {
      sonCount[pa[i]]--;
      if (sonCount[pa[i]] == 0 && pa[i] >= n) {
        del[pa[i]] = 1;
        que.push(pa[i]);
      }
    }
  }

  if (dell == 1) {
    while (true) {
      int flag = 0;
      srep(i, n, ansN)
      {
        if (del[i]) continue;
        if (sonCount[i] == 1) {
          int pai = pa[i];
          rep(j, ansN)
          {
            if (del[i]) continue;
            if (pa[j] == i) {
              pa[j]  = pai;
              del[i] = 1;
              flag   = 1;
              break;
            }
          }
        }
      }
      if (flag == 0) break;
    }
  }

  rep(i, ansN)
  {
    if (del[i]) continue;
    if (pa[i] == -1) {
      Ans(0, 0, a[i], b[i]);
      que.push(i);
    }
  }
  while (que.size()) {
    int i = que.front();
    que.pop();
    if (del[i]) continue;
    rep(j, ansN)
    {
      if (del[j]) continue;
      if (pa[j] == i) {
        Ans(a[i], b[i], a[j], b[j]);
        que.push(j);
      }
    }
  }

  int newAnsN = 0;
  rep(i, ansN)
  {
    if (del[i]) continue;
    a[newAnsN] = a[i];
    b[newAnsN] = b[i];
    newAnsN++;
  }
  ansN = newAnsN;
}

void Initialize6()
{
  Initialize32();

  ll beforeAnsCost = ansCost;

  while (true) {
    if (GetNowTime() > TL) break;
    int newAnsN = ansN;
    rep(i, ansN)
    {
      if (newAnsN == MAX_N) break;
      if (del[i]) continue;
      rep(j, ansN)
      {
        if (newAnsN == MAX_N) break;
        if (del[j]) continue;
        if (pa[i] == pa[j]) {
          a[newAnsN] = min(a[i], a[j]);
          b[newAnsN] = min(b[i], b[j]);
          newAnsN++;
        }
      }
    }

    ansN = newAnsN;
    EliminateDup();

    Initialize32();

    if (ansCost == beforeAnsCost) {
      break;
    }
    beforeAnsCost = ansCost;
  }

  Initialize32();
}

void Initialize8()
{
  Initialize32(0);

  ll beforeAnsCost = ansCost;
  int lastOne      = 0;

  while (true) {
    if (GetNowTime() > TL) break;

    set<P> se;
    rep(i, ansN)
    {
      se.insert(P(a[i], b[i]));
    }

    int newAnsN = ansN;

    rep(i, ansCount)
    {
      if (newAnsN == MAX_N) break;
      int aa = ans[i][0];
      int bb = ans[i][3];
      if (se.find(P(aa, bb)) == se.end()) {
        se.insert(P(aa, bb));
        a[newAnsN] = aa;
        b[newAnsN] = bb;
        newAnsN++;
      }

      if (newAnsN == MAX_N) break;
      bb = ans[i][1];
      aa = ans[i][2];
      if (se.find(P(aa, bb)) == se.end()) {
        se.insert(P(aa, bb));
        a[newAnsN] = aa;
        b[newAnsN] = bb;
        newAnsN++;
      }
    }

    ansN = newAnsN;
    EliminateDup();

    Initialize32(0);

    if (ansCost == beforeAnsCost) {
      lastOne++;
      if (lastOne >= 2) {
        break;
      }
    }
    beforeAnsCost = ansCost;
  }

  Initialize32(0);
}

void Initialize4()
{
  ansN = n;
  vector<P> vpA, vpB;
  rep(i, n)
  {
    vpA.push_back(P(a[i], b[i]));
    vpB.push_back(P(b[i], a[i]));
  }

  sort(vpA.begin(), vpA.end());
  sort(vpB.begin(), vpB.end());

  rep(i, n)
  {
    int ai  = vpA[i].first;
    int bi  = vpA[i].second;
    int cnt = 0;
    srep(j, i + 1, n)
    {
      if (cnt == 2) break;
      int bj = vpA[j].second;
      if (bj < bi) {
        a[ansN] = ai;
        b[ansN] = bj;
        ansN++;
        cnt++;
      }
    }
  }

  rep(i, n)
  {
    int ai  = vpB[i].second;
    int bi  = vpB[i].first;
    int cnt = 0;
    srep(j, i + 1, n)
    {
      if (cnt == 2) break;
      int aj = vpB[j].second;
      if (aj < ai) {
        a[ansN] = aj;
        b[ansN] = bi;
        ansN++;
        cnt++;
      }
    }
  }

  set<P> se;
  int f[5500] = {};
  rep(i, ansN)
  {
    P p = P(a[i], b[i]);
    if (se.find(p) == se.end()) {
      se.insert(p);
    }
  }

  ansN = 0;
  for (auto p : se) {
    a[ansN] = p.first;
    b[ansN] = p.second;
    ansN++;
  }

  if (mode != 0) {
    cout << "ansN = " << ansN << endl;
  }

  Initialize3();
}

void Initialize5()
{
  Initialize3();

  rep(i, ansCount)
  {
    a[ansN] = ans[i][0];
    b[ansN] = ans[i][3];
    ansN++;
    a[ansN] = ans[i][1];
    b[ansN] = ans[i][2];
    ansN++;
  }

  EliminateDup();

  ansCost  = 0;
  ansCount = 0;
  Initialize32();
}

ll Dist(ll x1, ll y1, ll x2, ll y2)
{
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

void Initialize7()
{
  priority_queue<pair<ll, P>, vector<pair<ll, P>>, greater<pair<ll, P>>> que;

  pair<ll, P> pp;
  rep(i, ansN)
  {
    srep(j, i + 1, ansN)
    {
      pp.first         = Dist(a[i], b[i], a[j], b[j]);
      pp.second.first  = i;
      pp.second.second = j;
      if ((double)pp.first > 1e16) continue;
      que.push(pp);
    }
  }

  set<P> se;
  rep(i, n)
  {
    se.insert(P(a[i], b[i]));
  }

  ansN = n;
  while (que.size() && ansN < MAX_N) {
    pp = que.top();
    que.pop();
    int i  = pp.second.first;
    int j  = pp.second.second;
    int aa = min(a[i], a[j]);
    int bb = min(b[i], b[j]);
    if (se.find(P(aa, bb)) != se.end()) continue;

    a[ansN] = aa;
    b[ansN] = bb;
    ansN++;
    se.insert(P(aa, bb));

    // rep(ii, n) {
    //   pp.first = Dist(a[ii], b[ii], a[ansN - 1], b[ansN - 1]);
    //   if ((double)pp.first > 1e16) continue;
    //   pp.second.first  = ii;
    //   pp.second.second = ansN - 1;
    //   que.push(pp);
    // }
  }

  Initialize3();
  Initialize32();
}

void Method1()
{
  while (true) {
    if (GetNowTime() > TL) {
      break;
    }

    int flag = 0;

    rep(i, ansN)
    {
      if (ansN == MAX_N) break;

      int ansIdx = -1;
      ll cost    = 0;
      rep(j, ansCount)
      {
        if (ans[j][2] == a[i] && ans[j][3] == b[i]) {
          ansIdx = j;
          cost   = abs(ans[j][2] - ans[j][0]) + abs(ans[j][3] - ans[j][1]);
          break;
        }
      }

      rep(j, ansCount)
      {
        if (j == ansIdx) continue;
        if (ans[j][0] < a[i]) {
          if (ans[j][1] < b[i] && b[i] < ans[j][3]) {
            if (a[i] - ans[j][0] < cost) {
              a[ansN] = ans[j][0];
              b[ansN] = b[i];

              ans[ansIdx][0] = a[ansN];
              ans[ansIdx][1] = b[ansN];
              ans[ansIdx][2] = a[i];
              ans[ansIdx][3] = b[i];

              int keepA = ans[j][2];
              int keepB = ans[j][3];
              ans[j][2] = a[ansN];
              ans[j][3] = b[ansN];

              ans[ansCount][0] = a[ansN];
              ans[ansCount][1] = b[ansN];
              ans[ansCount][2] = keepA;
              ans[ansCount][3] = keepB;
              ansCount++;

              ansCost -= cost;
              ansCost += (ans[ansIdx][2] + ans[ansIdx][3]) - (ans[ansIdx][0] + ans[ansIdx][1]);

              ansN++;
              flag = 1;

              break;
            }
          }
        }

        if (ans[j][1] < b[i]) {
          if (ans[j][0] < a[i] && a[i] < ans[j][2]) {
            if (b[i] - ans[j][1] < cost) {
              a[ansN] = a[i];
              b[ansN] = ans[j][1];

              ans[ansIdx][0] = a[ansN];
              ans[ansIdx][1] = b[ansN];
              ans[ansIdx][2] = a[i];
              ans[ansIdx][3] = b[i];

              int keepA = ans[j][2];
              int keepB = ans[j][3];
              ans[j][2] = a[ansN];
              ans[j][3] = b[ansN];

              ans[ansCount][0] = a[ansN];
              ans[ansCount][1] = b[ansN];
              ans[ansCount][2] = keepA;
              ans[ansCount][3] = keepB;
              ansCount++;

              ansCost -= cost;
              ansCost += (ans[ansIdx][2] + ans[ansIdx][3]) - (ans[ansIdx][0] + ans[ansIdx][1]);

              ansN++;
              flag = 1;

              break;
            }
          }
        }
      }
    }

    if (flag == 0) break;

    Initialize32(0);
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << ansCount << endl;
    rep(i, ansCount)
    {
      rep(j, 4)
      {
        cout << ans[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    ofs << ansCount << endl;
    rep(i, ansCount)
    {
      rep(j, 4)
      {
        ofs << ans[i][j] << ' ';
      }
      ofs << endl;
    }
  }
}

ll Solve(int probNum)
{
  startTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  ansN     = n;
  ansCost  = 0;
  ansCount = 0;
  Initialize7();
  Initialize6();

  Initialize8();

  Method1();

  CopyToReal();

  // int loop = 0;
  // while (true) {
  //   if (ansN == n * 5) break;

  //   int raA = randxor() % MA;
  //   int raB = randxor() % MA;

  //   if (GetNowTime() > TL) break;
  //   loop++;

  //   a[ansN] = raA;
  //   b[ansN] = raB;
  //   ansN++;
  //   ansCost  = 0;
  //   ansCount = 0;

  //   Initialize3();

  //   if (ansCost < real_ansCost) {
  //     CopyToReal();
  //   } else {
  //     ansN--;
  //   }
  // }

  // if (mode != 0) {
  //   cout << "loop = " << loop << endl;
  // }

  CopyToAns();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
