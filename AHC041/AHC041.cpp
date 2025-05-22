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

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

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
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
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

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
clock_t startTime, endTime;

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int n = 1000;
int m, h;
int a[n];
vector<int> G[n];
int x[n], y[n];

class Answer
{
public:
  vector<int> p;
  int score;

  Answer()
  {
    score = -1;
    p.resize(n);
    for (int i = 0; i < n; ++i) {
      p[i] = -1;
    }
  }

  void Init()
  {
    for (int i = 0; i < n; ++i) {
      p[i] = -1;
    }
    score = -1;
  }
};

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  for (int i = 0; i < n; ++i) G[i].clear();
}

// 入力受け取り
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  if (!ifs.is_open()) {
    int nn;
    cin >> nn >> m >> h;
    for (int i = 0; i < n; ++i) cin >> a[i];
    for (int i = 0; i < m; ++i) {
      int u, v;
      cin >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    for (int i = 0; i < n; ++i) cin >> x[i] >> y[i];
  }
  // ファイル入力する
  else {
    int nn;
    ifs >> nn >> m >> h;
    for (int i = 0; i < n; ++i) ifs >> a[i];
    for (int i = 0; i < m; ++i) {
      int u, v;
      ifs >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    for (int i = 0; i < n; ++i) ifs >> x[i] >> y[i];
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
int hhh[n];
int CalcScore(const Answer& ans)
{
  int res = 0;
  for (int i = 0; i < n; ++i) {
    int hh = 1;
    int x = ans.p[i];
    while (x != -1) {
      hh++;
      x = ans.p[x];
    }
    res += a[i] * hh;
    hhh[i] = hh - 1;
  }
  return res;
}

int f[n] = {};
int hei[n] = {};
int fcount;
int dfsLimit = 40;
void dfs(int x, Answer& ans)
{
  int raFlag = Rand() % 10;
  int raIte = -1;
  if (raFlag == 0) {
    if (G[x].size() >= 2) {
      raIte = Rand() % (G[x].size() - 1);
      swap(G[x][raIte], G[x][raIte + 1]);
    }
  }
  for (auto y : G[x]) {
    if (f[y]) continue;
    if (hei[x] == h - 1 && a[y] <= dfsLimit) continue;
    ans.p[y] = x;
    f[y] = 1;
    fcount++;
    hei[y] = hei[x] + 1;
    if (hei[y] < h) {
      dfs(y, ans);
    }
  }
  if (raIte >= 0) {
    swap(G[x][raIte], G[x][raIte + 1]);
  }
}

// 初期解生成
void Initialize(Answer& ans)
{
  for (int i = 0; i < n; ++i) {
    vector<P> vp;
    for (auto y : G[i]) {
      vp.push_back(P(a[y], y));
    }
    sort(vp.begin(), vp.end());
    G[i].clear();
    for (auto vpp : vp) {
      G[i].push_back(vpp.second);
    }
  }

  Answer best_ans;

  while (true) {
    if (GetNowTime() > TL / 3) {
      break;
    }

    for (int i = 0; i < n; ++i) {
      ans.p[i] = -2;
      f[i] = 0;
      hei[i] = -1;
    }

    fcount = 0;
    while (fcount < n) {
      int ra = Rand() % n;
      while (f[ra] == 1) {
        ra = (ra + 1) % n;
        ra = Rand() % n;
      }

      if (a[ra] >= 40) {
        ra = Rand() % n;
        while (f[ra] == 1) {
          ra = (ra + 1) % n;
          ra = Rand() % n;
        }
      }

      ans.p[ra] = -1;
      f[ra] = 1;
      hei[ra] = 0;
      fcount++;
      dfs(ra, ans);
    }

    if (fcount < n) {
      continue;
    }

    ans.score = CalcScore(ans);
    if (ans.score >= best_ans.score) {
      best_ans = ans;
    }
  }

  ans = best_ans;
}

vector<int> sons[n];
vector<int> roots;
void Method1(Answer& ans)
{
  int loop = 0;
  int flagCount = 0;

  Answer best_ans = ans;

  while (true) {
    if (GetNowTime() > TL) {
      break;
    }

    roots.clear();
    for (int i = 0; i < n; ++i) {
      sons[i].clear();
    }
    for (int i = 0; i < n; ++i) {
      ans.p[i] = best_ans.p[i];
      f[i] = 1;
      if (ans.p[i] == -1) {
        roots.push_back(i);
      }
      else {
        sons[ans.p[i]].push_back(i);
      }
    }

    if (Rand() % 20 == 0) {
      CalcScore(ans);

      int flag = 0;
      for (int i = 0; i < n; ++i) {
        if (sons[i].size() == 0 && hhh[i] < h) {
          for (auto y : G[i]) {
            if (hhh[i] <= hhh[y] + 1 && hhh[y] + 1 <= h) {
              hhh[i] = hhh[y] + 1;
              ans.p[i] = y;
              sons[y].push_back(i);
              flag = 1;
            }
          }
        }
      }

      if (flag) {
        flagCount++;
        ans.score = CalcScore(ans);
        if (ans.score >= best_ans.score) {
          best_ans = ans;
        }
        continue;
      }
    }

    std::shuffle(roots.begin(), roots.end(), engine);
    int raCount = Rand() % 5 + 1;
    fcount = n;
    for (int aespa = 0; aespa < raCount; ++aespa) {
      queue<int> que;
      que.push(roots[aespa]);
      f[roots[aespa]] = 0;
      fcount--;
      while (que.size()) {
        int x = que.front();
        que.pop();
        for (auto y : sons[x]) {
          f[y] = 0;
          fcount--;
          que.push(y);
        }
      }
    }

    while (fcount < n) {
      int ra = Rand() % n;
      while (f[ra] == 1) {
        ra = (ra + 1) % n;
        ra = Rand() % n;
      }

      ans.p[ra] = -1;
      f[ra] = 1;
      hei[ra] = 0;
      fcount++;
      dfs(ra, ans);
    }

    ans.score = CalcScore(ans);
    if (ans.score >= best_ans.score) {
      best_ans = ans;
    }

    loop++;
  }
  if (mode != 0) {
    cout << "loop = " << loop << ", flagCount = " << flagCount << endl;
  }

  ans = best_ans;
}

// 解答出力
void Output(ofstream& ofs, const Answer& ans)
{
  if (mode == 0) {
    for (int i = 0; i < n; ++i) cout << ans.p[i] << ' ';
    cout << endl;
  }
  else {
    for (int i = 0; i < n; ++i) ofs << ans.p[i] << ' ';
    ofs << endl;
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

  Answer ans;

  // 初期解生成
  Initialize(ans);
  Method1(ans);

  // 解答を出力
  Output(ofs, ans);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore(ans);
    int sum[11] = {};
    for (int i = 0; i < n; ++i) sum[hhh[i]]++;
  }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
