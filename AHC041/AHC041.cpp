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
  static uint32_t randxor() {
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
  static double rand01() {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
  }

  // 配列シャッフル
  void FisherYates(int* data, int n) {
    for (int i = n - 1; i >= 0; i--) {
      int j = randxor() % (i + 1);
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

double GetNowTime() {
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int n = 1000;
int m, h;
int a[n];
vector<int> G[n];
int x[n], y[n];

int p[n], best_p[n];
int score, best_score;

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp() {
  rep(i, n) G[i].clear();
  score = -1;
  best_score = -1;
}

// 入力受け取り
void Input(int problemNum) {
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4) {
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
    rep(i, n) cin >> a[i];
    rep(i, m) {
      int u, v;
      cin >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    rep(i, n) cin >> x[i] >> y[i];
  }
  // ファイル入力する
  else {
    int nn;
    ifs >> nn >> m >> h;
    rep(i, n) ifs >> a[i];
    rep(i, m) {
      int u, v;
      ifs >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    rep(i, n) ifs >> x[i] >> y[i];
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs) {
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4) {
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
int CalcScore() {
  int res = 0;
  rep(i, n) {
    int hh = 1;
    int x = p[i];
    while (x != -1) {
      hh++;
      x = p[x];
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
void dfs(int x) {
  int raFlag = randxor() % 10;
  int raIte = -1;
  if (raFlag == 0) {
    if (G[x].size() >= 2) {
      raIte = randxor() % (G[x].size() - 1);
      swap(G[x][raIte], G[x][raIte + 1]);
    }
  }
  for (auto y : G[x]) {
    if (f[y]) continue;
    if (hei[x] == h - 1 && a[y] <= dfsLimit) continue;
    p[y] = x;
    f[y] = 1;
    fcount++;
    hei[y] = hei[x] + 1;
    if (hei[y] < h) {
      dfs(y);
    }
  }
  if (raIte >= 0) {
    swap(G[x][raIte], G[x][raIte + 1]);
  }
}

// 初期解生成
void Initialize() {
  rep(i, n) {
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

  rep(i, n) {
    best_p[i] = -1;
  }

  while (true) {
    if (GetNowTime() > TL / 3) {
      break;
    }
    // dfsLimit = randxor() % 10 * 10;
    rep(i, n) {
      p[i] = -2;
      f[i] = 0;
      hei[i] = -1;
    }

    fcount = 0;
    while (fcount < n) {
      int ra = randxor() % n;
      while (f[ra] == 1) {
        ra = (ra + 1) % n;
        ra = randxor() % n;
      }

      if (a[ra] >= 40) {
        ra = randxor() % n;
        while (f[ra] == 1) {
          ra = (ra + 1) % n;
          ra = randxor() % n;
        }
      }

      p[ra] = -1;
      f[ra] = 1;
      hei[ra] = 0;
      fcount++;
      dfs(ra);
    }

    if (fcount < n) {
      continue;
    }

    score = CalcScore();
    if (score >= best_score) {
      best_score = score;
      rep(i, n) {
        best_p[i] = p[i];
      }
    }
  }

  score = best_score;
  rep(i, n) {
    p[i] = best_p[i];
  }
}

vector<int> sons[n];
vector<int> roots;
void Method1() {
  int loop = 0;
  int flagCount = 0;

  while (true) {
    if (GetNowTime() > TL) {
      break;
    }
    // dfsLimit = randxor() % 10 * 10;

    roots.clear();
    rep(i, n) {
      sons[i].clear();
    }
    rep(i, n) {
      p[i] = best_p[i];
      f[i] = 1;
      if (p[i] == -1) {
        roots.push_back(i);
      }
      else {
        sons[p[i]].push_back(i);
      }
    }

    if (randxor() % 20 == 0) {
      CalcScore();

      int flag = 0;
      rep(i, n) {
        if (sons[i].size() == 0 && hhh[i] < h) {
          for (auto y : G[i]) {
            if (hhh[i] <= hhh[y] + 1 && hhh[y] + 1 <= h) {
              hhh[i] = hhh[y] + 1;
              p[i] = y;
              sons[y].push_back(i);
              flag = 1;
            }
          }
        }
      }

      if (flag) {
        flagCount++;
        score = CalcScore();
        if (score >= best_score) {
          best_score = score;
          rep(i, n) {
            best_p[i] = p[i];
          }
        }
        continue;
      }
    }

    std::shuffle(roots.begin(), roots.end(), engine);
    int raCount = randxor() % 5 + 1;
    fcount = n;
    rep(aespa, raCount) {
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
      int ra = randxor() % n;
      while (f[ra] == 1) {
        ra = (ra + 1) % n;
        ra = randxor() % n;
      }

      // if (a[ra] >= 40) {
      //   ra = randxor() % n;
      //   while (f[ra] == 1) {
      //     ra = (ra + 1) % n;
      //     ra = randxor() % n;
      //   }
      // }

      p[ra] = -1;
      f[ra] = 1;
      hei[ra] = 0;
      fcount++;
      dfs(ra);
    }

    score = CalcScore();
    if (score >= best_score) {
      best_score = score;
      rep(i, n) {
        best_p[i] = p[i];
      }
    }

    loop++;
  }
  if (mode != 0) {
    cout << "loop = " << loop << ", flagCount = " << flagCount << endl;
  }

  score = best_score;
  rep(i, n) {
    p[i] = best_p[i];
  }
}

// 解答出力
void Output(ofstream& ofs) {
  if (mode == 0) {
    rep(i, n) cout << p[i] << ' ';
    cout << endl;
  }
  else {
    rep(i, n) ofs << p[i] << ' ';
    ofs << endl;
  }
}

ll Solve(int probNum) {
  startTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();
  Method1();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
    int sum[11] = {};
    rep(i, n) sum[hhh[i]]++;
    rep(i, 11) {
      // cout << i << " : " << sum[i] << endl;
    }
  }
  return score;
}

int main() {
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
    srep(i, 0, 15) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
