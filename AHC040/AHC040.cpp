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

// ループの簡略化マクロ
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t RandXor() {
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() {
  return (RandXor() + 0.5) * (1.0 / UINT_MAX);
}

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double rot) {
  return l + (rot - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n) {
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化（使用されていない）
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// 非常に大きな値
const int INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 2.8;                                       // 時間制限（Time Limit）
int mode;                                              // 実行モード
std::chrono::steady_clock::time_point startTimeClock;  // 時間計測用

// 時間計測をリセットする関数
void ResetTime() {
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime() {
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

// 区間の交差判定
bool isCrossing(double l1, double r1, double l2, double r2) {
  return (std::max(l1, l2) < std::min(r1, r2));
}

struct Column {
  int piece;
  int rot;
  int dir;
  int base;
};

struct Score {
  int score;
  int ww;
  int hh;
};

const int MAX_N = 100;
const int MAX_T = 400;

int n, t, sigma;
int w[MAX_N], h[MAX_N];

int W[MAX_N], H[MAX_N];
int dW[MAX_T], dH[MAX_T];

int queryCount;
Score tScore[MAX_T];

void CopyToBest() {
}

void CopyToAns() {
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp() {
  queryCount = 0;
}

// 入力を受け取る関数
void Input(int problemNum) {
  if (mode == 0) {
    // 標準入力
    cin >> n >> t >> sigma;
    rep(i, n) {
      cin >> w[i] >> h[i];
    }
  }
  else {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
    ifstream ifs(oss.str());
    // ファイル入力
    ifs >> n >> t >> sigma;
    rep(i, n) {
      ifs >> w[i] >> h[i];
    }
    rep(i, n) {
      ifs >> W[i] >> H[i];
    }
  }

  rep(i, n) {
    w[i] = max(10000, w[i]);
    w[i] = min(100000, w[i]);
    h[i] = max(10000, h[i]);
    h[i] = min(100000, h[i]);
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs) {
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

Score CalcTScore() {
  Score res;
  res.score = INF;
  rep(i, queryCount) {
    if (tScore[i].score < res.score) {
      res = tScore[i];
    }
  }
  return res;
}

int cs_use[MAX_N] = {};
int cs_up[MAX_N], cs_down[MAX_N], cs_left[MAX_N], cs_right[MAX_N];
// スコアを計算する関数
int CalcScore(const vector<Column>& columns, bool cheat, int& ww, int& hh) {
  rep(i, MAX_N) {
    cs_up[i] = -1;
    cs_down[i] = -1;
    cs_left[i] = -1;
    cs_right[i] = -1;
  }
  for (auto col : columns) {
    int num = col.piece;
    int wid = w[num];
    int hei = h[num];
    if (mode != 0 && cheat == true) {
      wid = W[num];
      hei = H[num];
    }
    if (col.rot == 1) swap(wid, hei);
    if (col.dir == 0) {
      cs_left[num] = 0;
      if (col.base != -1) {
        cs_left[num] = cs_right[col.base];
      }
      cs_right[num] = cs_left[num] + wid;

      cs_up[num] = 0;
      rep(i, num) {
        if (cs_use[i]) {
          if (isCrossing(cs_left[num], cs_right[num], cs_left[i], cs_right[i])) {
            cs_up[num] = max(cs_up[num], cs_down[i]);
          }
        }
      }
      cs_down[num] = cs_up[num] + hei;
    }
    else {
      cs_up[num] = 0;
      if (col.base != -1) {
        cs_up[num] = cs_down[col.base];
      }
      cs_down[num] = cs_up[num] + hei;

      cs_left[num] = 0;
      rep(i, num) {
        if (cs_use[i]) {
          if (isCrossing(cs_up[num], cs_down[num], cs_up[i], cs_down[i])) {
            cs_left[num] = max(cs_left[num], cs_right[i]);
          }
        }
      }
      cs_right[num] = cs_left[num] + wid;
    }

    cs_use[num] = 1;
  }

  int maxDown = 0, maxRight = 0;
  rep(i, n) {
    if (cs_use[i]) {
      maxDown = max(maxDown, cs_down[i]);
      maxRight = max(maxRight, cs_right[i]);
    }
  }

  ww = maxRight;
  hh = maxDown;

  int score = maxDown + maxRight;
  rep(i, n) {
    if (cs_use[i] == 0) {
      score += w[i] + h[i];
    }
  }

  return score;
}

void Print(const vector<Column>& columns, int& ww, int& hh, ofstream& ofs) {
  ww = 0;
  hh = 0;

  if (mode == 0) {
    cout << columns.size() << endl;
    rep(i, columns.size()) {
      cout << columns[i].piece << ' ' << columns[i].rot << ' ' << (columns[i].dir == 0 ? 'U' : 'L') << ' ' << columns[i].base << endl;
    }
    fflush(stdout);

    cin >> ww >> hh;
  }
  else {
    ofs << columns.size() << endl;
    rep(i, columns.size()) {
      ofs << columns[i].piece << ' ' << columns[i].rot << ' ' << (columns[i].dir == 0 ? 'U' : 'L') << ' ' << columns[i].base << endl;
    }

    CalcScore(columns, true, ww, hh);
    ww += dW[queryCount];
    hh += dH[queryCount];
  }


  tScore[queryCount].hh = hh;
  tScore[queryCount].ww = ww;
  tScore[queryCount].score = hh + ww;
  queryCount++;
}

vector<Column> Method2_Shoki2() {
  vector<Column> best(n);
  int bestScore = INF;

  vector<Column> tmp(n);
  rep(i, n) {
    tmp[i].piece = i;
    tmp[i].rot = 0;
    if (w[i] > h[i]) tmp[i].rot = 1;
    tmp[i].dir = 0;
    tmp[i].base = i - 1;
  }

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30) {
        break;
      }
    }

    loop++;

    int widSum = 0;
    int heiSum = 0;

    int widLimit = RandXor() % 1000000 + 200000;
    int now = 0;
    int maxHeight = 0;
    rep(i, n) {
      int wid = w[i];
      int hei = h[i];
      if (tmp[i].rot == 1) {
        wid = h[i];
        hei = w[i];
      }

      if (now + wid <= widLimit) {
        tmp[i].base = i - 1;
        now += wid;
        maxHeight = max(maxHeight, hei);
      }
      else {
        tmp[i].base = -1;
        widSum = max(widSum, now);
        now = wid;
        heiSum += maxHeight;
        maxHeight = hei;
      }
    }

    widSum = max(widSum, now);
    heiSum += maxHeight;

    if (widSum + heiSum < bestScore) {
      bestScore = widSum + heiSum;
      best = tmp;
    }
  }

  int keepRot[MAX_N];
  int keepRotCount = 0;
  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30 * 25) {
        break;
      }
    }

    loop++;

    int widSum = 0;
    int heiSum = 0;

    int widLimit = RandXor() % 1000000 + 200000;

    int now = 0;
    int last = -1;
    int maxHeight = 0;

    int beforeNow = INF;
    int beforeLast = -1;
    int beforeMaxHeight = 0;

    int isRandomRot = RandXor() % 10;

    keepRotCount = 0;
    rep(i, n) {
      if (isRandomRot >= 1 && RandXor() % n <= isRandomRot) {
        tmp[i].rot = 1 - tmp[i].rot;
        keepRot[keepRotCount] = i;
        keepRotCount++;
      }

      int wid = w[i];
      int hei = h[i];
      if (tmp[i].rot == 1) {
        wid = h[i];
        hei = w[i];
      }

      if (now < widLimit * 0.9 && beforeNow + wid <= widLimit) {
        tmp[i].base = beforeLast;
        beforeLast = i;
        beforeNow += wid;
        widSum = max(widSum, beforeNow);
        if (hei > beforeMaxHeight) {
          heiSum += hei - beforeMaxHeight;
          beforeMaxHeight = hei;
        }
      }
      else if (now + wid <= widLimit) {
        tmp[i].base = last;
        last = i;
        now += wid;
        maxHeight = max(maxHeight, hei);
      }
      else {
        beforeLast = last;
        beforeMaxHeight = maxHeight;
        beforeNow = now;

        tmp[i].base = -1;
        last = i;
        widSum = max(widSum, now);
        now = wid;
        heiSum += maxHeight;
        maxHeight = hei;
      }
    }

    widSum = max(widSum, now);
    heiSum += maxHeight;

    if (widSum + heiSum < bestScore) {
      bestScore = widSum + heiSum;
      best = tmp;
    }

    rep(i, keepRotCount) {
      int ii = keepRot[i];
      tmp[ii].rot = 1 - tmp[ii].rot;
    }
  }

  return best;
}

void Method2(ofstream& ofs) {
  vector<Column> best(n);
  int bestScore = INF;

  vector<Column> tmp(n);

  rep(aespa, t) {
    int raMode = RandXor() % 2;
    if (aespa == 0) {
      tmp = Method2_Shoki2();
    }
    else {
      tmp = best;
      if (raMode == 0) {
        int ra = RandXor() % n;
        tmp[ra].rot = 1 - tmp[ra].rot;
      }
      else {
        vector<vector<Column>> vvc;
        vector<Column> vc;
        for (auto col : best) {
          if (col.base == -1 && !vc.empty()) {
            vvc.push_back(vc);
            vc.clear();
          }
          vc.push_back(col);
        }
        vvc.push_back(vc);
        if (vvc.size() >= 2) {
          while (true) {
            int ra1 = RandXor() % (vvc.size() - 1);
            int ra2 = RandXor() % 2;
            if (ra2 == 0) {
              if (vvc[ra1].size() >= 2) {
                vvc[ra1 + 1].insert(vvc[ra1 + 1].begin(), vvc[ra1].back());
                vvc[ra1].pop_back();
                break;
              }
            }
            else {
              if (vvc[ra1 + 1].size() >= 2) {
                vvc[ra1].push_back(vvc[ra1 + 1][0]);
                vvc[ra1 + 1].erase(vvc[ra1 + 1].begin());
                break;
              }
            }
          }
        }
        tmp.clear();
        rep(i, vvc.size()) {
          rep(j, vvc[i].size()) {
            Column col = vvc[i][j];
            if (j == 0) {
              col.base = -1;
            }
            else {
              col.base = vvc[i][j - 1].piece;
            }
            tmp.push_back(col);
          }
        }
      }
    }

    int ww, hh;
    Print(tmp, ww, hh, ofs);

    int tmpScore = ww + hh;
    if (tmpScore < bestScore) {
      if (mode != 0 && aespa > 0) {
        //cout << "turn = " << aespa + 1 << ", raMode = " << raMode << endl;
      }
      bestScore = tmpScore;
      best = tmp;
    }
  }
}

// 問題を解く関数
ll Solve(int problem_num) {
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // 初期解生成
  Method2(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  int score = 0;
  if (mode != 0) {
    Score sc = CalcTScore();
    cout << "ww = " << sc.ww << ", hh = " << sc.hh << ", score = " << sc.score << endl;
    score = sc.score;
  }
  return score;
}

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main() {
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 3) {
    rep(_, 10) {
      ll sum = 0;
      srep(i, 0, 100) {
        ll score = Solve(i);
        sum += score;
      }
      cout << sum << endl;
    }
  }
  else {
    ll sum = 0;
    srep(i, 0, 100) {
      ll score = Solve(i);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
