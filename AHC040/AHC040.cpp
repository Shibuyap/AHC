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

const int MAX_N = 100;
const int MAX_T = 400;

int n, t, sigma;
int w[MAX_N], h[MAX_N];

int W[MAX_N], H[MAX_N];
int dW[MAX_T], dH[MAX_T];

int queryCount;
int tScore[MAX_T];

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

int CalcTScore() {
  int score = INF;
  rep(i, queryCount) {
    score = min(score, tScore[i]);
  }
  return score;
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

  tScore[queryCount] = ww + hh;
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

  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30 * 2) {
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
    rep(i, n) {
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

vector<Column> bestbest;
bool ff[MAX_N * 2 + MAX_T][MAX_N * 2];
int nowSum[MAX_N * 2 + MAX_T];
int keisoku[MAX_N * 2 + MAX_T];
void Yamanobori(ofstream& ofs) {
  vector<Column> best(n);
  int bestScore = INF;

  vector<Column> tmp(n);
  int tmpScore = INF;
  int loop = 0;

  double nowTime = GetNowTime();  // 現在の経過時間
  const double START_TEMP = 200.0;  // 焼きなまし法の開始温度
  const double END_TEMP = 0.1;      // 焼きなまし法の終了温度
  double temp = START_TEMP + (END_TEMP - START_TEMP) * nowTime / TL;  // 現在の温度

  int originH[MAX_N], originW[MAX_N];
  rep(i, n) {
    originH[i] = h[i];
    originW[i] = w[i];
  }


  rep(i, n * 2 + t) {
    rep(j, n * 2) {
      ff[i][j] = false;
    }
  }
  rep(i, n * 2) {
    ff[i][i] = true;
  }

  rep(i, n * 2) {
    nowSum[i * 2] = h[i];
    nowSum[i * 2 + 1] = w[i];
    keisoku[i * 2] = h[i];
    keisoku[i * 2 + 1] = w[i];
  }

  t = MAX_T;
  rep(aespa, t - 1) {
    vector<Column> cols;
    int pred = 0;
    rep(i, n) {
      if (RandXor() % 2) {
        Column col;
        col.piece = i;
        col.rot = RandXor() % 2;
        //col.rot = 0;
        col.dir = 0;
        col.base = -1;
        if (col.rot == 0) {
          pred += h[i];
        }
        else {
          pred += w[i];
        }
        cols.push_back(col);
      }
    }

    int ww, hh;
    Print(cols, ww, hh, ofs);

    keisoku[n * 2 + aespa] = hh;
    for (auto col : cols) {
      int num = col.piece * 2;
      if (col.rot == 1)num++;
      ff[n * 2 + aespa][num] = true;
      if (col.rot == 0) {
        nowSum[n * 2 + aespa] += h[col.piece];
      }
      else {
        nowSum[n * 2 + aespa] += w[col.piece];
      }
    }
  }

  while (true) {
    if (RandXor() % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL) {
        break;
      }
    }

    int raP = RandXor() % (n * 2);
    int raDiff = RandXor() % 201 - 100;
    if (raP % 2 == 0) {
      if (h[raP / 2] + raDiff < 10000)continue;
      if (h[raP / 2] + raDiff > 100000)continue;
    }
    else {
      if (w[raP / 2] + raDiff < 10000)continue;
      if (w[raP / 2] + raDiff > 100000)continue;
    }
    int beforeDiffSum = 0;
    int afterDiffSum = 0;

    beforeDiffSum += abs(keisoku[raP] - nowSum[raP]);
    nowSum[raP] += raDiff;
    afterDiffSum += abs(keisoku[raP] - nowSum[raP]);
    srep(aespa, n * 2, n * 2 + t - 1) {
      if (ff[aespa][raP]) {
        beforeDiffSum += abs(keisoku[aespa] - nowSum[aespa]);
        nowSum[aespa] += raDiff;
        afterDiffSum += abs(keisoku[aespa] - nowSum[aespa]);
      }
    }

    if (afterDiffSum <= beforeDiffSum) {
      if (raP % 2 == 0) {
        h[raP / 2] += raDiff;
      }
      else {
        w[raP / 2] += raDiff;
      }
    }
    else {
      nowSum[raP] -= raDiff;
      srep(aespa, n * 2, n * 2 + t - 1) {
        if (ff[aespa][raP]) {
          nowSum[aespa] -= raDiff;
        }
      }
    }
  }

  rep(i, n) {
    cout << originH[i] << ' ' << h[i] << ' ' << H[i] << endl;
  }

  while (loop < 100) {
    loop++;
    rep(i, n) {
      tmp[i].piece = i;
      tmp[i].rot = RandXor() % 2;
      tmp[i].dir = RandXor() % 2;
      tmp[i].base = RandXor() % (i + 1) - 1;
    }

    int width = sqrt(n * 2);
    rep(i, n) {
      tmp[i].piece = i;
      tmp[i].rot = 0;
      if (w[i] > h[i]) tmp[i].rot = 1;
      tmp[i].dir = 0;
      if (i % width == 0) {
        tmp[i].base = -1;
      }
      else {
        tmp[i].base = i - 1;
      }
    }

    int ww, hh;
    tmpScore = CalcScore(tmp, false, ww, hh);
    if (tmpScore <= bestScore) {
      bestScore = tmpScore;
      best = tmp;
    }
  }

  tmp = best;
  tmpScore = bestScore;

  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL * 3) {
        break;
      }
    }

    loop++;

    if (RandXor() % 100000 == 0) {
      tmp = best;
      tmpScore = bestScore;
    }

    int raP = RandXor() % n;

    Column keep = tmp[raP];

    tmp[raP].rot = RandXor() % 2;
    tmp[raP].dir = RandXor() % 2;
    tmp[raP].base = RandXor() % (raP + 1) - 1;

    int ww, hh;
    int tmp2Score = CalcScore(tmp, false, ww, hh);

    const double progressRatio = nowTime / TL;  // 進捗率（0.0～1.0）
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;  // 温度の更新

    double diff = tmpScore - tmp2Score;  // スコアの差分
    double prob = exp(diff / temp);     // 焼きなまし法の採用確率

    if (prob > Rand01()) {
      tmpScore = tmp2Score;
      if (tmpScore < bestScore) {
        bestScore = tmpScore;
        best = tmp;
      }
    }
    else {
      tmp[raP] = keep;
    }
  }

  int ww, hh;
  cout << "loop = " << loop << endl;
  cout << bestScore << ' ' << CalcScore(best, true, ww, hh) << endl;

  Print(best, ww, hh, ofs);

  bestbest = best;
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
  //Yamanobori(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    int ww, hh;
    score = CalcTScore();
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
