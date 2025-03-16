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
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t RandXor()
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() { return (RandXor() + 0.5) * (1.0 / UINT_MAX); }

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
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
const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.95; // 時間制限（Time Limit）
int mode;        // 実行モード
std::chrono::steady_clock::time_point startTimeClock; // 時間計測用

// 時間計測をリセットする関数
void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

//const int MAX_N = 30;

const int n = 100;
const int L = 500000;

int t[n];
int searchT[n];

int ansScore;
int ansScore2;
int a[n], b[n];
double cntDouble[n];
double cntDoubleAns[n];

int best_ansScore;
int best_ansScore2;
int best_a[n], best_b[n];

void CopyToBest()
{
  best_ansScore = ansScore;
  best_ansScore2 = ansScore2;
  rep(i, n) {
    best_a[i] = a[i];
    best_b[i] = b[i];
  }
}

void CopyToAns()
{
  ansScore = best_ansScore;
  ansScore2 = best_ansScore2;
  rep(i, n) {
    a[i] = best_a[i];
    b[i] = best_b[i];
  }
}

bool IsNG(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  ansScore2 = 0;
  best_ansScore = 0;
  best_ansScore2 = 0;
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nnn, lll;
    cin >> nnn >> lll;
    rep(i, n)cin >> t[i];
  }
  else {
    // ファイル入力
    int nnn, lll;
    ifs >> nnn >> lll;
    rep(i, n)ifs >> t[i];
  }

  vector<P> vp;
  rep(i, n) {
    vp.push_back(P(t[i], i));
  }
  sort(vp.begin(), vp.end());
  rep(i, n) {
    t[i] = vp[i].first;
    searchT[i] = vp[i].second;
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore()
{
  int cnt[n] = {};
  int now = 0;
  rep(i, L) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
  }

  int res = 1000000;
  rep(i, n) {
    res -= abs(cnt[i] - t[i]);
  }
  return res;
}

int CalcScoreEasy(int k)
{
  int cnt[n] = {};
  int now = 0;
  rep(i, k) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
  }

  double res = 1000000;
  double mul = (double)L / k;
  rep(i, n) {
    res -= abs(cnt[i] * mul - t[i]);
  }
  return max(0, (int)round(res));
}

int sortedA[n];
int sortedB[n];
vector<P> sortVec(n);
int convertArr[n];
int SortEasy(int k) {
  int cnt[n] = {};
  int now = 0;
  rep(i, k) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
    if (i == 10000) {
      rep(j, n) {
        if (t[j] > 100 && cnt[j] == 0) {
          return 0;
        }
      }
    }
  }

  rep(i, n) {
    sortVec[i].first = cnt[i];
    sortVec[i].second = i;
  }
  sort(sortVec.begin(), sortVec.end());
  rep(i, n) {
    convertArr[sortVec[i].second] = i;
  }

  rep(i, n) {
    sortedA[i] = convertArr[a[sortVec[i].second]];
    sortedB[i] = convertArr[b[sortVec[i].second]];
  }

  rep(i, n) {
    if (sortVec[i].first > 0)break;
    if (t[i] > 100) {
      return 0;
    }
  }

  double mul = (double)L / k;
  rep(i, n) {
    cntDouble[i] = sortVec[i].first * mul;
  }

  double res = 1000000;
  rep(i, n) {
    res -= abs(cntDouble[i] - t[i]);
  }
  return max(0, (int)round(res));
}

double earlyCheckArr[n];
int earlyCheckNumArr[1000];
int earlyCheckNumArr2[1000];
void ClearEarlyCheckArr() {
  rep(i, n) {
    earlyCheckArr[i] = cntDoubleAns[i];
  }
}
void UpdateEarlyCheck(int s, double diff) {
  earlyCheckNumArr[0] = s;
  int tail = 1;
  rep(dfs, 6) {
    int newTail = 0;
    rep(j, tail) {
      int num = earlyCheckNumArr[j];
      earlyCheckArr[num] += diff;
      earlyCheckNumArr2[newTail] = a[num];
      newTail++;
      earlyCheckNumArr2[newTail] = b[num];
      newTail++;
    }
    rep(j, newTail) {
      earlyCheckNumArr[j] = earlyCheckNumArr2[j];
    }
    tail = newTail;
    diff /= 2;
  }
}
double CalcEarlyCheck() {
  double res = 1000000;
  rep(i, n) {
    res -= abs(earlyCheckArr[i] - t[i]);
  }
  return max(0, (int)round(res));
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  int aaaa[n] = {};
  int bbbb[n] = {};
  rep(i, n) {
    aaaa[searchT[i]] = searchT[a[i]];
    bbbb[searchT[i]] = searchT[b[i]];
  }

  if (mode == 0) {
    // 標準出力
    rep(i, n)cout << aaaa[i] << ' ' << bbbb[i] << endl;
  }
  else {
    // ファイル出力
    rep(i, n)ofs << aaaa[i] << ' ' << bbbb[i] << endl;
  }
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition[10];
};

void SimulatedAnnealing(Hypers hypers)
{
  CopyToBest();

  // ランダムに初期解作成
  double nowTime = GetNowTime();
  int loop1 = 0;
  while (true) {
    loop1++;

    if (loop1 % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL / 10) break;
    }

    if (RandXor() % 2 == 0) {
      rep(i, n) {
        a[i] = RandXor() % n;
        b[i] = RandXor() % n;
      }
    }
    else {
      rep(i, n) {
        if (RandXor() % 2 == 0) {
          a[i] = RandXor() % n;
          b[i] = (i + 1) % n;
        }
        else {
          a[i] = (i + 1) % n;
          b[i] = RandXor() % n;
        }
      }
    }

    int tmpScore = SortEasy(10000);
    if (tmpScore > ansScore) {
      ansScore = tmpScore;
      rep(i, n) {
        a[i] = sortedA[i];
        b[i] = sortedB[i];
      }
      CopyToBest();
    }
  }


  CopyToAns();
  int k = 25000;
  ansScore = SortEasy(k);
  rep(i, n) {
    a[i] = sortedA[i];
    b[i] = sortedB[i];
    cntDoubleAns[i] = cntDouble[i];
  }
  CopyToBest();


  if (mode != 0 && mode != 3) {
    cout << loop1 << endl;
  }

  int saitakuCount[10][2];
  rep(i, 10) {
    rep(j, 2) {
      saitakuCount[i][j] = 0;
    }
  }

  nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;
  int loop2 = 0;
  while (true) {
    loop2++;

    if (loop2 % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    int ok = 1;

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
    int NEAR = 5;
    int raMode = RandXor() % hypers.Partition[9];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      saitakuCount[0][1]++;
      ra1 = RandXor() % n;
      if (RandXor() % 2 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }

      ClearEarlyCheckArr();
      UpdateEarlyCheck(a[ra1], -cntDoubleAns[ra1] / 2.0);

      keep1 = a[ra1];
      a[ra1] = ra2;

      UpdateEarlyCheck(ra2, cntDoubleAns[ra1] / 2.0);
      double earlyScore = CalcEarlyCheck();
      if (earlyScore < ansScore - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.Partition[1]) {
      saitakuCount[1][1]++;
      ra1 = RandXor() % n;
      if (RandXor() % 2 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }

      ClearEarlyCheckArr();
      UpdateEarlyCheck(b[ra1], -cntDoubleAns[ra1] / 2.0);

      keep1 = b[ra1];
      b[ra1] = ra2;

      UpdateEarlyCheck(ra2, cntDoubleAns[ra1] / 2.0);
      double earlyScore = CalcEarlyCheck();
      if (earlyScore < ansScore - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.Partition[2]) {
      saitakuCount[2][1]++;
      ra1 = RandXor() % n;
      if (RandXor() % 10 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }

      ClearEarlyCheckArr();
      UpdateEarlyCheck(a[ra1], -cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(a[ra2], -cntDoubleAns[ra2] / 2.0);

      swap(a[ra1], a[ra2]);

      UpdateEarlyCheck(a[ra1], cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(a[ra2], cntDoubleAns[ra2] / 2.0);
      double earlyScore = CalcEarlyCheck();
      if (earlyScore < ansScore - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.Partition[3]) {
      saitakuCount[3][1]++;
      ra1 = RandXor() % n;
      if (RandXor() % 10 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }

      ClearEarlyCheckArr();
      UpdateEarlyCheck(b[ra1], -cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(b[ra2], -cntDoubleAns[ra2] / 2.0);

      swap(b[ra1], b[ra2]);

      UpdateEarlyCheck(b[ra1], cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(b[ra2], cntDoubleAns[ra2] / 2.0);
      double earlyScore = CalcEarlyCheck();
      if (earlyScore < ansScore - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.Partition[4]) {
      saitakuCount[4][1]++;
      ra1 = RandXor() % n;
      swap(a[ra1], b[ra1]);
    }
    else if (raMode < hypers.Partition[5]) {
      saitakuCount[5][1]++;
      ra1 = RandXor() % n;
      if (RandXor() % 10 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }

      ClearEarlyCheckArr();
      UpdateEarlyCheck(a[ra1], -cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(b[ra2], -cntDoubleAns[ra2] / 2.0);

      swap(a[ra1], b[ra2]);

      UpdateEarlyCheck(a[ra1], cntDoubleAns[ra1] / 2.0);
      UpdateEarlyCheck(b[ra2], cntDoubleAns[ra2] / 2.0);
      double earlyScore = CalcEarlyCheck();
      if (earlyScore < ansScore - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.Partition[6]) {
      saitakuCount[6][1]++;
      ra1 = RandXor() % (n);
      if (RandXor() % 2 == 0) {
        ra2 = RandXor() % n;
      }
      else {
        ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + RandXor() % (NEAR * 2 + 1) - NEAR;
        }
      }
      rep(i, n) {
        if (a[i] == ra1) {
          a[i] = ra2;
        }
        else if (a[i] == ra2) {
          a[i] = ra1;
        }
        if (b[i] == ra1) {
          b[i] = ra2;
        }
        else if (b[i] == ra2) {
          b[i] = ra1;
        }
      }
    }

    // スコア計算
    double tmpScore2 = 0;
    if (ok) {
      tmpScore2 = SortEasy(k);

      // 焼きなまし
      double diffScore2 = (tmpScore2 - ansScore) * hypers.MultipleValue;
      double prob2 = exp(diffScore2 / temp);
      ok = prob2 > Rand01();
    }

    if (ok) {
      // 採用
      ansScore = tmpScore2;
      //ansScore2 = tmpScore2;
      rep(i, n) {
        a[i] = sortedA[i];
        b[i] = sortedB[i];
        cntDoubleAns[i] = cntDouble[i];
      }

      // Best解よりもいいか
      if (ansScore > best_ansScore) {
        //cout << nowTime << ' ' << ansScore << endl;
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {
        saitakuCount[0][0]++;
        a[ra1] = keep1;
      }
      else if (raMode < hypers.Partition[1]) {
        saitakuCount[1][0]++;
        b[ra1] = keep1;
      }
      else if (raMode < hypers.Partition[2]) {
        saitakuCount[2][0]++;
        swap(a[ra1], a[ra2]);
      }
      else if (raMode < hypers.Partition[3]) {
        saitakuCount[3][0]++;
        swap(b[ra1], b[ra2]);
      }
      else if (raMode < hypers.Partition[4]) {
        saitakuCount[4][0]++;
        swap(a[ra1], b[ra1]);
      }
      else if (raMode < hypers.Partition[5]) {
        saitakuCount[5][0]++;
        swap(a[ra1], b[ra2]);
      }
      else if (raMode < hypers.Partition[6]) {
        saitakuCount[6][0]++;
        rep(i, n) {
          if (a[i] == ra1) {
            a[i] = ra2;
          }
          else if (a[i] == ra2) {
            a[i] = ra1;
          }
          if (b[i] == ra1) {
            b[i] = ra2;
          }
          else if (b[i] == ra2) {
            b[i] = ra1;
          }
        }
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop2 << endl;
    rep(i, 10) {
      cout << saitakuCount[i][1] - saitakuCount[i][0] << " / " << saitakuCount[i][1] << endl;
    }
  }

  CopyToAns();
}


// 問題を解く関数
ll Solve(int problem_num, Hypers hypers)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // 焼きなまし
  SimulatedAnnealing(hypers);

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

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;

  Hypers HYPERS;
  HYPERS.StartTemp = 2000000.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 12345.0;
  HYPERS.Partition[0] = 100;
  HYPERS.Partition[1] = 200;
  HYPERS.Partition[2] = 300;
  HYPERS.Partition[3] = 400;
  HYPERS.Partition[4] = 440;
  HYPERS.Partition[5] = 700;
  HYPERS.Partition[6] = 700;
  HYPERS.Partition[7] = 700;
  HYPERS.Partition[8] = 700;
  HYPERS.Partition[9] = 700;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 15)
    {
      ll score = Solve(i, HYPERS);
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
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition[0] = RandXor() % 101;

      ll sum = 0;
      srep(i, 0, 4)
      {
        ll score = Solve(i, hypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition[0]
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
