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
typedef pair<int, int> PAIR;
const int INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

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


// 汎用変数
int mode;
std::mt19937 engine;
std::normal_distribution<> dist;

// 入力用変数
int L, N, S, LL, SS;
int Y[110], X[110];
int A[110];
const int F_SIZE = 10000;
int f[11000];

// 本番解答用変数
int keisokuCount;
int P[55][55];
int E[110];

// InputFile用変数
ll inputFileKeisokuCost;

// ハイパラ
const double TL = 3.8;
const int SetSize = 20;
int TansakuSize = 9;

// 焼きなまし用
int AA[110];
int f_SA[110][SetSize][11][11];

ll maxScore;
ll maxHaitiCost;
ll maxKeisokuCost;
int EE[SetSize][110];
int DIFFS[SetSize][110][110];

int KEEP_EE[SetSize][110];
int KEEP_DIFFS[SetSize][110][110];

ll real_maxScore;
ll real_maxHaitiCost;
ll real_maxKeisokuCost;
int real_P[55][55];
int real_EE[SetSize][110];
int real_DIFFS[SetSize][110][110];

// 情報
int MethodCount[20][10];

void InitDist()
{
  // 平均0.0、標準偏差Sで分布させる
  std::normal_distribution<>::param_type param(0.0, S);
  dist.param(param);
}

inline int Dist()
{
  return round(dist(engine));
}

void Init()
{
  keisokuCount = 0;
  rep(i, 20)
  {
    rep(j, 10)
    {
      MethodCount[i][j] = 0;
    }
  }
}

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
    cin >> L >> N >> S;
    rep(i, N) { cin >> Y[i] >> X[i]; }
  }
  // ファイル入力する
  else {
    ifs >> L >> N >> S;
    rep(i, N) { ifs >> Y[i] >> X[i]; }
    rep(i, N) { ifs >> A[i]; }
    rep(i, F_SIZE) { ifs >> f[i]; }
  }


  LL = L - 10;
  SS = 0;
  srep(i, 1, 31)
  {
    if (i * i == S) {
      SS = i - 1;
      break;
    }
  }

  InitDist();
}

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

ll CalcHaitiCost()
{
  ll haitiCost = 0;
  rep(i, L)
  {
    rep(j, L)
    {
      haitiCost += (P[i][j] - P[(i + 1) % L][j]) * (P[i][j] - P[(i + 1) % L][j]);
      haitiCost += (P[i][j] - P[i][(j + 1) % L]) * (P[i][j] - P[i][(j + 1) % L]);
    }
  }
  return haitiCost;
}

void InitHaiti()
{
  rep(i, L)
  {
    rep(j, L) { P[i][j] = 0; }
  }
}

// ランダムに0か1000
void InitHaiti2()
{
  rep(i, L)
  {
    rep(j, L) { P[i][j] = randxor() % 2 * 1000; }
  }
}

void InitHaiti3()
{
  rep(i, L)
  {
    rep(j, L) { P[i][j] = 500; }
  }
}

// 一様ランダム
void InitHaiti4()
{
  rep(i, L)
  {
    rep(j, L) { P[i][j] = randxor() % 1001; }
  }
}

void InitHaiti5()
{
  rep(i, L)
  {
    rep(j, L) { P[i][j] = 250 + randxor() % 2 * 500; }
  }
}

void PrintHaiti(ofstream& ofs)
{
  if (mode == 0) {
    rep(i, L)
    {
      rep(j, L) { cout << P[i][j] << ' '; }
      cout << endl;
      fflush(stdout);
    }
  }
  else {
    rep(i, L)
    {
      rep(j, L) { ofs << P[i][j] << ' '; }
      ofs << endl;
    }
  }
}

int Keisoku(int i, int y, int x, ofstream& ofs)
{
  if (mode == 0) {
    cout << i << ' ' << y << ' ' << x << endl;
    fflush(stdout);
  }
  else {
    ofs << i << ' ' << y << ' ' << x << endl;
  }

  int m = 0;
  if (mode == 0) {
    cin >> m;
  }
  else {
    int yy = (Y[A[i]] + y + L * 10) % L;
    int xx = (X[A[i]] + x + L * 10) % L;
    m = std::max(0, std::min(1000, (int)round(P[yy][xx] + f[keisokuCount])));
    keisokuCount++;
  }

  inputFileKeisokuCost += 100LL * (10LL + abs(y) + abs(x));

  return m;
}

int SA_Keisoku(int i, int y, int x, int se)
{
  int slide = (TansakuSize - 1) / 2;
  int m = 0;

  int yy = (Y[AA[i]] + y + L * 10) % L;
  int xx = (X[AA[i]] + x + L * 10) % L;
  m = std::max(0, std::min(1000, (int)round(P[yy][xx] + f_SA[i][se][y + slide][x + slide])));

  return m;
}

inline int SA_Keisoku_Value(int i, int y, int x, int se, int value)
{
  int slide = (TansakuSize - 1) / 2;
  return std::max(0, std::min(1000, (int)round(value + f_SA[i][se][y + slide][x + slide])));
}

void PrintKeisoku(ofstream& ofs)
{
  int slide = (TansakuSize - 1) / 2;
  rep(i, N)
  {
    int scores[31][31];
    srep(j, -slide, -slide + TansakuSize)
    {
      srep(k, -slide, -slide + TansakuSize)
      {
        scores[j + slide][k + slide] = Keisoku(i, j, k, ofs);
      }
    }

    int minDiff = INF;
    rep(j, N)
    {
      int sumDiff = 0;
      srep(k, -slide, -slide + TansakuSize)
      {
        srep(l, -slide, -slide + TansakuSize)
        {
          int y = (Y[j] + k) % L;
          int x = (X[j] + l) % L;
          sumDiff += abs(P[y][x] - scores[k + slide][l + slide]);
        }
      }

      if (sumDiff < minDiff) {
        minDiff = sumDiff;
        E[i] = j;
      }
    }
  }
}

void PrintKaitou(ofstream& ofs)
{
  if (mode == 0) {
    cout << "-1 -1 -1" << endl;
    fflush(stdout);
    rep(i, N)
    {
      cout << E[i] << endl;
      fflush(stdout);
    }
  }
  else {
    ofs << "-1 -1 -1" << endl;
    rep(i, N) { ofs << E[i] << endl; }
  }
}

ll CalcInputFileScore()
{
  double score = 1e14;
  rep(i, N)
  {
    if (E[i] != A[i]) {
      score *= 0.8;
    }
  }
  score = score / (1e5 + CalcHaitiCost() + inputFileKeisokuCost);
  return ceil(score);
}

void KeepReal()
{
  real_maxScore = maxScore;
  real_maxHaitiCost = maxHaitiCost;
  real_maxKeisokuCost = maxKeisokuCost;
  rep(i, L)
  {
    rep(j, L)
    {
      real_P[i][j] = P[i][j];
    }
  }
  rep(i, SetSize)
  {
    rep(j, N)
    {
      real_EE[i][j] = EE[i][j];
    }
    rep(j, N)
    {
      rep(k, N)
      {
        real_DIFFS[i][j][k] = DIFFS[i][j][k];
      }
    }
  }
}

void RollBackReal()
{
  maxScore = real_maxScore;
  maxHaitiCost = real_maxHaitiCost;
  maxKeisokuCost = real_maxKeisokuCost;
  rep(i, L)
  {
    rep(j, L)
    {
      P[i][j] = real_P[i][j];
    }
  }
  rep(i, SetSize)
  {
    rep(j, N)
    {
      EE[i][j] = real_EE[i][j];
    }
    rep(j, N)
    {
      rep(k, N)
      {
        DIFFS[i][j][k] = real_DIFFS[i][j][k];
      }
    }
  }
}

int scores_SA_Keisoku[SetSize][31][31];
void InitSA_Keisoku()
{
  maxKeisokuCost = 0;

  int slide = (TansakuSize - 1) / 2;
  rep(i, N)
  {
    srep(k, -slide, -slide + TansakuSize)
    {
      srep(l, -slide, -slide + TansakuSize)
      {
        rep(m, SetSize)
        {
          scores_SA_Keisoku[m][k + slide][l + slide] = SA_Keisoku(i, k, l, m);
        }
        maxKeisokuCost += 100LL * (10LL + abs(k) + abs(l));
      }
    }

    rep(m, SetSize)
    {
      int minDiff = INF;
      rep(j, N)
      {
        DIFFS[m][i][j] = 0;
        srep(k, -slide, -slide + TansakuSize)
        {
          srep(l, -slide, -slide + TansakuSize)
          {
            int y = (Y[j] + k + L) % L;
            int x = (X[j] + l + L) % L;
            DIFFS[m][i][j] += abs(P[y][x] - scores_SA_Keisoku[m][k + slide][l + slide]);
          }
        }

        if (DIFFS[m][i][j] < minDiff) {
          minDiff = DIFFS[m][i][j];
          EE[m][i] = j;
        }
      }
    }
  }
}

void InitSA()
{
  // 焼きなまし用データ作成
  rep(i, N)
  {
    AA[i] = i; // 答え
  }
  // ノイズ作成
  InitDist();
  rep(i, N)
  {
    rep(j, SetSize)
    {
      rep(k, TansakuSize)
      {
        rep(l, TansakuSize)
        {
          f_SA[i][j][k][l] = Dist();
        }
      }
    }
  }

  maxHaitiCost = CalcHaitiCost();

  maxKeisokuCost = 0;
  InitSA_Keisoku();

  maxScore = 0;
  rep(m, SetSize)
  {
    double score = 1e14;
    rep(i, N)
    {
      if (EE[m][i] != AA[i]) {
        score *= 0.8;
      }
    }
    score = score / (1e5 + maxHaitiCost + maxKeisokuCost);
    maxScore += ceil(score);
  }

  KeepReal();
}

int IsNear_iy;
int IsNear_ix;
inline bool IsNear(int i, int y, int x)
{
  int slide = (TansakuSize - 1) / 2;

  int iy = Y[i];
  int ix = X[i];
  int iU = Y[i] - slide;
  int iD = iU + TansakuSize;
  int iL = X[i] - slide;
  int iR = iL + TansakuSize;
  srep(j, -1, 2)
  {
    srep(k, -1, 2)
    {
      int yy = y + L * j;
      int xx = x * L * k;
      if (iU <= yy && yy < iD && iL <= xx && xx < iR) {
        IsNear_iy = yy - iU - slide;
        IsNear_ix = xx - iL - slide;
        return true;
      }
    }
  }
  return false;
}

ll CalcDiffHaitiCost(int y, int x, int beforeP, int afterP)
{
  ll ret = 0;
  rep(i, 4)
  {
    int ny = (y + dy[i] + L) % L;
    int nx = (x + dx[i] + L) % L;
    ret -= (beforeP - P[ny][nx]) * (beforeP - P[ny][nx]);
    ret += (afterP - P[ny][nx]) * (afterP - P[ny][nx]);
  }
  return ret;
}

// ランダムな1マスのPを少し変える
int Method1Vector[110];
void Method1(double temperature)
{
  MethodCount[1][0]++;

  int y = randxor() % L;
  int x = randxor() % L;
  int diff = rand() % 201 - 100;
  int beforeP = P[y][x];
  int afterP = max(0, min(1000, P[y][x] + diff));

  int cnt = 0;

  // 探索範囲に入っているワームホールの差分更新
  rep(i, N)
  {
    if (!IsNear(i, y, x)) {
      continue;
    }
    int iy = IsNear_iy;
    int ix = IsNear_ix;


    rep(m, SetSize)
    {
      KEEP_EE[m][i] = EE[m][i];
      rep(j, N)
      {
        KEEP_DIFFS[m][i][j] = DIFFS[m][i][j];
      }
    }

    rep(m, SetSize)
    {
      int beforeValue = SA_Keisoku_Value(i, iy, ix, m, beforeP);
      rep(j, N)
      {
        int beforeDiff = abs(P[(Y[j] + iy + L) % L][(X[j] + ix + L) % L] - beforeValue);
        DIFFS[m][i][j] -= beforeDiff;
      }
    }

    P[y][x] = afterP;

    rep(m, SetSize)
    {
      int afterValue = SA_Keisoku_Value(i, iy, ix, m, afterP);
      rep(j, N)
      {
        int afterDiff = abs(P[(Y[j] + iy + L) % L][(X[j] + ix + L) % L] - afterValue);
        DIFFS[m][i][j] += afterDiff;
      }
    }

    P[y][x] = beforeP;

    // EE[m][i]を更新
    rep(m, SetSize)
    {
      int minDiff = INF;
      rep(j, N)
      {
        if (DIFFS[m][i][j] < minDiff) {
          minDiff = DIFFS[m][i][j];
          EE[m][i] = j;
        }
      }
    }

    Method1Vector[cnt] = i;
    cnt++;
  }

  if (cnt == 0) {
    return;
  }
  MethodCount[1][1]++;

  P[y][x] = afterP;

  // スコア更新
  ll haitiDiff = CalcDiffHaitiCost(y, x, beforeP, afterP);
  ll keisokudiff = 0;
  ll tmpSumScore = 0;
  rep(m, SetSize)
  {
    double dTmpScore = 1e14;
    rep(i, N)
    {
      if (EE[m][i] != AA[i]) {
        dTmpScore *= 0.8;
      }
    }
    dTmpScore = dTmpScore / (1e5 + maxHaitiCost + haitiDiff + maxKeisokuCost + keisokudiff);
    ll tmpScore = ceil(dTmpScore);
    tmpSumScore += tmpScore;
  }

  ll diffScore = tmpSumScore - maxScore;
  double prob = exp((double)diffScore / temperature);
  if (prob > rand01()) {
    MethodCount[1][2]++;
    maxScore += diffScore;
    maxHaitiCost += haitiDiff;
    maxKeisokuCost += keisokudiff;

    if (maxScore > real_maxScore) {
      MethodCount[1][3]++;
      if (MethodCount[1][3] % 100 == 0) {
        cout << maxScore << ' ' << real_maxScore << endl;
        KeepReal();
      }
    }
  }
  else {
    // 元に戻す
    P[y][x] = beforeP;
    rep(ii, cnt)
    {
      int i = Method1Vector[ii];
      rep(m, SetSize)
      {
        EE[m][i] = KEEP_EE[m][i];
        rep(j, N)
        {
          DIFFS[m][i][j] = KEEP_DIFFS[m][i][j];
        }
      }
    }
  }
}

ll Solve(int probNum)
{
  Init();

  ofstream ofs;
  Input(probNum);
  OpenOfs(probNum, ofs);
  //InitHaiti();
  InitHaiti2();
  //InitHaiti3();
  //InitHaiti4();
  //InitHaiti5();

  // 焼きなまし
  InitSA();
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();
  double startTemperature = 2000;
  double endTemperature = 0;
  int loopCount = 0;
  double nowProgress = 0;
  while (true) {
    break;
    if (loopCount % 100 == 0) {
      endTime = clock();
      double nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
      if (nowTime > TL) {
        break;
      }
      nowProgress = nowTime / TL;
    }

    loopCount++;

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;


    Method1(temperature);
  }

  // 戻す
  RollBackReal();

  // 配置を回答
  PrintHaiti(ofs);

  // 計測
  PrintKeisoku(ofs);

  // 解答を回答
  PrintKaitou(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcInputFileScore();
  }
  return score;
}


int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }
  std::random_device rnd;
  engine.seed(rnd());

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      cout << "num = " << i << ", ";
      cout << "score = " << score << endl;
      //cout << "Method1 : ";
      //rep(j, 4) {
      //  cout << MethodCount[1][j] << ' ';
      //}
      //cout << endl;
      //rep(i, N) {
      //  rep(j, SetSize) {
      //    cout << EE[j][i] << ' ';
      //  }
      //  cout << endl;
      //}
    }
  }

  return 0;
}
