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
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

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

const ll INF      = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 2.8;
int mode;
clock_t startTime, endTime;

double GetNowTime()
{
  endTime        = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int MAX_D         = 55;
const int MAX_N         = 56;
const int MAX_LINECOUNT = 60;
int lineMaxLimit        = 100;

const int w = 1000;
int d, n;
double ee;
int a[MAX_D][MAX_N];
int maxA[MAX_N];
int sumA[MAX_D];
int daysDifficultySorted[MAX_D];

int mostVariableAsCount[MAX_D][w + 10];
int mostVariableAsValue[MAX_D][w + 10][10];
int mostVariableAs[MAX_D][w + 10][10][2];

int ans[MAX_D][MAX_N][4];
ll ansScore;
int ansLinePos[MAX_D][MAX_LINECOUNT] = {};
int ansLineCount[MAX_D];
int ansBaseLineCount;

int temp_ans[MAX_D][MAX_N][4];
ll temp_ansScore;
int temp_ansLinePos[MAX_D][MAX_LINECOUNT] = {};
int temp_ansLineCount[MAX_D];
int temp_ansBaseLineCount;

int real_ans[MAX_D][MAX_N][4];
ll real_ansScore;
int real_ansLinePos[MAX_D][MAX_LINECOUNT] = {};
int real_ansLineCount[MAX_D];
int real_ansBaseLineCount;

int real_real_ans[MAX_D][MAX_N][4];
ll real_real_ansScore;
int real_real_ansLinePos[MAX_D][MAX_LINECOUNT] = {};
int real_real_ansLineCount[MAX_D];
int real_real_ansBaseLineCount;

int keep31Count    = 0;
int keep31KeepSize = 2;
int keep31_ans[10][MAX_D][MAX_N][4];
ll keep31_ansScore[10];
int keep31_ansLinePos[10][MAX_D][MAX_LINECOUNT] = {};
int keep31_ansLineCount[10][MAX_D];
int keep31_ansBaseLineCount[10];

void CopyToRealRealAns()
{
  real_real_ansScore = real_ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        real_real_ans[i][j][k] = real_ans[i][j][k];
      }
    }
  }
  real_real_ansBaseLineCount = real_ansBaseLineCount;
  rep(i, d)
  {
    real_real_ansLineCount[i] = real_ansLineCount[i];
    rep(j, real_ansLineCount[i] + 1)
    {
      real_real_ansLinePos[i][j] = real_ansLinePos[i][j];
    }
  }
}

void CopyToRealAns()
{
  real_ansScore = ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        real_ans[i][j][k] = ans[i][j][k];
      }
    }
  }
  real_ansBaseLineCount = ansBaseLineCount;
  rep(i, d)
  {
    real_ansLineCount[i] = ansLineCount[i];
    rep(j, ansLineCount[i] + 1)
    {
      real_ansLinePos[i][j] = ansLinePos[i][j];
    }
  }
}

void CopyFromRealRealAns()
{
  real_ansScore = real_real_ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        real_ans[i][j][k] = real_real_ans[i][j][k];
      }
    }
  }
  real_ansBaseLineCount = real_real_ansBaseLineCount;
  rep(i, d)
  {
    real_ansLineCount[i] = real_real_ansLineCount[i];
    rep(j, real_ansLineCount[i] + 1)
    {
      real_ansLinePos[i][j] = real_real_ansLinePos[i][j];
    }
  }
}

void CopyFromRealAns()
{
  ansScore = real_ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        ans[i][j][k] = real_ans[i][j][k];
      }
    }
  }
  ansBaseLineCount = real_ansBaseLineCount;
  rep(i, d)
  {
    ansLineCount[i] = real_ansLineCount[i];
    rep(j, ansLineCount[i] + 1)
    {
      ansLinePos[i][j] = real_ansLinePos[i][j];
    }
  }
}

void CopyToKeep31(int idx)
{
  keep31_ansScore[idx] = ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        keep31_ans[idx][i][j][k] = ans[i][j][k];
      }
    }
  }
  keep31_ansBaseLineCount[idx] = ansBaseLineCount;
  rep(i, d)
  {
    keep31_ansLineCount[idx][i] = ansLineCount[i];
    rep(j, ansLineCount[i] + 1)
    {
      keep31_ansLinePos[idx][i][j] = ansLinePos[i][j];
    }
  }
}

void CopyFromKeep31(int idx)
{
  ansScore = keep31_ansScore[idx];
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        ans[i][j][k] = keep31_ans[idx][i][j][k];
      }
    }
  }
  ansBaseLineCount = keep31_ansBaseLineCount[idx];
  rep(i, d)
  {
    ansLineCount[i] = keep31_ansLineCount[idx][i];
    rep(j, ansLineCount[i] + 1)
    {
      ansLinePos[i][j] = keep31_ansLinePos[idx][i][j];
    }
  }
}

void CopyToTemp()
{
  temp_ansScore = ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        temp_ans[i][j][k] = ans[i][j][k];
      }
    }
  }
  temp_ansBaseLineCount = ansBaseLineCount;
  rep(i, d)
  {
    temp_ansLineCount[i] = ansLineCount[i];
    rep(j, ansLineCount[i] + 1)
    {
      temp_ansLinePos[i][j] = ansLinePos[i][j];
    }
  }
}

void CopyFromTemp()
{
  ansScore = temp_ansScore;
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, 4)
      {
        ans[i][j][k] = temp_ans[i][j][k];
      }
    }
  }
  ansBaseLineCount = temp_ansBaseLineCount;
  rep(i, d)
  {
    ansLineCount[i] = temp_ansLineCount[i];
    rep(j, ansLineCount[i] + 1)
    {
      ansLinePos[i][j] = temp_ansLinePos[i][j];
    }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ansScore      = INF;
  real_ansScore = INF;
}

void InitMostVariableAs()
{
  rep(i, d)
  {
    srep(j, 1, w + 1)
    {
      mostVariableAsCount[i][j] = 0;
      rep(k, 10)
      {
        mostVariableAsValue[i][j][k] = 0;
      }
      int n2     = 0;
      int n2Need = (a[i][n2] - 1) / j + 1;
      if (n2Need > w) continue;
      drep(n1, n)
      {
        if (n2 == n) break;
        int val    = a[i][n1];
        int n1Need = (val - 1) / j + 1;
        if (n1Need + n2Need > w) continue;
        while (n2 < n - 1) {
          int nxtNeed = (a[i][n2 + 1] - 1) / j + 1;
          if (n1Need + nxtNeed <= w) {
            n2++;
            n2Need = nxtNeed;
          }
          else {
            break;
          }
        }
        rep(k, 10)
        {
          if (n2 - k < 0) break;
          if (n2 - k == n1) continue;
          int valSum = val + a[i][n2 - k];
          int junni  = mostVariableAsCount[i][j];
          if (junni == 10) {
            if (valSum > mostVariableAsValue[i][j][9]) {
              mostVariableAsValue[i][j][9] = valSum;
              mostVariableAs[i][j][9][0]   = n1;
              mostVariableAs[i][j][9][1]   = n2 - k;
              junni                        = 9;
            }
          }
          if (junni == 10) break;
          while (junni > 0) {
            if (mostVariableAsValue[i][j][junni] > mostVariableAsValue[i][j][junni - 1]) {
              int swa                              = mostVariableAsValue[i][j][junni - 1];
              mostVariableAsValue[i][j][junni - 1] = mostVariableAsValue[i][j][junni];
              mostVariableAsValue[i][j][junni]     = swa;
              swa                                  = mostVariableAs[i][j][junni - 1][0];
              mostVariableAs[i][j][junni - 1][0]   = mostVariableAs[i][j][junni][0];
              mostVariableAs[i][j][junni][0]       = swa;
              swa                                  = mostVariableAs[i][j][junni - 1][1];
              mostVariableAs[i][j][junni - 1][1]   = mostVariableAs[i][j][junni][1];
              mostVariableAs[i][j][junni][1]       = swa;
              junni--;
            }
            else {
              break;
            }
          }
          if (mostVariableAsCount[i][j] < 10) mostVariableAsCount[i][j]++;
        }
      }
    }
  }
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
    int www;
    cin >> www >> d >> n;
    rep(i, d)
    {
      rep(j, n)
      {
        cin >> a[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    int www;
    ifs >> www >> d >> n;
    rep(i, d)
    {
      rep(j, n)
      {
        ifs >> a[i][j];
      }
    }
  }

  rep(j, n)
  {
    maxA[j] = 0;
  }
  rep(i, d)
  {
    rep(j, n)
    {
      maxA[j] = max(maxA[j], a[i][j]);
    }
  }
  rep(i, d)
  {
    sumA[i] = 0;
  }
  rep(i, d)
  {
    rep(j, n)
    {
      sumA[i] += a[i][j];
    }
  }
  vector<P> sumsPairs;
  rep(i, d)
  {
    sumsPairs.emplace_back(sumA[i], i);
  }
  sort(sumsPairs.begin(), sumsPairs.end());
  rep(i, d)
  {
    daysDifficultySorted[i] = sumsPairs[i].second;
  }

  int allSum = 0;
  rep(i, d)
  {
    allSum += sumA[i];
  }
  ee = (double)(w * w * d - allSum) / (w * w * d);

  // InitMostVariableAs();
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

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    rep(i, d)
    {
      rep(j, n)
      {
        rep(k, 4)
        {
          cout << ans[i][j][k] << ' ';
        }
        cout << endl;
      }
    }
  }
  else {
    rep(i, d)
    {
      rep(j, n)
      {
        rep(k, 4)
        {
          ofs << ans[i][j][k] << ' ';
        }
        ofs << endl;
      }
    }
  }
}

bool IsNGAns()
{
  rep(i, d)
  {
    rep(j, n)
    {
      if (ans[i][j][0] < 0) return true;
      if (ans[i][j][1] < 0) return true;
      if (ans[i][j][0] >= ans[i][j][2]) return true;
      if (ans[i][j][1] >= ans[i][j][3]) return true;
      if (ans[i][j][2] > w) return true;
      if (ans[i][j][3] > w) return true;
      if ((ans[i][j][2] - ans[i][j][0]) * (ans[i][j][3] - ans[i][j][1]) < a[i][j]) return true;
    }
  }
  rep(i, d)
  {
    rep(j, n)
    {
      srep(k, j + 1, n)
      {
        if (max(ans[i][j][0], ans[i][k][0]) < min(ans[i][j][2], ans[i][k][2])) {
          if (max(ans[i][j][1], ans[i][k][1]) < min(ans[i][j][3], ans[i][k][3])) { return true; }
        }
      }
    }
  }
  return false;
}

// スコア計算
ll CalcScore()
{
  ll res = 1;
  rep(i, d)
  {
    rep(j, n)
    {
      int sx = ans[i][j][0];
      int sy = ans[i][j][1];
      int tx = ans[i][j][2];
      int ty = ans[i][j][3];
      res += (ll)max(0, a[i][j] - (tx - sx) * (ty - sy)) * 100;
    }
  }
  {
    vector<P> before;
    rep(i, d)
    {
      vector<P> now;
      rep(j, n)
      {
        int sx = ans[i][j][0];
        int sy = ans[i][j][1];
        int tx = ans[i][j][2];
        int ty = ans[i][j][3];
        if (sy != 0 && sy != w) {
          now.emplace_back(sy * w + sx, 1);
          now.emplace_back(sy * w + tx, -1);
        }
        if (ty != 0 && ty != w) {
          now.emplace_back(ty * w + sx, 1);
          now.emplace_back(ty * w + tx, -1);
        }
      }
      sort(now.begin(), now.end());
      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < before.size() || itr2 < now.size()) {
          if (itr1 == before.size()) {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
          else if (itr2 == now.size()) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else if (before[itr1].first <= now[itr2].first) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
        }
      }
      before = now;
    }
  }
  {
    vector<P> before;
    rep(i, d)
    {
      vector<P> now;
      rep(j, n)
      {
        int sx = ans[i][j][0];
        int sy = ans[i][j][1];
        int tx = ans[i][j][2];
        int ty = ans[i][j][3];
        if (sx != 0 && sx != w) {
          now.emplace_back(sx * w + sy, 1);
          now.emplace_back(sx * w + ty, -1);
        }
        if (tx != 0 && tx != w) {
          now.emplace_back(tx * w + sy, 1);
          now.emplace_back(tx * w + ty, -1);
        }
      }
      sort(now.begin(), now.end());
      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < before.size() || itr2 < now.size()) {
          if (itr1 == before.size()) {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
          else if (itr2 == now.size()) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else if (before[itr1].first <= now[itr2].first) {
            int val = before[itr1].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += before[itr1].second;
            itr1++;
          }
          else {
            int val = now[itr2].first;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += now[itr2].second;
            itr2++;
          }
        }
      }
      before = now;
    }
  }
  return res;
}

int CalcScoreForMethod3BeforeArr[1000];
int CalcScoreForMethod3NowArr[1000];
int CalcScoreForMethod3Used[1000];
ll CalcScoreForMethod3()
{
  ll res = 1;
  rep(i, d)
  {
    rep(j, n)
    {
      int sx = ans[i][j][0];
      int sy = ans[i][j][1];
      int tx = ans[i][j][2];
      int ty = ans[i][j][3];
      res += (ll)max(0, a[i][j] - (tx - sx) * (ty - sy)) * 100;
    }
  }

  // 縦線
  {rep(i, d)
  {
    rep(j, ansLineCount[i] + 1)
    {
      if (j == 0 || j == ansLineCount[i]) {
        CalcScoreForMethod3Used[j] = 1;
      }
      else {
        CalcScoreForMethod3Used[j] = 0;
      }
    }
    rep(j, n)
    {
      rep(k, ansLineCount[i])
      {
        if (ans[i][j][1] == ansLinePos[i][k]) {
          CalcScoreForMethod3Used[k]     = 1;
          CalcScoreForMethod3Used[k + 1] = 1;
        }
      }
    }
    rep(j, ansLineCount[i] + 1)
    {
      if (CalcScoreForMethod3Used[j] == 0) { res += w * 2; }
    }
  }
  }

  // 横線
  {
    int beforeTail = 0;
    rep(i, d)
    {
      int nowTail = 0;
      rep(j, n)
      {
        int sx = ans[i][j][0];
        int sy = ans[i][j][1];
        int tx = ans[i][j][2];
        int ty = ans[i][j][3];
        if (sx != 0 && sx != w) {
          CalcScoreForMethod3NowArr[nowTail] = (sx * w + sy) * 10 + 1;
          nowTail++;
          CalcScoreForMethod3NowArr[nowTail] = (sx * w + ty) * 10 - 1;
          nowTail++;
        }
        if (tx != 0 && tx != w) {
          CalcScoreForMethod3NowArr[nowTail] = (tx * w + sy) * 10 + 1;
          nowTail++;
          CalcScoreForMethod3NowArr[nowTail] = (tx * w + ty) * 10 - 1;
          nowTail++;
        }
      }

      sort(CalcScoreForMethod3NowArr, CalcScoreForMethod3NowArr + nowTail);

      if (i > 0) {
        int itr1 = 0, itr2 = 0;
        int cnt1 = 0, cnt2 = 0;
        int pos = 0;
        while (itr1 < beforeTail || itr2 < nowTail) {
          if (itr1 == beforeTail) {
            int val = (CalcScoreForMethod3NowArr[itr2] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += CalcScoreForMethod3NowArr[itr2] - val * 10;
            itr2++;
          }
          else if (itr2 == nowTail) {
            int val = (CalcScoreForMethod3BeforeArr[itr1] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += CalcScoreForMethod3BeforeArr[itr1] - val * 10;
            itr1++;
          }
          else if ((CalcScoreForMethod3BeforeArr[itr1] + 1) / 10 <= (CalcScoreForMethod3NowArr[itr2] + 1) / 10) {
            int val = (CalcScoreForMethod3BeforeArr[itr1] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt1 += CalcScoreForMethod3BeforeArr[itr1] - val * 10;
            itr1++;
          }
          else {
            int val = (CalcScoreForMethod3NowArr[itr2] + 1) / 10;
            if ((cnt1 == 0 && cnt2 > 0) || (cnt1 > 0 && cnt2 == 0)) { res += (ll)val - pos; }
            pos = val;
            cnt2 += CalcScoreForMethod3NowArr[itr2] - val * 10;
            itr2++;
          }
        }
      }
      beforeTail = nowTail;
      rep(j, nowTail)
      {
        CalcScoreForMethod3BeforeArr[j] = CalcScoreForMethod3NowArr[j];
      }
    }
  }
  return res;
}

void Initialize()
{
  rep(i, d)
  {
    int now = 0;
    rep(j, n)
    {
      int need     = (a[i][j] - 1) / w + 1;
      int newY     = min(now + need, w - (n - 1 - j));
      ans[i][j][0] = 0;
      ans[i][j][1] = now;
      ans[i][j][2] = w;
      ans[i][j][3] = newY;
      now          = newY;
    }
  }

  ansScore = CalcScore();
  CopyToRealAns();
}

void Method1()
{
  int ok = 1;
  rep(i, d)
  {
    int x           = 0;
    int y           = 0;
    int used[MAX_N] = {};
    rep(j, n)
    {
      int minAmari = INT_INF;
      int tmpansjj = -1;
      int tmpans2  = -1;
      int tmpans3  = -1;
      int newx     = -1;
      int newy     = -1;
      rep(jj, n)
      {
        if (used[jj]) continue;
        int lenx = w - x;
        int leny = w - y;
        if (a[i][jj] > lenx * leny) continue;
        int need1  = (a[i][jj] - 1) / leny + 1;
        int need2  = (a[i][jj] - 1) / lenx + 1;
        int amari1 = need1 * leny - a[i][jj];
        int amari2 = need2 * lenx - a[i][jj];
        if (amari1 <= minAmari) {
          minAmari = amari1;
          tmpansjj = jj;
          tmpans2  = x + need1;
          tmpans3  = w;
          newx     = x + need1;
          newy     = y;
        }
        if (amari2 <= minAmari) {
          minAmari = amari2;
          tmpansjj = jj;
          tmpans2  = w;
          tmpans3  = y + need2;
          newx     = x;
          newy     = y + need2;
        }
      }
      if (tmpansjj == -1) {
        ok = 0;
        break;
      }
      ans[i][tmpansjj][0] = x;
      ans[i][tmpansjj][1] = y;
      ans[i][tmpansjj][2] = tmpans2;
      ans[i][tmpansjj][3] = tmpans3;
      used[tmpansjj]      = 1;
      x                   = newx;
      y                   = newy;
    }
    if (ok == 0) { break; }
  }
  if (ok) {
    ansScore = CalcScore();
    if (ansScore < real_ansScore) { CopyToRealAns(); }
  }
}

void MethodPerfect()
{
  ll maxASum = 0;
  rep(i, n)
  {
    maxASum += maxA[i];
  }
  if (maxASum > w * w) { return; }

  int ok = 1;
  rep(i, d)
  {
    int x           = 0;
    int y           = 0;
    int used[MAX_N] = {};
    rep(j, n)
    {
      int minAmari = INT_INF;
      int tmpansjj = -1;
      int tmpans2  = -1;
      int tmpans3  = -1;
      int newx     = -1;
      int newy     = -1;
      rep(jj, n)
      {
        if (used[jj]) continue;
        int lenx = w - x;
        int leny = w - y;
        if (maxA[jj] > lenx * leny) {
          ok = 0;
          break;
        }
        int need1  = (maxA[jj] - 1) / leny + 1;
        int need2  = (maxA[jj] - 1) / lenx + 1;
        int amari1 = need1 * leny - maxA[jj];
        int amari2 = need2 * lenx - maxA[jj];
        if (amari1 <= minAmari) {
          minAmari = amari1;
          tmpansjj = jj;
          tmpans2  = x + need1;
          tmpans3  = w;
          newx     = x + need1;
          newy     = y;
        }
        if (amari2 <= minAmari) {
          minAmari = amari2;
          tmpansjj = jj;
          tmpans2  = w;
          tmpans3  = y + need2;
          newx     = x;
          newy     = y + need2;
        }
      }
      if (ok == 0) { break; }
      if (tmpansjj >= 0) {
        ans[i][tmpansjj][0] = x;
        ans[i][tmpansjj][1] = y;
        ans[i][tmpansjj][2] = tmpans2;
        ans[i][tmpansjj][3] = tmpans3;
        used[tmpansjj]      = 1;
        x                   = newx;
        y                   = newy;
      }
    }
    if (ok == 0) { break; }
  }

  if (ok) {
    ansScore = CalcScore();
    CopyToRealAns();
  }
}

int preCalcScheduleSizes[MAX_D][MAX_N][MAX_LINECOUNT];
int widths[MAX_LINECOUNT] = {};

int ansColumnNum[MAX_D][MAX_N];
int ansColumnSchedules[MAX_D][MAX_LINECOUNT][MAX_N];
int ansColumnSchedulesCount[MAX_D][MAX_LINECOUNT];
int ansColumnSchedulesPosition[MAX_D][MAX_LINECOUNT][MAX_N];

int real_ansColumnNum[MAX_D][MAX_N];
int real_ansColumnSchedules[MAX_D][MAX_LINECOUNT][MAX_N];
int real_ansColumnSchedulesCount[MAX_D][MAX_LINECOUNT];
int real_ansColumnSchedulesPosition[MAX_D][MAX_LINECOUNT][MAX_N];

// realに格納
void CopyToReal_M42()
{
  rep(i, d)
  {
    rep(j, n)
    {
      real_ansColumnNum[i][j] = ansColumnNum[i][j];
    }
    rep(j, ansBaseLineCount)
    {
      real_ansColumnSchedulesCount[i][j] = ansColumnSchedulesCount[i][j];
      rep(k, real_ansColumnSchedulesCount[i][j] + 1)
      {
        real_ansColumnSchedules[i][j][k]         = ansColumnSchedules[i][j][k];
        real_ansColumnSchedulesPosition[i][j][k] = ansColumnSchedulesPosition[i][j][k];
      }
    }
  }
}

// realから戻す
void CoptToCurrent_M42()
{
  rep(i, d)
  {
    rep(j, n)
    {
      ansColumnNum[i][j] = real_ansColumnNum[i][j];
    }
    rep(j, ansBaseLineCount)
    {
      ansColumnSchedulesCount[i][j] = real_ansColumnSchedulesCount[i][j];
      rep(k, ansColumnSchedulesCount[i][j] + 1)
      {
        ansColumnSchedules[i][j][k]         = real_ansColumnSchedules[i][j][k];
        ansColumnSchedulesPosition[i][j][k] = real_ansColumnSchedulesPosition[i][j][k];
      }
    }
  }
}

int CalcDiffScore2(int day1, int day2, int lineNum)
{
  int diffScore2 = 0;
  int ite1 = 1, ite2 = 1;
  while (ite1 < ansColumnSchedulesCount[day1][lineNum] && ite2 < ansColumnSchedulesCount[day2][lineNum]) {
    int num1 = ansColumnSchedulesPosition[day1][lineNum][ite1];
    int num2 = ansColumnSchedulesPosition[day2][lineNum][ite2];
    if (num1 == num2) {
      diffScore2 -= widths[lineNum] * 2;
      ite1++;
      ite2++;
    }
    else if (num1 < num2) {
      ite1++;
    }
    else {
      ite2++;
    }
  }
  return diffScore2;
}

int CD3_Members[MAX_N];
int CD3_TmpPosition[MAX_N];
int CD3_NeighborPos[MAX_N * 2];
int shuffleThree[6][3] = { {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0} };
int CalcDiffScore3(int winter, int memberCount, int neighborPosCount, int day, int lineNum, int margin)
{
  int diffScore3 = 0;
  if (memberCount == 0) { return diffScore3; }
  if (winter >= 1) {
    if (memberCount == 1) {
      ;
    }
    else if (memberCount == 2) {
      int swa        = CD3_Members[0];
      CD3_Members[0] = CD3_Members[1];
      CD3_Members[1] = swa;
    }
    else if (memberCount == 3) {
      int tmp[3] = {};
      rep(j, 3)
      {
        tmp[shuffleThree[winter - 1][j]] = CD3_Members[j];
      }
      rep(j, 3)
      {
        CD3_Members[j] = tmp[shuffleThree[winter][j]];
      }
    }
    else {
      FisherYates(CD3_Members, memberCount);
    }
  }
  int tmpMargin         = margin;
  int ite1              = 0;
  int sum               = preCalcScheduleSizes[day][CD3_Members[0]][lineNum];
  CD3_TmpPosition[0]    = 0;
  int ite2              = 1;
  CD3_TmpPosition[ite2] = sum;
  while (ite1 < neighborPosCount && ite2 < memberCount) {
    int pos1 = CD3_NeighborPos[ite1];
    if (pos1 < sum) {
      ite1++;
    }
    else if (sum + tmpMargin < pos1) {
      sum += preCalcScheduleSizes[day][CD3_Members[ite2]][lineNum];
      ite2++;
      CD3_TmpPosition[ite2] = sum;
    }
    else {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      if (ite1 < neighborPosCount) {
        int nxtPos1 = CD3_NeighborPos[ite1];
        if (nxtPos1 == pos1) {
          diffScore3 += widths[lineNum] * 2;
          ite1++;
        }
      }
      tmpMargin -= pos1 - sum;
      sum                   = pos1;
      CD3_TmpPosition[ite2] = sum;
      sum += preCalcScheduleSizes[day][CD3_Members[ite2]][lineNum];
      ite2++;
      CD3_TmpPosition[ite2] = sum;
    }
  }
  while (ite2 < memberCount) {
    sum += preCalcScheduleSizes[day][CD3_Members[ite2]][lineNum];
    ite2++;
    CD3_TmpPosition[ite2] = sum;
  }
  CD3_TmpPosition[ite2] = w;
  return diffScore3;
}

int CD32_Members[MAX_N];
int CD32_TmpPosition[MAX_N];
int CD32_NeighborNewPos1[MAX_N];
int CD32_NeighborNewPos2[MAX_N];
int CalcDiffScore3_2(int winter, int memberCount, int day, int lineNum, int margin)
{
  int beforeCount = 0;
  if (day > 0) {
    beforeCount = ansColumnSchedulesCount[day - 1][lineNum];
    rep(k, beforeCount + 1)
    {
      CD32_NeighborNewPos1[k] = ansColumnSchedulesPosition[day - 1][lineNum][k];
    }
  }
  int afterCount = 0;
  if (day < d - 1) {
    afterCount = ansColumnSchedulesCount[day + 1][lineNum];
    rep(k, afterCount + 1)
    {
      CD32_NeighborNewPos2[k] = ansColumnSchedulesPosition[day + 1][lineNum][k];
    }
  }

  int diffScore3 = 0;

  if (memberCount == 0) { return diffScore3; }
  if (winter >= 1) {
    if (memberCount == 1) {
      ;
    }
    else if (memberCount == 2) {
      int swa         = CD32_Members[0];
      CD32_Members[0] = CD32_Members[1];
      CD32_Members[1] = swa;
    }
    else if (memberCount == 3) {
      int tmp[3] = {};
      rep(j, 3)
      {
        tmp[shuffleThree[winter - 1][j]] = CD32_Members[j];
      }
      rep(j, 3)
      {
        CD32_Members[j] = tmp[shuffleThree[winter][j]];
      }
    }
    else {
      FisherYates(CD32_Members, memberCount);
    }
  }
  int tmpMargin          = margin;
  int ite1               = 1;
  int sum                = preCalcScheduleSizes[day][CD32_Members[0]][lineNum];
  CD32_TmpPosition[0]    = 0;
  CD32_TmpPosition[ite1] = sum;

  int ite2 = 1;
  int ite3 = 1;

  while (ite1 < memberCount) {
    int num2 = 0, pos2 = 0, pos22 = 0;
    if (ite2 < beforeCount) {
      num2  = ansColumnSchedules[day - 1][lineNum][ite2 - 1];
      pos2  = ansColumnSchedulesPosition[day - 1][lineNum][ite2 - 1] + preCalcScheduleSizes[day - 1][num2][lineNum];
      pos22 = ansColumnSchedulesPosition[day - 1][lineNum][ite2];
      pos2  = min(pos22, max(pos2, sum + 1));
      pos2  = pos22;
      if (1 < day) {
        rep(k, ansColumnSchedulesCount[day - 2][lineNum])
        {
          if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos22) {
            pos2 = pos22;
            break;
          }
        }
      }
    }
    int num3 = 0, pos3 = 0, pos32 = 0;
    if (ite3 < afterCount) {
      num3  = ansColumnSchedules[day + 1][lineNum][ite3 - 1];
      pos3  = ansColumnSchedulesPosition[day + 1][lineNum][ite3 - 1] + preCalcScheduleSizes[day + 1][num3][lineNum];
      pos32 = ansColumnSchedulesPosition[day + 1][lineNum][ite3];
      pos3  = min(pos32, max(pos3, sum + 1));
      pos3  = pos32;
      if (day < n - 2) {
        rep(k, ansColumnSchedulesCount[day + 2][lineNum])
        {
          if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos32) {
            pos3 = pos32;
            break;
          }
        }
      }
    }
    if (ite2 < beforeCount && ite3 < afterCount) {
      if (pos22 < sum) {
        ite2++;
        continue;
      }
      if (pos32 < sum) {
        ite3++;
        continue;
      }
      if (sum + tmpMargin < pos2 && sum + tmpMargin < pos3) {
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
        continue;
      }
      if (pos2 <= sum + tmpMargin && pos3 <= sum + tmpMargin) {
        if (pos22 == pos32) {
          bool isConnect1 = false;
          if (1 < day) {
            rep(k, ansColumnSchedulesCount[day - 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos22) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect1 = true;
                break;
              }
            }
          }
          bool isConnect2 = false;
          if (day < n - 2) {
            rep(k, ansColumnSchedulesCount[day + 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos32) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect2 = true;
                break;
              }
            }
          }
          if (pos22 <= sum + tmpMargin && pos32 <= sum + tmpMargin && (isConnect1 || isConnect2)) {
            // nPos = pos22 = pos32
            if (isConnect1) diffScore3 += widths[lineNum] * 2;
            if (isConnect2) diffScore3 += widths[lineNum] * 2;
            tmpMargin -= pos22 - sum;
            sum                    = pos22;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
            ite2++;
            ite3++;
          }
          else {
            // nPos = min(pos2,pos3)
            if (pos2 <= pos3) {
              // nPos = pos2
              CD32_NeighborNewPos1[ite2] = pos2;
              if (1 < day) {
                rep(k, ansColumnSchedulesCount[day - 2][lineNum])
                {
                  if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos2) {
                    diffScore3 += widths[lineNum] * 2;
                    break;
                  }
                }
              }
              tmpMargin -= pos2 - sum;
              sum                    = pos2;
              CD32_TmpPosition[ite1] = sum;
              sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
              ite1++;
              CD32_TmpPosition[ite1] = sum;
              ite2++;
            }
            else {
              // nPos = pos3
              CD32_NeighborNewPos2[ite3] = pos3;
              if (day < n - 2) {
                rep(k, ansColumnSchedulesCount[day + 2][lineNum])
                {
                  if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos3) {
                    diffScore3 += widths[lineNum] * 2;
                    break;
                  }
                }
              }
              tmpMargin -= pos3 - sum;
              sum                    = pos3;
              CD32_TmpPosition[ite1] = sum;
              sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
              ite1++;
              CD32_TmpPosition[ite1] = sum;
              ite3++;
            }
          }
        }
        else if (pos22 < pos32) {
          bool isConnect = false;
          if (1 < day) {
            rep(k, ansColumnSchedulesCount[day - 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos22) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect = true;
                break;
              }
            }
          }
          if (pos22 <= sum + tmpMargin && isConnect) {
            // nPos = pos22
            if (isConnect) diffScore3 += widths[lineNum] * 2;
            tmpMargin -= pos22 - sum;
            sum                    = pos22;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          else {
            // nPos = pos2
            CD32_NeighborNewPos1[ite2] = pos2;
            if (1 < day) {
              rep(k, ansColumnSchedulesCount[day - 2][lineNum])
              {
                if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos2) {
                  diffScore3 += widths[lineNum] * 2;
                  break;
                }
              }
            }
            tmpMargin -= pos2 - sum;
            sum                    = pos2;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          ite2++;
        }
        else {
          bool isConnect = false;
          if (day < n - 2) {
            rep(k, ansColumnSchedulesCount[day + 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos32) {
                diffScore3 -= widths[lineNum] * 2;
                isConnect = true;
                break;
              }
            }
          }
          if (pos32 <= sum + tmpMargin && isConnect) {
            // nPos = pos32
            if (isConnect) diffScore3 += widths[lineNum] * 2;
            tmpMargin -= pos32 - sum;
            sum                    = pos32;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          else {
            // nPos = pos3
            CD32_NeighborNewPos2[ite3] = pos3;
            if (day < n - 2) {
              rep(k, ansColumnSchedulesCount[day + 2][lineNum])
              {
                if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos3) {
                  diffScore3 += widths[lineNum] * 2;
                  break;
                }
              }
            }
            tmpMargin -= pos3 - sum;
            sum                    = pos3;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          ite3++;
        }
      }
      else if (pos2 <= sum + tmpMargin) {
        if (pos22 < sum) {
          ite2++;
          continue;
        }
        if (sum + tmpMargin < pos2) {
          sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
          ite1++;
          CD32_TmpPosition[ite1] = sum;
          continue;
        }
        bool isConnect = false;
        if (1 < day) {
          rep(k, ansColumnSchedulesCount[day - 2][lineNum])
          {
            if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos22) {
              diffScore3 -= widths[lineNum] * 2;
              isConnect = true;
              break;
            }
          }
        }
        if (isConnect) {
          if (pos22 <= sum + tmpMargin) {
            diffScore3 += widths[lineNum] * 2;
            tmpMargin -= pos22 - sum;
            sum                    = pos22;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          else {
            CD32_NeighborNewPos1[ite2] = pos2;
            if (1 < day) {
              rep(k, ansColumnSchedulesCount[day - 2][lineNum])
              {
                if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos2) {
                  diffScore3 += widths[lineNum] * 2;
                  break;
                }
              }
            }
            tmpMargin -= pos2 - sum;
            sum                    = pos2;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
        }
        else {
          CD32_NeighborNewPos1[ite2] = pos2;
          if (1 < day) {
            rep(k, ansColumnSchedulesCount[day - 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos2) {
                diffScore3 += widths[lineNum] * 2;
                break;
              }
            }
          }
          tmpMargin -= pos2 - sum;
          sum                    = pos2;
          CD32_TmpPosition[ite1] = sum;
          sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
          ite1++;
          CD32_TmpPosition[ite1] = sum;
        }
        ite2++;
      }
      else {
        if (pos32 < sum) {
          ite3++;
          continue;
        }
        if (sum + tmpMargin < pos3) {
          sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
          ite1++;
          CD32_TmpPosition[ite1] = sum;
          continue;
        }
        bool isConnect = false;
        if (day < n - 2) {
          rep(k, ansColumnSchedulesCount[day + 2][lineNum])
          {
            if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos32) {
              diffScore3 -= widths[lineNum] * 2;
              isConnect = true;
              break;
            }
          }
        }
        if (isConnect) {
          if (pos32 <= sum + tmpMargin) {
            diffScore3 += widths[lineNum] * 2;
            tmpMargin -= pos32 - sum;
            sum                    = pos32;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
          else {
            CD32_NeighborNewPos2[ite3] = pos3;
            if (day < n - 2) {
              rep(k, ansColumnSchedulesCount[day + 2][lineNum])
              {
                if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos3) {
                  diffScore3 += widths[lineNum] * 2;
                  break;
                }
              }
            }
            tmpMargin -= pos3 - sum;
            sum                    = pos3;
            CD32_TmpPosition[ite1] = sum;
            sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
            ite1++;
            CD32_TmpPosition[ite1] = sum;
          }
        }
        else {
          CD32_NeighborNewPos2[ite3] = pos3;
          if (day < n - 2) {
            rep(k, ansColumnSchedulesCount[day + 2][lineNum])
            {
              if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos3) {
                diffScore3 += widths[lineNum] * 2;
                break;
              }
            }
          }
          tmpMargin -= pos3 - sum;
          sum                    = pos3;
          CD32_TmpPosition[ite1] = sum;
          sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
          ite1++;
          CD32_TmpPosition[ite1] = sum;
        }
        ite3++;
      }
    }
    else if (ite2 < beforeCount) {
      if (pos22 < sum) {
        ite2++;
        continue;
      }
      if (sum + tmpMargin < pos2) {
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
        continue;
      }
      bool isConnect = false;
      if (1 < day) {
        rep(k, ansColumnSchedulesCount[day - 2][lineNum])
        {
          if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos22) {
            diffScore3 -= widths[lineNum] * 2;
            isConnect = true;
            break;
          }
        }
      }
      if (isConnect && pos22 <= sum + tmpMargin) {
        diffScore3 += widths[lineNum] * 2;
        tmpMargin -= pos22 - sum;
        sum                    = pos22;
        CD32_TmpPosition[ite1] = sum;
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
      }
      else {
        CD32_NeighborNewPos1[ite2] = pos2;
        if (1 < day) {
          rep(k, ansColumnSchedulesCount[day - 2][lineNum])
          {
            if (ansColumnSchedulesPosition[day - 2][lineNum][k] == pos2) {
              diffScore3 += widths[lineNum] * 2;
              break;
            }
          }
        }
        tmpMargin -= pos2 - sum;
        sum                    = pos2;
        CD32_TmpPosition[ite1] = sum;
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
      }
      ite2++;
    }
    else if (ite3 < afterCount) {
      if (pos32 < sum) {
        ite3++;
        continue;
      }
      if (sum + tmpMargin < pos3) {
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
        continue;
      }

      bool isConnect = false;
      if (day < n - 2) {
        rep(k, ansColumnSchedulesCount[day + 2][lineNum])
        {
          if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos32) {
            diffScore3 -= widths[lineNum] * 2;
            isConnect = true;
            break;
          }
        }
      }

      if (sum + tmpMargin < pos32 || !isConnect) {
        CD32_NeighborNewPos2[ite3] = pos3;
        if (day < n - 2) {
          rep(k, ansColumnSchedulesCount[day + 2][lineNum])
          {
            if (ansColumnSchedulesPosition[day + 2][lineNum][k] == pos3) {
              diffScore3 += widths[lineNum] * 2;
              break;
            }
          }
        }
        tmpMargin -= pos3 - sum;
        sum                    = pos3;
        CD32_TmpPosition[ite1] = sum;
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
      }
      else {
        if (isConnect) diffScore3 += widths[lineNum] * 2;
        tmpMargin -= pos32 - sum;
        sum                    = pos32;
        CD32_TmpPosition[ite1] = sum;
        sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
        ite1++;
        CD32_TmpPosition[ite1] = sum;
      }

      ite3++;
    }
    else {
      sum += preCalcScheduleSizes[day][CD32_Members[ite1]][lineNum];
      ite1++;
      CD32_TmpPosition[ite1] = sum;
    }
  }

  CD32_TmpPosition[memberCount] = w;

  ite1 = 1;
  ite2 = 1;
  ite3 = 1;
  while (ite1 < memberCount && ite2 < beforeCount) {
    if (CD32_TmpPosition[ite1] == CD32_NeighborNewPos1[ite2]) {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      ite2++;
    }
    else if (CD32_TmpPosition[ite1] < CD32_NeighborNewPos1[ite2]) {
      ite1++;
    }
    else {
      ite2++;
    }
  }

  ite1 = 1;
  ite2 = 1;
  ite3 = 1;
  while (ite1 < memberCount && ite3 < afterCount) {
    if (CD32_TmpPosition[ite1] == CD32_NeighborNewPos2[ite3]) {
      diffScore3 += widths[lineNum] * 2;
      ite1++;
      ite3++;
    }
    else if (CD32_TmpPosition[ite1] < CD32_NeighborNewPos2[ite3]) {
      ite1++;
    }
    else {
      ite3++;
    }
  }

  return diffScore3;
}

double M43_start_temp = 10000.1;
double M43_end_temp   = 0.1;
double M43_startTime  = 0;
double M43_nowTime    = 0;
double M43_timeLimit  = 0;
int M43_karinas[32]   = {};

int M43_kouhos[MAX_N];
int M43_kouhoCount = 0;
int M43_nonKouhos[MAX_N];
int M43_nonKouhoCount    = 0;
int M43_neighborPosCount = 0;

int M43_tmpAnsNextLine[MAX_N];
int M43_tmpAnsNextLinePosition[MAX_N];
int M43_tmpAnsNextCount = 0;
int M43_tmpAnsCurrentLine[MAX_N];
int M43_tmpAnsCurrentLinePosition[MAX_N];
int M43_tmpAnsCurrentCount = 0;

int M43_tmpAnsNextLineBeforePosition[MAX_N];
int M43_tmpAnsNextLineAfterPosition[MAX_N];
int M43_tmpAnsCurrentLineBeforePosition[MAX_N];
int M43_tmpAnsCurrentLineAfterPosition[MAX_N];

// 1対多スワップ
int M431Count;
void Method4_3_1()
{
  int raD          = randxor() % d;
  int raN          = randxor() % n;
  int lineNum      = ansColumnNum[raD][raN];
  int lineCapacity = w;
  rep(k, ansColumnSchedulesCount[raD][lineNum])
  {
    int num = ansColumnSchedules[raD][lineNum][k];
    if (num != raN) { lineCapacity -= preCalcScheduleSizes[raD][num][lineNum]; }
  }

  int nextLine = randxor() % ansLineCount[raD];
  while (nextLine == lineNum) {
    nextLine = randxor() % ansLineCount[raD];
  }
  int beforeCount         = ansColumnSchedulesCount[raD][lineNum];
  int afterCount          = ansColumnSchedulesCount[raD][nextLine];
  M43_kouhoCount          = 0;
  M43_nonKouhoCount       = 0;
  int nextLineCapacity    = w;
  int maxNextLineCapacity = w;
  rep(k, ansColumnSchedulesCount[raD][nextLine])
  {
    int num = ansColumnSchedules[raD][nextLine][k];
    nextLineCapacity -= preCalcScheduleSizes[raD][num][nextLine];
    if (preCalcScheduleSizes[raD][num][lineNum] > lineCapacity) {
      maxNextLineCapacity -= preCalcScheduleSizes[raD][num][nextLine];
      M43_nonKouhos[M43_nonKouhoCount] = num;
      M43_nonKouhoCount++;
    }
    else {
      M43_kouhos[M43_kouhoCount] = num;
      M43_kouhoCount++;
    }
  }

  if (M43_kouhoCount > 30) return;
  // if (M43_kouhoCount == 0 && afterCount > 0) return;
  if (preCalcScheduleSizes[raD][raN][nextLine] > maxNextLineCapacity) return;

  int karinaCount = 0;
  if (M43_kouhoCount == 0) {
    M43_karinas[karinaCount] = 0;
    karinaCount++;
  }
  else if (M43_kouhoCount <= 5) {
    rep(karina, (1 << M43_kouhoCount))
    {
      // if (karina == 0) continue;
      M43_karinas[karinaCount] = karina;
      karinaCount++;
    }
  }
  else {
    rep(aespa, 32)
    {
      // int karina               = randxor() % ((1 << M43_kouhoCount) - 1) + 1;
      int karina               = randxor() % ((1 << M43_kouhoCount));
      M43_karinas[karinaCount] = karina;
      karinaCount++;
    }
  }
  FisherYates(M43_karinas, karinaCount);
  rep(my, karinaCount)
  {
    int karina    = M43_karinas[my];
    int nextSpace = 0;
    int needSpace = 0;
    rep(jj, M43_kouhoCount)
    {
      if (karina & (1 << jj)) {
        int j = M43_kouhos[jj];
        nextSpace += preCalcScheduleSizes[raD][j][nextLine];
        needSpace += preCalcScheduleSizes[raD][j][lineNum];
        if (needSpace > lineCapacity) { break; }
      }
    }
    if (needSpace <= lineCapacity && preCalcScheduleSizes[raD][raN][nextLine] <= nextLineCapacity + nextSpace) {
      // int widths2[MAX_LINECOUNT] = {};
      // rep(k, ansBaseLineCount) {
      //  widths2[k] = widths[k];
      //  widths[k] *= widths[k];
      //}

      // スワップ可能
      int moveCount = 0;
      rep(jj, M43_kouhoCount)
      {
        if (karina & (1 << jj)) { moveCount++; }
      }
      int diffScore1 = 0;
      if (M43_kouhoCount == 0 || karina == 0) {
        if (beforeCount > 1) {
          if (raD == 0 || raD == d - 1) {
            diffScore1 += widths[lineNum];
          }
          else {
            diffScore1 += widths[lineNum] * 2;
          }
        }
        if (afterCount > 0) {
          if (raD == 0 || raD == d - 1) {
            diffScore1 -= widths[nextLine];
          }
          else {
            diffScore1 -= widths[nextLine] * 2;
          }
        }
      }
      else {
        diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]) * 2;
        if (raD == 0 || raD == d - 1) { diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]); }
      }

      int diffScore2 = 0;
      if (raD > 0) {
        diffScore2 += CalcDiffScore2(raD, raD - 1, lineNum);
        diffScore2 += CalcDiffScore2(raD, raD - 1, nextLine);
      }
      if (raD < d - 1) {
        diffScore2 += CalcDiffScore2(raD, raD + 1, lineNum);
        diffScore2 += CalcDiffScore2(raD, raD + 1, nextLine);
      }

      // 10回シャッフルnextLine
      M43_neighborPosCount = 0;
      if (raD > 0) {
        srep(k, 1, ansColumnSchedulesCount[raD - 1][nextLine])
        {
          CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD - 1][nextLine][k];
          M43_neighborPosCount++;
        }
      }
      if (raD < d - 1) {
        srep(k, 1, ansColumnSchedulesCount[raD + 1][nextLine])
        {
          CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD + 1][nextLine][k];
          M43_neighborPosCount++;
        }
      }
      sort(CD3_NeighborPos, CD3_NeighborPos + M43_neighborPosCount);
      M43_nonKouhos[M43_nonKouhoCount] = raN;
      M43_nonKouhoCount++;
      rep(jj, M43_kouhoCount)
      {
        if ((karina & (1 << jj)) == 0) {
          int j                            = M43_kouhos[jj];
          M43_nonKouhos[M43_nonKouhoCount] = j;
          M43_nonKouhoCount++;
        }
      }
      rep(i, M43_nonKouhoCount)
      {
        CD3_Members[i] = M43_nonKouhos[i];
      }
      int margin = w;
      rep(i, M43_nonKouhoCount)
      {
        int num = M43_nonKouhos[i];
        margin -= preCalcScheduleSizes[raD][num][nextLine];
      }
      int diffScore3 = -INT_INF;
      rep(winter, 10)
      {
        if (M43_nonKouhoCount == 1 && winter >= 1) break;
        if (M43_nonKouhoCount == 2 && winter >= 2) break;
        if (M43_nonKouhoCount == 3 && winter >= 6) break;
        int tmpDiffScore3 = CalcDiffScore3(winter, M43_nonKouhoCount, M43_neighborPosCount, raD, nextLine, margin);
        if (tmpDiffScore3 > diffScore3) {
          diffScore3 = tmpDiffScore3;
          rep(i, M43_nonKouhoCount)
          {
            M43_tmpAnsNextLine[i]         = CD3_Members[i];
            M43_tmpAnsNextLinePosition[i] = CD3_TmpPosition[i];
          }
          M43_tmpAnsNextLinePosition[M43_nonKouhoCount] = CD3_TmpPosition[M43_nonKouhoCount];
          M43_tmpAnsNextCount                           = M43_nonKouhoCount;
        }
      }

      // 10回シャッフルlineNum
      M43_neighborPosCount = 0;
      if (raD > 0) {
        srep(k, 1, ansColumnSchedulesCount[raD - 1][lineNum])
        {
          CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD - 1][lineNum][k];
          M43_neighborPosCount++;
        }
      }
      if (raD < d - 1) {
        srep(k, 1, ansColumnSchedulesCount[raD + 1][lineNum])
        {
          CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD + 1][lineNum][k];
          M43_neighborPosCount++;
        }
      }
      sort(CD3_NeighborPos, CD3_NeighborPos + M43_neighborPosCount);
      int tmpKouhos[MAX_N];
      int tmpKouhoCount = 0;
      rep(jj, M43_kouhoCount)
      {
        if (karina & (1 << jj)) {
          int j                    = M43_kouhos[jj];
          tmpKouhos[tmpKouhoCount] = j;
          tmpKouhoCount++;
        }
      }
      rep(j, tmpKouhoCount)
      {
        M43_kouhos[j] = tmpKouhos[j];
      }
      M43_kouhoCount = tmpKouhoCount;
      rep(k, ansColumnSchedulesCount[raD][lineNum])
      {
        int num = ansColumnSchedules[raD][lineNum][k];
        if (num != raN) {
          M43_kouhos[M43_kouhoCount] = num;
          M43_kouhoCount++;
        }
      }
      rep(i, M43_kouhoCount)
      {
        CD3_Members[i] = M43_kouhos[i];
      }
      margin = w;
      rep(i, M43_kouhoCount)
      {
        int num = M43_kouhos[i];
        margin -= preCalcScheduleSizes[raD][num][lineNum];
      }
      int diffScore4 = -INT_INF;
      rep(winter, 10)
      {
        if (M43_kouhoCount == 1 && winter >= 1) break;
        if (M43_kouhoCount == 2 && winter >= 2) break;
        if (M43_kouhoCount == 3 && winter >= 6) break;
        int tmpDiffScore4 = CalcDiffScore3(winter, M43_kouhoCount, M43_neighborPosCount, raD, lineNum, margin);
        if (tmpDiffScore4 > diffScore4) {
          diffScore4 = tmpDiffScore4;
          rep(i, M43_kouhoCount)
          {
            M43_tmpAnsCurrentLine[i]         = CD3_Members[i];
            M43_tmpAnsCurrentLinePosition[i] = CD3_TmpPosition[i];
          }
          M43_tmpAnsCurrentLinePosition[M43_kouhoCount] = CD3_TmpPosition[M43_kouhoCount];
          M43_tmpAnsCurrentCount                        = M43_kouhoCount;
        }
      }

      int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

      double temp = (M43_start_temp + (M43_end_temp - M43_start_temp) * M43_nowTime / M43_timeLimit);
      // double timeRatio = (M43_nowTime - M43_startTime) / (M43_timeLimit - M43_startTime);
      // double temp      = pow(M43_start_temp, 1.0 - timeRatio) * pow(M43_end_temp, timeRatio);
      // cout << temp << endl;
      const double prob = exp((double)totalDiffScore * 100 / temp);

      // rep(k, ansBaseLineCount) {
      //  widths[k] = widths2[k];
      //}

      if (prob > rand01()) {
        M431Count++;
        ansColumnSchedulesCount[raD][nextLine] = M43_tmpAnsNextCount;
        rep(i, M43_tmpAnsNextCount)
        {
          int num                                      = M43_tmpAnsNextLine[i];
          ansColumnNum[raD][num]                       = nextLine;
          ansColumnSchedules[raD][nextLine][i]         = num;
          ansColumnSchedulesPosition[raD][nextLine][i] = M43_tmpAnsNextLinePosition[i];
        }
        ansColumnSchedulesPosition[raD][nextLine][M43_tmpAnsNextCount] = M43_tmpAnsNextLinePosition[M43_tmpAnsNextCount];

        ansColumnSchedulesCount[raD][lineNum] = M43_tmpAnsCurrentCount;
        rep(i, M43_tmpAnsCurrentCount)
        {
          int num                                     = M43_tmpAnsCurrentLine[i];
          ansColumnNum[raD][num]                      = lineNum;
          ansColumnSchedules[raD][lineNum][i]         = num;
          ansColumnSchedulesPosition[raD][lineNum][i] = M43_tmpAnsCurrentLinePosition[i];
        }
        ansColumnSchedulesPosition[raD][lineNum][M43_tmpAnsCurrentCount] = M43_tmpAnsCurrentLinePosition[M43_tmpAnsCurrentCount];

        ansScore -= totalDiffScore;
        if (ansScore < real_ansScore) {
          // cout << ansScore << ' ' << M43_nowTime << endl;
          CopyToRealAns();
          CopyToReal_M42();
        }
      }
      else {
        ;
      }

      break;
    }
  }
}

// 1対多スワップ（柔軟）
void Method4_3_2()
{
  int raD          = randxor() % d;
  int raN          = randxor() % n;
  int lineNum      = ansColumnNum[raD][raN];
  int lineCapacity = w;
  rep(k, ansColumnSchedulesCount[raD][lineNum])
  {
    int num = ansColumnSchedules[raD][lineNum][k];
    if (num != raN) { lineCapacity -= preCalcScheduleSizes[raD][num][lineNum]; }
  }

  int nextLine = randxor() % ansLineCount[raD];
  while (nextLine == lineNum) {
    nextLine = randxor() % ansLineCount[raD];
  }
  int beforeCount         = ansColumnSchedulesCount[raD][lineNum];
  int afterCount          = ansColumnSchedulesCount[raD][nextLine];
  M43_kouhoCount          = 0;
  M43_nonKouhoCount       = 0;
  int nextLineCapacity    = w;
  int maxNextLineCapacity = w;
  rep(k, ansColumnSchedulesCount[raD][nextLine])
  {
    int num = ansColumnSchedules[raD][nextLine][k];
    nextLineCapacity -= preCalcScheduleSizes[raD][num][nextLine];
    if (preCalcScheduleSizes[raD][num][lineNum] > lineCapacity) {
      maxNextLineCapacity -= preCalcScheduleSizes[raD][num][nextLine];
      M43_nonKouhos[M43_nonKouhoCount] = num;
      M43_nonKouhoCount++;
    }
    else {
      M43_kouhos[M43_kouhoCount] = num;
      M43_kouhoCount++;
    }
  }

  if (M43_kouhoCount > 30) return;
  // if (kouhoCount == 0 && afterCount > 0) return;
  if (M43_kouhoCount == 0) return;
  if (preCalcScheduleSizes[raD][raN][nextLine] > maxNextLineCapacity) return;

  int karinaCount = 0;
  if (M43_kouhoCount == 0) {
    M43_karinas[karinaCount] = 0;
    karinaCount++;
  }
  else if (M43_kouhoCount <= 5) {
    rep(karina, (1 << M43_kouhoCount))
    {
      if (karina == 0) continue;
      M43_karinas[karinaCount] = karina;
      karinaCount++;
    }
  }
  else {
    rep(aespa, 32)
    {
      int karina               = randxor() % ((1 << M43_kouhoCount) - 1) + 1;
      M43_karinas[karinaCount] = karina;
      karinaCount++;
    }
  }
  FisherYates(M43_karinas, karinaCount);
  rep(my, karinaCount)
  {
    int karina    = M43_karinas[my];
    int nextSpace = 0;
    int needSpace = 0;
    rep(jj, M43_kouhoCount)
    {
      if (karina & (1 << jj)) {
        int j = M43_kouhos[jj];
        nextSpace += preCalcScheduleSizes[raD][j][nextLine];
        needSpace += preCalcScheduleSizes[raD][j][lineNum];
        if (needSpace > lineCapacity) { break; }
      }
    }
    if (needSpace <= lineCapacity && preCalcScheduleSizes[raD][raN][nextLine] <= nextLineCapacity + nextSpace) {
      // スワップ可能
      int moveCount = 0;
      rep(jj, M43_kouhoCount)
      {
        if (karina & (1 << jj)) { moveCount++; }
      }
      int diffScore1 = 0;
      if (M43_kouhoCount == 0) {
        if (beforeCount > 1) {
          diffScore1 = widths[lineNum] * 2;
          if (raD == 0 || raD == d - 1) { diffScore1 = widths[lineNum]; }
        }
      }
      else {
        diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]) * 2;
        if (raD == 0 || raD == d - 1) { diffScore1 = (moveCount - 1) * (widths[nextLine] - widths[lineNum]); }
      }

      int diffScore2 = 0;
      if (raD > 0) {
        diffScore2 += CalcDiffScore2(raD, raD - 1, lineNum);
        diffScore2 += CalcDiffScore2(raD, raD - 1, nextLine);
      }
      if (raD < d - 1) {
        diffScore2 += CalcDiffScore2(raD, raD + 1, lineNum);
        diffScore2 += CalcDiffScore2(raD, raD + 1, nextLine);
      }

      // 10回シャッフルnextLine
      M43_nonKouhos[M43_nonKouhoCount] = raN;
      M43_nonKouhoCount++;
      rep(jj, M43_kouhoCount)
      {
        if ((karina & (1 << jj)) == 0) {
          int j                            = M43_kouhos[jj];
          M43_nonKouhos[M43_nonKouhoCount] = j;
          M43_nonKouhoCount++;
        }
      }
      rep(i, M43_nonKouhoCount)
      {
        CD32_Members[i] = M43_nonKouhos[i];
      }
      int margin = w;
      rep(i, M43_nonKouhoCount)
      {
        int num = M43_nonKouhos[i];
        margin -= preCalcScheduleSizes[raD][num][nextLine];
      }
      int diffScore3 = -INT_INF;
      rep(winter, 10)
      {
        if (M43_nonKouhoCount == 1 && winter >= 1) break;
        if (M43_nonKouhoCount == 2 && winter >= 2) break;
        if (M43_nonKouhoCount == 3 && winter >= 6) break;
        int tmpDiffScore3 = CalcDiffScore3_2(winter, M43_nonKouhoCount, raD, nextLine, 0);
        if (tmpDiffScore3 > diffScore3) {
          diffScore3 = tmpDiffScore3;
          rep(i, M43_nonKouhoCount)
          {
            M43_tmpAnsNextLine[i]         = CD32_Members[i];
            M43_tmpAnsNextLinePosition[i] = CD32_TmpPosition[i];
          }
          M43_tmpAnsNextLinePosition[M43_nonKouhoCount] = CD32_TmpPosition[M43_nonKouhoCount];
          M43_tmpAnsNextCount                           = M43_nonKouhoCount;
        }

        if (raD > 0) {
          beforeCount = ansColumnSchedulesCount[raD - 1][nextLine];
          rep(k, beforeCount + 1)
          {
            M43_tmpAnsNextLineBeforePosition[k] = CD32_NeighborNewPos1[k];
          }
        }
        if (raD < d - 1) {
          int afterCount = ansColumnSchedulesCount[raD + 1][nextLine];
          rep(k, afterCount + 1)
          {
            M43_tmpAnsNextLineAfterPosition[k] = CD32_NeighborNewPos2[k];
          }
        }
      }

      // 10回シャッフルlineNum
      int tmpKouhos[MAX_N];
      int tmpKouhoCount = 0;
      rep(jj, M43_kouhoCount)
      {
        if (karina & (1 << jj)) {
          int j                    = M43_kouhos[jj];
          tmpKouhos[tmpKouhoCount] = j;
          tmpKouhoCount++;
        }
      }
      rep(j, tmpKouhoCount)
      {
        M43_kouhos[j] = tmpKouhos[j];
      }
      M43_kouhoCount = tmpKouhoCount;
      rep(k, ansColumnSchedulesCount[raD][lineNum])
      {
        int num = ansColumnSchedules[raD][lineNum][k];
        if (num != raN) {
          M43_kouhos[M43_kouhoCount] = num;
          M43_kouhoCount++;
        }
      }
      rep(i, M43_kouhoCount)
      {
        CD32_Members[i] = M43_kouhos[i];
      }
      margin = w;
      rep(i, M43_kouhoCount)
      {
        int num = M43_kouhos[i];
        margin -= preCalcScheduleSizes[raD][num][lineNum];
      }
      int diffScore4 = -INT_INF;
      rep(winter, 10)
      {
        if (M43_kouhoCount == 1 && winter >= 1) break;
        if (M43_kouhoCount == 2 && winter >= 2) break;
        if (M43_kouhoCount == 3 && winter >= 6) break;
        int tmpDiffScore4 = CalcDiffScore3_2(winter, M43_kouhoCount, raD, lineNum, 0);
        if (tmpDiffScore4 > diffScore4) {
          diffScore4 = tmpDiffScore4;
          rep(i, M43_kouhoCount)
          {
            M43_tmpAnsCurrentLine[i]         = CD32_Members[i];
            M43_tmpAnsCurrentLinePosition[i] = CD32_TmpPosition[i];
          }
          M43_tmpAnsCurrentLinePosition[M43_kouhoCount] = CD32_TmpPosition[M43_kouhoCount];
          M43_tmpAnsCurrentCount                        = M43_kouhoCount;

          if (raD > 0) {
            beforeCount = ansColumnSchedulesCount[raD - 1][lineNum];
            rep(k, beforeCount + 1)
            {
              M43_tmpAnsCurrentLineBeforePosition[k] = CD32_NeighborNewPos1[k];
            }
          }
          if (raD < d - 1) {
            int afterCount = ansColumnSchedulesCount[raD + 1][lineNum];
            rep(k, afterCount + 1)
            {
              M43_tmpAnsCurrentLineAfterPosition[k] = CD32_NeighborNewPos2[k];
            }
          }
        }
      }

      int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

      double temp       = (M43_start_temp + (M43_end_temp - M43_start_temp) * M43_nowTime / M43_timeLimit);
      const double prob = exp((double)totalDiffScore * 100 / temp);

      if (prob > rand01()) {
        ansColumnSchedulesCount[raD][nextLine] = M43_tmpAnsNextCount;
        rep(i, M43_tmpAnsNextCount)
        {
          int num                                      = M43_tmpAnsNextLine[i];
          ansColumnNum[raD][num]                       = nextLine;
          ansColumnSchedules[raD][nextLine][i]         = num;
          ansColumnSchedulesPosition[raD][nextLine][i] = M43_tmpAnsNextLinePosition[i];
        }
        ansColumnSchedulesPosition[raD][nextLine][M43_tmpAnsNextCount] = M43_tmpAnsNextLinePosition[M43_tmpAnsNextCount];

        ansColumnSchedulesCount[raD][lineNum] = M43_tmpAnsCurrentCount;
        rep(i, M43_tmpAnsCurrentCount)
        {
          int num                                     = M43_tmpAnsCurrentLine[i];
          ansColumnNum[raD][num]                      = lineNum;
          ansColumnSchedules[raD][lineNum][i]         = num;
          ansColumnSchedulesPosition[raD][lineNum][i] = M43_tmpAnsCurrentLinePosition[i];
        }
        ansColumnSchedulesPosition[raD][lineNum][M43_tmpAnsCurrentCount] = M43_tmpAnsCurrentLinePosition[M43_tmpAnsCurrentCount];

        if (raD > 0) {
          beforeCount = ansColumnSchedulesCount[raD - 1][lineNum];
          rep(k, beforeCount + 1)
          {
            ansColumnSchedulesPosition[raD - 1][lineNum][k] = M43_tmpAnsCurrentLineBeforePosition[k];
          }
        }
        if (raD < d - 1) {
          int afterCount = ansColumnSchedulesCount[raD + 1][lineNum];
          rep(k, afterCount + 1)
          {
            ansColumnSchedulesPosition[raD + 1][lineNum][k] = M43_tmpAnsCurrentLineAfterPosition[k];
          }
        }
        if (raD > 0) {
          beforeCount = ansColumnSchedulesCount[raD - 1][nextLine];
          rep(k, beforeCount + 1)
          {
            ansColumnSchedulesPosition[raD - 1][nextLine][k] = M43_tmpAnsNextLineBeforePosition[k];
          }
        }
        if (raD < d - 1) {
          int afterCount = ansColumnSchedulesCount[raD + 1][nextLine];
          rep(k, afterCount + 1)
          {
            ansColumnSchedulesPosition[raD + 1][nextLine][k] = M43_tmpAnsNextLineAfterPosition[k];
          }
        }

        ansScore -= totalDiffScore;
        if (ansScore < real_ansScore) {
          CopyToRealAns();
          CopyToReal_M42();
        }
      }
      else {
        ;
      }

      break;
    }
  }
}

// 列スワップ
void Method4_3_3()
{
  int raD   = randxor() % d;
  int line1 = randxor() % ansLineCount[raD];
  int line2 = randxor() % ansLineCount[raD];
  while (line2 == line1) {
    line2 = randxor() % ansLineCount[raD];
  }
  if (ansColumnSchedulesCount[raD][line1] == 0 || ansColumnSchedulesCount[raD][line2] == 0) return;

  {
    int line1NextSum = 0;
    rep(k, ansColumnSchedulesCount[raD][line1])
    {
      int num = ansColumnSchedules[raD][line1][k];
      line1NextSum += preCalcScheduleSizes[raD][num][line2];
    }
    if (line1NextSum > w) return;
  }
  {
    int line2NextSum = 0;
    rep(k, ansColumnSchedulesCount[raD][line2])
    {
      int num = ansColumnSchedules[raD][line2][k];
      line2NextSum += preCalcScheduleSizes[raD][num][line1];
    }
    if (line2NextSum > w) return;
  }

  int diffScore1 = (ansColumnSchedulesCount[raD][line1] - ansColumnSchedulesCount[raD][line2]) * (widths[line1] - widths[line2]) * 2;
  if (raD == 0 || raD == d - 1) diffScore1 = (ansColumnSchedulesCount[raD][line1] - ansColumnSchedulesCount[raD][line2]) * (widths[line1] - widths[line2]);

  int diffScore2 = 0;
  if (raD > 0) {
    diffScore2 += CalcDiffScore2(raD, raD - 1, line1);
    diffScore2 += CalcDiffScore2(raD, raD - 1, line2);
  }
  if (raD < d - 1) {
    diffScore2 += CalcDiffScore2(raD, raD + 1, line1);
    diffScore2 += CalcDiffScore2(raD, raD + 1, line2);
  }

  // 10回シャッフルline2
  M43_neighborPosCount = 0;
  if (raD > 0) {
    srep(k, 1, ansColumnSchedulesCount[raD - 1][line2])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD - 1][line2][k];
      M43_neighborPosCount++;
    }
  }
  if (raD < d - 1) {
    srep(k, 1, ansColumnSchedulesCount[raD + 1][line2])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD + 1][line2][k];
      M43_neighborPosCount++;
    }
  }
  sort(CD3_NeighborPos, CD3_NeighborPos + M43_neighborPosCount);

  M43_kouhoCount = ansColumnSchedulesCount[raD][line1];
  int margin     = w;
  rep(i, M43_kouhoCount)
  {
    CD3_Members[i] = ansColumnSchedules[raD][line1][i];
    margin -= preCalcScheduleSizes[raD][CD3_Members[i]][line2];
  }

  int diffScore3 = -1;
  rep(winter, 10)
  {
    if (M43_kouhoCount == 1 && winter >= 1) break;
    if (M43_kouhoCount == 2 && winter >= 2) break;
    if (M43_kouhoCount == 3 && winter >= 6) break;
    int tmpDiffScore3 = CalcDiffScore3(winter, M43_kouhoCount, M43_neighborPosCount, raD, line2, margin);
    if (tmpDiffScore3 > diffScore3) {
      diffScore3 = tmpDiffScore3;
      rep(i, M43_kouhoCount)
      {
        M43_tmpAnsNextLine[i]         = CD3_Members[i];
        M43_tmpAnsNextLinePosition[i] = CD3_TmpPosition[i];
      }
      M43_tmpAnsNextLinePosition[M43_kouhoCount] = CD3_TmpPosition[M43_kouhoCount];
      M43_tmpAnsNextCount                        = M43_kouhoCount;
    }
  }

  // 10回シャッフルline1
  M43_neighborPosCount = 0;
  if (raD > 0) {
    srep(k, 1, ansColumnSchedulesCount[raD - 1][line1])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD - 1][line1][k];
      M43_neighborPosCount++;
    }
  }
  if (raD < d - 1) {
    srep(k, 1, ansColumnSchedulesCount[raD + 1][line1])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD + 1][line1][k];
      M43_neighborPosCount++;
    }
  }
  sort(CD3_NeighborPos, CD3_NeighborPos + M43_neighborPosCount);

  M43_kouhoCount = ansColumnSchedulesCount[raD][line2];
  margin         = w;
  rep(i, M43_kouhoCount)
  {
    CD3_Members[i] = ansColumnSchedules[raD][line2][i];
    margin -= preCalcScheduleSizes[raD][CD3_Members[i]][line1];
  }

  int diffScore4 = -1;
  rep(winter, 10)
  {
    if (M43_kouhoCount == 1 && winter >= 1) break;
    if (M43_kouhoCount == 2 && winter >= 2) break;
    if (M43_kouhoCount == 3 && winter >= 6) break;
    int tmpDiffScore4 = CalcDiffScore3(winter, M43_kouhoCount, M43_neighborPosCount, raD, line1, margin);
    if (tmpDiffScore4 > diffScore4) {
      diffScore4 = tmpDiffScore4;
      rep(i, M43_kouhoCount)
      {
        M43_tmpAnsCurrentLine[i]         = CD3_Members[i];
        M43_tmpAnsCurrentLinePosition[i] = CD3_TmpPosition[i];
      }
      M43_tmpAnsCurrentLinePosition[M43_kouhoCount] = CD3_TmpPosition[M43_kouhoCount];
      M43_tmpAnsCurrentCount                        = M43_kouhoCount;
    }
  }

  int totalDiffScore = diffScore1 + diffScore2 + diffScore3 + diffScore4;

  double temp       = (M43_start_temp + (M43_end_temp - M43_start_temp) * M43_nowTime / M43_timeLimit);
  const double prob = exp((double)totalDiffScore * 10000 / temp);

  if (totalDiffScore >= 0) {
    // if (prob > rand01()) {
    // cout << "OK";
    ansColumnSchedulesCount[raD][line2] = M43_tmpAnsNextCount;
    rep(i, M43_tmpAnsNextCount)
    {
      int num                                   = M43_tmpAnsNextLine[i];
      ansColumnNum[raD][num]                    = line2;
      ansColumnSchedules[raD][line2][i]         = num;
      ansColumnSchedulesPosition[raD][line2][i] = M43_tmpAnsNextLinePosition[i];
    }
    ansColumnSchedulesPosition[raD][line2][M43_tmpAnsNextCount] = M43_tmpAnsNextLinePosition[M43_tmpAnsNextCount];

    ansColumnSchedulesCount[raD][line1] = M43_tmpAnsCurrentCount;
    rep(i, M43_tmpAnsCurrentCount)
    {
      int num                                   = M43_tmpAnsCurrentLine[i];
      ansColumnNum[raD][num]                    = line1;
      ansColumnSchedules[raD][line1][i]         = num;
      ansColumnSchedulesPosition[raD][line1][i] = M43_tmpAnsCurrentLinePosition[i];
    }
    ansColumnSchedulesPosition[raD][line1][M43_tmpAnsCurrentCount] = M43_tmpAnsCurrentLinePosition[M43_tmpAnsCurrentCount];

    ansScore -= totalDiffScore;
    if (ansScore < real_ansScore) {
      // cout << "XXXXXXXXXXXX" << endl;
      CopyToRealAns();
      CopyToReal_M42();
    }
  }
}

// 横線を移動
void Method4_3_4()
{
  int raD     = randxor() % d;
  int lineNum = randxor() % ansBaseLineCount;
  if (ansColumnSchedulesCount[raD][lineNum] <= 1) return;
  int raIndex = randxor() % ansColumnSchedulesCount[raD][lineNum];
  int raDir   = randxor() % 2;
  if (raIndex == 0) { raDir = 1; }
  if (raIndex == ansColumnSchedulesCount[raD][lineNum] - 1) { raDir = 0; }
  int raN    = ansColumnSchedules[raD][lineNum][raIndex];
  int margin = (ansColumnSchedulesPosition[raD][lineNum][raIndex + 1] - ansColumnSchedulesPosition[raD][lineNum][raIndex]) - preCalcScheduleSizes[raD][raN][lineNum];
  if (margin == 0) return;
  int moveAmount = randxor() % margin + 1;
  int beforePos  = ansColumnSchedulesPosition[raD][lineNum][raIndex];
  int afterPos   = ansColumnSchedulesPosition[raD][lineNum][raIndex] + moveAmount;
  if (raDir == 1) {
    beforePos = ansColumnSchedulesPosition[raD][lineNum][raIndex + 1];
    afterPos  = ansColumnSchedulesPosition[raD][lineNum][raIndex + 1] - moveAmount;
  }
  int diffScore = 0;
  if (raD > 0) {
    rep(k, ansColumnSchedulesCount[raD - 1][lineNum])
    {
      if (ansColumnSchedulesPosition[raD - 1][lineNum][k] == beforePos) diffScore -= widths[lineNum] * 2;
      if (ansColumnSchedulesPosition[raD - 1][lineNum][k] == afterPos) diffScore += widths[lineNum] * 2;
    }
  }
  if (raD < d - 1) {
    rep(k, ansColumnSchedulesCount[raD + 1][lineNum])
    {
      if (ansColumnSchedulesPosition[raD + 1][lineNum][k] == beforePos) diffScore -= widths[lineNum] * 2;
      if (ansColumnSchedulesPosition[raD + 1][lineNum][k] == afterPos) diffScore += widths[lineNum] * 2;
    }
  }

  if (diffScore >= 0) {
    if (raDir == 0) {
      ansColumnSchedulesPosition[raD][lineNum][raIndex] = afterPos;
    }
    else {
      ansColumnSchedulesPosition[raD][lineNum][raIndex + 1] = afterPos;
    }
    ansScore -= diffScore;
    if (ansScore < real_ansScore) {
      CopyToRealAns();
      CopyToReal_M42();
    }
  }
}

// 横線を移動
void Method4_3_4_2()
{
  int raD     = randxor() % d;
  int lineNum = randxor() % ansBaseLineCount;
  if (ansColumnSchedulesCount[raD][lineNum] <= 1) return;
  int raIndex = randxor() % ansColumnSchedulesCount[raD][lineNum];
  int raDir   = randxor() % 2;
  if (raIndex == 0) { raDir = 1; }
  if (raIndex == ansColumnSchedulesCount[raD][lineNum] - 1) { raDir = 0; }
  int beforeLinePos = ansColumnSchedulesPosition[raD][lineNum][raIndex];
  if (raDir == 1) { beforeLinePos = ansColumnSchedulesPosition[raD][lineNum][raIndex + 1]; }

  int margin   = w;
  int startDay = raD;
  int endDay   = raD;
  drep(i, raD + 1)
  {
    int lineIndex = -1;
    srep(j, 1, ansColumnSchedulesCount[i][lineNum])
    {
      if (ansColumnSchedulesPosition[i][lineNum][j] == beforeLinePos) {
        lineIndex = j;
        break;
      }
    }
    if (lineIndex == -1) break;
    startDay = i;
    if (raDir == 0) {
      int raN = ansColumnSchedules[i][lineNum][lineIndex];
      margin  = min(margin, (ansColumnSchedulesPosition[i][lineNum][lineIndex + 1] - ansColumnSchedulesPosition[i][lineNum][lineIndex]) - preCalcScheduleSizes[i][raN][lineNum]);
    }
    else {
      int raN = ansColumnSchedules[i][lineNum][lineIndex - 1];
      margin  = min(margin, (ansColumnSchedulesPosition[i][lineNum][lineIndex] - ansColumnSchedulesPosition[i][lineNum][lineIndex - 1]) - preCalcScheduleSizes[i][raN][lineNum]);
    }
    if (margin == 0) return;
  }
  srep(i, raD + 1, d)
  {
    int lineIndex = -1;
    srep(j, 1, ansColumnSchedulesCount[i][lineNum])
    {
      if (ansColumnSchedulesPosition[i][lineNum][j] == beforeLinePos) {
        lineIndex = j;
        break;
      }
    }
    if (lineIndex == -1) break;
    endDay = i;
    if (raDir == 0) {
      int raN = ansColumnSchedules[i][lineNum][lineIndex];
      margin  = min(margin, (ansColumnSchedulesPosition[i][lineNum][lineIndex + 1] - ansColumnSchedulesPosition[i][lineNum][lineIndex]) - preCalcScheduleSizes[i][raN][lineNum]);
    }
    else {
      int raN = ansColumnSchedules[i][lineNum][lineIndex - 1];
      margin  = min(margin, (ansColumnSchedulesPosition[i][lineNum][lineIndex] - ansColumnSchedulesPosition[i][lineNum][lineIndex - 1]) - preCalcScheduleSizes[i][raN][lineNum]);
    }
    if (margin == 0) return;
  }

  int moveAmount   = randxor() % margin + 1;
  int afterLinePos = beforeLinePos + moveAmount;
  if (raDir == 1) { afterLinePos = beforeLinePos - moveAmount; }

  int diffScore = 0;
  if (startDay > 0) {
    rep(k, ansColumnSchedulesCount[startDay - 1][lineNum])
    {
      if (ansColumnSchedulesPosition[startDay - 1][lineNum][k] == beforeLinePos) diffScore -= widths[lineNum] * 2;
      if (ansColumnSchedulesPosition[startDay - 1][lineNum][k] == afterLinePos) diffScore += widths[lineNum] * 2;
    }
  }
  if (endDay < d - 1) {
    rep(k, ansColumnSchedulesCount[endDay + 1][lineNum])
    {
      if (ansColumnSchedulesPosition[endDay + 1][lineNum][k] == beforeLinePos) diffScore -= widths[lineNum] * 2;
      if (ansColumnSchedulesPosition[endDay + 1][lineNum][k] == afterLinePos) diffScore += widths[lineNum] * 2;
    }
  }

  if (diffScore >= 0) {
    srep(i, startDay, endDay + 1)
    {
      int ok = 0;
      srep(j, 1, ansColumnSchedulesCount[i][lineNum])
      {
        if (ansColumnSchedulesPosition[i][lineNum][j] == beforeLinePos) {
          ansColumnSchedulesPosition[i][lineNum][j] = afterLinePos;
          ok                                        = 1;
          break;
        }
      }
    }
    ansScore -= diffScore;
    if (ansScore < real_ansScore) {
      CopyToRealAns();
      CopyToReal_M42();
    }
  }
}

void Method4_3_5()
{
  int lineNum = randxor() % ansBaseLineCount;
  int raDir   = randxor() % 2;
  if (lineNum == 0) { raDir = 1; }
  if (lineNum == ansBaseLineCount - 1) { raDir = 0; }
  int margin = widths[lineNum] - 1;
  rep(i, d)
  {
    rep(k, ansColumnSchedulesCount[i][lineNum])
    {
      int num       = ansColumnSchedules[i][lineNum][k];
      int tmpMargin = widths[lineNum] - ((a[i][num] - 1) / (ansColumnSchedulesPosition[i][lineNum][k + 1] - ansColumnSchedulesPosition[i][lineNum][k]) + 1);
      if (tmpMargin < margin) margin = tmpMargin;
    }
    if (margin <= 0) break;
  }
  if (margin == 0) return;
  int moveAmount = randxor() % margin + 1;
  int diffScore  = 0;
  srep(i, 1, d)
  {
    int ite1 = 1;
    int ite2 = 1;
    while (ite1 < ansColumnSchedulesCount[i - 1][lineNum] || ite2 < ansColumnSchedulesCount[i][lineNum]) {
      if (ite1 == ansColumnSchedulesCount[i - 1][lineNum]) {
        diffScore += moveAmount;
        ite2++;
      }
      else if (ite2 == ansColumnSchedulesCount[i][lineNum]) {
        diffScore += moveAmount;
        ite1++;
      }
      else {
        if (ansColumnSchedulesPosition[i - 1][lineNum][ite1] == ansColumnSchedulesPosition[i][lineNum][ite2]) {
          ite1++;
          ite2++;
        }
        else if (ansColumnSchedulesPosition[i - 1][lineNum][ite1] < ansColumnSchedulesPosition[i][lineNum][ite2]) {
          diffScore += moveAmount;
          ite1++;
        }
        else {
          diffScore += moveAmount;
          ite2++;
        }
      }
    }
    if (raDir == 0) {
      ite1 = 1;
      ite2 = 1;
      while (ite1 < ansColumnSchedulesCount[i - 1][lineNum - 1] || ite2 < ansColumnSchedulesCount[i][lineNum - 1]) {
        if (ite1 == ansColumnSchedulesCount[i - 1][lineNum - 1]) {
          diffScore -= moveAmount;
          ite2++;
        }
        else if (ite2 == ansColumnSchedulesCount[i][lineNum - 1]) {
          diffScore -= moveAmount;
          ite1++;
        }
        else {
          if (ansColumnSchedulesPosition[i - 1][lineNum - 1][ite1] == ansColumnSchedulesPosition[i][lineNum - 1][ite2]) {
            ite1++;
            ite2++;
          }
          else if (ansColumnSchedulesPosition[i - 1][lineNum - 1][ite1] < ansColumnSchedulesPosition[i][lineNum - 1][ite2]) {
            diffScore -= moveAmount;
            ite1++;
          }
          else {
            diffScore -= moveAmount;
            ite2++;
          }
        }
      }
    }
    else {
      ite1 = 1;
      ite2 = 1;
      while (ite1 < ansColumnSchedulesCount[i - 1][lineNum + 1] || ite2 < ansColumnSchedulesCount[i][lineNum + 1]) {
        if (ite1 == ansColumnSchedulesCount[i - 1][lineNum + 1]) {
          diffScore -= moveAmount;
          ite2++;
        }
        else if (ite2 == ansColumnSchedulesCount[i][lineNum + 1]) {
          diffScore -= moveAmount;
          ite1++;
        }
        else {
          if (ansColumnSchedulesPosition[i - 1][lineNum + 1][ite1] == ansColumnSchedulesPosition[i][lineNum + 1][ite2]) {
            ite1++;
            ite2++;
          }
          else if (ansColumnSchedulesPosition[i - 1][lineNum + 1][ite1] < ansColumnSchedulesPosition[i][lineNum + 1][ite2]) {
            diffScore -= moveAmount;
            ite1++;
          }
          else {
            diffScore -= moveAmount;
            ite2++;
          }
        }
      }
    }
  }

  double temp       = (M43_start_temp + (M43_end_temp - M43_start_temp) * M43_nowTime / M43_timeLimit);
  const double prob = exp((double)diffScore * 1 / temp);

  // if (diffScore >= 0) {
  if (prob > rand01()) {
    if (raDir == 0) {
      widths[lineNum] -= moveAmount;
      widths[lineNum - 1] += moveAmount;
      rep(i, d)
      {
        rep(j, n)
        {
          preCalcScheduleSizes[i][j][lineNum]     = (a[i][j] - 1) / widths[lineNum] + 1;
          preCalcScheduleSizes[i][j][lineNum - 1] = (a[i][j] - 1) / widths[lineNum - 1] + 1;
        }
      }
    }
    else {
      widths[lineNum] -= moveAmount;
      widths[lineNum + 1] += moveAmount;
      rep(i, d)
      {
        rep(j, n)
        {
          preCalcScheduleSizes[i][j][lineNum]     = (a[i][j] - 1) / widths[lineNum] + 1;
          preCalcScheduleSizes[i][j][lineNum + 1] = (a[i][j] - 1) / widths[lineNum + 1] + 1;
        }
      }
    }
    ansScore -= diffScore;
    rep(i, d)
    {
      ansLinePos[i][0] = 0;
      rep(j, ansBaseLineCount)
      {
        ansLinePos[i][j + 1] = ansLinePos[i][j] + widths[j];
      }
    }
    if (ansScore < real_ansScore) {
      CopyToRealAns();
      CopyToReal_M42();
    }
  }
}

// 交換できるやつを交換
void Method4_3_6()
{
  int raD       = randxor() % d;
  int raN       = randxor() % n;
  int lineNum   = ansColumnNum[raD][raN];
  int lineIndex = 0;
  int raNSpace  = 0;
  rep(i, ansColumnSchedulesCount[raD][lineNum])
  {
    if (ansColumnSchedules[raD][lineNum][i] == raN) {
      lineIndex = i;
      raNSpace  = (ansColumnSchedulesPosition[raD][lineNum][i + 1] - ansColumnSchedulesPosition[raD][lineNum][i]) * widths[lineNum];
      break;
    }
  }
  int now = randxor() % n;
  rep(jisoo, n)
  {
    now = (now + 97) % n;
    if (now == raN) continue;
    int nextLine      = ansColumnNum[raD][now];
    int nextLineIndex = 0;
    int nowSpace      = 0;
    rep(i, ansColumnSchedulesCount[raD][nextLine])
    {
      if (ansColumnSchedules[raD][nextLine][i] == now) {
        nextLineIndex = i;
        nowSpace      = (ansColumnSchedulesPosition[raD][nextLine][i + 1] - ansColumnSchedulesPosition[raD][nextLine][i]) * widths[nextLine];
        break;
      }
    }
    if (a[raD][raN] <= nowSpace && a[raD][now] <= raNSpace) {
      ansColumnNum[raD][raN]                           = nextLine;
      ansColumnNum[raD][now]                           = lineNum;
      ansColumnSchedules[raD][lineNum][lineIndex]      = now;
      ansColumnSchedules[raD][nextLine][nextLineIndex] = raN;
      break;
    }
  }
}

// 列シャッフル
void Method4_3_7()
{
  int raD          = randxor() % d;
  int raN          = randxor() % n;
  int lineNum      = ansColumnNum[raD][raN];
  int lineCapacity = w;
  if (ansColumnSchedulesCount[raD][lineNum] == 1) return;

  int diffScore2 = 0;
  if (raD > 0) { diffScore2 += CalcDiffScore2(raD, raD - 1, lineNum); }
  if (raD < d - 1) { diffScore2 += CalcDiffScore2(raD, raD + 1, lineNum); }
  // 10回シャッフルlineNum
  M43_neighborPosCount = 0;
  if (raD > 0) {
    srep(k, 1, ansColumnSchedulesCount[raD - 1][lineNum])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD - 1][lineNum][k];
      M43_neighborPosCount++;
    }
  }
  if (raD < d - 1) {
    srep(k, 1, ansColumnSchedulesCount[raD + 1][lineNum])
    {
      CD3_NeighborPos[M43_neighborPosCount] = ansColumnSchedulesPosition[raD + 1][lineNum][k];
      M43_neighborPosCount++;
    }
  }

  M43_kouhoCount = ansColumnSchedulesCount[raD][lineNum];
  int margin     = w;
  rep(i, M43_kouhoCount)
  {
    CD3_Members[i] = ansColumnSchedules[raD][lineNum][i];
    margin -= preCalcScheduleSizes[raD][CD3_Members[i]][lineNum];
  }

  int diffScore4 = -1;
  rep(winter, 10)
  {
    if (M43_kouhoCount == 1 && winter >= 1) break;
    if (M43_kouhoCount == 2 && winter >= 2) break;
    if (M43_kouhoCount == 3 && winter >= 6) break;
    int tmpDiffScore4 = CalcDiffScore3(winter, M43_kouhoCount, M43_neighborPosCount, raD, lineNum, margin);
    if (tmpDiffScore4 > diffScore4) {
      diffScore4 = tmpDiffScore4;
      rep(i, M43_kouhoCount)
      {
        M43_tmpAnsCurrentLine[i]         = CD3_Members[i];
        M43_tmpAnsCurrentLinePosition[i] = CD3_TmpPosition[i];
      }
      M43_tmpAnsCurrentLinePosition[M43_kouhoCount] = CD3_TmpPosition[M43_kouhoCount];
      M43_tmpAnsCurrentCount                        = M43_kouhoCount;
    }
  }

  int totalDiffScore = diffScore2 + diffScore4;

  double temp       = (M43_start_temp + (M43_end_temp - M43_start_temp) * M43_nowTime / M43_timeLimit);
  const double prob = exp((double)totalDiffScore * 100 / temp);

  if (prob > rand01()) {
    ansColumnSchedulesCount[raD][lineNum] = M43_tmpAnsCurrentCount;
    rep(i, M43_tmpAnsCurrentCount)
    {
      int num                                     = M43_tmpAnsCurrentLine[i];
      ansColumnNum[raD][num]                      = lineNum;
      ansColumnSchedules[raD][lineNum][i]         = num;
      ansColumnSchedulesPosition[raD][lineNum][i] = M43_tmpAnsCurrentLinePosition[i];
    }
    ansColumnSchedulesPosition[raD][lineNum][M43_tmpAnsCurrentCount] = M43_tmpAnsCurrentLinePosition[M43_tmpAnsCurrentCount];

    ansScore -= totalDiffScore;
    if (ansScore < real_ansScore) {
      CopyToRealAns();
      CopyToReal_M42();
    }
  }
  else {
    ;
  }
}

// 幅を不要な列に移動
int M438_yokoLineCount[MAX_LINECOUNT];
int M438_LineNumbers[MAX_LINECOUNT];
void Method4_3_8()
{
  rep(j, ansBaseLineCount) M438_LineNumbers[j] = j;
  FisherYates(M438_LineNumbers, ansBaseLineCount);
  int minLineNum     = -1;
  int minLineCount   = INT_INF;
  int taisyouLineNum = M438_LineNumbers[0];
  rep(jj, 2)
  {
    int j                 = M438_LineNumbers[jj];
    M438_yokoLineCount[j] = 0;
    srep(i, 1, d)
    {
      M438_yokoLineCount[j] += max(0, ansColumnSchedulesCount[i - 1][j] - 1) + max(0, ansColumnSchedulesCount[i][j] - 1);
      // srep(k, 1, ansColumnSchedulesCount[i - 1][j]) {
      //  rep(l, ansColumnSchedulesCount[i][j]) {
      //    if (ansColumnSchedulesPosition[i - 1][j][k] == ansColumnSchedulesPosition[i][j][l]) { M438_yokoLineCount[j] -= 2; }
      //  }
      //}
      int ite1 = 1;
      int ite2 = 1;
      while (ite1 < ansColumnSchedulesCount[i - 1][j] && ite2 < ansColumnSchedulesCount[i][j]) {
        if (ansColumnSchedulesPosition[i - 1][j][ite1] == ansColumnSchedulesPosition[i][j][ite2]) {
          M438_yokoLineCount[j] -= 2;
          ite1++;
          ite2++;
        }
        else if (ansColumnSchedulesPosition[i - 1][j][ite1] < ansColumnSchedulesPosition[i][j][ite2]) {
          ite1++;
        }
        else {
          ite2++;
        }
      }
    }
    if (M438_yokoLineCount[j] < minLineCount) {
      minLineCount = M438_yokoLineCount[j];
      minLineNum   = j;
      if (minLineCount == 0) break;
    }
  }
  if (minLineNum == taisyouLineNum) return;

  // 対象列の中身をなるべく小さくする
  // marginを計算する
  int margin = widths[taisyouLineNum] - 1;
  rep(i, d)
  {
    rep(j, ansColumnSchedulesCount[i][taisyouLineNum])
    {
      int num    = ansColumnSchedules[i][taisyouLineNum][j];
      int height = ansColumnSchedulesPosition[i][taisyouLineNum][j + 1] - ansColumnSchedulesPosition[i][taisyouLineNum][j];
      rep(k, n)
      {
        if (k == num) break;
        if (ansColumnNum[i][k] == taisyouLineNum) continue;
        if (preCalcScheduleSizes[i][k][taisyouLineNum] <= height) {
          int nextLine      = ansColumnNum[i][k];
          int nextLineIndex = -1;
          rep(l, ansColumnSchedulesCount[i][nextLine])
          {
            if (ansColumnSchedules[i][nextLine][l] == k) {
              nextLineIndex = l;
              break;
            }
          }
          if (preCalcScheduleSizes[i][num][nextLine] > ansColumnSchedulesPosition[i][nextLine][nextLineIndex + 1] - ansColumnSchedulesPosition[i][nextLine][nextLineIndex]) continue;
          swap(ansColumnSchedules[i][taisyouLineNum][j], ansColumnSchedules[i][nextLine][nextLineIndex]);
          swap(ansColumnNum[i][num], ansColumnNum[i][k]);
          break;
        }
      }
      int newNum = ansColumnSchedules[i][taisyouLineNum][j];
      // margin計算
      margin = min(margin, widths[taisyouLineNum] - ((a[i][newNum] - 1) / height + 1));
    }
  }

  if (margin == 0) return;
  int diffScore = margin * (M438_yokoLineCount[taisyouLineNum] - minLineCount);

  // ansLinePosとpreCalcScheduleSizesとwidthsを更新
  widths[taisyouLineNum] -= margin;
  widths[minLineNum] += margin;
  rep(i, d)
  {
    rep(j, ansBaseLineCount)
    {
      ansLinePos[i][j + 1] = ansLinePos[i][j] + widths[j];
    }
  }
  rep(i, d)
  {
    rep(j, n)
    {
      preCalcScheduleSizes[i][j][taisyouLineNum] = (a[i][j] - 1) / widths[taisyouLineNum] + 1;
      preCalcScheduleSizes[i][j][minLineNum]     = (a[i][j] - 1) / widths[minLineNum] + 1;
    }
  }

  ansScore -= diffScore;
  if (ansScore < real_ansScore) {
    CopyToRealAns();
    CopyToReal_M42();
  }
}

// ソート
int M439Array[MAX_N];
int M439LineNum[MAX_N];
int M439LineIndex[MAX_N];
void Method4_3_9()
{
  int raD = randxor() % d;
  int cnt = 0;
  rep(j, ansBaseLineCount)
  {
    rep(k, ansColumnSchedulesCount[raD][j])
    {
      M439LineNum[cnt]   = j;
      M439LineIndex[cnt] = k;
      M439Array[cnt]     = widths[j] * (ansColumnSchedulesPosition[raD][j][k + 1] - ansColumnSchedulesPosition[raD][j][k]) * 100 + cnt;
      cnt++;
    }
  }
  sort(M439Array, M439Array + cnt);

  rep(j, ansBaseLineCount)
  {
    rep(k, ansColumnSchedulesCount[raD][j])
    {
      ansColumnSchedules[raD][j][k] = -1;
    }
  }

  int ite = 0;
  rep(j, n)
  {
    while (true) {
      if (a[raD][j] <= M439Array[ite] / 100) {
        int posNum                                  = M439Array[ite] % 100;
        int lineNum                                 = M439LineNum[posNum];
        int lineIndex                               = M439LineIndex[posNum];
        ansColumnNum[raD][j]                        = lineNum;
        ansColumnSchedules[raD][lineNum][lineIndex] = j;
        ite++;
        break;
      }
      else {
        ite++;
      }
    }
  }
}

void Method4_3(double timeLimit)
{
  M43_startTime = GetNowTime();
  M43_timeLimit = timeLimit;
  M431Count     = 0;

  rep(i, d)
  {
    rep(j, n)
    {
      int ng = 1;
      rep(k, ansLineCount[i])
      {
        if (ans[i][j][1] == ansLinePos[i][k] && ans[i][j][3] == ansLinePos[i][k + 1]) {
          ansColumnNum[i][j] = k;
          break;
        }
      }
    }
  }

  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, ansLineCount[i])
      {
        int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
        if (width == 0) {
          preCalcScheduleSizes[i][j][k] = 1001001;
        }
        else {
          preCalcScheduleSizes[i][j][k] = (a[i][j] - 1) / width + 1;
        }
      }
    }
  }

  // 初期解作成
  rep(i, d)
  {
    rep(j, ansLineCount[i])
    {
      ansColumnSchedulesCount[i][j] = 0;
    }
  }
  rep(i, d)
  {
    drep(j, n)
    {
      int lineNum                                                         = ansColumnNum[i][j];
      ansColumnSchedules[i][lineNum][ansColumnSchedulesCount[i][lineNum]] = j;
      ansColumnSchedulesCount[i][lineNum]++;
    }
  }
  rep(i, d)
  {
    rep(j, ansLineCount[i])
    {
      ansColumnSchedulesPosition[i][j][0] = 0;
      rep(k, ansColumnSchedulesCount[i][j])
      {
        int num                                 = ansColumnSchedules[i][j][k];
        ansColumnSchedulesPosition[i][j][k + 1] = ansColumnSchedulesPosition[i][j][k] + preCalcScheduleSizes[i][num][j];
        if (k == ansColumnSchedulesCount[i][j] - 1) { ansColumnSchedulesPosition[i][j][k + 1] = w; }
      }
    }
  }

  rep(j, ansBaseLineCount)
  {
    widths[j] = ansLinePos[0][j + 1] - ansLinePos[0][j];
  }

  // realに格納
  CopyToReal_M42();

  int loopCount = 0;
  endTime       = clock();
  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      M43_nowTime = GetNowTime();
      if (M43_nowTime > timeLimit) { break; }
    }

    int ra = randxor() % 101;
    if (ra < 0) {
      if (M43_nowTime < (M43_startTime + timeLimit) / 2) continue;
      Method4_3_2();
    }
    else if (ra < 69) {  // 69
      Method4_3_1();
    }
    else if (ra < 70) {  // 70
      Method4_3_3();
    }
    else if (ra < 90) {  // 90
      Method4_3_4_2();
    }
    else if (ra < 91) {  // 91
      Method4_3_5();
    }
    else if (ra < 100) {  // 100
      Method4_3_6();
    }
    else if (ra < 100) {
      Method4_3_8();
    }
    else if (ra < 100) {
      Method4_3_7();
    }
    else if (ra < 101) {
      if (randxor() % 100 != 0) continue;
      Method4_3_9();
    }

    // if (randxor() % 212121 == 0) {
    //  CopyFromRealAns();
    //  CoptToCurrent_M42();
    //}
  }

  CopyFromRealAns();
  CoptToCurrent_M42();

  // ansScheduleLineNumからansを作成
  rep(i, d)
  {
    rep(j, ansBaseLineCount)
    {
      rep(k, ansColumnSchedulesCount[i][j])
      {
        int num        = ansColumnSchedules[i][j][k];
        ans[i][num][0] = ansColumnSchedulesPosition[i][j][k];
        ans[i][num][2] = ansColumnSchedulesPosition[i][j][k + 1];
        int posNum     = ansColumnNum[i][num];
        ans[i][num][1] = ansLinePos[i][posNum];
        ans[i][num][3] = ansLinePos[i][posNum + 1];
      }
    }
  }

  // if (mode != 0) { cout << "M431Count = " << M431Count << ", loop = " << loopCount << ", ansScore = " << ansScore << endl; }
  ansScore = CalcScore();
  CopyToRealAns();

  CopyFromRealAns();
  CoptToCurrent_M42();
  // if (mode != 0) { cout << "loop = " << loopCount << ", ansScore = " << ansScore << endl; }
}

int oshiiLineCount[MAX_D];
int oshiiLinePos[MAX_D][MAX_LINECOUNT];
int M3_alreadyUsed[MAX_D][MAX_N];
int M3_alreadyCount = 0;
int Method3_Oshii()
{
  rep(i, d)
  {
    if (oshiiLineCount[i] == -1) { return 0; }
  }
  // 作り直し
  {
    ansBaseLineCount = oshiiLineCount[0];
    rep(i, d)
    {
      ansLineCount[i] = oshiiLineCount[i];
      rep(j, oshiiLineCount[i] + 1)
      {
        ansLinePos[i][j] = oshiiLinePos[i][j];
      }
    }

    int now[MAX_LINECOUNT] = {};
    int ng                 = 0;
    drep(ii, d)
    {
      int i = daysDifficultySorted[ii];
      rep(j, ansLineCount[i])
      {
        now[j] = 0;
      }

      drep(j, n)
      {
        if (M3_alreadyUsed[i][j]) continue;
        int minAmari = INT_INF;
        int minOver  = INT_INF;
        int posNum   = -1;
        int tmpNeed  = 0;
        dsrep(k, M3_alreadyCount, ansLineCount[i])
        {
          int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
          if (width * w < a[i][j]) break;
          int need  = (a[i][j] - 1) / width + 1;
          int over  = 0;
          int amari = need * width - a[i][j];
          if (now[k] + need > w) {
            over  = (now[k] + need - w) * width;
            amari = 0;
          }

          if (over < minOver || (over == minOver && amari <= minAmari)) {
            minAmari = amari;
            minOver  = over;
            posNum   = k;
            tmpNeed  = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans[i][j][0] = now[posNum];
        ans[i][j][2] = now[posNum] + tmpNeed;
        ans[i][j][1] = ansLinePos[i][posNum];
        ans[i][j][3] = ansLinePos[i][posNum + 1];
        now[posNum] += tmpNeed;
      }
      if (ng) { break; }
    }

    if (ng) {
      CopyFromRealAns();
      return 0;
    }
  }

  // 調整
  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, ansLineCount[i])
      {
        if (ans[i][j][1] == ansLinePos[i][k] && ans[i][j][3] == ansLinePos[i][k + 1]) {
          ansColumnNum[i][j] = k;
          break;
        }
      }
    }
  }

  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, ansLineCount[i])
      {
        int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
        if (width == 0) {
          preCalcScheduleSizes[i][j][k] = 1001001;
        }
        else {
          preCalcScheduleSizes[i][j][k] = (a[i][j] - 1) / width + 1;
        }
      }
    }
  }

  int ng = 0;
  int kouhos[MAX_N];
  int kouhoCount = 0;
  int loopCount  = 0;
  rep(ii, d)
  {
    int i  = daysDifficultySorted[ii];
    int ok = 1;
    rep(j, n)
    {
      if (ans[i][j][2] > w) {
        ok = 0;
        break;
      }
    }
    if (ok) continue;
    int nowSum[MAX_LINECOUNT] = {};
    rep(j, n)
    {
      nowSum[ansColumnNum[i][j]] = max(nowSum[ansColumnNum[i][j]], ans[i][j][2]);
    }
    int ngCount = 0;
    rep(j, n)
    {
      if (nowSum[j] > w) ngCount++;
    }

    // 逆引き作成
    rep(j, ansLineCount[i])
    {
      ansColumnSchedulesCount[i][j] = 0;
    }
    rep(j, n)
    {
      int lineNum                                                         = ansColumnNum[i][j];
      ansColumnSchedules[i][lineNum][ansColumnSchedulesCount[i][lineNum]] = j;
      ansColumnSchedulesCount[i][lineNum]++;
    }

    rep(leSserafim, 10000)
    {
      // if (randxor() % 5 != 0) {
      //  int raN1     = randxor() % n;
      //  int lineNum1 = ansColumnNum[i][raN1];
      //  int raN2     = randxor() % n;
      //  if (raN1 == raN2) continue;
      //  int lineNum2 = ansColumnNum[i][raN2];
      //  if (lineNum1 == lineNum2) continue;
      //  int diffOver = max(0, nowSum[lineNum1] - w) + max(0, nowSum[lineNum2] - w);
      //  diffOver -= max(0, nowSum[lineNum1] - preCalcScheduleSizes[i][raN1][lineNum1] + preCalcScheduleSizes[i][raN2][lineNum1] - w) + max(0, nowSum[lineNum2] - preCalcScheduleSizes[i][raN2][lineNum2] + preCalcScheduleSizes[i][raN1][lineNum2] - w);
      //  if (diffOver >= 0) {
      //    // スワップ
      //    if (nowSum[lineNum1] > w) ngCount--;
      //    if (nowSum[lineNum2] > w) ngCount--;
      //    nowSum[lineNum1] = nowSum[lineNum1] - preCalcScheduleSizes[i][raN1][lineNum1] + preCalcScheduleSizes[i][raN2][lineNum1];
      //    nowSum[lineNum2] = nowSum[lineNum2] - preCalcScheduleSizes[i][raN2][lineNum2] + preCalcScheduleSizes[i][raN1][lineNum2];
      //    if (nowSum[lineNum1] > w) ngCount++;
      //    if (nowSum[lineNum2] > w) ngCount++;
      //    ansColumnNum[i][raN1] = lineNum2;
      //    ansColumnNum[i][raN2] = lineNum1;

      //    // 逆引き更新
      //    rep(k, ansColumnSchedulesCount[i][lineNum1]) {
      //      if (ansColumnSchedules[i][lineNum1][k] == raN1) {
      //        ansColumnSchedules[i][lineNum1][k] = raN2;
      //        break;
      //      }
      //    }
      //    rep(k, ansColumnSchedulesCount[i][lineNum2]) {
      //      if (ansColumnSchedules[i][lineNum2][k] == raN2) {
      //        ansColumnSchedules[i][lineNum2][k] = raN1;
      //        break;
      //      }
      //    }
      //  }

      //  if (ngCount == 0) { break; }
      //  continue;
      //}

      int raN      = randxor() % n;
      int lineNum  = ansColumnNum[i][raN];
      int nextLine = randxor() % ansBaseLineCount;
      while (nextLine == lineNum) {
        nextLine = randxor() % ansBaseLineCount;
      }
      kouhoCount = 0;
      if (nowSum[lineNum] <= w && nowSum[nextLine] <= w && preCalcScheduleSizes[i][raN][nextLine] > w) continue;
      rep(j, ansColumnSchedulesCount[i][nextLine])
      {
        if (nowSum[lineNum] <= w && nowSum[nextLine] <= w && preCalcScheduleSizes[i][ansColumnSchedules[i][nextLine][j]][lineNum] > w) continue;
        kouhos[kouhoCount] = ansColumnSchedules[i][nextLine][j];
        kouhoCount++;
      }

      if (kouhoCount == 0 || kouhoCount > 30) continue;

      if (kouhoCount <= 5) {
        drep(chaewon, (1 << kouhoCount))
        {
          if (chaewon == 0) continue;
          int nextSpace = 0;
          int needSpace = 0;
          rep(jj, kouhoCount)
          {
            if (chaewon & (1 << jj)) {
              int j = kouhos[jj];
              nextSpace += preCalcScheduleSizes[i][j][nextLine];
              needSpace += preCalcScheduleSizes[i][j][lineNum];
            }
          }
          int diffOver = max(0, nowSum[lineNum] - w) + max(0, nowSum[nextLine] - w);
          diffOver -= max(0, nowSum[lineNum] - preCalcScheduleSizes[i][raN][lineNum] + needSpace - w) + max(0, nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][raN][nextLine] - w);
          if (diffOver >= 0) {
            // スワップ
            if (nowSum[lineNum] > w) ngCount--;
            if (nowSum[nextLine] > w) ngCount--;
            nowSum[lineNum]  = nowSum[lineNum] - preCalcScheduleSizes[i][raN][lineNum] + needSpace;
            nowSum[nextLine] = nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][raN][nextLine];
            if (nowSum[lineNum] > w) ngCount++;
            if (nowSum[nextLine] > w) ngCount++;
            ansColumnNum[i][raN] = nextLine;
            rep(jj, kouhoCount)
            {
              if (chaewon & (1 << jj)) {
                int j              = kouhos[jj];
                ansColumnNum[i][j] = lineNum;
              }
            }

            // 逆引き更新
            rep(j, ansLineCount[i])
            {
              ansColumnSchedulesCount[i][j] = 0;
            }
            rep(j, n)
            {
              int lineNum                                                         = ansColumnNum[i][j];
              ansColumnSchedules[i][lineNum][ansColumnSchedulesCount[i][lineNum]] = j;
              ansColumnSchedulesCount[i][lineNum]++;
            }

            break;
          }
        }
      }
      else {
        rep(eunchae, 32)
        {
          int chaewon = randxor() % ((1 << kouhoCount) - 1) + 1;
          if (chaewon == 0) continue;
          int nextSpace = 0;
          int needSpace = 0;
          rep(jj, kouhoCount)
          {
            if (chaewon & (1 << jj)) {
              int j = kouhos[jj];
              nextSpace += preCalcScheduleSizes[i][j][nextLine];
              needSpace += preCalcScheduleSizes[i][j][lineNum];
            }
          }
          int diffOver = max(0, nowSum[lineNum] - w) + max(0, nowSum[nextLine] - w);
          diffOver -= max(0, nowSum[lineNum] - preCalcScheduleSizes[i][raN][lineNum] + needSpace - w) + max(0, nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][raN][nextLine] - w);
          if (diffOver >= 0) {
            // スワップ
            if (nowSum[lineNum] > w) ngCount--;
            if (nowSum[nextLine] > w) ngCount--;
            nowSum[lineNum]  = nowSum[lineNum] - preCalcScheduleSizes[i][raN][lineNum] + needSpace;
            nowSum[nextLine] = nowSum[nextLine] - nextSpace + preCalcScheduleSizes[i][raN][nextLine];
            if (nowSum[lineNum] > w) ngCount++;
            if (nowSum[nextLine] > w) ngCount++;
            ansColumnNum[i][raN] = nextLine;
            rep(jj, kouhoCount)
            {
              if (chaewon & (1 << jj)) {
                int j              = kouhos[jj];
                ansColumnNum[i][j] = lineNum;
              }
            }

            // 逆引き更新
            rep(j, ansLineCount[i])
            {
              ansColumnSchedulesCount[i][j] = 0;
            }
            rep(j, n)
            {
              int lineNum                                                         = ansColumnNum[i][j];
              ansColumnSchedules[i][lineNum][ansColumnSchedulesCount[i][lineNum]] = j;
              ansColumnSchedulesCount[i][lineNum]++;
            }

            break;
          }
        }
      }

      if (ngCount == 0) { break; }
    }

    if (ngCount > 0) {
      ng = 1;
      break;
    }
  }

  if (ng) {
    CopyFromRealAns();
    return 0;
  }

  // ansScheduleLineNumからansを作成
  {
    int now[MAX_LINECOUNT]     = {};
    int lastNum[MAX_LINECOUNT] = {};
    rep(i, d)
    {
      rep(j, ansLineCount[i])
      {
        now[j]     = 0;
        lastNum[j] = -1;
      }
      drep(j, n)
      {
        int posNum   = ansColumnNum[i][j];
        int need     = preCalcScheduleSizes[i][j][posNum];
        ans[i][j][0] = now[posNum];
        ans[i][j][2] = now[posNum] + need;
        ans[i][j][1] = ansLinePos[i][posNum];
        ans[i][j][3] = ansLinePos[i][posNum + 1];
        now[posNum] += need;
        lastNum[posNum] = j;
      }
      rep(j, ansLineCount[i])
      {
        if (lastNum[j] == -1) continue;
        ans[i][lastNum[j]][2] = w;
      }
    }
    ansScore = CalcScore();
    int ret  = 1;
    if (ansScore < real_ansScore) {
      // cout << "Oshii" << endl;
      CopyToRealAns();
      ret = 2;
    }
    CopyFromRealAns();
    return ret;
  }
}

void Method3_Oshii2()
{
  rep(i, d)
  {
    if (oshiiLineCount[i] == -1) { return; }
  }

  // 作り直し
  {
    ansBaseLineCount = oshiiLineCount[0];
    rep(i, d)
    {
      ansLineCount[i] = oshiiLineCount[i];
      rep(j, oshiiLineCount[i] + 1)
      {
        ansLinePos[i][j] = oshiiLinePos[i][j];
      }
    }

    int now[MAX_LINECOUNT] = {};
    int ng                 = 0;
    drep(ii, d)
    {
      int i = daysDifficultySorted[ii];
      rep(j, ansLineCount[i])
      {
        now[j] = 0;
      }

      drep(j, n)
      {
        int minAmari = INT_INF;
        int minOver  = INT_INF;
        int posNum   = -1;
        int tmpNeed  = 0;
        drep(k, ansLineCount[i])
        {
          int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
          if (width * w < a[i][j]) break;
          int need  = (a[i][j] - 1) / width + 1;
          int over  = 0;
          int amari = need * width - a[i][j];
          if (now[k] + need > w) {
            over  = (now[k] + need - w) * width;
            amari = 0;
          }

          if (over < minOver || (over == minOver && amari <= minAmari)) {
            minAmari = amari;
            minOver  = over;
            posNum   = k;
            tmpNeed  = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans[i][j][0] = now[posNum];
        ans[i][j][2] = now[posNum] + tmpNeed;
        ans[i][j][1] = ansLinePos[i][posNum];
        ans[i][j][3] = ansLinePos[i][posNum + 1];
        now[posNum] += tmpNeed;
      }
      if (ng) { break; }
    }

    if (ng) {
      CopyFromRealAns();
      return;
    }
  }

  // 調整
  {
    int ng                     = 0;
    int now[MAX_LINECOUNT]     = {};
    int lastNum[MAX_LINECOUNT] = {};
    rep(ii, d)
    {
      int i  = daysDifficultySorted[ii];
      int ok = 1;
      rep(j, n)
      {
        if (ans[i][j][2] > w) {
          ok = 0;
          break;
        }
      }
      if (ok) continue;

      ok = 0;
      srep(newjeans, 1, ansLineCount[i])
      {
        int keepPos             = ansLinePos[i][newjeans];
        ansLinePos[i][newjeans] = ansLinePos[i][newjeans + 1];

        rep(j, ansLineCount[i])
        {
          now[j]     = 0;
          lastNum[j] = -1;
        }

        int hanni = 1;
        drep(j, n)
        {
          int minAmari = INT_INF;
          int posNum   = -1;
          int tmpNeed  = 0;
          drep(k, ansLineCount[i])
          {
            int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
            if (width * w < a[i][j]) continue;
            int need = (a[i][j] - 1) / width + 1;
            if (now[k] + need > w) continue;
            int amari = need * width - a[i][j];
            if (amari <= minAmari) {
              minAmari = amari;
              posNum   = k;
              tmpNeed  = need;
            }
          }

          if (posNum == -1) {
            hanni = 0;
            break;
          }

          ans[i][j][0] = now[posNum];
          ans[i][j][2] = now[posNum] + tmpNeed;
          ans[i][j][1] = ansLinePos[i][posNum];
          ans[i][j][3] = ansLinePos[i][posNum + 1];
          now[posNum] += tmpNeed;
          lastNum[posNum] = j;
        }
        if (hanni == 1) {
          rep(j, ansLineCount[i])
          {
            if (lastNum[j] == -1) continue;
            ans[i][lastNum[j]][2] = w;
          }
          ok = 1;
          break;
        }
        else {
          ansLinePos[i][newjeans] = keepPos;
        }
      }
      if (ok == 0) {
        ng = 1;
        break;
      }
    }

    if (ng) {
      CopyFromRealAns();
      return;
    }
  }

  // ansScheduleLineNumからansを作成
  {
    ansScore = CalcScore();
    if (ansScore < real_ansScore) {
      // if (mode != 0) { cout << "OK" << endl; }
      CopyToRealAns();
    }
    CopyFromRealAns();
  }
}

int oshiiDecideMethod = 0;
int oshiiMinMax       = INT_INF;
int oshiiMinNGCount   = INT_INF;

int Method3_Normal(int loopCount)
{
  const int MIN_LINECOUNT = 2;
  ansBaseLineCount        = randxor() % n + MIN_LINECOUNT;
  if (real_ansBaseLineCount >= MIN_LINECOUNT && randxor() % 2 == 0) {
    ansBaseLineCount = real_ansBaseLineCount + 1;
  }
  else if (loopCount >= 100000 && real_ansBaseLineCount != -1) {
    ansBaseLineCount = min(ansBaseLineCount, real_ansBaseLineCount + 1);
    if (randxor() % 10 != 0) { ansBaseLineCount = max(ansBaseLineCount, real_ansBaseLineCount - 1); }
  }
  if (loopCount >= 100000 && real_ansBaseLineCount == -1) { ansBaseLineCount = randxor() % 3 + MIN_LINECOUNT; }

  if (ansBaseLineCount <= M3_alreadyCount) return 0;
  // if (real_ansBaseLineCount != -1 && ansBaseLineCount < real_ansBaseLineCount) return 0;

  if (ansBaseLineCount > lineMaxLimit) return 0;

  int startW     = ansLinePos[0][M3_alreadyCount];
  int nokoriW    = w - startW;
  int startLine  = M3_alreadyCount;
  int nokoriLine = ansBaseLineCount - M3_alreadyCount;

  {
    int ra = randxor() % 100;
    if (ra < 25) {
      int hosyouWidth = randxor() % 21;
      if (startW + nokoriLine * hosyouWidth > 800) return 0;
      int ok2 = 0;
      rep(_, 10)
      {
        ansLinePos[0][ansBaseLineCount] = (startW + nokoriW - nokoriLine * hosyouWidth);
        srep(i, startLine + 1, ansBaseLineCount)
        {
          ansLinePos[0][i] = startW + randxor() % (nokoriW - nokoriLine * hosyouWidth);
        }
        sort(ansLinePos[0] + startLine, ansLinePos[0] + ansBaseLineCount);
        int ok = 1;
        rep(i, ansBaseLineCount)
        {
          if (ansLinePos[0][i] >= ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) return 0;
      srep(i, startLine, ansBaseLineCount)
      {
        ansLinePos[0][i + 1] += hosyouWidth * (i + 1 - startLine);
      }
    }
    else if (ra < 35) {
      int ok2 = 0;
      rep(_, 10)
      {
        ansLinePos[0][ansBaseLineCount]     = w;
        int maxNeed                         = (maxA[n - 1] - 1) / w + 1;
        ansLinePos[0][ansBaseLineCount - 1] = w - maxNeed;
        if (ansLinePos[0][ansBaseLineCount - 1] <= ansLinePos[0][startLine]) break;
        int ng = 0;
        srep(i, startLine + 1, ansBaseLineCount - 1)
        {
          ansLinePos[0][i] = startW + (w - maxNeed - startW) * i / (ansBaseLineCount - 1 - startLine);
          ansLinePos[0][i] += randxor() % 31 - 15;
          if (ansLinePos[0][i] <= ansLinePos[0][i - 1]) {
            ng = 1;
            break;
          }
        }
        if (ng) continue;
        sort(ansLinePos[0] + startLine, ansLinePos[0] + ansBaseLineCount);
        int ok = 1;
        rep(i, ansBaseLineCount)
        {
          if (ansLinePos[0][i] == ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) return 0;
    }
    else if (ra < 50) {
      int ok2 = 0;
      rep(_, 10)
      {
        ansLinePos[0][ansBaseLineCount]     = w;
        int maxNeed                         = (maxA[n - 1] - 1) / w + 1;
        ansLinePos[0][ansBaseLineCount - 1] = w - maxNeed;
        if (ansLinePos[0][ansBaseLineCount - 1] <= ansLinePos[0][startLine]) break;
        int ng = 0;
        srep(i, startLine + 1, ansBaseLineCount - 1)
        {
          ansLinePos[0][i] = startW + randxor() % (nokoriW - maxNeed);
        }
        sort(ansLinePos[0] + startLine, ansLinePos[0] + ansBaseLineCount);
        int ok = 1;
        rep(i, ansBaseLineCount)
        {
          if (ansLinePos[0][i] == ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) return 0;
    }
    else {
      int ok2 = 0;
      rep(_, 10)
      {
        ansLinePos[0][ansBaseLineCount] = w;
        srep(i, 1, ansBaseLineCount)
        {
          ansLinePos[0][i] = startW + randxor() % nokoriW;
        }
        sort(ansLinePos[0] + startLine, ansLinePos[0] + ansBaseLineCount);
        int ok = 1;
        rep(i, ansBaseLineCount)
        {
          if (ansLinePos[0][i] >= ansLinePos[0][i + 1]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          ok2 = 1;
          break;
        }
      }
      if (ok2 == 0) return 0;
    }
  }

  rep(i, ansBaseLineCount)
  {
    if (ansLinePos[0][i] >= ansLinePos[0][i + 1]) { return 0; }
  }
  if (ansBaseLineCount < real_ansBaseLineCount - 1) { return 0; }

  {
    int widths[MAX_LINECOUNT] = {};
    rep(i, nokoriLine)
    {
      widths[i] = ansLinePos[0][startLine + i + 1] - ansLinePos[0][startLine + i];
    }
    sort(widths, widths + nokoriLine);
    if (widths[nokoriLine - 1] * w < maxA[n - 1]) { return 0; }
    rep(i, nokoriLine)
    {
      ansLinePos[0][startLine + i + 1] = ansLinePos[0][startLine + i] + widths[i];
    }
    // if (ansBaseLineCount > 10) {
    //  ansBaseLineCount++;
    //  int ra = randxor() % 5 + 1;
    //  rep(i, ra) {
    //    ansLinePos[0][ansBaseLineCount - i] = ansLinePos[0][ansBaseLineCount - 1 - i];
    //  }
    //  ansLinePos[0][ansBaseLineCount - ra] = (ansLinePos[0][ansBaseLineCount - (ra + 1)] + ansLinePos[0][ansBaseLineCount - (ra - 1)]) / 2;
    //}
  }
  rep(i, d)
  {
    ansLineCount[i] = ansBaseLineCount;
  }
  srep(i, 1, d)
  {
    rep(j, ansBaseLineCount + 1)
    {
      ansLinePos[i][j] = ansLinePos[0][j];
    }
  }

  int ng                     = 0;
  int now[MAX_LINECOUNT]     = {};
  int lastNum[MAX_LINECOUNT] = {};
  int tmpOshiiMax            = 0;
  int tmpOshiiNGMax          = 0;
  int numsOrder[MAX_N];
  drep(ii, d)
  {
    int i = daysDifficultySorted[ii];
    rep(j, ansLineCount[i])
    {
      now[j]     = 0;
      lastNum[j] = -1;
    }

    int tmpOshiiSum        = 0;
    int tmpOshiiNGCount    = 0;
    rep(j, n) numsOrder[j] = j;
    if (randxor() % 5 == 0) {
      int cnt = randxor() % 30;
      rep(_, cnt)
      {
        int raN = randxor() % (n - 1);
        swap(numsOrder[raN], numsOrder[raN + 1]);
      }
    }
    drep(jj, n)
    {
      int j = numsOrder[jj];
      if (M3_alreadyUsed[i][j]) { continue; }
      int minAmari = INT_INF;
      int minOver  = INT_INF;
      int posNum   = -1;
      int tmpNeed  = 0;
      dsrep(k, startLine, ansLineCount[i])
      {
        int width = ansLinePos[i][k + 1] - ansLinePos[i][k];
        if (width * w < a[i][j]) break;
        int need  = (a[i][j] - 1) / width + 1;
        int over  = 0;
        int amari = need * width - a[i][j];
        if (now[k] + need > w) {
          over  = (now[k] + need - w) * width;
          amari = 0;
        }
        if (over < minOver || (over == minOver && amari <= minAmari)) {
          minAmari = amari;
          minOver  = over;
          posNum   = k;
          tmpNeed  = need;
        }
      }

      if (posNum == -1) {
        ng = 1;
        break;
      }
      if (minOver > 0) {
        if (real_ansBaseLineCount != -1 && ansBaseLineCount <= real_ansBaseLineCount) {
          ng = 1;
          break;
        }
        if (real_ansBaseLineCount == -1 && ansBaseLineCount < MIN_LINECOUNT) {
          ng = 1;
          break;
        }
      }

      ans[i][j][0] = now[posNum];
      ans[i][j][2] = now[posNum] + tmpNeed;
      ans[i][j][1] = ansLinePos[i][posNum];
      ans[i][j][3] = ansLinePos[i][posNum + 1];
      now[posNum] += tmpNeed;
      lastNum[posNum] = j;

      if (oshiiDecideMethod == 0) {
        tmpOshiiSum += minOver;
        if (tmpOshiiSum > oshiiMinMax) {
          srep(iii, ii, d - 1)
          {
            swap(daysDifficultySorted[iii], daysDifficultySorted[iii + 1]);
          }
          ng = 1;
          break;
        }
      }
      else if (oshiiDecideMethod == 1) {
        if (minOver > 0) { tmpOshiiNGCount++; }
        if (tmpOshiiNGCount > oshiiMinNGCount) {
          ng = 1;
          break;
        }
      }
    }
    if (ng == 1) { break; }

    srep(j, startLine, ansLineCount[i])
    {
      if (lastNum[j] == -1) continue;
      ans[i][lastNum[j]][2] = w;
    }

    tmpOshiiMax   = max(tmpOshiiMax, tmpOshiiSum);
    tmpOshiiNGMax = max(tmpOshiiNGMax, tmpOshiiNGCount);
  }

  if (oshiiDecideMethod == 0) {
    if (tmpOshiiMax > 0) { ng = 2; }
  }
  else if (oshiiDecideMethod == 1) {
    if (tmpOshiiNGMax > 0) { ng = 2; }
  }

  if (ng == 0) {
    int ret  = 1;
    ansScore = CalcScoreForMethod3();
    if (ansScore < real_ansScore) {
      ret = 2;
      CopyToRealAns();
      rep(i, d)
      {
        oshiiLineCount[i] = -1;
      }
      oshiiMinMax     = INT_INF;
      oshiiMinNGCount = INT_INF;
    }
    return ret;
  }
  else if (ng == 2) {
    if (oshiiDecideMethod == 0) {
      if (tmpOshiiMax < oshiiMinMax) {
        rep(i, d)
        {
          oshiiLineCount[i] = ansLineCount[i];
          rep(j, ansLineCount[i] + 1)
          {
            oshiiLinePos[i][j] = ansLinePos[i][j];
          }
        }
        oshiiMinMax = tmpOshiiMax;
      }
    }
    else if (oshiiDecideMethod == 1) {
      if (tmpOshiiNGMax < oshiiMinNGCount) {
        rep(i, d)
        {
          oshiiLineCount[i] = ansLineCount[i];
          rep(j, ansLineCount[i] + 1)
          {
            oshiiLinePos[i][j] = ansLinePos[i][j];
          }
        }
        oshiiMinNGCount = tmpOshiiNGMax;
      }
    }
    return 0;
  }
  return 0;
}

void Method3_1(double timeLimit)
{
  keep31Count           = 0;
  int loopCount         = 0;
  ansBaseLineCount      = -1;
  real_ansBaseLineCount = -1;
  rep(i, d)
  {
    rep(j, MAX_LINECOUNT)
    {
      ansLinePos[i][j] = 0;
    }
    ansLineCount[i]      = -1;
    real_ansLineCount[i] = -1;
    oshiiLineCount[i]    = -1;
  }
  oshiiMinMax     = INT_INF;
  oshiiMinNGCount = INT_INF;

  M3_alreadyCount = 0;
  rep(i, d)
  {
    rep(j, n)
    {
      M3_alreadyUsed[i][j] = 0;
    }
  }

  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      if (GetNowTime() > timeLimit) { break; }
    }

    // already更新
    if (false) {
      M3_alreadyCount = 0;
      rep(i, d)
      {
        rep(j, n)
        {
          M3_alreadyUsed[i][j] = 0;
        }
      }

      if (real_ansBaseLineCount == -1) {
        M3_alreadyCount = randxor() % 3;
      }
      else {
        if (real_ansBaseLineCount <= 5) {
          M3_alreadyCount = randxor() % 3;
        }
        else {
          M3_alreadyCount = randxor() % (n / 2 + 1);
        }
      }

      if (M3_alreadyCount == 0) continue;
      rep(blackpink, M3_alreadyCount)
      {
        int baseWidth = randxor() % 100 + 50;
        if (M3_alreadyCount >= 5) { baseWidth = randxor() % 100 + 20; }
        if (blackpink > 0 && ansLinePos[0][blackpink - 1] + baseWidth + 20 > w) {
          M3_alreadyCount = 0;
          rep(i, d)
          {
            rep(j, n)
            {
              M3_alreadyUsed[i][j] = 0;
            }
          }
          break;
        }

        int maxWidth    = -1;
        double maxValue = 0;
        rep(lisa, 20)
        {
          int width    = baseWidth + lisa;
          int minValue = INT_INF;
          int ng       = 0;
          rep(i, d)
          {
            int ok = 0;
            rep(k, mostVariableAsCount[i][width])
            {
              if (M3_alreadyUsed[i][mostVariableAs[i][width][k][0]] == 0 && M3_alreadyUsed[i][mostVariableAs[i][width][k][1]] == 0) {
                ok = 1;
                if (mostVariableAsValue[i][width][k] < minValue) { minValue = mostVariableAsValue[i][width][k]; }
                break;
              }
            }
            if (ok == 0) {
              ng = 1;
              break;
            }
          }
          if (ng) continue;
          double tmpValue = (double)minValue / (width * w);
          if (tmpValue > maxValue) {
            maxValue = tmpValue;
            maxWidth = width;
          }
        }
        if (maxWidth == -1) {
          M3_alreadyCount = 0;
          rep(i, d)
          {
            rep(j, n)
            {
              M3_alreadyUsed[i][j] = 0;
            }
          }
          break;
        }

        ansLinePos[0][blackpink + 1] = ansLinePos[0][blackpink] + maxWidth;
        rep(i, d)
        {
          rep(k, mostVariableAsCount[i][maxWidth])
          {
            if (M3_alreadyUsed[i][mostVariableAs[i][maxWidth][k][0]] == 0 && M3_alreadyUsed[i][mostVariableAs[i][maxWidth][k][1]] == 0) {
              int j1                = mostVariableAs[i][maxWidth][k][0];
              int j2                = mostVariableAs[i][maxWidth][k][1];
              M3_alreadyUsed[i][j1] = 1;
              M3_alreadyUsed[i][j2] = 1;

              int need1     = (a[i][j1] - 1) / maxWidth + 1;
              ans[i][j1][0] = 0;
              ans[i][j1][2] = need1;
              ans[i][j1][1] = ansLinePos[0][blackpink];
              ans[i][j1][3] = ansLinePos[0][blackpink + 1];
              ans[i][j2][0] = need1;
              ans[i][j2][2] = w;
              ans[i][j2][1] = ansLinePos[0][blackpink];
              ans[i][j2][3] = ansLinePos[0][blackpink + 1];
              break;
            }
          }
        }
      }

      continue;
    }

    int ra   = randxor() % 10000;
    int isOK = 0;
    if (ra < 9990) {
      isOK = Method3_Normal(loopCount);
    }
    else if (ra < 10000) {
      // CopyToTemp();
      isOK = Method3_Oshii();
      // CopyFromTemp();
      rep(i, d)
      {
        oshiiLineCount[i] = -1;
      }
      oshiiMinMax     = INT_INF;
      oshiiMinNGCount = INT_INF;
    }
    else if (ra < 10000) {
      // Method3_Oshii2();
      // rep(i, d) {
      //  oshiiLineCount[i] = -1;
      //}
      // oshiiMinMax     = INT_INF;
      // oshiiMinNGCount = INT_INF;
    }

    if (isOK == 2) {
      // cout << "aaa" << M3_alreadyCount << ' ' << real_ansBaseLineCount << endl;
      CopyToKeep31(keep31Count % keep31KeepSize);
      keep31Count++;
    }
  }

  CopyFromRealAns();
  // if (mode != 0) { cout << "loop = " << loopCount << ", keep31Count = " << keep31Count << ", ansScore = " << ansScore << endl; }
}

// 縦線をずらす
void Method3_2(double timeLimit)
{
  // 適用可能かチェック
  {
    rep(i, d - 1)
    {
      rep(j, ansLineCount[i])
      {
        if (ansLinePos[i][j] != ansLinePos[i + 1][j]) { return; }
      }
    }
  }

  int loopCount     = 0;
  endTime           = clock();
  double start_temp = 100.1;
  double end_temp   = 0.0;

  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      if (GetNowTime() > timeLimit) { break; }
    }

    int ra  = randxor() % (ansBaseLineCount - 1) + 1;
    int ra2 = 0;
    while (ra2 == 0) {
      ra2 = randxor() % 11 - 5;
    }
    ansLinePos[0][ra] += ra2;
    if (ansLinePos[0][ra - 1] >= ansLinePos[0][ra] || ansLinePos[0][ra] >= ansLinePos[0][ra + 1]) {
      ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ra >= 2 && ansLinePos[0][ra - 1] - ansLinePos[0][ra - 2] > ansLinePos[0][ra] - ansLinePos[0][ra - 1]) {
      ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ansLinePos[0][ra] - ansLinePos[0][ra - 1] > ansLinePos[0][ra + 1] - ansLinePos[0][ra]) {
      ansLinePos[0][ra] -= ra2;
      continue;
    }
    if (ra <= ansBaseLineCount - 2 && ansLinePos[0][ra + 1] - ansLinePos[0][ra] > ansLinePos[0][ra + 2] - ansLinePos[0][ra + 1]) {
      ansLinePos[0][ra] -= ra2;
      continue;
    }

    int ng                     = 0;
    int now[MAX_LINECOUNT]     = {};
    int lastNum[MAX_LINECOUNT] = {};
    rep(i, d)
    {
      rep(j, ansBaseLineCount)
      {
        now[j]     = 0;
        lastNum[j] = -1;
      }

      drep(j, n)
      {
        int minAmari = INT_INF;
        int posNum   = -1;
        int tmpNeed  = 0;
        drep(k, ansBaseLineCount)
        {
          int width = ansLinePos[0][k + 1] - ansLinePos[0][k];
          if (width * w < a[i][j]) break;
          int need = (a[i][j] - 1) / width + 1;
          if (now[k] + need > w) continue;
          int amari = need * width - a[i][j];
          if (amari <= minAmari) {
            minAmari = amari;
            posNum   = k;
            tmpNeed  = need;
          }
        }

        if (posNum == -1) {
          ng = 1;
          break;
        }

        ans[i][j][0] = now[posNum];
        ans[i][j][2] = now[posNum] + tmpNeed;
        ans[i][j][1] = ansLinePos[0][posNum];
        ans[i][j][3] = ansLinePos[0][posNum + 1];
        now[posNum] += tmpNeed;
        lastNum[posNum] = j;
      }
      if (ng) { break; }

      rep(j, ansBaseLineCount)
      {
        if (lastNum[j] == -1) continue;
        ans[i][lastNum[j]][2] = w;
      }
    }

    if (ng) {
      ansLinePos[0][ra] -= ra2;
      continue;
    }

    int beforeScore = ansScore;
    ansScore        = CalcScoreForMethod3();
    int diffScore   = beforeScore - ansScore;

    double temp       = (start_temp + (end_temp - start_temp) * GetNowTime() / timeLimit);
    const double prob = exp((double)diffScore / temp);

    if (prob > rand01()) {
      srep(i, 1, d)
      {
        ansLinePos[i][ra] = ansLinePos[0][ra];
      }
      if (ansScore <= real_ansScore) { CopyToRealAns(); }
    }
    else {
      ansScore = beforeScore;
      ansLinePos[0][ra] -= ra2;
    }
  }

  CopyFromRealAns();

  if (mode != 0) {
    // cout << "loop = " << loopCount << ", ansScore = " << ansScore << endl;
  }
}

void Method6_ColumnShuffle(double timeLimit)
{
  int widths[MAX_LINECOUNT] = {};
  rep(j, ansBaseLineCount)
  {
    widths[j] = ansLinePos[0][j + 1] - ansLinePos[0][j];
  }

  rep(i, d)
  {
    rep(j, n)
    {
      rep(k, ansBaseLineCount)
      {
        if (ans[i][j][1] == ansLinePos[i][k] && ans[i][j][3] == ansLinePos[i][k + 1]) {
          ansColumnNum[i][j] = k;
          break;
        }
      }
    }
  }

  int v[MAX_LINECOUNT] = {};
  rep(i, ansBaseLineCount)
  {
    v[i] = i;
  }
  rep(ningning, 100)
  {
    FisherYates(v, ansBaseLineCount);
    int nextLinePos[MAX_LINECOUNT] = {};
    int nextArgPos[MAX_LINECOUNT]  = {};
    rep(i, ansBaseLineCount)
    {
      nextLinePos[i + 1] = nextLinePos[i] + widths[v[i]];
      nextArgPos[v[i]]   = i;
    }
    rep(i, d)
    {
      rep(j, n)
      {
        int nextLine = nextArgPos[ansColumnNum[i][j]];
        ans[i][j][1] = nextLinePos[nextLine];
        ans[i][j][3] = nextLinePos[nextLine + 1];
      }
    }
    ansScore = CalcScore();
    if (ansScore < real_ansScore) { CopyToRealAns(); }
    if (ningning % 10 == 9 && GetNowTime() > timeLimit) break;
  }
  CopyFromRealAns();
}

void Method7()
{
  CopyFromRealAns();
  CoptToCurrent_M42();

  rep(i, ansBaseLineCount)
  {
    widths[i] = ansLinePos[0][i + 1] - ansLinePos[0][i];
  }

  int yokoLineCount[MAX_LINECOUNT] = {};
  rep(j, ansBaseLineCount)
  {
    srep(i, 1, d)
    {
      yokoLineCount[j] += max(0, ansColumnSchedulesCount[i - 1][j] - 1) + max(0, ansColumnSchedulesCount[i][j] - 1);
      int ite1 = 1;
      int ite2 = 1;
      while (ite1 < ansColumnSchedulesCount[i - 1][j] && ite2 < ansColumnSchedulesCount[i][j]) {
        if (ansColumnSchedulesPosition[i - 1][j][ite1] == ansColumnSchedulesPosition[i][j][ite2]) {
          yokoLineCount[j] -= 2;
          ite1++;
          ite2++;
        }
        else if (ansColumnSchedulesPosition[i - 1][j][ite1] < ansColumnSchedulesPosition[i][j][ite2]) {
          ite1++;
        }
        else {
          ite2++;
        }
      }
    }
  }

  rep(line1, ansBaseLineCount)
  {
    rep(line2, ansBaseLineCount)
    {
      if (line1 == line2) { continue; }
      if (yokoLineCount[line1] <= yokoLineCount[line2]) continue;
      int margin = widths[line1] - 1;
      rep(i, d)
      {
        rep(j, ansColumnSchedulesCount[i][line1])
        {
          int num    = ansColumnSchedules[i][line1][j];
          int height = ansColumnSchedulesPosition[i][line1][j + 1] - ansColumnSchedulesPosition[i][line1][j];
          rep(k, n)
          {
            if (k == num) break;
            if (ansColumnNum[i][k] == line1) continue;
            if (preCalcScheduleSizes[i][k][line1] <= height) {
              int nextLine      = ansColumnNum[i][k];
              int nextLineIndex = -1;
              rep(l, ansColumnSchedulesCount[i][nextLine])
              {
                if (ansColumnSchedules[i][nextLine][l] == k) {
                  nextLineIndex = l;
                  break;
                }
              }
              if (preCalcScheduleSizes[i][num][nextLine] > ansColumnSchedulesPosition[i][nextLine][nextLineIndex + 1] - ansColumnSchedulesPosition[i][nextLine][nextLineIndex]) continue;
              swap(ansColumnSchedules[i][line1][j], ansColumnSchedules[i][nextLine][nextLineIndex]);
              swap(ansColumnNum[i][num], ansColumnNum[i][k]);
              break;
            }
          }
          int newNum = ansColumnSchedules[i][line1][j];
          // margin計算
          margin = min(margin, widths[line1] - ((a[i][newNum] - 1) / height + 1));
        }
      }

      if (margin == 0) continue;
      int diffScore = margin * (yokoLineCount[line1] - yokoLineCount[line2]);

      // ansLinePosとpreCalcScheduleSizesとwidthsを更新
      widths[line1] -= margin;
      widths[line2] += margin;
      rep(i, d)
      {
        rep(j, ansBaseLineCount)
        {
          ansLinePos[i][j + 1] = ansLinePos[i][j] + widths[j];
        }
      }
      rep(i, d)
      {
        rep(j, n)
        {
          preCalcScheduleSizes[i][j][line1] = (a[i][j] - 1) / widths[line1] + 1;
          preCalcScheduleSizes[i][j][line2] = (a[i][j] - 1) / widths[line2] + 1;
        }
      }

      ansScore -= diffScore;
    }
  }

  if (ansScore < real_ansScore) {
    CopyToRealAns();
    CopyToReal_M42();
  }
}

int M8ans[1010][MAX_D][MAX_N][4];
int M8ansOK[1010][MAX_D];
int M8ansScore[1010][MAX_D];
int M8ansLineCount[1010];
int M8ansLinePos[1010][MAX_LINECOUNT];
int M8ansansCount;
int M8ansNum[MAX_D];
int M8ansansCountEachDay[MAX_D];
int M8ansansNumEachDay[MAX_D][1010];

void CopyM8ToAns()
{
  rep(i, d)
  {
    int num = M8ansNum[i];
    rep(j, n)
    {
      rep(k, 4)
      {
        ans[i][j][k] = M8ans[num][i][j][k];
      }
    }
  }
}

void Method8(double timeLimit)
{
  const int TMP_NUM = 1000;
  const int NG_NUM  = 1005;
  M8ansansCount     = 0;
  rep(i, d)
  {
    M8ansNum[i] = -1;
  }

  rep(i, d)
  {
    M8ansNum[i]           = NG_NUM;
    M8ansScore[NG_NUM][i] = 1001001;
    rep(j, n)
    {
      rep(k, 4)
      {
        M8ans[NG_NUM][i][j][k] = ans[i][j][k];
      }
    }
  }

  double M8_startTime = GetNowTime();
  double M8_NowTime   = GetNowTime();
  int loopCount       = 0;
  while (M8ansansCount < 1000) {
    loopCount++;
    if (loopCount % 25 == 0) {
      M8_NowTime = GetNowTime();
      if (M8_NowTime > M8_startTime + (timeLimit - M8_startTime) / 2) { break; }
    }

    // 縦線作成
    ansBaseLineCount = randxor() % 5 + 2;
    int ok2          = 0;
    rep(_2, 10)
    {
      ansLinePos[0][ansBaseLineCount] = w;
      srep(i, 1, ansBaseLineCount)
      {
        ansLinePos[0][i] = randxor() % w;
      }
      sort(ansLinePos[0], ansLinePos[0] + ansBaseLineCount);
      int ok = 1;
      rep(i, ansBaseLineCount)
      {
        if (ansLinePos[0][i] >= ansLinePos[0][i + 1]) {
          ok = 0;
          break;
        }
      }
      if (ok) {
        ok2 = 1;
        break;
      }
    }
    if (ok2 == 0) continue;

    {
      int widths[MAX_LINECOUNT] = {};
      rep(i, ansBaseLineCount)
      {
        widths[i] = ansLinePos[0][i + 1] - ansLinePos[0][i];
      }
      sort(widths, widths + ansBaseLineCount);
      rep(i, ansBaseLineCount)
      {
        ansLinePos[0][i + 1] = ansLinePos[0][i] + widths[i];
      }
    }

    M8ansLineCount[TMP_NUM] = ansBaseLineCount;
    srep(i, 0, ansBaseLineCount + 1)
    {
      M8ansLinePos[TMP_NUM][i] = ansLinePos[0][i];
    }

    int okCount                = 0;
    int now[MAX_LINECOUNT]     = {};
    int lastNum[MAX_LINECOUNT] = {};
    int numsOrder[MAX_N];
    drep(ii, d)
    {
      int i = daysDifficultySorted[ii];
      rep(j, ansBaseLineCount)
      {
        now[j]     = 0;
        lastNum[j] = -1;
      }

      rep(j, n) numsOrder[j] = j;
      if (randxor() % 5 == 0) {
        int cnt = randxor() % 30;
        rep(_, cnt)
        {
          int raN = randxor() % (n - 1);
          swap(numsOrder[raN], numsOrder[raN + 1]);
        }
      }
      int ok = 1;
      drep(jj, n)
      {
        int j        = numsOrder[jj];
        int minAmari = INT_INF;
        int posNum   = -1;
        int tmpNeed  = 0;
        dsrep(k, 0, ansBaseLineCount)
        {
          int width = ansLinePos[0][k + 1] - ansLinePos[0][k];
          if (width * w < a[i][j]) break;
          int need  = (a[i][j] - 1) / width + 1;
          int amari = need * width - a[i][j];
          if (now[k] + need > w) { continue; }
          if (amari < minAmari) {
            minAmari = amari;
            posNum   = k;
            tmpNeed  = need;
          }
        }

        if (posNum == -1) {
          ok = 0;
          break;
        }

        M8ans[TMP_NUM][i][j][0] = now[posNum];
        M8ans[TMP_NUM][i][j][2] = now[posNum] + tmpNeed;
        M8ans[TMP_NUM][i][j][1] = ansLinePos[0][posNum];
        M8ans[TMP_NUM][i][j][3] = ansLinePos[0][posNum + 1];
        now[posNum] += tmpNeed;
        lastNum[posNum] = j;
      }
      if (ok == 1) {
        okCount++;
        M8ansOK[TMP_NUM][i] = 1;
        srep(j, 0, ansBaseLineCount)
        {
          if (lastNum[j] == -1) continue;
          M8ans[TMP_NUM][i][lastNum[j]][2] = w;
        }
      }
      else {
        M8ansOK[TMP_NUM][i] = 0;
      }
    }

    if (okCount >= 2) {
      M8ansLineCount[M8ansansCount] = ansBaseLineCount;
      rep(i, ansBaseLineCount + 1)
      {
        M8ansLinePos[M8ansansCount][i] = M8ansLinePos[TMP_NUM][i];
      }
      rep(i, d)
      {
        M8ansOK[M8ansansCount][i] = M8ansOK[TMP_NUM][i];
        if (M8ansOK[M8ansansCount][i]) {
          int score = 0;
          rep(j, n)
          {
            rep(k, 4)
            {
              M8ans[M8ansansCount][i][j][k] = M8ans[TMP_NUM][i][j][k];
            }
            if (M8ans[M8ansansCount][i][j][0] != 0 && M8ans[M8ansansCount][i][j][0] != w) { score += M8ans[M8ansansCount][i][j][3] - M8ans[M8ansansCount][i][j][1]; }
            if (M8ans[M8ansansCount][i][j][2] != 0 && M8ans[M8ansansCount][i][j][2] != w) { score += M8ans[M8ansansCount][i][j][3] - M8ans[M8ansansCount][i][j][1]; }
          }
          score /= 2;
          score += (M8ansLineCount[M8ansansCount] - 1) * w;
          M8ansScore[M8ansansCount][i] = score;
          if (score < M8ansScore[M8ansNum[i]][i]) { M8ansNum[i] = M8ansansCount; }
        }
      }
      M8ansansCount++;
    }
  }

  CopyM8ToAns();
  ansScore = CalcScore();
  if (ansScore < real_ansScore) { CopyToRealAns(); }
  rep(i, d)
  {
    M8ansansCountEachDay[i] = 0;
    rep(j, M8ansansCount)
    {
      if (M8ansOK[j][i]) {
        M8ansansNumEachDay[i][M8ansansCountEachDay[i]] = j;
        M8ansansCountEachDay[i]++;
      }
    }
  }

  // cout << M8ansansCount << endl;
  // rep(i, d) {
  //  int cnt = 0;
  //  rep(j, M8ansansCount) {
  //    if (M8ansOK[j][i]) cnt++;
  //  }
  //  cout << cnt << ' ';
  //}
  // cout << endl;

  double M8_startTime2 = GetNowTime();
  loopCount            = 0;
  M8_NowTime           = GetNowTime();
  while (true) {
    loopCount++;
    if (loopCount % 100 == 0) {
      M8_NowTime = GetNowTime();
      if (M8_NowTime > timeLimit) { break; }
    }

    if (loopCount % 12345 == 0) {
      CopyM8ToAns();
      ansScore = CalcScore();
      if (ansScore < real_ansScore) { CopyToRealAns(); }
    }

    int raD = randxor() % d;
    if (M8ansansCountEachDay[raD] == 0) continue;
    int raNum = M8ansansNumEachDay[raD][randxor() % M8ansansCountEachDay[raD]];
    if (raNum == M8ansNum[raD]) continue;

    int beforeNum = M8ansNum[raD];
    int diffScore = (M8ansScore[beforeNum][raD] - M8ansScore[raNum][raD]) * 2;
    if (raD == 0 || raD == d - 1) { diffScore = (M8ansScore[beforeNum][raD] - M8ansScore[raNum][raD]); }

    if (raD > 0) {
      int beforeDayNum = M8ansNum[raD - 1];
      if (beforeNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[beforeNum] && ite2 < M8ansLineCount[beforeDayNum]) {
          if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            diffScore -= w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
      if (raNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[raNum] && ite2 < M8ansLineCount[beforeDayNum]) {
          if (M8ansLinePos[raNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            diffScore += w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[raNum][ite1] == M8ansLinePos[beforeDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
    }
    if (raD < d - 1) {
      int afterDayNum = M8ansNum[raD + 1];
      if (beforeNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[beforeNum] && ite2 < M8ansLineCount[afterDayNum]) {
          if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            diffScore -= w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[beforeNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
      if (raNum < 1000) {
        int ite1 = 1;
        int ite2 = 1;
        while (ite1 < M8ansLineCount[raNum] && ite2 < M8ansLineCount[afterDayNum]) {
          if (M8ansLinePos[raNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            diffScore += w * 2;
            ite1++;
            ite2++;
          }
          else if (M8ansLinePos[raNum][ite1] == M8ansLinePos[afterDayNum][ite2]) {
            ite1++;
          }
          else {
            ite2++;
          }
        }
      }
    }

    double timeRatio  = (M8_NowTime - M8_startTime2) / (timeLimit - M8_startTime2);
    double temp       = (M43_start_temp + (M43_end_temp - M43_start_temp) * timeRatio);
    const double prob = exp((double)diffScore * 100 / temp);

    if (prob > rand01()) { M8ansNum[raD] = raNum; }
  }

  CopyM8ToAns();
  ansScore = CalcScore();
  if (ansScore < real_ansScore) { CopyToRealAns(); }

  // cout << loopCount << endl;
}

int isFind = 0;
void Method4(int setCount)
{
  CopyToRealRealAns();

  CopyFromRealAns();
  CopyToTemp();

  double outerTL     = TL;
  double m4StartTime = GetNowTime();
  rep(twice, setCount)
  {
    CopyFromTemp();
    CopyToRealAns();
    double nowTime = GetNowTime();
    double innerTL = (outerTL - m4StartTime) * (twice + 1) / setCount + m4StartTime;
    Method3_1(nowTime + (innerTL - nowTime) * 0.5);

    CopyFromRealAns();

    if (ansBaseLineCount == -1) {
      isFind = 0;
    }
    else {
      isFind = 1;
    }

    if (ansScore == 1) { return; }

    if (ansBaseLineCount == -1) {
      Method8(nowTime + (innerTL - nowTime) * 0.98);
      if (real_ansScore < real_real_ansScore) { CopyToRealRealAns(); }
    }
    else {
      Method3_2(nowTime + (innerTL - nowTime) * 0.6);

      CopyFromRealAns();
      Method4_3(nowTime + (innerTL - nowTime) * 0.98);

      // Method7(); // バグってる

      // 列位置シャッフル
      // cout << GetNowTime() << endl;
      CopyFromRealAns();
      Method6_ColumnShuffle(nowTime + (innerTL - nowTime) * 1.0);
      // cout << GetNowTime() << endl;

      if (real_ansScore < real_real_ansScore) { CopyToRealRealAns(); }
    }
  }

  CopyFromRealRealAns();
  CopyFromRealAns();
}

void Method5()
{
  double TL31 = TL * 0.5;
  Method3_1(TL31);

  if (ansBaseLineCount == -1) return;

  CopyToRealRealAns();

  keep31Count = min(keep31KeepSize, keep31Count);
  double TL2  = TL * 0.7;
  rep(tzuyu, keep31Count)
  {
    double nowTime   = GetNowTime();
    double innerTL   = (TL2 - nowTime) * (tzuyu + 1) / keep31KeepSize + nowTime;
    double innerTL32 = (innerTL - nowTime) / 5.0 + nowTime;
    CopyFromKeep31(tzuyu);
    CopyToRealAns();
    Method3_2(innerTL32);
    CopyFromRealAns();
    Method4_3(innerTL);
    if (real_ansScore < real_real_ansScore) { CopyToRealRealAns(); }
  }

  CopyFromRealRealAns();

  // double TL32 = TL * 0.75;
  // CopyFromRealAns();
  // Method3_2(TL32);

  CopyFromRealAns();
  Method4_3(TL);
  if (real_ansScore < real_real_ansScore) { CopyToRealRealAns(); }

  CopyFromRealRealAns();
  CopyFromRealAns();
}

ll Solve(int probNum)
{
  startTime = clock();
  endTime   = clock();

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

  MethodPerfect();

  Method4(1);
  // Method5();

  // if (ee < 0.1) {
  //  cout << 0 << endl;
  //  return 0;
  //}

  // 解答を出力
  CopyFromRealAns();
  Output(ofs);

  if (ofs.is_open()) { ofs.close(); }

  if (mode != 0) {
    bool isNGAns = IsNGAns();
    // if (isNGAns) { cout << "!!!IsNGANS" << endl; }
  }

  ll score = 0;
  if (mode != 0) { score = CalcScore(); }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    rep(_, 1)
    {
      vector<ll> scores;
      ll sum = 0;
      srep(i, 0, 5000)
      {
        // lineMaxLimit = i;
        ll score = Solve(i);
        sum += score;
        // cout << "num = " << i << ", ";
        // cout << "score = " << score << ", ";
        // cout << "sum = " << sum << endl;
        int maxASum = 0;
        rep(j, n)
        {
          maxASum += maxA[j];
        }
        if (true || (isFind == 0 && ansScore > 1)) cout << i << ", " << d << ", " << n << ", " << 1 - ee << ", " << maxASum << ", " << isFind << ", " << ansBaseLineCount << ", " << score << ", " << sum << ' ' << GetNowTime() << endl;
        scores.push_back(score);
      }
      for (auto score : scores) {
        cout << score << endl;
      }
    }
  }

  return 0;
}
