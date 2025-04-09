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

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() { return (Rand() + 0.5) * (1.0 / UINT_MAX); }

// l以上r未満の実数をとる乱数
static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);


const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;


const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char dc[4] = { 'U','L','D','R' };

double TL = 1.9;
int mode;
std::chrono::steady_clock::time_point startTimeClock; // 時間計測用


void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}


double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}


const int MAX_N = 30;

const int n = 20;

int a[n + 2][n + 2];
int init_a[n + 2][n + 2];

int ansScore;
vector<P> ans;

int best_ansScore;
vector<P> best_ans;

void CopyToBest()
{
  best_ansScore = ansScore;
  best_ans = ans;
}

void CopyToAns()
{
  ansScore = best_ansScore;
  ans = best_ans;
}


bool IsNG(int x, int y)
{
  if (x < 1 || n < x || y < 1 || n < y)return true;
  return false;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  ans.clear();
  CopyToBest();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn;
    cin >> nn;
    rep(i, n)
    {
      string s;
      cin >> s;
      rep(j, n)
      {
        a[i + 1][j + 1] = 0;
        if (s[j] == 'x') {
          a[i + 1][j + 1] = 1;
        }
        else if (s[j] == 'o') {
          a[i + 1][j + 1] = 2;
        }
      }
    }
  }
  else {
    // ファイル入力
    int nn;
    ifs >> nn;
    rep(i, n)
    {
      string s;
      ifs >> s;
      rep(j, n)
      {
        a[i + 1][j + 1] = 0;
        if (s[j] == 'x') {
          a[i + 1][j + 1] = 1;
        }
        else if (s[j] == 'o') {
          a[i + 1][j + 1] = 2;
        }
      }
    }
  }

  rep(i, n + 2)
  {
    rep(j, n + 2)
    {
      init_a[i][j] = a[i][j];
    }
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
ll CalcScore()
{
  ll res = 8 * n * n - ans.size();
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    for (auto p : ans) {
      cout << dc[p.first] << ' ' << p.second - 1 << endl;
    }
  }
  else {
    // ファイル出力
    for (auto p : ans) {
      ofs << dc[p.first] << ' ' << p.second - 1 << endl;
    }
  }
}

void MoveU(int j, bool isClean)
{
  rep(i, n)
  {
    a[i][j] = a[i + 1][j];
  }
  if (isClean) {
    a[n][j] = 0;
  }
  else {
    a[n][j] = a[n + 1][j];
  }
  a[n + 1][j] = 0;
}

void MoveL(int i, bool isClean)
{
  rep(j, n)
  {
    a[i][j] = a[i][j + 1];
  }
  if (isClean) {
    a[i][n] = 0;
  }
  else {
    a[i][n] = a[i][n + 1];
  }
  a[i][n + 1] = 0;
}

void MoveD(int j, bool isClean)
{
  for (int i = n + 1; i >= 2; i--) {
    a[i][j] = a[i - 1][j];
  }
  if (isClean) {
    a[1][j] = 0;
  }
  else {
    a[1][j] = a[0][j];
  }
  a[0][j] = 0;
}

void MoveR(int i, bool isClean)
{
  for (int j = n + 1; j >= 2; j--) {
    a[i][j] = a[i][j - 1];
  }
  if (isClean) {
    a[i][1] = 0;
  }
  else {
    a[i][1] = a[i][0];
  }
  a[i][0] = 0;
}

int Method2_CalcScore(int dir, int num)
{
  int score = 0;
  queue<P> que;
  int f[n + 2][n + 2];
  srep(i, 1, n + 1)
  {
    srep(j, 1, n + 1)
    {
      if (a[i][j] == 1) {
        while (que.size()) {
          que.pop();
        }
        rep(k, n + 2)
        {
          rep(l, n + 2)
          {
            f[k][l] = 99999;
          }
        }
        int dist = 99999;
        f[i][j] = 0;
        que.push(P(i, j));
        while (que.size()) {
          int x = que.front().first;
          int y = que.front().second;
          que.pop();
          if (x == 1 || x == n || y == 1 || y == n) {
            dist = f[x][y] + 1;
            break;
          }
          rep(k, 4)
          {
            int nx = x + dx[k];
            int ny = y + dy[k];
            if (IsNG(nx, ny))continue;
            if (a[nx][ny] == 2)continue;
            if (f[nx][ny] > f[x][y] + 1) {
              f[nx][ny] = f[x][y] + 1;
              que.push(P(nx, ny));
            }
          }
        }
        score += dist + 100;
        if (dir == 0 || dir == 2) {
          if (j == num)score--;
        }
        else {
          if (i == num)score--;
        }
      }
    }
  }

  return score;
}

// BFS
void Method2()
{
  int outerLoop = 0;

  while (outerLoop == 0) {

    if (GetNowTime() > TL) {
      break;
    }
    outerLoop++;

    ans.clear();
    int cnt = 0;
    int loop = 0;

    int befDir = -1;
    int befNum = -1;

    while (cnt < n * 2 && loop < 500) {


      loop++;

      int dir = -1;
      int num = -1;
      int cost = 99999;

      // U
      srep(j, 1, n + 1)
      {
        if (a[1][j] == 2) {
          continue;
        }
        if (befDir == 2 && befNum == j)continue;
        MoveU(j, true);
        int score = Method2_CalcScore(0, j);
        if (score < cost) {
          dir = 0;
          num = j;
          cost = score;
        }
        MoveD(j, false);
      }

      // L
      srep(i, 1, n + 1)
      {
        if (a[i][1] == 2) {
          continue;
        }
        if (befDir == 3 && befNum == i)continue;
        MoveL(i, true);
        int score = Method2_CalcScore(1, i);
        if (score < cost) {
          dir = 1;
          num = i;
          cost = score;
        }
        MoveR(i, false);
      }

      // D
      srep(j, 1, n + 1)
      {
        if (a[n][j] == 2) {
          continue;
        }
        if (befDir == 0 && befNum == j)continue;
        MoveD(j, true);
        int score = Method2_CalcScore(2, j);
        if (score < cost) {
          dir = 2;
          num = j;
          cost = score;
        }
        MoveU(j, false);
      }

      // R
      srep(i, 1, n + 1)
      {
        if (a[i][n] == 2) {
          continue;
        }
        if (befDir == 1 && befNum == i)continue;
        MoveR(i, true);
        int score = Method2_CalcScore(3, i);
        if (score < cost) {
          dir = 3;
          num = i;
          cost = score;
        }
        MoveL(i, false);
      }

      ans.push_back(P(dir, num));
      if (dir == 0) {
        if (a[1][num] == 1)cnt++;
        MoveU(num, true);
      }
      else if (dir == 1) {
        if (a[num][1] == 1)cnt++;
        MoveL(num, true);
      }
      else if (dir == 2) {
        if (a[n][num] == 1)cnt++;
        MoveD(num, true);
      }
      else if (dir == 3) {
        if (a[num][n] == 1)cnt++;
        MoveR(num, true);
      }

      befDir = dir;
      befNum = num;
    }

    ansScore = CalcScore();

    if (ansScore > best_ansScore) {
      CopyToBest();
    }
  }

  CopyToAns();
}

void InitA()
{
  rep(i, n + 2)
  {
    rep(j, n + 2)
    {
      a[i][j] = init_a[i][j];
    }
  }
}

int Sim()
{
  InitA();
  int cnt = 0;
  rep(i, ans.size())
  {
    int dir = ans[i].first;
    int num = ans[i].second;
    if (dir == 0) {
      if (a[1][num] == 1) {
        cnt++;
      }
      else if (a[1][num] == 2) {
        return -1;
      }
      MoveU(num, true);
    }
    else if (dir == 1) {
      if (a[num][1] == 1) {
        cnt++;
      }
      else if (a[num][1] == 2) {
        return -1;
      }
      MoveL(num, true);
    }
    else if (dir == 2) {
      if (a[n][num] == 1) {
        cnt++;
      }
      else if (a[n][num] == 2) {
        return -1;
      }
      MoveD(num, true);
    }
    else if (dir == 3) {
      if (a[num][n] == 1) {
        cnt++;
      }
      else if (a[num][n] == 2) {
        return -1;
      }
      MoveR(num, true);
    }
    if (cnt == n * 2) {
      return i + 1;
    }
  }

  if (cnt < n * 2) {
    return -1;
  }

  return ans.size();
}

int Sim2(vector<P>& ans2, bool oniCheck)
{
  ans2.clear();
  InitA();
  int cnt = 0;
  rep(i, ans.size())
  {
    int dir = ans[i].first;
    int num = ans[i].second;
    if (dir == 0) {
      if (a[1][num] == 1) {
        cnt++;
      }
      else if (a[1][num] == 2) {
        return -1;
      }

      if (oniCheck) {
        int oni = 0;
        srep(i, 1, n + 1)
        {
          if (a[i][num] == 1) {
            oni = 1;
            break;
          }
        }
        if (oni == 0)continue;
      }

      MoveU(num, true);
      ans2.push_back(ans[i]);
    }
    else if (dir == 1) {
      if (a[num][1] == 1) {
        cnt++;
      }
      else if (a[num][1] == 2) {
        return -1;
      }

      if (oniCheck) {
        int oni = 0;
        srep(j, 1, n + 1)
        {
          if (a[num][j] == 1) {
            oni = 1;
            break;
          }
        }
        if (oni == 0)continue;
      }

      MoveL(num, true);
      ans2.push_back(ans[i]);
    }
    else if (dir == 2) {
      if (a[n][num] == 1) {
        cnt++;
      }
      else if (a[n][num] == 2) {
        return -1;
      }

      if (oniCheck) {
        int oni = 0;
        srep(i, 1, n + 1)
        {
          if (a[i][num] == 1) {
            oni = 1;
            break;
          }
        }
        if (oni == 0)continue;
      }

      MoveD(num, true);
      ans2.push_back(ans[i]);
    }
    else if (dir == 3) {
      if (a[num][n] == 1) {
        cnt++;
      }
      else if (a[num][n] == 2) {
        return -1;
      }

      if (oniCheck) {
        int oni = 0;
        srep(j, 1, n + 1)
        {
          if (a[num][j] == 1) {
            oni = 1;
            break;
          }
        }
        if (oni == 0)continue;
      }

      MoveR(num, true);
      ans2.push_back(ans[i]);
    }
    if (cnt == n * 2) {
      return ans2.size();
    }
  }

  srep(i, 1, n + 1)
  {
    //if (ans2.size() >= ansScore + 5)break;
    srep(j, 1, n + 1)
    {
      //if (ans2.size() >= ansScore + 5)break;

      if (a[i][j] == 1) {
        int dir = -1;
        int cost = 999;

        // U
        int ok = 1;
        srep(k, 1, i)
        {
          if (a[k][j] != 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          if (i < cost) {
            dir = 0;
            cost = i;
          }
        }

        // L
        ok = 1;
        srep(k, 1, j)
        {
          if (a[i][k] != 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          if (j < cost) {
            dir = 1;
            cost = j;
          }
        }

        // D
        ok = 1;
        srep(k, i + 1, n + 1)
        {
          if (a[k][j] != 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          if (n + 1 - i < cost) {
            dir = 2;
            cost = n + 1 - i;
          }
        }

        // R
        ok = 1;
        srep(k, j + 1, n + 1)
        {
          if (a[i][k] != 0) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          if (n + 1 - j < cost) {
            dir = 3;
            cost = n + 1 - j;
          }
        }

        if (dir == -1) {
          return -1;
        }

        cnt++;

        if (dir == 0) {
          rep(k, cost)
          {
            ans2.push_back(P(dir, j));
            MoveU(j, true);
          }
        }
        else if (dir == 1) {
          rep(k, cost)
          {
            ans2.push_back(P(dir, i));
            MoveL(i, true);
          }
        }
        else if (dir == 2) {
          rep(k, cost)
          {
            ans2.push_back(P(dir, j));
            MoveD(j, true);
          }
        }
        else if (dir == 3) {
          rep(k, cost)
          {
            ans2.push_back(P(dir, i));
            MoveR(i, true);
          }
        }
      }
    }
  }

  if (cnt < n * 2) {
    return -1;
  }

  return ans2.size();
}

struct Haiparas
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition1;
  double SimPartition;
};

void Mountain(Haiparas haiparas)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = haiparas.StartTemp;
  const double END_TEMP = haiparas.EndTemp;

  int loop = 0;
  int loop2 = 0;
  vector<P> ans2;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) break;
    }

    if (ans.size() >= best_ans.size() + 4) {
      CopyToAns();
    }

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    int raMode = Rand() % 100;
    if (raMode < haiparas.Partition1) {
      // swap
      int ra1 = Rand() % ans.size();
      int ra2 = Rand() % ans.size();
      while (ans[ra1] == ans[ra2]) {
        ra2 = Rand() % ans.size();
      }
      swap(ans[ra1], ans[ra2]);

      if (Rand01() < haiparas.SimPartition) {
        int res = Sim();
        if (res == -1) {
          swap(ans[ra1], ans[ra2]);
        }
        else {
          loop2++;
          while (ans.size() > res) {
            ans.pop_back();
          }
          ansScore = CalcScore();
          if (ansScore > best_ansScore) {
            CopyToBest();
          }
        }
      }
      else {
        bool oniCheck = !(Rand01() < haiparas.SimPartition);
        int res = Sim2(ans2, true);
        if (res == -1) {
          swap(ans[ra1], ans[ra2]);
        }
        else {
          double diffScore = ((double)ans.size() - ans2.size()) * haiparas.MultipleValue;
          double prob = exp(diffScore / temp);

          if (prob > Rand01()) {
            loop2++;
            // 採用
            ans = ans2;
            ansScore = CalcScore();
            if (ansScore > best_ansScore) {
              CopyToBest();
            }
          }
          else {
            // 元に戻す
            swap(ans[ra1], ans[ra2]);
          }

        }
      }
    }
    else {
      int ra1 = Rand() % (ans.size() - 1);
      int raDir = Rand() % 4;
      int raNum = Rand() % n + 1;
      P keep = ans[ra1];
      ans[ra1] = P(raDir, raNum);


      if (Rand01() < haiparas.SimPartition) {
        int res = Sim();
        if (res == -1) {
          ans[ra1] = keep;
        }
        else {
          loop2++;
          while (ans.size() > res) {
            ans.pop_back();
          }
          ansScore = CalcScore();
          if (ansScore > best_ansScore) {
            CopyToBest();
          }
        }
      }
      else {
        bool oniCheck = !(Rand01() < haiparas.SimPartition);
        int res = Sim2(ans2, true);
        if (res == -1) {
          ans[ra1] = keep;
        }
        else {
          double diffScore = ((double)ans.size() - ans2.size()) * haiparas.MultipleValue;
          double prob = exp(diffScore / temp);

          if (prob > Rand01()) {
            loop2++;
            // 採用
            ans = ans2;
            ansScore = CalcScore();
            if (ansScore > best_ansScore) {
              CopyToBest();
            }
          }
          else {
            // 元に戻す
            ans[ra1] = keep;
          }
        }
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << " " << loop2 << endl;
  }

  CopyToAns();
}

// 問題を解く関数
ll Solve(int problem_num, Haiparas haiparas)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // 初期解生成
  Method2();

  Mountain(haiparas);

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
    Rand();
  }

  mode = 2;

  Haiparas HAIPARAS;
  HAIPARAS.StartTemp = 302.535;
  HAIPARAS.EndTemp = 0.0;
  HAIPARAS.MultipleValue = 826.94;
  HAIPARAS.Partition1 = 94;
  HAIPARAS.SimPartition = 0.540721;

  if (mode == 0) {
    Solve(0, HAIPARAS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 15)
    {
      ll score = Solve(i, HAIPARAS);
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
    Haiparas bestHaiparas;
    ll bestSumScore = 0;

    while (true) {
      Haiparas haiparas;
      haiparas.StartTemp = pow(2.0, Rand01() * 10);
      haiparas.EndTemp = 0.0;
      haiparas.MultipleValue = pow(2.0, Rand01() * 10);
      haiparas.Partition1 = Rand() % 101;
      haiparas.SimPartition = Rand01();

      ll sum = 0;
      srep(i, 0, 15)
      {
        ll score = Solve(i, haiparas);
        sum += score;
        if (i == 0 && score < 3120) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << haiparas.StartTemp
        << ", EndTemp = " << haiparas.EndTemp
        << ", MultipleValue = " << haiparas.MultipleValue
        << ", Partition1 = " << haiparas.Partition1
        << ", SimPartition = " << haiparas.SimPartition
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHaiparas = haiparas;
      }

      loop++;
    }
  }

  return 0;
}
