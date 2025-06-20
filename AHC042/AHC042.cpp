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
typedef pair<P, P> PP;


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


static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
}


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
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

std::chrono::steady_clock::time_point startTimeClock;

void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTimeClock;
  return elapsed.count();
}


const int MAX_N = 30;

const int n = 20;

int a[n + 2][n + 2];
int original_a[n + 2][n + 2];

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
    for (int i = 0; i < n; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < n; ++j) {
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
    for (int i = 0; i < n; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < n; ++j) {
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

  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      original_a[i][j] = a[i][j];
    }
  }
}

// 出力ファイルストリームを開く関数
void OpenOutput(int probNum, ofstream& ofs)
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
  for (int i = 0; i < n; ++i) {
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
  for (int j = 0; j < n; ++j) {
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

int CalcMoveScore(int dir, int num)
{
  int score = 0;
  queue<P> q;
  int f[n + 2][n + 2];
  for (int i = 1; i < n + 1; ++i) {
    for (int j = 1; j < n + 1; ++j) {
      if (a[i][j] == 1) {
        while (q.size()) {
          q.pop();
        }
        for (int k = 0; k < n + 2; ++k) {
          for (int l = 0; l < n + 2; ++l) {
            f[k][l] = 99999;
          }
        }
        int dist = 99999;
        f[i][j] = 0;
        q.push(P(i, j));
        while (q.size()) {
          int x = q.front().first;
          int y = q.front().second;
          q.pop();
          if (x == 1 || x == n || y == 1 || y == n) {
            dist = f[x][y] + 1;
            break;
          }
          for (int k = 0; k < 4; ++k) {
            int nx = x + dx[k];
            int ny = y + dy[k];
            if (IsNG(nx, ny))continue;
            if (a[nx][ny] == 2)continue;
            if (f[nx][ny] > f[x][y] + 1) {
              f[nx][ny] = f[x][y] + 1;
              q.push(P(nx, ny));
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
void InitialSolve()
{
  int trial = 0;

  while (trial == 0) {

    if (GetNowTime() > TL) {
      break;
    }
    trial++;

    ans.clear();
    int cnt = 0;
    int loop = 0;

    int prev_dir = -1;
    int prev_num = -1;

    while (cnt < n * 2 && loop < 500) {


      loop++;

      int dir = -1;
      int num = -1;
      int cost = 99999;

      // U
      for (int j = 1; j < n + 1; ++j) {
        if (a[1][j] == 2) {
          continue;
        }
        if (prev_dir == 2 && prev_num == j)continue;
        MoveU(j, true);
        int score = CalcMoveScore(0, j);
        if (score < cost) {
          dir = 0;
          num = j;
          cost = score;
        }
        MoveD(j, false);
      }

      // L
      for (int i = 1; i < n + 1; ++i) {
        if (a[i][1] == 2) {
          continue;
        }
        if (prev_dir == 3 && prev_num == i)continue;
        MoveL(i, true);
        int score = CalcMoveScore(1, i);
        if (score < cost) {
          dir = 1;
          num = i;
          cost = score;
        }
        MoveR(i, false);
      }

      // D
      for (int j = 1; j < n + 1; ++j) {
        if (a[n][j] == 2) {
          continue;
        }
        if (prev_dir == 0 && prev_num == j)continue;
        MoveD(j, true);
        int score = CalcMoveScore(2, j);
        if (score < cost) {
          dir = 2;
          num = j;
          cost = score;
        }
        MoveU(j, false);
      }

      // R
      for (int i = 1; i < n + 1; ++i) {
        if (a[i][n] == 2) {
          continue;
        }
        if (prev_dir == 1 && prev_num == i)continue;
        MoveR(i, true);
        int score = CalcMoveScore(3, i);
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

      prev_dir = dir;
      prev_num = num;
    }

    ansScore = CalcScore();

    if (ansScore > best_ansScore) {
      CopyToBest();
    }
  }

  CopyToAns();
}

void ResetBoard()
{
  for (int i = 0; i < n + 2; ++i) {
    for (int j = 0; j < n + 2; ++j) {
      a[i][j] = original_a[i][j];
    }
  }
}

int Sim()
{
  ResetBoard();
  int cnt = 0;
  for (int i = 0; i < ans.size(); ++i) {
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

int Sim2(vector<P>& ans2, bool check_x)
{
  ans2.clear();
  ResetBoard();
  int cnt = 0;
  for (int i = 0; i < ans.size(); ++i) {
    int dir = ans[i].first;
    int num = ans[i].second;
    if (dir == 0) {
      if (a[1][num] == 1) {
        cnt++;
      }
      else if (a[1][num] == 2) {
        return -1;
      }

      if (check_x) {
        int has_x = 0;
        for (int i = 1; i < n + 1; ++i) {
          if (a[i][num] == 1) {
            has_x = 1;
            break;
          }
        }
        if (has_x == 0)continue;
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

      if (check_x) {
        int has_x = 0;
        for (int j = 1; j < n + 1; ++j) {
          if (a[num][j] == 1) {
            has_x = 1;
            break;
          }
        }
        if (has_x == 0)continue;
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

      if (check_x) {
        int has_x = 0;
        for (int i = 1; i < n + 1; ++i) {
          if (a[i][num] == 1) {
            has_x = 1;
            break;
          }
        }
        if (has_x == 0)continue;
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

      if (check_x) {
        int has_x = 0;
        for (int j = 1; j < n + 1; ++j) {
          if (a[num][j] == 1) {
            has_x = 1;
            break;
          }
        }
        if (has_x == 0)continue;
      }

      MoveR(num, true);
      ans2.push_back(ans[i]);
    }
    if (cnt == n * 2) {
      return ans2.size();
    }
  }

  for (int i = 1; i < n + 1; ++i) {
    //if (ans2.size() >= ansScore + 5)break;
    for (int j = 1; j < n + 1; ++j) {
      //if (ans2.size() >= ansScore + 5)break;

      if (a[i][j] == 1) {
        int dir = -1;
        int cost = 999;

        // U
        int ok = 1;
        for (int k = 1; k < i; ++k) {
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
        for (int k = 1; k < j; ++k) {
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
        for (int k = i + 1; k < n + 1; ++k) {
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
        for (int k = j + 1; k < n + 1; ++k) {
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
          for (int k = 0; k < cost; ++k) {
            ans2.push_back(P(dir, j));
            MoveU(j, true);
          }
        }
        else if (dir == 1) {
          for (int k = 0; k < cost; ++k) {
            ans2.push_back(P(dir, i));
            MoveL(i, true);
          }
        }
        else if (dir == 2) {
          for (int k = 0; k < cost; ++k) {
            ans2.push_back(P(dir, j));
            MoveD(j, true);
          }
        }
        else if (dir == 3) {
          for (int k = 0; k < cost; ++k) {
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

struct Params
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition1;
  double SimPartition;
};

void HillClimb(Params params)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = params.StartTemp;
  const double END_TEMP = params.EndTemp;

  int loop = 0;
  int loop2 = 0;
  vector<P> ans2;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) { break; }
    }

    if (ans.size() >= best_ans.size() + 4) {
      CopyToAns();
    }

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    int r_mode = Rand() % 100;
    if (r_mode < params.Partition1) {
      // swap
      int r1 = Rand() % ans.size();
      int r2 = Rand() % ans.size();
      while (ans[r1] == ans[r2]) {
        r2 = Rand() % ans.size();
      }
      swap(ans[r1], ans[r2]);

      if (Rand01() < params.SimPartition) {
        int res = Sim();
        if (res == -1) {
          swap(ans[r1], ans[r2]);
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
        bool check_x = !(Rand01() < params.SimPartition);
        int res = Sim2(ans2, check_x);
        if (res == -1) {
          swap(ans[r1], ans[r2]);
        }
        else {
          double diffScore = ((double)ans.size() - ans2.size()) * params.MultipleValue;
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
            swap(ans[r1], ans[r2]);
          }
        }
      }
    }
    else {
      int r1 = Rand() % (ans.size() - 1);
      int r_dir = Rand() % 4;
      int r_num = Rand() % n + 1;
      P keep = ans[r1];
      ans[r1] = P(r_dir, r_num);


      if (Rand01() < params.SimPartition) {
        int res = Sim();
        if (res == -1) {
          ans[r1] = keep;
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
        bool check_x = !(Rand01() < params.SimPartition);
        int res = Sim2(ans2, check_x);
        if (res == -1) {
          ans[r1] = keep;
        }
        else {
          double diffScore = ((double)ans.size() - ans2.size()) * params.MultipleValue;
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
            ans[r1] = keep;
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
ll Solve(int problem_num, Params params)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOutput(problem_num, ofs);

  // 初期解生成
  InitialSolve();

  HillClimb(params);

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
  mode = 2;

  Params PARAMS;
  PARAMS.StartTemp = 302.535;
  PARAMS.EndTemp = 0.0;
  PARAMS.MultipleValue = 826.94;
  PARAMS.Partition1 = 94;
  PARAMS.SimPartition = 0.540721;

  if (mode == 0) {
    Solve(0, PARAMS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = Solve(i, PARAMS);
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
    Params bestParams;
    ll bestSumScore = 0;

    while (true) {
      Params params;
      params.StartTemp = pow(2.0, Rand01() * 10);
      params.EndTemp = 0.0;
      params.MultipleValue = pow(2.0, Rand01() * 10);
      params.Partition1 = Rand() % 101;
      params.SimPartition = Rand01();

      ll sum = 0;
      for (int i = 0; i < 15; ++i) {
        ll score = Solve(i, params);
        sum += score;
        if (i == 0 && score < 3120) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << params.StartTemp
        << ", EndTemp = " << params.EndTemp
        << ", MultipleValue = " << params.MultipleValue
        << ", Partition1 = " << params.Partition1
        << ", SimPartition = " << params.SimPartition
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestParams = params;
      }

      loop++;
    }
  }

  return 0;
}
