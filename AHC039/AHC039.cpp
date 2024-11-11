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

double TL = 1.8; // 時間制限（Time Limit）
int mode;        // 実行モード
std::chrono::steady_clock::time_point startTimeClock, endTimeClock; // 時間計測用

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

// 二次元座標
struct Point
{
public:
  int x;
  int y;

  Point() { x = 0; y = 0; }
  Point(int _x, int _y) { x = _x; y = _y; }
};

const int MAX_N = 30;

const int n = 5000;

vector<Point> saba, iwashi;

vector<Point> ans;

int ansScore;

int best_ansScore;

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  ans.clear();
  saba.clear();
  iwashi.clear();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  saba.resize(n);
  iwashi.resize(n);

  if (!ifs.is_open()) {
    // 標準入力
    int _n;
    cin >> _n;
    rep(i, n)cin >> saba[i].x >> saba[i].y;
    rep(i, n)cin >> iwashi[i].x >> iwashi[i].y;
  }
  else {
    // ファイル入力
    int _n;
    ifs >> _n;
    rep(i, n)ifs >> saba[i].x >> saba[i].y;
    rep(i, n)ifs >> iwashi[i].x >> iwashi[i].y;
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
  ll res = ansScore + 1;
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    cout << ans.size() << endl;
    for (auto p : ans)cout << p.x << ' ' << p.y << endl;
  }
  else {
    // ファイル出力
    ofs << ans.size() << endl;
    for (auto p : ans)ofs << p.x << ' ' << p.y << endl;
  }
}

bool IsLengthOK(vector<Point> vp)
{
  int len = 0;
  rep(i, vp.size() - 1)
  {
    len += abs(vp[i + 1].x - vp[i].x);
    len += abs(vp[i + 1].y - vp[i].y);
  }

  len += abs(vp[0].x - vp.back().x);
  len += abs(vp[0].y - vp.back().y);

  return len <= 400000;
}

const int bSize = 20;
bool IsNG(int x, int y)
{
  if (x < 0 || bSize <= x || y < 0 || bSize <= y)return true;
  return false;
}

void Method3()
{
  int block[bSize][bSize];
  rep(i, bSize)rep(j, bSize)block[i][j] = 0;
  rep(i, n)
  {
    {
      int xx = saba[i].x / (100000 / bSize);
      xx = min(xx, bSize - 1);
      int yy = saba[i].y / (100000 / bSize);
      yy = min(yy, bSize - 1);

      block[xx][yy]++;
    }

    {
      int xx = iwashi[i].x / (100000 / bSize);
      xx = min(xx, bSize - 1);
      int yy = iwashi[i].y / (100000 / bSize);
      yy = min(yy, bSize - 1);

      block[xx][yy]--;
    }
  }

  int xx1, yy1, xx2, yy2;

  ansScore = 0;
  int loop1 = 0;
  while (true) {
    if (loop1 % 100 == 0) {
      if (GetNowTime() > TL / 2)break;
    }
    loop1++;
    int x1 = RandXor() % bSize;
    int x2 = RandXor() % bSize;
    int y1 = RandXor() % bSize;
    int y2 = RandXor() % bSize;
    if (x1 > x2)swap(x1, x2);
    if (y1 > y2)swap(y1, y2);

    int cnt = 0;
    srep(i, x1, x2 + 1)
    {
      srep(j, y1, y2 + 1)
      {
        cnt += block[i][j];
      }
    }

    if (cnt > ansScore) {
      ansScore = cnt;
      ans.clear();
      ans.emplace_back(x1 * (100000 / bSize), y1 * (100000 / bSize));
      ans.emplace_back((x2 + 1) * (100000 / bSize), y1 * (100000 / bSize));
      ans.emplace_back((x2 + 1) * (100000 / bSize), (y2 + 1) * (100000 / bSize));
      ans.emplace_back(x1 * (100000 / bSize), (y2 + 1) * (100000 / bSize));
      xx1 = x1;
      yy1 = y1;
      xx2 = x2;
      yy2 = y2;
    }
  }

  int f[bSize + 2][bSize + 2];
  rep(i, bSize + 2)
  {
    rep(j, bSize + 2)
    {
      f[i][j] = 0;
    }
  }

  srep(i, xx1, xx2 + 1)
  {
    srep(j, yy1, yy2 + 1)
    {
      f[i + 1][j + 1] = 1;
    }
  }

  int haba[bSize + 2][bSize + 2];
  queue<P> que;

  double nowTime = GetNowTime();
  const double START_TEMP = 100.0;
  const double END_TEMP = 0.1;
  double temp = START_TEMP + (END_TEMP - START_TEMP) * nowTime / TL;

  int loop2 = 0;
  while (true) {
    if (loop2 % 1 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL)break;
    }
    loop2++;

    int rax = RandXor() % bSize;
    int ray = RandXor() % bSize;

    int ng = 1;
    rep(i, 4)
    {
      int nx = rax + dx[i];
      int ny = ray + dy[i];
      if (IsNG(nx, ny))continue;
      if (f[nx + 1][ny + 1] != f[rax + 1][ray + 1]) ng = 0;
    }

    if (ng)continue;

    int tmpScore = ansScore;
    if (f[rax + 1][ray + 1] == 0) {
      tmpScore += block[rax][ray];
    }
    else {
      tmpScore += -block[rax][ray];
    }

    const double progressRatio = nowTime / TL;  // 進捗。開始時が0.0、終了時が1.0
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    double diff = tmpScore - ansScore;
    double prob = exp(diff / temp);
    f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];
    int upd = 0;
    if (prob > Rand01()) {
      // 幅優先
      rep(i, bSize + 2)rep(j, bSize + 2)haba[i][j] = 0;
      int now = 1;
      upd = 1;
      srep(i, 1, bSize + 1)
      {
        srep(j, 1, bSize + 1)
        {
          if (haba[i][j] != 0)continue;
          if (now == 3) {
            upd = 0;
            break;
          }

          haba[i][j] = now;
          que.push(P(i, j));
          while (que.size()) {
            int x = que.front().first;
            int y = que.front().second;
            que.pop();
            rep(k, 4)
            {
              int nx = x + dx[k];
              int ny = y + dy[k];
              if (IsNG(nx - 1, ny - 1))continue;
              if (haba[nx][ny] == 0 && f[nx][ny] == f[i][j]) {
                haba[nx][ny] = now;
                que.push(P(nx, ny));
              }
            }
          }

          now++;
        }
        if (upd == 0)break;
      }
      if (now > 3)upd = 0;
    }

    if (upd) {
      auto ans2 = ans;

      ans.clear();

      int sx = -1, sy = -1;
      int befx = -1, befy = -1;
      srep(i, 1, bSize + 1)
      {
        srep(j, 1, bSize + 1)
        {
          if (f[i][j] == 1 && f[i][j - 1] == 0) {
            sx = i;
            sy = j;
            befx = i + 1;
            befy = j;
          }
        }
      }

      if (sx == -1) {
        assert(false);
      }

      vector<P> vp;
      vp.emplace_back(sx, sy);
      while (true) {
        int x = -1, y = -1;

        if (x == -1) {
          int nx = sx - 1;
          int ny = sy;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy - 1] != f[sx - 1][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx + 1;
          int ny = sy;
          if (!(nx == befx && ny == befy)) {
            if (f[sx][sy - 1] != f[sx][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx;
          int ny = sy - 1;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy - 1] != f[sx][sy - 1]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx;
          int ny = sy + 1;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy] != f[sx][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == vp[0].first && y == vp[0].second)break;

        if (x == -1) {
          assert(false);
          for (auto p : vp) cout << p.first << ' ' << p.second << endl;
          cout << sx << ' ' << sy << ' ' << befx << ' ' << befy << endl;
          srep(i, 1, bSize + 1)
          {
            srep(j, 1, bSize + 1)
            {
              cout << f[i][j];
            }
            cout << endl;
          }
        }

        befx = sx;
        befy = sy;
        sx = x;
        sy = y;
        vp.emplace_back(sx, sy);
      }

      for (auto p : vp) {
        int x = (p.first - 1) * (100000 / bSize);
        int y = (p.second - 1) * (100000 / bSize);
        ans.emplace_back(x, y);
      }

      if (IsLengthOK(ans)) {
        upd = 1;
      }
      else {
        upd = 0;
        ans = ans2;
      }
    }

    if (upd) {
      ansScore = tmpScore;
    }
    else {
      f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];
    }
  }

  if (mode != 0) {
    cout << "loop1 = " << loop1 << ", ";
    cout << "loop2 = " << loop2 << ", ";
    cout << endl;
    srep(i, 1, bSize + 1)
    {
      srep(j, 1, bSize + 1)
      {
        cout << f[i][j];
      }
      cout << endl;
    }
  }
}

// 問題を解く関数
ll Solve(int problem_num)
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
  Method3();

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

  if (mode == 0) {
    Solve(0);
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
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
