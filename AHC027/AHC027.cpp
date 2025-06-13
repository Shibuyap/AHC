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
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const int INF = 1001001001;

double TL = 1.8;
int mode;

const int MAX_LENGTH = 100000;
int N;
vector<int> G[2000];
int d[2000];
vector<int> ans;
vector<int> V[2000];

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  ans.clear();
  rep(i, 2000)
  {
    G[i].clear();
    V[i].clear();
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
    cin >> N;
    rep(i, N - 1)
    {
      string s;
      cin >> s;
      rep(j, N)
      {
        if (s[j] == '0') {
          G[i * N + j].push_back((i + 1) * N + j);
          G[(i + 1) * N + j].push_back(i * N + j);
        }
      }
    }
    rep(i, N)
    {
      string s;
      cin >> s;
      rep(j, N - 1)
      {
        if (s[j] == '0') {
          G[i * N + j].push_back(i * N + j + 1);
          G[i * N + j + 1].push_back(i * N + j);
        }
      }
    }
    rep(i, N * N) { cin >> d[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N;
    rep(i, N - 1)
    {
      string s;
      ifs >> s;
      rep(j, N)
      {
        if (s[j] == '0') {
          G[i * N + j].push_back((i + 1) * N + j);
          G[(i + 1) * N + j].push_back(i * N + j);
        }
      }
    }
    rep(i, N)
    {
      string s;
      ifs >> s;
      rep(j, N - 1)
      {
        if (s[j] == '0') {
          G[i * N + j].push_back(i * N + j + 1);
          G[i * N + j + 1].push_back(i * N + j);
        }
      }
    }
    rep(i, N * N) { ifs >> d[i]; }
  }
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

// スコア計算
ll CalcScore()
{
  ll res = 0;
  return res;
}

// 解答出力
void Output(ofstream& ofs)
{
  string ansString;
  rep(i, ans.size() - 1)
  {
    int diff = ans[i + 1] - ans[i];
    if (diff == N) {
      ansString += 'D';
    }
    else if (diff == -N) {
      ansString += 'U';
    }
    else if (diff == 1) {
      ansString += 'R';
    }
    else if (diff == -1) {
      ansString += 'L';
    }
  }
  if (mode == 0) {
    cout << ansString << endl;
  }
  else {
    ofs << ansString << endl;
  }
}

int visited[2000];
void MakeLoopDfs(int now, vector<int>& vec)
{
  visited[now] = 1;
  vec.push_back(now);
  for (auto nxt : G[now]) {
    if (visited[nxt] == 0) {
      MakeLoopDfs(nxt, vec);
      vec.push_back(now);
    }
  }
}

vector<int> MakeLoop(int st)
{
  rep(i, N * N) { visited[i] = 0; }
  vector<int> res;
  MakeLoopDfs(st, res);
  return res;
}

ll Solve(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  vector<int> loop = MakeLoop(0);
  ans = loop;

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

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
  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
