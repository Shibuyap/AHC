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

namespace /* 乱数ライブラリ */
{
  static uint32_t rand_xorshift()
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


  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
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
  for (int i = 0; i < 2000; ++i) {
    G[i].clear();
    V[i].clear();
  }
}

// 入力受け取り
void Input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> N;
    for (int i = 0; i < N - 1; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < N; ++j) {
        if (s[j] == '0') {
          G[i * N + j].push_back((i + 1) * N + j);
          G[(i + 1) * N + j].push_back(i * N + j);
        }
      }
    }
    for (int i = 0; i < N; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < N - 1; ++j) {
        if (s[j] == '0') {
          G[i * N + j].push_back(i * N + j + 1);
          G[i * N + j + 1].push_back(i * N + j);
        }
      }
    }
    for (int i = 0; i < N * N; ++i) { cin >> d[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N;
    for (int i = 0; i < N - 1; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < N; ++j) {
        if (s[j] == '0') {
          G[i * N + j].push_back((i + 1) * N + j);
          G[(i + 1) * N + j].push_back(i * N + j);
        }
      }
    }
    for (int i = 0; i < N; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < N - 1; ++j) {
        if (s[j] == '0') {
          G[i * N + j].push_back(i * N + j + 1);
          G[i * N + j + 1].push_back(i * N + j);
        }
      }
    }
    for (int i = 0; i < N * N; ++i) { ifs >> d[i]; }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int case_num, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(case_num % 10 + '0');
      case_num /= 10;
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
  for (int i = 0; i < ans.size() - 1; ++i) {
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
  for (int i = 0; i < N * N; ++i) { visited[i] = 0; }
  vector<int> res;
  MakeLoopDfs(st, res);
  return res;
}

ll Solve(int case_num)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(case_num);

  vector<int> loop = MakeLoop(0);
  ans = loop;

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(case_num, ofs);

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
    for (int i = 0; i < 100; ++i) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
