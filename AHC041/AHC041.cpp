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

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;

const int N = 1000;
const int H = 10;
int A[N];
vector<int> G[N];
int X[N], Y[N];

class Answer
{
public:
  vector<int> parent;

  Answer()
  {
    parent.resize(N);
    Init();
  }

  void Init()
  {
    for (int i = 0; i < N; ++i) {
      parent[i] = -1;
    }
  }
};

class Heights
{
public:
  vector<int> height;

  Heights()
  {
    height.resize(N);
    Init();
  }

  void Init()
  {
    for (int i = 0; i < N; ++i) {
      height[i] = -1;
    }
  }

  int calc_score()
  {
    int res = 1;
    for (int i = 0; i < N; ++i) {
      res += A[i] * (height[i] + 1);
    }
    return res;
  }
};

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  for (int i = 0; i < N; ++i) G[i].clear();
}

// 入力受け取り
void Input(int num)
{
  // ファイルパス生成
  string filePath = "./in/" + to_string(num / 1000) + to_string((num / 100) % 10) +
    to_string((num / 10) % 10) + to_string(num % 10) + ".txt";

  // 入力ストリーム選択（ファイルか標準入力）
  ifstream ifs(filePath);
  istream& input = ifs.is_open() ? ifs : cin;

  // データ読み込み
  int m, _n, _h;
  input >> _n >> m >> _h;

  for (int i = 0; i < N; ++i) input >> A[i];

  for (int i = 0; i < m; ++i) {
    int u, v;
    input >> u >> v;
    G[u].push_back(v);
    G[v].push_back(u);
  }

  for (int i = 0; i < N; ++i) input >> X[i] >> Y[i];
}

void output(int case_num, const Answer& ans)
{
  auto output_answer = [&](ostream& out) {
    for (int i = 0; i < N; ++i) out << ans.parent[i] << ' ';
    out << endl;
    };

  if (mode == 0) {
    output_answer(cout);
  }
  else {
    ofstream ofs("./out/" + to_string(case_num / 1000) + to_string((case_num / 100) % 10) +
      to_string((case_num / 10) % 10) + to_string(case_num % 10) + ".txt");
    if (ofs.is_open()) {
      output_answer(ofs);
    }
  }
}

Answer heights_to_ans(const Heights& heights)
{
  Answer ans;
  for (int i = 0; i < N; ++i) {
    if (heights.height[i] == 0) {
      ans.parent[i] = -1;
    }
    for (auto y : G[i]) {
      if (heights.height[y] == heights.height[i] + 1) {
        ans.parent[y] = i;
      }
    }
  }
  return ans;
}

bool can_reach_all(Heights& heights, int start_h)
{
  queue<int> q;
  for (int i = 0; i < N; ++i) {
    if (heights.height[i] == start_h && start_h < H) {
      q.push(i);
    }
  }
  while (q.size()) {
    int x = q.front();
    q.pop();
    for (auto y : G[x]) {
      if (heights.height[y] == -1) {
        heights.height[y] = heights.height[x] + 1;
        if (heights.height[y] < H) {
          q.push(y);
        }
      }
    }
  }

  bool res = true;

  for (int i = 0; i < N; ++i) {
    if (heights.height[i] == -1) {
      res = false;
    }
  }

  // 元に戻す
  for (int i = 0; i < N; ++i) {
    if (heights.height[i] >= start_h + 1) {
      heights.height[i] = -1;
    }
  }

  return res;
}

void greedy_assign(Heights& heights)
{
  heights.Init();
  for (int i = 0; i <= H; i++) {
    set<int> candidates;
    if (i == 0) {
      for (int j = 0; j < N; ++j) {
        candidates.insert(j);
      }
    }
    else {
      for (int j = 0; j < N; ++j) {
        if (heights.height[j] == i - 1) {
          for (auto y : G[j]) {
            if (heights.height[y] == -1) {
              candidates.insert(y);
            }
          }
        }
      }
    }

    vector<pair<int, int>> nodes;
    for (auto num : candidates) {
      heights.height[num] = i;
      nodes.push_back(make_pair(A[num], num));
    }
    sort(nodes.begin(), nodes.end(), greater<pair<int, int>>());
    for (auto p : nodes) {
      int num = p.second;
      heights.height[num] = -1;
      if (!can_reach_all(heights, i)) {
        heights.height[num] = i;
      }
    }
  }
}

int Solve(int num)
{
  SetUp();

  Input(num);

  Heights heights;
  greedy_assign(heights);

  Answer ans = heights_to_ans(heights);

  output(num, ans);

  int score = 0;
  if (mode != 0) {
    score = heights.calc_score();
  }
  return score;
}

int main()
{
  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
