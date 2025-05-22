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
clock_t startTime, endTime;

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

const int N = 1000;
const int H = 10;
int A[N];
vector<int> G[N];
int X[N], Y[N];

class Answer
{
public:
  vector<int> p;

  Answer()
  {
    p.resize(N);
    for (int i = 0; i < N; ++i) {
      p[i] = -1;
    }
  }

  void Init()
  {
    for (int i = 0; i < N; ++i) {
      p[i] = -1;
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
    for (int i = 0; i < N; ++i) {
      height[i] = -1;
    }
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
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // 標準入力する
  int m, _n, _h;
  if (!ifs.is_open()) {
    cin >> _n >> m >> _h;
    for (int i = 0; i < N; ++i) cin >> A[i];
    for (int i = 0; i < m; ++i) {
      int u, v;
      cin >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    for (int i = 0; i < N; ++i) cin >> X[i] >> Y[i];
  }
  // ファイル入力する
  else {
    ifs >> _n >> m >> _h;
    for (int i = 0; i < N; ++i) ifs >> A[i];
    for (int i = 0; i < m; ++i) {
      int u, v;
      ifs >> u >> v;
      G[u].push_back(v);
      G[v].push_back(u);
    }
    for (int i = 0; i < N; ++i) ifs >> X[i] >> Y[i];
  }
}

void output_data(int case_num, const Answer& ans)
{
  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < N; ++i) cout << ans.p[i] << ' ';
    cout << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < N; ++i) ofs << ans.p[i] << ' ';
    ofs << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

Answer convert_heights_to_answer(const Heights& heights)
{
  Answer ans;
  for (int i = 0; i < N; ++i) {
    ans.p[i] = -1;
  }
  for (int i = 0; i < N; ++i) {
    if (heights.height[i] == 0) {
      ans.p[i] = -1;
    }
    for (auto y : G[i]) {
      if (heights.height[y] == heights.height[i] + 1) {
        ans.p[y] = i;
      }
    }
  }
  return ans;
}

bool check_if_valid(Heights& heights, int start_height)
{
  queue<int> que;
  for (int i = 0; i < N; ++i) {
    if (heights.height[i] == start_height && start_height < H) {
      que.push(i);
    }
  }
  while (que.size()) {
    int x = que.front();
    que.pop();
    for (auto y : G[x]) {
      if (heights.height[y] == -1) {
        heights.height[y] = heights.height[x] + 1;
        if (heights.height[y] < H) {
          que.push(y);
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
    if (heights.height[i] >= start_height + 1) {
      heights.height[i] = -1;
    }
  }

  return res;
}

void greedy_1(Heights& heights)
{
  heights.Init();
  for (int i = 0; i <= H; i++) {
    set<int> se;
    if (i == 0) {
      for (int j = 0; j < N; ++j) {
        se.insert(j);
      }
    }
    else {
      for (int j = 0; j < N; ++j) {
        if (heights.height[j] == i - 1) {
          for (auto y : G[j]) {
            if (heights.height[y] == -1) {
              se.insert(y);
            }
          }
        }
      }
    }

    vector<pair<int, int>> vp;
    for (auto num : se) {
      heights.height[num] = i;
      vp.push_back(make_pair(A[num], num));
    }
    sort(vp.begin(), vp.end(), greater<pair<int, int>>());
    for (auto p : vp) {
      int num = p.second;
      heights.height[num] = -1;
      if (!check_if_valid(heights, i)) {
        heights.height[num] = i;
      }
    }
  }
}

ll Solve(int probNum)
{
  startTime = clock();

  SetUp();

  Input(probNum);

  Heights heights;
  heights.Init();
  greedy_1(heights);

  Answer ans = convert_heights_to_answer(heights);

  output_data(probNum, ans);

  ll score = 0;
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
