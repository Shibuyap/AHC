#include <algorithm>
#include <array>
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

typedef pair<int, int> P;
typedef long long int ll;

// タイマー
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

Timer timer;

// 乱数
namespace
{
  static uint32_t rand_xorshift()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
  }

  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

// ユークリッド距離を計算する関数
double euclidean_distance(int x1, int y1, int x2, int y2)
{
  return sqrt(static_cast<double>((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

const double TIME_LIMIT = 1.9;
int exec_mode;

enum class PlaceType
{
  PROCESSOR = 0,
  SORTER = 1,
  INLET = 2,
};

class Place
{
public:
  PlaceType type;
  int id;
  int x, y;
  int k;
  int v1, v2;
};


// 入力データ
int n, m, k;
vector<Place> places;
vector<vector<double>> p;

// 出力データ
vector<int> d;
int s;

// 入力データの読み込み
void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> n >> m >> k;

    places.resize(n + m + 1);
    for (int i = 0; i < n; i++) {
      int x, y;
      cin >> x >> y;
      places[i] = { PlaceType::PROCESSOR, i, x, y, -1, 0, 0 };
    }
    for (int i = 0; i < m; i++) {
      int x, y;
      cin >> x >> y;
      places[n + i] = { PlaceType::SORTER, i, x, y, -1, 0, 0 };
    }
    places[n + m] = { PlaceType::INLET, 0, 0, 5000, -1, 0, 0 };

    p.assign(k, vector<double>(n));
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < n; j++) {
        cin >> p[i][j];
      }
    }
  }
  else {
    // ファイル入力
    ifs >> n >> m >> k;

    places.resize(n + m + 1);
    for (int i = 0; i < n; i++) {
      int x, y;
      ifs >> x >> y;
      places[i] = { PlaceType::PROCESSOR, i, x, y, -1, 0, 0 };
    }
    for (int i = 0; i < m; i++) {
      int x, y;
      ifs >> x >> y;
      places[n + i] = { PlaceType::SORTER, i, x, y, -1, 0, 0 };
    }
    places[n + m] = { PlaceType::INLET, 0, 0, 5000, -1, 0, 0 };

    p.assign(k, vector<double>(n));
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < n; j++) {
        ifs >> p[i][j];
      }
    }

    ifs.close();
  }
}

void output_data(int case_num)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < n; i++) {
      cout << d[i] << " ";
    }
    cout << endl;
    cout << s << endl;
    for (int i = n; i < n + m; i++) {
      if (places[i].k == -1) {
        cout << -1 << endl;
      }
      else {
        cout << places[i].k << " " << places[i].v1 << " " << places[i].v2 << endl;
      }
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    if (ofs.is_open()) {
      for (int i = 0; i < n; i++) {
        ofs << d[i] << " ";
      }
      ofs << endl;
      ofs << s << endl;
      for (int i = n; i < n + m; i++) {
        if (places[i].k == -1) {
          ofs << -1 << endl;
        }
        else {
          ofs << places[i].k << " " << places[i].v1 << " " << places[i].v2 << endl;
        }
      }

      ofs.close();
    }
  }
}

vector<int> GetNearOrder(int i)
{
  vector<pair<double, int>> distances;
  for (int j = 0; j < n + m; j++) {
    if (j == i) {
      continue;
    }
    double dist = euclidean_distance(places[i].x, places[i].y, places[j].x, places[j].y);
    distances.push_back({ dist, j });
  }
  sort(distances.begin(), distances.end());
  vector<int> near_order;
  for (const auto& p : distances) {
    near_order.push_back(p.second);
  }
  return near_order;
}

void initialize()
{
  d.resize(n);
  for (int i = 0; i < n; i++) {
    d[i] = i;
  }

  // s決定
  {
    vector<int> near_order = GetNearOrder(n + m);
    s = near_order[0]; // 最も近い場所をsに設定
  }

  vector<vector<int>> parents(n + m + 1);
  queue<int> q;
  if (s >= n) {
    q.push(s);
    parents[s].push_back(n + m);
  }

  while (!q.empty()) {
    int current = q.front();
    q.pop();
    auto near_order = GetNearOrder(current);
    places[current].k = 0;
    places[current].v1 = -1;
    places[current].v2 = -1;
    for (int next : near_order) {
      if (next < n) {
        if (places[current].v1 == -1) {
          places[current].v1 = next;
        }
        else {
          places[current].v2 = next;
        }
      }
      else {
        if (places[current].k != -1) {
          continue;
        }

      }

      if (places[current].v2 != -1) {
        break;
      }
    }
  }
}

ll calculate_score()
{
  ll res = 0;
  return res;
}

ll solve_case(int case_num)
{
  timer.start();

  input_data(case_num);

  initialize();

  ll score = 0;

  output_data(case_num);

  return score;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 100; i++) {
      ll score = solve_case(i);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << timer.get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
