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

#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

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

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

bool is_out_of_bounds(int x, int y, int n)
{
  if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

double TIME_LIMIT = 1.8;
int exec_mode;

const int BOARD_SIZE = 6;
const int ITEM_KIND = 15;
const int TURN_COUNT = 10;
const int EDGE_COUNT = 2 * BOARD_SIZE * (BOARD_SIZE - 1);
vector<vector<int>> item_value(EDGE_COUNT, vector<int>(ITEM_KIND, 0));
int item_value_max[ITEM_KIND];
int item_value_pow3[100][ITEM_KIND];
int sorted_order[EDGE_COUNT];

vector<vector<int>> placement(BOARD_SIZE, vector<int>(BOARD_SIZE, 0));

int horizontal_edge_block[TURN_COUNT][BOARD_SIZE][BOARD_SIZE][ITEM_KIND];
int vertical_edge_block[TURN_COUNT][BOARD_SIZE][BOARD_SIZE][ITEM_KIND];

// 複数ケース回すときに内部状態を初期値に戻す
void reset_state() {}

// 入力受け取り
void read_input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    int n, m, t;
    cin >> n >> m >> t;
    for (int i = 0; i < EDGE_COUNT; i++) {
      for (int j = 0; j < ITEM_KIND; j++) {
        cin >> item_value[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    int n, m, t;
    ifs >> n >> m >> t;
    for (int i = 0; i < EDGE_COUNT; i++) {
      for (int j = 0; j < ITEM_KIND; j++) {
        ifs >> item_value[i][j];
      }
    }
    rep(i, TURN_COUNT)
    {
      rep(j, BOARD_SIZE)
      {
        rep(k, BOARD_SIZE - 1)
        {
          string s;
          ifs >> s;
          rep(l, ITEM_KIND)
          {
            horizontal_edge_block[i][j][k][l] = s[l] - '0';
          }
        }
      }
      rep(j, BOARD_SIZE - 1)
      {
        rep(k, BOARD_SIZE)
        {
          string s;
          ifs >> s;
          rep(l, ITEM_KIND)
          {
            vertical_edge_block[i][j][k][l] = s[l] - '0';
          }
        }
      }
    }
  }

  rep(j, ITEM_KIND)
  {
    item_value_max[j] = 0;
  }
  rep(i, EDGE_COUNT)
  {
    rep(j, ITEM_KIND)
    {
      item_value_max[j] = max(item_value_max[j], item_value[i][j]);
    }
  }
}

// 出力ファイルストリームオープン
void open_output_file(int probNum, ofstream& ofs)
{
  if (exec_mode != 0) {
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
ll calc_final_score()
{
  ll sumX = 0;
  rep(j, ITEM_KIND)
  {
    sumX += item_value_max[j];
  }

  int ma = 0;
  rep(i, EDGE_COUNT)
  {
    int tmp = 0;
    rep(j, ITEM_KIND)
    {
      tmp += item_value[i][j];
    }
    ma = max(ma, tmp);
  }

  ll res = round(1000000.0 * ma / sumX);
  return res;
}

void apply_feedback(int turn)
{
  // placementに基づいて itemValue を再設定
  int idx = 0;
  // horizontal
  for (int row = 0; row < BOARD_SIZE; row++) {
    for (int col = 0; col < BOARD_SIZE - 1; col++) {
      for (int k = 0; k < ITEM_KIND; k++) {
        if (horizontal_edge_block[turn][row][col][k] == 0) {
          item_value[idx][k] = item_value_pow3[placement[row][col]][k];
        }
        else {
          item_value[idx][k] = item_value_pow3[placement[row][col + 1]][k];
        }
      }
      idx++;
    }
  }
  // vertical
  for (int row = 0; row < BOARD_SIZE - 1; row++) {
    for (int col = 0; col < BOARD_SIZE; col++) {
      for (int k = 0; k < ITEM_KIND; k++) {
        if (vertical_edge_block[turn][row][col][k] == 0) {
          item_value[idx][k] = item_value_pow3[placement[row][col]][k];
        }
        else {
          item_value[idx][k] = item_value_pow3[placement[row + 1][col]][k];
        }
      }
      idx++;
    }
  }
}

void precompute_values()
{
  vector<P> v[ITEM_KIND];
  rep(i, EDGE_COUNT)
  {
    rep(j, ITEM_KIND)
    {
      v[j].push_back(P(item_value[i][j], i));
    }
  }
  rep(i, ITEM_KIND)
  {
    sort(v[i].begin(), v[i].end(), greater<P>());
  }
  int val[EDGE_COUNT][ITEM_KIND];
  rep(i, EDGE_COUNT)
  {
    rep(j, ITEM_KIND)
    {
      val[v[j][i].second][j] = i;
    }
  }

  int ma[ITEM_KIND] = {};
  rep(i, EDGE_COUNT)
  {
    rep(j, ITEM_KIND)
    {
      ma[j] = max(ma[j], item_value[i][j]);
    }
  }

  rep(i, EDGE_COUNT)
  {
    rep(j, ITEM_KIND)
    {
      item_value_pow3[i][j] = item_value[i][j];
      item_value[i][j] = item_value[i][j] * item_value[i][j] * item_value[i][j];
    }
  }
}

void sort_by_total_value()
{
  vector<P> vp;
  rep(i, EDGE_COUNT)
  {
    int sum = 0;
    rep(j, ITEM_KIND)
    {
      sum += item_value[i][j];
    }
    vp.push_back(P(sum, i));
  }
  sort(vp.begin(), vp.end(), greater<P>());
  rep(i, vp.size())
  {
    sorted_order[i] = vp[i].second;
  }
}

const int SPIRAL_RANK[6][6] = {
    {21, 20, 19, 18, 17, 16},
    {22,  7,  6,  5, 16, 35},
    {23,  8,  1,  4, 15, 34},
    {24,  9,  2,  3, 14, 33},
    {25, 10, 11, 12, 13, 32},
    {26, 27, 28, 29, 30, 31},
};

vector<P> spiral_cells;
void init_spiral_order()
{
  spiral_cells.clear();
  vector<pair<int, P>> vp;
  rep(i, 6)
  {
    rep(j, 6)
    {
      vp.push_back(make_pair(SPIRAL_RANK[i][j], P(i, j)));
    }
  }
  sort(vp.begin(), vp.end());
  for (auto p : vp) {
    spiral_cells.push_back(p.second);
  }
}

ll local_neighbor_score(int x, int y)
{
  ll res = 0;
  int id = placement[x][y];
  rep(i, 4)
  {
    int nx = x + DX[i];
    int ny = y + DY[i];
    if (is_out_of_bounds(nx, ny, BOARD_SIZE)) {
      continue;
    }
    int nd = placement[nx][ny];
    rep(j, ITEM_KIND)
    {
      res += 10 * item_value[id][j] + 2 * abs(item_value[id][j] - item_value[nd][j]);
    }
  }
  return res;
}

void anneal_spiral_placement()
{
  rep(i, 36)
  {
    int x = spiral_cells[i].first;
    int y = spiral_cells[i].second;
    placement[x][y] = sorted_order[i];
  }

  const int LOOP_COUNT = 500000;
  rep(_, LOOP_COUNT)
  {
    int ra1 = rand_xorshift() % 6;
    int ra2 = rand_xorshift() % 6;
    int ra3 = rand_xorshift() % 6;
    int ra4 = rand_xorshift() % 6;

    int before = 0;
    before += local_neighbor_score(ra1, ra2);
    before += local_neighbor_score(ra3, ra4);

    swap(placement[ra1][ra2], placement[ra3][ra4]);

    int after = 0;
    after += local_neighbor_score(ra1, ra2);
    after += local_neighbor_score(ra3, ra4);

    double temp = (double)(LOOP_COUNT - _) / LOOP_COUNT * 500000;
    const double prob = exp((double)(after - before) / temp);

    if (prob > rand_01()) {
    }
    else {
      swap(placement[ra1][ra2], placement[ra3][ra4]);
    }
  }
}

ll solve_case(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  reset_state();

  init_spiral_order();

  // 入力受け取り
  read_input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  open_output_file(probNum, ofs);

  for (int t = 0; t < TURN_COUNT; t++) {
    for (int i = 0; i < BOARD_SIZE; i++) {
      for (int j = 0; j < BOARD_SIZE; j++) {
        placement[i][j] = i * BOARD_SIZE + j;
      }
    }

    precompute_values();
    sort_by_total_value();
    anneal_spiral_placement();

    for (int i = 0; i < BOARD_SIZE; i++) {
      for (int j = 0; j < BOARD_SIZE; j++) {
        if (exec_mode == 0) {
          cout << placement[i][j];
          if (j < BOARD_SIZE - 1) {
            cout << " ";
          }
          else {
            cout << endl;
          }
        }
        else {
          ofs << placement[i][j];
          if (j < BOARD_SIZE - 1) {
            ofs << " ";
          }
          else {
            ofs << endl;
          }
        }
      }
    }

    if (exec_mode == 0) {
      cout.flush();
      for (int i = 0; i < EDGE_COUNT; i++) {
        for (int j = 0; j < ITEM_KIND; j++) {
          cin >> item_value[i][j];
        }
      }
    }
    else {
      apply_feedback(t);
    }
  }

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = calc_final_score();
  }

  return score;
}

int main()
{
  exec_mode = 1;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = solve_case(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
