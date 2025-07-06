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

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

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

  double log_uniform(double low, double high)
  {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> dist(std::log(low), std::log(high));
    return std::exp(dist(rng));
  }
}

const double TIME_LIMIT = 1.9;
int exec_mode;

const int n = 40;
int m;
vector<vector<int>> init_grid(n, vector<int>(n, 0));

vector<vector<int>> input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  int _n;
  if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> m;
    for (int i = 0; i < _n; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < _n; ++j) {
        init_grid[i][j] = (s[j] == '#') ? 1 : 0;
      }
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> m;
    for (int i = 0; i < _n; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < _n; ++j) {
        init_grid[i][j] = (s[j] == '#') ? 1 : 0;
      }
    }
    ifs.close();
  }

  auto grid = init_grid;
  return grid;
}

void output_data(int case_num, vector<P>& ans)
{
  if (exec_mode == 0) {
    // 標準出力
    for (auto p : ans) {
      cout << p.first << ' ' << p.second << endl;
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (auto p : ans) {
      ofs << p.first << ' ' << p.second << endl;
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

ll calculate_score()
{
  ll res = 0;
  return res;
}

void Move(const vector<vector<int>>& g, const vector<vector<double>>& cnt, vector<vector<double>>& cnt3, int dir)
{
  // 左→右
  if (dir == 0) {
    for (int i = 0; i < n; i++) {
      double now = 0;
      for (int j = 0; j < n; j++) {
        now += cnt[i][j];
        if (j == n - 1) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i][j + 1] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 右→左
  if (dir == 1) {
    for (int i = 0; i < n; i++) {
      double now = 0;
      for (int j = n - 1; j >= 0; j--) {
        now += cnt[i][j];
        if (j == 0) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i][j - 1] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 左→右
  if (dir == 2) {
    for (int j = 0; j < n; j++) {
      double now = 0;
      for (int i = 0; i < n; i++) {
        now += cnt[i][j];
        if (i == n - 1) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i + 1][j] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 左→右
  if (dir == 3) {
    for (int j = 0; j < n; j++) {
      double now = 0;
      for (int i = n - 1; i >= 0; i--) {
        now += cnt[i][j];
        if (i == 0) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i - 1][j] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
}

int tttThres = 500;
vector<vector<int>> canErase(n, vector<int>(n, 0));

void Search2(const vector<vector<int>>& g, const vector<vector<double>>& cnt2, vector<vector<double>>& cnt4, vector<vector<double>>& cnt5, int& x, int& y, double thre, int turn)
{
  if (turn < 50) {
    //thre = 0.5;
  }
  x = -1;
  y = -1;
  double mi = thre;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cnt4[i][j] = 0;
      cnt5[i][j] = 0;
    }
  }

  for (int i = 0; i < n; i++) {
    double now = 0;
    int si = 0;
    for (int j = 0; j < n; j++) {
      now += cnt2[i][j];
      if (j == n - 1) {
        for (int k = si; k <= j; k++) {
          cnt4[i][k] = now;
        }
        now = 0;
        si = j + 1;
      }
      else if (g[i][j + 1] == 1) {
        for (int k = si; k <= j; k++) {
          cnt4[i][k] = now;
        }
        now = 0;
        si = j + 1;
      }
    }
  }

  for (int j = 0; j < n; j++) {
    double now = 0;
    int si = 0;
    for (int i = 0; i < n; i++) {
      now += cnt2[i][j];
      if (i == n - 1) {
        for (int k = si; k <= i; k++) {
          cnt5[k][j] = now;
        }
        now = 0;
        si = i + 1;
      }
      else if (g[i][j + 1] == 1) {
        for (int k = si; k <= i; k++) {
          cnt5[k][j] = now;
        }
        now = 0;
        si = i + 1;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (g[i][j] == 1) {
        continue;
      }

      double tmp = cnt4[i][j] + cnt5[i][j];
      if (turn < tttThres && canErase[i][j] == 0) {
        tmp += 1e4;
      }
      if (tmp < mi) {
        mi = tmp;
        x = i;
        y = j;
      }
    }
  }
}

int aaa[3][3];
bool CanErase(vector<vector<int>>& g, int x, int y)
{
  int cnt = 0;
  for (int i = 0; i < 4; i++) {
    int nx = x + dx[i];
    int ny = y + dy[i];
    if (nx < 0 || n <= nx || ny < 0 || n <= ny) {
    }
    else if (g[nx][ny] == 1) {
    }
    else {
      cnt++;
    }
  }

  if (cnt == 0) {
    return true;
  }

  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      int nx = x + i;
      int ny = y + j;
      if (nx < 0 || n <= nx || ny < 0 || n <= ny) {
        aaa[i + 1][j + 1] = 0;
      }
      else if (g[nx][ny] == 1) {
        aaa[i + 1][j + 1] = 0;
      }
      else {
        aaa[i + 1][j + 1] = 1;
      }
    }
  }

  cnt -= aaa[0][0] * aaa[0][1] * aaa[1][0];
  cnt -= aaa[2][0] * aaa[1][0] * aaa[2][1];
  cnt -= aaa[0][2] * aaa[0][1] * aaa[1][2];
  cnt -= aaa[2][2] * aaa[2][1] * aaa[1][2];

  return cnt == 1;
}

int Method1(vector<vector<int>>& g, vector<P>& ans, const vector<double>& thres, vector<int>& vec, int turn)
{
  double score = 0;
  double survive = n * n - m;

  vector<vector<double>> cnt(n, vector<double>(n, 0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (g[i][j] == 0) {
        cnt[i][j] = 1;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      canErase[i][j] = CanErase(g, i, j);
    }
  }

  vector<vector<double>> cnt2(n, vector<double>(n, 0));
  vector<vector<double>> cnt3(n, vector<double>(n, 0));
  vector<vector<double>> cnt4(n, vector<double>(n, 0));
  vector<vector<double>> cnt5(n, vector<double>(n, 0));

  int q = n * n - m;
  for (int ttt = 0; ttt < turn; ttt++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cnt2[i][j] = 0;
      }
    }

    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          cnt3[i][j] = 0;
        }
      }

      Move(g, cnt, cnt3, k);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          cnt2[i][j] += cnt3[i][j] / 4.0;
        }
      }
    }

    int x = -1;
    int y = -1;

    Search2(g, cnt2, cnt4, cnt5, x, y, thres[1], ttt);

    if (x == -1) {
      double mi1 = 1e9;
      double ma2 = -1e9;

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (g[i][j] == 1) {
            continue;
          }
          double keep = cnt2[i][j];
          if (cnt2[i][j] < thres[0]) {
            cnt2[i][j] = 0.0;
          }
          if (cnt2[i][j] < mi1) {
            x = i;
            y = j;
            mi1 = cnt2[i][j];
            ma2 = 0;
            int cntb = 0;
            for (int k = 0; k < 4; k++) {
              int nx = i + dx[k];
              int ny = j + dy[k];
              if (nx < 0 || n <= nx || ny < 0 || n <= ny) {
                cntb++;
              }
              else if (g[nx][ny] == 1) {
                cntb++;
              }
              else {
                ma2 += cnt2[nx][ny];
              }
            }
            //ma2 += cntb * 1024;
            for (int k = 0; k < 5; k++) {
              if (cntb == k && vec[k] == 1) {
                ma2 += 1e5;
              }
            }
            if (ttt < tttThres && canErase[i][j] == 0) {
              ma2 -= 1e6;
            }
          }
          else if (cnt2[i][j] == mi1) {
            double tmp = 0;
            int cntb = 0;
            for (int k = 0; k < 4; k++) {
              int nx = i + dx[k];
              int ny = j + dy[k];
              if (nx < 0 || n <= nx || ny < 0 || n <= ny) {
                cntb++;
              }
              else if (g[nx][ny] == 1) {
                cntb++;
              }
              else {
                tmp += cnt2[nx][ny];
              }
            }
            //tmp += cntb * 1024;
            for (int k = 0; k < 5; k++) {
              if (cntb == k && vec[k] == 1) {
                tmp += 1e5;
              }
            }
            if (ttt < tttThres && canErase[i][j] == 0) {
              tmp -= 1e6;
            }
            if (tmp > ma2) {
              x = i;
              y = j;
              mi1 = cnt2[i][j];
              ma2 = tmp;
            }
          }
          cnt2[i][j] = keep;
        }
      }
      //cout << ma2 << endl;
      //cout << x << ' ' << y << endl;
    }

    ans.push_back(P(x, y));
    survive -= cnt2[x][y];
    score += survive;
    cnt2[x][y] = 0;
    g[x][y] = 1;
    cnt = cnt2;

    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        int nx = x + i;
        int ny = y + j;
        if (nx < 0 || n <= nx || ny < 0 || n <= ny) {
          continue;
        }
        canErase[nx][ny] = CanErase(g, nx, ny);
      }
    }
  }

  vector<double> nokori;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (g[i][j] == 0) {
        nokori.push_back(cnt[i][j]);
      }
    }
  }
  sort(nokori.begin(), nokori.end());

  for (int i = 0; i < nokori.size(); i++) {
    survive -= nokori[i];
    score += survive;
  }

  score /= n * n - m;
  int scoreInt = round(1e6 * score / (n * n - m - 1));
  return scoreInt;
}

int Method2(vector<vector<int>>& g, vector<P>& ans, double thre, vector<int>& vec, int turn)
{
  int q = n * n - m;
  for (int _ = 0; _ < q; _++) {
    int ok = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (ok == 0 && g[i][j] == 0) {
          ans.push_back(P(i, j));
          ok = 1;
          g[i][j] = 1;
        }
      }
    }
  }
  return 11001001;
}

vector<vector<int>> transpose(const vector<vector<int>>& mat)
{
  if (mat.empty() || mat[0].empty()) return {};

  int H = mat.size();       // 元の行数
  int W = mat[0].size();    // 元の列数
  vector<vector<int>> res(W, vector<int>(H));

  for (int i = 0; i < H; ++i)
    for (int j = 0; j < W; ++j)
      res[j][i] = mat[i][j];

  return res;
}

ll solve_case(int case_num)
{
  timer.start();

  vector<P> best_ans;
  int bestScore = 0;
  vector<double> bestThres;
  vector<int> bestVec;

  auto g = input_data(case_num);

  int loop = 0;
  while (true) {
    if (timer.get_elapsed_time() > TIME_LIMIT) {
      break;
    }
    loop++;
    g = init_grid;
    vector<double> thres;
    thres.push_back(log_uniform(1e-10, 1.0));
    thres.push_back(log_uniform(1e-10, 1.0));
    if (timer.get_elapsed_time() > TIME_LIMIT > 2) {
      thres.clear();
      for (int k = 0; k < 2; k++) {
        thres.push_back(log_uniform(bestThres[k] / 10, bestThres[k] * 10));
      }
    }
    vector<int> vec;
    for (int i = 0; i < 5; i++) {
      vec.push_back(rand_xorshift() % 2);
      //if (i == 0) {
      //  vec[i] = 1;
      //}
    }
    int tenti = rand_xorshift() % 2;
    if (tenti) {
      g = transpose(g);
    }
    vector<P> ans;
    int score = Method1(g, ans, thres, vec, n * n - m);
    //int score = Method2(g, ans, thre, vec, 300);
    if (score > bestScore) {
      bestScore = score;
      best_ans = ans;
      bestThres = thres;
      bestVec = vec;
      if (tenti) {
        for (int i = 0; i < best_ans.size(); i++) {
          swap(best_ans[i].first, best_ans[i].second);
        }
      }
    }
  }

  //g = init_grid;
  //best_ans.clear();
  //bestScore = Method1(g, best_ans, bestThre, bestVec, n * n - m);

  cerr << loop << ' ';
  for (int k = 0; k < 2; k++) {
    cerr << bestThres[k] << ' ';
  }
  for (int k = 0; k < 5; k++) {
    cerr << bestVec[k];
  }
  cerr << endl;

  output_data(case_num, best_ans);
  return bestScore;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
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
      if (score > 990000) {
        //break;
      }
    }
  }

  return 0;
}
