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

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

// 乱数
namespace
{
  static uint32_t rand_xorshift() {
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

  static double rand_01() {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r) {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r) {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n) {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

const double TIME_LIMIT = 1.8;
int exec_mode;

const int n = 36;
const int m = 12;
const int L = 1000000;

string s[n];
vector<int> sv[n];
int p[n];

int c[m];

int a[m][m];
int current_score;

int best_a[m][m];
int best_score;

void store_best_score() {
  best_score = current_score;
  rep(i, m) {
    rep(j, m) {
      best_a[i][j] = a[i][j];
    }
  }
}

void restore_best_score() {
  current_score = best_score;
  rep(i, m) {
    rep(j, m) {
      a[i][j] = best_a[i][j];
    }
  }
}

bool is_out_of_range(int x, int y) {
  //if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

void initialize_state() {
  current_score = 0;
}

void input_data(int case_num) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (exec_mode == 0 || !ifs.is_open()) {
    // 標準入力
    int _n, _m, _L;
    cin >> _n >> _m >> _L;
    rep(i, n) {
      cin >> s[i];
      cin >> p[i];
    }
  }
  else {
    // ファイル入力
    int _n, _m, _L;
    ifs >> _n >> _m >> _L;
    rep(i, n) {
      ifs >> s[i];
      ifs >> p[i];
    }
  }

  rep(i, n) {
    vector<int> v;
    rep(j, s[i].size()) {
      if (s[i][j] == ' ') continue;
      v.push_back(s[i][j] - 'a');
    }
    sv[i] = v;
  }

  vector<pair<int, vector<int>>> tmp;
  rep(i, n) {
    tmp.push_back(make_pair(p[i], sv[i]));
  }
  sort(tmp.begin(), tmp.end());
  rep(i, n) {
    p[i] = tmp[i].first;
    sv[i] = tmp[i].second;
  }
}

void open_ofs(int case_num, ofstream& ofs) {
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

inline long double prob_one_or_more(long double k) {
  constexpr long double ONE_MILLION = 1'000'000.0L;
  if (k <= 0.0L)                  return 0.0L;   // 起こらない
  if (k >= ONE_MILLION - 1)       return 1.0L;   // ほぼ確実に起こる

  long double p = k / ONE_MILLION;               // 1 試行あたり成功確率
  // Q = 1 - (1-p)^1e6  = 1 - exp(1e6 * log(1-p))
  long double ln_term = ONE_MILLION * log1pl(-p);
  return 1.0L - expl(ln_term);
}

int calculate_score(bool all = false) {
  int visited[m] = {};

  int next[m][100];
  rep(i, m) {
    rep(j, 100) {
      next[i][j] = 0;
    }
    int sum = 0;
    rep(j, m) {
      rep(k, a[i][j]) {
        next[i][sum] = j;
        sum++;
      }
    }
  }

  int now = 0;
  rep(i, 10000) {
    visited[now]++;
    int ra = rand_xorshift() % 100;
    now = next[now][ra];
  }
  rep(i, m) {
    visited[i] = visited[i] * 100;
  }

  double res = 0;

  rep(i, n) {
    double kitai = 0;

    rep(k, 2) {
      double prob = 0;
      rep(j, sv[i].size()) {
        if (j == 0) {
          prob = visited[sv[i][j] + k * 6];
        }
        else {
          prob *= a[sv[i][j - 1] + k * 6][sv[i][j] + k * 6] / 100.0;
        }
      }
      kitai += prob;
    }

    double kaku = prob_one_or_more(kitai);
    res += kaku * p[i];
    if (all) {
      cout << i << " " << kitai << " " << kaku << " " << p[i] << endl;
    }
  }

  return round(res);
}

int calculate_score_12(bool all = false, int siki = 0) {
  int visited[m] = {};

  int next[m][100];
  rep(i, m) {
    rep(j, 100) {
      next[i][j] = 0;
    }
    int sum = 0;
    rep(j, m) {
      rep(k, a[i][j]) {
        next[i][sum] = j;
        sum++;
      }
    }
  }

  int now = 0;
  rep(i, 1000) {
    visited[now]++;
    int ra = rand_xorshift() % 100;
    now = next[now][ra];
  }
  rep(i, m) {
    visited[i] = visited[i] * 1000;
  }

  double res = 0;

  rep(i, n) {
    double dp[2] = {};
    double dp2[2] = {};

    rep(j, sv[i].size()) {
      int num = sv[i][j];
      if (j == 0) {
        dp[0] = visited[num];
        dp[1] = visited[num + 6];
      }
      else {
        dp2[0] = 0;
        dp2[1] = 0;

        dp2[0] += dp[0] * a[sv[i][j - 1]][sv[i][j]] / 100.0;
        dp2[1] += dp[0] * a[sv[i][j - 1]][sv[i][j] + 6] / 100.0;
        dp2[0] += dp[1] * a[sv[i][j - 1] + 6][sv[i][j]] / 100.0;
        dp2[1] += dp[1] * a[sv[i][j - 1] + 6][sv[i][j] + 6] / 100.0;

        dp[0] = dp2[0];
        dp[1] = dp2[1];
      }
    }

    double kaku = prob_one_or_more(dp[0] + dp[1]);
    if (siki == 0) {
      res += kaku * pow(max(p[i] - 1, 0), 1.0);
    }
    else {
      res += kaku * p[i];
    }
    if (all) {
      cout << i << " " << dp[0] + dp[1] << " " << kaku << " " << p[i] << endl;
    }
  }

  return round(res);
}

void output_data(ofstream& ofs) {
  if (exec_mode == 0) {
    // 標準出力
    rep(i, m) {
      cout << (char)('a' + c[i]);
      rep(j, m) {
        cout << ' ' << a[i][j];
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    rep(i, m) {
      ofs << (char)('a' + c[i]);
      rep(j, m) {
        ofs << ' ' << a[i][j];
      }
      ofs << endl;
    }
  }
}

void build_initial_solution() {
  rep(i, m) {
    c[i] = i % 6;
  }

  rep(i, m) {
    rep(j, m) {
      a[i][j] = 0;
    }
  }

  rep(i, m) {
    rep(j, 100) {
      a[i][j % m]++;
    }
  }
}

void build_initial_solution_2() {
  rep(i, m) {
    c[i] = i % 6;
  }

  rep(i, m) {
    rep(j, m) {
      a[i][j] = 0;
    }
  }

  int cnt[6][6];
  rep(i, 6) {
    rep(j, 6) {
      cnt[i][j] = 0;
    }
  }

  rep(i, m) {
    rep(j, s[i].size() - 1) {
      int x = s[i][j] - 'a';
      int y = s[i][j + 1] - 'a';
      cnt[x][y]++;
    }
  }

  rep(i, m) {
    double sum = 0;
    rep(j, 6) {
      sum += cnt[i % 6][j];
    }
    int sum2 = 0;
    rep(j, m) {
      a[i][j] = round(cnt[i % 6][j % 6] / sum * 100.0 / 2);
      sum2 += a[i][j];
    }
    if (sum2 > 100) {
      int diff = sum2 - 100;
      while (diff > 0) {
        int x = rand_xorshift() % m;
        if (a[i][x] > 0) {
          a[i][x]--;
          diff--;
        }
      }
    }
    else if (sum2 < 100) {
      int diff = 100 - sum2;
      while (diff > 0) {
        int x = rand_xorshift() % m;
        a[i][x]++;
        diff--;
      }
    }
  }
}

void build_initial_solution_3() {
  rep(i, m) {
    c[i] = i % 6;
  }

  rep(i, m) {
    rep(j, m) {
      a[i][j] = 0;
    }
  }

  rep(i, sv[n - 1].size()) {
    c[i] = sv[n - 1][i];
  }
  rep(i, m) {
    rep(j, m) {
      a[i][j] = 0;
    }
    a[i][(i + 1) % m] = 40;
    rep(j, 60) {
      a[i][(i + j + 1) % m]++;
    }
  }
}

vector<int> make_99(vector<double> v) {
  double sum = 0;
  rep(i, 6) {
    sum += v[i];
  }
  vector<int> res(6);
  if (sum < 0.1) {
    rep(i, 99) {
      res[i % 6]++;
    }
    return res;
  }

  rep(i, 6) {
    res[i] = round(v[i] / sum * 99.0);
  }
  int sum2 = 0;
  rep(i, 6) {
    sum2 += res[i];
  }
  if (sum2 > 99) {
    int diff = sum2 - 99;
    while (diff > 0) {
      int x = rand_xorshift() % 6;
      if (res[x] > 0) {
        res[x]--;
        diff--;
      }
    }
  }
  else if (sum2 < 99) {
    int diff = 99 - sum2;
    while (diff > 0) {
      int x = rand_xorshift() % 6;
      res[x]++;
      diff--;
    }
  }
  return res;
}

vector<int> make_100_12(vector<double> v) {
  double sum = 0;
  rep(i, 12) {
    sum += v[i];
  }
  vector<int> res(12);
  if (sum < 0.1) {
    rep(i, 100) {
      res[i % 12]++;
    }
    return res;
  }

  rep(i, 12) {
    res[i] = round(v[i] / sum * 100.0);
  }
  int sum2 = 0;
  rep(i, 12) {
    sum2 += res[i];
  }
  if (sum2 > 100) {
    int diff = sum2 - 100;
    while (diff > 0) {
      int x = rand_xorshift() % 12;
      if (res[x] > 0) {
        res[x]--;
        diff--;
      }
    }
  }
  else if (sum2 < 100) {
    int diff = 100 - sum2;
    while (diff > 0) {
      int x = rand_xorshift() % 12;
      res[x]++;
      diff--;
    }
  }
  return res;
}


void build_initial_solution_4() {
  int max_a[m][m];
  int max_score = -1;

  rep(i, m) {
    c[i] = i % 6;
  }

  rep(sss, 10) {
    rep(i, m) {
      rep(j, m) {
        a[i][j] = 0;
      }
      a[i][(i + 6) % m] = 1;
    }

    // 一番大きいやつ
    rep(kk, 2) {
      vector<double> cnt_sum[6];
      rep(i, 6) {
        cnt_sum[i] = vector<double>(6);
      }
      rep(l, sss + 1) {
        int lll = n - 1 - (kk + l * 2);
        rep(i, sv[lll].size() - 1) {
          int x = sv[lll][i];
          int y = sv[lll][i + 1];
          cnt_sum[x][y] += p[lll];
        }
      }
      int ma = 99;
      rep(i, 6) {
        auto v = make_99(cnt_sum[i]);
        rep(j, 6) {
          a[kk * 6 + i][kk * 6 + j] = v[j];
        }
      }
    }
    int score = calculate_score();
    if (score > max_score) {
      max_score = score;
      rep(i, m) {
        rep(j, m) {
          max_a[i][j] = a[i][j];
        }
      }
    }
  }

  int iter = 0;
  while (true) {
    int ra = rand_xorshift() % 10 + 2;
    vector<int> v[2];
    rep(i, ra) {
      if (rand_xorshift() % 2 == 0) {
        v[0].push_back(n - 1 - i);
      }
      else {
        v[1].push_back(n - 1 - i);
      }
    }
    if (v[0].empty() || v[1].empty()) continue;


    if (get_elapsed_time() > 0.5) {
      break;
    }
    iter++;

    rep(i, m) {
      rep(j, m) {
        a[i][j] = 0;
      }
      a[i][(i + 6) % m] = 1;
    }

    rep(kk, 2) {
      vector<double> cnt_sum[6];
      rep(i, 6) {
        cnt_sum[i] = vector<double>(6);
      }
      rep(l, v[kk].size()) {
        int lll = v[kk][l];
        rep(i, sv[lll].size() - 1) {
          int x = sv[lll][i];
          int y = sv[lll][i + 1];
          cnt_sum[x][y] += p[lll];
        }
      }
      int ma = 99;
      rep(i, 6) {
        auto v = make_99(cnt_sum[i]);
        rep(j, 6) {
          a[kk * 6 + i][kk * 6 + j] = v[j];
        }
      }
    }
    int score = calculate_score();
    if (score > max_score) {
      max_score = score;
      rep(i, m) {
        rep(j, m) {
          max_a[i][j] = a[i][j];
        }
      }
    }
  }

  cerr << iter << ' ' << max_score << ' ' << get_elapsed_time() << endl;

  rep(i, m) {
    rep(j, m) {
      a[i][j] = max_a[i][j];
    }
  }
}

void build_initial_solution_5(double time_limit) {
  int max_a[m][m];
  int max_score = -1;

  rep(i, m) {
    c[i] = i % 6;
  }

  rep(sss, 10) {
    rep(i, m) {
      rep(j, m) {
        a[i][j] = 0;
      }
      a[i][(i + 6) % m] = 1;
    }

    // 一番大きいやつ
    rep(kk, 2) {
      vector<double> cnt_sum[6];
      rep(i, 6) {
        cnt_sum[i] = vector<double>(6);
      }
      rep(l, sss + 1) {
        int lll = n - 1 - (kk + l * 2);
        rep(i, sv[lll].size() - 1) {
          int x = sv[lll][i];
          int y = sv[lll][i + 1];
          cnt_sum[x][y] += p[lll];
        }
      }
      int ma = 99;
      rep(i, 6) {
        auto v = make_99(cnt_sum[i]);
        rep(j, 6) {
          a[kk * 6 + i][kk * 6 + j] = v[j];
        }
      }
    }
    int score = calculate_score();
    if (score > max_score) {
      max_score = score;
      rep(i, m) {
        rep(j, m) {
          max_a[i][j] = a[i][j];
        }
      }
    }
  }

  int iter = 0;
  while (true) {

    if (get_elapsed_time() > time_limit) {
      break;
    }
    iter++;

    rep(i, m) {
      rep(j, m) {
        a[i][j] = 0;
      }
    }

    int ra = rand_xorshift() % 10 + 2;
    vector<int> ras;
    rep(i, ra) {
      int lll = n - 1 - i;
      if (rand_xorshift() % 100 < 90) {
        ras.push_back(lll);
      }
    }

    vector<double> cnt_sum[12];
    rep(i, 12) {
      cnt_sum[i] = vector<double>(12);
    }
    rep(l, ras.size()) {
      int lll = ras[l];
      vector<int> vv;
      rep(i, sv[lll].size()) {
        vv.push_back(sv[lll][i] + rand_xorshift() % 2 * 6);
      }
      rep(i, sv[lll].size() - 1) {
        int x = vv[i];
        int y = vv[i + 1];
        cnt_sum[x][y] += p[lll];
      }
    }
    int ma = 99;
    rep(i, 12) {
      auto v = make_100_12(cnt_sum[i]);
      rep(j, 12) {
        a[i][j] = v[j];
      }
    }

    int score = calculate_score_12();
    if (score > max_score) {
      max_score = score;
      rep(i, m) {
        rep(j, m) {
          max_a[i][j] = a[i][j];
        }
      }
    }
  }

  cerr << iter << ' ' << max_score << ' ' << get_elapsed_time() << endl;

  rep(i, m) {
    rep(j, m) {
      a[i][j] = max_a[i][j];
    }
  }
}


struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams, int siki, double time_limit) {
  int sum = 0;
  rep(i, 10) {
    sum += calculate_score_12();
  }
  current_score = sum / 10;
  store_best_score();

  double start_time = get_elapsed_time();
  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) break;
    }

    //if (current_score * 1.2 < best_score) {
    //  restore_best_score();
    //  continue;
    //}

    double progress_ratio = (now_time - start_time) / (time_limit - start_time);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[2];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra5 = rand_xorshift() % 15 + 1;
      ra1 = rand_xorshift() % 2;
      ra2 = rand_xorshift() % 6;
      ra3 = rand_xorshift() % 6;
      int ng = 0;
      while (a[ra2 + ra1 * 6][ra3 + ra1 * 6] < ra5) {
        ra2 = rand_xorshift() % 6;
        ra3 = rand_xorshift() % 6;
        ng++;
        if (ng > 100)break;
      }
      if (ng > 100)continue;
      ra4 = rand_xorshift() % 6;
      while (ra3 == ra4) {
        ra4 = rand_xorshift() % 6;
      }

      a[ra2 + ra1 * 6][ra3 + ra1 * 6] -= ra5;
      a[ra2 + ra1 * 6][ra4 + ra1 * 6] += ra5;
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
      ra1 = 0;
      ra5 = rand_xorshift() % 15 + 1;
      ra2 = rand_xorshift() % 12;
      ra3 = rand_xorshift() % 12;
      int ng = 0;
      while (a[ra2][ra3] < ra5) {
        ra2 = rand_xorshift() % 12;
        ra3 = rand_xorshift() % 12;
        ng++;
        if (ng > 100)break;
      }
      if (ng > 100)continue;
      ra4 = rand_xorshift() % 12;
      while (ra3 == ra4) {
        ra4 = rand_xorshift() % 12;
      }

      a[ra2][ra3] -= ra5;
      a[ra2][ra4] += ra5;
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
      // 近傍操作2
      ra1 = rand_xorshift() % 6;
      ra2 = rand_xorshift() % 6;
      swap(a[ra1][ra2], a[ra1][ra2 + 6]);
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[3]) {
      // 近傍操作2
      ra1 = rand_xorshift() % 6;
      rep(j, 12) {
        swap(a[ra1][j], a[ra1 + 6][j]);
      }
    }

    // スコア計算
    double tmp_score = calculate_score_12(false, siki);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score > best_score) {
        int sum = 0;
        rep(i, 10) {
          sum += calculate_score_12(false, siki);
        }
        current_score = sum / 10;
        if (current_score > best_score) {
          store_best_score();
        }
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        a[ra2 + ra1 * 6][ra3 + ra1 * 6] += ra5;
        a[ra2 + ra1 * 6][ra4 + ra1 * 6] -= ra5;
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
        a[ra2][ra3] += ra5;
        a[ra2][ra4] -= ra5;
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
        // 近傍操作2 の巻き戻し
        swap(a[ra1][ra2], a[ra1][ra2 + 6]);
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[3]) {
        // 近傍操作2 の巻き戻し
        rep(j, 12) {
          swap(a[ra1][j], a[ra1 + 6][j]);
        }
      }
    }
  }

  cerr << loop << endl;

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams) {
  start_timer();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  //build_initial_solution();
  //build_initial_solution_2();
  //build_initial_solution_3();
  //build_initial_solution_4();
  build_initial_solution_5(1.0);

  // 焼きなまし実行
  run_simulated_annealing(annealingParams, 0, 1.5);

  annealingParams.start_temperature[0] = 5000048.0;
  run_simulated_annealing(annealingParams, 1, 1.9);

  // 解答を出力
  output_data(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (true) {
    score = calculate_score_12(false, 1);
  }
  return score;
}

int main() {
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 2000048.0;
  annealingParams.start_temperature[1] = 2048.0;
  annealingParams.start_temperature[2] = 2048.0;
  annealingParams.start_temperature[3] = 2048.0;
  annealingParams.start_temperature[4] = 2048.0;
  annealingParams.start_temperature[5] = 2048.0;
  annealingParams.start_temperature[6] = 2048.0;
  annealingParams.start_temperature[7] = 2048.0;
  annealingParams.start_temperature[8] = 2048.0;
  annealingParams.start_temperature[9] = 2048.0;
  annealingParams.end_temperature = 0.0;
  annealingParams.score_scale = 12345.0;
  annealingParams.operation_thresholds[0] = 100;
  annealingParams.operation_thresholds[1] = 200;
  annealingParams.operation_thresholds[2] = 300;
  annealingParams.operation_thresholds[3] = 400;
  annealingParams.operation_thresholds[4] = 500;
  annealingParams.operation_thresholds[5] = 600;
  annealingParams.operation_thresholds[6] = 700;
  annealingParams.operation_thresholds[7] = 800;
  annealingParams.operation_thresholds[8] = 900;
  annealingParams.operation_thresholds[9] = 1000;

  if (exec_mode == 0) {
    cerr << solve_case(0, annealingParams) << endl;
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    srep(i, 0, 150) {
      ll score = solve_case(i, annealingParams);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << get_elapsed_time() << ", "
          << endl;
      }
    }
  }
  else if (exec_mode == 3) {
    int loop_count = 0;
    AnnealingParams best_annealingParams;
    ll best_sum_score = 0;

    while (true) {
      AnnealingParams new_annealingParams;
      new_annealingParams.start_temperature[0] = pow(2.0, rand_01() * 20);
      new_annealingParams.end_temperature = 0.0;
      new_annealingParams.score_scale = pow(2.0, rand_01() * 20);
      new_annealingParams.operation_thresholds[0] = rand() % 101;

      ll sum_score = 0;
      srep(i, 0, 150) {
        ll score = solve_case(i, new_annealingParams);
        sum_score += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cerr << "loop_count = " << loop_count
        << ", sum_score = " << sum_score
        << ", start_temperature = " << new_annealingParams.start_temperature[0]
        << ", end_temperature = " << new_annealingParams.end_temperature
        << ", score_scale = " << new_annealingParams.score_scale
        << ", operation_thresholds = " << new_annealingParams.operation_thresholds[0]
        << endl;

      if (sum_score > best_sum_score) {
        best_sum_score = sum_score;
        best_annealingParams = new_annealingParams;
      }

      loop_count++;
    }
  }

  return 0;
}
