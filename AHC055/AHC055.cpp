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

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

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

// 1次元キューのクラス
class Queue1D
{
private:
  static const int MAX_SIZE = 500;
  int arr[MAX_SIZE];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue1D() : head(0), tail(0) {}

  void clear_queue()
  {
    head = 0;
    tail = 0;
  }

  int front() const
  {
    return arr[head];
  }

  void push(int val)
  {
    arr[tail] = val;
    tail++;
  }

  void pop()
  {
    head++;
  }

  int size() const
  {
    return tail - head;
  }

  bool empty() const
  {
    return head == tail;
  }
};

Queue1D que;

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

struct CoolingParams
{
  double start_temp;   // START_TEMP
  double end_temp;     // END_TEMP
  double poly_p;       // 多項式冷却のp (例:2.0)
  double exp_k;        // 指数冷却のk (例:5.0)
  double phase_split;  // フェーズ切り替え割合 (例:0.7)
};

enum class CoolingType
{
  Linear,
  Poly,        // 多項式 (pow)
  Exponential, // 指数で後半ガッと冷やす
  Cosine,
  TwoPhase     // 前半ほぼ一定→後半で一気に落とす
};

// progress_ratio は 0.0 ～ 1.0
inline double get_temp(double progress_ratio,
  CoolingType type,
  const CoolingParams& p)
{
  // 安全のためクランプ
  progress_ratio = std::clamp(progress_ratio, 0.0, 1.0);

  double t01 = 1.0; // 1.0→0.0 に落ちる係数をここで作って最後に補間

  switch (type) {
  case CoolingType::Linear:
    // 線形: t01 = 1 - ratio
    t01 = 1.0 - progress_ratio;
    break;

  case CoolingType::Poly:
    // 多項式: t01 = (1 - ratio)^p
    // p.poly_p >1 で後半急降下、<1で序盤急降下
    t01 = std::pow(1.0 - progress_ratio, p.poly_p);
    break;

  case CoolingType::Exponential:
    // 指数: t01 = exp(-k * ratio)
    // kが大きいと最後に一気に冷える
    t01 = std::exp(-p.exp_k * progress_ratio);
    break;

  case CoolingType::Cosine:
    // コサイン: t01 = (1 + cos(pi * ratio)) / 2
    // なめらかに 1→0
    t01 = (1.0 + std::cos(3.1415926535 * progress_ratio)) * 0.5;
    break;

  case CoolingType::TwoPhase:
    // 2フェーズ:
    // ratio < split までは一定(=ほぼstart_temp相当の高さ)
    // ratio >= split から線形に 1→0
  {
    double split = p.phase_split; // 例:0.7
    if (progress_ratio < split) {
      t01 = 1.0; // ずっと高温
    }
    else {
      double r2 = (progress_ratio - split) / (1.0 - split); // 0→1
      // r2=0 のとき t01=1, r2=1 のとき t01=0
      t01 = 1.0 - r2;
    }
  }
  break;
  }

  // t01 = 1.0 のとき start_temp
  // t01 = 0.0 のとき end_temp
  double temp = p.end_temp + (p.start_temp - p.end_temp) * t01;
  return temp;
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

int exec_mode;

const int n = 200;
const int MAX_C = 6;
const int MAX_H = 500;

int h[n];
int initial_h[n];
int c[n];
int initial_c[n];
int a[n][n];
int sort_a_idx[n][n];
int sum_h;

int score;
int ans_count;
int w[200000], b[200000];

int best_score;
int best_ans_count;
int best_w[200000], best_b[200000];

int best_best_score;
int best_best_ans_count;
int best_best_w[200000], best_best_b[200000];

void store_best_score()
{
  best_score = score;
  best_ans_count = ans_count;
  for (int i = 0; i < ans_count; i++) {
    best_w[i] = w[i];
    best_b[i] = b[i];
  }
}

void restore_best_score()
{
  score = best_score;
  ans_count = best_ans_count;
  for (int i = 0; i < ans_count; i++) {
    w[i] = best_w[i];
    b[i] = best_b[i];
  }
}

void store_best_best_score()
{
  best_best_score = best_score;
  best_best_ans_count = best_ans_count;
  for (int i = 0; i < best_ans_count; i++) {
    best_best_w[i] = best_w[i];
    best_best_b[i] = best_b[i];
  }
}

void restore_best_best_score()
{
  best_score = best_best_score;
  best_ans_count = best_best_ans_count;
  for (int i = 0; i < best_ans_count; i++) {
    best_w[i] = best_best_w[i];
    best_b[i] = best_best_b[i];
  }
}

void initialize_state()
{
  for (int i = 0; i < n; i++) {
    initial_h[i] = h[i];
    initial_c[i] = c[i];
  }
}

void reset_state()
{
  for (int i = 0; i < n; i++) {
    h[i] = initial_h[i];
    c[i] = initial_c[i];
  }
}

std::istream& open_input_stream(int case_num, std::ifstream& fin, int exec_mode)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    fin.open(oss.str());
    if (fin.is_open()) {
      return fin;
    }
  }
  return std::cin;
}

std::ostream& open_output_stream(int case_num, std::ofstream& fout, int exec_mode)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    fout.open(oss.str());
    if (fout.is_open()) {
      return fout;
    }
  }
  return std::cout;
}

void read_case(std::istream& is)
{
  int _n;
  is >> _n;
  for (int i = 0; i < n; i++) {
    is >> h[i];
  }
  for (int i = 0; i < n; i++) {
    is >> c[i];
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      is >> a[i][j];
    }
  }

  sum_h = 0;
  for (int i = 0; i < n; i++) {
    sum_h += h[i];
  }

  for (int i = 0; i < n; i++) {
    initial_h[i] = h[i];
    initial_c[i] = c[i];
  }

  // sort_a_idxの初期化
  for (int i = 0; i < n; i++) {
    vector<P> vec;
    for (int j = 0; j < n; j++) {
      vec.push_back(P(a[i][j], j));
    }
    sort(vec.begin(), vec.end(), greater<P>());
    for (int j = 0; j < n; j++) {
      sort_a_idx[i][j] = vec[j].second;
    }
  }
}

void write_solution(std::ostream& os)
{
  for (int i = 0; i < ans_count; ++i) {
    os << w[i] << ' ' << b[i] << '\n';
  }
}

int calculate_score()
{
  int res = sum_h - ans_count + 1;
  return res;
}

void build_initial_solution()
{
  ans_count = 0;

  // 残りを武器で殴る
  while (true) {
    int max_effect = -1;
    int use_w = -1;
    int use_b = -1;
    for (int i = 0; i < n; i++) {
      if (h[i] > 0) {
        continue;
      }
      if (c[i] == 0) {
        continue;
      }
      for (int j = 0; j < n; j++) {
        if (h[j] <= 0) {
          continue;
        }
        int tmp_effect = a[i][j];
        if (a[i][j] * c[i] >= h[j]) {
          tmp_effect = h[j] + 10000;
        }
        if (tmp_effect > max_effect) {
          max_effect = tmp_effect;
          use_w = i;
          use_b = j;
        }
      }
    }

    if (max_effect < 0) {
      int min_h = INT_INF;
      int min_i = -1;
      for (int i = 0; i < n; i++) {
        if (h[i] <= 0) {
          continue;
        }
        if (h[i] < min_h) {
          min_h = h[i];
          min_i = i;
        }
      }
      use_w = -1;
      use_b = min_i;
    }

    if (use_b == -1) {
      break;
    }

    w[ans_count] = use_w;
    b[ans_count] = use_b;
    if (use_w == -1) {
      h[use_b]--;
    }
    else {
      h[use_b] -= a[use_w][use_b];
      c[use_w]--;
    }
    ans_count++;
  }

  score = calculate_score();
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

int best_sum;
int best_hukuzatsu;
int best_target[n][MAX_C];
int sum_attack[n];

int order[n];

int Sim(int sim_mode)
{
  ans_count = 0;
  reset_state();
  int mottainai = 0;

  // best_targetを解に変換
  ans_count = 0;

  int got_count = 0;
  que.clear_queue();
  int order_itr = 0;

  // まずは素手攻撃をすべて行う
  for (int i = 0; i < n; i++) {
    if (h[i] - sum_attack[i] > 0) {
      if (sim_mode == 1) {
        while (h[i] - sum_attack[i] > 0) {
          w[ans_count] = -1;
          b[ans_count] = i;
          h[i]--;
          ans_count++;
        }
      }
      else {
        ans_count += h[i] - sum_attack[i];
        h[i] = sum_attack[i];
      }
    }
    if (h[i] <= 0) {
      got_count++;
      que.push(i);
    }
  }

  // 武器攻撃を行う
  // 残りを武器で殴る
  while (got_count < n) {
    if (que.empty()) {
      // 仕方がないので素手攻撃
      while (true) {
        int x = order[order_itr];
        if (h[x] <= 0) {
          order_itr++;
          continue;
        }
        if (h[x] > 0) {
          if (sim_mode == 1) {
            while (h[x] > 0) {
              w[ans_count] = -1;
              b[ans_count] = x;
              h[x]--;
              ans_count++;
            }
          }
          else {
            ans_count += h[x];
            h[x] = 0;
          }
          got_count++;
          que.push(x);
          break;
        }
      }
      continue;
    }

    int cur = que.front();
    que.pop();
    while (c[cur] > 0) {
      int target = best_target[cur][c[cur] - 1];
      if (h[target] <= 0) {
        if (sim_mode == 1) {
          mottainai += a[cur][target];
        }
        // 次のターゲットへ
        // sorted_a_idxを使う
        for (int i = 0; i < n; i++) {
          int nxt = sort_a_idx[cur][i];
          if (h[nxt] > 0) {
            target = nxt;
            break;
          }
        }
      }
      if (h[target] <= 0) {
        break;
      }
      if (sim_mode == 1) {
        w[ans_count] = cur;
        b[ans_count] = target;
        if (a[cur][target] == 1) {
          w[ans_count] = -1;
        }
      }
      h[target] -= a[cur][target];
      c[cur]--;
      ans_count++;
      if (h[target] <= 0) {
        got_count++;
        que.push(target);
      }

      if (1 <= h[target] && h[target] <= 1) {
        if (sim_mode == 1) {
          while (h[target] > 0) {
            w[ans_count] = -1;
            b[ans_count] = target;
            h[target]--;
            ans_count++;
          }
        }
        else {
          ans_count += h[target];
          h[target] = 0;
        }
        got_count++;
        que.push(target);
        break;
      }
    }
  }

  if (sim_mode == 1) {
    cerr << "mottainai: " << mottainai << endl;
  }
  return  calculate_score();
}

inline bool sa_accept_fast(double deltaScore, double temp)
{
  if (deltaScore >= 0) return true;
  // u ∈ (0,1) → -log(u) は Exp(1) 乱数
  const double y = -std::log((rand_xorshift() + 0.5) * (1.0 / UINT_MAX));
  // 受理条件: -deltaScore <= temp * y  ≡  deltaScore > -temp*y
  return (-deltaScore) <= (temp * y);
}

void Method1(AnnealingParams annealingParams, AnnealingParams annealingParams2, double timeLimit1, double timeLimit2)
{
  reset_state();

  int target[n][MAX_C];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < c[i]; j++) {
      target[i][j] = rand_xorshift() % n;
      //target[i][j] = sort_a_idx[i][n-1];
      while (target[i][j] == i) {
        target[i][j] = rand_xorshift() % n;
      }
    }
  }

  // sum_attack初期化
  for (int i = 0; i < n; i++) {
    sum_attack[i] = 0;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < c[i]; j++) {
      sum_attack[target[i][j]] += a[i][target[i][j]];
    }
  }

  int current_sum = 0;
  for (int i = 0; i < n; i++) {
    current_sum += min(h[i], sum_attack[i]);
  }

  int current_hukuzatsu = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < c[i]; j++) {
      if (target[i][j] != target[i][(j + 1) % c[i]]) {
        current_hukuzatsu++;
      }
    }
  }

  best_sum = current_sum;
  best_hukuzatsu = current_hukuzatsu;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < c[i]; j++) {
      best_target[i][j] = target[i][j];
    }
  }

  // targetを
  {
    double now_time = timer.get_elapsed_time();
    const double START_TEMP = annealingParams.start_temperature[0];
    const double END_TEMP = annealingParams.end_temperature;

    CoolingParams params;
    params.start_temp = START_TEMP;
    params.end_temp = END_TEMP;
    params.poly_p = 2.0;
    params.exp_k = 5.0;
    params.phase_split = 0.7;

    int loop = 0;
    while (true) {
      loop++;

      if (loop % 100 == 0) {
        now_time = timer.get_elapsed_time();
        if (now_time > timeLimit1) { break; }
      }

      if (rand_xorshift() % 123456 == 0) {
        // bestに戻す
        current_sum = best_sum;
        current_hukuzatsu = best_hukuzatsu;
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < c[i]; j++) {
            target[i][j] = best_target[i][j];
          }
        }
        // sum_attack初期化
        for (int i = 0; i < n; i++) {
          sum_attack[i] = 0;
        }
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < c[i]; j++) {
            sum_attack[target[i][j]] += a[i][target[i][j]];
          }
        }
      }

      double progress_ratio = now_time / timeLimit1;
      double temp = get_temp(progress_ratio, CoolingType::Cosine, params);

      // 近傍解作成
      int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[1];
      int ra1, ra2, ra3, ra4, ra5;
      int keep1, keep2, keep3, keep4, keep5;

      int tmp_current_sum = current_sum;
      int tmp_current_hukuzatsu = current_hukuzatsu;

      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1
        ra1 = rand_xorshift() % n;
        ra2 = rand_xorshift() % c[ra1];
        keep1 = target[ra1][ra2];
        target[ra1][ra2] = sort_a_idx[ra1][rand_xorshift() % (n - 180)];
        while (target[ra1][ra2] == ra1) {
          target[ra1][ra2] = sort_a_idx[ra1][rand_xorshift() % (n - 180)];
        }

        // current_sum更新
        tmp_current_sum -= min(h[keep1], sum_attack[keep1]);
        tmp_current_sum -= min(h[target[ra1][ra2]], sum_attack[target[ra1][ra2]]);

        // sum_attack更新
        sum_attack[keep1] -= a[ra1][keep1];
        sum_attack[target[ra1][ra2]] += a[ra1][target[ra1][ra2]];

        // current_sum更新
        tmp_current_sum += min(h[keep1], sum_attack[keep1]);
        tmp_current_sum += min(h[target[ra1][ra2]], sum_attack[target[ra1][ra2]]);

        // hukuzatsu更新
        if (c[ra1] >= 2) {
          int before_idx = (ra2 - 1 + c[ra1]) % c[ra1];
          int after_idx = (ra2 + 1) % c[ra1];
          if (target[ra1][before_idx] != keep1) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra1][before_idx] != target[ra1][ra2]) {
            tmp_current_hukuzatsu++;
          }
          if (target[ra1][after_idx] != keep1) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra1][after_idx] != target[ra1][ra2]) {
            tmp_current_hukuzatsu++;
          }
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 スワップ
        ra1 = rand_xorshift() % n;
        ra2 = rand_xorshift() % c[ra1];
        ra3 = rand_xorshift() % n;
        ra4 = rand_xorshift() % c[ra3];
        keep1 = target[ra1][ra2];
        keep2 = target[ra3][ra4];
        target[ra1][ra2] = keep2;
        target[ra3][ra4] = keep1;

        // current_sum更新
        tmp_current_sum -= min(h[keep1], sum_attack[keep1]);
        tmp_current_sum -= min(h[keep2], sum_attack[keep2]);

        // sum_attack更新
        sum_attack[keep1] += a[ra3][keep1];
        sum_attack[keep1] -= a[ra1][keep1];
        sum_attack[keep2] += a[ra1][keep2];
        sum_attack[keep2] -= a[ra3][keep2];
        // current_sum更新
        tmp_current_sum += min(h[keep1], sum_attack[keep1]);
        tmp_current_sum += min(h[keep2], sum_attack[keep2]);

        // hukuzatsu更新
        // ra1側
        if (c[ra1] >= 2) {
          int before_idx = (ra2 - 1 + c[ra1]) % c[ra1];
          int after_idx = (ra2 + 1) % c[ra1];
          if (target[ra1][before_idx] != keep1) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra1][before_idx] != keep2) {
            tmp_current_hukuzatsu++;
          }
          if (target[ra1][after_idx] != keep1) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra1][after_idx] != keep2) {
            tmp_current_hukuzatsu++;
          }
        }
        // ra3側
        if (c[ra3] >= 2) {
          int before_idx = (ra4 - 1 + c[ra3]) % c[ra3];
          int after_idx = (ra4 + 1) % c[ra3];
          if (target[ra3][before_idx] != keep2) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra3][before_idx] != keep1) {
            tmp_current_hukuzatsu++;
          }
          if (target[ra3][after_idx] != keep2) {
            tmp_current_hukuzatsu--;
          }
          if (target[ra3][after_idx] != keep1) {
            tmp_current_hukuzatsu++;
          }
        }
      }

      //tmp_current_sum = min(tmp_current_sum, sum_h);

      // 焼きなましで採用判定
      const double hukuzatsu_penalty = 13.0;
      double diff_score = (tmp_current_sum - tmp_current_hukuzatsu * hukuzatsu_penalty - (current_sum - current_hukuzatsu * hukuzatsu_penalty)) * annealingParams.score_scale;
      if (sa_accept_fast(diff_score, temp)) {
        // 採用
        current_sum = tmp_current_sum;
        current_hukuzatsu = tmp_current_hukuzatsu;

        // ベスト更新
        if (current_sum > best_sum) {
          best_sum = current_sum;
          best_hukuzatsu = current_hukuzatsu;
          for (int i = 0; i < n; i++) {
            for (int j = 0; j < c[i]; j++) {
              best_target[i][j] = target[i][j];
            }
          }
        }
      }
      else {
        // 元に戻す
        if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
          // 近傍操作1 の巻き戻し

          // sum_attack更新
          sum_attack[keep1] += a[ra1][keep1];
          sum_attack[target[ra1][ra2]] -= a[ra1][target[ra1][ra2]];

          target[ra1][ra2] = keep1;
        }
        else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
          // 近傍操作2 の巻き戻し

          // sum_attack更新
          sum_attack[keep1] += a[ra1][keep1];
          sum_attack[keep1] -= a[ra3][keep1];
          sum_attack[keep2] += a[ra3][keep2];
          sum_attack[keep2] -= a[ra1][keep2];
          target[ra1][ra2] = keep1;
          target[ra3][ra4] = keep2;
        }
      }
    }

    if (exec_mode != 0 && exec_mode != 3) {
      cerr << loop << endl;
    }

    cerr << "Best Sum: " << best_sum << ", Best Hukuzatsu: " << best_hukuzatsu << ' ' << sum_h << endl;

    // best_targetをtargetに戻す
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < c[i]; j++) {
        target[i][j] = best_target[i][j];
      }
    }

    // sum_attack初期化
    for (int i = 0; i < n; i++) {
      sum_attack[i] = 0;
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < c[i]; j++) {
        sum_attack[target[i][j]] += a[i][target[i][j]];
      }
    }
  }

  // order初期化
  for (int i = 0; i < n; i++) {
    order[i] = i;
  }

  {
    double now_time = timer.get_elapsed_time();
    const double START_TEMP = annealingParams2.start_temperature[0];
    const double END_TEMP = annealingParams2.end_temperature;

    CoolingParams params;
    params.start_temp = START_TEMP;
    params.end_temp = END_TEMP;
    params.poly_p = 2.0;
    params.exp_k = 5.0;
    params.phase_split = 0.7;

    best_score = Sim(1);
    restore_best_score();
    int best_order[n];
    for (int i = 0; i < n; i++) {
      best_order[i] = order[i];
    }

    int keep_order[n];

    int loop2 = 0;

    now_time = timer.get_elapsed_time();
    double progress_ratio = now_time / timeLimit2;
    double temp = get_temp(progress_ratio, CoolingType::Poly, params);

    while (true) {
      //reset_state();
      loop2++;

      if (loop2 % 100 == 0) {
        now_time = timer.get_elapsed_time();
        progress_ratio = now_time / timeLimit2;
        temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;
        if (now_time > timeLimit2) {
          break;
        }
      }

      // 近傍解作成
      int ra_exec_mode = rand_xorshift() % annealingParams2.operation_thresholds[0];
      if (progress_ratio > 0.9) {
        ra_exec_mode = rand_xorshift() % annealingParams2.operation_thresholds[2];
      }
      int ra1, ra2, ra3, ra4, ra5;
      int keep1, keep2, keep3, keep4, keep5;

      for (int i = 0; i < n; i++) {
        keep_order[i] = order[i];
      }

      if (ra_exec_mode < annealingParams2.operation_thresholds[0]) {
        // 近傍操作1
        ra1 = rand_xorshift() % n;
        ra2 = rand_xorshift() % n;
        while (ra1 == ra2) {
          ra2 = rand_xorshift() % n;
        }
        swap(order[ra1], order[ra2]);
      }
      else if (ra_exec_mode < annealingParams2.operation_thresholds[1]) {
        // 近傍操作2 
        int seg_len = rand_xorshift() % (n / 10) + 1;
        ra1 = rand_xorshift() % (n - seg_len);
        ra2 = rand_xorshift() % (n - seg_len);
        while (ra1 == ra2) {
          ra2 = rand_xorshift() % (n - seg_len);
        }
        for (int i = 0; i < seg_len; i++) {
          swap(order[ra1 + i], order[ra2 + i]);
        }
      }
      else if (ra_exec_mode < annealingParams2.operation_thresholds[2]) {
        ra1 = rand_xorshift() % n;
        ra2 = rand_xorshift() % initial_c[ra1];
        keep1 = best_target[ra1][ra2];
        best_target[ra1][ra2] = sort_a_idx[ra1][rand_xorshift() % (n - 70)];
        while (best_target[ra1][ra2] == ra1) {
          best_target[ra1][ra2] = sort_a_idx[ra1][rand_xorshift() % (n - 70)];
        }

        // sum_attack更新
        sum_attack[keep1] -= a[ra1][keep1];
        sum_attack[best_target[ra1][ra2]] += a[ra1][best_target[ra1][ra2]];
      }

      int tmp_score = Sim(0);

      // 焼きなましで採用判定
      double diff_score = (tmp_score - score) * annealingParams2.score_scale;
      double prob = exp(diff_score / temp);
      bool update = false;
      if (tmp_score >= score) {
        update = true;
      }
      else if (ra_exec_mode < annealingParams2.operation_thresholds[0] && rand_01() < prob) {
        update = true;
      }
      if (update) {
        // 採用
        score = tmp_score;

        // ベスト更新
        if (score > best_score) {
          best_score = score;
          for (int i = 0; i < n; i++) {
            best_order[i] = order[i];
          }
        }
      }
      else {
        // 元に戻す
        if (ra_exec_mode < annealingParams2.operation_thresholds[0]) {
          // 近傍操作1 の巻き戻し
          swap(order[ra1], order[ra2]);
        }
        else if (ra_exec_mode < annealingParams2.operation_thresholds[1]) {
          // 近傍操作2
          for (int i = 0; i < n; i++) {
            order[i] = keep_order[i];
          }
        }
        else if (ra_exec_mode < annealingParams2.operation_thresholds[2]) {
          // 近傍操作3 の巻き戻し

          // sum_attack更新
          sum_attack[keep1] += a[ra1][keep1];
          sum_attack[best_target[ra1][ra2]] -= a[ra1][best_target[ra1][ra2]];

          best_target[ra1][ra2] = keep1;
        }
      }
    }

    if (exec_mode != 0 && exec_mode != 3) {
      cerr << loop2 << endl;
    }

    restore_best_score();
    for (int i = 0; i < n; i++) {
      order[i] = best_order[i];
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < c[i]; j++) {
        target[i][j] = best_target[i][j];
      }
    }
    score = Sim(1);

    store_best_score();
  }

  cerr << "Final Score: " << best_score << endl;
}

ll solve_case(int case_num, AnnealingParams annealingParams, AnnealingParams annealingParams2)
{
  timer.start();

  initialize_state();

  std::ifstream fin;
  std::istream& is = open_input_stream(case_num, fin, exec_mode);
  read_case(is);

  std::ofstream fout;
  std::ostream& os = open_output_stream(case_num, fout, exec_mode);

  best_score = -INF;

  build_initial_solution();

  store_best_score();

  // 焼きなまし実行
  // run_simulated_annealing(annealingParams);
  const double TIME_LIMIT = 1.8;
  best_best_score = 0;
  int LOOP_COUNT = 1;
  for (int loop = 0; loop < LOOP_COUNT; loop++) {
    timer.start();
    double timeLimit2 = TIME_LIMIT / LOOP_COUNT;
    double timeLimit1 = (timeLimit2) * 0.5;
    Method1(annealingParams, annealingParams2, timeLimit1, timeLimit2);

    if (best_score > best_best_score) {
      store_best_best_score();
    }

    cerr << "Loop: " << loop << ", Best Score: " << best_score << ", Best Best Score: " << best_best_score << ", Time: " << timer.get_elapsed_time() << endl;
  }

  restore_best_best_score();
  restore_best_score();

  // 解答を出力
  write_solution(os);

  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main()
{
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 200048.0;
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

  AnnealingParams annealingParams2 = annealingParams;
  annealingParams2.start_temperature[0] = 50000.0;
  annealingParams2.operation_thresholds[0] = 100;
  annealingParams2.operation_thresholds[1] = 100;
  annealingParams2.operation_thresholds[2] = 150;

  if (exec_mode == 0) {
    solve_case(0, annealingParams, annealingParams2);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 3; ++i) {
      ll score = solve_case(i, annealingParams, annealingParams2);
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
