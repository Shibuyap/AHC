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

const double TIME_LIMIT = 1.8;
int exec_mode;

int n;

int current_score;

int best_score;

void store_best_score()
{
  best_score = current_score;
}

void restore_best_score()
{
  current_score = best_score;
}

bool is_out_of_range(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

void initialize_state()
{
  current_score = 0;
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
}

void write_solution(std::ostream& os)
{
  os << 0 << endl;
}

ll calculate_score()
{
  ll res = 0;
  return res;
}

void build_initial_solution()
{
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams)
{
  store_best_score();

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
      if (now_time > TIME_LIMIT) { break; }
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = get_temp(progress_ratio, CoolingType::Linear, params);

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[1];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算
    double tmp_score = calculate_score();

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score > best_score) {
        store_best_score();
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  timer.start();

  initialize_state();

  std::ifstream fin;
  std::istream& is = open_input_stream(case_num, fin, exec_mode);
  read_case(is);

  std::ofstream fout;
  std::ostream& os = open_output_stream(case_num, fout, exec_mode);

  build_initial_solution();

  // 焼きなまし実行
  // run_simulated_annealing(annealingParams);

  write_solution(os);

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main()
{
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 2048.0;
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
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = solve_case(i, annealingParams);
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
      for (int i = 0; i < 15; ++i) {
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
