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

int exec_mode;

class Board
{
public:
  static constexpr int MIN_LINE = 5;
  static constexpr int MAX_LINE = 50;

  int n;
  int a[14];
  int a_sum;
  vector<vector<int>> x, y;

  Board() : n(0), x(20000), y(20000)
  {
  }

  void init()
  {
    for (int i = 0; i < 20000; i++) {
      sort(x[i].begin(), x[i].end());
      sort(y[i].begin(), y[i].end());
    }

    a_sum = 0;
    for (int i = 1; i <= 10; i++) {
      a_sum += a[i];
    }
  }
};

class Answer
{
public:
  static constexpr int MIN = 0;
  static constexpr int MAX = 20000;

  vector<int> xs, ys; // 半開区間 : [xs[i], xs[i+1])

  vector<vector<int>> counts;
  int b[14];

  void initialize(int v_num, int h_num)
  {
    xs.clear();
    xs.push_back(MIN);
    for (int i = 1; i < v_num + 1; i++) {
      xs.push_back(MIN + (MAX - MIN) * i / (v_num + 1));
    }
    xs.push_back(MAX);

    ys.clear();
    ys.push_back(MIN);
    for (int i = 1; i < h_num + 1; i++) {
      ys.push_back(MIN + (MAX - MIN) * i / (h_num + 1));
    }
    ys.push_back(MAX);
  }

  void calc_counts(const Board& board)
  {
    for (int i = 0; i < 14; i++) {
      b[i] = 0;
    }

    counts.resize(xs.size() - 1, vector<int>(ys.size() - 1));
    for (int i = 0; i < xs.size() - 1; i++) {
      int sx = xs[i];
      int gx = xs[i + 1];
      for (int j = 0; j < ys.size() - 1; j++) {
        int sy = ys[j];
        int gy = ys[j + 1];

        counts[i][j] = 0;
        for (int k = sx; k < gx; k++) {
          for (auto l : board.x[k]) {
            if (sy <= l && l < gy) {
              counts[i][j]++;
            }
          }
        }

        if (counts[i][j] <= 12) {
          b[counts[i][j]]++;
        }
        else {
          b[13]++;
        }
      }
    }
  }

  void update_xs(const Board& board, int num, int diff)
  {
    if (diff < 0) {
      for (int i = xs[num] + diff; i < xs[num]; ++i) {
        for (auto y : board.x[i]) {
          //for (int j = 0; j < ys.size() - 1; j++) {
          {
            int j = upper_bound(ys.begin(), ys.end(), y) - ys.begin() - 1;
            if (ys[j] <= y && y < ys[j + 1]) {
              // 上側を減らす
              if (counts[num - 1][j] <= 12) {
                b[counts[num - 1][j]]--;
              }
              else {
                b[13]--;
              }
              counts[num - 1][j]--;
              if (counts[num - 1][j] <= 12) {
                b[counts[num - 1][j]]++;
              }
              else {
                b[13]++;
              }

              // 下側を増やす
              if (counts[num][j] <= 12) {
                b[counts[num][j]]--;
              }
              else {
                b[13]--;
              }
              counts[num][j]++;
              if (counts[num][j] <= 12) {
                b[counts[num][j]]++;
              }
              else {
                b[13]++;
              }
            }
          }
        }
      }
    }
    else {
      for (int i = xs[num]; i <= xs[num] + diff - 1; ++i) {
        for (auto y : board.x[i]) {
          //for (int j = 0; j < ys.size() - 1; j++) {
          {
            int j = upper_bound(ys.begin(), ys.end(), y) - ys.begin() - 1;
            if (ys[j] <= y && y < ys[j + 1]) {
              // 上側を増やす
              if (counts[num - 1][j] <= 12) {
                b[counts[num - 1][j]]--;
              }
              else {
                b[13]++;
              }
              counts[num - 1][j]++;
              if (counts[num - 1][j] <= 12) {
                b[counts[num - 1][j]]++;
              }
              else {
                b[13]++;
              }

              // 下側を減らす
              if (counts[num][j] <= 12) {
                b[counts[num][j]]--;
              }
              else {
                b[13]++;
              }
              counts[num][j]--;
              if (counts[num][j] <= 12) {
                b[counts[num][j]]++;
              }
              else {
                b[13]++;
              }
            }
          }
        }
      }
    }

    xs[num] += diff;
  }

  void update_ys(const Board& board, int num, int diff)
  {
    if (diff < 0) {
      for (int i = ys[num] + diff; i < ys[num]; ++i) {
        for (auto x : board.y[i]) {
          //for (int j = 0; j < xs.size() - 1; j++) {
          {
            int j = upper_bound(xs.begin(), xs.end(), x) - xs.begin() - 1;
            if (xs[j] <= x && x < xs[j + 1]) {
              // 左側を減らす
              if (counts[j][num - 1] <= 12) {
                b[counts[j][num - 1]]--;
              }
              else {
                b[13]++;
              }
              counts[j][num - 1]--;
              if (counts[j][num - 1] <= 12) {
                b[counts[j][num - 1]]++;
              }
              else {
                b[13]++;
              }
              // 右側を増やす
              if (counts[j][num] <= 12) {
                b[counts[j][num]]--;
              }
              else {
                b[13]++;
              }
              counts[j][num]++;
              if (counts[j][num] <= 12) {
                b[counts[j][num]]++;
              }
              else {
                b[13]++;
              }
            }
          }
        }
      }
    }
    else {
      for (int i = ys[num]; i <= ys[num] + diff - 1; ++i) {
        for (auto x : board.y[i]) {
          //for (int j = 0; j < xs.size() - 1; j++) {
          {
            int j = upper_bound(xs.begin(), xs.end(), x) - xs.begin() - 1;
            if (xs[j] <= x && x < xs[j + 1]) {
              // 左側を増やす
              if (counts[j][num - 1] <= 12) {
                b[counts[j][num - 1]]--;
              }
              else {
                b[13]++;
              }
              counts[j][num - 1]++;
              if (counts[j][num - 1] <= 12) {
                b[counts[j][num - 1]]++;
              }
              else {
                b[13]++;
              }
              // 右側を減らす
              if (counts[j][num] <= 12) {
                b[counts[j][num]]--;
              }
              else {
                b[13]++;
              }
              counts[j][num]--;
              if (counts[j][num] <= 12) {
                b[counts[j][num]]++;
              }
              else {
                b[13]++;
              }
            }
          }
        }
      }
    }

    ys[num] += diff;
  }
};

Board input_data(int case_num)
{
  Board board;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int  _k;
    cin >> board.n >> _k;
    for (int i = 0; i < 14; i++) {
      board.a[i] = 0;
    }
    for (int i = 1; i <= 10; i++) {
      cin >> board.a[i];
    }
    for (int i = 0; i < board.n; i++) {
      int x, y;
      cin >> x >> y;
      x += 10000;
      y += 10000;
      board.x[x].push_back(y);
      board.y[y].push_back(x);
    }
  }
  else {
    // ファイル入力
    int _k;
    ifs >> board.n >> _k;
    for (int i = 1; i <= 10; i++) {
      ifs >> board.a[i];
    }
    for (int i = 0; i < board.n; i++) {
      int x, y;
      ifs >> x >> y;
      x += 10000;
      y += 10000;
      board.x[x].push_back(y);
      board.y[y].push_back(x);
    }
  }

  board.init();

  return board;
}

int calculate_score(const Board& board, const Answer& answer)
{
  int ok_cnt = 0;
  for (int i = 1; i <= 10; i++) {
    ok_cnt += min(answer.b[i], board.a[i]);
  }

  double penalty = 0;
  for (int i = 1; i <= 10; i++) {
    penalty += max(0, answer.b[i] - board.a[i]) * (11 - i) * 0.01;
  }

  int res = round(1e6 * ok_cnt / board.a_sum - penalty);

  return res;
}

void output_data(int case_num, const Answer& answer)
{
  if (exec_mode == 0) {
    // 標準出力
    cout << answer.xs.size() + answer.ys.size() - 4 << endl;
    for (int i = 1; i < answer.xs.size() - 1; i++) {
      cout << answer.xs[i] - 10000 - 1 << " " << -10000000 << " " << answer.xs[i] - 10000 << " " << 10000000 << endl;
    }
    for (int i = 1; i < answer.ys.size() - 1; i++) {
      cout << -10000000 << " " << answer.ys[i] - 10000 - 1 << " " << 10000000 << " " << answer.ys[i] - 10000 << endl;
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << answer.xs.size() + answer.ys.size() - 4 << endl;
    for (int i = 1; i < answer.xs.size() - 1; i++) {
      ofs << answer.xs[i] - 10000 - 1 << " " << -10000000 << " " << answer.xs[i] - 10000 << " " << 10000000 << endl;
    }
    for (int i = 1; i < answer.ys.size() - 1; i++) {
      ofs << -10000000 << " " << answer.ys[i] - 10000 - 1 << " " << 10000000 << " " << answer.ys[i] - 10000 << endl;
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

Answer build_one_solution(const Board& board)
{
  Answer answer;
  int sum = rand_xorshift() % 61 + 40;
  int v_num = rand_xorshift() % (sum - 1) + 1;
  int h_num = sum - v_num;
  answer.initialize(v_num, h_num);
  answer.calc_counts(board);
  return answer;
}

Answer build_one_solution_2(const Board& board)
{
  double best_pieces = board.a_sum * 4 / 3.141592;

  Answer answer;
  int sum = rand_xorshift() % 61 + 40;
  int v_num = rand_xorshift() % (sum - 1) + 1;
  int h_num = sum - v_num;
  while (true) {
    double pieces = (v_num + 1) * (h_num + 1);
    double ratio = pieces / best_pieces;
    if (0.90 < ratio && ratio < 1.1) {
      break;
    }
    sum = rand_xorshift() % 61 + 40;
    v_num = rand_xorshift() % (sum - 1) + 1;
    h_num = sum - v_num;
  }
  answer.initialize(v_num, h_num);
  answer.calc_counts(board);
  return answer;
}

struct AnnealingParams
{
  double start_temperature;
  double end_temperature;
  double score_scale;
};

void annealing(const Board& board, Answer& answer, const AnnealingParams& params, double time_limit, bool print)
{
  Answer best_answer = answer;
  answer.calc_counts(board);
  int best_score = calculate_score(board, answer);

  double start_time = get_elapsed_time();
  double now_time = start_time;
  int iter = 0;
  while (true) {
    iter++;

    if (iter % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) {
        break;
      }
    }

    int choice = rand_xorshift() % 2;
    int num = 1;
    if (choice == 0) {
      num = rand_xorshift() % (answer.xs.size() - 2) + 1;
    }
    else {
      num = rand_xorshift() % (answer.ys.size() - 2) + 1;
    }
    int diff = rand_xorshift() % 61 - 30;
    while (diff == 0) {
      diff = rand_xorshift() % 61 - 30;
    }

    int current_score = calculate_score(board, answer);
    if (choice == 0) {
      if (answer.xs[num] + diff <= answer.xs[num - 1]) {
        continue;
      }
      if (answer.xs[num] + diff >= answer.xs[num + 1]) {
        continue;
      }
      answer.update_xs(board, num, diff);
    }
    else {
      if (answer.ys[num] + diff <= answer.ys[num - 1]) {
        continue;
      }
      if (answer.ys[num] + diff >= answer.ys[num + 1]) {
        continue;
      }
      answer.update_ys(board, num, diff);
    }

    double progress_ratio = (now_time - start_time) / (time_limit - start_time);
    double temp = params.start_temperature + (params.end_temperature - params.start_temperature) * progress_ratio;

    int new_score = calculate_score(board, answer);
    double diff_score = (new_score - current_score) * params.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用

      if (new_score > best_score) {
        best_score = new_score;
        best_answer = answer;
        if (best_score == 1000000) {
          break;
        }
      }
    }
    else {
      if (choice == 0) {
        answer.update_xs(board, num, -diff);
      }
      else {
        answer.update_ys(board, num, -diff);
      }
    }
  }

  if (print) {
    cerr << "iter = " << iter << ", score = " << best_score << endl;
  }

  answer = best_answer;
}

Answer build_initial_solution(const Board& board, double time_limit)
{
  Answer best_answer;
  int best_score = -1;

  int iter = 0;
  while (get_elapsed_time() < time_limit) {
    iter++;
    //Answer answer = build_one_solution(board);
    Answer answer = build_one_solution_2(board);
    int score = calculate_score(board, answer);
    if (score > best_score) {
      best_score = score;
      best_answer = answer;
    }
  }
  cerr << "iter = " << iter << ", score = " << best_score << endl;

  return best_answer;
}

Answer build_initial_solution_2(const Board& board, double time_limit)
{
  Answer best_answer;
  int best_score = -1;

  int iter = 0;
  while (get_elapsed_time() < time_limit) {
    iter++;
    Answer answer = build_one_solution_2(board);

    AnnealingParams params;
    params.start_temperature = 10048.0;
    params.end_temperature = 0.0;
    params.score_scale = 12345.0;
    annealing(board, answer, params, get_elapsed_time() + 0.002, false);

    int score = calculate_score(board, answer);
    if (score > best_score) {
      best_score = score;
      best_answer = answer;
    }
  }
  cerr << "iter = " << iter << ", score = " << best_score << endl;

  return best_answer;
}

const int VEC_SIZE = 2;
vector<Answer> build_initial_solution_3(const Board& board, double time_limit)
{
  vector<Answer> best_answers(VEC_SIZE);
  int best_scores[VEC_SIZE];
  for (int i = 0; i < VEC_SIZE; i++) {
    best_scores[i] = -1;
  }

  int iter = 0;
  while (get_elapsed_time() < time_limit) {
    iter++;
    Answer answer = build_one_solution_2(board);

    AnnealingParams params;
    params.start_temperature = 10048.0;
    params.end_temperature = 0.0;
    params.score_scale = 12345.0;
    annealing(board, answer, params, get_elapsed_time() + 0.002, false);

    int score = calculate_score(board, answer);
    if (score > best_scores[VEC_SIZE - 1]) {
      best_scores[VEC_SIZE - 1] = score;
      best_answers[VEC_SIZE - 1] = answer;

      int idx = VEC_SIZE - 1;
      while (idx > 0 && best_scores[idx] > best_scores[idx - 1]) {
        swap(best_scores[idx], best_scores[idx - 1]);
        swap(best_answers[idx], best_answers[idx - 1]);
        idx--;
      }
    }
  }
  cerr << "iter = " << iter << ", score = " << best_scores[0] << endl;

  return best_answers;
}

ll solve_case(int case_num)
{
  start_timer();

  Board board = input_data(case_num);

  const double TIME_LIMIT = 2.9;

  //Answer answer = build_initial_solution_2(board, TIME_LIMIT / 3);
  vector<Answer> answers = build_initial_solution_3(board, TIME_LIMIT / 3);

  Answer best_answer = answers[0];
  int best_score = calculate_score(board, best_answer);

  AnnealingParams params;
  params.start_temperature = 10000048.0;
  params.end_temperature = 0.0;
  params.score_scale = 12345.0;

  double now_time = get_elapsed_time();
  for (int i = 0; i < VEC_SIZE; i++) {
    double time_limit = now_time + (TIME_LIMIT - now_time) / VEC_SIZE * (i + 1);
    annealing(board, answers[i], params, time_limit, true);
    int score = calculate_score(board, answers[i]);
    if (score > best_score) {
      best_score = score;
      best_answer = answers[i];
    }

    cerr << i << ' ' << best_score << ' ' << score << endl;
    if (best_score == 1000000) {
      break;
    }
  }

  output_data(case_num, best_answer);

  int score = 0;
  if (exec_mode != 0) {
    score = calculate_score(board, best_answer);
    cerr << best_answer.xs.size() + best_answer.ys.size() - 4 << endl;
  }
  return score;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    int sum_score = 0;
    for (int i = 0; i < 150; i++) {
      int score = solve_case(i);
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

  return 0;
}
