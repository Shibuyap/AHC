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

//------------------------------------------------------------------------------
// 定数 (UPPER_SNAKE_CASE)
//------------------------------------------------------------------------------
const ll INF_LL = 1001001001001001001LL;  // 旧: INF
const int DX[6] = { -1, 0, 0, 1, 0, 0 };
const int DY[6] = { 0, -1, 0, 0, 1, 0 };
const int DZ[6] = { 0, 0, -1, 0, 0, 1 };

//------------------------------------------------------------------------------
// 関数プロトタイプ (先に宣言しておく)
//------------------------------------------------------------------------------
int get_direction(int num);
void clear_all_multicase();
void init_state();
void normal_clear();
void real_clear();
void seed_clear();
void copy_to_real();
void copy_to_seed();
void roll_back_from_real();
void roll_back_from_seed();
double calc_score();
bool can_delete_block(int i, int x, int y, int z);
bool is_invalid_coord(int x, int y, int z);
void method_1(double temperature);
void method_2(double temperature);
void method_3(double temperature);
void method_4(double temperature);
double solve_problem(int run_mode, int problem_num = 0);
double solve_outer(int run_mode, int problem_num = 0);
void read_input(int problem_num);
void write_output(int run_mode, int problem_num);

//------------------------------------------------------------------------------
// グローバル変数
//------------------------------------------------------------------------------
namespace
{
  //=== 入力関連 (旧: 変数) ===//
  int dimension;             // 旧: D
  bool f_matrix[2][20][20];  // 旧: F
  bool r_matrix[2][20][20];  // 旧: R

  //=== 解答用変数 ===//
  double min_score;               // 旧: minScore
  int answer_grid[2][15][15][15]; // 旧: ans
  int block_count[2][100];        // 旧: bcount

  double real_min_score;               // 旧: real_minScore
  int real_answer_grid[2][15][15][15]; // 旧: real_ans
  int real_block_count[2][100];        // 旧: real_bcount

  double seed_min_score;               // 旧: seed_minScore
  int seed_answer_grid[2][15][15][15]; // 旧: seed_ans
  int seed_block_count[2][100];        // 旧: seed_bcount

  int method_count[20][2];  // 旧: methodCount

} // unnamed namespace

//------------------------------------------------------------------------------
// (1) ユーティリティ系 関数
//------------------------------------------------------------------------------

// 旧: GetDir
int get_direction(int num)
{
  for (int i = 0; i < 6; ++i) {
    if (num & (1 << i)) return i;
  }
  return -1;
}

//------------------------------------------------------------------------------
// (2) 乱数ライブラリ
//------------------------------------------------------------------------------
static uint32_t rand32()
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
  return (rand32() + 0.5) * (1.0 / UINT_MAX);
}

//------------------------------------------------------------------------------
// (3) 変数の初期化リセット系 関数
//------------------------------------------------------------------------------
void normal_clear()
{
  min_score = INF_LL;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    answer_grid[i][j][k][l] = -1;
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    block_count[i][j] = 0;
  }
}

void real_clear()
{
  real_min_score = INF_LL;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    real_answer_grid[i][j][k][l] = -1;
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    real_block_count[i][j] = 0;
  }
}

void seed_clear()
{
  seed_min_score = INF_LL;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    seed_answer_grid[i][j][k][l] = -1;
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    seed_block_count[i][j] = 0;
  }
}

// ローカルで複数ケース試すための全て消す関数 (旧: AllClear_MultiCase)
void clear_all_multicase()
{
  normal_clear();
  real_clear();
  seed_clear();
}

//------------------------------------------------------------------------------
// (4) 初期状態作成 (旧: Init)
//------------------------------------------------------------------------------
void init_state()
{
  // block_countの初期化
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 100; ++j) {
      block_count[i][j] = 0;
    }
  }
  // answer_grid の初期化 (FやRがtrueなら最初からブロックID=0)
  for (int i = 0; i < 2; ++i) {
    for (int x = 0; x < dimension; ++x) {
      for (int y = 0; y < dimension; ++y) {
        for (int z = 0; z < dimension; ++z) {
          if (f_matrix[i][z][x] && r_matrix[i][z][y]) {
            answer_grid[i][x][y][z] = 0;
            block_count[i][1]++;
          }
          else {
            answer_grid[i][x][y][z] = -1;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// (5) 入力受け取り (旧: Input)
//------------------------------------------------------------------------------
void read_input(int problem_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problem_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> dimension;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < dimension; ++j) {
        string s;
        cin >> s;
        for (int k = 0; k < dimension; ++k) {
          f_matrix[i][j][k] = (s[k] - '0');
        }
      }
      for (int j = 0; j < dimension; ++j) {
        string s;
        cin >> s;
        for (int k = 0; k < dimension; ++k) {
          r_matrix[i][j][k] = (s[k] - '0');
        }
      }
    }
  }
  else {
    // ファイル入力
    ifs >> dimension;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < dimension; ++j) {
        string s;
        ifs >> s;
        for (int k = 0; k < dimension; ++k) {
          f_matrix[i][j][k] = (s[k] - '0');
        }
      }
      for (int j = 0; j < dimension; ++j) {
        string s;
        ifs >> s;
        for (int k = 0; k < dimension; ++k) {
          r_matrix[i][j][k] = (s[k] - '0');
        }
      }
    }
  }

  // 初期化
  init_state();
}

//------------------------------------------------------------------------------
// (6) 解答出力 (旧: Output)
//------------------------------------------------------------------------------
void write_output(int run_mode, int problem_num)
{
  int ans_n[100] = {};
  int ans_sum[100] = {};
  // 大きい方を取って合計するロジック
  for (int j = 0; j < 100; ++j) {
    ans_n[j] = max(block_count[0][j], block_count[1][j]);
    ans_sum[j] = ans_n[j];
    if (j > 0) ans_sum[j] += ans_sum[j - 1];
  }

  int ans_print[2][15][15][15];
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < dimension; ++j) {
      for (int k = 0; k < dimension; ++k) {
        for (int l = 0; l < dimension; ++l) {
          ans_print[i][j][k][l] = 0;
        }
      }
    }
  }

  // 各ブロックのID割り振り
  for (int i = 0; i < 2; ++i) {
    int count_sum[100] = {};
    for (int j = 1; j < 100; ++j) {
      count_sum[j] = ans_sum[j - 1];
    }
    for (int j = 0; j < dimension; ++j) {
      for (int k = 0; k < dimension; ++k) {
        for (int l = 0; l < dimension; ++l) {
          if (ans_print[i][j][k][l] != 0) { continue; }
          if (answer_grid[i][j][k][l] == -1) {
            ans_print[i][j][k][l] = 0;
          }
          else if (answer_grid[i][j][k][l] == 0) {
            count_sum[1]++;
            ans_print[i][j][k][l] = count_sum[1];
          }
          else {
            // 2ブロックつなぎ
            count_sum[2]++;
            ans_print[i][j][k][l] = count_sum[2];
            int dir = get_direction(answer_grid[i][j][k][l]);
            int nx = j + DX[dir];
            int ny = k + DY[dir];
            int nz = l + DZ[dir];
            ans_print[i][nx][ny][nz] = count_sum[2];
          }
        }
      }
    }
  }

  if (run_mode == 0) {
    // 標準出力
    cout << ans_sum[99] << endl;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < dimension; ++j) {
        for (int k = 0; k < dimension; ++k) {
          for (int l = 0; l < dimension; ++l) {
            cout << ans_print[i][j][k][l] << ' ';
          }
        }
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    string file_name_ofs = "./out/";
    {
      string str_num;
      for (int i = 0; i < 4; ++i) {
        str_num += char((problem_num % 10) + '0');
        problem_num /= 10;
      }
      reverse(str_num.begin(), str_num.end());
      file_name_ofs += str_num + ".txt";
    }

    ofstream ofs(file_name_ofs);

    ofs << ans_sum[99] << endl;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < dimension; ++j) {
        for (int k = 0; k < dimension; ++k) {
          for (int l = 0; l < dimension; ++l) {
            ofs << ans_print[i][j][k][l] << ' ';
          }
        }
      }
      ofs << endl;
    }
    ofs.close();
  }
}

//------------------------------------------------------------------------------
// (7) ロールバック用 関数
//------------------------------------------------------------------------------
void copy_to_real()
{
  real_min_score = min_score;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    real_answer_grid[i][j][k][l] = answer_grid[i][j][k][l];
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    real_block_count[i][j] = block_count[i][j];
  }
}

void copy_to_seed()
{
  seed_min_score = min_score;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    seed_answer_grid[i][j][k][l] = answer_grid[i][j][k][l];
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    seed_block_count[i][j] = block_count[i][j];
  }
}

void roll_back_from_real()
{
  min_score = real_min_score;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    answer_grid[i][j][k][l] = real_answer_grid[i][j][k][l];
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    block_count[i][j] = real_block_count[i][j];
  }
}

void roll_back_from_seed()
{
  min_score = seed_min_score;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < dimension; ++j) for (int k = 0; k < dimension; ++k) for (int l = 0; l < dimension; ++l) {
    answer_grid[i][j][k][l] = seed_answer_grid[i][j][k][l];
  }
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 100; ++j) {
    block_count[i][j] = seed_block_count[i][j];
  }
}

//------------------------------------------------------------------------------
// (8) スコア計算 (旧: CalcScore)
//------------------------------------------------------------------------------
double calc_score()
{
  double resd = 0;
  double ma, mi;
  // srep(j, 1, 3)  => j=1,2
  for (int j = 1; j < 3; ++j) {
    if (block_count[0][j] >= block_count[1][j]) {
      ma = block_count[0][j];
      mi = block_count[1][j];
    }
    else {
      ma = block_count[1][j];
      mi = block_count[0][j];
    }
    // (ma - mi)*j + mi/j
    resd += (ma - mi) * j;
    resd += mi / j;
  }
  return resd;
}

//------------------------------------------------------------------------------
// (9) ブロック消去の可否判定 (旧: CanDelete)
//------------------------------------------------------------------------------
bool can_delete_block(int i, int x, int y, int z)
{
  // 条件: F[i][z][x] が true の場合
  if (f_matrix[i][z][x]) {
    bool ok = false;
    for (int k = 0; k < dimension; ++k) {
      if (k == y) { continue; }
      if (answer_grid[i][x][k][z] != -1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  // 条件: R[i][z][y] が true の場合
  if (r_matrix[i][z][y]) {
    bool ok = false;
    for (int j = 0; j < dimension; ++j) {
      if (j == x) { continue; }
      if (answer_grid[i][j][y][z] != -1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// (10) 座標チェック (旧: IsNG)
//------------------------------------------------------------------------------
bool is_invalid_coord(int x, int y, int z)
{
  if (x < 0 || x >= dimension) return true;
  if (y < 0 || y >= dimension) return true;
  if (z < 0 || z >= dimension) return true;
  return false;
}

//------------------------------------------------------------------------------
// (11) メソッド群 (焼きなまし操作) - 旧: Method1/2/3/4
//------------------------------------------------------------------------------
void method_1(double temperature)
{
  int i = rand32() % 2;
  int x = rand32() % dimension;
  int y = rand32() % dimension;
  int z = rand32() % dimension;
  if (answer_grid[i][x][y][z] != -1) { return; }
  if (!f_matrix[i][z][x] || !r_matrix[i][z][y]) { return; }

  method_count[1][1]++;
  answer_grid[i][x][y][z] = 0;
  block_count[i][1]++;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_01()) {
    // 受け入れる
    method_count[1][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // 元に戻す
    block_count[i][1]--;
    answer_grid[i][x][y][z] = -1;
  }
}

void method_2(double temperature)
{
  int i = rand32() % 2;
  int x = rand32() % dimension;
  int y = rand32() % dimension;
  int z = rand32() % dimension;
  if (answer_grid[i][x][y][z] != 0) { return; }
  if (!can_delete_block(i, x, y, z)) { return; }

  method_count[2][1]++;
  answer_grid[i][x][y][z] = -1;
  block_count[i][1]--;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_01()) {
    method_count[2][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // 元に戻す
    block_count[i][1]++;
    answer_grid[i][x][y][z] = 0;
  }
}

void method_3(double temperature)
{
  int i = rand32() % 2;
  int x = rand32() % dimension;
  int y = rand32() % dimension;
  int z = rand32() % dimension;
  int dir = rand32() % 6;
  if (answer_grid[i][x][y][z] != 0) { return; }
  int nx = x + DX[dir];
  int ny = y + DY[dir];
  int nz = z + DZ[dir];
  if (is_invalid_coord(nx, ny, nz)) { return; }
  if (answer_grid[i][nx][ny][nz] != 0) { return; }

  method_count[3][1]++;
  // つなげる
  answer_grid[i][x][y][z] = (1 << dir);
  answer_grid[i][nx][ny][nz] = (1 << ((dir + 3) % 6));
  block_count[i][1] -= 2;
  block_count[i][2]++;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_01()) {
    method_count[3][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ロールバック
    block_count[i][2]--;
    block_count[i][1] += 2;
    answer_grid[i][x][y][z] = 0;
    answer_grid[i][nx][ny][nz] = 0;
  }
}

void method_4(double temperature)
{
  int i = rand32() % 2;
  int x = rand32() % dimension;
  int y = rand32() % dimension;
  int z = rand32() % dimension;
  if (answer_grid[i][x][y][z] <= 0) { return; }
  int dir = get_direction(answer_grid[i][x][y][z]);
  int nx = x + DX[dir];
  int ny = y + DY[dir];
  int nz = z + DZ[dir];

  method_count[4][1]++;
  // 分離
  answer_grid[i][x][y][z] = 0;
  answer_grid[i][nx][ny][nz] = 0;
  block_count[i][1] += 2;
  block_count[i][2]--;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_01()) {
    method_count[4][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ロールバック
    block_count[i][2]++;
    block_count[i][1] -= 2;
    answer_grid[i][x][y][z] = (1 << dir);
    answer_grid[i][nx][ny][nz] = (1 << ((dir + 3) % 6));
  }
}

//------------------------------------------------------------------------------
// (12) メインの焼きなまし処理 (旧: Solve)
//------------------------------------------------------------------------------
double solve_problem(int run_mode, int problem_num)
{
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();

  // 初期状態 (init_state は read_input内で呼ばれた後にも呼ばれるが念のため)
  init_state();

  // 最初のスコア
  min_score = calc_score();
  copy_to_real();
  copy_to_seed();

  //--- シード作成用 (初段階) ---
  int seed_count = 100;
  for (int tei = 0; tei < seed_count; ++tei) {
    start_time = clock();
    // 初期化し直す
    init_state();
    min_score = calc_score();

    // 焼きなまし
    end_time = clock();
    double now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
    double time_limit = 4.0 / seed_count; // 旧: TL
    double now_progress = now_time / time_limit;
    double start_temperature = 20;
    double end_temperature = 0;
    int loop_count = 0;
    int rollback_count = 0;

    while (true) {
      loop_count++;
      if (loop_count % 100 == 1) {
        end_time = clock();
        now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
        now_progress = now_time / time_limit;
      }
      if (now_progress > 1.0) { break; }

      // もしスコアが大きく悪化していればロールバック (コメントアウト原文通り)
      if (min_score > real_min_score * 10) {
        // roll_back_from_real();
        // rollback_count++;
      }

      double temperature =
        start_temperature + (end_temperature - start_temperature) * now_progress;

      int ra = rand32() % 100;
      if (ra < 20) {
        method_1(temperature);
      }
      else if (ra < 40) {
        method_2(temperature);
      }
      else if (ra < 70) {
        method_3(temperature);
      }
      else {
        method_4(temperature);
      }
    } // while(シード作成)

    // 最良を戻す
    roll_back_from_real();
    // シード更新
    if (min_score <= seed_min_score) {
      copy_to_seed();
    }
  }

  // シードからロールバック
  roll_back_from_seed();
  copy_to_real();

  //--- メイン焼きなまし ---
  start_time = clock();
  end_time = clock();
  double now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
  double time_limit = 0.9;  // 旧: TL
  double now_progress = now_time / time_limit;
  double start_temperature = 2;
  double end_temperature = 0;
  int loop_count = 0;
  int rollback_count = 0;

  while (true) {
    loop_count++;
    if (loop_count % 100 == 1) {
      end_time = clock();
      now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
      now_progress = now_time / time_limit;
    }
    if (now_progress > 1.0) { break; }

    double temperature =
      start_temperature + (end_temperature - start_temperature) * now_progress;

    int ra = rand32() % 100;
    if (ra < 20) {
      method_1(temperature);
    }
    else if (ra < 40) {
      method_2(temperature);
    }
    else if (ra < 70) {
      method_3(temperature);
    }
    else {
      method_4(temperature);
    }
  }

  // 最良解を適用
  roll_back_from_real();

  // デバッグ出力
  if (run_mode != 0) {
    cout << "problem_num = " << problem_num << ", D = " << dimension << endl;
    cout << "min_score = " << min_score << endl;
    cout << "loop_count = " << loop_count
      << ", rollback_count = " << rollback_count << endl;
    for (int i = 1; i < 5; ++i) {
      cout << "method_" << i << " = "
        << method_count[i][0] << " / " << method_count[i][1] << endl;
    }
    cout << endl;
  }

  return min_score;
}

//------------------------------------------------------------------------------
// (13) solve_outer (旧: SolveOuter)
//------------------------------------------------------------------------------
double solve_outer(int run_mode, int problem_num)
{
  // 入力受け取り
  read_input(problem_num);

  // 実行
  double score = solve_problem(run_mode, problem_num);

  // 解答の出力
  write_output(run_mode, problem_num);

  return score;
}

//------------------------------------------------------------------------------
// (14) main
//------------------------------------------------------------------------------
int main()
{
  int mode = 2;  // 0: 提出用, 1: 1ケースのみ, 2: 複数ケース

  if (mode == 0) {
    solve_outer(mode);
  }
  else if (mode == 1) {
    solve_outer(mode, 19);
  }
  else if (mode == 2) {
    double score_sum = 0;
    for (int i = 0; i < 100; ++i) {
      score_sum += solve_outer(mode, i);
      clear_all_multicase();
    }
    cout << "scoreSum = " << score_sum << endl;
  }

  return 0;
}
