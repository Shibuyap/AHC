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

//=== Macros (�K�v�ɉ����Ē���) ===//
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;

//------------------------------------------------------------------------------
// �萔 (UPPER_SNAKE_CASE)
//------------------------------------------------------------------------------
const ll INF_LL = 1001001001001001001LL;  // ��: INF
const int DX[6] = { -1, 0, 0, 1, 0, 0 };
const int DY[6] = { 0, -1, 0, 0, 1, 0 };
const int DZ[6] = { 0, 0, -1, 0, 0, 1 };

//------------------------------------------------------------------------------
// �֐��v���g�^�C�v (��ɐ錾���Ă���)
//------------------------------------------------------------------------------
int get_direction(int num);
double rand_double_01();
uint32_t rand_uint32();
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
// �O���[�o���ϐ�
//------------------------------------------------------------------------------
namespace {
  //=== ���͊֘A (��: �ϐ�) ===//
  int dimension;             // ��: D
  bool f_matrix[2][20][20];  // ��: F
  bool r_matrix[2][20][20];  // ��: R

  //=== �𓚗p�ϐ� ===//
  double min_score;               // ��: minScore
  int answer_grid[2][15][15][15]; // ��: ans
  int block_count[2][100];        // ��: bcount

  double real_min_score;               // ��: real_minScore
  int real_answer_grid[2][15][15][15]; // ��: real_ans
  int real_block_count[2][100];        // ��: real_bcount

  double seed_min_score;               // ��: seed_minScore
  int seed_answer_grid[2][15][15][15]; // ��: seed_ans
  int seed_block_count[2][100];        // ��: seed_bcount

  int method_count[20][2];  // ��: methodCount

} // unnamed namespace

//------------------------------------------------------------------------------
// (1) ���[�e�B���e�B�n �֐�
//------------------------------------------------------------------------------

// ��: GetDir
int get_direction(int num) {
  rep(i, 6) {
    if (num & (1 << i)) return i;
  }
  return -1;
}

//------------------------------------------------------------------------------
// (2) �������C�u����
//------------------------------------------------------------------------------
static uint32_t rand_uint32() {
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

static double rand_double_01() {
  return (rand_uint32() + 0.5) * (1.0 / UINT_MAX);
}

//------------------------------------------------------------------------------
// (3) �ϐ��̏��������Z�b�g�n �֐�
//------------------------------------------------------------------------------
void normal_clear() {
  min_score = INF_LL;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    answer_grid[i][j][k][l] = -1;
  }
  rep(i, 2) rep(j, 100) {
    block_count[i][j] = 0;
  }
}

void real_clear() {
  real_min_score = INF_LL;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    real_answer_grid[i][j][k][l] = -1;
  }
  rep(i, 2) rep(j, 100) {
    real_block_count[i][j] = 0;
  }
}

void seed_clear() {
  seed_min_score = INF_LL;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    seed_answer_grid[i][j][k][l] = -1;
  }
  rep(i, 2) rep(j, 100) {
    seed_block_count[i][j] = 0;
  }
}

// ���[�J���ŕ����P�[�X�������߂̑S�ď����֐� (��: AllClear_MultiCase)
void clear_all_multicase() {
  normal_clear();
  real_clear();
  seed_clear();
}

//------------------------------------------------------------------------------
// (4) ������ԍ쐬 (��: Init)
//------------------------------------------------------------------------------
void init_state() {
  // block_count�̏�����
  rep(i, 2) {
    rep(j, 100) {
      block_count[i][j] = 0;
    }
  }
  // answer_grid �̏����� (F��R��true�Ȃ�ŏ�����u���b�NID=0)
  rep(i, 2) {
    rep(x, dimension) {
      rep(y, dimension) {
        rep(z, dimension) {
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
// (5) ���͎󂯎�� (��: Input)
//------------------------------------------------------------------------------
void read_input(int problem_num) {
  string file_name_ifs = "./in/";
  {
    string str_num;
    rep(i, 4) {
      str_num += char((problem_num % 10) + '0');
      problem_num /= 10;
    }
    reverse(str_num.begin(), str_num.end());
    file_name_ifs += str_num + ".txt";
  }

  ifstream ifs(file_name_ifs);

  if (!ifs.is_open()) {
    // �W������
    cin >> dimension;
    rep(i, 2) {
      rep(j, dimension) {
        string s;
        cin >> s;
        rep(k, dimension) {
          f_matrix[i][j][k] = (s[k] - '0');
        }
      }
      rep(j, dimension) {
        string s;
        cin >> s;
        rep(k, dimension) {
          r_matrix[i][j][k] = (s[k] - '0');
        }
      }
    }
  }
  else {
    // �t�@�C������
    ifs >> dimension;
    rep(i, 2) {
      rep(j, dimension) {
        string s;
        ifs >> s;
        rep(k, dimension) {
          f_matrix[i][j][k] = (s[k] - '0');
        }
      }
      rep(j, dimension) {
        string s;
        ifs >> s;
        rep(k, dimension) {
          r_matrix[i][j][k] = (s[k] - '0');
        }
      }
    }
  }

  // ������
  init_state();
}

//------------------------------------------------------------------------------
// (6) �𓚏o�� (��: Output)
//------------------------------------------------------------------------------
void write_output(int run_mode, int problem_num) {
  int ans_n[100] = {};
  int ans_sum[100] = {};
  // �傫����������č��v���郍�W�b�N
  rep(j, 100) {
    ans_n[j] = max(block_count[0][j], block_count[1][j]);
    ans_sum[j] = ans_n[j];
    if (j > 0) ans_sum[j] += ans_sum[j - 1];
  }

  int ans_print[2][15][15][15];
  rep(i, 2) {
    rep(j, dimension) {
      rep(k, dimension) {
        rep(l, dimension) {
          ans_print[i][j][k][l] = 0;
        }
      }
    }
  }

  // �e�u���b�N��ID����U��
  rep(i, 2) {
    int count_sum[100] = {};
    srep(j, 1, 100) {
      count_sum[j] = ans_sum[j - 1];
    }
    rep(j, dimension) {
      rep(k, dimension) {
        rep(l, dimension) {
          if (ans_print[i][j][k][l] != 0) continue;
          if (answer_grid[i][j][k][l] == -1) {
            ans_print[i][j][k][l] = 0;
          }
          else if (answer_grid[i][j][k][l] == 0) {
            count_sum[1]++;
            ans_print[i][j][k][l] = count_sum[1];
          }
          else {
            // 2�u���b�N�Ȃ�
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
    // �W���o��
    cout << ans_sum[99] << endl;
    rep(i, 2) {
      rep(j, dimension) {
        rep(k, dimension) {
          rep(l, dimension) {
            cout << ans_print[i][j][k][l] << ' ';
          }
        }
      }
      cout << endl;
    }
  }
  else {
    // �t�@�C���o��
    string file_name_ofs = "./out/";
    {
      string str_num;
      rep(i, 4) {
        str_num += char((problem_num % 10) + '0');
        problem_num /= 10;
      }
      reverse(str_num.begin(), str_num.end());
      file_name_ofs += str_num + ".txt";
    }

    ofstream ofs(file_name_ofs);

    ofs << ans_sum[99] << endl;
    rep(i, 2) {
      rep(j, dimension) {
        rep(k, dimension) {
          rep(l, dimension) {
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
// (7) ���[���o�b�N�p �֐�
//------------------------------------------------------------------------------
void copy_to_real() {
  real_min_score = min_score;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    real_answer_grid[i][j][k][l] = answer_grid[i][j][k][l];
  }
  rep(i, 2) rep(j, 100) {
    real_block_count[i][j] = block_count[i][j];
  }
}

void copy_to_seed() {
  seed_min_score = min_score;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    seed_answer_grid[i][j][k][l] = answer_grid[i][j][k][l];
  }
  rep(i, 2) rep(j, 100) {
    seed_block_count[i][j] = block_count[i][j];
  }
}

void roll_back_from_real() {
  min_score = real_min_score;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    answer_grid[i][j][k][l] = real_answer_grid[i][j][k][l];
  }
  rep(i, 2) rep(j, 100) {
    block_count[i][j] = real_block_count[i][j];
  }
}

void roll_back_from_seed() {
  min_score = seed_min_score;
  rep(i, 2) rep(j, dimension) rep(k, dimension) rep(l, dimension) {
    answer_grid[i][j][k][l] = seed_answer_grid[i][j][k][l];
  }
  rep(i, 2) rep(j, 100) {
    block_count[i][j] = seed_block_count[i][j];
  }
}

//------------------------------------------------------------------------------
// (8) �X�R�A�v�Z (��: CalcScore)
//------------------------------------------------------------------------------
double calc_score() {
  double resd = 0;
  double ma, mi;
  // srep(j, 1, 3)  => j=1,2
  srep(j, 1, 3) {
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
// (9) �u���b�N�����̉۔��� (��: CanDelete)
//------------------------------------------------------------------------------
bool can_delete_block(int i, int x, int y, int z) {
  // ����: F[i][z][x] �� true �̏ꍇ
  if (f_matrix[i][z][x]) {
    bool ok = false;
    rep(k, dimension) {
      if (k == y) continue;
      if (answer_grid[i][x][k][z] != -1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  // ����: R[i][z][y] �� true �̏ꍇ
  if (r_matrix[i][z][y]) {
    bool ok = false;
    rep(j, dimension) {
      if (j == x) continue;
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
// (10) ���W�`�F�b�N (��: IsNG)
//------------------------------------------------------------------------------
bool is_invalid_coord(int x, int y, int z) {
  if (x < 0 || x >= dimension) return true;
  if (y < 0 || y >= dimension) return true;
  if (z < 0 || z >= dimension) return true;
  return false;
}

//------------------------------------------------------------------------------
// (11) ���\�b�h�Q (�Ă��Ȃ܂�����) - ��: Method1/2/3/4
//------------------------------------------------------------------------------
void method_1(double temperature) {
  int i = rand_uint32() % 2;
  int x = rand_uint32() % dimension;
  int y = rand_uint32() % dimension;
  int z = rand_uint32() % dimension;
  if (answer_grid[i][x][y][z] != -1) return;
  if (!f_matrix[i][z][x] || !r_matrix[i][z][y]) return;

  method_count[1][1]++;
  answer_grid[i][x][y][z] = 0;
  block_count[i][1]++;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_double_01()) {
    // �󂯓����
    method_count[1][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ���ɖ߂�
    block_count[i][1]--;
    answer_grid[i][x][y][z] = -1;
  }
}

void method_2(double temperature) {
  int i = rand_uint32() % 2;
  int x = rand_uint32() % dimension;
  int y = rand_uint32() % dimension;
  int z = rand_uint32() % dimension;
  if (answer_grid[i][x][y][z] != 0) return;
  if (!can_delete_block(i, x, y, z)) return;

  method_count[2][1]++;
  answer_grid[i][x][y][z] = -1;
  block_count[i][1]--;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_double_01()) {
    method_count[2][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ���ɖ߂�
    block_count[i][1]++;
    answer_grid[i][x][y][z] = 0;
  }
}

void method_3(double temperature) {
  int i = rand_uint32() % 2;
  int x = rand_uint32() % dimension;
  int y = rand_uint32() % dimension;
  int z = rand_uint32() % dimension;
  int dir = rand_uint32() % 6;
  if (answer_grid[i][x][y][z] != 0) return;
  int nx = x + DX[dir];
  int ny = y + DY[dir];
  int nz = z + DZ[dir];
  if (is_invalid_coord(nx, ny, nz)) return;
  if (answer_grid[i][nx][ny][nz] != 0) return;

  method_count[3][1]++;
  // �Ȃ���
  answer_grid[i][x][y][z] = (1 << dir);
  answer_grid[i][nx][ny][nz] = (1 << ((dir + 3) % 6));
  block_count[i][1] -= 2;
  block_count[i][2]++;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_double_01()) {
    method_count[3][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ���[���o�b�N
    block_count[i][2]--;
    block_count[i][1] += 2;
    answer_grid[i][x][y][z] = 0;
    answer_grid[i][nx][ny][nz] = 0;
  }
}

void method_4(double temperature) {
  int i = rand_uint32() % 2;
  int x = rand_uint32() % dimension;
  int y = rand_uint32() % dimension;
  int z = rand_uint32() % dimension;
  if (answer_grid[i][x][y][z] <= 0) return;
  int dir = get_direction(answer_grid[i][x][y][z]);
  int nx = x + DX[dir];
  int ny = y + DY[dir];
  int nz = z + DZ[dir];

  method_count[4][1]++;
  // ����
  answer_grid[i][x][y][z] = 0;
  answer_grid[i][nx][ny][nz] = 0;
  block_count[i][1] += 2;
  block_count[i][2]--;

  double tmp_score = calc_score();
  double diff_score = min_score - tmp_score;
  double prob = exp(diff_score / temperature);

  if (prob > rand_double_01()) {
    method_count[4][0]++;
    min_score = tmp_score;
    if (min_score < real_min_score) {
      copy_to_real();
    }
  }
  else {
    // ���[���o�b�N
    block_count[i][2]++;
    block_count[i][1] -= 2;
    answer_grid[i][x][y][z] = (1 << dir);
    answer_grid[i][nx][ny][nz] = (1 << ((dir + 3) % 6));
  }
}

//------------------------------------------------------------------------------
// (12) ���C���̏Ă��Ȃ܂����� (��: Solve)
//------------------------------------------------------------------------------
double solve_problem(int run_mode, int problem_num) {
  clock_t start_time, end_time;
  start_time = clock();
  end_time = clock();

  // ������� (init_state �� read_input���ŌĂ΂ꂽ��ɂ��Ă΂�邪�O�̂���)
  init_state();

  // �ŏ��̃X�R�A
  min_score = calc_score();
  copy_to_real();
  copy_to_seed();

  //--- �V�[�h�쐬�p (���i�K) ---
  int seed_count = 100;
  rep(tei, seed_count) {
    start_time = clock();
    // ������������
    init_state();
    min_score = calc_score();

    // �Ă��Ȃ܂�
    end_time = clock();
    double now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
    double time_limit = 4.0 / seed_count; // ��: TL
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
      if (now_progress > 1.0) break;

      // �����X�R�A���傫���������Ă���΃��[���o�b�N (�R�����g�A�E�g�����ʂ�)
      if (min_score > real_min_score * 10) {
        // roll_back_from_real();
        // rollback_count++;
      }

      double temperature =
        start_temperature + (end_temperature - start_temperature) * now_progress;

      int ra = rand_uint32() % 100;
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
    } // while(�V�[�h�쐬)

    // �ŗǂ�߂�
    roll_back_from_real();
    // �V�[�h�X�V
    if (min_score <= seed_min_score) {
      copy_to_seed();
    }
  }

  // �V�[�h���烍�[���o�b�N
  roll_back_from_seed();
  copy_to_real();

  //--- ���C���Ă��Ȃ܂� ---
  start_time = clock();
  end_time = clock();
  double now_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
  double time_limit = 0.9;  // ��: TL
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
    if (now_progress > 1.0) break;

    if (min_score > real_min_score * 10) {
      // roll_back_from_real();
      // rollback_count++;
    }

    double temperature =
      start_temperature + (end_temperature - start_temperature) * now_progress;

    int ra = rand_uint32() % 100;
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

  // �ŗǉ���K�p
  roll_back_from_real();

  // �f�o�b�O�o��
  if (run_mode != 0) {
    cout << "problem_num = " << problem_num << ", D = " << dimension << endl;
    cout << "min_score = " << min_score << endl;
    cout << "loop_count = " << loop_count
      << ", rollback_count = " << rollback_count << endl;
    srep(i, 1, 5) {
      cout << "method_" << i << " = "
        << method_count[i][0] << " / " << method_count[i][1] << endl;
    }
    cout << endl;
  }

  return min_score;
}

//------------------------------------------------------------------------------
// (13) solve_outer (��: SolveOuter)
//------------------------------------------------------------------------------
double solve_outer(int run_mode, int problem_num) {
  // ���͎󂯎��
  read_input(problem_num);

  // ���s
  double score = solve_problem(run_mode, problem_num);

  // �𓚂̏o��
  write_output(run_mode, problem_num);

  return score;
}

//------------------------------------------------------------------------------
// (14) main
//------------------------------------------------------------------------------
int main() {
  // ��������
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    rand_uint32();
  }

  int mode = 2;  // 0: ��o�p, 1: 1�P�[�X�̂�, 2: �����P�[�X

  if (mode == 0) {
    solve_outer(mode);
  }
  else if (mode == 1) {
    solve_outer(mode, 19);
  }
  else if (mode == 2) {
    double score_sum = 0;
    rep(i, 100) {
      score_sum += solve_outer(mode, i);
      clear_all_multicase();
    }
    cout << "scoreSum = " << score_sum << endl;
  }

  return 0;
}
