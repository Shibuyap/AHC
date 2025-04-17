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

//=== Macros (�K�v�ɉ����ďC��) ===//
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;

//------------------------------------------------------------------------------
// (1) Union-Find�֘A�Funnamed namespace
//------------------------------------------------------------------------------
namespace
{
  // �萔 (UPPER_SNAKE_CASE)
  const int UNION_FIND_MAX = 50000;

  // Union-Find�Ŏg�p����z�� (snake_case)
  int uf_parent[UNION_FIND_MAX];
  int uf_rank[UNION_FIND_MAX];
  int uf_count[UNION_FIND_MAX];

  // n�v�f�ŏ�����
  void uf_init()
  {
    for (int i = 0; i < UNION_FIND_MAX; i++) {
      uf_parent[i] = i;
      uf_rank[i] = 0;
      uf_count[i] = 1;
    }
  }

  // �؂̍������߂�
  int uf_find(int x)
  {
    if (uf_parent[x] == x) {
      return x;
    }
    else {
      return uf_parent[x] = uf_find(uf_parent[x]);
    }
  }

  // x��y�̑�����W���𕹍�
  void uf_unite(int x, int y)
  {
    x = uf_find(x);
    y = uf_find(y);
    if (x == y) return;

    if (uf_rank[x] < uf_rank[y]) {
      uf_parent[x] = y;
      uf_count[y] += uf_count[x];
    }
    else {
      uf_parent[y] = x;
      uf_count[x] += uf_count[y];
      if (uf_rank[x] == uf_rank[y]) uf_rank[x]++;
    }
  }

  // x��y�������W���ɑ����邩�ۂ�
  bool uf_same(int x, int y)
  {
    return uf_find(x) == uf_find(y);
  }
}

//------------------------------------------------------------------------------
// (2) �֗��n�E�萔�֘A�Funnamed namespace
//------------------------------------------------------------------------------
namespace
{
  const int INF = 1001001001;

  // �ړ��p�̍��� (UPPER_SNAKE_CASE)
  const int DELTA_X[4] = { -1, 0, 1, 0 };
  const int DELTA_Y[4] = { 0, -1, 0, 1 };

  // (��) �ړ����������������� (�g���Ă��Ȃ��\������)
  const char MOVE_CHARS[4] = { 'U', 'L', 'D', 'R' };
}

//------------------------------------------------------------------------------
// (3) �������C�u�����Funnamed namespace
//------------------------------------------------------------------------------
namespace
{
  static uint32_t rand_uint32()
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

  static double rand_double_01()
  {
    return (rand_uint32() + 0.5) * (1.0 / UINT_MAX);
  }
}

//------------------------------------------------------------------------------
// (4) ��v�O���[�o���ϐ��Funnamed namespace
//------------------------------------------------------------------------------
namespace
{
  // ���Ƃ̃R�[�h�� MIN_ROCK, MAX_ROCK �͎g���Ă��Ȃ����ߒu���Ă����܂�
  const int MIN_ROCK = 10;
  const int MAX_ROCK = 5000;

  // ��: CC (�U���R�X�g�̌��)
  const int C_OPTIONS[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  // ���s���[�h (0: ��o�p / ����ȊO: ���[�J��)
  int execution_mode = 0;
  ofstream output_file_stream;

  // ���͊֘A (��N, W, K, C)
  int n_size;          // N
  int w_count;         // W
  int k_count;         // K
  int attack_cost;     // C

  // �����̍��W (��: a[i], b[i])
  int water_x[20];
  int water_y[20];

  // �Ƃ̍��W (��: c[i], d[i])
  int house_x[20];
  int house_y[20];

  // ��̋��x�}�b�v (��: S)
  int rock_strength[210][210];

  // �e�}�X�����Ă��邩�ǂ����Ȃ� (��: f, minS, maxS)
  int is_broken[210][210];
  int min_strength[210][210];
  int max_strength[210][210];

  // �e�Ƃ��ł��߂��������L�^ (��: nearX, nearY)
  int nearest_water_x[20];
  int nearest_water_y[20];

  // �X�R�A�v�Z�p: HP (��: hp)
  int health_points;

  // �U���͂̌�� (��: attack_power[8])
  const int ATTACK_POWER_VALUES[8] = { 50, 50, 50, 50, 50, 100, 100, 100 };

  // ���ۂɎg���U���� (��: ATTACK_POWER)
  int attack_power_global = 100;
}

//------------------------------------------------------------------------------
// (5) �w���p�[�֐�
//------------------------------------------------------------------------------
inline bool is_out_of_bounds(int x, int y)
{
  return (x < 0 || x >= n_size || y < 0 || y >= n_size);
}

// x,y�����ƂȂ����Ă��邩 (��: IsUniteWater)
inline bool is_united_with_water(int x, int y)
{
  // UF_MAX - 1 = UNION_FIND_MAX - 1
  // x * n_size + y �� UNION_FIND_MAX - 1 �Ɠ����A��������
  return uf_same(x * n_size + y, UNION_FIND_MAX - 1);
}

// �}���n�b�^������ (��: Manhattan)
inline int manhattan_distance(int x1, int y1, int x2, int y2)
{
  return abs(x1 - x2) + abs(y1 - y2);
}

//------------------------------------------------------------------------------
// (6) �X�R�A�v�Z (��: CalcScore)
//------------------------------------------------------------------------------
int calc_score()
{
  return health_points;
}

//------------------------------------------------------------------------------
// (7) �U������ (��: Attack)
//------------------------------------------------------------------------------
int attack_cell(int x, int y, int power)
{
  // min_strength������max_strength�ɍX�V
  min_strength[x][y] = max_strength[x][y];
  // max_strength�ɍ���̒ǉ��U���ʂ����Z
  max_strength[x][y] += power;
  // HP(�X�R�A)�����Z
  health_points += attack_cost + power;

  int res = 0;
  if (execution_mode == 0) {
    // ��o�p: �W���o�͂ɍU�����e���o���A����(��ꂽ���ǂ���)���󂯎��
    cout << x << ' ' << y << ' ' << power << endl;
    fflush(stdout);
    cin >> res;
  }
  else {
    // ���[�J�����s: rock_strength�����炵�ĉ�ꂽ������
    output_file_stream << x << ' ' << y << ' ' << power << endl;
    rock_strength[x][y] -= power;
    if (rock_strength[x][y] <= 0) {
      res = 1; // ��ꂽ
    }
  }

  // ��ꂽ�ꍇ
  if (res != 0) {
    is_broken[x][y] = 1;

    // ���̓}�X�����łɉ��Ă���� Union-Find �łȂ�
    rep(i, 4) {
      int nx = x + DELTA_X[i];
      int ny = y + DELTA_Y[i];
      if (is_out_of_bounds(nx, ny)) continue;
      if (is_broken[nx][ny] != 0) {
        uf_unite(x * n_size + y, nx * n_size + ny);
      }
    }

    // ���������}�X��������AUF_MAX-1 �ƂȂ�
    rep(i, w_count) {
      if (x == water_x[i] && y == water_y[i]) {
        uf_unite(x * n_size + y, UNION_FIND_MAX - 1);
      }
    }

    // �����ōēx�u�S�ẲƂ��q���������ǂ����v�𔻒肵�Ă���炵��
    res = 2;
    rep(i, k_count) {
      if (!is_united_with_water(house_x[i], house_y[i])) {
        res = 1; // �܂����ڑ��̉Ƃ�����
      }
    }
  }

  return res;
}

//------------------------------------------------------------------------------
// (8) �A���U���ŉ���܂Ō@�葱���� (��: Challenge)
//------------------------------------------------------------------------------
int challenge_cell(int x, int y, int power)
{
  int res = is_broken[x][y];
  while (res == 0) {
    res = attack_cell(x, y, power);
  }
  return res;
}

//------------------------------------------------------------------------------
// (9) �����P�[�X�p�ɂ��ׂăN���A (��: AllClear_MultiCase)
//------------------------------------------------------------------------------
void clear_all_multicase()
{
  // ��������(�T���v��)
}

//------------------------------------------------------------------------------
// (10) ������ԍ쐬 (��: Init)
//------------------------------------------------------------------------------
void init(int problem_num)
{
  uf_init();
  health_points = 0;

  // ��ԃN���A
  rep(i, n_size) {
    rep(j, n_size) {
      is_broken[i][j] = 0;
      min_strength[i][j] = 0;
      max_strength[i][j] = 0;
    }
  }

  // �U���R�X�gC�ɉ����čU���͂�I��
  rep(i, 8) {
    if (attack_cost == C_OPTIONS[i]) {
      attack_power_global = ATTACK_POWER_VALUES[i];
    }
  }

  // ���[�J�����s���̓t�@�C���o�͂��I�[�v��
  if (execution_mode != 0) {
    string file_name_ofs = "./out/";
    {
      // problem_num��4�������ăt�@�C���������
      int tmp = problem_num;
      string str_num;
      rep(i, 4) {
        str_num += char((tmp % 10) + '0');
        tmp /= 10;
      }
      reverse(str_num.begin(), str_num.end());
      file_name_ofs += str_num + ".txt";
    }

    output_file_stream.open(file_name_ofs);
  }
}

//------------------------------------------------------------------------------
// (11) ���͎󂯎�� (��: Input)
//------------------------------------------------------------------------------
void read_input(int problem_num)
{
  // �t�@�C��������
  string file_name_ifs = "./in/";
  {
    int tmp = problem_num;
    string str_num;
    rep(i, 4) {
      str_num += char((tmp % 10) + '0');
      tmp /= 10;
    }
    reverse(str_num.begin(), str_num.end());
    file_name_ifs += str_num + ".txt";
  }

  ifstream ifs(file_name_ifs);

  // ���s���[�h=0�܂��̓t�@�C�����J���Ȃ�������W�����͂���ǂ�
  if (execution_mode == 0 || !ifs.is_open()) {
    cin >> n_size >> w_count >> k_count >> attack_cost;
    rep(i, w_count) {
      cin >> water_x[i] >> water_y[i];
    }
    rep(i, k_count) {
      cin >> house_x[i] >> house_y[i];
    }
  }
  else {
    // ���[�J�����̓t�@�C������
    ifs >> n_size >> w_count >> k_count >> attack_cost;
    rep(i, n_size) {
      rep(j, n_size) {
        ifs >> rock_strength[i][j];
      }
    }
    rep(i, w_count) {
      ifs >> water_x[i] >> water_y[i];
    }
    rep(i, k_count) {
      ifs >> house_x[i] >> house_y[i];
    }
  }
}

//------------------------------------------------------------------------------
// (12) �𓚏o�� (��: Output)
//------------------------------------------------------------------------------
void write_output(int /*problem_num*/)
{
  // ���[�J�����s���̓t�@�C���o�̓N���[�Y
  if (execution_mode != 0) {
    output_file_stream.close();
  }
}

//------------------------------------------------------------------------------
// (13) ���ۂ̏��� (��: Solve)
//------------------------------------------------------------------------------
int solve_problem(int problem_num = 0)
{
  clock_t start_time = clock();
  clock_t end_time = clock();
  (void)start_time; // �g��Ȃ��Ȃ�void�L���X�g�ŏ���
  (void)end_time;

  // ������
  init(problem_num);

  // �e�Ƃ̈�ԋ߂�������T��
  rep(i, k_count) {
    int dist = INF;
    rep(j, w_count) {
      int mdist = manhattan_distance(house_x[i], house_y[i], water_x[j], water_y[j]);
      if (mdist < dist) {
        dist = mdist;
        nearest_water_x[i] = water_x[j];
        nearest_water_y[i] = water_y[j];
      }
    }
  }

  // �Ƃ𐅌��ɂȂ���
  rep(i, k_count) {
    int phase = 0;
    int now_x = house_x[i], now_y = house_y[i];
    int next_x = -1, next_y = -1;
    int dir = -1;

    // ���̉Ƃ��܂����ƌq�����Ă��Ȃ��Ԃ͌@�葱����
    while (!is_united_with_water(house_x[i], house_y[i])) {
      if (phase == 0) {
        // �Ƃ̃}�X���@��
        challenge_cell(now_x, now_y, attack_power_global);
        phase = 1;
      }
      else if (phase == 1) {
        // ���́u�^�[�Q�b�g�n�_(next_x, next_y)�v�����߂�
        int diff_x = nearest_water_x[i] - now_x;
        if (abs(diff_x) > 20) {
          diff_x = (diff_x > 0) ? 20 : -20;
        }
        int diff_y = nearest_water_y[i] - now_y;
        if (abs(diff_y) > 20) {
          diff_y = (diff_y > 0) ? 20 : -20;
        }

        // ���ɉ��Đ��H���ʂ��Ă���}�X���߂��ɂ���΁A�������D��
        int diff_sum = abs(diff_x) + abs(diff_y);
        srep(k, 1, diff_sum) {
          bool found = false;
          rep(j, 4) {
            int nx = now_x + DELTA_X[j] * k;
            int ny = now_y + DELTA_Y[j] * k;
            if (!is_out_of_bounds(nx, ny) && is_united_with_water(nx, ny)) {
              diff_x = nx - now_x;
              diff_y = ny - now_y;
              found = true;
              break;
            }
          }
          if (found) break;
        }

        // x���� or y���������̏ꍇ
        if (diff_x == 0) {
          next_x = now_x;
          next_y = now_y + diff_y;
          challenge_cell(next_x, next_y, attack_power_global);
        }
        else if (diff_y == 0) {
          next_x = now_x + diff_x;
          next_y = now_y;
          challenge_cell(next_x, next_y, attack_power_global);
        }
        else {
          // �΂߂ɐi�ނƂ��͂܂� x�����Ɍ@���āA���ꂩ�� y�����Ɍ@�� or �t��
          int nx1 = now_x + diff_x;
          int ny1 = now_y;
          challenge_cell(nx1, ny1, attack_power_global);

          int nx2 = now_x;
          int ny2 = now_y + diff_y;
          challenge_cell(nx2, ny2, attack_power_global);

          // �ǂ���̃}�X�̂ق��� max_strength �����������Ŏ��̌��ݒn�����߂�
          if (max_strength[nx1][ny1] <= max_strength[nx2][ny2]) {
            next_x = nx1;
            next_y = ny1;
          }
          else {
            next_x = nx2;
            next_y = ny2;
          }
        }

        // �ړ�����dir�����߂� (DELTA_X/DELTA_Y�̓Y����Ή�)
        if (next_x - now_x < 0) {
          dir = 0; // �����(DELTA_X[0],DELTA_Y[0])
        }
        else if (next_x - now_x > 0) {
          dir = 2; // ������(DELTA_X[2],DELTA_Y[2])
        }
        else if (next_y - now_y < 0) {
          dir = 1; // ������(DELTA_X[1],DELTA_Y[1])
        }
        else if (next_y - now_y > 0) {
          dir = 3; // �E����(DELTA_X[3],DELTA_Y[3])
        }

        phase = 2;
      }
      else if (phase == 2) {
        // next_x, next_y �ɓ��B�ςȂ�ēx�^�[�Q�b�g���Čv�Z
        if (now_x == next_x && now_y == next_y) {
          phase = 1;
        }
        else {
          // next_x, next_y �Ɍ������Ĉ���i��
          int nx = now_x + DELTA_X[dir];
          int ny = now_y + DELTA_Y[dir];
          if (max_strength[nx][ny] == 0) {
            // �܂��U������Ă��Ȃ��}�X�Ȃ�C���͂̏󋵂ɉ������U���ʌv�Z
            int d1 = manhattan_distance(now_x, now_y, nx, ny);
            int d2 = manhattan_distance(nx, ny, next_x, next_y);
            int p1 = min_strength[now_x][now_y];
            int p2 = min_strength[next_x][next_y];
            int power = (d2 * p1 + d1 * p2) / (d1 + d2);
            power = max(power, 10);
            attack_cell(nx, ny, power);
          }
          else {
            // ���łɍU�����ꂽ���Ƃ̂���}�X�Ȃ�A����܂Ō@��
            challenge_cell(nx, ny, attack_power_global);
            now_x = nx;
            now_y = ny;
          }
        }
      }
    }
  }

  // ���[�J�����s���̓f�o�b�O�o��
  if (execution_mode != 0) {
    cerr << "problem_num = " << problem_num
      << ", health_points = " << health_points << endl;
  }

  return calc_score();
}

//------------------------------------------------------------------------------
// (14) ���b�p�֐� (��: SolveOuter)
//------------------------------------------------------------------------------
int solve_outer(int problem_num = 0)
{
  // ���͎󂯎��
  read_input(problem_num);

  // ���s
  int score = solve_problem(problem_num);

  // �o��
  write_output(problem_num);

  return score;
}

//------------------------------------------------------------------------------
// (15) main
//------------------------------------------------------------------------------
int main()
{
  srand((unsigned)time(NULL));
  // ��������
  while (rand() % 100) {
    rand_uint32();
  }

  // execution_mode = 0 => ��o�p
  execution_mode = 2;

  if (execution_mode == 0) {
    // ��o�p: �P��P�[�X
    solve_outer();
  }
  else if (execution_mode == 1) {
    // 1�P�[�X�̂݃t�@�C���ǂ݁E����
    solve_outer(0);
  }
  else if (execution_mode == 2) {
    // �����P�[�X
    ll score_sum = 0;
    rep(i, 100) {
      score_sum += solve_outer(i);
      clear_all_multicase();
    }
    cout << "scoreSum = " << score_sum << endl;
  }

  return 0;
}
