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

// ���[�v�̊ȗ����}�N��
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

// �^��`�̃G�C���A�X
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// ���������iXorShift�@�ɂ��[������������j
static uint32_t randxor()
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

// 0�ȏ�1�����̏�����Ԃ������֐�
static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }

// �z����V���b�t������֐��iFisher-Yates�A���S���Y���j
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = randxor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// �����_���f�o�C�X�ƃ����Z���k�E�c�C�X�^�̏������i�g�p����Ă��Ȃ��j
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;   // ���ɑ傫�Ȓl�i�I�[�o�[�t���[�ɒ��Ӂj
const int INT_INF = 1001001001;       // int�^�̔��ɑ傫�Ȓl

// �ړ������̔z��i�g�p����Ă��Ȃ��j
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;  // ���Ԑ����iTime Limit�j
int mode;         // ���s���[�h
std::chrono::steady_clock::time_point startTime, endTime;  // ���Ԍv���p

// ���Ԍv�������Z�b�g����֐�
void ResetTime()
{
  startTime = std::chrono::steady_clock::now();
}

// ���݂̌o�ߎ��Ԃ��擾����֐�
double GetNowTime()
{
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTime - startTime;
  return elapsed.count();
}

// ���̈ʒu��\���\����
struct Point {
  int x;  // �R�̔ԍ�
  int y;  // �R�̒��ł̍����i������̈ʒu�j
};

const int MAX_N = 30;  // �g�p����Ă��Ȃ�

const int n = 200;  // ���̑���
const int m = 10;   // �R�̑���

vector<int> init_stacks[m];   // ������Ԃ̊e�R�ɐς܂ꂽ���̔ԍ����X�g
vector<Point> init_positions; // ������Ԃ̊e���̈ʒu���

// ���̏�Ԃ�\���\����
struct Problem {
  vector<int> stacks[m];      // ���݂̊e�R�̏��
  vector<Point> positions;    // ���݂̊e���̈ʒu
  vector<P> ans;              // �����̋L�^�i�𓚁j
};

// �����̃P�[�X����������ۂɁA������Ԃ�����������֐�
void SetUp()
{
  rep(i, m) {
    init_stacks[i].clear();
  }
  init_positions.clear();
}

// ���͂��󂯎��֐�
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  rep(i, m) {
    init_stacks[i].resize(n / m);  // �e�R�̃T�C�Y��ݒ�
  }
  init_positions.resize(n);  // ���̈ʒu���̃T�C�Y��ݒ�

  // �W�����͂���󂯎��ꍇ
  if (!ifs.is_open()) {
    int _n, _m;
    cin >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> init_stacks[i][j];
        init_stacks[i][j]--;  // 0-indexed�ɕϊ�
      }
    }
  }
  // �t�@�C������󂯎��ꍇ
  else {
    int _n, _m;
    ifs >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        ifs >> init_stacks[i][j];
        init_stacks[i][j]--;  // 0-indexed�ɕϊ�
      }
    }
  }

  // �e���̈ʒu����ݒ�
  rep(i, m) {
    rep(j, n / m) {
      init_positions[init_stacks[i][j]].x = i;  // �R�̔ԍ�
      init_positions[init_stacks[i][j]].y = j;  // ����
    }
  }
}

// �o�̓t�@�C���X�g���[�����J���֐�
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// �X�R�A���v�Z����֐�
int CalcScore(const vector<P>& ans)
{
  vector<int> tmp_stacks[m];          // �ꎞ�I�ȎR�̏��
  vector<Point> tmp_positions = init_positions;  // �ꎞ�I�Ȕ��̈ʒu���

  int cnt = 0;  // �^�яo�������̐�
  rep(i, m)
  {
    tmp_stacks[i] = init_stacks[i];  // ������Ԃ��R�s�[
  }

  int res = 10000;  // �����X�R�A
  rep(i, ans.size())
  {
    int num = ans[i].first;       // ���삷�锠�̔ԍ�
    int x = tmp_positions[num].x; // ���̌��݂̎R�̔ԍ�
    int y = tmp_positions[num].y; // ���̌��݂̍���
    int nx = ans[i].second;       // �ړ���̎R�̔ԍ�
    if (nx == -1) {
      // ����2�F�����^�яo��
      tmp_stacks[x].pop_back();
      cnt++;
    }
    else {
      // ����1�F���Ƃ��̏�̔����ړ�
      int k = tmp_stacks[x].size() - y;
      res -= (k + 1); // ����̗͂��v�Z

      rep(j, tmp_stacks[x].size() - y)
      {
        int num2 = tmp_stacks[x][y + j];
        tmp_positions[num2].x = nx;                 // ���̐V�����R�̔ԍ���ݒ�
        tmp_positions[num2].y = tmp_stacks[nx].size(); // �V����������ݒ�
        tmp_stacks[nx].push_back(num2);             // �ړ���̎R�ɔ���ǉ�
      }
      tmp_stacks[x].resize(y);  // ���̎R����ړ����������폜
    }
  }
  if (cnt != n) return -1;  // �S�Ă̔����^�яo���Ȃ������ꍇ
  return res;
}

// �𓚂��o�͂���֐�
void Output(ofstream& ofs, const vector<P>& ans)
{
  if (mode == 0) {
    // �W���o�͂ɏo��
    rep(i, ans.size()) { cout << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
  else {
    // �t�@�C���ɏo��
    rep(i, ans.size()) { ofs << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
}

// 1�^�[���ł̑�������s����֐�
void ExecuteTurn(Problem& problem, int current_box, int from_stack, const vector<int>& move_targets) {
  int current_y = problem.positions[current_box].y;  // ���݂̃^�[���ŉ^�яo�����̍���

  // �^�яo�����̏�ɂ��锠���ړ�
  for (int k = problem.stacks[from_stack].size() - 1; k > current_y;) {
    int to_stack = move_targets[k];  // �ړ���̎R
    while (k - 1 >= 0 && move_targets[k - 1] == to_stack) k--;  // �����ړ���̔����܂Ƃ߂�

    // �����ړ�
    srep(l, k, problem.stacks[from_stack].size()) {
      int num = problem.stacks[from_stack][l];
      problem.positions[num].x = to_stack;                   // �V�����R�̔ԍ�
      problem.positions[num].y = problem.stacks[to_stack].size(); // �V��������
      problem.stacks[to_stack].push_back(num);               // �ړ���̎R�ɒǉ�
    }

    // ������L�^
    problem.ans.emplace_back(problem.stacks[from_stack][k], to_stack);
    problem.stacks[from_stack].resize(k);  // ���̎R����폜

    k--;
  }

  // ���݂̃^�[���̔����^�яo��
  problem.stacks[from_stack].pop_back();
  problem.ans.emplace_back(current_box, -1);  // ������L�^
}

// ���̈ړ�������肷��֐�
void DecideMoveDestination(const Problem& problem, vector<int>& move_targets, const vector<P>& min_box_per_stack, int from_stack, int position) {
  int bb = problem.stacks[from_stack][position];  // �ړ����锠�̔ԍ�
  int id = 1;
  while (id + 1 < m && min_box_per_stack[id].first < bb) id++;
  move_targets[position] = min_box_per_stack[id].second;  // �ړ���̎R�̔ԍ���ݒ�
}

// �c��̑�����V�~�����[�V��������֐��i�v���C�A�E�g�j
int SimulateRemainingMoves(Problem problem, int current_box, int position, vector<int> move_targets, vector<P> min_box_per_stack, int from_stack) {
  int current_y = problem.positions[current_box].y;
  // �ړ��������
  srep(i, position + 1, problem.stacks[from_stack].size()) {
    DecideMoveDestination(problem, move_targets, min_box_per_stack, from_stack, i);
  }
  // 1�^�[�����̑�������s
  ExecuteTurn(problem, current_box, from_stack, move_targets);

  // �c��̃^�[�����������s
  srep(turn, current_box + 1, n) {
    vector<P> min_boxes_in_stacks(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, problem.stacks[i].size()) {
        minI = min(minI, problem.stacks[i][j]);  // �e�R�̍ŏ��̔��̔ԍ����擾
      }
      min_boxes_in_stacks[i] = P(minI, i);
    }

    sort(min_boxes_in_stacks.begin(), min_boxes_in_stacks.end());  // �ŏ��̔��̔ԍ��Ń\�[�g

    int next_from_stack = min_boxes_in_stacks[0].second;  // ���ɑ��삷��R
    int next_y = problem.positions[turn].y;

    vector<int> next_move_targets(problem.stacks[next_from_stack].size(), -1);
    // ���̈ړ��������
    srep(k, next_y + 1, problem.stacks[next_from_stack].size()) {
      DecideMoveDestination(problem, next_move_targets, min_boxes_in_stacks, next_from_stack, k);
    }

    // 1�^�[�����̑�������s
    ExecuteTurn(problem, turn, next_from_stack, next_move_targets);
  }

  return CalcScore(problem.ans);  // �X�R�A���v�Z���ĕԂ�
}

// Greedy�ȉ�@�����������֐�
Problem GreedySolution() {
  Problem problem;
  rep(i, m) problem.stacks[i] = init_stacks[i];  // ������Ԃ��R�s�[
  problem.positions = init_positions;

  rep(current_box, n) {  // �e�^�[���i�e���j�ɂ���
    vector<P> min_box_per_stack(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, problem.stacks[i].size()) {
        minI = min(minI, problem.stacks[i][j]);  // �e�R�̍ŏ��̔��̔ԍ����擾
      }
      min_box_per_stack[i] = P(minI, i);
    }

    sort(min_box_per_stack.begin(), min_box_per_stack.end());  // �ŏ��̔��̔ԍ��Ń\�[�g

    int from_stack = min_box_per_stack[0].second;  // ���݂̃^�[���ő��삷��R
    int current_y = problem.positions[current_box].y;

    vector<int> move_targets(problem.stacks[from_stack].size(), -1);

    // �e���ɂ��čœK�Ȉړ����T��
    srep(position, current_y + 1, problem.stacks[from_stack].size()) {
      int maxScore = -1;
      int maxId = -1;

      // �e�\�Ȉړ���ɂ��ăv���C�A�E�g
      srep(l, 1, m) {
        move_targets[position] = min_box_per_stack[l].second;
        int tmpScore = SimulateRemainingMoves(problem, current_box, position, move_targets, min_box_per_stack, from_stack);
        if (tmpScore > maxScore) {
          maxScore = tmpScore;
          maxId = min_box_per_stack[l].second;
        }
      }

      move_targets[position] = maxId;  // �ŗǂ̈ړ����ݒ�
    }

    // ���݂̃X�R�A���v�Z
    int score = SimulateRemainingMoves(problem, current_box, problem.stacks[from_stack].size() - 1, move_targets, min_box_per_stack, from_stack);

    // �����_���Ɉړ����ύX���ĒT���i�Ă��Ȃ܂��I�Ȏ�@�j
    if (problem.stacks[from_stack].size() - (current_y + 1) >= 2 && GetNowTime() < TL) {
      rep(iteration, 500) {
        int randomPosition = randxor() % (problem.stacks[from_stack].size() - (current_y + 1)) + current_y + 1;
        int keep = move_targets[randomPosition];

        if (randxor() % 2 == 0) {
          // �����_���Ɉړ����ύX
          int randomIdx = randxor() % (m - 1) + 1;
          move_targets[randomPosition] = min_box_per_stack[randomIdx].second;
        }
        else {
          // �㉺�̔��̈ړ���ɍ��킹��
          int randomDir = -1;
          if (randomPosition == current_y + 1) {
            randomDir = 1;
          }
          else if (randomPosition != current_y + 1 && randomPosition != problem.stacks[from_stack].size() - 1) {
            if (randxor() % 2 == 0) {
              randomDir = 1;
            }
          }
          move_targets[randomPosition] = move_targets[randomPosition + randomDir];
        }

        int tmpScore = SimulateRemainingMoves(problem, current_box, problem.stacks[from_stack].size() - 1, move_targets, min_box_per_stack, from_stack);
        if (tmpScore >= score) {
          score = tmpScore;
        }
        else {
          move_targets[randomPosition] = keep;  // ���P���Ȃ���Ό��ɖ߂�
        }
      }
    }

    // ���肵���ړ���ő�������s
    ExecuteTurn(problem, current_box, from_stack, move_targets);
  }

  return problem;
}

// ���������֐�
ll Solve(int probNum)
{
  ResetTime();  // ���Ԍv�������Z�b�g

  // ������Ԃ�������
  SetUp();

  // ���͂��󂯎��
  Input(probNum);

  // �o�̓t�@�C���X�g���[�����J��
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // �������𐶐�
  auto problem = GreedySolution();

  // �𓚂��o��
  Output(ofs, problem.ans);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore(problem.ans);  // �X�R�A���v�Z
  }
  return score;
}

////////////////////////////////////////////////////////////////////
/*
����

*/
////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));  // �����̎��ݒ�
  while (rand() % 100) {
    randxor();  // ������i�߂Ă���
  }

  mode = 2;  // ���s���[�h�̐ݒ�

  if (mode == 0) {
    Solve(0);  // �P��̃P�[�X������
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);  // �����̃P�[�X������
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
