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
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

// �^��`�̃G�C���A�X
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// ���������iXorShift�@�ɂ��[������������j
static uint32_t RandXor()
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

// 0�ȏ�1�����̎�����Ԃ������֐�
static double Rand01() { return (RandXor() + 0.5) * (1.0 / UINT_MAX); }

// l�ȏ�r�����̎������Ƃ闐��
static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

// �z����V���b�t������֐��iFisher-Yates�A���S���Y���j
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// �����_���f�o�C�X�ƃ����Z���k�E�c�C�X�^�̏������i�g�p����Ă��Ȃ��j
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// ���ɑ傫�Ȓl
const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

// �ړ������̔z��
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8; // ���Ԑ����iTime Limit�j
int mode;        // ���s���[�h
std::chrono::steady_clock::time_point startTimeClock; // ���Ԍv���p

// ���Ԍv�������Z�b�g����֐�
void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

// ���݂̌o�ߎ��Ԃ��擾����֐�
double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}


const int MAX_N = 30;

const int n = 50;
int m;
int k;
const int T = 800;

vector<int> sx, sy, tx, ty;

int ansScore;

int best_ansScore;

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

// �����̃P�[�X����������ۂɁA������Ԃ�����������֐�
void SetUp()
{
  ansScore = 0;
  sx.clear();
  sy.clear();
  tx.clear();
  ty.clear();
}

// ���͂��󂯎��֐�
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // �W������
    int nn, tt;
    cin >> nn >> m >> k >> tt;

  }
  else {
    // �t�@�C������
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
ll CalcScore()
{
  ll res = 0;
  return res;
}

// �𓚂��o�͂���֐�
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // �W���o��
  }
  else {
    // �t�@�C���o��
  }
}

// �i�C�[�u�ȉ�@
void Method1()
{

}

// ���������֐�
ll Solve(int problem_num)
{
  ResetTime();

  // �����P�[�X�񂷂Ƃ��ɓ�����Ԃ������l�ɖ߂�
  SetUp();

  // ���͎󂯎��
  Input(problem_num);

  // �o�̓t�@�C���X�g���[���I�[�v��
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // �����𐶐�
  Method1();

  // �𓚂��o��
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

/////////////////////////////////////////////////////////////////////////
/*
����

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
