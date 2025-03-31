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

namespace
{
  const int MAX_UF = 10005;
  int UF_par[MAX_UF];   // �e
  int UF_rank[MAX_UF];  // �؂̐[��
  int UF_cnt[MAX_UF];   // �����钸�_�̌�(�e�̂ݐ�����)

  // ������
  void UF_Init(int n)
  {
    for (int i = 0; i < n; i++) {
      UF_par[i]  = i;
      UF_rank[i] = 0;
      UF_cnt[i]  = 1;
    }
  }

  // ������
  void UF_Init()
  {
    UF_Init(MAX_UF);
  }

  // �؂̍������߂�
  int UF_Find(int x)
  {
    if (UF_par[x] == x) {
      return x;
    }
    else {
      return UF_par[x] = UF_Find(UF_par[x]);
    }
  }

  // x��y�̑�����W���𕹍�
  void UF_Unite(int x, int y)
  {
    x = UF_Find(x);
    y = UF_Find(y);
    if (x == y) return;

    if (UF_rank[x] < UF_rank[y]) {
      UF_par[x] = y;
      UF_cnt[y] += UF_cnt[x];
    }
    else {
      UF_par[y] = x;
      UF_cnt[x] += UF_cnt[y];
      if (UF_rank[x] == UF_rank[y]) UF_rank[x]++;
    }
  }

  // x��y�������W���ɑ����邩�ۂ�
  bool UF_Same(int x, int y)
  {
    return UF_Find(x) == UF_Find(y);
  }

  // x�̑�����W���̃T�C�Y
  int UF_Count(int x)
  {
    return UF_cnt[UF_Find(x)];
  }
}  // namespace



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

const int n = 800;
const int MAX_M = 400;
const int q = 400;
const int MAX_L = 15;

int m, l, w;

int g[MAX_M];
int sort_g[MAX_M];
int arg_g[MAX_M];
int arg_rev_g[MAX_M];

int lx[n], rx[n], ly[n], ry[n];

int true_x[n], true_y[n];

int ansScore;
int ans[n];

int best_ansScore;
int best_ans[n];

struct Edge
{
  int dist;
  int u;
  int v;
};

P GetPoint(int num)
{
  return P((lx[num] + rx[num]) / 2, (ly[num] + ry[num]) / 2);
}

int Distance(const P& p1, const P& p2)
{
  return (p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second);
}

vector<P> MakeZenikigi(const vector<int>& nums)
{
  vector<P> points;
  for (auto num : nums) {
    points.push_back(GetPoint(num));
  }

  vector<Edge> edges;
  rep(i, nums.size())
  {
    srep(j, i + 1, points.size())
    {
      Edge e;
      e.dist = Distance(points[i], points[j]);
      e.u = nums[i];
      e.v = nums[j];
      edges.push_back(e);
    }
  }

  sort(edges.begin(), edges.end());

  vector<P> res;



  return res;
}

void CopyToBest()
{
  best_ansScore = ansScore;
}

void CopyToAns()
{
  ansScore = best_ansScore;
}

bool IsNG(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// �����̃P�[�X����������ۂɁA������Ԃ�����������֐�
void SetUp()
{
  ansScore = 0;
}

// ���͂��󂯎��֐�
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // �W������
    int _n, _q;
    cin >> _n >> m >> _q >> l >> w;
    rep(i, m)cin >> g[i];
    rep(i, n)
    {
      cin >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
  }
  else {
    // �t�@�C������
    int _n, _q;
    ifs >> _n >> m >> _q >> l >> w;
    rep(i, m)ifs >> g[i];
    rep(i, n)
    {
      ifs >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
    rep(i, n)cin >> true_x[i] >> true_y[i];
  }

  vector<P> vp;
  rep(i, m)
  {
    vp.emplace_back(g[i], i);
  }
  sort(vp.begin(), vp.end());
  rep(i, m)
  {
    sort_g[i] = vp[i].first;
    arg_g[i] = vp[i].second;
    arg_rev_g[vp[i].second] = i;
  }
}

vector<P> QueryLocal()
{
  return {};
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
  UF_Init(n);

  if (mode == 0) {
    // �W���o��
    cout << '!' << endl;
    rep(i, m)
    {
      int rev_i = arg_rev_g[i];
      vector<int> points;
      rep(j, n)
      {
        if (ans[j] == rev_i) {
          points.push_back(j);
        }
      }
      for (auto po : points) {
        cout << po << ' ';
      }
      cout << endl;
      vector<P> edges = MakeZenikigi(points);
      for (auto edge : edges) {
        cout << edge.first << ' ' << edge.second << endl;
      }
    }
  }
  else {
    // �t�@�C���o��
    ofs << '!' << endl;
    rep(i, m)
    {
      int rev_i = arg_rev_g[i];
      vector<int> points;
      rep(j, n)
      {
        if (ans[j] == rev_i) {
          points.push_back(j);
        }
      }
      for (auto po : points) {
        ofs << po << ' ';
      }
      ofs << endl;
      vector<P> edges = MakeZenikigi(points);
      for (auto edge : edges) {
        ofs << edge.first << ' ' << edge.second << endl;
      }
    }
  }
}

// �i�C�[�u�ȉ�@
void Method1()
{

}

// �n�C�p�[�p�����[�^
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition;
};

void SimulatedAnnealing(Hypers hypers)
{
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;


  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) break;
    }

    // �߂�
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    int raMode = RandXor() % 100;
    if (raMode < hypers.Partition) {
      // �ߖT���쐬

      // �X�R�A�v�Z
      double tmpScore = CalcScore();

      // �Ă��Ȃ܂�
      double diffScore = (tmpScore - ansScore) * hypers.MultipleValue;
      double prob = exp(diffScore / temp);
      if (prob > Rand01()) {
        // �̗p
        ansScore = tmpScore;

        // Best������������
        if (ansScore > best_ansScore) {
          CopyToBest();
        }
      }
      else {
        // ���ɖ߂�
      }
    }
    else if (raMode < 100) {

    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
  }

  CopyToAns();
}


// ���������֐�
ll Solve(int problem_num, Hypers hypers)
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

  // �Ă��Ȃ܂�
  //SimulatedAnnealing(hypers);

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

  Hypers HYPERS;
  HYPERS.StartTemp = 2048.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 1.0;
  HYPERS.Partition = 50;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 15)
    {
      ll score = Solve(i, HYPERS);
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
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition = RandXor() % 101;

      ll sum = 0;
      srep(i, 0, 15)
      {
        ll score = Solve(i, hypers);
        sum += score;

        // �V�[�h0��������Αł��؂�
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
