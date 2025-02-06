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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)
using namespace std;
typedef long long int ll;
typedef pair<ll, ll> P;
typedef pair<P, P> PP;

namespace /* �������C�u���� */
{
  static uint32_t randxor()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t        = x ^ (x << 11);
    x        = y;
    y        = z;
    z        = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0�ȏ�1�����̏������Ƃ闐��
  static double rand01()
  {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
  }

  // �z��V���b�t��
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j   = randxor() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

// �z��V���b�t��
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF      = 1001001001001001001;
const int INT_INF = 1001001001;
const int MA      = 1000000000;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
std::chrono::steady_clock::time_point startTime, endTime;

void ResetTime()
{
  startTime = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTime - startTime;
  return elapsed.count();
}

const int n     = 1000;
const int MAX_M = 5000;

ll a[n], b[n];
ll L;

ll ans[5500][4];
int ansCount;
ll ansCost;

ll real_ans[5500][4];
int real_ansCount;
ll real_ansCost;

void CopyToReal()
{
  real_ansCount = ansCount;
  real_ansCost  = ansCost;

  rep(i, ansCount)
  {
    rep(j, 4)
    {
      real_ans[i][j] = ans[i][j];
    }
  }
}

void CopyToAns()
{
  ansCount = real_ansCount;
  ansCost  = real_ansCost;

  rep(i, ansCount)
  {
    rep(j, 4)
    {
      ans[i][j] = real_ans[i][j];
    }
  }
}

// �����P�[�X�񂷂Ƃ��ɓ�����Ԃ������l�ɖ߂�
void SetUp()
{
  ansCount = 0;
  ansCost  = 0;
}

// ���͎󂯎��
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // �W�����͂���
  if (!ifs.is_open()) {
    int nnn;
    cin >> nnn;
    rep(i, n)
    {
      cin >> a[i] >> b[i];
    }
  }
  // �t�@�C�����͂���
  else {
    int nnn;
    ifs >> nnn;
    rep(i, n)
    {
      ifs >> a[i] >> b[i];
    }
  }

  L = 0;
  rep(i, n)
  {
    L = max(L, a[i]);
    L = max(L, b[i]);
  }

  vector<P> vp;
  rep(i, n)
  {
    vp.push_back(P(a[i], b[i]));
  }
  sort(vp.begin(), vp.end());
  rep(i, n)
  {
    a[i] = vp[i].first;
    b[i] = vp[i].second;
  }
}

// �o�̓t�@�C���X�g���[���I�[�v��
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// �X�R�A�v�Z
ll CalcScore()
{
  ll res = round(1.0 * n * L / (ansCost + 1) * 1000000.0);
  return res;
}

// �𓚏o��
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << ansCount << endl;
    rep(i, ansCount)
    {
      rep(j, 4)
      {
        cout << ans[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    ofs << ansCount << endl;
    rep(i, ansCount)
    {
      rep(j, 4)
      {
        ofs << ans[i][j] << ' ';
      }
      ofs << endl;
    }
  }
}

void AddOneStepToAns(int x1, int y1, int x2, int y2)
{
  ans[ansCount][0] = x1;
  ans[ansCount][1] = y1;
  ans[ansCount][2] = x2;
  ans[ansCount][3] = y2;
  ansCount++;
  ansCost += (x2 - x1) + (y2 - y1);
}

// 23�ʉ�@
void Method23()
{
  vector<P> points;
  vector<int> used;
  priority_queue<PP> que;
  rep(i, n)
  {
    points.emplace_back(a[i], b[i]);
    used.push_back(0);
  }
  rep(i, n)
  {
    srep(j, i + 1, n)
    {
      PP pp;
      pp.first.first = points[i].first + points[i].second + points[j].first + points[j].second
        - abs(points[i].first - points[j].first) - abs(points[i].second - points[j].second);
      pp.second.first = i;
      pp.second.second = j;
      que.push(pp);
    }
  }

  ansCount = 0;
  ansCost = 0;

  while (que.size()) {
    PP pp = que.top();
    que.pop();

    int idx1 = pp.second.first;
    int idx2 = pp.second.second;
    if (used[idx1] || used[idx2]) continue;

    int newX = min(points[idx1].first, points[idx2].first);
    int newY = min(points[idx1].second, points[idx2].second);
    AddOneStepToAns(newX, newY, points[idx1].first, points[idx1].second);
    AddOneStepToAns(newX, newY, points[idx2].first, points[idx2].second);

    used[idx1] = 1;
    used[idx2] = 1;

    rep(i, points.size())
    {
      if (used[i])continue;

      PP ppp;
      ppp.first.first = points[i].first + points[i].second + newX + newY
        - abs(points[i].first - newX) - abs(points[i].second - newY);
      ppp.second.first = i;
      ppp.second.second = points.size();
      que.push(ppp);
    }

    points.emplace_back(newX, newY);
    used.push_back(0);
  }

  rep(i, ansCount / 2)
  {
    rep(j, 4)
    {
      swap(ans[i][j], ans[ansCount - 1 - i][j]);
    }
  }
}

ll Solve(int probNum)
{
  ResetTime();

  // �����P�[�X�񂷂Ƃ��ɓ�����Ԃ������l�ɖ߂�
  SetUp();

  // ���͎󂯎��
  Input(probNum);

  // �o�̓t�@�C���X�g���[���I�[�v��
  ofstream ofs;
  OpenOfs(probNum, ofs);

  Method23();

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
    randxor();
  }

  mode = 0;

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
        cout << endl;
      }
    }
  }

  return 0;
}
