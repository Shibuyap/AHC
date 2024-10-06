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
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

namespace /* �������C�u���� */
{
  static uint32_t randxor() {
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

  // 0�ȏ�1�����̏������Ƃ闐��
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }

  // �z��V���b�t��
  void FisherYates(int* data, int n) {
    for (int i = n - 1; i >= 0; i--) {
      int j = randxor() % (i + 1);
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

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;
clock_t startTime, endTime;

const int MAX_N = 30;
const int MAX_M = 450;
const int MAX_V = 15;

int n, m, v;
int a[MAX_N][MAX_N];
int b[MAX_N][MAX_N];

double GetNowTime() {
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

// �����P�[�X�񂷂Ƃ��ɓ�����Ԃ������l�ɖ߂�
void SetUp() {}

// ���͎󂯎��
void Input(int problemNum) {
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // �W�����͂���
  if (!ifs.is_open()) {
    cin >> n >> m >> v;
    rep(i, n) {
      string s;
      cin >> s;
      rep(j, n) {
        a[i][j] = s[j] - '0';
      }
    }
    rep(i, n) {
      string t;
      cin >> t;
      rep(j, n) {
        b[i][j] = t[j] - '0';
      }
    }
  }
  else {
    // �t�@�C�����͂���
    ifs >> n >> m >> v;
    rep(i, n) {
      string s;
      ifs >> s;
      rep(j, n) {
        a[i][j] = s[j] - '0';
      }
    }
    rep(i, n) {
      string t;
      ifs >> t;
      rep(j, n) {
        b[i][j] = t[j] - '0';
      }
    }
  }
}

// �o�̓t�@�C���X�g���[���I�[�v��
void OpenOfs(int probNum, ofstream& ofs) {
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4) {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// �X�R�A�v�Z
ll CalcScore() {
  ll res = 0;
  return res;
}

// �����𐶐�
void Initialize() {}

// �𓚏o��
void Output(ofstream& ofs) {
  if (mode == 0) {
  }
  else {
  }
}

ll Solve(int probNum) {
  // �����P�[�X�񂷂Ƃ��ɓ�����Ԃ������l�ɖ߂�
  SetUp();

  // ���͎󂯎��
  Input(probNum);

  // �o�̓t�@�C���X�g���[���I�[�v��
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // �����𐶐�
  Initialize();

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

int main() {
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
