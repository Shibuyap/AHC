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

vector<int> init_b[m];   // ������Ԃ̊e�R�ɐς܂ꂽ���̔ԍ����X�g
vector<Point> init_c;    // ������Ԃ̊e���̈ʒu���

// ���̏�Ԃ�\���\����
struct Problem {
  vector<int> b[m];      // ���݂̊e�R�̏��
  vector<Point> c;       // ���݂̊e���̈ʒu
  vector<P> ans;         // �����̋L�^�i�𓚁j
};

// �s�g�p�̊֐��i�����炭�����I�ȍŗǉ��̕ۑ��p�j
void CopyToBest()
{
}

// �s�g�p�̊֐��i�����炭���݂̉𓚂��R�s�[���邽�߂̂��́j
void CopyToAns()
{
}

// �����̃P�[�X����������ۂɁA������Ԃ�����������֐�
void SetUp()
{
  rep(i, m) {
    init_b[i].clear();
  }
  init_c.clear();
}

// ���͂��󂯎��֐�
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  rep(i, m) {
    init_b[i].resize(n / m);  // �e�R�̃T�C�Y��ݒ�
  }
  init_c.resize(n);  // ���̈ʒu���̃T�C�Y��ݒ�

  // �W�����͂���󂯎��ꍇ
  if (!ifs.is_open()) {
    int _n, _m;
    cin >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> init_b[i][j];
        init_b[i][j]--;  // 0-indexed�ɕϊ�
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
        ifs >> init_b[i][j];
        init_b[i][j]--;  // 0-indexed�ɕϊ�
      }
    }
  }

  // �e���̈ʒu����ݒ�
  rep(i, m) {
    rep(j, n / m) {
      init_c[init_b[i][j]].x = i;  // �R�̔ԍ�
      init_c[init_b[i][j]].y = j;  // ����
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
  vector<int> tmp_b[m];          // �ꎞ�I�ȎR�̏��
  vector<Point> tmp_c = init_c;  // �ꎞ�I�Ȕ��̈ʒu���

  int cnt = 0;  // �^�яo�������̐�
  rep(i, m)
  {
    tmp_b[i] = init_b[i];  // ������Ԃ��R�s�[
  }

  int res = 10000;  // �����X�R�A
  rep(i, ans.size())
  {
    int num = ans[i].first;  // ���삷�锠�̔ԍ�
    int x = tmp_c[num].x;    // ���̌��݂̎R�̔ԍ�
    int y = tmp_c[num].y;    // ���̌��݂̍���
    int nx = ans[i].second;  // �ړ���̎R�̔ԍ�
    if (nx == -1) {
      // ����2�F�����^�яo��
      tmp_b[x].pop_back();
      cnt++;
    }
    else {
      // ����1�F���Ƃ��̏�̔����ړ�
      res--;
      rep(j, tmp_b[x].size() - y)
      {
        int num2 = tmp_b[x][y + j];
        tmp_c[num2].x = nx;           // ���̐V�����R�̔ԍ���ݒ�
        tmp_c[num2].y = tmp_b[nx].size();  // �V����������ݒ�
        tmp_b[nx].push_back(num2);    // �ړ���̎R�ɔ���ǉ�
        res--;
      }
      tmp_b[x].resize(y);  // ���̎R����ړ����������폜
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
void MoveOneTurn(Problem& prob, int turn, int from, const vector<int>& idx) {
  int turnY = prob.c[turn].y;  // ���݂̃^�[���ŉ^�яo�����̍���

  // �^�яo�����̏�ɂ��锠���ړ�
  for (int k = prob.b[from].size() - 1; k > turnY;) {
    int to = idx[k];  // �ړ���̎R
    while (idx[k - 1] == to) k--;  // �����ړ���̔����܂Ƃ߂�

    // �����ړ�
    srep(l, k, prob.b[from].size()) {
      int num = prob.b[from][l];
      prob.c[num].x = to;           // �V�����R�̔ԍ�
      prob.c[num].y = prob.b[to].size();  // �V��������
      prob.b[to].push_back(num);    // �ړ���̎R�ɒǉ�
    }

    // ������L�^
    prob.ans.emplace_back(prob.b[from][k], to);
    prob.b[from].resize(k);  // ���̎R����폜

    k--;
  }

  // ���݂̃^�[���̔����^�яo��
  prob.b[from].pop_back();
  prob.ans.emplace_back(turn, -1);  // ������L�^
}

// ���̈ړ�������肷��֐�
void DecideIdx(const Problem& prob, vector<int>& idx, const vector<P>& minIs, int from, int k) {
  int bb = prob.b[from][k];  // �ړ����锠�̔ԍ�
  int id = 1;
  while (id + 1 < m && minIs[id].first < bb) id++;
  idx[k] = minIs[id].second;  // �ړ���̎R�̔ԍ���ݒ�
}

// �c��̑�����V�~�����[�V��������֐��i�v���C�A�E�g�j
int PlayOut(Problem prob, int _turn, int _k, vector<int> _idx, vector<P> _minIs, int _from) {
  int turnY = prob.c[_turn].y;
  // �ړ��������
  srep(k, _k + 1, prob.b[_from].size()) {
    DecideIdx(prob, _idx, _minIs, _from, k);
  }
  // 1�^�[�����̑�������s
  MoveOneTurn(prob, _turn, _from, _idx);

  // �c��̃^�[�����������s
  srep(turn, _turn + 1, n) {
    vector<P> minIs(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, prob.b[i].size()) {
        minI = min(minI, prob.b[i][j]);  // �e�R�̍ŏ��̔��̔ԍ����擾
      }
      minIs[i] = P(minI, i);
    }

    sort(minIs.begin(), minIs.end());  // �ŏ��̔��̔ԍ��Ń\�[�g

    int from = minIs[0].second;  // ���ɑ��삷��R
    int turnY = prob.c[turn].y;

    vector<int> idx(prob.b[from].size(), -1);
    // ���̈ړ��������
    srep(k, turnY + 1, prob.b[from].size()) {
      DecideIdx(prob, idx, minIs, from, k);
    }

    // 1�^�[�����̑�������s
    MoveOneTurn(prob, turn, from, idx);
  }

  return CalcScore(prob.ans);  // �X�R�A���v�Z���ĕԂ�
}

// kotatsugame������×~�@�����������֐�
Problem Method2() {
  Problem prob;
  rep(i, m) prob.b[i] = init_b[i];  // ������Ԃ��R�s�[
  prob.c = init_c;

  rep(turn, n) {  // �e�^�[���i�e���j�ɂ���
    vector<P> minIs(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, prob.b[i].size()) {
        minI = min(minI, prob.b[i][j]);  // �e�R�̍ŏ��̔��̔ԍ����擾
      }
      minIs[i] = P(minI, i);
    }

    sort(minIs.begin(), minIs.end());  // �ŏ��̔��̔ԍ��Ń\�[�g

    int from = minIs[0].second;  // ���݂̃^�[���ő��삷��R
    int turnY = prob.c[turn].y;

    vector<int> keepIdx(prob.b[from].size(), -1);
    // ���̈ړ��������
    srep(k, turnY + 1, prob.b[from].size()) {
      DecideIdx(prob, keepIdx, minIs, from, k);
    }

    vector<int> idx(prob.b[from].size(), -1);
    // �e���ɂ��čœK�Ȉړ����T��
    srep(k, turnY + 1, prob.b[from].size()) {
      int maxScore = -1;
      int maxId = -1;

      // �e�\�Ȉړ���ɂ��ăv���C�A�E�g
      srep(l, 1, m) {
        idx[k] = minIs[l].second;
        int tmpScore = PlayOut(prob, turn, k, idx, minIs, from);
        if (tmpScore > maxScore) {
          maxScore = tmpScore;
          maxId = minIs[l].second;
        }
      }

      idx[k] = maxId;  // �ŗǂ̈ړ����ݒ�
    }

    // ���݂̃X�R�A���v�Z
    int score = PlayOut(prob, turn, prob.b[from].size() - 1, idx, minIs, from);

    // �����_���Ɉړ����ύX���ĒT���i�Ă��Ȃ܂��I�Ȏ�@�j
    if (prob.b[from].size() - (turnY + 1) >= 2 && GetNowTime() < TL) {
      rep(aespa, 500) {
        int randomK = randxor() % (prob.b[from].size() - (turnY + 1)) + turnY + 1;
        int keep = idx[randomK];

        if (randxor() % 2 == 0) {
          // �����_���Ɉړ����ύX
          int randomIdx = randxor() % (m - 1) + 1;
          idx[randomK] = minIs[randomIdx].second;
        }
        else {
          // �㉺�̔��̈ړ���ɍ��킹��
          int randomDir = -1;
          if (randomK == turnY + 1) {
            randomDir = 1;
          }
          else if (randomK != turnY + 1 && randomK != prob.b[from].size() - 1) {
            if (randxor() % 2 == 0) {
              randomDir = 1;
            }
          }
          idx[randomK] = idx[randomK + randomDir];
        }

        int tmpScore = PlayOut(prob, turn, prob.b[from].size() - 1, idx, minIs, from);
        if (tmpScore >= score) {
          score = tmpScore;
        }
        else {
          idx[randomK] = keep;  // ���P���Ȃ���Ό��ɖ߂�
        }
      }
    }

    // ���肵���ړ���ő�������s
    MoveOneTurn(prob, turn, from, idx);
  }

  return prob;
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
  auto problem = Method2();

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

/////////////////////////////////////////////////////////////////////
/*
����

*/
/////////////////////////////////////////////////////////////////////
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
