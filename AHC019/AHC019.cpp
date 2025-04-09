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
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

const ll INF = 1001001001001001001;
const int dx[6] = { -1, 0, 0, 1, 0, 0 };
const int dy[6] = { 0, -1, 0, 0, 1, 0 };
const int dz[6] = { 0, 0, -1, 0, 0, 1 };

int GetDir(int num)
{
  rep(i, 6)
  {
    if (num & (1 << i)) return i;
  }
  return -1;
}

namespace /* �������C�u���� */
{
  static uint32_t Rand()
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


  static double Rand01() {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace /* �ϐ� */
{
  // ���͗p�ϐ�
  int D;
  bool F[2][20][20];
  bool R[2][20][20];

  // �𓚗p�ϐ�
  double minScore;
  int ans[2][15][15][15];
  int bcount[2][100];

  double real_minScore;
  int real_ans[2][15][15][15];
  int real_bcount[2][100];

  double seed_minScore;
  int seed_ans[2][15][15][15];
  int seed_bcount[2][100];

  // ���̑�
  int methodCount[20][2];

}  // namespace

void MethodCountReset()
{
  rep(i, 20)
  {
    rep(j, 2) { methodCount[i][j] = 0; }
  }
}

// �X�R�A�v�Z
double CalcScore()
{
  double resd = 0;
  double ma, mi;

  srep(j, 1, 3)
  {
    if (bcount[0][j] >= bcount[1][j]) {
      ma = bcount[0][j];
      mi = bcount[1][j];
    }
    else {
      ma = bcount[1][j];
      mi = bcount[0][j];
    }
    resd += (ma - mi) * j;
    resd += mi / j;
  }

  return resd;
}

void NormalClear()
{
  minScore = INF;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) ans[i][j][k][l] = -1;
  rep(i, 2) rep(j, 100) bcount[i][j] = 0;
}

void RealClear()
{
  real_minScore = INF;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) real_ans[i][j][k][l] = -1;
  rep(i, 2) rep(j, 100) real_bcount[i][j] = 0;
}

void SeedClear()
{
  seed_minScore = INF;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) seed_ans[i][j][k][l] = -1;
  rep(i, 2) rep(j, 100) seed_bcount[i][j] = 0;
}

// ���[�J���ŕ����P�[�X�������߂̑S�ď����֐�
void AllClear_MultiCase()
{
  NormalClear();
  RealClear();
  SeedClear();
}

// ������ԍ쐬�i������Ăׂ΃X�^�[�g�ʒu�ɖ߂�邱�Ƃ�z��Areal_minScore���͖߂��Ȃ��j
void Init()
{
  rep(i, 2)
  {
    rep(j, 100) { bcount[i][j] = 0; }
  }
  rep(i, 2)
  {
    rep(x, D)
    {
      rep(y, D)
      {
        rep(z, D)
        {
          if (F[i][z][x] && R[i][z][y]) {
            ans[i][x][y][z] = 0;
            bcount[i][1]++;
          }
          else {
            ans[i][x][y][z] = -1;
          }
        }
      }
    }
  }
}

// ���͎󂯎��i���s����x�����Ă΂�Ȃ����Ƃ�z��j
void Input(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // �W�����͂���
  if (!ifs.is_open()) {
    cin >> D;
    rep(i, 2)
    {
      rep(j, D)
      {
        string s;
        cin >> s;
        rep(k, D) { F[i][j][k] = s[k] - '0'; }
      }
      rep(j, D)
      {
        string s;
        cin >> s;
        rep(k, D) { R[i][j][k] = s[k] - '0'; }
      }
    }
  }
  // �t�@�C�����͂���
  else {
    ifs >> D;
    rep(i, 2)
    {
      rep(j, D)
      {
        string s;
        ifs >> s;
        rep(k, D) { F[i][j][k] = s[k] - '0'; }
      }
      rep(j, D)
      {
        string s;
        ifs >> s;
        rep(k, D) { R[i][j][k] = s[k] - '0'; }
      }
    }
  }

  Init();
}

// �𓚏o��
void Output(int mode, int problemNum)
{
  int ansN[100] = {};
  int ansSum[100] = {};
  rep(j, 100)
  {
    ansN[j] = max(bcount[0][j], bcount[1][j]);
    ansSum[j] = ansN[j];
    if (j > 0) ansSum[j] += ansSum[j - 1];
  }

  int ansPrint[2][15][15][15];
  rep(i, 2)
  {
    rep(j, D)
    {
      rep(k, D)
      {
        rep(l, D) { ansPrint[i][j][k][l] = 0; }
      }
    }
  }
  rep(i, 2)
  {
    int countSum[100] = {};
    srep(j, 1, 100) countSum[j] = ansSum[j - 1];
    rep(j, D)
    {
      rep(k, D)
      {
        rep(l, D)
        {
          if (ansPrint[i][j][k][l] != 0) continue;
          if (ans[i][j][k][l] == -1) {
            ansPrint[i][j][k][l] = 0;
          }
          else if (ans[i][j][k][l] == 0) {
            countSum[1]++;
            ansPrint[i][j][k][l] = countSum[1];
          }
          else {
            countSum[2]++;
            ansPrint[i][j][k][l] = countSum[2];
            int dir = GetDir(ans[i][j][k][l]);
            int nx = j + dx[dir];
            int ny = k + dy[dir];
            int nz = l + dz[dir];
            ansPrint[i][nx][ny][nz] = countSum[2];
          }
        }
      }
    }
  }

  if (mode == 0) {
    cout << ansSum[99] << endl;
    rep(i, 2)
    {
      rep(j, D)
      {
        rep(k, D)
        {
          rep(l, D) { cout << ansPrint[i][j][k][l] << ' '; }
        }
      }
      cout << endl;
    }
  }
  else {
    // �t�@�C���o��
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    ofs << ansSum[99] << endl;
    rep(i, 2)
    {
      rep(j, D)
      {
        rep(k, D)
        {
          rep(l, D) { ofs << ansPrint[i][j][k][l] << ' '; }
        }
      }
      ofs << endl;
    }

    ofs.close();
  }
}

void CopyToReal()
{
  real_minScore = minScore;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) real_ans[i][j][k][l] =
    ans[i][j][k][l];
  rep(i, 2) rep(j, 100) real_bcount[i][j] = bcount[i][j];
}

void CopyToSeed()
{
  seed_minScore = minScore;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) seed_ans[i][j][k][l] =
    ans[i][j][k][l];
  rep(i, 2) rep(j, 100) seed_bcount[i][j] = bcount[i][j];
}

void RollBackFromReal()
{
  minScore = real_minScore;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) ans[i][j][k][l] =
    real_ans[i][j][k][l];
  rep(i, 2) rep(j, 100) bcount[i][j] = real_bcount[i][j];
}

void RollBackFromSeed()
{
  minScore = seed_minScore;
  rep(i, 2) rep(j, D) rep(k, D) rep(l, D) ans[i][j][k][l] =
    seed_ans[i][j][k][l];
  rep(i, 2) rep(j, 100) bcount[i][j] = seed_bcount[i][j];
}

// 1*1*1�̃u���b�N�����
void Method1(double temperature)
{
  int i = Rand() % 2;
  int x = Rand() % D;
  int y = Rand() % D;
  int z = Rand() % D;
  if (ans[i][x][y][z] != -1) return;
  if (!F[i][z][x] || !R[i][z][y]) return;

  methodCount[1][1]++;
  ans[i][x][y][z] = 0;
  bcount[i][1]++;

  double tmpScore = CalcScore();
  double diffScore = minScore - tmpScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[1][0]++;
    minScore = tmpScore;

    if (minScore < real_minScore) {
      CopyToReal();
    }
  }
  else {
    // ���ɖ߂�
    bcount[i][1]--;
    ans[i][x][y][z] = -1;
  }
}

bool CanDelete(int i, int x, int y, int z)
{
  if (F[i][z][x]) {
    bool ok = false;
    rep(k, D)
    {
      if (k == y) continue;
      if (ans[i][x][k][z] != -1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  if (R[i][z][y]) {
    bool ok = false;
    rep(j, D)
    {
      if (j == x) continue;
      if (ans[i][j][y][z] != -1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  return true;
}

// 1*1*1�̃u���b�N������
void Method2(double temperature)
{
  int i = Rand() % 2;
  int x = Rand() % D;
  int y = Rand() % D;
  int z = Rand() % D;
  if (ans[i][x][y][z] != 0) return;
  if (!CanDelete(i, x, y, z)) return;

  methodCount[2][1]++;
  ans[i][x][y][z] = -1;
  bcount[i][1]--;

  double tmpScore = CalcScore();
  double diffScore = minScore - tmpScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[2][0]++;
    minScore = tmpScore;

    if (minScore < real_minScore) {
      CopyToReal();
    }
  }
  else {
    // ���ɖ߂�
    bcount[i][1]++;
    ans[i][x][y][z] = 0;
  }
}

bool IsNG(int x, int y, int z)
{
  if (x < 0 || D <= x || y < 0 || D <= y || z < 0 || D <= z) {
    return true;
  }
  return false;
}

// 1�̃u���b�N������
void Method3(double temperature)
{
  int i = Rand() % 2;
  int x = Rand() % D;
  int y = Rand() % D;
  int z = Rand() % D;
  int dir = Rand() % 6;
  if (ans[i][x][y][z] != 0) return;
  int nx = x + dx[dir];
  int ny = y + dy[dir];
  int nz = z + dz[dir];
  if (IsNG(nx, ny, nz)) return;
  if (ans[i][nx][ny][nz] != 0) return;

  methodCount[3][1]++;
  ans[i][x][y][z] = 1 << dir;
  ans[i][nx][ny][nz] = 1 << ((dir + 3) % 6);
  bcount[i][1] -= 2;
  bcount[i][2]++;

  double tmpScore = CalcScore();
  double diffScore = minScore - tmpScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[3][0]++;
    minScore = tmpScore;

    if (minScore < real_minScore) {
      CopyToReal();
    }
  }
  else {
    // ���ɖ߂�
    bcount[i][2]--;
    bcount[i][1] += 2;
    ans[i][x][y][z] = 0;
    ans[i][nx][ny][nz] = 0;
  }
}

// 2�̃u���b�N�𕪗�
void Method4(double temperature)
{
  int i = Rand() % 2;
  int x = Rand() % D;
  int y = Rand() % D;
  int z = Rand() % D;
  if (ans[i][x][y][z] <= 0) return;
  int dir = GetDir(ans[i][x][y][z]);
  int nx = x + dx[dir];
  int ny = y + dy[dir];
  int nz = z + dz[dir];

  methodCount[4][1]++;
  ans[i][x][y][z] = 0;
  ans[i][nx][ny][nz] = 0;
  bcount[i][1] += 2;
  bcount[i][2]--;

  double tmpScore = CalcScore();
  double diffScore = minScore - tmpScore;

  double prob = exp((double)diffScore / temperature);
  if (prob > Rand01()) {
    methodCount[4][0]++;
    minScore = tmpScore;

    if (minScore < real_minScore) {
      CopyToReal();
    }
  }
  else {
    // ���ɖ߂�
    bcount[i][2]++;
    bcount[i][1] -= 2;
    ans[i][x][y][z] = 1 << dir;
    ans[i][nx][ny][nz] = 1 << ((dir + 3) % 6);
  }
}

double Solve(int mode, int problemNum = 0)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // ������ԍ쐬
  Init();

  // �𒼉��쐬
  minScore = CalcScore();
  CopyToReal();
  CopyToSeed();

  // �V�[�h���
  int seedCount = 100;  // 0�ɂ���ƃV�[�h�쐬���s��Ȃ�
  rep(tei, seedCount)
  {
    startTime = clock();

    // ������Ԃɖ߂�
    Init();
    minScore = CalcScore();

    // �Ă��Ȃ܂�
    endTime = clock();
    double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
    double TL = 4.0 / seedCount;
    double nowProgress = nowTime / TL;
    double startTemperature = 20;
    double endTemperature = 0;
    int loop = 0;
    int rollbackCount = 0;
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        endTime = clock();
        nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        nowProgress = nowTime / TL;
      }
      if (nowProgress > 1.0) break;

      // ���݂̃X�R�A�������Ƃ��͌��ɖ߂�
      if (minScore > real_minScore * 10) {
        // RollBackFromReal();
        // rollbackCount++;
      }

      double temperature =
        startTemperature + (endTemperature - startTemperature) * nowProgress;

      // ���\�b�h�I��
      int ra = Rand() % 100;

      // �e���\�b�h����
      if (ra < 20) {
        Method1(temperature);
      }
      else if (ra < 40) {
        Method2(temperature);
      }
      else if (ra < 70) {
        Method3(temperature);
      }
      else if (ra < 100) {
        Method4(temperature);
      }

    }  // while�������܂Łi�V�[�h�쐬�j

    // �X�R�A���ǂ���΃V�[�h���X�V
    RollBackFromReal();
    if (minScore <= seed_minScore) {
      CopyToSeed();
    }

    // �����ŏ���������̂�����Ώ�������
  }

  // �V�[�h����߂�
  RollBackFromSeed();
  CopyToReal();
  // �Ă��Ȃ܂�
  startTime = clock();
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  double TL = 0.9;
  double nowProgress = nowTime / TL;
  double startTemperature = 2;
  double endTemperature = 0;
  int loop = 0;
  int rollbackCount = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      endTime = clock();
      nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
    }
    if (nowProgress > 1.0) break;

    // ���݂̃X�R�A�������Ƃ��͌��ɖ߂�
    if (minScore > real_minScore * 10) {
      // RollBackFromReal();
      // rollbackCount++;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    // ���\�b�h�I��
    int ra = Rand() % 100;

    // �e���\�b�h����
    if (ra < 20) {
      Method1(temperature);
    }
    else if (ra < 40) {
      Method2(temperature);
    }
    else if (ra < 70) {
      Method3(temperature);
    }
    else if (ra < 100) {
      Method4(temperature);
    }
  }  // while�������܂Łi���C�����[�v�j

  // ��ԃX�R�A�̗ǂ���
  RollBackFromReal();

  // �f�o�b�O���O
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", D = " << D << endl;
    cout << "minScore = " << minScore << endl;
    cout << "loop = " << loop << ", rollbackCount = " << rollbackCount << endl;
    srep(i, 1, 5)
    {
      cout << "Method" << i << " = " << methodCount[i][0] << " / "
        << methodCount[i][1] << endl;
    }
    cout << endl;
  }

  return minScore;
}

double SolveOuter(int mode, int problemNum = 0)
{
  // ���͎󂯎��
  Input(problemNum);

  double score = Solve(mode, problemNum);

  // �𓚂̏o��
  Output(mode, problemNum);

  return score;
}

int main()
{
  // ��������
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

  int mode = 0;

  // ��o�p
  if (mode == 0) {
    SolveOuter(mode);
  }
  // 1�P�[�X����
  else if (mode == 1) {
    SolveOuter(mode, 19);
  }
  // �����P�[�X����
  else if (mode == 2) {
    double scoreSum = 0;
    rep(i, 100)
    {
      scoreSum += SolveOuter(mode, i);
      AllClear_MultiCase();
    }
    cout << "scoreSum = " << scoreSum << endl;
  }

  return 0;
}
