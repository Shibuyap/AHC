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

namespace /* Union Find */
{
  const int UF_MAX = 50000;
  int UF_par[UF_MAX];   // �e
  int UF_rank[UF_MAX];  // �؂̐[��
  int UF_cnt[UF_MAX];   // �����钸�_�̌�(�e�̂ݐ�����)

  // n�v�f�ŏ�����
  void UFinit()
  {
    for (int i = 0; i < UF_MAX; i++) {
      UF_par[i] = i;
      UF_rank[i] = 0;
      UF_cnt[i] = 1;
    }
  }

  // �؂̍������߂�
  int find(int x)
  {
    if (UF_par[x] == x) {
      return x;
    }
    else {
      return UF_par[x] = find(UF_par[x]);
    }
  }

  // x��y�̑�����W���𕹍�
  void unite(int x, int y)
  {
    x = find(x);
    y = find(y);
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
  bool same(int x, int y) { return find(x) == find(y); }
}  // namespace

namespace /* ���낢�� */
{
  const int INF = 1001001001;
  const int dx[4] = { -1, 0, 1, 0 };
  const int dy[4] = { 0, -1, 0, 1 };
  const char cc[4] = { 'U', 'L', 'D', 'R' };
}  // namespace

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

  // 0�ȏ�1�����̏������Ƃ闐��
  static double rand01() { return (Rand() + 0.5) * (1.0 / UINT_MAX); }
}  // namespace

namespace /* �ϐ� */
{
  int mode = 0;
  ofstream ofs;

  // ���͗p�ϐ�
  const int MIN_ROCK = 10;
  const int MAX_ROCK = 5000;
  const int CC[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  int N, W, K, C;
  int a[20], b[20], c[20], d[20];
  int S[210][210];

  int f[210][210];
  int minS[210][210], maxS[210][210];
  int nearX[20], nearY[20];

  // �ϐ�
  int hp;
  int attack_power[8] = { 50, 50, 50, 50, 50, 100, 100, 100 };
  int ATTACK_POWER = 100;

}  // namespace

inline bool IsNG(int x, int y)
{
  if (x < 0 || N <= x || y < 0 || N <= y) {
    return true;
  }
  return false;
}

inline bool IsUniteWater(int x, int y) { return same(x * N + y, UF_MAX - 1); }

inline int Manhattan(int x1, int y1, int x2, int y2)
{
  return abs(x1 - x2) + abs(y1 - y2);
}

// �X�R�A�v�Z
int CalcScore() { return hp; }

int Attack(int x, int y, int power)
{
  minS[x][y] = maxS[x][y];
  maxS[x][y] += power;
  hp += C + power;
  int res = 0;
  if (mode == 0) {
    cout << x << ' ' << y << ' ' << power << endl;
    fflush(stdout);
    cin >> res;
  }
  else {
    ofs << x << ' ' << y << ' ' << power << endl;
    S[x][y] -= power;
    if (S[x][y] <= 0) {
      res = 1;
    }
  }

  if (res != 0) {
    f[x][y] = 1;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (IsNG(nx, ny)) continue;
      if (f[nx][ny] != 0) {
        unite(x * N + y, nx * N + ny);
      }
    }
    rep(i, W)
    {
      if (x == a[i] && y == b[i]) {
        unite(x * N + y, UF_MAX - 1);
      }
    }
    res = 2;
    rep(i, K)
    {
      if (!IsUniteWater(c[i], d[i])) {
        res = 1;
      }
    }
  }

  return res;
}

int Challenge(int x, int y, int power)
{
  int res = f[x][y];
  while (res == 0) {
    res = Attack(x, y, power);
  }
  return res;
}

// ���[�J���ŕ����P�[�X�������߂̑S�ď����֐�
void AllClear_MultiCase() {}

// ������ԍ쐬�i������Ăׂ΃X�^�[�g�ʒu�ɖ߂�邱�Ƃ�z��Areal_maxScore���͖߂��Ȃ��j
void Init(int problemNum)
{
  UFinit();
  hp = 0;

  rep(i, N)
  {
    rep(j, N)
    {
      f[i][j] = 0;
      minS[i][j] = 0;
      maxS[i][j] = 0;
    }
  }

  rep(i, 8)
  {
    if (C == CC[i]) {
      ATTACK_POWER = attack_power[i];
    }
  }

  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
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
  if (mode == 0 || !ifs.is_open()) {
    cin >> N >> W >> K >> C;
    rep(i, W) { cin >> a[i] >> b[i]; }
    rep(i, K) { cin >> c[i] >> d[i]; }
  }
  // �t�@�C�����͂���
  else {
    ifs >> N >> W >> K >> C;
    rep(i, N)
    {
      rep(j, N) { ifs >> S[i][j]; }
    }
    rep(i, W) { ifs >> a[i] >> b[i]; }
    rep(i, K) { ifs >> c[i] >> d[i]; }
  }
}

// �𓚏o��
void Output(int problemNum)
{
  if (mode != 0) {
    ofs.close();
  }
}

int Solve(int problemNum = 0)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // ������ԍ쐬
  Init(problemNum);

  // �e�Ƃ̈�ԋ߂�������T��
  rep(i, K)
  {
    int dist = INF;
    rep(j, W)
    {
      if (Manhattan(c[i], d[i], a[j], b[j]) < dist) {
        dist = Manhattan(c[i], d[i], a[j], b[j]);
        nearX[i] = a[j];
        nearY[i] = b[j];
      }
    }
  }

  // �Ƃ𐅌��ɂȂ��Ă���
  rep(i, K)
  {
    int phase = 0;
    int nowX = c[i], nowY = d[i];
    int nextX = -1, nextY = -1;
    int dir = -1;
    while (!IsUniteWater(c[i], d[i])) {
      if (phase == 0) {
        // �ƂɌ���������
        Challenge(nowX, nowY, ATTACK_POWER);
        phase = 1;
      }
      else if (phase == 1) {
        // ���̃`�F�b�N�|�C���g�����߂�
        int diffX = nearX[i] - nowX;
        if (abs(diffX) > 20) {
          if (diffX > 0) {
            diffX = 20;
          }
          else {
            diffX = -20;
          }
        }
        int diffY = nearY[i] - nowY;
        if (abs(diffY) > 20) {
          if (diffY > 0) {
            diffY = 20;
          }
          else {
            diffY = -20;
          }
        }

        // �߂��ɐ��H������΂������Ɍ�����
        int diffSum = abs(diffX) + abs(diffY);
        srep(k, 1, diffSum)
        {
          int ok = 0;
          rep(j, 4)
          {
            int nx = nowX + dx[j] * k;
            int ny = nowY + dy[j] * k;
            if (!IsNG(nx, ny) && IsUniteWater(nx, ny)) {
              diffX = nx - nowX;
              diffY = ny - nowY;
              ok = 1;
              break;
            }
          }
          if (ok) break;
        }

        if (diffX == 0) {
          // Y�����ɐi��
          nextX = nowX;
          nextY = nowY + diffY;
          Challenge(nextX, nextY, ATTACK_POWER);
        }
        else if (diffY == 0) {
          // X�����ɐi��
          nextX = nowX + diffX;
          nextY = nowY;
          Challenge(nextX, nextY, ATTACK_POWER);
        }
        else {
          int nextX1 = nowX + diffX;
          int nextY1 = nowY;
          Challenge(nextX1, nextY1, ATTACK_POWER);
          int nextX2 = nowX;
          int nextY2 = nowY + diffY;
          Challenge(nextX2, nextY2, ATTACK_POWER);
          if (maxS[nextX1][nextY1] <= maxS[nextX2][nextY2]) {
            nextX = nextX1;
            nextY = nextY1;
          }
          else {
            nextX = nextX2;
            nextY = nextY2;
          }
        }

        if (nextX - nowX < 0) {
          dir = 0;
        }
        else if (nextX - nowX > 0) {
          dir = 2;
        }
        else if (nextY - nowY < 0) {
          dir = 1;
        }
        else if (nextY - nowY > 0) {
          dir = 3;
        }

        phase = 2;
      }
      else if (phase == 2) {
        if (nowX == nextX && nowY == nextY) {
          phase = 1;
        }
        else {
          // �`�F�b�N�|�C���g�ɐi��
          int nx = nowX + dx[dir];
          int ny = nowY + dy[dir];
          if (maxS[nx][ny] == 0) {
            int d1 = Manhattan(nowX, nowY, nx, ny);
            int d2 = Manhattan(nx, ny, nextX, nextY);
            int p1 = minS[nowX][nowY];
            int p2 = minS[nextX][nextY];
            int power = (d2 * p1 + d1 * p2) / (d1 + d2);
            power = max(power, 10);
            Attack(nx, ny, power);
          }
          else {
            Challenge(nx, ny, ATTACK_POWER);
            nowX = nx;
            nowY = ny;
          }
        }
      }
    }
  }

  // �f�o�b�O���O
  if (mode != 0) {
    cout << "problemNum = " << problemNum << ", hp = " << hp << endl;
  }

  return CalcScore();
}

int SolveOuter(int problemNum = 0)
{
  // ���͎󂯎��
  Input(problemNum);

  int score = Solve(problemNum);

  // �𓚂̏o��
  Output(problemNum);

  return score;
}

int main()
{
  // ��������
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

  mode = 0;

  // ��o�p
  if (mode == 0) {
    SolveOuter();
  }
  // 1�P�[�X����
  else if (mode == 1) {
    SolveOuter(0);
  }
  // �����P�[�X����
  else if (mode == 2) {
    ll scoreSum = 0;
    rep(i, 100)
    {
      scoreSum += SolveOuter(i);
      AllClear_MultiCase();
    }
    cout << "scoreSum = " << scoreSum << endl;
  }

  return 0;
}
