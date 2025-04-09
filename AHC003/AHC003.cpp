#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <climits>
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <string>
#include <cstring>
#include <ctime>
#include<queue>
// #include <atcoder/all>
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
// using namespace atcoder;
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<double, P> PDP;
#define yn {puts("Yes");}else{puts("No");}
#define MAX_N 200005
const int INF = 1001001001;

namespace /* ���� */
{
  static uint32_t Rand()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0�ȏ�1�����̏������Ƃ闐��
  static double rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}

namespace /* �O���b�h�p */
{
  inline int ManhattanDistance(int x1, int y1, int x2, int y2)
  {
    return std::abs(x1 - x2) + std::abs(y1 - y2);
  }

  const int dx[4] = { -1, 0, 1, 0 };
  const int dy[4] = { 0, -1, 0, 1 };
  const string ULDR = "ULDR";
}

namespace /* ���ݒ� */
{
  const int n = 30;
  const int Q = 1000;
}

namespace /* �n�C�p�[�p�����[�^ */
{
  double initialD = 3585;
  double start_temp = 103.48026;
  double end_temp = 6.74495e-07;
  int loopTimes = 3007;
  int addTimes = 5007;

  double SabunCostMultiple = 1.5;
}

namespace /* ���[�J���p */
{
  double dReal[2][n + 1][n + 1]; // �����̒���
  double scoreSumGlobal; // �X�R�A�Ǘ�
}

namespace /* �ϐ� */
{
  PDP p;
}

namespace /* �Ă��Ȃ܂��p�ϐ� */
{
  double DistResponce[Q] = {};				// �󂯎��������
}

namespace /* �s�A��A�؂�ڂ̍\�� */
{
  double Upper[n] = {}, Downer[n] = {}, Lefter[n] = {}, Righter[n] = {};
  int CutUD[n] = {}, CutLR[n] = {}; // 0~30���Ƃ锼�J���
  int VerticalSum[Q][n + 1][n + 1];
  int HorizontalSum[Q][n + 1][n + 1];
  vector<int> TurnUD[n], TurnLR[n];

  double BestUpper[n], BestDowner[n], BestLefter[n], BestRighter[n];
  int BestCutUD[n], BestCutLR[n];

  double DifferenceSum = 0;
  double DistEstimate[Q] = {};

  vector<pair<double, double>> vec(1100);

  double dUD[n + 1][n + 1], dLR[n + 1][n + 1];
  vector<int> PathIDVectorUD[n][n], PathIDVectorLR[n][n];
}

void ClearGlobalVariables()
{
  scoreSumGlobal = 0;
  rep(i, 2)rep(j, n + 1)rep(k, n + 1) dReal[i][j][k] = 0;
  rep(i, Q)
  {
    rep(j, n + 1)
    {
      rep(k, n + 1)
      {
        VerticalSum[i][j][k] = 0;
        HorizontalSum[i][j][k] = 0;
      }
    }
  }
  rep(i, n)
  {
    TurnUD[i].clear();
    TurnLR[i].clear();
  }
  rep(i, n)
  {
    rep(j, n)
    {
      PathIDVectorLR[i][j].clear();
      PathIDVectorUD[i][j].clear();
      dUD[i][j] = 0;
      dLR[i][j] = 0;
    }
  }
}

bool IsNgNxNy(int nx, int ny)
{
  if (nx < 0 || n <= nx || ny < 0 || n <= ny) return 1;
  return 0;
}

int CalcScore(const int sx, const int sy, const string& ans, const double aValue, const double eValue, int k)
{
  int x = sx, y = sy;
  int res = 0;
  rep(i, ans.size())
  {
    if (ans[i] == 'U') {
      res += dReal[0][x][y];
      x--;
    }
    if (ans[i] == 'L') {
      res += dReal[1][x][y];
      y--;
    }
    if (ans[i] == 'D') {
      res += dReal[0][x + 1][y];
      x++;
    }
    if (ans[i] == 'R') {
      res += dReal[1][x][y + 1];
      y++;
    }
  }
  scoreSumGlobal += pow(0.998, 1000 - k) * aValue / res;
  return res * eValue;
}

double dp[n][n];
int nxt[n][n];

void Dijkstra(int sx, int sy, int gx, int gy)
{
  rep(i, n) rep(j, n) dp[i][j] = INF;
  dp[sx][sy] = 0;
  priority_queue<PDP, vector<PDP>, greater<PDP>> que;
  p.first = 0;
  p.second.first = sx;
  p.second.second = sy;
  que.push(p);
  while (que.size()) {
    p = que.top();
    que.pop();
    int x = p.second.first;
    int y = p.second.second;
    double val = p.first;
    if (x == gx && y == gy) break;
    if (val > dp[x][y]) continue;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (IsNgNxNy(nx, ny)) continue;
      double dd = 0;
      if (i == 0) {
        if (x <= CutUD[y]) dd = Upper[y];
        else dd = Downer[y];
      }
      if (i == 1) {
        if (y <= CutLR[x]) dd = Lefter[x];
        else dd = Righter[x];
      }
      if (i == 2) {
        if (x + 1 <= CutUD[y]) dd = Upper[y];
        else dd = Downer[y];
      }
      if (i == 3) {
        if (y + 1 <= CutLR[x]) dd = Lefter[x];
        else dd = Righter[x];
      }
      if (dp[nx][ny] > dp[x][y] + dd) {
        dp[nx][ny] = dp[x][y] + dd;
        nxt[nx][ny] = i;
        p.first = dp[nx][ny];
        p.second.first = nx;
        p.second.second = ny;
        que.push(p);
      }
    }
  }
}

void Dijkstra2(int sx, int sy, int gx, int gy)
{
  rep(i, n) rep(j, n) dp[i][j] = INF;
  dp[sx][sy] = 0;
  priority_queue<PDP, vector<PDP>, greater<PDP>> que;
  p.first = 0;
  p.second.first = sx;
  p.second.second = sy;
  que.push(p);
  while (que.size()) {
    p = que.top();
    que.pop();
    int x = p.second.first;
    int y = p.second.second;
    double val = p.first;
    if (x == gx && y == gy) break;
    if (val > dp[x][y]) continue;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (IsNgNxNy(nx, ny)) continue;
      double dd = 0;
      if (i == 0) {
        if (x <= CutUD[y]) dd = Upper[y];
        else dd = Downer[y];
        dd += dUD[x][y];
      }
      if (i == 1) {
        if (y <= CutLR[x]) dd = Lefter[x];
        else dd = Righter[x];
        dd += dLR[x][y];
      }
      if (i == 2) {
        if (x + 1 <= CutUD[y]) dd = Upper[y];
        else dd = Downer[y];
        dd += dUD[x + 1][y];
      }
      if (i == 3) {
        if (y + 1 <= CutLR[x]) dd = Lefter[x];
        else dd = Righter[x];
        dd += dLR[x][y + 1];
      }

      dd = max(dd, 1.0);

      if (dp[nx][ny] > dp[x][y] + dd) {
        dp[nx][ny] = dp[x][y] + dd;
        nxt[nx][ny] = i;
        p.first = dp[nx][ny];
        p.second.first = nx;
        p.second.second = ny;
        que.push(p);
      }
    }
  }
}


/* ����
  �ӂ̔ԍ���1~29(1-index)
  �c��(i-1,j)->(i,j)��d[0][i][j];
  ����(i,j-1)->(i,j)��d[1][i][j];
*/

int Solve(string iunputFileNum)
{
  // ���Ԍv��
  clock_t start, end;
  start = clock();


  int inputMode = 1;
  string fileNameIfs = (string)"./in/" + iunputFileNum + ".txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs) { // �W�����͂���
    inputMode = 0;
  }
  else {
    rep(i, n) srep(j, 1, n) ifs >> dReal[1][i][j];
    srep(i, 1, n) rep(j, n) ifs >> dReal[0][i][j];
  }

  string fileName = (string)"./out/" + iunputFileNum + "_out.txt";
  const char* cstr = fileName.c_str();
  ofstream ofs(cstr);

  rep(i, n)
  {
    Upper[i] = initialD;
    Downer[i] = initialD;
    Lefter[i] = initialD;
    Righter[i] = initialD;
    CutUD[i] = 15;
    CutLR[i] = 15;
  }
  int timeOverFlag = 0;
  int ALLloop = 0;
  rep(turn, Q)
  {
    int sx, sy, gx, gy;

    // ���͎󂯎��
    if (inputMode == 0) {
      cin >> sx >> sy >> gx >> gy;
    }
    else {
      ifs >> sx >> sy >> gx >> gy;
    }

    // �L�[�v
    int keepSx = sx, keepSy = sy, keepGx = gx, keepGy = gy;

    // �v���C�I���e�B�[�L���[�ōŒZ�H�����߂�i�_�C�N�X�g���j
    Dijkstra2(sx, sy, gx, gy);

    // �񓚕����񐶐�
    string ans;
    vector<int> v;
    while (gx != sx || gy != sy) {
      v.push_back(nxt[gx][gy]);
      ans += ULDR[nxt[gx][gy]];
      int ggx = gx - dx[nxt[gx][gy]];
      int ggy = gy - dy[nxt[gx][gy]];
      gx = ggx;
      gy = ggy;
    }
    gx = keepGx; gy = keepGy;
    reverse(ans.begin(), ans.end());
    reverse(v.begin(), v.end());
    if (inputMode == 0) {
      cout << ans << endl;
      fflush(stdout);
    }
    else {
      ofs << ans << endl;
    }

    // �p�X�̒����̎󂯎��
    double dist = 0;
    if (inputMode == 0) {
      cin >> dist;
    }
    else {
      double aValue, eValue;
      ifs >> aValue >> eValue;
      dist = CalcScore(sx, sy, ans, aValue, eValue, turn + 1);
    }

    // �R�o��p
    DistResponce[turn] = dist;

    int length = ans.size();
    double myDistance = dp[gx][gy];

    if (turn < 0) continue;
    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 1.8) {
      if (timeOverFlag == 0 && ifs) {
        cout << "TIME OVER" << ' ' << turn << endl;
        timeOverFlag = 1;
      }
      continue;
    }

    // �V�����f�[�^�\���X�V
    int tatePath[n + 1][n + 1];
    int yokoPath[n + 1][n + 1];
    rep(i, n + 1)
    {
      rep(j, n + 1)
      {
        tatePath[i][j] = 0;
        yokoPath[i][j] = 0;
      }
    }
    rep(i, length)
    {
      int nsx = sx;
      int nsy = sy;
      if (v[i] == 2) nsx += 1;
      if (v[i] == 3) nsy += 1;

      if (v[i] % 2 == 0) {
        tatePath[nsx][nsy] = 1;
        PathIDVectorUD[nsx][nsy].push_back(turn);
      }
      else {
        yokoPath[nsx][nsy] = 1;
        PathIDVectorLR[nsx][nsy].push_back(turn);
      }

      sx += dx[v[i]];
      sy += dy[v[i]];
    }
    sx = keepSx; sy = keepSy;
    rep(j, n)
    {
      srep(i, 1, n)
      {
        VerticalSum[turn][i][j] = VerticalSum[turn][i - 1][j] + tatePath[i][j];
      }
      if (VerticalSum[turn][n - 1][j]) TurnUD[j].push_back(turn);
    }
    rep(i, n)
    {
      srep(j, 1, n)
      {
        HorizontalSum[turn][i][j] = HorizontalSum[turn][i][j - 1] + yokoPath[i][j];
      }
      if (HorizontalSum[turn][i][n - 1]) TurnLR[i].push_back(turn);
    }

    DifferenceSum = 0;
    rep(i, Q) DistEstimate[i] = 0;
    rep(j, n)
    {
      rep(i, TurnUD[j].size())
      {
        int turnID = TurnUD[j][i];
        DistEstimate[turnID] += Upper[j] * VerticalSum[turnID][CutUD[j]][j] + Downer[j] * (VerticalSum[turnID][n - 1][j] - VerticalSum[turnID][CutUD[j]][j]);
      }
    }
    rep(j, n)
    {
      rep(i, TurnLR[j].size())
      {
        int turnID = TurnLR[j][i];
        DistEstimate[turnID] += Lefter[j] * HorizontalSum[turnID][j][CutLR[j]] + Righter[j] * (HorizontalSum[turnID][j][n - 1] - HorizontalSum[turnID][j][CutLR[j]]);
      }
    }
    rep(i, turn + 1)
    {
      DifferenceSum += std::abs(DistResponce[i] - DistEstimate[i]) * (40000.0 / DistResponce[i]);
    }

    rep(i, n)
    {
      DifferenceSum += std::abs(Upper[i] - Downer[i]) * SabunCostMultiple;
      DifferenceSum += std::abs(Lefter[i] - Righter[i]) * SabunCostMultiple;
    }

    double BestDifferenceSum = DifferenceSum;
    rep(i, n)
    {
      BestUpper[i] = Upper[i];
      BestDowner[i] = Downer[i];
      BestLefter[i] = Lefter[i];
      BestRighter[i] = Righter[i];
      BestCutUD[i] = CutUD[i];
      BestCutLR[i] = CutLR[i];
    }

    // �Ă��Ȃ܂��ŗ\���l������ɒ���
    end = clock();
    double now_time = (double)(end - start) / CLOCKS_PER_SEC;
    double TL = 1.8;

    rep(loop, loopTimes)
    {
      ALLloop++;
      if (loop % 100 == 1) {
        end = clock();
        now_time = (double)(end - start) / CLOCKS_PER_SEC;
      }
      if (now_time > TL) {
        break;
      }

      double temp = start_temp + (end_temp - start_temp) * ((double)loop / (loopTimes + addTimes));


      int MODE = loop % 3;
      if (turn >= 200 && ALLloop % 5247 == 1853) MODE = 3;

      // Upper, Downer, Lefter, Righter�̂����ꂩ���ύX
      if (MODE == 0) {
        int rULDR = Rand() % 4;

        int rN = Rand() % n;

        double rA = rand01() * 20 - 10;

        if (rULDR == 0) if (Upper[rN] + rA < 1000 || 9000 < Upper[rN] + rA) continue;
        if (rULDR == 1) if (Lefter[rN] + rA < 1000 || 9000 < Lefter[rN] + rA) continue;
        if (rULDR == 2) if (Downer[rN] + rA < 1000 || 9000 < Downer[rN] + rA) continue;
        if (rULDR == 3) if (Righter[rN] + rA < 1000 || 9000 < Righter[rN] + rA) continue;

        double diff = 0;

        if (rULDR == 0) {
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = rA * VerticalSum[turnID][CutUD[rN]][rN];
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
          diff += std::abs(Upper[rN] - Downer[rN]) * SabunCostMultiple - std::abs((Upper[rN] + rA) - Downer[rN]) * SabunCostMultiple;
        }
        if (rULDR == 1) {
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = rA * HorizontalSum[turnID][rN][CutLR[rN]];
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
          diff += std::abs(Lefter[rN] - Righter[rN]) * SabunCostMultiple - std::abs((Lefter[rN] + rA) - Righter[rN]) * SabunCostMultiple;
        }
        if (rULDR == 2) {
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = rA * (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
          diff += std::abs(Upper[rN] - Downer[rN]) * SabunCostMultiple - std::abs(Upper[rN] - (Downer[rN] + rA)) * SabunCostMultiple;
        }
        if (rULDR == 3) {
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = rA * (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][CutLR[rN]]);
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
          diff += std::abs(Lefter[rN] - Righter[rN]) * SabunCostMultiple - std::abs(Lefter[rN] - (Righter[rN] + rA)) * SabunCostMultiple;
        }



        double prob = exp((double)diff / temp);

        if (prob > rand01()) {
          if (rULDR == 0) {
            rep(i, TurnUD[rN].size())
            {
              int turnID = TurnUD[rN][i];
              double tmpDiff = rA * VerticalSum[turnID][CutUD[rN]][rN];
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rULDR == 1) {
            rep(i, TurnLR[rN].size())
            {
              int turnID = TurnLR[rN][i];
              double tmpDiff = rA * HorizontalSum[turnID][rN][CutLR[rN]];
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rULDR == 2) {
            rep(i, TurnUD[rN].size())
            {
              int turnID = TurnUD[rN][i];
              double tmpDiff = rA * (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rULDR == 3) {
            rep(i, TurnLR[rN].size())
            {
              int turnID = TurnLR[rN][i];
              double tmpDiff = rA * (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][CutLR[rN]]);
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rULDR == 0) Upper[rN] += rA;
          if (rULDR == 1) Lefter[rN] += rA;
          if (rULDR == 2) Downer[rN] += rA;
          if (rULDR == 3) Righter[rN] += rA;
          DifferenceSum -= diff;
        }
      }

      // Upper��Downer�A��������Lefter��Righter��ύX
      if (MODE == 1) {
        int rUL = Rand() % 2;

        int rN = Rand() % n;

        double rA = rand01() * 20 - 10;

        if (rUL == 0) if (Upper[rN] + rA < 1000 || 9000 < Upper[rN] + rA) continue;
        if (rUL == 1) if (Lefter[rN] + rA < 1000 || 9000 < Lefter[rN] + rA) continue;
        if (rUL == 0) if (Downer[rN] + rA < 1000 || 9000 < Downer[rN] + rA) continue;
        if (rUL == 1) if (Righter[rN] + rA < 1000 || 9000 < Righter[rN] + rA) continue;

        double diff = 0;

        if (rUL == 0) {
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = rA * VerticalSum[turnID][n - 1][rN];
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
        }
        if (rUL == 1) {
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = rA * HorizontalSum[turnID][rN][n - 1];
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
        }



        double prob = exp((double)diff / temp);

        if (prob > rand01()) {
          if (rUL == 0) {
            rep(i, TurnUD[rN].size())
            {
              int turnID = TurnUD[rN][i];
              double tmpDiff = rA * VerticalSum[turnID][n - 1][rN];
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rUL == 1) {
            rep(i, TurnLR[rN].size())
            {
              int turnID = TurnLR[rN][i];
              double tmpDiff = rA * HorizontalSum[turnID][rN][n - 1];
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rUL == 0) Upper[rN] += rA;
          if (rUL == 1) Lefter[rN] += rA;
          if (rUL == 0) Downer[rN] += rA;
          if (rUL == 1) Righter[rN] += rA;
          DifferenceSum -= diff;
        }
      }

      // ��؂���ꂩ���ύX
      if (MODE == 2) {
        int rUL = Rand() % 2;
        int rN = Rand() % 30;
        int rA = Rand() % 30 + 1;
        if (Rand() % 2 == 1) rA *= -1;
        if (rUL == 0) if (CutUD[rN] + rA < 3 || 27 <= CutUD[rN] + rA) continue;
        if (rUL == 1) if (CutLR[rN] + rA < 3 || 27 <= CutLR[rN] + rA) continue;

        double diff = 0;

        if (rUL == 0) {
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = 0;
            if (rA > 0) {
              tmpDiff = (Upper[rN] - Downer[rN]) * (VerticalSum[turnID][CutUD[rN] + rA][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
            }
            else {
              tmpDiff = (Downer[rN] - Upper[rN]) * (VerticalSum[turnID][CutUD[rN]][rN] - VerticalSum[turnID][CutUD[rN] + rA][rN]);
            }
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
        }
        if (rUL == 1) {
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = 0;
            if (rA > 0) {
              tmpDiff = (Lefter[rN] - Righter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN] + rA] - HorizontalSum[turnID][rN][CutLR[rN]]);
            }
            else {
              tmpDiff = (Righter[rN] - Lefter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN]] - HorizontalSum[turnID][rN][CutLR[rN] + rA]);
            }
            diff += (std::abs(DistResponce[turnID] - DistEstimate[turnID]) - std::abs(DistResponce[turnID] - (DistEstimate[turnID] + tmpDiff))) * (40000.0 / DistResponce[turnID]);
          }
        }




        double prob = exp((double)diff / temp);

        if (prob > rand01()) {
          if (rUL == 0) {
            rep(i, TurnUD[rN].size())
            {
              int turnID = TurnUD[rN][i];
              double tmpDiff = 0;
              if (rA > 0) {
                tmpDiff = (Upper[rN] - Downer[rN]) * (VerticalSum[turnID][CutUD[rN] + rA][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
              }
              else {
                tmpDiff = (Downer[rN] - Upper[rN]) * (VerticalSum[turnID][CutUD[rN]][rN] - VerticalSum[turnID][CutUD[rN] + rA][rN]);
              }
              DistEstimate[turnID] += tmpDiff;
            }
          }
          if (rUL == 1) {
            rep(i, TurnLR[rN].size())
            {
              int turnID = TurnLR[rN][i];
              double tmpDiff = 0;
              if (rA > 0) {
                tmpDiff = (Lefter[rN] - Righter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN] + rA] - HorizontalSum[turnID][rN][CutLR[rN]]);
              }
              else {
                tmpDiff = (Righter[rN] - Lefter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN]] - HorizontalSum[turnID][rN][CutLR[rN] + rA]);
              }
              DistEstimate[turnID] += tmpDiff;
            }
          }

          if (rUL == 0) CutUD[rN] += rA;
          if (rUL == 1) CutLR[rN] += rA;

          DifferenceSum -= diff;
        }
      }

      // 1���C������
      if (MODE == 3) {
        int rUL = Rand() % 2;
        int rN = Rand() % n;
        if (rUL == 0) {
          vector<double> amari;
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            amari.push_back(DistResponce[turnID] - (DistEstimate[turnID] - (Upper[rN] * VerticalSum[turnID][CutUD[rN]][rN] + Downer[rN] * (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][CutUD[rN]][rN]))));
          }
          double argmaxUpper = Upper[rN], argmaxDowner = Downer[rN];
          int argmaxCutUD = CutUD[rN];
          double maxDiff = 0;
          int changeFlag = 0;
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            maxDiff += (abs(amari[i]) - std::abs(DistResponce[turnID] - DistEstimate[turnID])) * (40000.0 / DistResponce[turnID]);
          }
          double keepMaxDiff = maxDiff;
          srep(kugiri, 3, 27)
          {
            rep(randomChallenge, 20)
            {
              double rA = Rand() % 8001 + 1000;

              double countSum = 0;
              rep(i, TurnUD[rN].size())
              {
                int turnID = TurnUD[rN][i];
                if (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][kugiri][rN] == 0) vec[i].first = 1001001;
                else vec[i].first = (amari[i] - rA * VerticalSum[turnID][kugiri][rN]) / (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][kugiri][rN]);
                vec[i].second = (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][kugiri][rN]) * (40000.0 / DistResponce[turnID]);
                countSum += vec[i].second;
              }
              if (countSum == 0) continue;
              vec[TurnUD[rN].size()].first = rA;
              vec[TurnUD[rN].size()].second = SabunCostMultiple;
              countSum += SabunCostMultiple;
              sort(vec.begin(), vec.begin() + TurnUD[rN].size() + 1);
              int ite = 0;
              double countNow = vec[ite].second;
              while (ite < TurnUD[rN].size() + 1) {
                if (countNow < countSum - countNow) {
                  ite++;
                  countNow += vec[ite].second;
                }
                else {
                  break;
                }
              }
              double tmpDiff = 0;
              double DownValue = vec[ite].first;
              DownValue = max(DownValue, 1000.0);
              DownValue = min(DownValue, 9000.0);
              rep(i, TurnUD[rN].size())
              {
                int turnID = TurnUD[rN][i];
                tmpDiff += (abs(amari[i]) - std::abs(amari[i] - (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][kugiri][rN]) * DownValue)) * (40000.0 / DistResponce[turnID]);
              }
              tmpDiff += std::abs(Upper[rN] - Downer[rN]) * SabunCostMultiple - std::abs(rA - DownValue) * SabunCostMultiple;

              if (tmpDiff > maxDiff) {
                maxDiff = tmpDiff;
                argmaxUpper = rA;
                argmaxDowner = DownValue;
                argmaxCutUD = kugiri;
                changeFlag = 1;
              }
            }
          }
          if (changeFlag == 0) continue;
          int kugiriDiff = argmaxCutUD - CutUD[rN];
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = 0;
            if (kugiriDiff > 0) {
              tmpDiff = (Upper[rN] - Downer[rN]) * (VerticalSum[turnID][CutUD[rN] + kugiriDiff][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
            }
            else {
              tmpDiff = (Downer[rN] - Upper[rN]) * (VerticalSum[turnID][CutUD[rN]][rN] - VerticalSum[turnID][CutUD[rN] + kugiriDiff][rN]);
            }
            DistEstimate[turnID] += tmpDiff;
          }
          CutUD[rN] = argmaxCutUD;
          double UpperDiff = argmaxUpper - Upper[rN];
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = UpperDiff * VerticalSum[turnID][CutUD[rN]][rN];
            DistEstimate[turnID] += tmpDiff;
          }
          Upper[rN] = argmaxUpper;
          double DownerDiff = argmaxDowner - Downer[rN];
          rep(i, TurnUD[rN].size())
          {
            int turnID = TurnUD[rN][i];
            double tmpDiff = DownerDiff * (VerticalSum[turnID][n - 1][rN] - VerticalSum[turnID][CutUD[rN]][rN]);
            DistEstimate[turnID] += tmpDiff;
          }
          Downer[rN] = argmaxDowner;
          DifferenceSum -= (maxDiff - keepMaxDiff);
        }
        else {
          vector<double> amari;
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            amari.push_back(DistResponce[turnID] - (DistEstimate[turnID] - (Lefter[rN] * HorizontalSum[turnID][rN][CutLR[rN]] + Righter[rN] * (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][CutLR[rN]]))));
          }
          double argmaxLefter = Lefter[rN], argmaxRighter = Righter[rN];
          int argmaxCutLR = CutLR[rN];
          double maxDiff = 0;
          int changeFlag = 0;
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            maxDiff += (abs(amari[i]) - std::abs(DistResponce[turnID] - DistEstimate[turnID])) * (40000.0 / DistResponce[turnID]);
          }
          double keepMaxDiff = maxDiff;
          srep(kugiri, 3, 27)
          {
            rep(randomChallenge, 20)
            {
              double rA = Rand() % 8001 + 1000;
              double countSum = 0;
              rep(i, TurnLR[rN].size())
              {
                int turnID = TurnLR[rN][i];
                if (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][kugiri] == 0) vec[i].first = 1001001;
                else vec[i].first = (amari[i] - rA * HorizontalSum[turnID][rN][kugiri]) / (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][kugiri]);
                vec[i].second = (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][kugiri]) * (40000.0 / DistResponce[turnID]);
                countSum += vec[i].second;
              }
              if (countSum == 0) continue;
              vec[TurnLR[rN].size()].first = rA;
              vec[TurnLR[rN].size()].second = SabunCostMultiple;
              countSum += SabunCostMultiple;
              sort(vec.begin(), vec.begin() + TurnLR[rN].size() + 1);
              int ite = 0;
              double countNow = vec[ite].second;
              while (ite < TurnLR[rN].size() + 1) {
                if (countNow < countSum - countNow) {
                  ite++;
                  countNow += vec[ite].second;
                }
                else {
                  break;
                }
              }
              double tmpDiff = 0;
              double RightValue = vec[ite].first;
              RightValue = max(RightValue, 1000.0);
              RightValue = min(RightValue, 9000.0);
              rep(i, TurnLR[rN].size())
              {
                int turnID = TurnLR[rN][i];
                tmpDiff += (abs(amari[i]) - std::abs(amari[i] - (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][kugiri]) * RightValue)) * (40000.0 / DistResponce[turnID]);
              }
              tmpDiff += std::abs(Lefter[rN] - Righter[rN]) * SabunCostMultiple - std::abs(rA - RightValue) * SabunCostMultiple;

              if (tmpDiff > maxDiff) {
                maxDiff = tmpDiff;
                argmaxLefter = rA;
                argmaxRighter = RightValue;
                argmaxCutLR = kugiri;
                changeFlag = 1;
              }
            }
          }
          if (changeFlag == 0) continue;
          int kugiriDiff = argmaxCutLR - CutLR[rN];
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = 0;
            if (kugiriDiff > 0) {
              tmpDiff = (Lefter[rN] - Righter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN] + kugiriDiff] - HorizontalSum[turnID][rN][CutLR[rN]]);
            }
            else {
              tmpDiff = (Righter[rN] - Lefter[rN]) * (HorizontalSum[turnID][rN][CutLR[rN]] - HorizontalSum[turnID][rN][CutLR[rN] + kugiriDiff]);
            }
            DistEstimate[turnID] += tmpDiff;
          }
          CutLR[rN] = argmaxCutLR;
          double LefterDiff = argmaxLefter - Lefter[rN];
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = LefterDiff * HorizontalSum[turnID][rN][CutLR[rN]];
            DistEstimate[turnID] += tmpDiff;
          }
          Lefter[rN] = argmaxLefter;
          double RighterDiff = argmaxRighter - Righter[rN];
          rep(i, TurnLR[rN].size())
          {
            int turnID = TurnLR[rN][i];
            double tmpDiff = RighterDiff * (HorizontalSum[turnID][rN][n - 1] - HorizontalSum[turnID][rN][CutLR[rN]]);
            DistEstimate[turnID] += tmpDiff;
          }
          Righter[rN] = argmaxRighter;
          DifferenceSum -= (maxDiff - keepMaxDiff);
        }
      }
    }


    // �X�R�A�����������烍�[���o�b�N

    if (BestDifferenceSum < DifferenceSum) {

      DifferenceSum = BestDifferenceSum;
      rep(i, n)
      {
        Upper[i] = BestUpper[i];
        Downer[i] = BestDowner[i];
        Lefter[i] = BestLefter[i];
        Righter[i] = BestRighter[i];
        CutUD[i] = BestCutUD[i];
        CutLR[i] = BestCutLR[i];
      }

    }

    rep(i, n)
    {
      rep(j, n)
      {
        dUD[i][j] = 0;
        dLR[i][j] = 0;
      }
    }

    rep(_, 20000)
    {
      int rv = Rand() % 2;
      int rx = 0, ry = 0;
      if (rv == 0) {
        rx = Rand() % 29 + 1;
        ry = Rand() % 30;
      }
      else {
        rx = Rand() % 30;
        ry = Rand() % 29 + 1;
      }
      double ra = rand01() * 200.0 - 100.0;
      if (rv == 0) {
        if (rx <= CutUD[ry]) {
          if (Upper[ry] + ra < 1000 || 9000 < Upper[ry] + ra) continue;
        }
        else {
          if (Downer[ry] + ra < 1000 || 9000 < Downer[ry] + ra) continue;
        }
      }
      else {
        if (ry <= CutLR[rx]) {
          if (Lefter[rx] + ra < 1000 || 9000 < Lefter[rx] + ra) continue;
        }
        else {
          if (Righter[rx] + ra < 1000 || 9000 < Righter[rx] + ra) continue;
        }
      }

      double diff = 0;

      if (rv == 0) {
        rep(i, PathIDVectorUD[rx][ry].size())
        {
          int pathID = PathIDVectorUD[rx][ry][i];
          diff += (std::abs(DistResponce[pathID] - DistEstimate[pathID]) - std::abs(DistResponce[pathID] - (DistEstimate[pathID] + ra))) * (40000.0 / DistResponce[pathID]);
        }
      }
      else {
        rep(i, PathIDVectorLR[rx][ry].size())
        {
          int pathID = PathIDVectorLR[rx][ry][i];
          diff += (std::abs(DistResponce[pathID] - DistEstimate[pathID]) - std::abs(DistResponce[pathID] - (DistEstimate[pathID] + ra))) * (40000.0 / DistResponce[pathID]);
        }
      }

      if (diff > 0) {
        if (rv == 0) {
          rep(i, PathIDVectorUD[rx][ry].size())
          {
            int pathID = PathIDVectorUD[rx][ry][i];
            DistEstimate[pathID] += ra;
          }
        }
        else {
          rep(i, PathIDVectorLR[rx][ry].size())
          {
            int pathID = PathIDVectorLR[rx][ry][i];
            DistEstimate[pathID] += ra;
          }
        }
        if (rv == 0) dUD[rx][ry] += ra;
        else dLR[rx][ry] += ra;
        DifferenceSum -= diff;
      }
    }

  }

  scoreSumGlobal *= 2312311;

  if (inputMode == 1) {
    string fileNameParam = (string)"./parameter/" + iunputFileNum + "_param.txt";
    const char* cstrParam = fileNameParam.c_str();
    ofstream ofsParam(cstrParam);
    ofsParam << "Upper / Downer" << endl;
    rep(i, n) ofsParam << i << ' ' << CutUD[i] << ' ' << Upper[i] << ' ' << Downer[i] << endl;
    ofsParam << "Lefter / Righter" << endl;
    rep(i, n) ofsParam << i << ' ' << CutLR[i] << ' ' << Lefter[i] << ' ' << Righter[i] << endl;

    ofsParam << "dUD" << endl;
    rep(i, n)
    {
      rep(j, n)
      {
        ofsParam << dUD[i][j] << ' ';
      }
      ofsParam << endl;
    }
    ofsParam << "dLR" << endl;
    rep(i, n)
    {
      rep(j, n)
      {
        ofsParam << dLR[i][j] << ' ';
      }
      ofsParam << endl;
    }

    ofsParam.close();
  }

  ofs.close();

  return scoreSumGlobal;
}


int main()
{
  clock_t start, end;
  start = clock();

  srand((unsigned)time(NULL));
  while (rand() % 128) Rand();

  int mode = 0;

  if (mode == 0) { // ��o�p
    Solve("noinput");
  }
  else if (mode == 1) { // �T���v��0~99�Ń`�F�b�N
    vector<P> ranking;
    ll allScore = 0;
    srep(i, 0, 100)
    {
      string inputFileNum;
      inputFileNum += (char)((i % 10000) / 1000 + '0');
      inputFileNum += (char)((i % 1000) / 100 + '0');
      inputFileNum += (char)((i % 100) / 10 + '0');
      inputFileNum += (char)((i % 10) / 1 + '0');
      Solve(inputFileNum);
      cout << (int)scoreSumGlobal << " " << i << endl;
      ranking.push_back(P(scoreSumGlobal, i));
      allScore += scoreSumGlobal;
      ClearGlobalVariables();
    }
    sort(ranking.begin(), ranking.end());
    rep(i, 10) cout << ranking[i].first << ' ' << ranking[i].second << endl;
    cout << allScore << endl;
  }

  end = clock();
  if (mode != 0) cout << (double)(end - start) / CLOCKS_PER_SEC << endl;
}

