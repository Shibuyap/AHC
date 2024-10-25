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

namespace /* 乱数ライブラリ */
{
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

  // 0以上1未満の小数をとる乱数
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const int INF = 1001001001;
double TL = 1.8;
int mode;

int n, m;
int ans[5100][2];
int ansSize = 0;
int max_ans[5100][2];
int max_ansSize;
int maxScore;
int shoki_b[20][210];
int shoki_c[210][2];
int shoki_bCnt[20];
int b[20][210];
int c[210][2];
int bCnt[20];
int tmp_b[20][210];
int tmp_c[210][2];
int tmp_bCnt[20];
int gomi1;
int gomi2;
int gomi3;

clock_t s_time, e_time;

void CopyAll()
{
  rep(i, m)
  {
    rep(j, n / m) { b[i][j] = shoki_b[i][j]; }
  }
  rep(i, n)
  {
    rep(j, 2) { c[i][j] = shoki_c[i][j]; }
  }
  rep(i, m) { bCnt[i] = shoki_bCnt[i]; }
}

void CopyAns()
{
  rep(i, ansSize)
  {
    rep(j, 2) { max_ans[i][j] = ans[i][j]; }
  }
  max_ansSize = ansSize;
}

void RollbackAns()
{
  rep(i, max_ansSize)
  {
    rep(j, 2) { ans[i][j] = max_ans[i][j]; }
  }
  ansSize = max_ansSize;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  maxScore = -INF;
  ansSize = 0;
  gomi1 = -1;
  gomi2 = -1;
  gomi3 = -1;
}

void Clean()
{
  CopyAll();
  ansSize = 0;
  gomi1 = -1;
  gomi2 = -1;
  gomi3 = -1;
}

// スコア計算
int CalcScore()
{
  int cnt = 0;
  rep(i, m)
  {
    rep(j, n / m) { tmp_b[i][j] = shoki_b[i][j]; }
  }
  rep(i, n)
  {
    rep(j, 2) { tmp_c[i][j] = shoki_c[i][j]; }
  }
  rep(i, m) { tmp_bCnt[i] = shoki_bCnt[i]; }
  int res = 10000;
  rep(i, ansSize)
  {
    int ii = ans[i][0];
    int x = tmp_c[ii][0];
    int z = tmp_c[ii][1];
    int y = ans[i][1];
    if (ans[i][1] == -1) {
      tmp_bCnt[x]--;
      cnt++;
    }
    else {
      res--;
      rep(j, tmp_bCnt[x] - z)
      {
        int jj = tmp_b[x][z + j];
        tmp_c[jj][0] = y;
        tmp_c[jj][1] = tmp_bCnt[y];
        tmp_b[y][tmp_bCnt[y]] = tmp_b[x][z + j];
        tmp_bCnt[y]++;
        res--;
      }
      tmp_bCnt[x] = z;
    }
  }
  if (cnt != n) return -1;
  return res;
}

void UpdateAns()
{
  int score = CalcScore();
  if (score > maxScore) {
    maxScore = score;
    CopyAns();
  }
}

// 入力受け取り
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

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> n >> m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> b[i][j];
        b[i][j]--;
        c[b[i][j]][0] = i;
        c[b[i][j]][1] = j;
      }
    }
  }
  // ファイル入力する
  else {
    ifs >> n >> m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        ifs >> b[i][j];
        b[i][j]--;
        c[b[i][j]][0] = i;
        c[b[i][j]][1] = j;
      }
    }
  }

  rep(i, m)
  {
    rep(j, n / m) { shoki_b[i][j] = b[i][j]; }
  }
  rep(i, n)
  {
    rep(j, 2) { shoki_c[i][j] = c[i][j]; }
  }
  rep(i, m)
  {
    bCnt[i] = n / m;
    shoki_bCnt[i] = bCnt[i];
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

void Move(int x, int y, int z)
{
  rep(j, bCnt[x] - z)
  {
    int jj = b[x][z + j];
    b[y][bCnt[y]] = b[x][z + j];
    c[jj][0] = y;
    c[jj][1] = bCnt[y];
    bCnt[y]++;
  }
  bCnt[x] = z;
}

// 初期解生成
void Initialize()
{
  Clean();
  rep(i, n)
  {
    int x = c[i][0];
    int z = c[i][1];
    int y = randxor() % m;
    while (y == x) {
      y = randxor() % m;
    }
    if (z < bCnt[x] - 1) {
      int ii = b[x][z + 1];
      ans[ansSize][0] = ii;
      ans[ansSize][1] = y;
      ansSize++;
      rep(j, bCnt[x] - (z + 1))
      {
        int jj = b[x][z + 1 + j];
        b[y][bCnt[y]] = b[x][z + 1 + j];
        c[jj][0] = y;
        c[jj][1] = bCnt[y];
        bCnt[y]++;
      }
      bCnt[x] = z + 1;
    }
    ans[ansSize][0] = i;
    ans[ansSize][1] = -1;
    ansSize++;
    bCnt[x]--;
  }
  UpdateAns();
  Clean();
}

void Initialize2()
{
  Clean();
  rep(i, n)
  {
    int x = c[i][0];
    int z = c[i][1];
    int y = randxor() % m;
    while (y == x) {
      y = randxor() % m;
    }
    int minY = -1001;
    rep(j, m)
    {
      if (j == x) continue;
      int tmpMin = 1001001;
      rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
      if (tmpMin > minY) {
        minY = tmpMin;
        y = j;
      }
    }
    if (z < bCnt[x] - 1) {
      int ii = b[x][z + 1];
      ans[ansSize][0] = ii;
      ans[ansSize][1] = y;
      ansSize++;
      rep(j, bCnt[x] - (z + 1))
      {
        int jj = b[x][z + 1 + j];
        b[y][bCnt[y]] = b[x][z + 1 + j];
        c[jj][0] = y;
        c[jj][1] = bCnt[y];
        bCnt[y]++;
      }
      bCnt[x] = z + 1;
    }
    ans[ansSize][0] = i;
    ans[ansSize][1] = -1;
    ansSize++;
    bCnt[x]--;
  }
  UpdateAns();
  Clean();
}

void Initialize3()
{
  Clean();
  rep(i, n)
  {
    int x = c[i][0];
    int z = c[i][1];
    int y = randxor() % m;
    while (y == x) {
      y = randxor() % m;
    }

    int yy = -1;
    int zz = -1;
    drep(j, bCnt[x])
    {
      if (j == z) break;
      int jj = b[x][j];
      rep(k, m)
      {
        if (k == x) continue;
        if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5) {
          yy = k;
          zz = j;
          break;
        }
      }
      if (yy != -1) {
        break;
      }
    }
    if (yy != -1) {
      ans[ansSize][0] = b[x][zz];
      ans[ansSize][1] = yy;
      ansSize++;
      Move(x, yy, zz);
      i--;
      continue;
    }

    int minY = -1001;
    rep(j, m)
    {
      if (j == x) continue;
      int tmpMin = 1001001;
      rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
      if (tmpMin > minY) {
        minY = tmpMin;
        y = j;
      }
    }
    if (z < bCnt[x] - 1) {
      int ii = b[x][z + 1];
      ans[ansSize][0] = ii;
      ans[ansSize][1] = y;
      ansSize++;
      rep(j, bCnt[x] - (z + 1))
      {
        int jj = b[x][z + 1 + j];
        b[y][bCnt[y]] = b[x][z + 1 + j];
        c[jj][0] = y;
        c[jj][1] = bCnt[y];
        bCnt[y]++;
      }
      bCnt[x] = z + 1;
    }
    ans[ansSize][0] = i;
    ans[ansSize][1] = -1;
    ansSize++;
    bCnt[x]--;
  }
  UpdateAns();
  Clean();
}

int CalcBestY(int x, int z = -1)
{
  int y = -1;
  int minY = -1001;
  rep(j, m)
  {
    if (j == x) continue;
    if (j == gomi1 || j == gomi2 || j == gomi3) continue;
    int tmpMin = 1001001;
    rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
    if (tmpMin > minY) {
      minY = tmpMin;
      y = j;
    }
  }
  return y;
}

int CalcBestY2(int x, int z)
{
  int y = -1;
  int minCnt = 1001001;
  rep(j, m)
  {
    if (j == x) continue;
    if (j == gomi1 || j == gomi2 || j == gomi3) continue;
    int cnt = 0;
    rep(k, bCnt[j])
    {
      srep(l, z, bCnt[x])
      {
        if (b[j][k] < b[x][l]) {
          cnt++;
        }
      }
    }
    if (cnt < minCnt) {
      minCnt = cnt;
      y = j;
    }
  }
  return y;
}

void Initialize4()
{
  Clean();
  for (int karina = 30; karina < 190; karina += 3) {
    for (int winter = karina + 10; winter < 195; winter += 3) {
      rep(i, n)
      {
        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          // rep(k, m) {
          //   if (k == x) continue;
          //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5) {
          //     yy = k;
          //     zz = j;
          //     break;
          //   }
          // }
          // if (yy != -1) {
          //   break;
          // }
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              break;
            }
          }
        }
        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }
        y = CalcBestY(x);
        if (z < bCnt[x] - 1) {
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize5()
{
  Clean();
  for (int ningning = 30; ningning < 190; ningning += 5) {
    for (int karina = ningning + 10; karina < 190; karina += 5) {
      for (int winter = karina + 10; winter < 195; winter += 5) {
        rep(i, n)
        {
          int x = c[i][0];
          int z = c[i][1];
          int y = randxor() % m;
          while (y == x) {
            y = randxor() % m;
          }

          rep(j, m)
          {
            if (bCnt[j] == 0) {
              if (gomi1 == -1) {
                gomi1 = j;
              }
              else if (gomi2 == -1 && j != gomi1) {
                gomi2 = j;
              }
              else if (gomi3 == -1 && j != gomi1 && j != gomi2) {
                gomi3 = j;
              }
            }
          }

          int yy = -1;
          int zz = -1;
          drep(j, bCnt[x])
          {
            if (j == z) break;
            int jj = b[x][j];
            // rep(k, m) {
            //   if (k == x) continue;
            //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5)
            //   {
            //     yy = k;
            //     zz = j;
            //     break;
            //   }
            // }
            // if (yy != -1) {
            //   break;
            // }
            if (i < ningning) {
              if (jj >= winter && gomi1 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi1;
                zz = j;
                break;
              }
              else if (jj >= karina && gomi2 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi2;
                zz = j;
                break;
              }
              else if (jj >= ningning && gomi3 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi3;
                zz = j;
                break;
              }
            }
          }
          if (yy != -1) {
            ans[ansSize][0] = b[x][zz];
            ans[ansSize][1] = yy;
            ansSize++;
            Move(x, yy, zz);
            i--;
            continue;
          }

          // int minY = -1001;
          // rep(j, m) {
          //   if (j == x) continue;
          //   int tmpMin = 1001001;
          //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
          //   if (tmpMin > minY) {
          //     minY = tmpMin;
          //     y = j;
          //   }
          // }
          y = CalcBestY(x);
          if (z < bCnt[x] - 1) {
            int ii = b[x][z + 1];
            ans[ansSize][0] = ii;
            ans[ansSize][1] = y;
            ansSize++;
            rep(j, bCnt[x] - (z + 1))
            {
              int jj = b[x][z + 1 + j];
              b[y][bCnt[y]] = b[x][z + 1 + j];
              c[jj][0] = y;
              c[jj][1] = bCnt[y];
              bCnt[y]++;
            }
            bCnt[x] = z + 1;
          }
          ans[ansSize][0] = i;
          ans[ansSize][1] = -1;
          ansSize++;
          bCnt[x]--;
        }
        UpdateAns();
        Clean();
      }
    }
  }
}

void Initialize6()
{
  Clean();
  for (int karina = 30; karina < 190; karina += 3) {
    for (int winter = karina + 10; winter < 195; winter += 3) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (gomi1 == -1) {
          rep(j, m)
          {
            if (bCnt[j] == 1) {
              int yy = CalcBestY(j);
              if (yy == -1) continue;
              ans[ansSize][0] = b[j][0];
              ans[ansSize][1] = yy;
              ansSize++;
              Move(j, yy, 0);
              gomi1 = j;
              break;
            }
          }
        }
        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          // rep(k, m) {
          //   if (k == x) continue;
          //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5) {
          //     yy = k;
          //     zz = j;
          //     break;
          //   }
          // }
          // if (yy != -1) {
          //   break;
          // }
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              break;
            }
          }
        }
        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }
        y = CalcBestY(x);
        if (z < bCnt[x] - 1) {
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize7()
{
  Clean();
  for (int karina = 30; karina < 190; karina += 3) {
    for (int winter = karina + 10; winter < 195; winter += 3) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (gomi1 == -1) {
          rep(j, m)
          {
            if (bCnt[j] == 1) {
              int yy = CalcBestY(j);
              if (yy == -1) continue;
              ans[ansSize][0] = b[j][0];
              ans[ansSize][1] = yy;
              ansSize++;
              Move(j, yy, 0);
              gomi1 = j;
              break;
            }
          }
        }

        if (gomi1 != -1 && gomi2 == -1) {
          rep(j, m)
          {
            if (bCnt[j] == 1) {
              int yy = CalcBestY(j);
              if (yy == -1) continue;
              ans[ansSize][0] = b[j][0];
              ans[ansSize][1] = yy;
              ansSize++;
              Move(j, yy, 0);
              gomi2 = j;
              break;
            }
          }
        }
        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          // rep(k, m) {
          //   if (k == x) continue;
          //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5) {
          //     yy = k;
          //     zz = j;
          //     break;
          //   }
          // }
          // if (yy != -1) {
          //   break;
          // }
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              break;
            }
          }
        }
        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }
        y = CalcBestY(x);
        if (z < bCnt[x] - 1) {
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize8()
{
  Clean();
  for (int ningning = 31; ningning < 190; ningning += 5) {
    for (int karina = ningning + 1; karina < 190; karina += 5) {
      for (int winter = karina + 1; winter < 195; winter += 5) {
        rep(i, n)
        {
          rep(j, m)
          {
            if (bCnt[j] == 0) {
              if (gomi1 == -1) {
                gomi1 = j;
              }
              else if (gomi2 == -1 && j != gomi1) {
                gomi2 = j;
              }
              else if (gomi3 == -1 && j != gomi1 && j != gomi2) {
                gomi3 = j;
              }
            }
          }

          if (gomi1 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY(j);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi1 = j;
                break;
              }
            }
          }

          if (gomi1 != -1 && gomi2 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY(j);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi2 = j;
                break;
              }
            }
          }

          int x = c[i][0];
          int z = c[i][1];
          int y = randxor() % m;
          while (y == x) {
            y = randxor() % m;
          }

          int yy = -1;
          int zz = -1;
          drep(j, bCnt[x])
          {
            if (j == z) break;
            int jj = b[x][j];
            // rep(k, m) {
            //   if (k == x) continue;
            //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5)
            //   {
            //     yy = k;
            //     zz = j;
            //     break;
            //   }
            // }
            // if (yy != -1) {
            //   break;
            // }
            if (i < ningning) {
              if (jj >= winter && gomi1 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi1;
                zz = j;
                break;
              }
              else if (jj >= karina && gomi2 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi2;
                zz = j;
                break;
              }
              else if (jj >= ningning && gomi3 != -1) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi3;
                zz = j;
                break;
              }
            }
          }
          if (yy != -1) {
            ans[ansSize][0] = b[x][zz];
            ans[ansSize][1] = yy;
            ansSize++;
            Move(x, yy, zz);
            i--;
            continue;
          }

          // int minY = -1001;
          // rep(j, m) {
          //   if (j == x) continue;
          //   int tmpMin = 1001001;
          //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
          //   if (tmpMin > minY) {
          //     minY = tmpMin;
          //     y = j;
          //   }
          // }
          y = CalcBestY(x);
          if (z < bCnt[x] - 1) {
            int ii = b[x][z + 1];
            ans[ansSize][0] = ii;
            ans[ansSize][1] = y;
            ansSize++;
            rep(j, bCnt[x] - (z + 1))
            {
              int jj = b[x][z + 1 + j];
              b[y][bCnt[y]] = b[x][z + 1 + j];
              c[jj][0] = y;
              c[jj][1] = bCnt[y];
              bCnt[y]++;
            }
            bCnt[x] = z + 1;
          }
          ans[ansSize][0] = i;
          ans[ansSize][1] = -1;
          ansSize++;
          bCnt[x]--;
        }
        UpdateAns();
        Clean();
      }
    }
  }
}

void Initialize9()
{
  Clean();
  for (int karina = 30; karina < 190; karina += 3) {
    for (int winter = karina + 10; winter < 195; winter += 3) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (i < winter) {
          if (gomi1 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY(j);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi1 = j;
                break;
              }
            }
          }

          if (gomi1 != -1 && gomi2 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY(j);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi2 = j;
                break;
              }
            }
          }
        }

        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          // rep(k, m) {
          //   if (k == x) continue;
          //   if (jj + 1 <= b[k][bCnt[k] - 1] && b[k][bCnt[k] - 1] <= jj + 5) {
          //     yy = k;
          //     zz = j;
          //     break;
          //   }
          // }
          // if (yy != -1) {
          //   break;
          // }
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              while (zz > z + 1) {
                if (b[x][zz - 1] >= winter) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              break;
            }
          }
        }

        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }
        y = CalcBestY(x);
        if (z < bCnt[x] - 1) {
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize10()
{
  Clean();
  for (int karina = 30; karina < 190; karina += 3) {
    for (int winter = karina + 10; winter < 195; winter += 3) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (i < winter) {
          if (gomi1 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                if (b[j][0] >= winter) {
                  gomi1 = j;
                }
                else {
                  int yy = CalcBestY2(j, 0);
                  if (yy == -1) continue;
                  ans[ansSize][0] = b[j][0];
                  ans[ansSize][1] = yy;
                  ansSize++;
                  Move(j, yy, 0);
                  gomi1 = j;
                }

                break;
              }
            }
          }

          if (gomi1 != -1 && gomi2 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY2(j, 0);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi2 = j;
                break;
              }
            }
          }
        }

        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY2(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              while (zz > z + 1) {
                if (b[x][zz - 1] >= winter) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY2(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              break;
            }
          }
        }

        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }

        if (z < bCnt[x] - 1) {
          y = CalcBestY2(x, z + 1);
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize11()
{
  Clean();
  for (int karina = 31; karina < 190; karina += 1) {
    for (int winter = karina + 1; winter < 195; winter += 1) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (i < winter) {
          if (gomi1 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                if (b[j][0] >= winter) {
                  gomi1 = j;
                }
                else {
                  int yy = CalcBestY(j, 0);
                  if (yy == -1) continue;
                  ans[ansSize][0] = b[j][0];
                  ans[ansSize][1] = yy;
                  ansSize++;
                  Move(j, yy, 0);
                  gomi1 = j;
                }

                break;
              }
            }
          }

          if (gomi1 != -1 && gomi2 == -1) {
            rep(j, m)
            {
              if (bCnt[j] == 1) {
                int yy = CalcBestY(j, 0);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi2 = j;
                break;
              }
            }
          }
        }

        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              while (zz > z + 1) {
                if (b[x][zz - 1] >= winter) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              while (zz > z + 1) {
                if (winter > b[x][zz - 1] && b[x][zz - 1] >= karina) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
          }
        }

        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }

        if (z < bCnt[x] - 1) {
          y = CalcBestY(x, z + 1);
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize12()
{
  Clean();
  for (int karina = 31; karina < 190; karina += 1) {
    for (int winter = karina + 1; winter < 195; winter += 1) {
      rep(i, n)
      {
        rep(j, m)
        {
          if (bCnt[j] == 0) {
            if (gomi1 == -1) {
              gomi1 = j;
            }
            else if (gomi2 == -1 && j != gomi1) {
              gomi2 = j;
            }
          }
        }

        if (i < winter) {
          if (gomi1 == -1) {
            rep(j, m)
            {
              if (0 < bCnt[j] && bCnt[j] <= 3) {
                int yy = CalcBestY(j, 0);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi1 = j;

                break;
              }
            }
          }

          if (gomi1 != -1 && gomi2 == -1) {
            rep(j, m)
            {
              if (0 < bCnt[j] && bCnt[j] <= 3) {
                int yy = CalcBestY(j, 0);
                if (yy == -1) continue;
                ans[ansSize][0] = b[j][0];
                ans[ansSize][1] = yy;
                ansSize++;
                Move(j, yy, 0);
                gomi2 = j;
                break;
              }
            }
          }
        }

        int x = c[i][0];
        int z = c[i][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        int yy = -1;
        int zz = -1;
        drep(j, bCnt[x])
        {
          if (j == z) break;
          int jj = b[x][j];
          if (i < karina) {
            if (jj >= winter && gomi1 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi1;
              zz = j;
              while (zz > z + 1) {
                if (b[x][zz - 1] >= winter) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
            else if (jj >= karina && gomi2 != -1) {
              if (j < bCnt[x] - 1) {
                int yyy = CalcBestY(x, j + 1);
                if (yyy == -1) continue;
                ans[ansSize][0] = b[x][j + 1];
                ans[ansSize][1] = yyy;
                ansSize++;
                Move(x, yyy, j + 1);
              }

              yy = gomi2;
              zz = j;
              while (zz > z + 1) {
                if (winter > b[x][zz - 1] && b[x][zz - 1] >= karina) {
                  zz--;
                }
                else {
                  break;
                }
              }
              break;
            }
          }
        }

        if (yy != -1) {
          ans[ansSize][0] = b[x][zz];
          ans[ansSize][1] = yy;
          ansSize++;
          Move(x, yy, zz);
          i--;
          continue;
        }

        // int minY = -1001;
        // rep(j, m) {
        //   if (j == x) continue;
        //   int tmpMin = 1001001;
        //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
        //   if (tmpMin > minY) {
        //     minY = tmpMin;
        //     y = j;
        //   }
        // }

        if (z < bCnt[x] - 1) {
          y = CalcBestY(x, z + 1);
          int ii = b[x][z + 1];
          ans[ansSize][0] = ii;
          ans[ansSize][1] = y;
          ansSize++;
          rep(j, bCnt[x] - (z + 1))
          {
            int jj = b[x][z + 1 + j];
            b[y][bCnt[y]] = b[x][z + 1 + j];
            c[jj][0] = y;
            c[jj][1] = bCnt[y];
            bCnt[y]++;
          }
          bCnt[x] = z + 1;
        }
        ans[ansSize][0] = i;
        ans[ansSize][1] = -1;
        ansSize++;
        bCnt[x]--;
      }
      UpdateAns();
      Clean();
    }
  }
}

void Initialize13()
{
  Clean();
  for (int karina = 31; karina < 190; karina += 3) {
    for (int winter = karina + 1; winter < 195; winter += 3) {
      for (int eri = 10; eri <= 100; eri += 10) {
        rep(i, n)
        {
          rep(j, m)
          {
            if (bCnt[j] == 0) {
              if (gomi1 == -1) {
                gomi1 = j;
              }
              else if (gomi2 == -1 && j != gomi1) {
                gomi2 = j;
              }
            }
          }

          if (i < winter) {
            if (gomi1 == -1) {
              rep(j, m)
              {
                if (0 < bCnt[j] && bCnt[j] <= 3) {
                  int yy = CalcBestY(j, 0);
                  if (yy == -1) continue;
                  ans[ansSize][0] = b[j][0];
                  ans[ansSize][1] = yy;
                  ansSize++;
                  Move(j, yy, 0);
                  gomi1 = j;

                  break;
                }
              }
            }

            if (gomi1 != -1 && gomi2 == -1) {
              rep(j, m)
              {
                if (0 < bCnt[j] && bCnt[j] <= 3) {
                  int yy = CalcBestY(j, 0);
                  if (yy == -1) continue;
                  ans[ansSize][0] = b[j][0];
                  ans[ansSize][1] = yy;
                  ansSize++;
                  Move(j, yy, 0);
                  gomi2 = j;
                  break;
                }
              }
            }
          }

          int x = c[i][0];
          int z = c[i][1];
          int y = randxor() % m;
          while (y == x) {
            y = randxor() % m;
          }

          int yy = -1;
          int zz = -1;
          drep(j, bCnt[x])
          {
            if (j == z) break;
            int jj = b[x][j];
            if (i < karina) {
              if (jj >= winter && gomi1 != -1 && bCnt[gomi1] < eri) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x, j + 1);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi1;
                zz = j;
                while (zz > z + 1) {
                  if (b[x][zz - 1] >= winter) {
                    zz--;
                  }
                  else {
                    break;
                  }
                }
                break;
              }
              else if (jj >= karina && gomi2 != -1 && bCnt[gomi2] < eri) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x, j + 1);
                  if (yyy == -1) continue;
                  ans[ansSize][0] = b[x][j + 1];
                  ans[ansSize][1] = yyy;
                  ansSize++;
                  Move(x, yyy, j + 1);
                }

                yy = gomi2;
                zz = j;
                while (zz > z + 1) {
                  if (winter > b[x][zz - 1] && b[x][zz - 1] >= karina) {
                    zz--;
                  }
                  else {
                    break;
                  }
                }
                break;
              }
            }
          }

          if (yy != -1) {
            ans[ansSize][0] = b[x][zz];
            ans[ansSize][1] = yy;
            ansSize++;
            Move(x, yy, zz);
            i--;
            continue;
          }

          // int minY = -1001;
          // rep(j, m) {
          //   if (j == x) continue;
          //   int tmpMin = 1001001;
          //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
          //   if (tmpMin > minY) {
          //     minY = tmpMin;
          //     y = j;
          //   }
          // }

          if (z < bCnt[x] - 1) {
            y = CalcBestY(x, z + 1);
            int ii = b[x][z + 1];
            ans[ansSize][0] = ii;
            ans[ansSize][1] = y;
            ansSize++;
            rep(j, bCnt[x] - (z + 1))
            {
              int jj = b[x][z + 1 + j];
              b[y][bCnt[y]] = b[x][z + 1 + j];
              c[jj][0] = y;
              c[jj][1] = bCnt[y];
              bCnt[y]++;
            }
            bCnt[x] = z + 1;
          }
          ans[ansSize][0] = i;
          ans[ansSize][1] = -1;
          ansSize++;
          bCnt[x]--;
        }
        UpdateAns();
        Clean();
      }
    }
  }
}

void Shuffle(int x, int y, int z)
{
  int kijun = 0;
  int arg = -1;
  srep(i, z, bCnt[x])
  {
    srep(j, i + 1, bCnt[x])
    {
      if (b[x][j] > b[x][i]) kijun++;
    }
  }
  srep(i, z + 1, bCnt[x])
  {
    vector<int> v;
    srep(j, i, bCnt[x]) { v.push_back(b[x][j]); }
    srep(j, z, i) { v.push_back(b[x][j]); }
    int cnt = 0;
    rep(j, v.size())
    {
      srep(k, j + 1, v.size())
      {
        if (v[k] > v[j]) cnt++;
      }
    }
    if (cnt < kijun) {
      kijun = cnt;
      arg = i;
    }
  }
  if (arg != -1) {
    ans[ansSize][0] = b[x][arg];
    ans[ansSize][1] = y;
    ansSize++;
    Move(x, y, arg);
    ans[ansSize][0] = b[x][z];
    ans[ansSize][1] = y;
    ansSize++;
    Move(x, y, z);
  }
  else {
    ans[ansSize][0] = b[x][z];
    ans[ansSize][1] = y;
    ansSize++;
    Move(x, y, z);
  }
}

void Initialize14()
{
  Clean();
  for (int karina = 31; karina < 190; karina += 3) {
    for (int winter = karina + 1; winter < 195; winter += 3) {
      for (int eri = 10; eri <= 100; eri += 10) {
        e_time = clock();
        double now_time = (double)(e_time - s_time) / CLOCKS_PER_SEC;
        if (now_time > TL) {
          break;
        }
        rep(i, n)
        {
          rep(j, m)
          {
            if (bCnt[j] == 0) {
              if (gomi1 == -1) {
                gomi1 = j;
              }
              else if (gomi2 == -1 && j != gomi1) {
                gomi2 = j;
              }
            }
          }

          if (i < winter) {
            if (gomi1 == -1) {
              rep(j, m)
              {
                if (0 < bCnt[j] && bCnt[j] <= 3) {
                  int yy = CalcBestY(j, 0);
                  if (yy == -1) continue;
                  Shuffle(j, yy, 0);
                  gomi1 = j;

                  break;
                }
              }
            }

            if (gomi1 != -1 && gomi2 == -1) {
              rep(j, m)
              {
                if (0 < bCnt[j] && bCnt[j] <= 3) {
                  int yy = CalcBestY(j, 0);
                  if (yy == -1) continue;
                  Shuffle(j, yy, 0);
                  gomi2 = j;
                  break;
                }
              }
            }
          }

          int x = c[i][0];
          int z = c[i][1];
          int y = randxor() % m;
          while (y == x) {
            y = randxor() % m;
          }

          int yy = -1;
          int zz = -1;
          drep(j, bCnt[x])
          {
            if (j == z) break;
            int jj = b[x][j];
            if (i < karina) {
              if (jj >= winter && gomi1 != -1 && bCnt[gomi1] < eri) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x, j + 1);
                  if (yyy == -1) continue;
                  Shuffle(x, yyy, j + 1);
                }

                yy = gomi1;
                zz = j;
                while (zz > z + 1) {
                  if (b[x][zz - 1] >= winter) {
                    zz--;
                  }
                  else {
                    break;
                  }
                }
                break;
              }
              else if (jj >= karina && gomi2 != -1 && bCnt[gomi2] < eri) {
                if (j < bCnt[x] - 1) {
                  int yyy = CalcBestY(x, j + 1);
                  if (yyy == -1) continue;
                  Shuffle(x, yyy, j + 1);
                }

                yy = gomi2;
                zz = j;
                while (zz > z + 1) {
                  if (winter > b[x][zz - 1] && b[x][zz - 1] >= karina) {
                    zz--;
                  }
                  else {
                    break;
                  }
                }
                break;
              }
            }
          }

          if (yy != -1) {
            Shuffle(x, yy, zz);
            i--;
            continue;
          }

          // int minY = -1001;
          // rep(j, m) {
          //   if (j == x) continue;
          //   int tmpMin = 1001001;
          //   rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
          //   if (tmpMin > minY) {
          //     minY = tmpMin;
          //     y = j;
          //   }
          // }

          if (z < bCnt[x] - 1) {
            y = CalcBestY(x, z + 1);
            Shuffle(x, y, z + 1);
          }
          ans[ansSize][0] = i;
          ans[ansSize][1] = -1;
          ansSize++;
          bCnt[x]--;
        }
        UpdateAns();
        Clean();
      }
    }
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    rep(i, ansSize) { cout << ans[i][0] + 1 << ' ' << ans[i][1] + 1 << endl; }
  }
  else {
    rep(i, ansSize) { ofs << ans[i][0] + 1 << ' ' << ans[i][1] + 1 << endl; }
  }
}

void Method1()
{
  while (true) {
    e_time = clock();
    double now_time = (double)(e_time - s_time) / CLOCKS_PER_SEC;
    if (now_time > TL) {
      break;
    }

    int startX = randxor() % n;
    CopyAll();
    RollbackAns();

    int i = 0;
    int now = 0;
    while (now < n) {
      if (now >= startX) {
        int x = c[now][0];
        int z = c[now][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }
        if (z < bCnt[x] - 1) {
          Move(x, y, z + 1);
          int ii = b[x][z + 1];
          ans[i][0] = ii;
          ans[i][1] = y;
          i++;
        }
        ans[i][0] = now;
        ans[i][1] = -1;
        bCnt[x]--;
        i++;
        now++;
      }
      else {
        if (ans[i][1] == -1) {
          int x = c[ans[i][0]][0];
          bCnt[x]--;
          i++;
          now++;
        }
        else {
          int x = c[ans[i][0]][0];
          int z = c[ans[i][0]][1];
          int y = ans[i][1];
          Move(x, y, z);
          i++;
        }
      }
    }
    ansSize = i;
    int tmpScore = CalcScore();
    if (tmpScore >= maxScore) {
      CopyAns();
    }
  }
}

void Method2()
{
  int loop = 0;
  while (true) {
    e_time = clock();
    double now_time = (double)(e_time - s_time) / CLOCKS_PER_SEC;
    if (now_time > TL) {
      break;
    }
    loop++;
    int startX = randxor() % n;
    CopyAll();
    RollbackAns();

    int i = 0;
    int now = 0;
    while (now < n) {
      if (now >= startX) {
        int x = c[now][0];
        int z = c[now][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }
        if (now != startX) {
          int minY = -1001;
          rep(j, m)
          {
            if (j == x) continue;
            int tmpMin = 1001001;
            rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
            if (tmpMin > minY) {
              minY = tmpMin;
              y = j;
            }
          }
        }

        if (z < bCnt[x] - 1) {
          Move(x, y, z + 1);
          int ii = b[x][z + 1];
          ans[i][0] = ii;
          ans[i][1] = y;
          i++;
        }
        ans[i][0] = now;
        ans[i][1] = -1;
        bCnt[x]--;
        i++;
        now++;
      }
      else {
        if (ans[i][1] == -1) {
          int x = c[ans[i][0]][0];
          bCnt[x]--;
          i++;
          now++;
        }
        else {
          int x = c[ans[i][0]][0];
          int z = c[ans[i][0]][1];
          int y = ans[i][1];
          Move(x, y, z);
          i++;
        }
      }
    }
    ansSize = i;
    int tmpScore = CalcScore();
    if (tmpScore >= maxScore) {
      CopyAns();
    }
  }
  // if (mode != 0) {
  //   cout << loop << endl;
  // }
}

void Method3()
{
  int loop = 0;
  while (true) {
    e_time = clock();
    double now_time = (double)(e_time - s_time) / CLOCKS_PER_SEC;
    if (now_time > TL) {
      break;
    }
    loop++;
    int startX = randxor() % n;
    CopyAll();
    RollbackAns();

    int i = 0;
    int now = 0;
    int use = 0;
    while (now < n) {
      if (now >= startX) {
        int x = c[now][0];
        int z = c[now][1];
        int y = randxor() % m;
        while (y == x) {
          y = randxor() % m;
        }

        if (randxor() % 10 == 0 && use == 0) {
          int ma = -1001;
          int pos = -1;
          srep(j, z + 1, bCnt[x])
          {
            if (b[x][j] > ma) {
              ma = b[x][j];
              pos = j;
            }
          }
          if (pos != -1) {
            int minY = -1001;
            rep(j, m)
            {
              if (j == x) continue;
              int tmpMin = 1001001;
              rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
              if (tmpMin > minY) {
                minY = tmpMin;
                y = j;
              }
            }
            ans[i][0] = ma;
            ans[i][1] = y;
            Move(x, y, pos);
            i++;
            use = 1;
          }

          continue;
        }
        else {
          int minY = -1001;
          rep(j, m)
          {
            if (j == x) continue;
            int tmpMin = 1001001;
            rep(k, bCnt[j]) { tmpMin = min(tmpMin, b[j][k]); }
            if (tmpMin > minY) {
              minY = tmpMin;
              y = j;
            }
          }
        }

        if (z < bCnt[x] - 1) {
          int ii = b[x][z + 1];
          ans[i][0] = ii;
          ans[i][1] = y;
          Move(x, y, z + 1);
          i++;
        }
        ans[i][0] = now;
        ans[i][1] = -1;
        bCnt[x]--;
        i++;
        now++;
      }
      else {
        if (ans[i][1] == -1) {
          int x = c[ans[i][0]][0];
          bCnt[x]--;
          i++;
          now++;
        }
        else {
          int x = c[ans[i][0]][0];
          int z = c[ans[i][0]][1];
          int y = ans[i][1];
          Move(x, y, z);
          i++;
        }
      }
    }
    ansSize = i;
    int tmpScore = CalcScore();
    if (tmpScore >= maxScore) {
      CopyAns();
    }
  }
  // if (mode != 0) {
  //   cout << loop << endl;
  // }
}

ll Solve(int probNum)
{
  s_time = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  // Initialize();
  // Initialize2();
  // Initialize3();
  // Initialize4();
  // Initialize5();
  // Initialize6();
  // Initialize7();
  // Initialize8();
  // Initialize9();
  // Initialize10();
  // Initialize11();
  // Initialize12();
  // Initialize13();
  Initialize14();
  UpdateAns();

  // Method1();
  // Method2();
  // Method3();

  // 解答を出力
  RollbackAns();
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
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 10)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
