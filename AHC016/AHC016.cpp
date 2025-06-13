#include "Hypers.h"

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
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

const int dx[4] = { -1, 1, 0, 0 };
const int dy[4] = { 0, 0, -1, 1 };

namespace /* 乱数ライブラリ */
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


  static double Rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

namespace
{
  // 共通変数
  double com[105][105];
  int MODE = 0;
  const int TURN = 100;
  int m, iEps;
  double eps;

  int n;
  int a[100][100][100];
  int b[100][100];

  int numSingleArr[1000];
  int numPairArr[1000][2];
  int numPairArrOK;

  int numThreeArr[1000][3];
  int numFourArr[1000][4];

  // ハイパラ調整用
  int hyperSolverNum;
  int hyperMinDiff = 10;
  int hyperMaxRound = 7;
  int hyperStep1 = 1;
  int hyperStep2 = 1;
}  // namespace

// ハイパラ
namespace
{
  int maxNumArray[100];
  int real_real_maxNumArray[11][100];
  int real_maxNumArray[100];
  int hyperN[101][41];
  double hyperMaxScore[101][41];
  int hyperSolver[101][41];
  int hyperMinDiffArr[101][41];
  int hyperMaxRoundArr[101][41];
  int hyperStep1Arr[101][41];
  int hyperStep2Arr[101][41];
}  // namespace

// Input（m, eps, iEpsの設定）
namespace
{
  int judgeArr[100];
  void Input(int mode, int problemNum = 0)
  {
    if (mode == 0) {
      cin >> m;
      string sEps;
      cin >> sEps;
      iEps = (sEps[2] - '0') * 10 + (sEps[3] - '0');
      eps = (double)iEps / 100.0;
    }
    else if (mode == 1000) {
      string fileNameIfs = "./in/";
      string strNum;
      for (int i = 0; i < (4); ++i) {
        strNum += (char)(problemNum % 10 + '0');
        problemNum /= 10;
      }
      reverse(strNum.begin(), strNum.end());
      fileNameIfs += strNum + ".txt";

      ifstream ifs(fileNameIfs);

      ifs >> m;
      string sEps;
      ifs >> sEps;
      iEps = (sEps[2] - '0') * 10 + (sEps[3] - '0');
      eps = (double)iEps / 100.0;

      for (int i = 0; i < (100); ++i) ifs >> judgeArr[i];
    }
    else if (mode == 100 || mode == 110) {
      m = Rand() % 91 + 10;
      // m = Rand() % 41 + 10;
      iEps = Rand() % 41;
      // iEps = Rand() % 11 + 30;
      eps = iEps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = Rand() % m;
    }
    else if (mode == 333) {
      m = problemNum * 10;
      iEps = 40;
      eps = iEps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = Rand() % m;
    }
  }
}  // namespace

// Output
namespace
{
  ofstream ofs1000Out;
  void OpenOfs1000Out(int problemNum)
  {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < (4); ++i) {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs1000Out.open(fileNameOfs);
  }

  void CloseOfs1000Out() { ofs1000Out.close(); }

  void OutputArrayAsString(int mode)
  {
    if (mode == 0) {
      cout << n << endl;
      for (int i = 0; i < (m); ++i) {
        string s;
        for (int j = 0; j < (n); ++j) {
          for (int k = j + 1; k < n; ++k) { s += (char)(a[i][j][k] + '0'); }
        }
        cout << s << endl;
      }
      fflush(stdout);
    }
    else if (mode == 1000) {
      ofs1000Out << n << endl;
      ofs1000Out << "# n = " << n << endl;
      ofs1000Out << "# hyperSolverNum = " << hyperSolverNum << endl;
      ofs1000Out << "# hyperMinDiff = " << hyperMinDiff << endl;
      ofs1000Out << "# hyperMaxRound = " << hyperMaxRound << endl;
      ofs1000Out << "# hyperStep1 = " << hyperStep1 << endl;
      ofs1000Out << "# hyperStep2 = " << hyperStep2 << endl;
      for (int i = 0; i < (m); ++i) {
        string s;
        for (int j = 0; j < (n); ++j) {
          for (int k = j + 1; k < n; ++k) { s += (char)(a[i][j][k] + '0'); }
        }
        ofs1000Out << s << endl;
        ofs1000Out << "# " << numPairArr[i][0] << ' ' << numPairArr[i][1] << endl;
      }
    }
  }

  int answersFor1000Out[TURN];
  void OutputAnsToOfs1000Out()
  {
    for (int i = 0; i < (100); ++i) { ofs1000Out << answersFor1000Out[i] << endl; }
  }

  // ハイパーパラメータ配列を出力する共通関数
  template<typename T>
  void OutputHyperArray(ofstream& ofs, const string& typeName, const string& arrayName, T arr[101][41])
  {
    ofs << typeName << " " << arrayName << "[101][41] = {" << endl;
    for (int i = 0; i < (101); ++i) {
      ofs << "{";
      for (int j = 0; j < (41); ++j) {
        ofs << arr[i][j];
        if (j == 40)
          ofs << "}";
        else
          ofs << ",";
      }
      if (i == 100) {
        ofs << "};" << endl;
      }
      else {
        ofs << ',' << endl;
      }
    }
    ofs << endl;
  }

  void OutputHaipara()
  {
    for (int i = 0; i < (101); ++i) {
      for (int j = 0; j < (41); ++j) {
        if (hyperSolver[i][j] / 10 == 10 || hyperSolver[i][j] / 10 == 12) {
          hyperStep2Arr[i][j] = hyperStep1Arr[i][j];
        }
      }
    }
    ofstream ofs("Haipara.txt");
    OutputHyperArray(ofs, "int", "hyperN", hyperN);
    OutputHyperArray(ofs, "double", "hyperMaxScore", hyperMaxScore);
    OutputHyperArray(ofs, "int", "hyperSolver", hyperSolver);
    OutputHyperArray(ofs, "int", "hyperMinDiffArr", hyperMinDiffArr);
    OutputHyperArray(ofs, "int", "hyperMaxRoundArr", hyperMaxRoundArr);
    OutputHyperArray(ofs, "int", "hyperStep1Arr", hyperStep1Arr);
    OutputHyperArray(ofs, "int", "hyperStep2Arr", hyperStep2Arr);
    ofs.close();
  }
}  // namespace

void RandmizeGraph(int x)
{
  for (int j = 0; j < (n); ++j) {
    for (int k = j + 1; k < n; ++k) {
      b[j][k] = a[x][j][k];
      if (Rand01() < eps) {
        b[j][k] = 1 - b[j][k];
      }
      b[k][j] = b[j][k];
    }
  }
}

int judgeNum;
void InitB(int mode, int turn = 0)
{
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { b[i][j] = 0; }
  }
  if (mode == 0) {
    string s;
    cin >> s;
    int ite = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        b[i][j] = s[ite] - '0';
        b[j][i] = b[i][j];
        ite++;
      }
    }
  }
  else {
    int judge = judgeArr[turn];
    judgeNum = judge;
    RandmizeGraph(judge);
  }
}

// numArray
namespace
{
  int numArr[100];

  // 共通のグラフ初期化関数
  void SetGraphFromNumArray(int arraySize = 100)
  {
    for (int i = 0; i < (m); ++i) {
      int num = numArr[i];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
        }
      }
    }
  }

  // 条件付きグラフ初期化関数（2つの数値用）
  void SetGraphFromPairArray()
  {
    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 <= j && k < num1 + num2) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
        }
      }
    }
  }

  // 条件付きグラフ初期化関数（3つの数値用）
  void SetGraphFromThreeArray()
  {
    for (int i = 0; i < (m); ++i) {
      int num1 = numThreeArr[i][0];
      int num2 = numThreeArr[i][1];
      int num3 = numThreeArr[i][2];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 <= j && k < num1 + num2) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 + num2 <= j && k < num1 + num2 + num3) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
        }
      }
    }
  }

  // 条件付きグラフ初期化関数（4つの数値用）
  void SetGraphFromFourArray()
  {
    for (int i = 0; i < (m); ++i) {
      int num1 = numFourArr[i][0];
      int num2 = numFourArr[i][1];
      int num3 = numFourArr[i][2];
      int num4 = numFourArr[i][3];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 <= j && k < num1 + num2) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 + num2 <= j && k < num1 + num2 + num3) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 + num2 + num3 <= j && k < num1 + num2 + num3 + num4) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
        }
      }
    }
  }

  void InitNumArray1()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
    }
    else {
      int cnt = 0;
      numArr[cnt] = n / 2 + 1;
      cnt++;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
    }

    SetGraphFromNumArray();
  }

  void InitNumArray2()
  {
    for (int i = 0; i < (100); ++i) { numArr[i] = maxNumArray[i]; }

    SetGraphFromNumArray();
  }

  void InitNumArray3()
  {
    for (int i = 0; i < (100); ++i) { numArr[i] = real_real_maxNumArray[(m + 9) / 10][i]; }

    SetGraphFromNumArray();
  }

  void InitNumArray4()
  {
    for (int i = 0; i < (100); ++i) numArr[i] = 0;
    for (int i = 0; i < (20); ++i) { numArr[i] = (i + 1) * 5; }
    SetGraphFromNumArray();
  }

  void InitNumArray5()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 1) continue;
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 0) continue;
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
    }
    else {
      int cnt = 0;
      numArr[cnt] = n / 2 + 1;
      cnt++;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 0) continue;
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 1) continue;
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
    }

    SetGraphFromNumArray();
  }

  void InitNumArray6()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int j = 0; j < (3); ++j) {
        for (int i = (n / 2) - 1; i >= 0; --i) {
          if (i % 3 != j) continue;
          numArr[cnt] = i + 1;
          cnt++;
          numArr[cnt] = n - i;
          cnt++;
        }
      }
    }
    else {
      int cnt = 0;
      numArr[cnt] = n / 2 + 1;
      cnt++;
      for (int j = 2; j < 5; ++j) {
        for (int i = (n / 2) - 1; i >= 0; --i) {
          if (i % 3 != j % 3) continue;
          numArr[cnt] = i + 1;
          cnt++;
          numArr[cnt] = n - i;
          cnt++;
        }
      }
    }

    SetGraphFromNumArray();
  }

  void InitNumArray7()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = i + 1;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = n - i;
        cnt++;
      }
    }
    else {
      int cnt = 0;
      numArr[cnt] = n / 2 + 1;
      cnt++;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = i + 1;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        numArr[cnt] = n - i;
        cnt++;
      }
    }

    for (int i = 0; i < (m); ++i) {
      int num = numArr[i];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num <= j) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
        }
      }
    }
  }

  void InitNumArray8()
  {
    int cnt = 0;
    for (int i = (n)-1; i >= 0; --i) {
      numArr[cnt] = i + 1;
      cnt++;
    }

    SetGraphFromNumArray();
  }

  void InitNumArray9()
  {
    int cnt = 0;
    for (int j = 1; j < 3; ++j) {
      for (int i = (n)-1; i >= 0; --i) {
        if (i % 2 != j % 2) continue;
        numArr[cnt] = i + 1;
        cnt++;
      }
    }

    SetGraphFromNumArray();
  }

  void InitNumArray10()
  {
    numPairArrOK = 1;
    int cnt = 0;
    for (int i = n; i > 0; i -= hyperStep1) {
      numPairArr[cnt][0] = i;
      numPairArr[cnt][1] = 0;
      cnt++;
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        cnt++;
        j -= hyperStep1;
      }
      if (cnt > 200) break;
    }

    SetGraphFromPairArray();
    if (cnt < m) numPairArrOK = 0;
  }

  void InitNumArray11()
  {
    numPairArrOK = 1;
    int cnt = 0;
    for (int i = n; i > 0; i -= hyperStep1) {
      numPairArr[cnt][0] = i;
      numPairArr[cnt][1] = 0;
      cnt++;
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        cnt++;
        j -= hyperStep2;
      }
      if (cnt > 200) break;
    }

    SetGraphFromPairArray();
    if (cnt < m) numPairArrOK = 0;
  }

  void InitNumArray12()
  {
    numPairArrOK = 1;
    int cnt = 0;
    int one = 1;
    for (int i = n; i > 0; i -= hyperStep1) {
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        one = 0;
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        cnt++;
        j -= hyperStep1;
      }
      if (one) {
        numPairArr[cnt][0] = i;
        numPairArr[cnt][1] = 0;
        cnt++;
      }
      if (cnt > 200) break;
    }

    SetGraphFromPairArray();
    if (cnt < m) numPairArrOK = 0;
  }

  void InitNumArray13()
  {
    numPairArrOK = 1;
    int cnt = 0;
    int one = 1;
    for (int i = n; i > 0; i -= hyperStep1) {
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        one = 0;
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        cnt++;
        j -= hyperStep2;
      }
      if (one) {
        numPairArr[cnt][0] = i;
        numPairArr[cnt][1] = 0;
        cnt++;
      }
      if (cnt > 200) break;
    }

    SetGraphFromPairArray();
    if (cnt < m) numPairArrOK = 0;
  }

  void InitNumArray14()
  {
    for (int i = 0; i < (m); ++i) {
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (cnt < i) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
          cnt++;
        }
      }
    }
  }

  void InitNumArray15()
  {
    numPairArrOK = 0;
    int minn = 4;
    int cnt = 0;
    set<int> iset;
    set<P> jset;
    while (true) {
      cnt = 0;
      int k = minn;
      while (k <= n) {
        int j = k + hyperMinDiff;
        while (j <= n) {
          int i = j + hyperMinDiff;
          while (i <= n) {
            // 1クラ
            if (iset.find(i) == iset.end()) {
              iset.insert(i);
              cnt++;
            }
            // 2クラ
            if (i + j <= n) {
              if (jset.find(P(i, j)) == jset.end()) {
                jset.insert(P(i, j));
                cnt++;
              }
            }

            // 3クラ
            if (i + j + k <= n) {
              cnt++;
            }
            i += hyperStep1;
          }
          j += hyperStep2;
        }
        k += hyperStep2;
        if (cnt >= m) break;
      }
      iset.clear();
      jset.clear();
      if (cnt >= m) {
        minn++;
        numPairArrOK = 1;
      }
      else {
        break;
      }
    }
    minn--;

    {
      cnt = 0;
      int k = minn;
      while (k <= n) {
        int j = k + hyperMinDiff;
        while (j <= n) {
          int i = j + hyperMinDiff;
          while (i <= n) {
            // 1クラ
            if (iset.find(i) == iset.end()) {
              iset.insert(i);
              numThreeArr[cnt][0] = i;
              numThreeArr[cnt][1] = 0;
              numThreeArr[cnt][2] = 0;
              cnt++;
            }
            // 2クラ
            if (i + j <= n) {
              if (jset.find(P(i, j)) == jset.end()) {
                jset.insert(P(i, j));
                numThreeArr[cnt][0] = i;
                numThreeArr[cnt][1] = j;
                numThreeArr[cnt][2] = 0;
                cnt++;
              }
            }

            // 3クラ
            if (i + j + k <= n) {
              numThreeArr[cnt][0] = i;
              numThreeArr[cnt][1] = j;
              numThreeArr[cnt][2] = k;
              cnt++;
            }
            i += hyperStep1;
            if (cnt >= 500) break;
          }
          j += hyperStep2;
          if (cnt >= 500) break;
        }
        k += hyperStep2;
        if (cnt >= 500) break;
      }
    }

    for (int i = 0; i < (cnt / 2); ++i) {
      swap(numThreeArr[i][0], numThreeArr[cnt - 1 - i][0]);
      swap(numThreeArr[i][1], numThreeArr[cnt - 1 - i][1]);
      swap(numThreeArr[i][2], numThreeArr[cnt - 1 - i][2]);
    }

    SetGraphFromThreeArray();
  }

  void InitNumArray16()
  {
    numPairArrOK = 0;
    int minn = 4;
    int cnt = 0;
    set<int> iset;
    set<P> jset;
    while (true) {
      cnt = 0;
      int k = minn;
      while (k <= n) {
        int j = k + hyperMinDiff;
        while (j <= n) {
          int i = j + hyperMinDiff;
          while (i <= n) {
            // 1クラ
            // if (iset.find(i) == iset.end()) {
            //  iset.insert(i);
            //  cnt++;
            //}
            // 2クラ
            // if (i + j <= n) {
            //  if (jset.find(P(i, j)) == jset.end()) {
            //    jset.insert(P(i, j));
            //    cnt++;
            //  }
            //}

            // 3クラ
            if (i + j + k <= n) {
              cnt++;
            }
            i += hyperStep1;
          }
          j += hyperStep2;
        }
        k += hyperStep2;
        if (cnt >= m) break;
      }
      iset.clear();
      jset.clear();
      if (cnt >= m) {
        minn++;
        numPairArrOK = 1;
      }
      else {
        break;
      }
    }
    minn--;

    {
      cnt = 0;
      int k = minn;
      while (k <= n) {
        int j = k + hyperMinDiff;
        while (j <= n) {
          int i = j + hyperMinDiff;
          while (i <= n) {
            //// 1クラ
            // if (iset.find(i) == iset.end()) {
            //   iset.insert(i);
            //   numThreeArr[cnt][0] = i;
            //   numThreeArr[cnt][1] = 0;
            //   numThreeArr[cnt][2] = 0;
            //   cnt++;
            // }
            //// 2クラ
            // if (i + j <= n) {
            //   if (jset.find(P(i, j)) == jset.end()) {
            //     jset.insert(P(i, j));
            //     numThreeArr[cnt][0] = i;
            //     numThreeArr[cnt][1] = j;
            //     numThreeArr[cnt][2] = 0;
            //     cnt++;
            //   }
            // }

            // 3クラ
            if (i + j + k <= n) {
              numThreeArr[cnt][0] = i;
              numThreeArr[cnt][1] = j;
              numThreeArr[cnt][2] = k;
              cnt++;
            }
            i += hyperStep1;
            if (cnt >= 500) break;
          }
          j += hyperStep2;
          if (cnt >= 500) break;
        }
        k += hyperStep2;
        if (cnt >= 500) break;
      }
    }

    for (int i = 0; i < (cnt / 2); ++i) {
      swap(numThreeArr[i][0], numThreeArr[cnt - 1 - i][0]);
      swap(numThreeArr[i][1], numThreeArr[cnt - 1 - i][1]);
      swap(numThreeArr[i][2], numThreeArr[cnt - 1 - i][2]);
    }

    SetGraphFromThreeArray();
  }

  // 4コア
  void InitNumArray17()
  {
    numPairArrOK = 0;
    int minn = 4;
    int cnt = 0;
    set<int> iset;
    set<P> jset;
    set<pair<P, P>> kset;
    pair<P, P> p;
    p.second.second = 0;
    while (true) {
      cnt = 0;
      int l = minn;
      while (l <= n) {
        int k = l + hyperMinDiff;
        while (k <= n) {
          int j = k + hyperMinDiff;
          while (j <= n) {
            int i = j + hyperMinDiff;
            while (i <= n) {
              // 1クラ
              if (iset.find(i) == iset.end()) {
                iset.insert(i);
                cnt++;
              }
              // 2クラ
              if (i + j <= n) {
                if (jset.find(P(i, j)) == jset.end()) {
                  jset.insert(P(i, j));
                  cnt++;
                }
              }

              // 3クラ
              if (i + j + k <= n) {
                p.first.first = i;
                p.first.second = j;
                p.second.first = k;
                if (kset.find(p) == kset.end()) {
                  kset.insert(p);
                  cnt++;
                }
              }
              // 4クラ
              if (i + j + k + l <= n) {
                cnt++;
              }
              i += hyperStep1;
            }
            j += hyperStep1;
          }
          k += hyperStep2;
          if (cnt >= m) break;
        }
        l += hyperStep2;
        if (cnt >= m) break;
      }

      iset.clear();
      jset.clear();
      kset.clear();
      if (cnt >= m) {
        minn++;
        numPairArrOK = 1;
      }
      else {
        break;
      }
    }
    minn--;

    if (numPairArrOK == 0) return;

    {
      cnt = 0;
      int l = minn;
      while (l <= n) {
        int k = l + hyperMinDiff;
        while (k <= n) {
          int j = k + hyperMinDiff;
          while (j <= n) {
            int i = j + hyperMinDiff;
            while (i <= n) {
              // 1クラ
              if (iset.find(i) == iset.end()) {
                iset.insert(i);
                numFourArr[cnt][0] = i;
                numFourArr[cnt][1] = 0;
                numFourArr[cnt][2] = 0;
                numFourArr[cnt][3] = 0;
                cnt++;
              }
              // 2クラ
              if (i + j <= n) {
                if (jset.find(P(i, j)) == jset.end()) {
                  jset.insert(P(i, j));
                  numFourArr[cnt][0] = i;
                  numFourArr[cnt][1] = j;
                  numFourArr[cnt][2] = 0;
                  numFourArr[cnt][3] = 0;
                  cnt++;
                }
              }

              // 3クラ
              if (i + j + k <= n) {
                p.first.first = i;
                p.first.second = j;
                p.second.first = k;
                if (kset.find(p) == kset.end()) {
                  kset.insert(p);
                  numFourArr[cnt][0] = i;
                  numFourArr[cnt][1] = j;
                  numFourArr[cnt][2] = k;
                  numFourArr[cnt][3] = 0;
                  cnt++;
                }
              }

              // 4クラ
              if (i + j + k + l <= n) {
                numFourArr[cnt][0] = i;
                numFourArr[cnt][1] = j;
                numFourArr[cnt][2] = k;
                numFourArr[cnt][3] = l;
                cnt++;
              }
              i += hyperStep1;
              if (cnt >= 500) break;
            }
            j += hyperStep1;
            if (cnt >= 500) break;
          }
          k += hyperStep2;
          if (cnt >= 500) break;
        }
        l += hyperStep2;
        if (cnt >= 500) break;
      }
    }

    for (int i = 0; i < (cnt / 2); ++i) {
      swap(numFourArr[i][0], numFourArr[cnt - 1 - i][0]);
      swap(numFourArr[i][1], numFourArr[cnt - 1 - i][1]);
      swap(numFourArr[i][2], numFourArr[cnt - 1 - i][2]);
      swap(numFourArr[i][3], numFourArr[cnt - 1 - i][3]);
    }

    SetGraphFromFourArray();
  }

  // 13表裏
  int omoteArr[1000];
  void InitNumArray18()
  {
    for (int i = 0; i < (1000); ++i)omoteArr[i] = 0;
    numPairArrOK = 1;
    int cnt = 0;
    int one = 1;
    for (int i = n; i > 0; i -= hyperStep1) {
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        one = 0;
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        omoteArr[cnt] = 0;
        cnt++;
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        omoteArr[cnt] = 1;
        cnt++;
        j -= hyperStep2;
      }
      if (one) {
        numPairArr[cnt][0] = i;
        numPairArr[cnt][1] = 0;
        omoteArr[cnt] = 0;
        cnt++;
        numPairArr[cnt][0] = i;
        numPairArr[cnt][1] = 0;
        omoteArr[cnt] = 1;
        cnt++;
      }
      if (cnt > 200) break;
    }

    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      int omote = omoteArr[i];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 <= j && k < num1 + num2) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
          if (omote) {
            a[i][j][k] = 1 - a[i][j][k];
            a[i][k][j] = 1 - a[i][k][j];
          }
        }
      }
    }
    if (cnt < m) numPairArrOK = 0;
  }

  // 10表裏
  void InitNumArray19()
  {
    for (int i = 0; i < (1000); ++i)omoteArr[i] = 0;
    numPairArrOK = 1;
    int cnt = 0;
    for (int i = n; i > 0; i -= hyperStep1) {
      numPairArr[cnt][0] = i;
      numPairArr[cnt][1] = 0;
      omoteArr[cnt] = 0;
      cnt++;
      numPairArr[cnt][0] = i;
      numPairArr[cnt][1] = 0;
      omoteArr[cnt] = 1;
      cnt++;
      int j = n - i;
      while (j >= i + hyperMinDiff) {
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        omoteArr[cnt] = 0;
        cnt++;
        numPairArr[cnt][0] = j;
        numPairArr[cnt][1] = i;
        omoteArr[cnt] = 1;
        cnt++;
        j -= hyperStep1;
      }
      if (cnt > 200) break;
    }

    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      int omote = omoteArr[i];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else if (num1 <= j && k < num1 + num2) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
          if (omote) {
            a[i][j][k] = 1 - a[i][j][k];
            a[i][k][j] = 1 - a[i][k][j];
          }
        }
      }
    }
    if (cnt < m) numPairArrOK = 0;
  }

  // 0.0用
  vector<vector<P>> zeroPairs;
  vector<int> dfsvec;
  void dfs(int sum, int x)
  {
    if (dfsvec.size()) {

      if (dfsvec.size() == 1) {
        int num1 = dfsvec[0];
        for (int i = num1 - 1; i < num1 * (num1 - 1) / 2 + 1; ++i) {
          vector<P> vp;
          vp.push_back(P(num1, i));
          zeroPairs.push_back(vp);
        }
      }
      if (dfsvec.size() == 2) {
        int num1 = dfsvec[0];
        int num2 = dfsvec[1];
        for (int i = num1 - 1; i < num1 * (num1 - 1) / 2 + 1; ++i) {
          for (int j = num2 - 1; j < num2 * (num2 - 1) / 2 + 1; ++j) {
            vector<P> vp;
            vp.push_back(P(num1, i));
            vp.push_back(P(num2, j));
            zeroPairs.push_back(vp);
          }
        }
      }
    }

    for (int i = x + 1; i < 10; ++i) {
      if (sum + i <= n) {
        dfsvec.push_back(i);
        dfs(sum + i, i);
        dfsvec.pop_back();
      }
    }
  }
  void InitNumArray20()
  {
    zeroPairs.clear();
    numPairArrOK = 0;
    dfs(0, 1);
    if (zeroPairs.size() >= m) {
      numPairArrOK = 1;
      for (int i = 0; i < (m); ++i) {
        for (int j = 0; j < (n); ++j) {
          for (int k = 0; k < (n); ++k) {
            a[i][j][k] = 0;
          }
        }
      }
      for (int i = 0; i < (m); ++i) {
        vector<P> vp = zeroPairs[i];
        int sz = vp.size();
        int sum = 0;
        for (int j = 0; j < (sz); ++j) {
          int sq = vp[j].first;
          int hon = vp[j].second;
          int cnt = 0;
          for (int k = sum; k < sum + sq; ++k) {
            for (int l = k + 1; l < sum + sq; ++l) {
              if (cnt < hon) {
                a[i][k][l] = 1;
                a[i][l][k] = 1;
                cnt++;
              }
            }
          }
          sum += sq;
        }
      }
    }
  }

  // 1コア
  void InitNumArray21()
  {
    for (int i = 0; i < (1000); ++i)omoteArr[i] = 0;
    numPairArrOK = 1;
    int cnt = 0;
    for (int i = n; i > 0; i -= hyperStep1) {
      numSingleArr[cnt] = i;
      omoteArr[cnt] = 0;
      cnt++;
      numSingleArr[cnt] = i;
      omoteArr[cnt] = 1;
      cnt++;
      if (cnt > 200) break;
    }

    for (int i = 0; i < (m); ++i) {
      int num1 = numSingleArr[i];
      int omote = omoteArr[i];
      for (int j = 0; j < (n); ++j) {
        for (int k = j + 1; k < n; ++k) {
          if (k < num1) {
            a[i][j][k] = 1;
            a[i][k][j] = 1;
          }
          else {
            a[i][j][k] = 0;
            a[i][k][j] = 0;
          }
          if (omote) {
            a[i][j][k] = 1 - a[i][j][k];
            a[i][k][j] = 1 - a[i][k][j];
          }
        }
      }
    }
    if (cnt < m) numPairArrOK = 0;
  }


  void InitNumArray(int mode)
  {
    numPairArrOK = 1;
    int ra = hyperSolverNum % 1000 / 10;

    if (ra == 1) {
      InitNumArray1();
    }
    else if (ra == 2) {
      n = 100;
      InitNumArray2();
    }
    else if (ra == 3) {
      n = 100;
      InitNumArray3();
    }
    else if (ra == 4) {
      n = 100;
      InitNumArray4();
    }
    else if (ra == 5) {
      InitNumArray5();
    }
    else if (ra == 6) {
      InitNumArray6();
    }
    else if (ra == 7) {
      InitNumArray7();
    }
    else if (ra == 8) {
      InitNumArray8();
    }
    else if (ra == 9) {
      InitNumArray9();
    }
    else if (ra == 10) {
      InitNumArray10();
    }
    else if (ra == 11) {
      InitNumArray11();
    }
    else if (ra == 12) {
      InitNumArray12();
    }
    else if (ra == 13) {
      InitNumArray13();
    }
    else if (ra == 14) {
      InitNumArray14();
    }
    else if (ra == 15) {
      InitNumArray15();
    }
    else if (ra == 16) {
      InitNumArray16();
    }
    else if (ra == 17) {
      InitNumArray17();
    }
    else if (ra == 18) {
      InitNumArray18();
    }
    else if (ra == 19) {
      InitNumArray19();
    }
    else if (ra == 20) {
      InitNumArray20();
    }
    else if (ra == 21) {
      InitNumArray21();
    }
  }

}  // namespace

int Solver1()
{
  vector<int> keep[110];

  int cnt[100] = {};
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
  }

  for (int i = 0; i < (m); ++i) {
    int num = numArr[i];
    double e0 = ((double)n - 1) * eps;
    double e1 = ((double)num - 1) + ((double)n - num) * eps - ((double)num - 1) * eps;
    int kijun = round(e0 + (e1 - e0) / 2.0);
    int count = 0;
    for (int j = 0; j < (100); ++j) {
      if (cnt[j] >= kijun) {
        count++;
      }
    }
    keep[abs(num - count)].push_back(i);
  }

  int res = -1;
  for (int i = 0; i < (110); ++i) {
    if (keep[i].size()) {
      int sz = keep[i].size();
      res = keep[i][sz / 2];
      break;
    }
  }

  return res;
}

int Solver2()
{
  int cnt[100] = {};
  int f[100] = {};
  int res = n;
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
    f[i] = 1;
  }

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num = numArr[i];
    if (abs(num - res) < diff) {
      diff = abs(num - res);
      argRes = i;
    }
  }

  return argRes;
}

int Solver3()
{
  vector<int> vec[2];
  int f[100];
  for (int i = 0; i < (n); ++i) {
    if (i == 0) {
      vec[0].push_back(i);
      f[i] = 0;
    }
    else {
      if (b[0][i]) {
        vec[0].push_back(i);
        f[i] = 0;
      }
      else {
        vec[1].push_back(i);
        f[i] = 1;
      }
    }
  }
  for (int _ = 0; _ < (100); ++_) {
    if (vec[0].empty() || vec[1].empty()) break;
    double cnt[110][2];
    for (int i = 0; i < (n); ++i) for (int j = 0; j < (2); ++j) cnt[i][j] = 0;
    vector<int> nxt[2];
    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        if (i == j) {
          cnt[i][0]++;
        }
        else if (b[i][j]) {
          if (f[i] == f[j]) {
            cnt[i][0]++;
          }
          else {
            cnt[i][1]++;
          }
        }
      }
      int sz = vec[f[i]].size();
      if (cnt[i][0] / sz >= cnt[i][1] / ((double)n - sz)) {
        nxt[f[i]].push_back(i);
      }
      else {
        nxt[1 - f[i]].push_back(i);
      }
    }

    for (int i = 0; i < (2); ++i) vec[i] = nxt[i];
  }
  int res = min(vec[0].size(), vec[1].size());
  // cout << res << endl;
  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num = numArr[i];
    if (abs(num - res) < diff) {
      diff = abs(num - res);
      argRes = i;
    }
  }

  return argRes;
}

int Solver4()
{
  int cnt[100] = {};
  int f[100] = {};
  int res = n;
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
    f[i] = 1;
  }

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }

  if (res >= 20) {
    int res2 = 0;
    vector<int> vec;
    for (int i = 0; i < (n); ++i) {
      if (f[i]) {
        vec.push_back(i);
      }
    }
    int ff[100] = {};
    double kijun = (eps * eps + (1.0 - eps) * (1.0 - eps)) / 2.0;
    for (int i = 0; i < (n); ++i) {
      int tri[2] = {};
      for (int j = 0; j < (res); ++j) {
        int jj = vec[j];
        if (jj == i) continue;
        for (int k = j + 1; k < res; ++k) {
          int kk = vec[k];
          if (kk == i) continue;
          tri[1]++;
          if (b[i][jj] && b[i][kk]) {
            tri[0]++;
          }
        }
      }
      if ((double)tri[0] / tri[1] >= kijun) {
        ff[i] = 1;
        res2++;
      }
    }
    res = res2;
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num = numArr[i];
    if (abs(num - res) < diff) {
      diff = abs(num - res);
      argRes = i;
    }
  }

  return argRes;
}

int Solver5()
{
  int cnt[100] = {};
  int f[100] = {};
  int res = n;
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
    f[i] = 1;
  }

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }

  int res1 = res;
  res = n - res;
  for (int i = 0; i < (n); ++i) { f[i] = 1 - f[i]; }
  for (int i = 0; i < (n); ++i) { cnt[i] = 0; }
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      if (f[i] && f[j]) cnt[i] += b[i][j];
    }
  }
  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }
  int res2 = res;
  if (res2 <= hyperMaxRound) res2 = 0;
  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

int Solver6()
{
  int cnt[100] = {};
  int f[100] = {};
  int ff[100] = {};
  int res = n;
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
    f[i] = 1;
  }

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }

  int res1 = res;
  res = n - res;
  for (int i = 0; i < (n); ++i) { f[i] = 1 - f[i]; }
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 0) ff[i] = 1;
  }
  for (int i = 0; i < (n); ++i) { cnt[i] = 0; }
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      if (f[i] && f[j]) cnt[i] += b[i][j];
    }
  }
  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }
  int res2 = res;
  if (res2 <= hyperMaxRound) {
    res2 = 0;
    for (int i = 0; i < (n); ++i) f[i] = 0;
  }
  for (int i = 0; i < (n); ++i) {
    if (f[i]) ff[i] = 2;
  }

  for (int i = 0; i < (n); ++i) f[i] = ff[i];

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[3] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 3;
    while (ra == f[x]) {
      ra = Rand() % 3;
    }
    if (res2 = 0) {
      ra = 1 - f[x];
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

int Solver7()
{
  int cnt[100] = {};
  int f[100] = {};
  int ff[100] = {};
  int res = n;
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { cnt[i] += b[i][j]; }
    f[i] = 1;
  }

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }

  int res1 = res;
  res = n - res;
  for (int i = 0; i < (n); ++i) { f[i] = 1 - f[i]; }
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 0) ff[i] = 1;
  }
  for (int i = 0; i < (n); ++i) { cnt[i] = 0; }
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      if (f[i] && f[j]) cnt[i] += b[i][j];
    }
  }
  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) break;
    for (int i = 0; i < (n); ++i) {
      if (i == arg) continue;
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }
    res--;
    f[arg] = 0;
  }
  int res2 = res;
  if (res2 <= hyperMaxRound) {
    res2 = 0;
    for (int i = 0; i < (n); ++i) f[i] = 0;
  }
  for (int i = 0; i < (n); ++i) {
    if (f[i]) ff[i] = 2;
  }

  for (int i = 0; i < (n); ++i) f[i] = ff[i];

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[3] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 3;
    while (ra == f[x]) {
      ra = Rand() % 3;
    }
    if (res2 == 0) {
      ra = 1 - f[x];
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

int Solver8()
{
  int cnt = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) { cnt += b[i][j]; }
  }
  cnt = min(cnt, m - 1);
  return cnt;
}

int Solver9()
{
  map<P, int> mp;
  int fff[10][100];
  P kp[10];

  int kcnt[100] = {};
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { kcnt[i] += b[i][j]; }
  }

  for (int _ = 0; _ < (10); ++_) {
    int cnt[100] = {};
    int f[100] = {};
    int ff[100] = {};
    int res = n;

    for (int i = 0; i < (n); ++i) {
      cnt[i] = kcnt[i];
      f[i] = 1;
    }

    while (res > 1) {
      int mi = 1000;
      vector<int> arv;
      for (int i = 0; i < (n); ++i) {
        if (f[i] && cnt[i] <= mi && cnt[i] < (res + 1) / 2) {
          if (cnt[i] == mi) {
            arv.push_back(i);
          }
          else {
            arv.clear();
            arv.push_back(i);
          }
          mi = cnt[i];
        }
      }
      if (arv.empty()) break;
      int arg = arv[Rand() % arv.size()];
      for (int i = 0; i < (n); ++i) {
        if (i == arg) continue;
        if (f[i] && b[i][arg]) {
          cnt[i]--;
        }
      }
      res--;
      f[arg] = 0;
    }

    int res1 = res;
    res = n - res;
    for (int i = 0; i < (n); ++i) { f[i] = 1 - f[i]; }
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 0) ff[i] = 1;
    }
    for (int i = 0; i < (n); ++i) { cnt[i] = 0; }
    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        if (f[i] && f[j]) cnt[i] += b[i][j];
      }
    }
    while (res > 1) {
      int mi = 1000;
      vector<int> arv;
      for (int i = 0; i < (n); ++i) {
        if (f[i] && cnt[i] <= mi && cnt[i] < (res + 1) / 2) {
          if (cnt[i] == mi) {
            arv.push_back(i);
          }
          else {
            arv.clear();
            arv.push_back(i);
          }
          mi = cnt[i];
        }
      }
      if (arv.empty()) break;
      int arg = arv[Rand() % arv.size()];
      for (int i = 0; i < (n); ++i) {
        if (i == arg) continue;
        if (f[i] && b[i][arg]) {
          cnt[i]--;
        }
      }
      res--;
      f[arg] = 0;
    }
    int res2 = res;
    if (res2 <= hyperMaxRound) {
      res2 = 0;
      for (int i = 0; i < (n); ++i) f[i] = 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (f[i]) ff[i] = 2;
    }

    for (int i = 0; i < (n); ++i) f[i] = ff[i];

    for (int i = 0; i < (n); ++i) { fff[_][i] = f[i]; }
    mp[P(res1, res2)]++;
    kp[_] = P(res1, res2);
  }

  int ma = 0;
  P maxP;
  for (auto elem : mp) {
    if (elem.second > ma) {
      ma = elem.second;
      maxP = elem.first;
    }
  }

  int f[100] = {};
  int res1, res2;
  res1 = maxP.first;
  res2 = maxP.second;
  for (int i = 0; i < (10); ++i) {
    if (maxP == kp[i]) {
      for (int j = 0; j < (100); ++j) { f[j] = fff[i][j]; }
    }
  }

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[3] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 3;
    while (ra == f[x]) {
      ra = Rand() % 3;
    }
    if (res2 == 0) {
      ra = 1 - f[x];
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

int Solver10()
{
  map<P, int> mp;
  int fff[31][100];
  P kp[31];

  int kcnt[100] = {};
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { kcnt[i] += b[i][j]; }
  }

  for (int _ = 0; _ < (31); ++_) {
    int cnt[100] = {};
    int f[100] = {};
    int ff[100] = {};
    int res = n;

    for (int i = 0; i < (n); ++i) {
      cnt[i] = kcnt[i];
      f[i] = 1;
    }

    while (res > 1) {
      int mi = 1000;
      vector<int> arv;
      for (int i = 0; i < (n); ++i) {
        if (f[i] && cnt[i] <= mi && cnt[i] < (res + 1) / 2) {
          if (cnt[i] == mi) {
            arv.push_back(i);
          }
          else {
            arv.clear();
            arv.push_back(i);
          }
          mi = cnt[i];
        }
      }
      if (arv.empty()) break;
      int arg = arv[Rand() % arv.size()];
      for (int i = 0; i < (n); ++i) {
        if (i == arg) continue;
        if (f[i] && b[i][arg]) {
          cnt[i]--;
        }
      }
      res--;
      f[arg] = 0;
    }

    int res1 = res;
    res = n - res;
    for (int i = 0; i < (n); ++i) { f[i] = 1 - f[i]; }
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 0) ff[i] = 1;
    }
    for (int i = 0; i < (n); ++i) { cnt[i] = 0; }
    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        if (f[i] && f[j]) cnt[i] += b[i][j];
      }
    }
    while (res > 1) {
      int mi = 1000;
      vector<int> arv;
      for (int i = 0; i < (n); ++i) {
        if (f[i] && cnt[i] <= mi && cnt[i] < (res + 1) / 2) {
          if (cnt[i] == mi) {
            arv.push_back(i);
          }
          else {
            arv.clear();
            arv.push_back(i);
          }
          mi = cnt[i];
        }
      }
      if (arv.empty()) break;
      int arg = arv[Rand() % arv.size()];
      for (int i = 0; i < (n); ++i) {
        if (i == arg) continue;
        if (f[i] && b[i][arg]) {
          cnt[i]--;
        }
      }
      res--;
      f[arg] = 0;
    }
    int res2 = res;
    if (res2 <= hyperMaxRound) {
      res2 = 0;
      for (int i = 0; i < (n); ++i) f[i] = 0;
    }
    for (int i = 0; i < (n); ++i) {
      if (f[i]) ff[i] = 2;
    }

    for (int i = 0; i < (n); ++i) f[i] = ff[i];

    for (int i = 0; i < (n); ++i) { fff[_][i] = f[i]; }
    mp[P(res1, res2)]++;
    kp[_] = P(res1, res2);
  }

  int ma = 0;
  P maxP;
  for (auto elem : mp) {
    if (elem.second > ma) {
      ma = elem.second;
      maxP = elem.first;
    }
  }

  int f[100] = {};
  int res1, res2;
  res1 = maxP.first;
  res2 = maxP.second;
  for (int i = 0; i < (10); ++i) {
    if (maxP == kp[i]) {
      for (int j = 0; j < (100); ++j) { f[j] = fff[i][j]; }
    }
  }

  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

// 4-クリークを見つける共通関数
bool findClique4(const vector<int>& kouho, int f[], vector<int>& cores, int markValue)
{
  for (int loop1 = 0; loop1 < (5000); ++loop1) {
    int core[4] = {};
    for (int i = 0; i < (4); ++i) {
      while (true) {
        core[i] = kouho[Rand() % kouho.size()];
        for (int j = 0; j < (i); ++j) {
          if (core[j] == core[i]) core[i] = -1;
        }
        if (core[i] != -1) break;
      }
    }
    int mitu = 1;
    for (int i = 0; i < (4); ++i) {
      for (int j = i + 1; j < 4; ++j) {
        if (!b[core[i]][core[j]]) mitu = 0;
      }
    }
    if (mitu) {
      for (int i = 0; i < (4); ++i) {
        f[core[i]] = markValue;
        cores.push_back(core[i]);
      }
      return true;
    }
  }
  return false;
}

int Solver11()
{
  int f[110] = {};

  // コア1を作る
  vector<int> cores1;
  vector<int> kouho;
  for (int i = 0; i < (n); ++i) kouho.push_back(i);
  if (kouho.size() < 4) return 0;
  if (!findClique4(kouho, f, cores1, 1)) return 0;

  // コア1を大きくしていく
  while (true) {
    int sz = cores1.size();
    int arg = -1;
    int ma = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (sz); ++j) {
        if (b[i][cores1[j]]) cnt++;
      }
      if (ma < cnt && cnt >= (sz + 2) / 2) {
        ma = cnt;
        arg = i;
      }
    }
    if (arg == -1) break;
    f[arg] = 1;
    cores1.push_back(arg);
  }

  // コア2を作る
  vector<int> cores2;
  kouho.clear();
  for (int i = 0; i < (n); ++i) {
    if (f[i] != 0) continue;
    int cnt = 0;
    for (int j = 0; j < (n); ++j) {
      if (f[j] == 0) cnt += b[i][j];
    }

    if (cnt <= 4) continue;
    kouho.push_back(i);
  }
  if (kouho.size() >= 4) {
    findClique4(kouho, f, cores2, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      while (true) {
        int sz = cores2.size();
        int arg = -1;
        int ma = -1;
        for (int i = 0; i < (n); ++i) {
          if (f[i] != 0) continue;
          int cnt = 0;
          for (int j = 0; j < (sz); ++j) {
            if (b[i][cores2[j]]) cnt++;
          }
          if (ma < cnt && cnt >= (sz + 2) / 2) {
            ma = cnt;
            arg = i;
          }
        }
        if (arg == -1) break;
        f[arg] = 2;
        cores2.push_back(arg);
      }
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[3] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 3;
    while (ra == f[x]) {
      ra = Rand() % 3;
    }
    if (res2 == 0) {
      ra = 1 - f[x];
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
  }
  if (res2 > res1) swap(res1, res2);
  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }

  return argRes;
}

int Solver12()
{
  int f[110] = {};

  // コア1を作る
  vector<int> cores1;
  vector<int> kouho;
  for (int i = 0; i < (n); ++i) kouho.push_back(i);
  if (kouho.size() < 4) return 0;
  findClique4(kouho, f, cores1, 1);
  if (cores1.size() == 0) return 0;

  // コア1を大きくしていく
  while (true) {
    int sz = cores1.size();
    int arg = -1;
    int ma = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (sz); ++j) {
        if (b[i][cores1[j]]) cnt++;
      }
      if (ma < cnt && cnt >= (sz + 2) / 2) {
        ma = cnt;
        arg = i;
      }
    }
    if (arg == -1) break;
    f[arg] = 1;
    cores1.push_back(arg);
  }

  // コア2を作る
  vector<int> cores2;
  kouho.clear();
  for (int i = 0; i < (n); ++i) {
    if (f[i] != 0) continue;
    int cnt = 0;
    for (int j = 0; j < (n); ++j) {
      if (f[j] == 0) cnt += b[i][j];
    }

    if (cnt <= 4) continue;
    kouho.push_back(i);
  }
  if (kouho.size() >= 4) {
    findClique4(kouho, f, cores2, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      while (true) {
        int sz = cores2.size();
        int arg = -1;
        int ma = -1;
        for (int i = 0; i < (n); ++i) {
          if (f[i] != 0) continue;
          int cnt = 0;
          for (int j = 0; j < (sz); ++j) {
            if (b[i][cores2[j]]) cnt++;
          }
          if (ma < cnt && cnt >= (sz + 2) / 2) {
            ma = cnt;
            arg = i;
          }
        }
        if (arg == -1) break;
        f[arg] = 2;
        cores2.push_back(arg);
      }
    }
  }

  // コア3を作る
  vector<int> cores3;
  if (cores2.size() > 0) {
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores3, 3);
      if (cores3.size() > 0) {
        // コア3を大きくしていく
        while (true) {
          int sz = cores3.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores3[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 3;
          cores3.push_back(arg);
        }
      }
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();
  int res3 = cores3.size();

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[4] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 4;
    while (ra == f[x]) {
      ra = Rand() % 4;
    }
    if (res2 == 0) {
      ra = 1 - f[x];
    }
    else if (res3 == 0) {
      ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  res3 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
    if (f[i] == 3) res3++;
  }
  vector<int> resv;
  // cout << res1 << ' ' << res2 << ' ' << res3 << endl;
  resv.push_back(res1);
  resv.push_back(res2);
  resv.push_back(res3);
  sort(resv.begin(), resv.end());
  res1 = resv[2];
  res2 = resv[1];
  res3 = resv[0];
  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numThreeArr[i][0];
    int num2 = numThreeArr[i][1];
    int num3 = numThreeArr[i][2];
    if (abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3);
      argRes = i;
    }
  }

  // if (argRes != judgeNum) {
  //   for (int i = 0; i < (n); ++i) cout << f[i];
  //   cout << "   " << diff << endl;
  //   cout << endl;
  // }

  return argRes;
}

// 4コア
int Solver13()
{
  int f[110] = {};

  // コア1を作る
  vector<int> cores1;
  vector<int> kouho;
  for (int i = 0; i < (n); ++i) kouho.push_back(i);
  if (kouho.size() < 4) return 0;
  findClique4(kouho, f, cores1, 1);
  if (cores1.size() == 0) return 0;

  // コア1を大きくしていく
  while (true) {
    int sz = cores1.size();
    int arg = -1;
    int ma = -1;
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (sz); ++j) {
        if (b[i][cores1[j]]) cnt++;
      }
      if (ma < cnt && cnt >= (sz + 2) / 2) {
        ma = cnt;
        arg = i;
      }
    }
    if (arg == -1) break;
    f[arg] = 1;
    cores1.push_back(arg);
  }

  // コア2を作る
  vector<int> cores2;
  kouho.clear();
  for (int i = 0; i < (n); ++i) {
    if (f[i] != 0) continue;
    int cnt = 0;
    for (int j = 0; j < (n); ++j) {
      if (f[j] == 0) cnt += b[i][j];
    }

    if (cnt <= 4) continue;
    kouho.push_back(i);
  }
  if (kouho.size() >= 4) {
    findClique4(kouho, f, cores2, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      while (true) {
        int sz = cores2.size();
        int arg = -1;
        int ma = -1;
        for (int i = 0; i < (n); ++i) {
          if (f[i] != 0) continue;
          int cnt = 0;
          for (int j = 0; j < (sz); ++j) {
            if (b[i][cores2[j]]) cnt++;
          }
          if (ma < cnt && cnt >= (sz + 2) / 2) {
            ma = cnt;
            arg = i;
          }
        }
        if (arg == -1) break;
        f[arg] = 2;
        cores2.push_back(arg);
      }
    }
  }

  // コア3を作る
  vector<int> cores3;
  if (cores2.size() > 0) {
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores3, 3);
      if (cores3.size() > 0) {
        // コア3を大きくしていく
        while (true) {
          int sz = cores3.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores3[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 3;
          cores3.push_back(arg);
        }
      }
    }
  }

  // コア4を作る
  vector<int> cores4;
  if (cores3.size() > 0) {
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores4, 4);
      if (cores4.size() > 0) {
        // コア4を大きくしていく
        while (true) {
          int sz = cores4.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores4[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 4;
          cores4.push_back(arg);
        }
      }
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();
  int res3 = cores3.size();
  int res4 = cores4.size();

  int score = 0;
  for (int i = 0; i < (n); ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (f[i] == 0 && f[j] == 0) {
        score += 1 - b[i][j];
      }
      else {
        if ((f[i] == f[j]) == (b[i][j])) {
          score++;
        }
      }
    }
  }

  bitset<100> bif[5] = {}, bib[100] = {};
  for (int i = 0; i < (n); ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
  }

  bitset<100> bione(0);
  for (int i = 0; i < (n); ++i) { bione[i] = 1; }

  int flipLoop = 1000;
  if (MODE == 0) flipLoop = 10000;
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = Rand() % n;
    int ra = Rand() % 5;
    while (ra == f[x]) {
      ra = Rand() % 5;
    }
    if (res2 == 0) {
      ra = 1 - f[x];
    }
    else if (res3 == 0) {
      ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
    }
    else if (res4 == 0) {
      ra = Rand() % 4;
      while (ra == f[x]) {
        ra = Rand() % 4;
      }
    }
    int tmp = score;
    int keep = f[x];

    // tmpからxの点を引く
    if (f[x] == 0) {
      tmp -= n - bib[x].count();
      tmp++;  // 自分の分
    }
    else {
      tmp -= (bif[f[x]] & bib[x]).count();
      tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    bif[f[x]][x] = 0;
    f[x] = ra;
    bif[f[x]][x] = 1;
    // tmpからraの点を足す
    if (f[x] == 0) {
      tmp += n - bib[x].count();
      tmp--;  // 自分の分
    }
    else {
      tmp += (bif[f[x]] & bib[x]).count();
      tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
    }

    if (tmp >= score) {
      score = tmp;
    }
    else {
      bif[f[x]][x] = 0;
      f[x] = keep;
      bif[f[x]][x] = 1;
    }
  }
  res1 = 0;
  res2 = 0;
  res3 = 0;
  res4 = 0;
  for (int i = 0; i < (n); ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
    if (f[i] == 3) res3++;
    if (f[i] == 4) res4++;
  }
  vector<int> resv;
  // cout << res1 << ' ' << res2 << ' ' << res3 << endl;
  resv.push_back(res1);
  resv.push_back(res2);
  resv.push_back(res3);
  resv.push_back(res4);
  sort(resv.begin(), resv.end());
  res1 = resv[3];
  res2 = resv[2];
  res3 = resv[1];
  res4 = resv[0];
  int diff = 1000;
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    int num1 = numFourArr[i][0];
    int num2 = numFourArr[i][1];
    int num3 = numFourArr[i][2];
    int num4 = numFourArr[i][3];
    if (abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) +
      abs(num4 - res4) <
      diff) {
      diff = abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) +
        abs(num4 - res4);
      argRes = i;
    }
  }

  return argRes;
}

int Solver14()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;
    ;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores2, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    if (diff < real_minDiff) {
      real_minDiff = diff;
      real_argRes = argRes;
      if (real_minDiff == 0) break;
    }
  }

  return real_argRes;
}

int Solver15()
{
  map<int, int> argMap;

  for (int wataruoop = 0; wataruoop < (5); ++wataruoop) {
    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;
    ;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores2, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    argMap[argRes]++;

    // if (argRes != judgeNum) {
    //   for (int i = 0; i < (n); ++i) cout << f[i];
    //   cout << endl;
    //   for (int i = 0; i < (n); ++i) {
    //     if (i < numPairArr[judgeNum][0])
    //       cout << 1;
    //     else if (numPairArr[judgeNum][0] <= i &&
    //              i <= numPairArr[judgeNum][0] + numPairArr[judgeNum][1])
    //       cout << 2;
    //     else
    //       cout << 0;
    //   }
    //   cout << "   " << diff << endl;
    // }
  }

  int maxCnt = 0;
  int real_argRes = 0;
  for (auto p : argMap) {
    if (p.second > maxCnt) {
      maxCnt = p.second;
      real_argRes = p.first;
    }
  }

  return real_argRes;
}

// 11表裏
int Solver16()
{
  int real_argRes = 0;
  int real_score = 0;
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }
  for (int tei = 0; tei < (2); ++tei) {

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (tei % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }


    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores2, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != tei % 2) continue;
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == f[j] && f[i] != 0) {
          if (b[i][j]) tmpScore++;
        }
        else {
          if (!b[i][j]) tmpScore++;
        }
      }
    }
    if (tmpScore > real_score) {
      real_score = tmpScore;
      real_argRes = argRes;
    }
  }

  return real_argRes;
}

// 11表裏5セット
int Solver17()
{
  int real_argRes = 0;
  int real_score = 0;
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }
  for (int tei = 0; tei < (10); ++tei) {

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (tei % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }


    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores2, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != tei % 2) continue;
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == f[j] && f[i] != 0) {
          if (b[i][j]) tmpScore++;
        }
        else {
          if (!b[i][j]) tmpScore++;
        }
      }
    }
    if (tmpScore > real_score) {
      real_score = tmpScore;
      real_argRes = argRes;
    }
  }

  return real_argRes;
}

// 0.0用
int Solver18()
{
  int visit[110] = {};
  int bb[110][110];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      bb[i][j] = b[i][j];
    }
  }

  vector<P> vp;
  for (int i = 0; i < (n); ++i) {
    if (visit[i]) {
      continue;
    }
    visit[i] = 1;
    int cnt1 = 1;
    int cnt2 = 0;
    queue<int> que;
    que.push(i);
    while (que.size()) {
      int x = que.front();
      que.pop();
      for (int j = 0; j < (n); ++j) {
        if (bb[x][j]) {
          cnt2++;
          bb[x][j] = 0;
          if (!visit[j]) {
            cnt1++;
            visit[j] = 1;
            que.push(j);
          }
        }
      }
    }
    cnt2 /= 2;
    if (cnt1 == 1) { continue; }
    vp.push_back(P(cnt1, cnt2));
  }
  sort(vp.begin(), vp.end());
  int argRes = 0;
  for (int i = 0; i < (m); ++i) {
    if (zeroPairs[i] == vp) {
      argRes = i;
      break;
    }
  }
  return argRes;
}

// 14表裏
int Solver19()
{
  int real_argRes[2] = {};
  int real_minDiff[2] = { 1000,1000 };
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }


  for (int wataruoop = 0; wataruoop < (10); ++wataruoop) {
    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (wataruoop % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }

    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;
    ;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 4) {
      findClique4(kouho, f, cores2, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != wataruoop % 2) continue;
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    if (diff < real_minDiff[wataruoop % 2]) {
      real_minDiff[wataruoop % 2] = diff;
      real_argRes[wataruoop % 2] = argRes;
    }
  }

  int real_real_argRes = real_argRes[0];
  if (real_minDiff[1] < real_minDiff[0]) {
    real_real_argRes = real_argRes[1];
  }
  return real_real_argRes;
}

// 21表裏5セット
int Solver20()
{
  int real_argRes = 0;
  int real_score = 0;
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }
  for (int tei = 0; tei < (2); ++tei) {

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (tei % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }


    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;;
    findClique4(kouho, f, cores1, 1);
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    int res1 = cores1.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 2;
      while (ra == f[x]) {
        ra = Rand() % 2;
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
    }
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != tei % 2) continue;
      int num1 = numSingleArr[i];
      if (abs(num1 - res1) < diff) {
        diff = abs(num1 - res1);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == f[j] && f[i] != 0) {
          if (b[i][j]) tmpScore++;
        }
        else {
          if (!b[i][j]) tmpScore++;
        }
      }
    }
    if (tmpScore > real_score) {
      real_score = tmpScore;
      real_argRes = argRes;
    }
  }

  return real_argRes;
}

// 14亜種
int Solver21()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 4) continue;
    ;
    for (int loop1 = 0; loop1 < (5000); ++loop1) {
      int core[5] = {};
      for (int i = 0; i < (5); ++i) {
        while (true) {
          core[i] = kouho[Rand() % kouho.size()];
          for (int j = 0; j < (i); ++j) {
            if (core[j] == core[i]) core[i] = -1;
          }
          if (core[i] != -1) break;
        }
      }
      int mitu = 1;
      for (int i = 0; i < (5); ++i) {
        for (int j = i + 1; j < 5; ++j) {
          if (!b[core[i]][core[j]]) mitu = 0;
        }
      }
      if (mitu) {
        for (int i = 0; i < (5); ++i) {
          f[core[i]] = 1;
          cores1.push_back(core[i]);
        }
        break;
      }
    }
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 5) {
      for (int loop1 = 0; loop1 < (5000); ++loop1) {
        int core[5] = {};
        for (int i = 0; i < (5); ++i) {
          while (true) {
            core[i] = kouho[Rand() % kouho.size()];
            for (int j = 0; j < (i); ++j) {
              if (core[j] == core[i]) core[i] = -1;
            }
            if (core[i] != -1) break;
          }
        }
        int mitu = 1;
        for (int i = 0; i < (5); ++i) {
          for (int j = i + 1; j < 5; ++j) {
            if (!b[core[i]][core[j]]) mitu = 0;
          }
        }
        if (mitu) {
          for (int i = 0; i < (5); ++i) {
            f[core[i]] = 2;
            cores2.push_back(core[i]);
          }
          break;
        }
      }
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    if (diff < real_minDiff) {
      real_minDiff = diff;
      real_argRes = argRes;
      if (real_minDiff == 0) break;
    }
  }

  return real_argRes;
}

// 14亜種
int Solver22()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 3) continue;
    ;
    for (int loop1 = 0; loop1 < (5000); ++loop1) {
      int core[5] = {};
      for (int i = 0; i < (3); ++i) {
        while (true) {
          core[i] = kouho[Rand() % kouho.size()];
          for (int j = 0; j < (i); ++j) {
            if (core[j] == core[i]) core[i] = -1;
          }
          if (core[i] != -1) break;
        }
      }
      int mitu = 1;
      for (int i = 0; i < (3); ++i) {
        for (int j = i + 1; j < 3; ++j) {
          if (!b[core[i]][core[j]]) mitu = 0;
        }
      }
      if (mitu) {
        for (int i = 0; i < (3); ++i) {
          f[core[i]] = 1;
          cores1.push_back(core[i]);
        }
        break;
      }
    }
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 2) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 3) {
      for (int loop1 = 0; loop1 < (5000); ++loop1) {
        int core[5] = {};
        for (int i = 0; i < (3); ++i) {
          while (true) {
            core[i] = kouho[Rand() % kouho.size()];
            for (int j = 0; j < (i); ++j) {
              if (core[j] == core[i]) core[i] = -1;
            }
            if (core[i] != -1) break;
          }
        }
        int mitu = 1;
        for (int i = 0; i < (3); ++i) {
          for (int j = i + 1; j < 3; ++j) {
            if (!b[core[i]][core[j]]) mitu = 0;
          }
        }
        if (mitu) {
          for (int i = 0; i < (3); ++i) {
            f[core[i]] = 2;
            cores2.push_back(core[i]);
          }
          break;
        }
      }
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    if (diff < real_minDiff) {
      real_minDiff = diff;
      real_argRes = argRes;
      if (real_minDiff == 0) break;
    }
  }

  return real_argRes;
}

// 16亜種
int Solver23()
{
  int real_argRes = 0;
  int real_score = 0;
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }
  for (int tei = 0; tei < (2); ++tei) {

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (tei % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }


    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 3) continue;;
    for (int loop1 = 0; loop1 < (5000); ++loop1) {
      int core[4] = {};
      for (int i = 0; i < (3); ++i) {
        while (true) {
          core[i] = kouho[Rand() % kouho.size()];
          for (int j = 0; j < (i); ++j) {
            if (core[j] == core[i]) core[i] = -1;
          }
          if (core[i] != -1) break;
        }
      }
      int mitu = 1;
      for (int i = 0; i < (3); ++i) {
        for (int j = i + 1; j < 3; ++j) {
          if (!b[core[i]][core[j]]) mitu = 0;
        }
      }
      if (mitu) {
        for (int i = 0; i < (3); ++i) {
          f[core[i]] = 1;
          cores1.push_back(core[i]);
        }
        break;
      }
    }
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 2) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 3) {
      for (int loop1 = 0; loop1 < (5000); ++loop1) {
        int core[3] = {};
        for (int i = 0; i < (3); ++i) {
          while (true) {
            core[i] = kouho[Rand() % kouho.size()];
            for (int j = 0; j < (i); ++j) {
              if (core[j] == core[i]) core[i] = -1;
            }
            if (core[i] != -1) break;
          }
        }
        int mitu = 1;
        for (int i = 0; i < (3); ++i) {
          for (int j = i + 1; j < 3; ++j) {
            if (!b[core[i]][core[j]]) mitu = 0;
          }
        }
        if (mitu) {
          for (int i = 0; i < (3); ++i) {
            f[core[i]] = 2;
            cores2.push_back(core[i]);
          }
          break;
        }
      }
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != tei % 2) continue;
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == f[j] && f[i] != 0) {
          if (b[i][j]) tmpScore++;
        }
        else {
          if (!b[i][j]) tmpScore++;
        }
      }
    }
    if (tmpScore > real_score) {
      real_score = tmpScore;
      real_argRes = argRes;
    }
  }

  return real_argRes;
}

// 16亜種
int Solver24()
{
  int real_argRes = 0;
  int real_score = 0;
  int keepB[100][100];
  for (int i = 0; i < (n); ++i) {
    for (int j = 0; j < (n); ++j) {
      keepB[i][j] = b[i][j];
    }
  }
  for (int tei = 0; tei < (2); ++tei) {

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) {
        b[i][j] = keepB[i][j];
      }
    }
    if (tei % 2 == 1) {
      for (int i = 0; i < (n); ++i) {
        for (int j = 0; j < (n); ++j) {
          if (i == j) continue;
          b[i][j] = 1 - b[i][j];
        }
      }
    }


    int f[110] = {};

    // コア1を作る
    vector<int> cores1;
    vector<int> kouho;
    for (int i = 0; i < (n); ++i) kouho.push_back(i);
    if (kouho.size() < 5) continue;;
    for (int loop1 = 0; loop1 < (5000); ++loop1) {
      int core[5] = {};
      for (int i = 0; i < (5); ++i) {
        while (true) {
          core[i] = kouho[Rand() % kouho.size()];
          for (int j = 0; j < (i); ++j) {
            if (core[j] == core[i]) core[i] = -1;
          }
          if (core[i] != -1) break;
        }
      }
      int mitu = 1;
      for (int i = 0; i < (5); ++i) {
        for (int j = i + 1; j < 5; ++j) {
          if (!b[core[i]][core[j]]) mitu = 0;
        }
      }
      if (mitu) {
        for (int i = 0; i < (5); ++i) {
          f[core[i]] = 1;
          cores1.push_back(core[i]);
        }
        break;
      }
    }
    if (cores1.size() == 0) continue;

    // コア1を大きくしていく
    while (true) {
      int sz = cores1.size();
      int arg = -1;
      int ma = -1;
      for (int i = 0; i < (n); ++i) {
        if (f[i] != 0) continue;
        int cnt = 0;
        for (int j = 0; j < (sz); ++j) {
          if (b[i][cores1[j]]) cnt++;
        }
        if (ma < cnt && cnt >= (sz + 2) / 2) {
          ma = cnt;
          arg = i;
        }
      }
      if (arg == -1) break;
      f[arg] = 1;
      cores1.push_back(arg);
    }

    // コア2を作る
    vector<int> cores2;
    kouho.clear();
    for (int i = 0; i < (n); ++i) {
      if (f[i] != 0) continue;
      int cnt = 0;
      for (int j = 0; j < (n); ++j) {
        if (f[j] == 0) cnt += b[i][j];
      }

      if (cnt <= 4) continue;
      kouho.push_back(i);
    }
    if (kouho.size() >= 5) {
      for (int loop1 = 0; loop1 < (5000); ++loop1) {
        int core[5] = {};
        for (int i = 0; i < (5); ++i) {
          while (true) {
            core[i] = kouho[Rand() % kouho.size()];
            for (int j = 0; j < (i); ++j) {
              if (core[j] == core[i]) core[i] = -1;
            }
            if (core[i] != -1) break;
          }
        }
        int mitu = 1;
        for (int i = 0; i < (5); ++i) {
          for (int j = i + 1; j < 5; ++j) {
            if (!b[core[i]][core[j]]) mitu = 0;
          }
        }
        if (mitu) {
          for (int i = 0; i < (5); ++i) {
            f[core[i]] = 2;
            cores2.push_back(core[i]);
          }
          break;
        }
      }
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        while (true) {
          int sz = cores2.size();
          int arg = -1;
          int ma = -1;
          for (int i = 0; i < (n); ++i) {
            if (f[i] != 0) continue;
            int cnt = 0;
            for (int j = 0; j < (sz); ++j) {
              if (b[i][cores2[j]]) cnt++;
            }
            if (ma < cnt && cnt >= (sz + 2) / 2) {
              ma = cnt;
              arg = i;
            }
          }
          if (arg == -1) break;
          f[arg] = 2;
          cores2.push_back(arg);
        }
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == 0 && f[j] == 0) {
          score += 1 - b[i][j];
        }
        else {
          if ((f[i] == f[j]) == (b[i][j])) {
            score++;
          }
        }
      }
    }

    bitset<100> bif[3] = {}, bib[100] = {};
    for (int i = 0; i < (n); ++i) {
      if (f[i] > 0) {
        bif[f[i]][i] = 1;
      }
    }

    for (int i = 0; i < (n); ++i) {
      for (int j = 0; j < (n); ++j) { bib[i][j] = b[i][j]; }
    }

    bitset<100> bione(0);
    for (int i = 0; i < (n); ++i) { bione[i] = 1; }

    int flipLoop = 1000;
    if (MODE == 0) flipLoop = 10000;
    for (int _ = 0; _ < (flipLoop); ++_) {
      int x = Rand() % n;
      int ra = Rand() % 3;
      while (ra == f[x]) {
        ra = Rand() % 3;
      }
      if (res2 == 0) {
        ra = 1 - f[x];
      }
      int tmp = score;
      int keep = f[x];

      // tmpからxの点を引く
      if (f[x] == 0) {
        tmp -= n - bib[x].count();
        tmp++;  // 自分の分
      }
      else {
        tmp -= (bif[f[x]] & bib[x]).count();
        tmp -= ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      bif[f[x]][x] = 0;
      f[x] = ra;
      bif[f[x]][x] = 1;
      // tmpからraの点を足す
      if (f[x] == 0) {
        tmp += n - bib[x].count();
        tmp--;  // 自分の分
      }
      else {
        tmp += (bif[f[x]] & bib[x]).count();
        tmp += ((bif[f[x]] ^ bione) & (bib[x] ^ bione)).count();
      }

      if (tmp >= score) {
        score = tmp;
      }
      else {
        bif[f[x]][x] = 0;
        f[x] = keep;
        bif[f[x]][x] = 1;
      }
    }
    res1 = 0;
    res2 = 0;
    for (int i = 0; i < (n); ++i) {
      if (f[i] == 1) res1++;
      if (f[i] == 2) res2++;
    }
    if (res2 > res1) swap(res1, res2);
    int diff = 1000;
    int argRes = 0;
    for (int i = 0; i < (m); ++i) {
      if (omoteArr[i] != tei % 2) continue;
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < (n); ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (f[i] == f[j] && f[i] != 0) {
          if (b[i][j]) tmpScore++;
        }
        else {
          if (!b[i][j]) tmpScore++;
        }
      }
    }
    if (tmpScore > real_score) {
      real_score = tmpScore;
      real_argRes = argRes;
    }
  }

  return real_argRes;
}

int ComputeAnswer(int mode, int turn = 0)
{
  int res = 0;
  int num = num = hyperSolverNum % 10 + hyperSolverNum / 1000 * 10;

  if (num == 1) {
    res = Solver1();
  }
  else if (num == 2) {
    res = Solver2();
  }
  else if (num == 3) {
    res = Solver3();
  }
  else if (num == 4) {
    res = Solver4();
  }
  else if (num == 5) {
    res = Solver5();
  }
  else if (num == 6) {
    res = Solver6();
  }
  else if (num == 7) {
    res = Solver7();
  }
  else if (num == 8) {
    res = Solver8();
  }
  else if (num == 9) {
    res = Solver9();
  }
  else if (num == 10) {
    res = Solver10();
  }
  else if (num == 11) {
    res = Solver11();
  }
  else if (num == 12) {
    res = Solver12();
  }
  else if (num == 13) {
    res = Solver13();
  }
  else if (num == 14) {
    res = Solver14();
  }
  else if (num == 15) {
    res = Solver15();
  }
  else if (num == 16) {
    res = Solver16();
  }
  else if (num == 17) {
    res = Solver17();
  }
  else if (num == 18) {
    res = Solver18();
  }
  else if (num == 19) {
    res = Solver19();
  }
  else if (num == 20) {
    res = Solver20();
  }
  else if (num == 21) {
    res = Solver21();
  }
  else if (num == 22) {
    res = Solver22();
  }
  else if (num == 23) {
    res = Solver23();
  }
  else if (num == 24) {
    res = Solver24();
  }

  return res;
}

double Simulate(int mode)
{
  double res = 1e9 / n;
  for (int turn = 0; turn < (100); ++turn) {
    InitB(mode, turn);
    int ans = ComputeAnswer(mode, turn);
    answersFor1000Out[turn] = ans;
    int judge = judgeArr[turn];
    if (ans != judge) res *= 0.9;

    if (mode == 0) {
      cout << ans << endl;
      fflush(stdout);
    }
  }
  return res;
}

void solve(int mode)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // 提出用
  if (mode == 0) {
    Input(mode);

    n = hyperN[m][iEps];
    hyperSolverNum = hyperSolver[m][iEps];
    hyperMinDiff = hyperMinDiffArr[m][iEps];
    hyperMaxRound = hyperMaxRoundArr[m][iEps];
    hyperStep1 = hyperStep1Arr[m][iEps];
    hyperStep2 = hyperStep2Arr[m][iEps];

    InitNumArray(mode);

    OutputArrayAsString(mode);

    Simulate(mode);
  }

  // ハイパラ調整
  if (mode == 100) {
    int loop = 0;
    struct winner
    {
      int winM = -1;
      int winEps = -1;
      int winLife = 0;
    };
    stack<winner> winners;
    ofstream changeOfs("changelog.txt");
    while (true) {
      loop++;
      if (loop % 10 == 1) {
        endTime = clock();
        double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        if (nowTime > 360000.0) break;
      }

      if (loop % 200 == 77) {
        OutputHaipara();
        cout << "Updated" << endl;
      }

      if (loop % 2000 == 1833) {
        mode = 1000;
        MODE = 1000;
        // 1000ケース実行
        double sumScore = 0;
        ofstream ofsScore("Score.txt");
        double hi = 0;
        for (int _ = 0; _ < (1000); ++_) {
          if (_ % 100 == 0) {
            // cout << _ << endl;
          }
          if (_ < 100) {
            OpenOfs1000Out(_);
          }

          Input(mode, _);
          n = hyperN[m][iEps];
          hyperSolverNum = hyperSolver[m][iEps];
          hyperMinDiff = hyperMinDiffArr[m][iEps];
          hyperMaxRound = hyperMaxRoundArr[m][iEps];
          hyperStep1 = hyperStep1Arr[m][iEps];
          hyperStep2 = hyperStep2Arr[m][iEps];

          InitNumArray(mode);
          if (!numPairArrOK) continue;

          if (_ < 100) {
            OutputArrayAsString(mode);
          }

          double score = Simulate(mode);
          sumScore += score;

          if (_ < 100) {
            OutputAnsToOfs1000Out();
            CloseOfs1000Out();
          }

          hi += score / hyperMaxScore[m][iEps];
          ofsScore << score / hyperMaxScore[m][iEps] << ' ' << fixed
            << setprecision(6) << score << ' ' << hyperMaxScore[m][iEps]
            << endl;
            // cout << _ << ' ' << score << endl;
            // if (_ < 100)
            //   ofsScore << setw(4) << setfill('0') << _ << " : " << score <<
            //   endl;
        }
        changeOfs << hi << endl;
        changeOfs << "sumScore = " << sumScore << endl;
        ofsScore << "sumScore = " << sumScore << endl;
        ofsScore.close();
        mode = 100;
        MODE = 100;
      }

      Input(mode);
      m = Rand() % 91 + 10;
      iEps = Rand() % 41;
      eps = iEps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = Rand() % m;

      int initMode = loop % 2;

      // if (initMode == 99) {
      //   iEps = Rand() % 6 + 35;
      //   m = Rand() % 11 + 90;
      //   eps = iEps / 100.0;
      //   n = hyperN[m][iEps];
      //   hyperSolverNum = hyperSolver[m][iEps];
      //   hyperMinDiff = hyperMinDiffArr[m][iEps];
      //   hyperMaxRound = hyperMaxRoundArr[m][iEps];
      //   hyperStep1 = hyperStep1Arr[m][iEps];
      //   hyperStep2 = hyperStep2Arr[m][iEps];
      //   hyperSolverNum = 148;
      //   n = 4;
      // }
      // initMode = 0;
      if (initMode == 1 || !winners.empty()) {
        // 上下左右の丸コピー
        int nm = m;
        int niEps = iEps;
        while (true) {
          int ra = Rand() % 4;
          nm = m + dx[ra];
          niEps = iEps + dy[ra];
          if (10 <= nm && nm <= 100 && 0 <= niEps && niEps <= 40) break;
        }
        if (!winners.empty()) {
          winners.top().winLife--;
          m = winners.top().winM;
          iEps = winners.top().winEps;
          if (winners.top().winLife == 0) winners.pop();
          nm = m;
          niEps = iEps;
          while (true) {
            int ra = Rand() % 4;
            nm = m + dx[ra];
            niEps = iEps + dy[ra];
            if (10 <= nm && nm <= 100 && 0 <= niEps && niEps <= 40) break;
          }
          swap(m, nm);
          swap(iEps, niEps);
          eps = (double)iEps / 100.0;
        }
        n = hyperN[nm][niEps];
        hyperSolverNum = hyperSolver[nm][niEps];
        hyperMinDiff = hyperMinDiffArr[nm][niEps];
        hyperMaxRound = hyperMaxRoundArr[nm][niEps];
        hyperStep1 = hyperStep1Arr[nm][niEps];
        hyperStep2 = hyperStep2Arr[nm][niEps];

        // 隣を改変
        if (winners.empty() && Rand() % 2 == 0) {
          //vector<int> selection = { 1114, 1134, 1134, 1134, 1152,
          //                         1152, 1135, 1135, 1135 };
          //vector<int> selection = { 1186,1196 ,1187,1197};
          // vector<int> selection = { 1186,1196 ,1187,1197};
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };
          hyperSolverNum = selection[Rand() % selection.size()];

          n = n + Rand() % 5 - 2;
          n = max(n, 4);
          n = min(n, 100);
          hyperStep1 = hyperStep1 + Rand() % 3 - 1;
          hyperStep1 = max(1, hyperStep1);
          hyperStep2 = hyperStep2 + Rand() % 3 - 1;
          hyperStep2 = max(1, hyperStep2);
          if (Rand() % 2 == 0) {
            hyperStep2 = hyperStep1;
          }

          if (Rand() % 2 == 0) {
            hyperMinDiff = hyperMinDiff + Rand() % 3 - 1;
            hyperMinDiff = max(0, hyperMinDiff);
          }
          if (Rand() % 2 == 0) {
            hyperMinDiff = 0;
          }
          hyperMaxRound = hyperMaxRound + Rand() % 3 - 1;
          hyperMaxRound = max(1, hyperMaxRound);
        }
      }
      // ランダム生成
      else {
        n = hyperN[m][iEps];
        hyperSolverNum = hyperSolver[m][iEps];
        hyperMinDiff = hyperMinDiffArr[m][iEps];
        hyperMaxRound = hyperMaxRoundArr[m][iEps];
        hyperStep1 = hyperStep1Arr[m][iEps];
        hyperStep2 = hyperStep2Arr[m][iEps];

        if (false && Rand() % 2 == 0) {
          //vector<int> selection = { 1186,1196 ,1187,1197};
          //vector<int> selection = { 1189,1199 };
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };
          hyperSolverNum = selection[Rand() % selection.size()];
        }
        else {
          // vector<int> selection = {
          //     105,  106,  107,  115,  116,  117,  125,  126,  127,  135,
          //     136,  137,  109,  119,  129,  139,  1100, 1110, 1120, 1130,
          //     1101, 1111, 1121, 1131, 1104, 1114, 1124, 1134, 1152, 1152,
          //     1152, 1152, 1152, 1134, 1134, 1134, 1134};
          //vector<int> selection = { 1114, 1134, 1134, 1134, 1152,
          //                         1152, 1135, 1135, 1135 };
          //vector<int> selection = { 1186,1196 ,1187,1197};
          // vector<int> selection = { 1189,1199 };
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };

          hyperSolverNum = selection[Rand() % selection.size()];

          n = hyperN[m][iEps] + Rand() % 21 - 10;
          if (Rand() % 2 == 0) {
            n = hyperN[m][iEps] - 1;
          }
          n = max(n, 4);
          n = min(n, 100);
          if (Rand() % 2 == 0) {
            hyperStep1 = hyperStep1Arr[m][iEps] + Rand() % 3 - 1;
            hyperStep1 = max(1, hyperStep1);
            hyperStep2 = hyperStep2Arr[m][iEps] + Rand() % 3 - 1;
            hyperStep2 = max(1, hyperStep2);
            if (Rand() % 2 == 0) {
              hyperStep2 = hyperStep1;
            }
          }


          if (Rand() % 2 == 0) {
            hyperMinDiff = hyperMinDiffArr[m][iEps] + Rand() % 15 - 7;
            hyperMinDiff = max(0, hyperMinDiff);
          }
          if (Rand() % 2 == 0) {
            hyperMinDiff = 0;
          }
          if (Rand() % 2 == 0) {
            hyperMaxRound = hyperMaxRoundArr[m][iEps] + Rand() % 5 - 2;
            hyperMaxRound = max(1, hyperMaxRound);
          }

        }
      }

      InitNumArray(mode);
      if (!numPairArrOK) continue;

      int nown = hyperN[m][iEps];
      int nowhyperSolverNum = hyperSolver[m][iEps];
      int nowhyperMinDiff = hyperMinDiffArr[m][iEps];
      int nowhyperMaxRound = hyperMaxRoundArr[m][iEps];
      int nowhyperStep1 = hyperStep1Arr[m][iEps];
      int nowhyperStep2 = hyperStep2Arr[m][iEps];

      int chan = n;
      int chahyperSolverNum = hyperSolverNum;
      int chahyperMinDiff = hyperMinDiff;
      int chahyperMaxRound = hyperMaxRound;
      int chahyperStep1 = hyperStep1;
      int chahyperStep2 = hyperStep2;

      int same = 1;
      if (nown != chan) same = 0;
      if (nowhyperSolverNum != chahyperSolverNum) same = 0;
      if (nowhyperMinDiff != chahyperMinDiff) same = 0;
      if (nowhyperMaxRound != chahyperMaxRound) same = 0;
      if (nowhyperStep1 != chahyperStep1) same = 0;
      if (nowhyperStep2 != chahyperStep2) same = 0;
      if (same) continue;

      int winCount = 0;
      double score = 0;
      double matchCount = 0;
      int CHAMP = 17;
      int LOSE = 3;
      int loseCount = 0;
      for (int i = 0; i < (CHAMP + LOSE + 100); ++i) {
        matchCount++;
        for (int j = 0; j < (100); ++j) judgeArr[j] = Rand() % m;

        n = nown;
        hyperSolverNum = nowhyperSolverNum;
        hyperMinDiff = nowhyperMinDiff;
        hyperMaxRound = nowhyperMaxRound;
        hyperStep1 = nowhyperStep1;
        hyperStep2 = nowhyperStep2;
        InitNumArray(mode);
        double nowscore = Simulate(mode);

        n = chan;
        hyperSolverNum = chahyperSolverNum;
        hyperMinDiff = chahyperMinDiff;
        hyperMaxRound = chahyperMaxRound;
        hyperStep1 = chahyperStep1;
        hyperStep2 = chahyperStep2;
        InitNumArray(mode);
        double chascore = Simulate(mode);
        score += chascore;
        // cout << chascore << ' ' << nowscore << endl;
        if (chascore > nowscore) {
          winCount++;
          // cout << "WIN" << endl;
          // cout << hyperSolverNum << ' ' << chascore << ' ' << nowscore << endl;
        }
        else if (chascore < nowscore) {
          // cout << "LOSE" << endl;
          loseCount++;
        }
        if (loseCount > LOSE) {
          winCount = 0;
          break;
        }
        if (winCount == CHAMP) break;
      }
      score /= matchCount;

      // double score = 0;
      // for (int i = 0; i < (30); ++i) score += Simulate(mode);
      // score /= 30;
      // if (score >= hyperMaxScore[m][iEps]) {

      if (winCount == CHAMP && (hyperMaxScore[m][iEps] < score)) {
        changeOfs << loop << ' ' << hyperSolver[m][iEps] << ' '
          << hyperSolverNum << ' ' << m << ' ' << eps << ' '
          << hyperN[m][iEps] << ' ' << n << ' '
          << hyperMaxScore[m][iEps] << ' ' << score << ' ' << score * n
          << ' ' << matchCount << endl;
        cout << loop << ' ' << hyperSolver[m][iEps] << ' ' << hyperSolverNum
          << ' ' << m << ' ' << eps << ' ' << hyperN[m][iEps] << ' ' << n
          << ' ' << hyperMaxScore[m][iEps] << ' ' << score << ' '
          << score * n << ' ' << matchCount << endl;
        hyperMaxScore[m][iEps] = score;
        hyperN[m][iEps] = n;
        hyperSolver[m][iEps] = hyperSolverNum;
        hyperMinDiffArr[m][iEps] = hyperMinDiff;
        hyperMaxRoundArr[m][iEps] = hyperMaxRound;
        hyperStep1Arr[m][iEps] = hyperStep1;
        hyperStep2Arr[m][iEps] = hyperStep2;

        winner er;
        er.winM = m;
        er.winEps = iEps;
        er.winLife = 10;
        winners.push(er);
      }
    }

    cout << "loop = " << loop << endl;

    OutputHaipara();
  }

  // 1000ケース実行
  if (mode == 1000) {
    double sumScore = 0;
    ofstream ofsScore("Score.txt");
    double hi = 0;
    for (int _ = 0; _ < (1000); ++_) {
      if (_ % 100 == 0) {
        cout << _ << endl;
      }
      if (_ < 100) {
        OpenOfs1000Out(_);
      }

      Input(mode, _);
      n = hyperN[m][iEps];
      hyperSolverNum = hyperSolver[m][iEps];
      hyperMinDiff = hyperMinDiffArr[m][iEps];
      hyperMaxRound = hyperMaxRoundArr[m][iEps];
      hyperStep1 = hyperStep1Arr[m][iEps];
      hyperStep2 = hyperStep2Arr[m][iEps];

      InitNumArray(mode);
      if (!numPairArrOK) continue;

      if (_ < 100) {
        OutputArrayAsString(mode);
      }

      double score = Simulate(mode);
      sumScore += score;

      if (_ < 100) {
        OutputAnsToOfs1000Out();
        CloseOfs1000Out();
      }

      hi += score / hyperMaxScore[m][iEps];
      ofsScore << score / hyperMaxScore[m][iEps] << ' ' << fixed
        << setprecision(6) << score << ' ' << hyperMaxScore[m][iEps]
        << endl;
        // cout << _ << ' ' << score << endl;
        // if (_ < 100)
        //   ofsScore << setw(4) << setfill('0') << _ << " : " << score << endl;
    }
    cout << hi << endl;
    cout << "sumScore = " << sumScore << endl;
    ofsScore << "sumScore = " << sumScore << endl;
    ofsScore.close();
  }

  // ケース19を100回実行
  if (mode == 1001) {
    mode = 1000;
    // 1000ケース実行
    double sumScore = 0;
    ofstream ofsScore("Score.txt");
    double hi = 0;
    for (int looop = 0; looop < (100); ++looop) {
      for (int _ = 19; _ < 20; ++_) {
        if (_ % 100 == 0) {
          cout << _ << endl;
        }
        if (_ < 100) {
          OpenOfs1000Out(_);
        }

        Input(mode, _);
        n = hyperN[m][iEps];
        hyperSolverNum = hyperSolver[m][iEps];
        hyperMinDiff = hyperMinDiffArr[m][iEps];
        hyperMaxRound = hyperMaxRoundArr[m][iEps];
        hyperStep1 = hyperStep1Arr[m][iEps];
        hyperStep2 = hyperStep2Arr[m][iEps];
        hyperSolverNum = 1154;

        InitNumArray(mode);
        if (!numPairArrOK) continue;

        if (_ < 100) {
          OutputArrayAsString(mode);
        }

        double score = Simulate(mode);
        sumScore += score;

        if (_ < 100) {
          OutputAnsToOfs1000Out();
          CloseOfs1000Out();
        }

        hi += score / hyperMaxScore[m][iEps];
        cout << hyperSolverNum << ' ' << score << endl;
      }
    }

    cout << hi << endl;
    cout << "sumScore = " << sumScore << endl;
    ofsScore << "sumScore = " << sumScore << endl;
    ofsScore.close();
  }
}

int main()
{
  for (int i = 0; i < (105); ++i) {
    for (int j = 0; j < (105); ++j) {
      com[i][j] = 0;
      if (j > i) continue;
      for (int k = 0; k < (j); ++k) {
        com[i][j] += log(i - k);
        com[i][j] -= log(j - k);
      }
    }
  }

  for (int i = 0; i < 100; i++) {
    maxNumArray[i] = MaxNumArray[i];
  }
  for (int i = 0; i < 11; i++) {
    for (int j = 0; j < 100; j++) {
      real_real_maxNumArray[i][j] = Real_real_maxNumArray[i][j];
    }
  }
  for (int i = 0; i < 100; i++) {
    real_maxNumArray[i] = Real_maxNumArray[i];
  }
  for (int i = 0; i < 101; i++) {
    for (int j = 0; j < 41; j++) {
      hyperN[i][j] = HyperN[i][j];
      hyperMaxScore[i][j] = HyperMaxScore[i][j];
      hyperSolver[i][j] = HyperSolver[i][j];
      hyperMinDiffArr[i][j] = HyperMinDiffArr[i][j];
      hyperMaxRoundArr[i][j] = HyperMaxRoundArr[i][j];
      hyperStep1Arr[i][j] = HyperStep1Arr[i][j];
      hyperStep2Arr[i][j] = HyperStep2Arr[i][j];
    }
  }

  MODE = 1000;
  solve(MODE);

  return 0;
}
