#include "Hypers.h"

#include <algorithm>
#include <array>
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

const std::array<int, 4> dx = { -1, 1, 0, 0 };
const std::array<int, 4> dy = { 0, 0, -1, 1 };

class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

namespace /* 乱数ライブラリ */
{
  static uint32_t rand32()
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


  static double rand_01()
  {
    return (rand32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

// 定数定義
const int MAX_N = 100;
const int MAX_M = 101;
const int MIN_M = 10;
const int MAX_TURN = 100;
const int MAX_IEPS = 41;
const int MAX_ATTEMPTS = 5000;
const int MAX_KOUHO_LIMIT = 200;
const int FILE_NUM_DIGITS = 4;
const int ARRAY_SIZE = 1000;
const int SPECIAL_MODE = 1000;
const int INITIAL_DIFF = 1000;
const int FLIP_ITERATIONS = 1000;
const double TIME_LIMIT_MS = 360000.0;
const double MISMATCH_PENALTY = 0.9;

namespace
{
  // 共通変数
  std::array<std::array<double, 105>, 105> com;
  int MODE = 0;
  const int TURN = 100;
  int m, i_eps;  // m: 符号化するグラフの数(10-100), i_eps: エラー率(整数%, 0-40)
  double eps;    // eps: エラー率(0.00-0.40)

  int n;  // n: グラフの頂点数(4-100)
  std::array<std::array<std::array<int, 100>, 100>, 100> a;  // a[i][j][k]: i番目の符号化グラフの隣接行列
  std::array<std::array<int, 100>, 100> b;  // b[i][j]: 受信したノイズ付きグラフ

  std::vector<int> numSingleArr(ARRAY_SIZE);
  std::vector<std::array<int, 2>> numPairArr(ARRAY_SIZE);
  int numPairArrOK;

  std::vector<std::array<int, 3>> numThreeArr(ARRAY_SIZE);
  std::vector<std::array<int, 4>> numFourArr(ARRAY_SIZE);

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
  std::array<int, 100> maxNumArray;
  std::array<std::array<int, 100>, 11> real_real_maxNumArray;
  std::array<int, 100> real_maxNumArray;
  std::array<std::array<int, 41>, 101> hyperN;
  std::array<std::array<double, 41>, 101> hyperMaxScore;
  std::array<std::array<int, 41>, 101> hyperSolver;
  std::array<std::array<int, 41>, 101> hyperMinDiffArr;
  std::array<std::array<int, 41>, 101> hyperMaxRoundArr;
  std::array<std::array<int, 41>, 101> hyperStep1Arr;
  std::array<std::array<int, 41>, 101> hyperStep2Arr;
}  // namespace

// 共通ユーティリティ関数
namespace
{
  // 数値をゼロ埋めした文字列に変換
  string numberToString(int num, int digits = FILE_NUM_DIGITS)
  {
    string result;
    for (int i = 0; i < digits; ++i) {
      result += (char)(num % 10 + '0');
      num /= 10;
    }
    reverse(result.begin(), result.end());
    return result;
  }

  // イプシロン文字列をパース
  int parseEpsilon(const string& sEps)
  {
    return (sEps[2] - '0') * 10 + (sEps[3] - '0');
  }

  // 対称的な値を設定（a[i][j][k]とa[i][k][j]に同じ値を設定）
  void setSymmetricValue(int i, int j, int k, int value)
  {
    a[i][j][k] = value;
    a[i][k][j] = value;
  }
}

// Input（m, eps, i_epsの設定）
namespace
{
  std::array<int, 100> judgeArr;
  void Input(int mode, int case_num = 0)
  {
    if (mode == 0) {
      cin >> m;
      string sEps;
      cin >> sEps;
      i_eps = parseEpsilon(sEps);
      eps = (double)i_eps / 100.0;
    }
    else if (mode == SPECIAL_MODE) {
      string fileNameIfs = "./in/" + numberToString(case_num) + ".txt";
      ifstream ifs(fileNameIfs);

      ifs >> m;
      string sEps;
      ifs >> sEps;
      i_eps = parseEpsilon(sEps);
      eps = (double)i_eps / 100.0;

      for (int i = 0; i < (100); ++i) ifs >> judgeArr[i];
    }
    else if (mode == 100 || mode == 110) {
      m = rand32() % 91 + 10;
      i_eps = rand32() % 41;
      eps = i_eps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = rand32() % m;
    }
    else if (mode == 333) {
      m = case_num * 10;
      i_eps = 40;
      eps = i_eps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = rand32() % m;
    }
  }
}  // namespace

// Output
namespace
{
  ofstream ofs1000Out;
  void OpenOfs1000Out(int case_num)
  {
    string fileNameOfs = "./out/" + numberToString(case_num) + ".txt";
    ofs1000Out.open(fileNameOfs);
  }

  void close_ofs_1000_out() { ofs1000Out.close(); }

  void output_array_as_string(int mode)
  {
    if (mode == 0) {
      cout << n << endl;
      for (int i = 0; i < m; ++i) {
        string s;
        for (int j = 0; j < n; ++j) {
          for (int k = j + 1; k < n; ++k) { s += (char)(a[i][j][k] + '0'); }
        }
        cout << s << endl;
      }
      fflush(stdout);
    }
    else if (mode == SPECIAL_MODE) {
      ofs1000Out << n << endl;
      ofs1000Out << "# n = " << n << endl;
      ofs1000Out << "# hyperSolverNum = " << hyperSolverNum << endl;
      ofs1000Out << "# hyperMinDiff = " << hyperMinDiff << endl;
      ofs1000Out << "# hyperMaxRound = " << hyperMaxRound << endl;
      ofs1000Out << "# hyperStep1 = " << hyperStep1 << endl;
      ofs1000Out << "# hyperStep2 = " << hyperStep2 << endl;
      for (int i = 0; i < m; ++i) {
        string s;
        for (int j = 0; j < n; ++j) {
          for (int k = j + 1; k < n; ++k) { s += (char)(a[i][j][k] + '0'); }
        }
        ofs1000Out << s << endl;
        ofs1000Out << "# " << numPairArr[i][0] << ' ' << numPairArr[i][1] << endl;
      }
    }
  }

  std::array<int, TURN> answersFor1000Out;
  void output_ans_to_ofs_1000_out()
  {
    for (int i = 0; i < (100); ++i) { ofs1000Out << answersFor1000Out[i] << endl; }
  }

  // ハイパーパラメータ配列を出力する共通関数
  template<typename T>
  void OutputHyperArray(ofstream& ofs, const string& typeName, const string& arrayName, const std::array<std::array<T, 41>, 101>& arr)
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

  void output_hyper_params()
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

// 前方宣言
void setSymmetricValue(int i, int j, int k, int value);
int calculateScore(std::array<int, 100>& f);
int performGreedyElimination(std::array<int, 100>& f, std::array<int, 100>& cnt, int initialRes);
int performRandomizedGreedyElimination(std::array<int, 100>& f, std::array<int, 100>& cnt, int initialRes);
void collectCandidates(vector<int>& kouho, const std::array<int, 100>& f, int minConnectivity);
void sortResultsDescending(int& res1, int& res2, int& res3);
void sortResultsDescending(int& res1, int& res2, int& res3, int& res4);
int findBestMatchingPairCores(int res1, int res2, const std::vector<std::array<int, 2>>& numPairArr);
int findBestMatchingPairCoresWithFilter(int res1, int res2, const std::vector<std::array<int, 2>>& numPairArr, const std::vector<int>& omoteArr, int tei);
template<size_t N>
void initializeBitsets(std::array<int, 100>& f, std::array<bitset<100>, N>& bif, std::array<bitset<100>, 100>& bib, bitset<100>& bione);
int createAndExpandCore(std::array<int, 100>& f, vector<int>& cores, int cliqueSize, int markValue, bool useReturn = true);
template<size_t N>
void flipOptimization(std::array<int, 100>& f, std::array<bitset<100>, N>& bif, std::array<bitset<100>, 100>& bib, bitset<100>& bione,
  int& score, int& res1, int& res2, int flipLoop);
int evaluateAndOptimizeCores(std::array<int, 100>& f, vector<int>& cores1, vector<int>& cores2,
  const std::vector<int>& omoteArr, int tei, int& real_score, int& real_argRes);

// グラフにノイズを適用（確率epsで各辺を反転）
void apply_noise_to_graph(int x)
{
  for (int j = 0; j < n; ++j) {
    for (int k = j + 1; k < n; ++k) {
      b[j][k] = a[x][j][k];
      if (rand_01() < eps) {
        b[j][k] = 1 - b[j][k];
      }
      b[k][j] = b[j][k];
    }
  }
}

int judgeNum;
// ノイズ付きグラフを受信してb配列に格納
void receive_noisy_graph(int mode, int turn = 0)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      b[i][j] = 0;
    }
  }
  if (mode == 0) {
    string s;
    cin >> s;
    int ite = 0;
    for (int i = 0; i < n; ++i) {
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
    apply_noise_to_graph(judge);
  }
}

// numArray
namespace
{
  std::array<int, 100> numArr;

  // 共通のグラフ初期化関数（完全グラフのサイズで符号化）
  void set_graph_from_num_array(int arraySize = 100)
  {
    for (int i = 0; i < m; ++i) {
      int num = numArr[i];
      for (int j = 0; j < n; ++j) {
        for (int k = j + 1; k < n; ++k) {
          setSymmetricValue(i, j, k, (k < num) ? 1 : 0);
        }
      }
    }
  }

  // 条件付きグラフ初期化関数（2つのクリークで符号化）
  void set_graph_from_pair_array()
  {
    for (int i = 0; i < m; ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      for (int j = 0; j < n; ++j) {
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

  // 条件付きグラフ初期化関数（3つのクリークで符号化）
  void set_graph_from_three_array()
  {
    for (int i = 0; i < m; ++i) {
      int num1 = numThreeArr[i][0];
      int num2 = numThreeArr[i][1];
      int num3 = numThreeArr[i][2];
      for (int j = 0; j < n; ++j) {
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

  // 条件付きグラフ初期化関数（4つのクリークで符号化）
  void SetGraphFromFourArray()
  {
    for (int i = 0; i < m; ++i) {
      int num1 = numFourArr[i][0];
      int num2 = numFourArr[i][1];
      int num3 = numFourArr[i][2];
      int num4 = numFourArr[i][3];
      for (int j = 0; j < n; ++j) {
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

    set_graph_from_num_array();
  }

  void InitNumArray2()
  {
    for (int i = 0; i < (100); ++i) { numArr[i] = maxNumArray[i]; }

    set_graph_from_num_array();
  }

  void InitNumArray3()
  {
    for (int i = 0; i < (100); ++i) { numArr[i] = real_real_maxNumArray[(m + 9) / 10][i]; }

    set_graph_from_num_array();
  }

  void InitNumArray4()
  {
    for (int i = 0; i < (100); ++i) numArr[i] = 0;
    for (int i = 0; i < (20); ++i) { numArr[i] = (i + 1) * 5; }
    set_graph_from_num_array();
  }

  void InitNumArray5()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 1) { continue; }
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 0) { continue; }
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
        if (i % 2 == 0) { continue; }
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
      for (int i = (n / 2) - 1; i >= 0; --i) {
        if (i % 2 == 1) { continue; }
        numArr[cnt] = i + 1;
        cnt++;
        numArr[cnt] = n - i;
        cnt++;
      }
    }

    set_graph_from_num_array();
  }

  void InitNumArray6()
  {
    if (n % 2 == 0) {
      int cnt = 0;
      for (int j = 0; j < (3); ++j) {
        for (int i = (n / 2) - 1; i >= 0; --i) {
          if (i % 3 != j) { continue; }
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
          if (i % 3 != j % 3) { continue; }
          numArr[cnt] = i + 1;
          cnt++;
          numArr[cnt] = n - i;
          cnt++;
        }
      }
    }

    set_graph_from_num_array();
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

    for (int i = 0; i < m; ++i) {
      int num = numArr[i];
      for (int j = 0; j < n; ++j) {
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

    set_graph_from_num_array();
  }

  void InitNumArray9()
  {
    int cnt = 0;
    for (int j = 1; j < 3; ++j) {
      for (int i = (n)-1; i >= 0; --i) {
        if (i % 2 != j % 2) { continue; }
        numArr[cnt] = i + 1;
        cnt++;
      }
    }

    set_graph_from_num_array();
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    set_graph_from_pair_array();
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    set_graph_from_pair_array();
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    set_graph_from_pair_array();
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    set_graph_from_pair_array();
    if (cnt < m) numPairArrOK = 0;
  }

  void InitNumArray14()
  {
    for (int i = 0; i < m; ++i) {
      int cnt = 0;
      for (int j = 0; j < n; ++j) {
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
        if (cnt >= m) { break; }
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
            if (cnt >= 500) { break; }
          }
          j += hyperStep2;
          if (cnt >= 500) { break; }
        }
        k += hyperStep2;
        if (cnt >= 500) { break; }
      }
    }

    for (int i = 0; i < (cnt / 2); ++i) {
      swap(numThreeArr[i][0], numThreeArr[cnt - 1 - i][0]);
      swap(numThreeArr[i][1], numThreeArr[cnt - 1 - i][1]);
      swap(numThreeArr[i][2], numThreeArr[cnt - 1 - i][2]);
    }

    set_graph_from_three_array();
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
            // 3クラ
            if (i + j + k <= n) {
              cnt++;
            }
            i += hyperStep1;
          }
          j += hyperStep2;
        }
        k += hyperStep2;
        if (cnt >= m) { break; }
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
            // 3クラ
            if (i + j + k <= n) {
              numThreeArr[cnt][0] = i;
              numThreeArr[cnt][1] = j;
              numThreeArr[cnt][2] = k;
              cnt++;
            }
            i += hyperStep1;
            if (cnt >= 500) { break; }
          }
          j += hyperStep2;
          if (cnt >= 500) { break; }
        }
        k += hyperStep2;
        if (cnt >= 500) { break; }
      }
    }

    for (int i = 0; i < (cnt / 2); ++i) {
      swap(numThreeArr[i][0], numThreeArr[cnt - 1 - i][0]);
      swap(numThreeArr[i][1], numThreeArr[cnt - 1 - i][1]);
      swap(numThreeArr[i][2], numThreeArr[cnt - 1 - i][2]);
    }

    set_graph_from_three_array();
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
          if (cnt >= m) { break; }
        }
        l += hyperStep2;
        if (cnt >= m) { break; }
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

    if (numPairArrOK == 0) { return; }

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
              if (cnt >= 500) { break; }
            }
            j += hyperStep1;
            if (cnt >= 500) { break; }
          }
          k += hyperStep2;
          if (cnt >= 500) { break; }
        }
        l += hyperStep2;
        if (cnt >= 500) { break; }
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
  std::vector<int> omoteArr(ARRAY_SIZE);
  void InitNumArray18()
  {
    for (int i = 0; i < (ARRAY_SIZE); ++i)omoteArr[i] = 0;
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    for (int i = 0; i < m; ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      int omote = omoteArr[i];
      for (int j = 0; j < n; ++j) {
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
    for (int i = 0; i < (ARRAY_SIZE); ++i)omoteArr[i] = 0;
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
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    for (int i = 0; i < m; ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      int omote = omoteArr[i];
      for (int j = 0; j < n; ++j) {
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
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
          for (int k = 0; k < (n); ++k) {
            a[i][j][k] = 0;
          }
        }
      }
      for (int i = 0; i < m; ++i) {
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
    for (int i = 0; i < (ARRAY_SIZE); ++i)omoteArr[i] = 0;
    numPairArrOK = 1;
    int cnt = 0;
    for (int i = n; i > 0; i -= hyperStep1) {
      numSingleArr[cnt] = i;
      omoteArr[cnt] = 0;
      cnt++;
      numSingleArr[cnt] = i;
      omoteArr[cnt] = 1;
      cnt++;
      if (cnt > MAX_KOUHO_LIMIT) { break; }
    }

    for (int i = 0; i < m; ++i) {
      int num1 = numSingleArr[i];
      int omote = omoteArr[i];
      for (int j = 0; j < n; ++j) {
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
    int ra = hyperSolverNum % ARRAY_SIZE / 10;

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

// Solver共通関数群

// cnt[]とf[]を初期化し、全頂点の次数をカウント
void initializeCountAndFlags(std::array<int, 100>& cnt, std::array<int, 100>& f, int& res)
{
  cnt = {};
  f = {};
  res = n;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cnt[i] += b[i][j];
    }
    f[i] = 1;
  }
}

// 単一コアのベストマッチを見つける
int findBestMatchingSingleCore(int res, const std::array<int, 100>& numArr)
{
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    int num = numArr[i];
    if (abs(num - res) < diff) {
      diff = abs(num - res);
      argRes = i;
    }
  }
  return argRes;
}

// 2つのコアのベストマッチを見つける
int findBestMatchingPairCores(int res1, int res2, const std::vector<std::array<int, 2>>& numPairArr)
{
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }
  return argRes;
}

// 2つのコアのベストマッチを見つける（omoteArrフィルタ付き）
int findBestMatchingPairCoresWithFilter(int res1, int res2, const std::vector<std::array<int, 2>>& numPairArr, const std::vector<int>& omoteArr, int tei)
{
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    if (omoteArr[i] != tei % 2) {
      continue;
    }
    int num1 = numPairArr[i][0];
    int num2 = numPairArr[i][1];
    if (abs(num1 - res1) + abs(num2 - res2) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2);
      argRes = i;
    }
  }
  return argRes;
}

// f[]で指定された頂点集合に対してcnt[]を再計算
void recountForSubset(std::array<int, 100>& cnt, const std::array<int, 100>& f)
{
  cnt = {};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (f[i] && f[j]) {
        cnt[i] += b[i][j];
      }
    }
  }
}

// 第2コアの探索（第1コアの補集合から）
int performSecondCoreSearch(std::array<int, 100>& f, std::array<int, 100>& cnt,
  std::array<int, 100>& ff, int res1)
{
  // 第1コアの補集合を第2コアの候補とする
  int res = n - res1;
  for (int i = 0; i < n; ++i) {
    f[i] = 1 - f[i];
    if (f[i] == 0) ff[i] = 1;
  }

  // 第2コア候補の次数を再計算
  recountForSubset(cnt, f);

  // 貪欲除去を実行
  res = performGreedyElimination(f, cnt, res);

  return res;
}

// ビットセット最適化を実行する共通関数
template<size_t N>
void performBitsetOptimization(std::array<int, 100>& f, int& res1, int& res2)
{
  int score = calculateScore(f);

  std::array<bitset<100>, N> bif = {};
  std::array<bitset<100>, 100> bib = {};
  bitset<100> bione(0);
  initializeBitsets(f, bif, bib, bione);

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
}

int solver_1()
{
  std::array<vector<int>, 100> keep;

  std::array<int, 100> cnt = {};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cnt[i] += b[i][j];
    }
  }

  for (int i = 0; i < m; ++i) {
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
  for (int i = 0; i < 100; ++i) {
    if (keep[i].size()) {
      int sz = keep[i].size();
      res = keep[i][sz / 2];
      break;
    }
  }

  return res;
}

int solver_2()
{
  std::array<int, 100> cnt, f;
  int res;

  // 初期化
  initializeCountAndFlags(cnt, f, res);

  // 貪欲除去
  res = performGreedyElimination(f, cnt, res);

  // ベストマッチを見つける
  return findBestMatchingSingleCore(res, numArr);
}

int solver_3()
{
  std::array<vector<int>, 2> vec;
  std::array<int, 100> f = {};
  for (int i = 0; i < n; ++i) {
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
    if (vec[0].empty() || vec[1].empty()) { break; }
    std::array<std::array<double, 2>, 100> cnt;
    for (int i = 0; i < n; ++i) for (int j = 0; j < (2); ++j) cnt[i][j] = 0;
    std::array<vector<int>, 2> nxt;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
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
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    int num = numArr[i];
    if (abs(num - res) < diff) {
      diff = abs(num - res);
      argRes = i;
    }
  }

  return argRes;
}

int solver_4()
{
  std::array<int, 100> cnt, f;
  int res;

  // 初期化
  initializeCountAndFlags(cnt, f, res);

  // 貪欲除去
  res = performGreedyElimination(f, cnt, res);

  if (res >= 20) {
    int res2 = 0;
    vector<int> vec;
    for (int i = 0; i < n; ++i) {
      if (f[i]) {
        vec.push_back(i);
      }
    }
    std::array<int, 100> ff = {};
    double kijun = (eps * eps + (1.0 - eps) * (1.0 - eps)) / 2.0;
    for (int i = 0; i < n; ++i) {
      std::array<int, 2> tri = {};
      for (int j = 0; j < (res); ++j) {
        int jj = vec[j];
        if (jj == i) { continue; }
        for (int k = j + 1; k < res; ++k) {
          int kk = vec[k];
          if (kk == i) { continue; }
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

  // ベストマッチを見つける
  return findBestMatchingSingleCore(res, numArr);
}

int solver_5()
{
  std::array<int, 100> cnt, f, ff = {};
  int res;

  // 初期化
  initializeCountAndFlags(cnt, f, res);

  // 第1コアを見つける
  res = performGreedyElimination(f, cnt, res);
  int res1 = res;

  // 第2コアを探索
  int res2 = performSecondCoreSearch(f, cnt, ff, res1);

  // hyperMaxRoundチェック
  if (res2 <= hyperMaxRound) res2 = 0;

  // ベストマッチを見つける
  return findBestMatchingPairCores(res1, res2, numPairArr);
}

int solver_6()
{
  std::array<int, 100> cnt, f;
  std::array<int, 100> ff = {};
  int res;

  // 初期化
  initializeCountAndFlags(cnt, f, res);

  // 第1コアを見つける
  res = performGreedyElimination(f, cnt, res);
  int res1 = res;

  // 第2コアを探索
  int res2 = performSecondCoreSearch(f, cnt, ff, res1);

  // hyperMaxRoundチェック
  if (res2 <= hyperMaxRound) {
    res2 = 0;
    for (int i = 0; i < n; ++i) f[i] = 0;
  }

  // ffにマージ
  for (int i = 0; i < n; ++i) {
    if (f[i]) ff[i] = 2;
  }
  for (int i = 0; i < n; ++i) f[i] = ff[i];

  // ビットセット最適化
  performBitsetOptimization<3>(f, res1, res2);

  // ベストマッチを見つける
  return findBestMatchingPairCores(res1, res2, numPairArr);
}

int solver_7()
{
  std::array<int, 100> cnt = {};
  std::array<int, 100> f = {};
  int ff[100] = {};
  int res = n;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cnt[i] += b[i][j];
    }
    f[i] = 1;
  }

  res = performGreedyElimination(f, cnt, res);

  int res1 = res;
  res = n - res;
  for (int i = 0; i < n; ++i) { f[i] = 1 - f[i]; }
  for (int i = 0; i < n; ++i) {
    if (f[i] == 0) ff[i] = 1;
  }
  for (int i = 0; i < n; ++i) { cnt[i] = 0; }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (f[i] && f[j]) cnt[i] += b[i][j];
    }
  }
  res = performGreedyElimination(f, cnt, res);
  int res2 = res;
  if (res2 <= hyperMaxRound) {
    res2 = 0;
    for (int i = 0; i < n; ++i) f[i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    if (f[i]) ff[i] = 2;
  }

  for (int i = 0; i < n; ++i) f[i] = ff[i];

  int score = calculateScore(f);

  std::array<bitset<100>, 3> bif = {};
  std::array<bitset<100>, 100> bib = {};
  bitset<100> bione(0);
  initializeBitsets(f, bif, bib, bione);

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);

  return findBestMatchingPairCores(res1, res2, numPairArr);
}

int solver_8()
{
  int cnt = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) { cnt += b[i][j]; }
  }
  cnt = min(cnt, m - 1);
  return cnt;
}

int solver_9()
{
  map<P, int> mp;
  int fff[10][100];
  P kp[10];

  int kcnt[100] = {};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      kcnt[i] += b[i][j];
    }
  }

  for (int _ = 0; _ < (10); ++_) {
    std::array<int, 100> cnt = {};
    std::array<int, 100> f = {};
    std::array<int, 100> ff = {};
    int res = n;

    for (int i = 0; i < n; ++i) {
      cnt[i] = kcnt[i];
      f[i] = 1;
    }

    res = performRandomizedGreedyElimination(f, cnt, res);

    int res1 = res;
    res = n - res;
    for (int i = 0; i < n; ++i) { f[i] = 1 - f[i]; }
    for (int i = 0; i < n; ++i) {
      if (f[i] == 0) ff[i] = 1;
    }
    for (int i = 0; i < n; ++i) { cnt[i] = 0; }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (f[i] && f[j]) cnt[i] += b[i][j];
      }
    }
    res = performRandomizedGreedyElimination(f, cnt, res);
    int res2 = res;
    if (res2 <= hyperMaxRound) {
      res2 = 0;
      for (int i = 0; i < n; ++i) f[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
      if (f[i]) ff[i] = 2;
    }

    for (int i = 0; i < n; ++i) f[i] = ff[i];

    for (int i = 0; i < n; ++i) { fff[_][i] = f[i]; }
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

  std::array<int, 100> f = {};
  int res1, res2;
  res1 = maxP.first;
  res2 = maxP.second;
  for (int i = 0; i < (10); ++i) {
    if (maxP == kp[i]) {
      for (int j = 0; j < (100); ++j) { f[j] = fff[i][j]; }
    }
  }

  int score = calculateScore(f);

  std::array<bitset<100>, 3> bif = {};
  std::array<bitset<100>, 100> bib = {};
  bitset<100> bione(0);
  initializeBitsets(f, bif, bib, bione);

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);

  return findBestMatchingPairCores(res1, res2, numPairArr);
}

int solver_10()
{
  map<P, int> mp;
  int fff[31][100];
  P kp[31];

  int kcnt[100] = {};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      kcnt[i] += b[i][j];
    }
  }

  for (int _ = 0; _ < (31); ++_) {
    std::array<int, 100> cnt = {};
    std::array<int, 100> f = {};
    std::array<int, 100> ff = {};
    int res = n;

    for (int i = 0; i < n; ++i) {
      cnt[i] = kcnt[i];
      f[i] = 1;
    }

    res = performRandomizedGreedyElimination(f, cnt, res);

    int res1 = res;
    res = n - res;
    for (int i = 0; i < n; ++i) { f[i] = 1 - f[i]; }
    for (int i = 0; i < n; ++i) {
      if (f[i] == 0) ff[i] = 1;
    }
    for (int i = 0; i < n; ++i) { cnt[i] = 0; }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (f[i] && f[j]) cnt[i] += b[i][j];
      }
    }
    res = performRandomizedGreedyElimination(f, cnt, res);
    int res2 = res;
    if (res2 <= hyperMaxRound) {
      res2 = 0;
      for (int i = 0; i < n; ++i) f[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
      if (f[i]) ff[i] = 2;
    }

    for (int i = 0; i < n; ++i) f[i] = ff[i];

    for (int i = 0; i < n; ++i) { fff[_][i] = f[i]; }
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

  std::array<int, 100> f = {};
  int res1, res2;
  res1 = maxP.first;
  res2 = maxP.second;
  for (int i = 0; i < (10); ++i) {
    if (maxP == kp[i]) {
      for (int j = 0; j < (100); ++j) { f[j] = fff[i][j]; }
    }
  }

  return findBestMatchingPairCores(res1, res2, numPairArr);
}

// 指定サイズのクリークを見つける汎用関数
bool findClique(const vector<int>& kouho, std::array<int, 100>& f, vector<int>& cores, int cliqueSize, int markValue)
{
  for (int loop1 = 0; loop1 < MAX_ATTEMPTS; ++loop1) {
    vector<int> core(cliqueSize);
    for (int i = 0; i < cliqueSize; ++i) {
      while (true) {
        core[i] = kouho[rand32() % kouho.size()];
        bool duplicate = false;
        for (int j = 0; j < i; ++j) {
          if (core[j] == core[i]) {
            duplicate = true;
            break;
          }
        }
        if (!duplicate) {
          break;
        }
      }
    }

    bool isClique = true;
    for (int i = 0; i < cliqueSize; ++i) {
      for (int j = i + 1; j < cliqueSize; ++j) {
        if (!b[core[i]][core[j]]) {
          isClique = false;
          break;
        }
      }
      if (!isClique) {
        break;
      }
    }

    if (isClique) {
      for (int i = 0; i < cliqueSize; ++i) {
        f[core[i]] = markValue;
        cores.push_back(core[i]);
      }
      return true;
    }
  }
  return false;
}


// コアを大きくしていく共通関数
void expandCore(vector<int>& cores, std::array<int, 100>& f, int markValue)
{
  while (true) {
    int sz = cores.size();
    int arg = -1;
    int ma = -1;
    for (int i = 0; i < n; ++i) {
      if (f[i] != 0) { continue; }
      int cnt = 0;
      for (int j = 0; j < (sz); ++j) {
        if (b[i][cores[j]]) cnt++;
      }
      if (ma < cnt && cnt >= (sz + 2) / 2) {
        ma = cnt;
        arg = i;
      }
    }
    if (arg == -1) { break; }
    f[arg] = markValue;
    cores.push_back(arg);
  }
}

// スコア計算の共通関数
int calculateScore(std::array<int, 100>& f)
{
  int score = 0;
  for (int i = 0; i < n; ++i) {
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
  return score;
}

// 貪欲除去アルゴリズムの共通関数
// 最小次数の頂点を繰り返し除去し、残った頂点数を返す
int performGreedyElimination(std::array<int, 100>& f, std::array<int, 100>& cnt, int initialRes)
{
  int res = initialRes;

  while (res > 1) {
    int mi = 1000;
    int arg = -1;
    for (int i = 0; i < n; ++i) {
      if (f[i] && cnt[i] < mi && cnt[i] < (res + 1) / 2) {
        mi = cnt[i];
        arg = i;
      }
    }
    if (arg == -1) { break; }

    for (int i = 0; i < n; ++i) {
      if (i == arg) { continue; }
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }

    res--;
    f[arg] = 0;
  }

  return res;
}

// ランダム選択を含む貪欲除去アルゴリズムの共通関数
// 同じ最小次数を持つ頂点からランダムに選択して除去
int performRandomizedGreedyElimination(std::array<int, 100>& f, std::array<int, 100>& cnt, int initialRes)
{
  int res = initialRes;

  while (res > 1) {
    int mi = 1000;
    vector<int> arv;

    // 最小次数の頂点を収集
    for (int i = 0; i < n; ++i) {
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

    if (arv.empty()) { break; }

    // ランダムに選択
    int arg = arv[rand32() % arv.size()];

    // 隣接頂点の次数を更新
    for (int i = 0; i < n; ++i) {
      if (i == arg) { continue; }
      if (f[i] && b[i][arg]) {
        cnt[i]--;
      }
    }

    res--;
    f[arg] = 0;
  }

  return res;
}

// 候補収集の共通関数
void collectCandidates(vector<int>& kouho, const std::array<int, 100>& f, int minConnectivity)
{
  kouho.clear();
  for (int i = 0; i < n; ++i) {
    if (f[i] != 0) {
      continue;
    }
    int cnt = 0;
    for (int j = 0; j < n; ++j) {
      if (f[j] == 0) {
        cnt += b[i][j];
      }
    }
    if (cnt <= minConnectivity) {
      continue;
    }
    kouho.push_back(i);
  }
}

// 結果ソートの共通関数 (降順にソートして返す)
void sortResultsDescending(int& res1, int& res2, int& res3)
{
  vector<int> resv = { res1, res2, res3 };
  sort(resv.begin(), resv.end());
  res1 = resv[2];
  res2 = resv[1];
  res3 = resv[0];
}

void sortResultsDescending(int& res1, int& res2, int& res3, int& res4)
{
  vector<int> resv = { res1, res2, res3, res4 };
  sort(resv.begin(), resv.end());
  res1 = resv[3];
  res2 = resv[2];
  res3 = resv[1];
  res4 = resv[0];
}

// ビットセット初期化の共通関数
template<size_t N>
void initializeBitsets(std::array<int, 100>& f, std::array<bitset<100>, N>& bif, std::array<bitset<100>, 100>& bib, bitset<100>& bione)
{
  // bifの初期化
  for (int i = 0; i < N; ++i) {
    bif[i].reset();
  }
  for (int i = 0; i < n; ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  // bibの初期化
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      bib[i][j] = b[i][j];
    }
  }

  // bioneの初期化
  bione.reset();
  for (int i = 0; i < n; ++i) {
    bione[i] = 1;
  }
}

// コアを作成して拡張する汎用関数
// 戻り値: 成功時はcores1のサイズ、失敗時は-1（returnの場合）または0（continueの場合）
int createAndExpandCore(std::array<int, 100>& f, vector<int>& cores, int cliqueSize, int markValue, bool useReturn)
{
  vector<int> kouho;
  for (int i = 0; i < n; ++i) kouho.push_back(i);

  if (kouho.size() < cliqueSize) {
    return useReturn ? -1 : 0;
  }

  if (!findClique(kouho, f, cores, cliqueSize, markValue)) {
    return useReturn ? -1 : 0;
  }

  if (cores.size() == 0) {
    return useReturn ? -1 : 0;
  }

  // コアを大きくしていく
  expandCore(cores, f, markValue);

  return cores.size();
}

// 最適化とスコア評価の共通関数
int evaluateAndOptimizeCores(std::array<int, 100>& f, vector<int>& cores1, vector<int>& cores2,
  const std::vector<int>& omoteArr, int tei, int& real_score, int& real_argRes)
{
  int res1 = cores1.size();
  int res2 = cores2.size();

  int score = calculateScore(f);

  std::array<bitset<100>, 3> bif = {};
  std::array<bitset<100>, 100> bib = {};
  bitset<100> bione(0);
  initializeBitsets(f, bif, bib, bione);

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
  if (res2 > res1) swap(res1, res2);
  int argRes = findBestMatchingPairCoresWithFilter(res1, res2, numPairArr, omoteArr, tei);

  // スコア計算
  int tmpScore = 0;
  for (int i = 0; i < n; ++i) {
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

  return tmpScore;
}

// フリップ最適化の共通関数
template<size_t N>
void flipOptimization(std::array<int, 100>& f, std::array<bitset<100>, N>& bif, std::array<bitset<100>, 100>& bib, bitset<100>& bione,
  int& score, int& res1, int& res2, int flipLoop)
{
  for (int _ = 0; _ < (flipLoop); ++_) {
    int x = rand32() % n;
    int ra = rand32() % 3;
    while (ra == f[x]) {
      ra = rand32() % 3;
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
  for (int i = 0; i < n; ++i) {
    if (f[i] == 1) res1++;
    if (f[i] == 2) res2++;
  }
}

int solver_11()
{
  std::array<int, 100> f = {};

  // コア1を作る
  vector<int> cores1;
  if (createAndExpandCore(f, cores1, 4, 1, true) == -1) {
    return 0;
  }

  // コア2を作る
  vector<int> cores2;
  vector<int> kouho;
  collectCandidates(kouho, f, 4);
  if (kouho.size() >= 4) {
    findClique(kouho, f, cores2, 4, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      expandCore(cores2, f, 2);
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();

  int score = calculateScore(f);

  std::array<bitset<100>, 3> bif = {};
  std::array<bitset<100>, 100> bib = {};
  bitset<100> bione(0);
  initializeBitsets(f, bif, bib, bione);

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
  if (res2 > res1) swap(res1, res2);
  return findBestMatchingPairCores(res1, res2, numPairArr);
}

int solver_12()
{
  std::array<int, 100> f = {};

  // コア1を作る
  vector<int> cores1;
  if (createAndExpandCore(f, cores1, 4, 1, true) == -1) {
    return 0;
  }

  // コア2を作る
  vector<int> cores2;
  vector<int> kouho;
  collectCandidates(kouho, f, 4);
  if (kouho.size() >= 4) {
    findClique(kouho, f, cores2, 4, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      expandCore(cores2, f, 2);
    }
  }

  // コア3を作る
  vector<int> cores3;
  if (cores2.size() > 0) {
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores3, 4, 3);
      if (cores3.size() > 0) {
        // コア3を大きくしていく
        expandCore(cores3, f, 3);
      }
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();
  int res3 = cores3.size();

  int score = calculateScore(f);

  std::array<bitset<100>, 4> bif = {};
  std::array<bitset<100>, 100> bib = {};
  for (int i = 0; i < n; ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      bib[i][j] = b[i][j];
    }
  }

  bitset<100> bione(0);
  for (int i = 0; i < n; ++i) { bione[i] = 1; }

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
  res3 = 0;
  for (int i = 0; i < n; ++i) {
    if (f[i] == 3) res3++;
  }
  sortResultsDescending(res1, res2, res3);
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    int num1 = numThreeArr[i][0];
    int num2 = numThreeArr[i][1];
    int num3 = numThreeArr[i][2];
    if (abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3);
      argRes = i;
    }
  }

  return argRes;
}

// 4コア
int solver_13()
{
  std::array<int, 100> f = {};

  // コア1を作る
  vector<int> cores1;
  if (createAndExpandCore(f, cores1, 4, 1, true) == -1) {
    return 0;
  }

  // コア2を作る
  vector<int> cores2;
  vector<int> kouho;
  collectCandidates(kouho, f, 4);
  if (kouho.size() >= 4) {
    findClique(kouho, f, cores2, 4, 2);
    if (cores2.size() > 0) {
      // コア2を大きくしていく
      expandCore(cores2, f, 2);
    }
  }

  // コア3を作る
  vector<int> cores3;
  if (cores2.size() > 0) {
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores3, 4, 3);
      if (cores3.size() > 0) {
        // コア3を大きくしていく
        expandCore(cores3, f, 3);
      }
    }
  }

  // コア4を作る
  vector<int> cores4;
  if (cores3.size() > 0) {
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores4, 4, 4);
      if (cores4.size() > 0) {
        // コア4を大きくしていく
        expandCore(cores4, f, 4);
      }
    }
  }

  int res1 = cores1.size();
  int res2 = cores2.size();
  int res3 = cores3.size();
  int res4 = cores4.size();

  int score = calculateScore(f);

  std::array<bitset<100>, 5> bif = {};
  std::array<bitset<100>, 100> bib = {};
  for (int i = 0; i < n; ++i) {
    if (f[i] > 0) {
      bif[f[i]][i] = 1;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      bib[i][j] = b[i][j];
    }
  }

  bitset<100> bione(0);
  for (int i = 0; i < n; ++i) { bione[i] = 1; }

  int flipLoop = FLIP_ITERATIONS;
  if (MODE == 0) flipLoop = 10000;
  flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
  res3 = 0;
  res4 = 0;
  for (int i = 0; i < n; ++i) {
    if (f[i] == 3) res3++;
    if (f[i] == 4) res4++;
  }
  sortResultsDescending(res1, res2, res3, res4);
  int diff = INITIAL_DIFF;
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    int num1 = numFourArr[i][0];
    int num2 = numFourArr[i][1];
    int num3 = numFourArr[i][2];
    int num4 = numFourArr[i][3];
    if (abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) + abs(num4 - res4) < diff) {
      diff = abs(num1 - res1) + abs(num2 - res2) + abs(num3 - res3) + abs(num4 - res4);
      argRes = i;
    }
  }

  return argRes;
}

int solver_14()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores2, 4, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    if (res2 > res1) swap(res1, res2);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
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
      if (real_minDiff == 0) {
        break;
      }
    }
  }

  return real_argRes;
}

int solver_15()
{
  map<int, int> argMap;

  for (int wataruoop = 0; wataruoop < (5); ++wataruoop) {
    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores2, 4, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    if (res2 > res1) swap(res1, res2);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
      int num1 = numPairArr[i][0];
      int num2 = numPairArr[i][1];
      if (abs(num1 - res1) + abs(num2 - res2) < diff) {
        diff = abs(num1 - res1) + abs(num2 - res2);
        argRes = i;
      }
    }

    argMap[argRes]++;
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

void reverse_b()
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      b[i][j] = 1 - b[i][j];
    }
  }
}

// 11表裏
int solver_16()
{
  int real_argRes = 0;
  int real_score = 0;
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int tei = 0; tei < 2; ++tei) {
    b = keepB;

    if (tei % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores2, 4, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    evaluateAndOptimizeCores(f, cores1, cores2, omoteArr, tei, real_score, real_argRes);
  }

  return real_argRes;
}

// 11表裏5セット
int solver_17()
{
  int real_argRes = 0;
  int real_score = 0;
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int tei = 0; tei < 10; ++tei) {
    b = keepB;

    if (tei % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores2, 4, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    evaluateAndOptimizeCores(f, cores1, cores2, omoteArr, tei, real_score, real_argRes);
  }

  return real_argRes;
}

// 0.0用
int solver_18()
{
  int visit[100] = {};
  int bb[100][100];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      bb[i][j] = b[i][j];
    }
  }

  vector<P> vp;
  for (int i = 0; i < n; ++i) {
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
      for (int j = 0; j < n; ++j) {
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
    if (cnt1 == 1) {
      continue;
    }
    vp.push_back(P(cnt1, cnt2));
  }
  sort(vp.begin(), vp.end());
  int argRes = 0;
  for (int i = 0; i < m; ++i) {
    if (zeroPairs[i] == vp) {
      argRes = i;
      break;
    }
  }
  return argRes;
}

// 14表裏
int solver_19()
{
  int real_argRes[2] = {};
  int real_minDiff[2] = { 1000,1000 };
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int wataruoop = 0; wataruoop < (10); ++wataruoop) {
    b = keepB;

    if (wataruoop % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 4) {
      findClique(kouho, f, cores2, 4, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    if (res2 > res1) swap(res1, res2);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
      if (omoteArr[i] != wataruoop % 2) {
        continue;
      }
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
int solver_20()
{
  int real_argRes = 0;
  int real_score = 0;
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int tei = 0; tei < 2; ++tei) {
    b = keepB;

    if (tei % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 4, 1, false) == 0) {
      continue;
    }

    int res1 = cores1.size();
    int res2 = 0;

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
      if (omoteArr[i] != tei % 2) {
        continue;
      }
      int num1 = numSingleArr[i];
      if (abs(num1 - res1) < diff) {
        diff = abs(num1 - res1);
        argRes = i;
      }
    }

    // スコア計算
    int tmpScore = 0;
    for (int i = 0; i < n; ++i) {
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
int solver_21()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 5, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 5) {
      findClique(kouho, f, cores2, 5, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    if (res2 > res1) swap(res1, res2);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
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
      if (real_minDiff == 0) {
        break;
      }
    }
  }

  return real_argRes;
}

// 14亜種
int solver_22()
{
  int real_argRes = 0;
  int real_minDiff = 1000;

  for (int wataruoop = 0; wataruoop < (15); ++wataruoop) {
    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 3, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 2);
    if (kouho.size() >= 3) {
      findClique(kouho, f, cores2, 3, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    int res1 = cores1.size();
    int res2 = cores2.size();

    int score = calculateScore(f);

    std::array<bitset<100>, 3> bif = {};
    std::array<bitset<100>, 100> bib = {};
    bitset<100> bione(0);
    initializeBitsets(f, bif, bib, bione);

    int flipLoop = FLIP_ITERATIONS;
    if (MODE == 0) flipLoop = 10000;
    flipOptimization(f, bif, bib, bione, score, res1, res2, flipLoop);
    if (res2 > res1) swap(res1, res2);
    int diff = INITIAL_DIFF;
    int argRes = 0;
    for (int i = 0; i < m; ++i) {
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
      if (real_minDiff == 0) {
        break;
      }
    }
  }

  return real_argRes;
}

// 16亜種
int solver_23()
{
  int real_argRes = 0;
  int real_score = 0;
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int tei = 0; tei < 2; ++tei) {
    b = keepB;

    if (tei % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 3, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 2);
    if (kouho.size() >= 3) {
      findClique(kouho, f, cores2, 3, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    evaluateAndOptimizeCores(f, cores1, cores2, omoteArr, tei, real_score, real_argRes);
  }

  return real_argRes;
}

// 16亜種
int solver_24()
{
  int real_argRes = 0;
  int real_score = 0;
  std::array<std::array<int, 100>, 100> keepB = b;

  for (int tei = 0; tei < 2; ++tei) {
    b = keepB;

    if (tei % 2 == 1) {
      reverse_b();
    }

    std::array<int, 100> f = {};

    // コア1を作る
    vector<int> cores1;
    if (createAndExpandCore(f, cores1, 5, 1, false) == 0) {
      continue;
    }

    // コア2を作る
    vector<int> cores2;
    vector<int> kouho;
    collectCandidates(kouho, f, 4);
    if (kouho.size() >= 5) {
      findClique(kouho, f, cores2, 5, 2);
      if (cores2.size() > 0) {
        // コア2を大きくしていく
        expandCore(cores2, f, 2);
      }
    }

    evaluateAndOptimizeCores(f, cores1, cores2, omoteArr, tei, real_score, real_argRes);
  }

  return real_argRes;
}

// Solver関数ポインタ型の定義
typedef int (*SolverFunction)();

// Solver関数ポインタ配列
const SolverFunction solverFunctions[] = {
    nullptr,  // index 0 (not used)
    solver_1, solver_2, solver_3, solver_4, solver_5,
    solver_6, solver_7, solver_8, solver_9, solver_10,
    solver_11, solver_12, solver_13, solver_14, solver_15,
    solver_16, solver_17, solver_18, solver_19, solver_20,
    solver_21, solver_22, solver_23, solver_24
};

// 受信したノイズ付きグラフから元のグラフ番号を推定
int solver_use_count[25] = { 0 };  // 各Solver関数の使用回数をカウント
int DecodeGraphIndex(int mode, int turn = 0)
{
  int res = 0;
  int num = hyperSolverNum % 10 + hyperSolverNum / ARRAY_SIZE * 10;
  solver_use_count[num]++;  // 使用回数をカウント

  // 範囲チェックとSolver関数の呼び出し
  if (num >= 1 && num <= 24) {
    res = solverFunctions[num]();
  }

  return res;
}

double Simulate(int mode)
{
  double res = 1e9 / n;
  for (int turn = 0; turn < (100); ++turn) {
    receive_noisy_graph(mode, turn);
    int ans = DecodeGraphIndex(mode, turn);
    answersFor1000Out[turn] = ans;
    int judge = judgeArr[turn];
    if (ans != judge) res *= MISMATCH_PENALTY;

    if (mode == 0) {
      cout << ans << endl;
      fflush(stdout);
    }
  }
  return res;
}

void solve(int mode)
{
  Timer timer;
  timer.start();

  // 提出用
  if (mode == 0) {
    Input(mode);

    n = hyperN[m][i_eps];
    hyperSolverNum = hyperSolver[m][i_eps];
    hyperMinDiff = hyperMinDiffArr[m][i_eps];
    hyperMaxRound = hyperMaxRoundArr[m][i_eps];
    hyperStep1 = hyperStep1Arr[m][i_eps];
    hyperStep2 = hyperStep2Arr[m][i_eps];

    InitNumArray(mode);

    output_array_as_string(mode);

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
        double nowTime = timer.get_elapsed_time();
        if (nowTime > TIME_LIMIT_MS) {
          break;
        }
      }

      if (loop % 200 == 77) {
        output_hyper_params();
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
          if (_ < 100) {
            OpenOfs1000Out(_);
          }

          Input(mode, _);
          n = hyperN[m][i_eps];
          hyperSolverNum = hyperSolver[m][i_eps];
          hyperMinDiff = hyperMinDiffArr[m][i_eps];
          hyperMaxRound = hyperMaxRoundArr[m][i_eps];
          hyperStep1 = hyperStep1Arr[m][i_eps];
          hyperStep2 = hyperStep2Arr[m][i_eps];

          InitNumArray(mode);
          if (!numPairArrOK) {
            continue;
          }

          if (_ < 100) {
            output_array_as_string(mode);
          }

          double score = Simulate(mode);
          sumScore += score;

          if (_ < 100) {
            output_ans_to_ofs_1000_out();
            close_ofs_1000_out();
          }

          hi += score / hyperMaxScore[m][i_eps];
          ofsScore << score / hyperMaxScore[m][i_eps] << ' ' << fixed
            << setprecision(6) << score << ' ' << hyperMaxScore[m][i_eps]
            << endl;
        }
        changeOfs << hi << endl;
        changeOfs << "sumScore = " << sumScore << endl;
        ofsScore << "sumScore = " << sumScore << endl;
        ofsScore.close();
        mode = 100;
        MODE = 100;
      }

      Input(mode);
      m = rand32() % 91 + 10;
      i_eps = rand32() % 41;
      eps = i_eps / 100.0;
      for (int i = 0; i < (100); ++i) judgeArr[i] = rand32() % m;

      int initMode = loop % 2;

      if (initMode == 1 || !winners.empty()) {
        // 上下左右の丸コピー
        int nm = m;
        int ni_eps = i_eps;
        while (true) {
          int ra = rand32() % 4;
          nm = m + dx[ra];
          ni_eps = i_eps + dy[ra];
          if (10 <= nm && nm <= 100 && 0 <= ni_eps && ni_eps <= 40) {
            break;
          }
        }
        if (!winners.empty()) {
          winners.top().winLife--;
          m = winners.top().winM;
          i_eps = winners.top().winEps;
          if (winners.top().winLife == 0) winners.pop();
          nm = m;
          ni_eps = i_eps;
          while (true) {
            int ra = rand32() % 4;
            nm = m + dx[ra];
            ni_eps = i_eps + dy[ra];
            if (10 <= nm && nm <= 100 && 0 <= ni_eps && ni_eps <= 40) {
              break;
            }
          }
          swap(m, nm);
          swap(i_eps, ni_eps);
          eps = (double)i_eps / 100.0;
        }
        n = hyperN[nm][ni_eps];
        hyperSolverNum = hyperSolver[nm][ni_eps];
        hyperMinDiff = hyperMinDiffArr[nm][ni_eps];
        hyperMaxRound = hyperMaxRoundArr[nm][ni_eps];
        hyperStep1 = hyperStep1Arr[nm][ni_eps];
        hyperStep2 = hyperStep2Arr[nm][ni_eps];

        // 隣を改変
        if (winners.empty() && rand32() % 2 == 0) {
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };
          hyperSolverNum = selection[rand32() % selection.size()];

          n = n + rand32() % 5 - 2;
          n = max(n, 4);
          n = min(n, 100);
          hyperStep1 = hyperStep1 + rand32() % 3 - 1;
          hyperStep1 = max(1, hyperStep1);
          hyperStep2 = hyperStep2 + rand32() % 3 - 1;
          hyperStep2 = max(1, hyperStep2);
          if (rand32() % 2 == 0) {
            hyperStep2 = hyperStep1;
          }

          if (rand32() % 2 == 0) {
            hyperMinDiff = hyperMinDiff + rand32() % 3 - 1;
            hyperMinDiff = max(0, hyperMinDiff);
          }
          if (rand32() % 2 == 0) {
            hyperMinDiff = 0;
          }
          hyperMaxRound = hyperMaxRound + rand32() % 3 - 1;
          hyperMaxRound = max(1, hyperMaxRound);
        }
      }
      // ランダム生成
      else {
        n = hyperN[m][i_eps];
        hyperSolverNum = hyperSolver[m][i_eps];
        hyperMinDiff = hyperMinDiffArr[m][i_eps];
        hyperMaxRound = hyperMaxRoundArr[m][i_eps];
        hyperStep1 = hyperStep1Arr[m][i_eps];
        hyperStep2 = hyperStep2Arr[m][i_eps];

        if (false && rand32() % 2 == 0) {
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };
          hyperSolverNum = selection[rand32() % selection.size()];
        }
        else {
          vector<int> selection = { 2183, 2184, 2193,2194,1186,1196 ,1187,1197,1134, 2131, 2132,1134, 2131, 2132 };

          hyperSolverNum = selection[rand32() % selection.size()];

          n = hyperN[m][i_eps] + rand32() % 21 - 10;
          if (rand32() % 2 == 0) {
            n = hyperN[m][i_eps] - 1;
          }
          n = max(n, 4);
          n = min(n, 100);
          if (rand32() % 2 == 0) {
            hyperStep1 = hyperStep1Arr[m][i_eps] + rand32() % 3 - 1;
            hyperStep1 = max(1, hyperStep1);
            hyperStep2 = hyperStep2Arr[m][i_eps] + rand32() % 3 - 1;
            hyperStep2 = max(1, hyperStep2);
            if (rand32() % 2 == 0) {
              hyperStep2 = hyperStep1;
            }
          }


          if (rand32() % 2 == 0) {
            hyperMinDiff = hyperMinDiffArr[m][i_eps] + rand32() % 15 - 7;
            hyperMinDiff = max(0, hyperMinDiff);
          }
          if (rand32() % 2 == 0) {
            hyperMinDiff = 0;
          }
          if (rand32() % 2 == 0) {
            hyperMaxRound = hyperMaxRoundArr[m][i_eps] + rand32() % 5 - 2;
            hyperMaxRound = max(1, hyperMaxRound);
          }

        }
      }

      InitNumArray(mode);
      if (!numPairArrOK) {
        continue;
      }

      int nown = hyperN[m][i_eps];
      int nowhyperSolverNum = hyperSolver[m][i_eps];
      int nowhyperMinDiff = hyperMinDiffArr[m][i_eps];
      int nowhyperMaxRound = hyperMaxRoundArr[m][i_eps];
      int nowhyperStep1 = hyperStep1Arr[m][i_eps];
      int nowhyperStep2 = hyperStep2Arr[m][i_eps];

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
      if (same) {
        continue;
      }

      int winCount = 0;
      double score = 0;
      double matchCount = 0;
      int CHAMP = 17;
      int LOSE = 3;
      int loseCount = 0;
      for (int i = 0; i < (CHAMP + LOSE + 100); ++i) {
        matchCount++;
        for (int j = 0; j < (100); ++j) judgeArr[j] = rand32() % m;

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

        if (chascore > nowscore) {
          winCount++;
        }
        else if (chascore < nowscore) {
          loseCount++;
        }
        if (loseCount > LOSE) {
          winCount = 0;
          break;
        }
        if (winCount == CHAMP) {
          break;
        }
      }
      score /= matchCount;

      if (winCount == CHAMP && (hyperMaxScore[m][i_eps] < score)) {
        changeOfs << loop << ' ' << hyperSolver[m][i_eps] << ' '
          << hyperSolverNum << ' ' << m << ' ' << eps << ' '
          << hyperN[m][i_eps] << ' ' << n << ' '
          << hyperMaxScore[m][i_eps] << ' ' << score << ' ' << score * n
          << ' ' << matchCount << endl;
        cout << loop << ' ' << hyperSolver[m][i_eps] << ' ' << hyperSolverNum
          << ' ' << m << ' ' << eps << ' ' << hyperN[m][i_eps] << ' ' << n
          << ' ' << hyperMaxScore[m][i_eps] << ' ' << score << ' '
          << score * n << ' ' << matchCount << endl;
        hyperMaxScore[m][i_eps] = score;
        hyperN[m][i_eps] = n;
        hyperSolver[m][i_eps] = hyperSolverNum;
        hyperMinDiffArr[m][i_eps] = hyperMinDiff;
        hyperMaxRoundArr[m][i_eps] = hyperMaxRound;
        hyperStep1Arr[m][i_eps] = hyperStep1;
        hyperStep2Arr[m][i_eps] = hyperStep2;

        winner er;
        er.winM = m;
        er.winEps = i_eps;
        er.winLife = 10;
        winners.push(er);
      }
    }

    cout << "loop = " << loop << endl;

    output_hyper_params();
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
      n = hyperN[m][i_eps];
      hyperSolverNum = hyperSolver[m][i_eps];
      hyperMinDiff = hyperMinDiffArr[m][i_eps];
      hyperMaxRound = hyperMaxRoundArr[m][i_eps];
      hyperStep1 = hyperStep1Arr[m][i_eps];
      hyperStep2 = hyperStep2Arr[m][i_eps];

      InitNumArray(mode);
      if (!numPairArrOK) { continue; }

      if (_ < 100) {
        output_array_as_string(mode);
      }

      double score = Simulate(mode);
      sumScore += score;

      if (_ < 100) {
        output_ans_to_ofs_1000_out();
        close_ofs_1000_out();
      }

      hi += score / hyperMaxScore[m][i_eps];
      ofsScore << score / hyperMaxScore[m][i_eps] << ' ' << fixed
        << setprecision(6) << score << ' ' << hyperMaxScore[m][i_eps]
        << endl;
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
        n = hyperN[m][i_eps];
        hyperSolverNum = hyperSolver[m][i_eps];
        hyperMinDiff = hyperMinDiffArr[m][i_eps];
        hyperMaxRound = hyperMaxRoundArr[m][i_eps];
        hyperStep1 = hyperStep1Arr[m][i_eps];
        hyperStep2 = hyperStep2Arr[m][i_eps];
        hyperSolverNum = 1154;

        InitNumArray(mode);
        if (!numPairArrOK) { continue; }

        if (_ < 100) {
          output_array_as_string(mode);
        }

        double score = Simulate(mode);
        sumScore += score;

        if (_ < 100) {
          output_ans_to_ofs_1000_out();
          close_ofs_1000_out();
        }

        hi += score / hyperMaxScore[m][i_eps];
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
  for (int i = 0; i < 105; ++i) {
    for (int j = 0; j < 105; ++j) {
      com[i][j] = 0;
      if (j > i) { continue; }
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

  for (int i = 0; i < 25; i++) {
    cout << "solver_use_count[" << i << "] = " << solver_use_count[i] << endl;
  }

  return 0;
}
