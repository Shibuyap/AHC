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

// 10のx乗
const ll D1 = 10LL;
const ll D2 = 100LL;
const ll D3 = 1000LL;
const ll D4 = 10000LL;
const ll D5 = 100000LL;
const ll D6 = 1000000LL;
const ll D7 = 10000000LL;
const ll D8 = 100000000LL;
const ll D9 = 1000000000LL;
const ll D10 = 10000000000LL;
const ll D11 = 100000000000LL;
const ll D12 = 1000000000000LL;
const ll D13 = 10000000000000LL;
const ll D14 = 100000000000000LL;
const ll D15 = 1000000000000000LL;
const ll D16 = 10000000000000000LL;
const ll D17 = 100000000000000000LL;
const ll D18 = 1000000000000000000LL;

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

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
std::normal_distribution<> normalDistribution(0, 1);

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int INF = 1001001001;
double TL = 2.5;
int mode;
double g_cost;
int g_query2Count;
clock_t g_startTime, g_endTime;

const int MAX_N = 20;
const int MAX_M = 20;
const int MAX_NNNN = MAX_N * MAX_N * MAX_N * MAX_N;

int n, m;
int nn, nnn, nnnn;
double eps;
int intEps;
int d[MAX_M];
int ii[MAX_M][250], jj[MAX_M][250];
int minIi[MAX_M], maxIi[MAX_M], minJj[MAX_M], maxJj[MAX_M];
int sum;

// ハイパラ調整用
const int M_INDEX_SIZE = 21;
const int EPS_INDEX_SIZE = 21;
int mIndex, epsIndex;

ll hyperParams[21][21] = { {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 10000000006, 10000000005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

// ローカルテスト用
int di[MAX_M], dj[MAX_M];
int localAnsGrid[MAX_N][MAX_N];
vector<P> localAns;
double sampleEs[3 * MAX_N * MAX_N];

double GetNowTime()
{
  g_endTime = clock();
  double nowTime = (double)(g_endTime - g_startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

bool IsNG(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y) { return true; }
  return false;
}

void InitIndices()
{
  mIndex = m;
  epsIndex = intEps;
}

void TestCaseGenerator(bool isMakeN, bool isMakeM, bool isMakeEps)
{
  for (int i = 0; i < MAX_N; ++i) {
    for (int j = 0; j < MAX_N; ++j) {
      localAnsGrid[i][j] = 0;
    }
  }
  if (isMakeN && isMakeM) {
    n = Rand() % (20 - 10 + 1) + 10;
    int ma = n * n / 20;
    m = Rand() % (ma - 2 + 1) + 2;
  }
  else if (isMakeN && !isMakeM) {
    while (true) {
      n = Rand() % (20 - 10 + 1) + 10;
      int ma = n * n / 20;
      if (m <= ma) break;
    }
  }
  else if (!isMakeN && !isMakeM) {
    int ma = n * n / 20;
    m = Rand() % (ma - 2 + 1) + 2;
  }

  if (isMakeEps) {
    eps = (double)(Rand() % (20 - 1 + 1) + 1) * 0.01;
    intEps = static_cast<int>(round(eps * 100 + 0.2));
  }
  {
    int paraA = (Rand() % (n * n / 2 - n * n / 5 + 1) + n * n / 5) / m;
    int paraD = Rand() % ((paraA - 4) - 0 + 1) + 0;
    for (int i = 0; i < m; ++i) {
      d[i] = Rand() % ((paraA + paraD) - (paraA - paraD) + 1) + (paraA - paraD);
      int grid[MAX_N][MAX_N];
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          grid[j][k] = 0;
        }
      }
      grid[n / 2][n / 2] = 1;
      for (int j = 0; j < d[i] - 1; ++j) {
        set<P> se;
        for (int k = 0; k < n; ++k) {
          for (int l = 0; l < n; ++l) {
            if (grid[k][l] == 1) {
              for (int z = 0; z < 4; ++z) {
                int nk = k + dx[z];
                int nl = l + dy[z];
                if (!IsNG(nk, nl) && grid[nk][nl] == 0) { se.emplace(nk, nl); }
              }
            }
          }
        }
        vector<P> v(se.begin(), se.end());
        P p = v[Rand() % v.size()];
        grid[p.first][p.second] = 1;
      }
      int miX = INF, miY = INF;
      int maX = -1, maY = -1;
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (grid[j][k] == 1) {
            miX = min(miX, j);
            miY = min(miY, k);
            maX = max(maX, j);
            maY = max(maY, k);
          }
        }
      }
      int now = 0;
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (grid[j][k] == 1) {
            ii[i][now] = j - miX;
            jj[i][now] = k - miY;
            now++;
          }
        }
      }

      // ローカルテスタ変数
      int sx = Rand() % ((n - 1 - (maX - miX)) - 0 + 1);
      int sy = Rand() % ((n - 1 - (maY - miY)) - 0 + 1);
      di[i] = sx;
      dj[i] = sy;
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (grid[j][k] == 1) { localAnsGrid[j + sx - miX][k + sy - miY]++; }
        }
      }
    }
    for (int i = 0; i < n * n * 2; ++i) {
      sampleEs[i] = normalDistribution(engine);
    }
  }

  InitIndices();
}

bool QueryAns(const vector<P>& ans, ofstream& ofs)
{
  bool res = true;
  if (mode == 0) {
    cout << "a " << ans.size();
    for (auto p : ans) {
      cout << " " << p.first << " " << p.second;
    }
    cout << endl;
    fflush(stdout);
    int xxx;
    cin >> xxx;
    if (xxx == 0) res = false;
  }
  else {
    ofs << "a " << ans.size();
    for (auto p : ans) {
      ofs << " " << p.first << " " << p.second;
    }
    ofs << endl;
    auto tmpAns = ans;
    sort(tmpAns.begin(), tmpAns.end());
    if (tmpAns != localAns) res = false;
  }
  if (!res) { g_cost += 1.0; }
  return res;
}

int Query1(int x, int y, ofstream& ofs)
{
  int val;
  if (mode == 0) {
    cout << "q 1 " << x << ' ' << y << endl;
    fflush(stdout);
    cin >> val;
  }
  else {
    ofs << "q 1 " << x << ' ' << y << endl;
    val = localAnsGrid[x][y];
  }
  g_cost += 1.0;
  return val;
}

int Query2Local(const vector<P>& points)
{
  int k = static_cast<int>(points.size());
  int vSum = 0;
  for (const auto& p : points) {
    vSum += localAnsGrid[p.first][p.second];
  }
  double mean = (static_cast<double>(k) - vSum) * eps + vSum * (1.0 - eps);
  double variance = k * eps * (1 - eps);
  double std_dev = std::sqrt(variance);  // 分散から標準偏差を計算

  double transformed_sample = sampleEs[g_query2Count] * std_dev + mean;
  g_query2Count++;

  return max(0, (int)round(transformed_sample));
}

int Query2(const vector<P>& points, ofstream& ofs)
{
  int res = 0;
  if (mode == 0) {
    cout << "q " << points.size();
    for (auto point : points) {
      cout << ' ' << point.first << ' ' << point.second;
    }
    cout << endl;
    fflush(stdout);
    cin >> res;
  }
  else {
    if (points.size() == 1) {
      res = Query1(points[0].first, points[0].second, ofs);
    }
    else {
      ofs << "q " << points.size();
      for (auto point : points) {
        ofs << ' ' << point.first << ' ' << point.second;
      }
      ofs << endl;
      res = Query2Local(points);
    }
  }
  g_cost += 1.0 / sqrt(points.size());
  return res;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  g_cost = 0;
  g_query2Count = 0;
  localAns.clear();
}

// 入力受け取り
void Input(int problemNum)
{
  if (mode == 0 || mode == 1) {
    string fileNameIfs = "./in/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameIfs += strNum + ".txt";

    ifstream ifs(fileNameIfs);

    // 標準入力する
    if (!ifs.is_open()) {
      cin >> n >> m >> eps;
      for (int i = 0; i < m; ++i) {
        cin >> d[i];
        for (int j = 0; j < d[i]; ++j) {
          cin >> ii[i][j] >> jj[i][j];
        }
      }
    }
    // ファイル入力する
    else {
      ifs >> n >> m >> eps;
      for (int i = 0; i < m; ++i) {
        ifs >> d[i];
        for (int j = 0; j < d[i]; ++j) {
          ifs >> ii[i][j] >> jj[i][j];
        }
      }
      for (int i = 0; i < m; ++i) {
        ifs >> di[i] >> dj[i];
      }
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          ifs >> localAnsGrid[i][j];
        }
      }
      for (int i = 0; i < n * n * 2; ++i) {
        ifs >> sampleEs[i];
      }
    }
  }
  else if (mode == 2) {
    ;
  }
}

void Initialize()
{
  nn = n * n;
  nnn = n * n * n;
  nnnn = n * n * n * n;

  sum = 0;
  for (int i = 0; i < m; ++i) {
    sum += d[i];
  }
  if (mode != 0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (localAnsGrid[i][j] > 0) { localAns.emplace_back(i, j); }
      }
    }
    sort(localAns.begin(), localAns.end());
  }

  for (int i = 0; i < m; ++i) {
    minIi[i] = 0;
    minJj[i] = 0;
    maxIi[i] = 0;
    maxJj[i] = 0;
    for (int j = 0; j < d[i]; ++j) {
      maxIi[i] = max(maxIi[i], ii[i][j]);
      maxJj[i] = max(maxJj[i], jj[i][j]);
    }
  }

  intEps = static_cast<int>(round(eps * 100 + 0.2));

  InitIndices();
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
ll CalcScore()
{
  ll res = static_cast<ll>(round(g_cost * 1000000LL));
  return res;
}

int m2_kouho[MAX_NNNN][MAX_N][MAX_N];
bool m2_ok[MAX_NNNN];
int m2_kakuteiValue1[MAX_N][MAX_N];
int m2_kakuteiValue2[MAX_N][MAX_N];
int m2_aimaiValue[MAX_N][MAX_N];
int m2_nums[MAX_N][MAX_N][MAX_M + 1];

// Method2のヘルパー関数：インデックスを指定して解答を出力
bool Method2_PrintAnsWithIndex(int idx, ofstream& ofs)
{
  vector<P> ans;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (m2_kouho[idx][i][j] >= 1 || m2_kakuteiValue2[i][j] >= 1) { ans.emplace_back(i, j); }
    }
  }
  bool isCorrect = QueryAns(ans, ofs);
  return isCorrect;
}

// Method2のヘルパー関数：候補数のカウント
int Method2_CountCandidates()
{
  int kouhoCount = 0;

  // カウント配列を初期化
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < m + 1; ++k) {
        m2_nums[i][j][k] = 0;
      }
    }
  }

  // 各候補パターンでの値をカウント
  for (int i = 0; i < nnnn; ++i) {
    if (!m2_ok[i]) continue;
    kouhoCount++;
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        m2_nums[j][k][m2_kouho[i][j][k]]++;
      }
    }
  }

  return kouhoCount;
}

// Method2のヘルパー関数：候補パターンの列挙
void Method2_EnumerateCandidates(int isTotyuu, int mNum0, int mNum1)
{
  if (isTotyuu == 0 || isTotyuu == 2) {
    // 2つのパターンの全組み合わせを列挙
    for (int i0 = 0; i0 < n; ++i0) {
      for (int j0 = 0; j0 < n; ++j0) {
        for (int i1 = 0; i1 < n; ++i1) {
          for (int j1 = 0; j1 < n; ++j1) {
            int idx = i0 * nnn + j0 * nn + i1 * n + j1;

            // パターン0の配置チェック
            for (int k = 0; k < d[mNum0]; ++k) {
              int x = i0 + ii[mNum0][k];
              int y = j0 + jj[mNum0][k];
              if (IsNG(x, y)) {
                m2_ok[idx] = false;
                break;
              }
              m2_kouho[idx][x][y]++;
            }

            // パターン1の配置チェック
            for (int k = 0; k < d[mNum1]; ++k) {
              int x = i1 + ii[mNum1][k];
              int y = j1 + jj[mNum1][k];
              if (IsNG(x, y)) {
                m2_ok[idx] = false;
                break;
              }
              m2_kouho[idx][x][y]++;
            }
          }
        }
      }
    }
  }
  else if (isTotyuu == 1) {
    // 1つのパターンのみ列挙
    for (int i = 0; i < nnnn; ++i) {
      m2_ok[i] = false;
    }
    for (int i0 = 0; i0 < n; ++i0) {
      for (int j0 = 0; j0 < n; ++j0) {
        int idx = i0 * n + j0;
        m2_ok[idx] = true;
        for (int k = 0; k < d[mNum0]; ++k) {
          int x = i0 + ii[mNum0][k];
          int y = j0 + jj[mNum0][k];
          if (IsNG(x, y)) {
            m2_ok[idx] = false;
            break;
          }
          m2_kouho[idx][x][y]++;
        }
      }
    }
  }
}

// Method2のヘルパー関数：クエリ位置の選択
P Method2_SelectQueryPosition()
{
  P result = { -1, -1 };
  bool findOne = false;

  // 既に確定した1があるかチェック
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (m2_kakuteiValue1[i][j] >= 1) { findOne = true; }
    }
  }

  if (findOne) {
    // 最も候補が少ない位置を選択
    int minMax = INF;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (m2_kakuteiValue1[i][j] != -1) continue;
        int ma = max(m2_nums[i][j][0], max(m2_nums[i][j][1], m2_nums[i][j][2]));
        if (ma < minMax) {
          minMax = ma;
          result.first = i;
          result.second = j;
        }
      }
    }
  }
  else {
    // 最も多くの候補がある位置を選択
    int maxNums = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        int numsSum = 0;
        for (int k = 1; k < m + 1; ++k) numsSum += m2_nums[i][j][k];
        if (numsSum > maxNums && m2_kakuteiValue1[i][j] == -1) {
          maxNums = numsSum;
          result.first = i;
          result.second = j;
        }
      }
    }
  }

  return result;
}

// Method2のヘルパー関数：単一候補の処理
void Method2_HandleSingleCandidate(ofstream& ofs)
{
  int idx = -1;
  for (int i = 0; i < nnnn; ++i) {
    if (m2_ok[i]) {
      idx = i;
      break;
    }
  }
  Method2_PrintAnsWithIndex(idx, ofs);
}

// Method2のヘルパー関数：クエリ実行と更新
void Method2_QueryAndUpdate(ofstream& ofs, int& cnt)
{
  P pos = Method2_SelectQueryPosition();
  int x = pos.first;
  int y = pos.second;

  int val = Query1(x, y, ofs);
  if (m2_kakuteiValue2[x][y] >= 0) val -= m2_kakuteiValue2[x][y];
  m2_kakuteiValue1[x][y] = val;
  cnt += val;

  // 候補を更新
  for (int i = 0; i < nnnn; ++i) {
    if (!m2_ok[i]) continue;
    if (m2_kouho[i][x][y] != val) m2_ok[i] = false;
  }
}

void Method2_Kakutei(ofstream& ofs)
{
  vector<P> ans;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (m2_kakuteiValue1[i][j] >= 1 || m2_kakuteiValue2[i][j] >= 1) { ans.emplace_back(i, j); }
    }
  }
  bool isCorrect = QueryAns(ans, ofs);
}

// Method2のヘルパー関数：初期化処理
void Method2_Initialize(int isTotyuu, const vector<vector<int>>& kakuteiValue1, const vector<vector<int>>& kakuteiValue2, int& cnt, int& cnt2)
{
  // 候補配列の初期化
  for (int i = 0; i < nnnn; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        m2_kouho[i][j][k] = 0;
      }
    }
    m2_ok[i] = true;
  }

  // 確定値の初期化
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m2_kakuteiValue1[i][j] = -1;
      m2_kakuteiValue2[i][j] = -1;
      m2_aimaiValue[i][j] = -1;
    }
  }

  // 途中からの場合は確定値をコピー
  cnt = 0;
  cnt2 = 0;
  if (isTotyuu != 0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        m2_kakuteiValue1[i][j] = kakuteiValue1[i][j];
        m2_kakuteiValue2[i][j] = kakuteiValue2[i][j];
        if (m2_kakuteiValue1[i][j] >= 0) cnt += m2_kakuteiValue1[i][j];
        if (m2_kakuteiValue2[i][j] >= 0) cnt2 += m2_kakuteiValue2[i][j];
      }
    }
  }
}

void Method2(ofstream& ofs, int isTotyuu = 0, int mNum0 = 0, int mNum1 = 1, const vector<vector<int>>& kakuteiValue1 = {}, const vector<vector<int>>& kakuteiValue2 = {})
{
  int cnt, cnt2;
  Method2_Initialize(isTotyuu, kakuteiValue1, kakuteiValue2, cnt, cnt2);

  // 候補パターンを列挙
  Method2_EnumerateCandidates(isTotyuu, mNum0, mNum1);

  while (true) {
    if (cnt + cnt2 == sum) {
      Method2_Kakutei(ofs);
      break;
    }
    else {
      // 候補数をカウント
      int kouhoCount = Method2_CountCandidates();

      if (kouhoCount == 1) {
        Method2_HandleSingleCandidate(ofs);
        break;
      }
      else {
        Method2_QueryAndUpdate(ofs, cnt);
      }
    }
  }
}

// M番号、D番号の配置がN*Nの枠内に収まっているか（障害物無し）
bool CanHaiti(int sx, int sy, int mNum)
{
  for (int i = 0; i < d[mNum]; ++i) {
    int x = sx + ii[mNum][i];
    int y = sy + jj[mNum][i];
    if (IsNG(x, y)) { return false; }
  }
  return true;
}

double m2_2_penalty[MAX_NNNN];

void Method2_2_Initialize()
{
  for (int i = 0; i < nnnn; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        m2_kouho[i][j][k] = 0;
      }
    }
    m2_ok[i] = true;
  }
  for (int i0 = 0; i0 < n; ++i0) {
    for (int j0 = 0; j0 < n; ++j0) {
      for (int i1 = 0; i1 < n; ++i1) {
        for (int j1 = 0; j1 < n; ++j1) {
          int idx = i0 * nnn + j0 * nn + i1 * n + j1;
          if (!CanHaiti(i0, j0, 0)) { m2_ok[idx] = false; }
          if (!CanHaiti(i1, j1, 1)) { m2_ok[idx] = false; }
        }
      }
    }
  }
  for (int i0 = 0; i0 < n; ++i0) {
    for (int j0 = 0; j0 < n; ++j0) {
      for (int i1 = 0; i1 < n; ++i1) {
        for (int j1 = 0; j1 < n; ++j1) {
          int idx = i0 * nnn + j0 * nn + i1 * n + j1;
          if (!m2_ok[idx]) continue;
          for (int k = 0; k < d[0]; ++k) {
            int x = i0 + ii[0][k];
            int y = j0 + jj[0][k];
            m2_kouho[idx][x][y]++;
          }
          for (int k = 0; k < d[1]; ++k) {
            int x = i1 + ii[1][k];
            int y = j1 + jj[1][k];
            m2_kouho[idx][x][y]++;
          }
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m2_kakuteiValue1[i][j] = -1;
      m2_kakuteiValue2[i][j] = -1;
      m2_aimaiValue[i][j] = -1;
    }
  }
}

void Method2_2_SearchKakuteiMasu(vector<pair<P, int>>& kakuteiValues)
{
  int valCount[MAX_N][MAX_N];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      valCount[i][j] = 0;
    }
  }
  for (int i = 0; i < nnnn; ++i) {
    if (!m2_ok[i]) continue;
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        valCount[j][k] += m2_kouho[i][j][k];
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (valCount[i][j] == 0) {
        kakuteiValues.emplace_back(P(i, j), 0);
        m2_kakuteiValue1[i][j] = 0;
        m2_kakuteiValue2[i][j] = 0;
      }
    }
  }
}

pair<double, int> m2_2_penaltySortedArray[MAX_NNNN];
int m2_2_kouhoCount;

int Method2_2_ChallengeAns(ofstream& ofs)
{
  for (int i = 0; i < nnnn; ++i) {
    if (!m2_ok[i]) {
      m2_2_penalty[i] = INF;
    }
    else {
      m2_2_penalty[i] = 0;
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (m2_kakuteiValue1[j][k] != -1) {
            if (m2_kouho[i][j][k] != m2_kakuteiValue1[j][k]) { m2_ok[i] = false; }
          }
          else if (m2_aimaiValue[j][k] != -1) {
            if (m2_kouho[i][j][k] != m2_aimaiValue[j][k]) { m2_2_penalty[i] += 100.0; }
          }
        }
      }
    }
  }

  pair<double, int> tmp;
  m2_2_kouhoCount = 0;
  for (int i = 0; i < nnnn; ++i) {
    if (m2_ok[i]) {
      tmp.first = m2_2_penalty[i];
      tmp.second = i;
      m2_2_penaltySortedArray[m2_2_kouhoCount] = tmp;
      m2_2_kouhoCount++;
    }
  }
  sort(m2_2_penaltySortedArray, m2_2_penaltySortedArray + m2_2_kouhoCount);
  if (m2_2_kouhoCount == 1) {
    int idx = m2_2_penaltySortedArray[0].second;
    Method2_PrintAnsWithIndex(idx, ofs);
    return 1;
  }
  else {
    if (m2_2_penaltySortedArray[0].first + 180.0 < m2_2_penaltySortedArray[1].first) {
      int idx = m2_2_penaltySortedArray[0].second;
      bool isCorect = Method2_PrintAnsWithIndex(idx, ofs);
      if (isCorect) {
        return 1;
      }
      else {
        m2_ok[idx] = false;
        return 2;
      }
    }
  }

  return 0;
}

void Method2_2(int query2Size, ofstream& ofs)
{
  Method2_2_Initialize();
  vector<pair<P, int>> kakuteiValues;
  // 確定マスを探す
  Method2_2_SearchKakuteiMasu(kakuteiValues);

  int isAC = 0;
  while (true) {
    if (g_cost > 20) { break; }
    // 回答するかどうか
    {
      int ret = Method2_2_ChallengeAns(ofs);
      if (ret == 0) {         // 回答せず
      }
      else if (ret == 1) {  // 正解
        isAC = 1;
        break;
      }
      else if (ret == 2) {  // 不正解
        continue;
      }
    }

    // 次のクエリ作成
    {
      int nums[MAX_N][MAX_N][3];
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          for (int k = 0; k < 3; ++k) {
            nums[i][j][k] = 0;
          }
        }
      }

      if (kakuteiValues.size() >= query2Size - 1) {
        // あいまいクエリ
        double minPenalty = m2_2_penaltySortedArray[0].first;
        for (int i = 0; i < m2_2_kouhoCount; ++i) {
          if (m2_2_penaltySortedArray[i].first > minPenalty + 180.0) { break; }
          int idx = m2_2_penaltySortedArray[i].second;
          for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
              if (m2_kakuteiValue1[j][k] != -1 || m2_aimaiValue[j][k] != -1) continue;
              nums[j][k][m2_kouho[idx][j][k]]++;
            }
          }
        }
        int x = -1, y = -1, minMax = INF;
        for (int j = 0; j < n; ++j) {
          for (int k = 0; k < n; ++k) {
            int ma = max(nums[j][k][0], max(nums[j][k][1], nums[j][k][2]));
            if (ma == 0) continue;
            if (ma < minMax) {
              x = j;
              y = k;
              minMax = ma;
            }
          }
        }

        if (x != -1) {
          vector<P> points;
          int valSum = 0;
          for (int i = 0; i < query2Size - 1; ++i) {
            points.push_back(kakuteiValues[i].first);
            valSum += kakuteiValues[i].second;
          }
          points.emplace_back(x, y);
          int val = Query2(points, ofs) - valSum;
          m2_aimaiValue[x][y] = val;
          continue;
        }
      }

      // 確定クエリ
      {
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
              nums[i][j][k] = 0;
            }
          }
        }
        double minPenalty = m2_2_penaltySortedArray[0].first;
        for (int i = 0; i < m2_2_kouhoCount; ++i) {
          if (m2_2_penaltySortedArray[i].first > minPenalty + 180.0) { break; }
          int idx = m2_2_penaltySortedArray[i].second;
          for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
              if (m2_kakuteiValue1[j][k] != -1) continue;
              nums[j][k][m2_kouho[idx][j][k]]++;
            }
          }
        }

        int x = -1, y = -1, minMax = INF;
        for (int j = 0; j < n; ++j) {
          for (int k = 0; k < n; ++k) {
            int ma = max(nums[j][k][0], max(nums[j][k][1], nums[j][k][2]));
            if (ma == 0) continue;
            if (ma < minMax) {
              x = j;
              y = k;
              minMax = ma;
            }
          }
        }

        if (x != -1) {
          int val = Query1(x, y, ofs);
          m2_kakuteiValue1[x][y] = val;
        }
      }
    }
  }

  if (!isAC) { Method2(ofs); }
}

int m3_used[MAX_N][MAX_N];
int m3_used2[MAX_N][MAX_N];
int m3_usedM[MAX_M];
int m3_findMCount;
void Method3_PrintAns(ofstream& ofs)
{
  vector<P> ans;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (m3_used[i][j] >= 1) {
        ans.emplace_back(i, j);
      }
      else if (m3_used2[i][j] >= 1) {
        ans.emplace_back(i, j);
      }
    }
  }
  bool isCorrect = QueryAns(ans, ofs);
}

vector<P> Method3_GetOnes()
{
  vector<P> ones;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (m3_used[i][j] >= 1) { ones.emplace_back(i, j); }
    }
  }
  return ones;
}

// {M番号、D番号}
struct Kouho
{
  int mNum;
  int dNum;
  int xPos;
  int yPos;
};
vector<Kouho> Method3_GetAllKouho()
{
  vector<Kouho> res;
  for (int i = 0; i < m; ++i) {
    if (m3_usedM[i]) continue;
    // チェック
    for (int j = 0; j < d[i]; ++j) {
      for (int xxx = 0; xxx < n; ++xxx) {
        for (int yyy = 0; yyy < n; ++yyy) {
          int isOk = 1;
          for (int k = 0; k < d[i]; ++k) {
            int xx = xxx + ii[i][k] - ii[i][j];
            int yy = yyy + jj[i][k] - jj[i][j];
            if (IsNG(xx, yy)) {
              isOk = 0;
              break;
            }
            if (m3_used[xx][yy] == 0) {
              isOk = 0;
              break;
            }
          }
          if (isOk) { res.emplace_back(Kouho({ i, j, xxx, yyy })); }
        }
      }
    }
  }
  return res;
}

// {M番号,D番号}
vector<P> Method3_GetAllKouhoSpecificPoint1(int xxx, int yyy)
{
  vector<P> res;
  for (int i = 0; i < m; ++i) {
    if (m3_usedM[i]) continue;
    // チェック
    for (int j = 0; j < d[i]; ++j) {
      int isOk = 1;
      for (int k = 0; k < d[i]; ++k) {
        int xx = xxx + ii[i][k] - ii[i][j];
        int yy = yyy + jj[i][k] - jj[i][j];
        if (IsNG(xx, yy)) {
          isOk = 0;
          break;
        }
        if (m3_used[xx][yy] == 0) {
          isOk = 0;
          break;
        }
      }
      if (!isOk) { continue; }
      res.emplace_back(i, j);
    }
  }
  return res;
}

// {M番号,D番号}
vector<P> Method3_GetAllKouhoSpecificPoint2(int xxx, int yyy)
{
  vector<P> res;
  for (int i = 0; i < m; ++i) {
    if (m3_usedM[i]) continue;
    // チェック
    const int shuuiSize = 3;
    for (int j = 0; j < d[i]; ++j) {
      int isOk = 1;
      int unused[shuuiSize][shuuiSize];
      for (int k = 0; k < shuuiSize; ++k) {
        for (int l = 0; l < shuuiSize; ++l) {
          unused[k][l] = 0;
          int xx = xxx + k - (shuuiSize / 2);
          int yy = yyy + l - (shuuiSize / 2);
          if (!IsNG(xx, yy) && m3_used[xx][yy] >= 1) { unused[k][l] = m3_used[xx][yy]; }
        }
      }
      for (int k = 0; k < d[i]; ++k) {
        int xx = xxx + ii[i][k] - ii[i][j];
        int yy = yyy + jj[i][k] - jj[i][j];
        if (IsNG(xx, yy)) {
          isOk = 0;
          break;
        }
        if (m3_used[xx][yy] == 0) {
          isOk = 0;
          break;
        }
        int xxxx = xx - xxx + (shuuiSize / 2);
        int yyyy = yy - yyy + (shuuiSize / 2);
        if (0 <= xxxx && xxxx < shuuiSize && 0 <= yyyy && yyyy < shuuiSize) { unused[xxxx][yyyy]--; }
      }
      if (!isOk) { continue; }

      // 周囲9マスで使っていないマスのチェック
      vector<P> unusedPoints;
      for (int k = 0; k < shuuiSize; ++k) {
        for (int l = 0; l < shuuiSize; ++l) {
          if (unused[k][l] > 0) { unusedPoints.emplace_back(xxx + k - (shuuiSize / 2), yyy + l - (shuuiSize / 2)); }
        }
      }

      for (auto unusedPoint : unusedPoints) {
        int ux = unusedPoint.first;
        int uy = unusedPoint.second;
        m3_usedM[i] = 1;
        for (int k = 0; k < d[i]; ++k) {
          int xx = xxx + ii[i][k] - ii[i][j];
          int yy = yyy + jj[i][k] - jj[i][j];
          if (m3_used[xx][yy] >= 1) { m3_used[xx][yy]--; }
        }
        auto kouho2 = Method3_GetAllKouhoSpecificPoint1(ux, uy);
        for (int k = 0; k < d[i]; ++k) {
          int xx = xxx + ii[i][k] - ii[i][j];
          int yy = yyy + jj[i][k] - jj[i][j];
          if (m3_used[xx][yy] >= 0) { m3_used[xx][yy]++; }
        }
        m3_usedM[i] = 0;

        if (kouho2.empty()) {
          isOk = 0;
          break;
        }
      }

      if (!isOk) { continue; }

      res.emplace_back(i, j);
    }
  }
  return res;
}

int onesEmptyCount = 0;
void Method3_Initialize()
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m3_used[i][j] = -1;
      m3_used2[i][j] = 0;
    }
  }
  for (int i = 0; i < m; ++i) {
    m3_usedM[i] = 0;
  }
  m3_findMCount = 0;
}

// Method3: 空白マスを探索してクエリを投げ、パターンをマッチングする

// パターンが一意に定まるかチェックし、定まる場合は確定させる
bool Method3_CheckUniquePattern(const vector<P>& ones, ofstream& ofs, int& cnt, int& cnt2)
{
  for (auto p : ones) {
    int x = p.first;
    int y = p.second;
    vector<P> allKouho;
    if (GetNowTime() < 1.0) {
      allKouho = Method3_GetAllKouhoSpecificPoint2(x, y);
    }
    else {
      allKouho = Method3_GetAllKouhoSpecificPoint1(x, y);
    }

    int realCount = 0;
    vector<P> v1;
    int onem = -1;

    // 候補パターンをチェック
    for (auto md : allKouho) {
      int i = md.first;
      int j = md.second;

      vector<P> v2;
      for (int k = 0; k < d[i]; ++k) {
        int xx = x + ii[i][k] - ii[i][j];
        int yy = y + jj[i][k] - jj[i][j];
        v2.emplace_back(xx, yy);
      }
      sort(v2.begin(), v2.end());

      if (realCount == 0) {
        v1 = v2;
        onem = i;
        realCount++;
      }
      else {
        if (v1 != v2) {
          realCount++;
          break;
        }
      }
    }

    // パターンが一意に定まる場合
    if (realCount == 1) {
      for (auto p : v1) {
        int xx = p.first;
        int yy = p.second;
        if (m3_used[xx][yy] > 0) {
          cnt--;
          cnt2++;
          m3_used[xx][yy]--;
          m3_used2[xx][yy]++;
        }
        else {
          cnt2++;
          m3_used2[xx][yy]++;
        }
      }
      m3_usedM[onem] = 1;
      m3_findMCount++;
      return true;
    }
  }
  return false;
}

// パターン候補から最適なクエリ位置を選択
void Method3_QueryBestPosition(const vector<P>& ones, ofstream& ofs, int& cnt, int& cnt2)
{
  int minKouhoCount = INF;
  int minNums[MAX_N][MAX_N];
  int minx = -1, miny = -1;

  // 各点での候補数を計算
  for (auto p : ones) {
    int x = p.first;
    int y = p.second;
    int kouhoCount = 0;
    int nums[MAX_N][MAX_N];

    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) nums[i][j] = 0;

    vector<P> allKouho;
    if (GetNowTime() < 1.0) {
      allKouho = Method3_GetAllKouhoSpecificPoint2(x, y);
    }
    else {
      allKouho = Method3_GetAllKouhoSpecificPoint1(x, y);
    }

    for (auto md : allKouho) {
      int i = md.first;
      int j = md.second;
      kouhoCount++;
      for (int k = 0; k < d[i]; ++k) {
        int xx = x + ii[i][k] - ii[i][j];
        int yy = y + jj[i][k] - jj[i][j];
        nums[xx][yy]++;
      }
    }

    // 最も候補数が少ない点を選択
    if (kouhoCount < minKouhoCount) {
      // 適切な探索点があるか確認
      int nx = -1;
      int ny = -1;
      int minDiff = INF;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        if (m3_used[i][j] != -1 || nums[i][j] == 0) continue;
        int diff = abs(kouhoCount / 2 - nums[i][j]);
        if (diff < minDiff) {
          minDiff = diff;
          nx = i;
          ny = j;
        }
      }
      if (nx == -1) continue;

      minKouhoCount = kouhoCount;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) minNums[i][j] = nums[i][j];
      minx = x;
      miny = y;
    }
  }

  if (minKouhoCount > 1) {
    // クエリ位置を決定
    int nx = -1;
    int ny = -1;
    int minDiff = INF;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
      if (m3_used[i][j] != -1 || minNums[i][j] == 0) continue;
      int diff = abs(minKouhoCount / 2 - minNums[i][j]);
      if (diff < minDiff) {
        minDiff = diff;
        nx = i;
        ny = j;
      }
    }

    // 適切な位置が見つからない場合はランダム選択
    if (nx == -1) {
      if (mode != 0) {
        ofs << "# random choice, kouhoCount = " << minKouhoCount << endl;
      }
      while (nx == -1) {
        int nnx = Rand() % n;
        int nny = Rand() % n;
        if (m3_used[nnx][nny] == -1) {
          nx = nnx;
          ny = nny;
        }
      }
    }

    // 現在のグリッド状態を出力（デバッグ用）
    if (mode != 0) {
      for (int i = 0; i < n; ++i) {
        ofs << "# ";
        for (int j = 0; j < n; ++j) {
          if (m3_used[i][j] == -1) {
            ofs << 'X';
          }
          else {
            ofs << m3_used[i][j];
          }
        }
        ofs << endl;
      }
    }

    int val = Query1(nx, ny, ofs);
    m3_used[nx][ny] = val - m3_used2[nx][ny];
    cnt += m3_used[nx][ny];
  }
}

// 空白マスからパターン候補を生成して最適な場所にクエリを投げる
void Method3_QueryEmptyCell(ofstream& ofs, int& cnt, int& cnt2)
{
  int nums[MAX_N][MAX_N];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nums[i][j] = 0;
    }
  }

  // 全ての候補パターンを取得
  auto allKouho = Method3_GetAllKouho();
  for (const auto& kouho : allKouho) {
    int i = kouho.mNum;
    int j = kouho.dNum;
    for (int k = 0; k < d[i]; ++k) {
      int xx = kouho.xPos + ii[i][k] - ii[i][j];
      int yy = kouho.yPos + jj[i][k] - jj[i][j];
      nums[xx][yy]++;
    }
  }

  // 最も多くの候補がある場所を選択
  int maxNums = 0;
  int x = -1, y = -1;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (nums[i][j] > maxNums && m3_used[i][j] == -1) {
        // 初期段階では離れた場所を優先
        if (m3_findMCount == 0 && m >= 10 && g_cost < 5) {
          int hasNeighbor = false;
          for (int k = -3; k < 4; ++k) {
            for (int l = -3; l < 4; ++l) {
              int ni = i + k;
              int nj = j + l;
              if (!IsNG(ni, nj) && m3_used[ni][nj] >= 1) {
                hasNeighbor = true;
                break;
              }
            }
          }
          if (hasNeighbor) continue;
        }
        maxNums = nums[i][j];
        x = i;
        y = j;
      }
    }
  }

  int val = Query1(x, y, ofs);
  m3_used[x][y] = val - m3_used2[x][y];
  cnt += m3_used[x][y];
}

// Method3: パターンマッチングによる油田探索
void Method3(ofstream& ofs)
{
  Method3_Initialize();
  onesEmptyCount = 0;

  int cnt = 0;   // 確定した油田セル数（m3_used）
  int cnt2 = 0;  // 確定した油田セル数（m3_used2）

  while (true) {
    // 全ての油田を発見した場合
    if (cnt + cnt2 == sum) {
      Method3_PrintAns(ofs);
      break;
    }

    vector<P> ones = Method3_GetOnes();

    // 油田があると分かっているマスがない場合
    if (ones.empty() || (m3_findMCount == 0 && m >= 5 && g_cost < 10)) {
      if (mode != 0) {
        onesEmptyCount++;
        ofs << "# ones.empty() : count = " << onesEmptyCount << endl;
      }
      Method3_QueryEmptyCell(ofs, cnt, cnt2);
    }
    // 油田があるマスからパターンを探索
    else {
      std::shuffle(ones.begin(), ones.end(), engine);

      // パターンが一意に定まるかチェック
      if (Method3_CheckUniquePattern(ones, ofs, cnt, cnt2)) {
        continue;
      }

      // 最適な位置にクエリを投げる
      Method3_QueryBestPosition(ones, ofs, cnt, cnt2);
    }
  }
}

int m4_values1[MAX_N][MAX_N];
int m4_values2[MAX_N][MAX_N];
vector<P> m4_kouhos[MAX_N][MAX_N];  // {M番号,D番号}
vector<int> m4_kouhos_OK[MAX_N][MAX_N];
void Method4(int kakuteiCount, ofstream& ofs)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m4_values1[i][j] = -1;
      m4_values2[i][j] = -1;
      m4_kouhos[i][j].clear();
    }
  }

  vector<pair<P, int>> kakuteiValues;

  // 全てのマスの候補を列挙する
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < m; ++k) {
        for (int l = 0; l < d[k]; ++l) {
          int miI = i - ii[k][l];
          int maI = i + (maxIi[k] - ii[k][l]);
          int miJ = j - jj[k][l];
          int maJ = j + (maxJj[k] - jj[k][l]);
          if (IsNG(miI, miJ) || IsNG(maI, maJ)) continue;
          m4_kouhos[i][j].emplace_back(k, l);
          m4_kouhos_OK[i][j].push_back(1);
        }
      }
    }
  }

  // 失敗した解答を覚えておく
  int wrongAnswerCount = 0;
  vector<vector<P>> wrongAnswers;

  while (true) {
    // 解答できるかチェック
    {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          if (m4_values1[i][j] <= 0 && m4_values2[i][j] <= 0) continue;
          for (int k = 0; k < m4_kouhos[i][j].size(); ++k) {
            if (m4_kouhos_OK[i][j][k] == 0) continue;
          }
        }
      }
    }
    // 確定マスの数が足りなければ確定操作を行う
  }
}

ll Solve(int probNum)
{
  g_startTime = clock();
  g_endTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  Initialize();

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 解答を出力
  ll para = hyperParams[mIndex][epsIndex];
  if (m == 2) {
    int methodNum = para % D12 / D10;
    if (methodNum == 0) {
      Method2(ofs);
    }
    else {
      int query2Size = para % D2;
      Method2_2(query2Size, ofs);
    }
  }
  else {
    Method3(ofs);
  }

  if (ofs.is_open()) { ofs.close(); }

  ll score = 0;
  if (mode != 0) { score = CalcScore(); }
  return score;
}

// hyperParamsをテキストファイルに出力
void PrintHyperParams(int loop)
{
  string str = "./hyperParams/hyperParams" + to_string(loop) + ".txt";
  std::ofstream file(str);

  // 出力したい文字列
  file << "ll hyperParams[21][21] = {" << endl;
  for (int _mIndex = 0; _mIndex < M_INDEX_SIZE; _mIndex++) {
    file << "{";
    for (int _epsIndex = 0; _epsIndex < EPS_INDEX_SIZE; _epsIndex++) {
      file << setw(19) << hyperParams[_mIndex][_epsIndex];
      if (_epsIndex < EPS_INDEX_SIZE - 1) { file << ","; }
    }
    if (_mIndex < M_INDEX_SIZE - 1) {
      file << "}," << endl;
    }
    else {
      file << "}" << endl;
    }
  }
  file << "};" << endl;
  file << endl;

  file.close();
}

void Battle()
{
  int loop = 0;  // ループ回数のカウンタ
  queue<pair<int, ll>> winQueue;
  while (true) {
    ll newPara = 0;
    if (!winQueue.empty()) {
      // winQueueに中身がある場合はそちらからシミュレーション
      // 次に試したい場所
      int nextIndex = winQueue.front().first;
      mIndex = nextIndex / D2;
      epsIndex = nextIndex % D2;
      // 試したいパラメータ
      newPara = winQueue.front().second;
      winQueue.pop();
    }
    else {
      // ランダムに1箇所選びシミュレーション
      // 1ケースジェネレート
      TestCaseGenerator(true, true, true);

      // 新しいハイパラ（挑戦者）を作成
      int solverNum = rand() % 2;
      if (solverNum == 0) {
        newPara = solverNum * D10;
      }
      else if (solverNum == 1) {
        int query2Size = Rand() % 9 + 1;
        newPara = solverNum * D10 + query2Size;
      }
    }

    // 同じ値の場合はcontinue
    if (hyperParams[mIndex][epsIndex] == newPara) continue;

    int win = 0;
    int lose = 0;
    int draw = 0;

    int judgeMode = 0;
    bool isWin = false;
    if (m != 2 || eps > 0.1) continue;
    cout << loop << ' ' << m << ' ' << eps << ", oldPara = " << hyperParams[mIndex][epsIndex] << ", newPara = " << newPara;
    if (judgeMode == 0) {  // 20戦18勝
      for (int _ = 0; _ < 20; _++) {
        // ここで毎回テストケースを1ケース作成
        TestCaseGenerator(true, false, false);
        ll oldScore = Solve(0);
        ll oldPara = hyperParams[mIndex][epsIndex];
        hyperParams[mIndex][epsIndex] = newPara;
        ll newScore = Solve(0);
        hyperParams[mIndex][epsIndex] = oldPara;

        if (newScore < oldScore) {
          // AHC030はスコアが小さいほど良い
          win++;
        }
        else if (newScore == oldScore) {
          draw++;
        }
        else {
          lose++;
        }

        // 時間短縮(結果が明らかなときはbreak)
        if (lose >= 3) break;
      }

      if (win >= 3 && (double)win / (win + lose) >= 0.9) {  // 1勝0敗19分などは避ける
        isWin = true;
      }
    }
    else if (judgeMode == 1) {  // 200戦120勝
      for (int _ = 0; _ < 200; _++) {
        // ここで毎回テストケースを1ケース作成
        TestCaseGenerator(true, false, false);
        ll oldScore = Solve(0);
        ll oldPara = hyperParams[mIndex][epsIndex];
        hyperParams[mIndex][epsIndex] = newPara;
        ll newScore = Solve(0);
        hyperParams[mIndex][epsIndex] = oldPara;

        if (newScore < oldScore) {
          // AHC030はスコアが小さいほど良い
          win++;
        }
        else if (newScore == oldScore) {
          draw++;
        }
        else {
          lose++;
        }

        // 時間短縮(結果が明らかなときはbreak)
        if (win >= 120) break;
        if (lose >= 80) break;
        if (win <= lose - 10) break;
        if (lose >= 20 && win <= lose) break;
        if (win - lose >= 40) break;
      }

      if (win >= 10 && (double)win / (win + lose) >= 0.6) {  // 3勝1敗116分などは避ける
        isWin = true;
      }
    }
    else if (judgeMode == 2) {  // 1000戦550勝
      for (int _ = 0; _ < 1000; _++) {
        // ここで毎回テストケースを1ケース作成
        TestCaseGenerator(true, false, false);
        ll oldScore = Solve(0);
        ll oldPara = hyperParams[mIndex][epsIndex];
        hyperParams[mIndex][epsIndex] = newPara;
        ll newScore = Solve(0);
        hyperParams[mIndex][epsIndex] = oldPara;

        if (newScore < oldScore) {
          // AHC030はスコアが小さいほど良い
          win++;
        }
        else if (newScore == oldScore) {
          draw++;
        }
        else {
          lose++;
        }

        // 時間短縮(結果が明らかなときはbreak)
        if (win >= 550) break;
        if (lose >= 450) break;
        if (win - lose >= 100) break;
        if (win <= lose - 10) break;
        if (lose >= 20 && win <= lose) break;
      }

      if (win >= 30 && (double)win / (win + lose) >= 0.55) { isWin = true; }
    }
    cout << " : win = " << win << ", lose = " << lose << ", draw = " << draw << endl;

    if (isWin) {
      cout << "win! loop = " << loop << " : " << endl;
      // hyperParamsを新パラメータで上書き
      hyperParams[mIndex][epsIndex] = newPara;

      // 付近(各インデックス±1)をキューに入れる
      for (int i = 0; i < 4; i++) {
        int nextMIndex = mIndex;
        int nextEpsIndex = epsIndex;
        if (i == 0) {
          nextMIndex--;
        }
        else if (i == 1) {
          nextMIndex++;
        }
        else if (i == 2) {
          nextEpsIndex--;
        }
        else if (i == 3) {
          nextEpsIndex++;
        }

        // 範囲チェック
        if (nextMIndex < 2 || 20 < nextMIndex) continue;
        if (nextEpsIndex <= 0 || 20 < nextEpsIndex) continue;

        // キューに入れる
        winQueue.push(make_pair(nextEpsIndex, newPara));
      }
    }

    loop++;
    if (loop % 100 == 0) {  // 適当なタイミングで出力
      PrintHyperParams(loop);
    }
  }
}

int main()
{
  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    ll onesEmptySum = 0;
    for (int i = 0; i < 10; ++i) {
      onesEmptyCount = 0;
      ll score = Solve(i);
      if (score == 0) continue;
      sum += score;
      onesEmptySum += onesEmptyCount * 1000000LL;
      cout << "num = " << setw(3) << i << ", ";
      cout << "score = " << setw(9) << score << ", ";
      cout << "sum = " << setw(11) << sum << ", ";
      cout << "ones.empty() : count = " << setw(3) << onesEmptyCount << ", ";
      cout << "onesEmptySum = " << setw(11) << onesEmptySum << endl;
    }
  }
  else if (mode == 2) {
    Battle();
  }

  return 0;
}
