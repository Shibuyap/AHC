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

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

const int INF = 1001001001;
double TL = 1.8;
int mode;

int n, m;
int si, sj;
string S[20];
int a[20][20];
string t[210];
int b[210][5];
vector<P> v[26];
int vcnt[26];

int ans[210];
int id[210];
vector<P> answer;
int real_ans[210];
int real_id[210];
int real_score;

struct Path
{
  int si;
  int sj;
  int ti;
  int tj;
  int val;
  vector<P> path;

  bool operator<(const Path& right) const { return val < right.val; }
};
vector<Path> V[210];
int baseVal[210];

int Distance(int i1, int j1, int i2, int j2)
{
  return abs(i1 - i2) + abs(j1 - j2);
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  rep(i, 26) { v[i].clear(); }
  rep(i, 210) { V[i].clear(); }
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
    cin >> si >> sj;
    rep(i, n) cin >> S[i];
    rep(i, m) cin >> t[i];
  }
  // ファイル入力する
  else {
    ifs >> n >> m;
    ifs >> si >> sj;
    rep(i, n) ifs >> S[i];
    rep(i, m) ifs >> t[i];
  }

  rep(i, n)
  {
    rep(j, n) { a[i][j] = S[i][j] - 'A'; }
  }
  rep(i, m)
  {
    rep(j, 5) { b[i][j] = t[i][j] - 'A'; }
  }

  rep(i, n)
  {
    rep(j, n) { v[a[i][j]].push_back(P(i, j)); }
  }
  rep(i, 26) { vcnt[i] = v[i].size(); }
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

int Same3(const Path& x, const Path& y)
{
  if (x.path[2] == y.path[0] && x.path[3] == y.path[1] && x.path[4] == y.path[2]) {
    int res = Distance(x.path[4].first, x.path[4].second, y.path[0].first, y.path[0].second)
      + Distance(y.path[0].first, y.path[0].second, y.path[1].first, y.path[1].second)
      + Distance(y.path[1].first, y.path[1].second, y.path[2].first, y.path[2].second) + 3;
    return res;
  }
  return 0;
}

int Same2(const Path& x, const Path& y)
{
  if (x.path[3] == y.path[0] && x.path[4] == y.path[1]) {
    int res = Distance(x.path[4].first, x.path[4].second, y.path[0].first, y.path[0].second)
      + Distance(y.path[0].first, y.path[0].second, y.path[1].first, y.path[1].second) + 2;
    return res;
  }
  return 0;
}

// スコア計算
ll CalcScore()
{
  int score = 10000 - m * 5;
  rep(i, m)
  {
    if (i == 0) {
      score -= Distance(V[ans[i]][id[i]].si, V[ans[i]][id[i]].sj, si, sj);
    }
    else {
      int dist =
        Distance(V[ans[i]][id[i]].si, V[ans[i]][id[i]].sj,
          V[ans[i - 1]][id[i - 1]].ti, V[ans[i - 1]][id[i - 1]].tj);
      score -= dist;
      int same3 = Same3(V[ans[i - 1]][id[i - 1]], V[ans[i]][id[i]]);
      int same2 = Same2(V[ans[i - 1]][id[i - 1]], V[ans[i]][id[i]]);
      if (same3 > 0) {
        score += same3;
      }
      else if (same2 > 0) {
        score += same2;
      }
      else {
        if (dist == 0) score++;
      }
    }

    score -= V[ans[i]][id[i]].val;
  }
  return score;
}

// 初期解生成
void Initialize()
{
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());

  real_score = -1;

  rep(aespa, 400)
  {
    // 貪欲に作る
    int f[210] = {};
    int ti = si;
    int tj = sj;

    vector<int> nums;
    rep(j, m)nums.push_back(j);
    rep(i, m)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, m)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, 20)
        {
          if (V[j].size() <= k)break;
          int dist1 = Distance(V[j][k].si, V[j][k].sj, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = Same3(V[ans[i - 1]][id[i - 1]], V[j][k]);
            same2 = Same2(V[ans[i - 1]][id[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].val - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].val - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].val - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && Rand() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      ans[i] = miAns;
      id[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].ti;
      tj = V[miAns][miId].tj;
    }

    int score = CalcScore();
    if (score > real_score) {
      real_score = score;
      rep(i, m)
      {
        real_ans[i] = ans[i];
        real_id[i] = id[i];
      }
    }
  }

  rep(i, m)
  {
    ans[i] = real_ans[i];
    id[i] = real_id[i];
  }


  rep(aespa, 400)
  {
    // 貪欲に作る
    int f[210] = {};
    int ti = si;
    int tj = sj;

    vector<int> nums;
    rep(j, m)nums.push_back(j);
    rep(i, 150)
    {
      f[ans[i]] = 1;
    }
    srep(i, 150, m)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, m)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, 20)
        {
          if (V[j].size() <= k)break;
          int dist1 = Distance(V[j][k].si, V[j][k].sj, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = Same3(V[ans[i - 1]][id[i - 1]], V[j][k]);
            same2 = Same2(V[ans[i - 1]][id[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].val - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].val - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].val - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && Rand() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      ans[i] = miAns;
      id[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].ti;
      tj = V[miAns][miId].tj;
    }

    int score = CalcScore();
    if (score > real_score) {
      real_score = score;
      rep(i, m)
      {
        real_ans[i] = ans[i];
        real_id[i] = id[i];
      }
    }
  }

  rep(i, m)
  {
    ans[i] = real_ans[i];
    id[i] = real_id[i];
  }

  rep(aespa, 1000)
  {
    // 貪欲に作る
    int f[210] = {};
    int ti = si;
    int tj = sj;

    vector<int> nums;
    rep(j, m)nums.push_back(j);
    rep(i, 180)
    {
      f[ans[i]] = 1;
    }
    srep(i, 180, m)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, m)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, 20)
        {
          if (V[j].size() <= k)break;
          int dist1 = Distance(V[j][k].si, V[j][k].sj, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = Same3(V[ans[i - 1]][id[i - 1]], V[j][k]);
            same2 = Same2(V[ans[i - 1]][id[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].val - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].val - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].val - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && Rand() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      ans[i] = miAns;
      id[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].ti;
      tj = V[miAns][miId].tj;
    }

    int score = CalcScore();
    if (score > real_score) {
      real_score = score;
      rep(i, m)
      {
        real_ans[i] = ans[i];
        real_id[i] = id[i];
      }
    }
  }

  rep(i, m)
  {
    ans[i] = real_ans[i];
    id[i] = real_id[i];
  }

  rep(aespa, 1000)
  {
    // 貪欲に作る
    int f[210] = {};
    int ti = si;
    int tj = sj;

    vector<int> nums;
    rep(j, m)nums.push_back(j);
    rep(i, 190)
    {
      f[ans[i]] = 1;
    }
    srep(i, 190, m)
    {
      int mi = INF;
      int miAns = -1;
      int miId = -1;

      std::shuffle(nums.begin(), nums.end(), engine);
      rep(jj, m)
      {
        int j = nums[jj];
        if (f[j])continue;
        rep(k, 20)
        {
          if (V[j].size() <= k)break;
          int dist1 = Distance(V[j][k].si, V[j][k].sj, ti, tj);
          int same3 = 0;
          int same2 = 0;
          if (i > 0) {
            same3 = Same3(V[ans[i - 1]][id[i - 1]], V[j][k]);
            same2 = Same2(V[ans[i - 1]][id[i - 1]], V[j][k]);
          }
          int tmp = 0;
          if (same3 > 0) {
            tmp = dist1 - same3 + V[j][k].val - baseVal[j];
          }
          else if (same2 > 0) {
            tmp = dist1 - same2 + V[j][k].val - baseVal[j];
          }
          else {
            int same = 0;
            if (i > 0 && dist1 == 0) same++;
            tmp = dist1 - same + V[j][k].val - baseVal[j];
          }

          if (tmp < mi) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
          else if (tmp == mi && Rand() % 2) {
            mi = tmp;
            miAns = j;
            miId = k;
          }
        }
      }

      ans[i] = miAns;
      id[i] = miId;
      f[miAns] = 1;
      ti = V[miAns][miId].ti;
      tj = V[miAns][miId].tj;
    }

    int score = CalcScore();
    if (score > real_score) {
      real_score = score;
      rep(i, m)
      {
        real_ans[i] = ans[i];
        real_id[i] = id[i];
      }
    }
  }

  rep(i, m)
  {
    ans[i] = real_ans[i];
    id[i] = real_id[i];
  }


  //cout << real_score << endl;


  //// 先頭から順番に一番短いやつを入れていく
  //rep(i, m) {
  //  int mi = INF;
  //  int mid = INF;
  //  ans[i] = i;
  //  rep(j, V[i].size()) {
  //    if (V[i][j].val < mi) {
  //      mi = V[i][j].val;
  //      id[i] = j;
  //      mid = abs(V[i][j].si - 7) + abs(V[i][j].sj - 7) + abs(V[i][j].ti - 7) + abs(V[i][j].tj - 7);
  //    }
  //    else if (V[i][j].val == mi) {
  //      int tmp = abs(V[i][j].si - 7) + abs(V[i][j].sj - 7) + abs(V[i][j].ti - 7) + abs(V[i][j].tj - 7);
  //      if (tmp < mid) {
  //        mi = V[i][j].val;
  //        id[i] = j;
  //        mid = tmp;
  //      }
  //    }
  //  }
  //}
}

// 解答出力
void Output(ofstream& ofs)
{
  answer.clear();
  rep(i, m)
  {
    Path path = V[ans[i]][id[i]];
    rep(j, 5)
    {
      if (j <= 2 && i > 0 && Same3(V[ans[i - 1]][id[i - 1]], path) > 0) {
        //cout << answer.size() << endl;
        continue;
      }
      else if (j <= 1 && i > 0 && Same2(V[ans[i - 1]][id[i - 1]], path) > 0) {
        //cout << answer.size() << endl;
        continue;
      }
      else if (j == 0 && !answer.empty()) {
        if (answer.back() == path.path[0]) {
          continue;
        }
      }
      answer.push_back(path.path[j]);
    }
  }
  if (mode == 0) {
    for (auto p : answer) {
      cout << p.first << ' ' << p.second << endl;
    }
  }
  else {
    for (auto p : answer) {
      ofs << p.first << ' ' << p.second << endl;
    }
  }
}

void Method11(double temperature)
{
  // 2点スワップ
  int x1 = Rand() % m;
  int x2 = Rand() % m;
  if (x1 > x2) swap(x1, x2);
  if (x2 - x1 <= 1) return;
  int x11 = x1 - 1;
  int x12 = x1 + 1;
  int x21 = x2 - 1;
  int x22 = x2 + 1;

  int beforeScore = 0;
  int afterScore = 0;
  if (x1 != 0 && x2 != m - 1) {
    {
      int dist1 = Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj,
        V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
      int dist2 = Distance(V[ans[x1]][id[x1]].ti, V[ans[x1]][id[x1]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);
      int dist3 = Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int dist4 = Distance(V[ans[x2]][id[x2]].ti, V[ans[x2]][id[x2]].tj,
        V[ans[x22]][id[x22]].si, V[ans[x22]][id[x22]].sj);
      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 + dist4 - same;
    }
    {
      int dist1 = Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int dist2 = Distance(V[ans[x1]][id[x1]].ti, V[ans[x1]][id[x1]].tj,
        V[ans[x22]][id[x22]].si, V[ans[x22]][id[x22]].sj);
      int dist3 = Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj,
        V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
      int dist4 = Distance(V[ans[x2]][id[x2]].ti, V[ans[x2]][id[x2]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist2 + dist3 + dist4 - same;
    }
  }
  else if (x1 == 0 && x2 != m - 1) {
    {
      int dist1 =
        Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj, si, sj);
      int dist2 = Distance(V[ans[x1]][id[x1]].ti, V[ans[x1]][id[x1]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);
      int dist3 = Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int dist4 = Distance(V[ans[x2]][id[x2]].ti, V[ans[x2]][id[x2]].tj,
        V[ans[x22]][id[x22]].si, V[ans[x22]][id[x22]].sj);
      int same = 0;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 + dist4 - same;
    }
    {
      int dist1 = Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int dist2 = Distance(V[ans[x1]][id[x1]].ti, V[ans[x1]][id[x1]].tj,
        V[ans[x22]][id[x22]].si, V[ans[x22]][id[x22]].sj);
      int dist3 =
        Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj, si, sj);
      int dist4 = Distance(V[ans[x2]][id[x2]].ti, V[ans[x2]][id[x2]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist2 + dist3 + dist4 - same;
    }
  }
  else if (x1 != 0 && x2 == m - 1) {
    {
      int dist1 = Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj,
        V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
      int dist2 = Distance(V[ans[x1]][id[x1]].ti, V[ans[x1]][id[x1]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);
      int dist3 = Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int same = 0;
      if (dist1 == 0) same++;
      if (dist2 == 0) same++;
      if (dist3 == 0) same++;
      beforeScore = dist1 + dist2 + dist3 - same;
    }
    {
      int dist1 = Distance(V[ans[x1]][id[x1]].si, V[ans[x1]][id[x1]].sj,
        V[ans[x21]][id[x21]].ti, V[ans[x21]][id[x21]].tj);
      int dist3 = Distance(V[ans[x2]][id[x2]].si, V[ans[x2]][id[x2]].sj,
        V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
      int dist4 = Distance(V[ans[x2]][id[x2]].ti, V[ans[x2]][id[x2]].tj,
        V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);

      int same = 0;
      if (dist1 == 0) same++;
      if (dist3 == 0) same++;
      if (dist4 == 0) same++;
      afterScore = dist1 + dist3 + dist4 - same;
    }
  }
  else {
    return;
  }

  int diffScore = beforeScore - afterScore;
  double prob = exp((double)diffScore / temperature);
  //if (prob > Rand01()) {
  if (diffScore >= 0) {
    swap(ans[x1], ans[x2]);
    swap(id[x1], id[x2]);
  }
}

// ID変更
void Method12(double temperature)
{
  int x = Rand() % m;
  int y = Rand() % 10;
  if (V[ans[x]].size() <= y) return;
  if (x == 0 || x == m - 1) return;
  if (y == id[x]) return;
  int x11 = x - 1;
  int x12 = x + 1;

  int beforeScore = 0;
  int afterScore = 0;
  {
    int dist1 = Distance(V[ans[x]][id[x]].si, V[ans[x]][id[x]].sj,
      V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
    int dist2 = Distance(V[ans[x]][id[x]].ti, V[ans[x]][id[x]].tj,
      V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    beforeScore = dist1 + dist2 - same + V[ans[x]][id[x]].val;
  }
  {
    int dist1 = Distance(V[ans[x]][y].si, V[ans[x]][y].sj,
      V[ans[x11]][id[x11]].ti, V[ans[x11]][id[x11]].tj);
    int dist2 = Distance(V[ans[x]][y].ti, V[ans[x]][y].tj,
      V[ans[x12]][id[x12]].si, V[ans[x12]][id[x12]].sj);
    int same = 0;
    if (dist1 == 0) same++;
    if (dist2 == 0) same++;
    afterScore = dist1 + dist2 - same + V[ans[x]][y].val;
  }

  int diffScore = beforeScore - afterScore;
  double prob = exp((double)diffScore / temperature);
  //if (prob > Rand01()) {
  if (diffScore >= 0) {
    id[x] = y;
    // cout << y << ' ' << diffScore << endl;
  }
}

void Method1(clock_t start)
{
  clock_t end;
  end = clock();

  int loop = 0;
  double startTemperature = 2;
  double endTemperature = 0;
  double nowProgress = 0;
  while (false) {
    loop++;
    if (loop % 100 == 0) {
      end = clock();
      double nowTime = (double)(end - start) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
      if (nowProgress > 1.0) break;
    }

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    int ra = Rand() % 100;
    if (ra < 50) {
      // 2点スワップ
      Method11(temperature);
    }
    else {
      Method12(temperature);
    }
  }

  if (mode != 0) {
    cout << "loop = " << loop << endl;
  }
}

ll Solve(int probNum)
{
  clock_t start, end;
  start = clock();
  end = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // dp
  {
    int dp[5][40][40];
    int dp2[5][40][40];
    rep(i, m)
    {
      rep(j, 5)
      {
        int x = b[i][0];
        if (j == 0) {
          rep(k, vcnt[x])
          {
            rep(l, vcnt[x]) { dp[j][k][l] = INF; }
          }
        }
        else {
          int y = b[i][j];
          rep(k, vcnt[x])
          {
            rep(l, vcnt[y]) { dp[j][k][l] = INF; }
          }
        }
      }

      int x = b[i][0];
      rep(j, 5)
      {
        if (j == 0) {
          rep(k, vcnt[x])
          {
            rep(l, vcnt[x])
            {
              if (k == l)
                dp[j][k][l] = 0;
              else
                dp[j][k][l] = INF;
            }
          }
        }
        else {
          int y = b[i][j - 1];
          int z = b[i][j];
          rep(k, vcnt[x])
          {
            rep(l, vcnt[y])
            {
              rep(o, vcnt[z])
              {
                int dist = abs(v[z][o].first - v[y][l].first) +
                  abs(v[z][o].second - v[y][l].second);
                if (dp[j][k][o] > dp[j - 1][k][l] + dist) {
                  dp[j][k][o] = dp[j - 1][k][l] + dist;
                  dp2[j][k][o] = l;
                }
              }
            }
          }
        }
      }

      int z = b[i][4];
      rep(j, vcnt[x])
      {
        rep(k, vcnt[z])
        {
          Path path;
          path.si = v[x][j].first;
          path.sj = v[x][j].second;
          path.ti = v[z][k].first;
          path.tj = v[z][k].second;
          path.val = dp[4][j][k];
          vector<P> tmp;
          int now = k;
          drep(l, 5)
          {
            int y = b[i][l];
            tmp.push_back(v[y][now]);
            if (l == 0) break;
            now = dp2[l][j][now];
          }
          reverse(tmp.begin(), tmp.end());
          path.path = tmp;
          V[i].push_back(path);
        }
      }

      sort(V[i].begin(), V[i].end());
      baseVal[i] = V[i][0].val;
    }
  }


  //rep(i, m) {
  //  rep(j, m) {
  //    if (b[i][2] == b[j][0] && b[i][3] == b[j][1] && b[i][4] == b[j][2]) {
  //      cout << t[i] << ' ' << t[j] << endl;
  //    }
  //  }
  //}

  // 初期解生成
  Initialize();

  Method1(start);

  // 解答を出力
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
    Rand();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
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
