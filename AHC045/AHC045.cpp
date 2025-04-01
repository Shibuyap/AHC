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

// ループの簡略化マクロ
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t RandXor()
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

// 0以上1未満の実数を返す乱数関数
static double Rand01() { return (RandXor() + 0.5) * (1.0 / UINT_MAX); }

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化（使用されていない）
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

namespace
{
  const int MAX_UF = 10005;
  int UF_par[MAX_UF];   // 親
  int UF_rank[MAX_UF];  // 木の深さ
  int UF_cnt[MAX_UF];   // 属する頂点の個数(親のみ正しい)

  // 初期化
  void UF_Init(int n)
  {
    for (int i = 0; i < n; i++) {
      UF_par[i] = i;
      UF_rank[i] = 0;
      UF_cnt[i] = 1;
    }
  }

  // 初期化
  void UF_Init()
  {
    UF_Init(MAX_UF);
  }

  // 木の根を求める
  int UF_Find(int x)
  {
    if (UF_par[x] == x) {
      return x;
    }
    else {
      return UF_par[x] = UF_Find(UF_par[x]);
    }
  }

  // xとyの属する集合を併合
  void UF_Unite(int x, int y)
  {
    x = UF_Find(x);
    y = UF_Find(y);
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

  // xとyが同じ集合に属するか否か
  bool UF_Same(int x, int y)
  {
    return UF_Find(x) == UF_Find(y);
  }

  // xの属する集合のサイズ
  int UF_Count(int x)
  {
    return UF_cnt[UF_Find(x)];
  }
}  // namespace



// 非常に大きな値
const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8; // 時間制限（Time Limit）
int mode;        // 実行モード
std::chrono::steady_clock::time_point startTimeClock; // 時間計測用

// 時間計測をリセットする関数
void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime()
{
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

ofstream ofs;

const int n = 800;
const int MAX_M = 400;
const int q = 400;
const int MAX_L = 15;

int m, l, w;

int g[MAX_M];
int sort_g[MAX_M];
int arg_g[MAX_M];
int arg_rev_g[MAX_M];

int lx[n], rx[n], ly[n], ry[n];
int pred_x[n], pred_y[n];

int queryCount;
vector<int> queries[q];
vector<P> queryAnswers[q];
vector<int> num_queries[n];
set<vector<int>> querySet;

int true_x[n], true_y[n];

ll ansScore;
double dScore;
int ans[n];
vector<int> ans_nums[MAX_M];
vector<P> ans_edges[MAX_M];

ll best_ansScore;
double best_dScore;
int best_ans[n];
vector<int> best_ans_nums[MAX_M];
vector<P> best_ans_edges[MAX_M];

struct Edge
{
  int dist;
  int u;
  int v;

  bool operator<(const Edge& other) const
  {
    return dist < other.dist;
  }
};

P GetMIddlePoint(int num)
{
  return P((lx[num] + rx[num]) / 2, (ly[num] + ry[num]) / 2);
}

P GetPoint(int num)
{
  return P(pred_x[num], pred_y[num]);
}

P GetTruePoint(int num)
{
  return P(true_x[num], true_y[num]);
}

int Distance(const P& p1, const P& p2)
{
  return (p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second);
}

int Distance(int num1, int num2)
{
  return Distance(GetPoint(num1), GetPoint(num2));
}

double buildMST_sum;
P buildMST_points[n];
Edge buildMST_edges[n * n];
int buildMST_idxDict[n];
vector<P> BuildMST(const vector<int>& nums, bool isTrue = false)
{
  buildMST_sum = 0;

  UF_Init(nums.size());

  rep(i, nums.size())
  {
    if (isTrue) {
      buildMST_points[i] = GetTruePoint(nums[i]);
    }
    else {
      buildMST_points[i] = GetPoint(nums[i]);
    }
  }

  int edgeCount = 0;
  rep(i, nums.size())
  {
    srep(j, i + 1, nums.size())
    {
      Edge e;
      buildMST_edges[edgeCount].dist = Distance(buildMST_points[i], buildMST_points[j]);
      buildMST_edges[edgeCount].u = i;
      buildMST_edges[edgeCount].v = j;
      edgeCount++;
    }
  }

  sort(buildMST_edges, buildMST_edges + edgeCount);

  vector<P> res(nums.size() - 1);
  int resCount = 0;

  rep(i, edgeCount)
  {
    if (!UF_Same(buildMST_edges[i].u, buildMST_edges[i].v)) {
      buildMST_sum += sqrt(buildMST_edges[i].dist);
      UF_Unite(buildMST_edges[i].u, buildMST_edges[i].v);
      res[resCount] = P(nums[buildMST_edges[i].u], nums[buildMST_edges[i].v]);
      resCount++;
      if (UF_Count(buildMST_edges[i].u) == nums.size()) {
        break;
      }
    }
  }

  return res;
}

vector<P> BuildMST(const int g_num1, const int g_num2, bool isTrue = false)
{
  const auto& nums1 = ans_nums[g_num1];
  const auto& nums2 = ans_nums[g_num2];
  int n1 = nums1.size();
  int n2 = nums2.size();

  buildMST_sum = 0;

  UF_Init(n1 + n2);

  rep(i, n1)
  {
    if (isTrue) {
      buildMST_points[i] = GetTruePoint(nums1[i]);
    }
    else {
      buildMST_points[i] = GetPoint(nums1[i]);
    }

    buildMST_idxDict[nums1[i]] = i;
  }
  rep(i, n2)
  {
    if (isTrue) {
      buildMST_points[i + n1] = GetTruePoint(nums2[i]);
    }
    else {
      buildMST_points[i + n1] = GetPoint(nums2[i]);
    }

    buildMST_idxDict[nums2[i]] = n1 + i;
  }

  int maxDist = -1;
  int edgeCount = 0;
  for (auto p : ans_edges[g_num1]) {
    buildMST_edges[edgeCount].dist = Distance(buildMST_points[buildMST_idxDict[p.first]], buildMST_points[buildMST_idxDict[p.second]]);
    maxDist = max(maxDist, buildMST_edges[edgeCount].dist);
    buildMST_edges[edgeCount].u = buildMST_idxDict[p.first];
    buildMST_edges[edgeCount].v = buildMST_idxDict[p.second];
    edgeCount++;
  }
  for (auto p : ans_edges[g_num2]) {
    buildMST_edges[edgeCount].dist = Distance(buildMST_points[buildMST_idxDict[p.first]], buildMST_points[buildMST_idxDict[p.second]]);
    maxDist = max(maxDist, buildMST_edges[edgeCount].dist);
    buildMST_edges[edgeCount].u = buildMST_idxDict[p.first];
    buildMST_edges[edgeCount].v = buildMST_idxDict[p.second];
    edgeCount++;
  }

  Edge minE;
  minE.dist = INT_INF;
  bool contains = false;
  rep(ii, nums1.size())
  {
    rep(jj, nums2.size())
    {
      int i = ii;
      int j = jj + n1;
      buildMST_edges[edgeCount].dist = Distance(buildMST_points[i], buildMST_points[j]);
      buildMST_edges[edgeCount].u = i;
      buildMST_edges[edgeCount].v = j;
      if (buildMST_edges[edgeCount].dist < minE.dist) {
        minE = buildMST_edges[edgeCount];
      }
      if (buildMST_edges[edgeCount].dist <= maxDist) {
        edgeCount++;
        contains = true;
      }
    }
  }

  if (!contains) {
    buildMST_edges[edgeCount] = minE;
    edgeCount++;
  }

  sort(buildMST_edges, buildMST_edges + edgeCount);

  vector<P> res(n1 + n2 - 1);
  int resCount = 0;

  rep(i, edgeCount)
  {
    if (!UF_Same(buildMST_edges[i].u, buildMST_edges[i].v)) {
      buildMST_sum += sqrt(buildMST_edges[i].dist);
      UF_Unite(buildMST_edges[i].u, buildMST_edges[i].v);
      int uu, vv;
      if (buildMST_edges[i].u < n1) {
        uu = nums1[buildMST_edges[i].u];
      }
      else {
        uu = nums2[buildMST_edges[i].u - n1];
      }
      if (buildMST_edges[i].v < n1) {
        vv = nums1[buildMST_edges[i].v];
      }
      else {
        vv = nums2[buildMST_edges[i].v - n1];
      }
      res[resCount] = P(uu, vv);
      resCount++;
      if (UF_Count(buildMST_edges[i].u) == n1 + n2) {
        break;
      }
    }
  }

  return res;
}

void CopyToBest()
{
  best_ansScore = ansScore;
  best_dScore = dScore;
  rep(i, n)
  {
    best_ans[i] = ans[i];
  }
  rep(i, m)
  {
    best_ans_nums[i] = ans_nums[i];
    best_ans_edges[i] = ans_edges[i];
  }
}

void CopyToAns()
{
  ansScore = best_ansScore;
  dScore = best_dScore;
  rep(i, n)
  {
    ans[i] = best_ans[i];
  }
  rep(i, m)
  {
    ans_nums[i] = best_ans_nums[i];
    ans_edges[i] = best_ans_edges[i];
  }
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = INT_INF;
  rep(i, MAX_M)
  {
    ans_nums[i].clear();
    ans_edges[i].clear();
  }
  CopyToBest();

  queryCount = 0;
  rep(i, q)
  {
    queries[i].clear();
    queryAnswers[i].clear();
  }
  rep(i, n)
  {
    num_queries[i].clear();
  }
  querySet.clear();
}

bool isSimulateTruePoint = false;
// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int _n, _q;
    cin >> _n >> m >> _q >> l >> w;
    rep(i, m)cin >> g[i];
    rep(i, n)
    {
      cin >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
  }
  else {
    // ファイル入力
    int _n, _q;
    ifs >> _n >> m >> _q >> l >> w;
    rep(i, m)ifs >> g[i];
    rep(i, n)
    {
      ifs >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
    rep(i, n)ifs >> true_x[i] >> true_y[i];
  }

  rep(i, n)
  {
    pred_x[i] = GetMIddlePoint(i).first;
    pred_y[i] = GetMIddlePoint(i).second;
  }

  vector<P> vp;
  rep(i, m)
  {
    vp.emplace_back(g[i], i);
  }
  sort(vp.begin(), vp.end());
  rep(i, m)
  {
    sort_g[i] = vp[i].first;
    arg_g[i] = vp[i].second;
    arg_rev_g[vp[i].second] = i;
  }

  if (isSimulateTruePoint) {
    rep(i, n)
    {
      lx[i] = true_x[i];
      rx[i] = true_x[i];
      ly[i] = true_y[i];
      ry[i] = true_y[i];
    }
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (ofs.is_open()) {
    ofs.close();
  }

  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore(int g_num)
{
  int res = 0;

  for (auto e : ans_edges[g_num]) {
    res += sqrt(Distance(GetPoint(e.first), GetPoint(e.second)));
  }

  return res;
}

int CalcScoreAll()
{
  int res = 0;

  rep(i, m)
  {
    res += CalcScore(i);
  }

  return res;
}

int CalcScoreLocal()
{
  int res = 0;

  rep(i, m)
  {
    for (auto e : ans_edges[i]) {
      res += sqrt(Distance(P(true_x[e.first], true_y[e.first]), P(true_x[e.second], true_y[e.second])));
    }
  }

  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    cout << '!' << endl;
    rep(i, m)
    {
      int rev_i = arg_rev_g[i];
      for (auto po : ans_nums[rev_i]) {
        cout << po << ' ';
      }
      cout << endl;
      for (auto edge : ans_edges[rev_i]) {
        cout << edge.first << ' ' << edge.second << endl;
      }
    }
  }
  else {
    // ファイル出力
    ofs << '!' << endl;
    rep(i, m)
    {
      int rev_i = arg_rev_g[i];
      for (auto po : ans_nums[rev_i]) {
        ofs << po << ' ';
      }
      ofs << endl;
      for (auto edge : ans_edges[rev_i]) {
        ofs << edge.first << ' ' << edge.second << endl;
      }
    }
  }
}

void Query()
{
  if (mode == 0) {
    cout << "? " << l;
    rep(j, l)
    {
      cout << ' ' << queries[queryCount][j];
    }
    cout << endl;
    fflush(stdout);

    rep(j, l - 1)
    {
      int a, b;
      cin >> a >> b;
      queryAnswers[queryCount].push_back(P(a, b));
    }
  }
  else {
    ofs << "? " << l;
    rep(j, l)
    {
      ofs << ' ' << queries[queryCount][j];
    }
    ofs << endl;

    queryAnswers[queryCount] = BuildMST(queries[queryCount], true);
  }

  queryCount++;
}

// 面積大きい順にQ個クエリを投げる
// Lは最大まで使う
// 残りL-1頂点は近い順
void Method1_Query(int start, int queryEnd)
{
  vector<P> vp;
  rep(i, n)
  {
    vp.emplace_back((rx[i] - lx[i] + 1) * (ry[i] - ly[i] + 1), i);
  }
  sort(vp.begin(), vp.end(), greater<P>());

  srep(i, start, n)
  {
    int num = vp[i].second;
    auto po = GetPoint(num);
    vector<P> dists;
    rep(j, n)
    {
      dists.emplace_back(Distance(po, GetPoint(j)), j);
    }
    sort(dists.begin(), dists.end());
    queries[queryCount].clear();
    rep(j, l)
    {
      queries[queryCount].push_back(dists[j].second);
    }
    sort(queries[queryCount].begin(), queries[queryCount].end());
    if (querySet.find(queries[queryCount]) != querySet.end()) {
      continue;
    }
    rep(j, l)
    {
      num_queries[dists[j].second].push_back(queryCount);
    }
    querySet.insert(queries[queryCount]);
    Query();
    if (queryCount == queryEnd) {
      if (mode != 0) {
        cout << i + 1 << endl;
      }
      break;
    }
  }
}

// ナイーブな解法
void Method1()
{
  int now = 0;
  rep(i, m)
  {
    rep(j, sort_g[i])
    {
      ans[now] = i;
      ans_nums[i].push_back(now);
      now++;
    }
  }

  rep(i, m)
  {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }

  ansScore = CalcScoreAll();
  CopyToBest();
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition[10];
};

int BoundingScore(int g_num)
{
  int minx = INT_INF;
  int maxx = -INT_INF;
  int miny = INT_INF;
  int maxy = -INT_INF;
  for (auto num : ans_nums[g_num]) {
    auto po = GetPoint(num);
    minx = min(minx, po.first);
    maxx = max(maxx, po.first);
    miny = min(miny, po.second);
    maxy = max(maxy, po.second);
  }

  return (maxx - minx) * (maxx - minx) + (maxy - miny) * (maxy - miny);
}

void SimulatedAnnealing0(Hypers hypers, double timeLimit)
{
  double startTime = GetNowTime();

  dScore = 0;
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    double tmpScore = dScore;

    // 近傍解作成
    int raMode = RandXor() % hypers.Partition[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      if (m == 1)continue;
      ra1 = RandXor() % n;
      while (num_queries[ra1].empty()) {
        ra1 = RandXor() % n;
      }
      if (RandXor() % 2 == 0) {
        ra2 = RandXor() % (rx[ra1] - lx[ra1] + 1) + lx[ra1];
        ra3 = RandXor() % (ry[ra1] - ly[ra1] + 1) + ly[ra1];
      }
      else {
        int len = (rx[ra1] - lx[ra1] + 1) / 10 + 1;
        ra2 = RandXor() % (len * 2 + 1) - len + pred_x[ra1];
        ra3 = RandXor() % (len * 2 + 1) - len + pred_y[ra1];
      }
      ra2 = max(ra2, lx[ra1]);
      ra2 = min(ra2, rx[ra1]);
      ra3 = max(ra3, ly[ra1]);
      ra3 = min(ra3, ry[ra1]);

      for (auto q_num : num_queries[ra1]) {
        auto zen = BuildMST(queries[q_num]);
        double zenSum = buildMST_sum;
        double querySum = 0;
        for (auto p : queryAnswers[q_num]) {
          querySum += sqrt(Distance(GetPoint(p.first), GetPoint(p.second)));
        }
        tmpScore -= querySum / zenSum;
      }

      keep2 = pred_x[ra1];
      keep3 = pred_y[ra1];
      pred_x[ra1] = ra2;
      pred_y[ra1] = ra3;

      for (auto q_num : num_queries[ra1]) {
        auto zen = BuildMST(queries[q_num]);
        double zenSum = buildMST_sum;
        double querySum = 0;
        for (auto p : queryAnswers[q_num]) {
          querySum += sqrt(Distance(GetPoint(p.first), GetPoint(p.second)));
        }
        tmpScore += querySum / zenSum;
      }
    }
    else if (raMode < hypers.Partition[1]) {

    }


    // 焼きなまし
    double diffScore = (dScore - tmpScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > Rand01()) {
      // 採用
      dScore = tmpScore;

      // Best解よりもいいか
      if (dScore < best_dScore) {
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {
        pred_x[ra1] = keep2;
        pred_y[ra1] = keep3;
      }
      else if (raMode < hypers.Partition[1]) {
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << "sa0 : " << loop << endl;
  }

  CopyToAns();

  rep(i, m)
  {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }
  ansScore = CalcScoreAll();

  CopyToBest();

  //rep(i, 100)
  //{
  //  cout << pred_x[i] << ' ' << pred_y[i] << endl;
  //}
}

void SimulatedAnnealing1(Hypers hypers)
{
  double timeLimit = TL * 3 / 4;
  double startTime = GetNowTime();

  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    ll tmpScore = ansScore;

    // 近傍解作成
    int raMode = RandXor() % hypers.Partition[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      if (m == 1)continue;
      ra1 = RandXor() % n;
      ra2 = RandXor() % n;
      while (ans[ra1] == ans[ra2]) {
        ra1 = RandXor() % n;
        ra2 = RandXor() % n;
      }

      int g_num1 = ans[ra1];
      int g_num2 = ans[ra2];

      tmpScore -= BoundingScore(g_num1);
      tmpScore -= BoundingScore(g_num2);

      ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), ra1));
      ans_nums[g_num1].push_back(ra2);
      ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), ra2));
      ans_nums[g_num2].push_back(ra1);

      tmpScore += BoundingScore(g_num1);
      tmpScore += BoundingScore(g_num2);

      swap(ans[ra1], ans[ra2]);
    }
    else if (raMode < hypers.Partition[1]) {

    }

    // 焼きなまし
    double diffScore = (ansScore - tmpScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > Rand01()) {
      // 採用
      ansScore = tmpScore;

      // Best解よりもいいか
      if (ansScore < best_ansScore) {
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {
        int g_num1 = ans[ra1];
        int g_num2 = ans[ra2];

        ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), ra1));
        ans_nums[g_num1].push_back(ra2);
        ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), ra2));
        ans_nums[g_num2].push_back(ra1);

        swap(ans[ra1], ans[ra2]);
      }
      else if (raMode < hypers.Partition[1]) {
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << "sa1 : " << loop << endl;
  }

  CopyToAns();

  rep(i, m)
  {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }
  ansScore = CalcScoreAll();

  CopyToBest();
}

void SimulatedAnnealing2(Hypers hypers)
{
  double timeLimit = TL;
  double startTime = GetNowTime();

  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0 || sort_g[m - 1] >= 100) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    double tmpScore = INT_INF;

    // 近傍解作成
    int raMode = RandXor() % hypers.Partition[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      if (m == 1)continue;
      ra1 = RandXor() % n;
      ra2 = RandXor() % n;
      while (ans[ra1] == ans[ra2]) {
        ra1 = RandXor() % n;
        ra2 = RandXor() % n;
      }

      int g_num1 = ans[ra1];
      int g_num2 = ans[ra2];

      tmpScore = ansScore;
      tmpScore -= CalcScore(g_num1);
      tmpScore -= CalcScore(g_num2);

      ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), ra1));
      ans_nums[g_num1].push_back(ra2);
      ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), ra2));
      ans_nums[g_num2].push_back(ra1);

      ans_edges[g_num1] = BuildMST(ans_nums[g_num1]);
      ans_edges[g_num2] = BuildMST(ans_nums[g_num2]);

      tmpScore += CalcScore(g_num1);
      tmpScore += CalcScore(g_num2);

      swap(ans[ra1], ans[ra2]);
    }
    else if (raMode < hypers.Partition[1]) {

    }

    // 焼きなまし
    double diffScore = (ansScore - tmpScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > Rand01()) {
      // 採用
      ansScore = tmpScore;

      // Best解よりもいいか
      if (ansScore < best_ansScore) {
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {
        int g_num1 = ans[ra1];
        int g_num2 = ans[ra2];

        ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), ra1));
        ans_nums[g_num1].push_back(ra2);
        ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), ra2));
        ans_nums[g_num2].push_back(ra1);

        ans_edges[g_num1] = BuildMST(ans_nums[g_num1]);
        ans_edges[g_num2] = BuildMST(ans_nums[g_num2]);

        swap(ans[ra1], ans[ra2]);
      }
      else if (raMode < hypers.Partition[1]) {
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << "sa2 : " << loop << endl;
  }

  CopyToAns();
  ansScore = CalcScoreAll();
  CopyToBest();
}

int sa3_graph[n][n];
int sa3_graph_count[n];

int sa3_subtree_size[n];
bool sa3_visited[n];
int sa3_parent[n];
int sa3_dfs(int u, int p)
{
  sa3_parent[u] = p;
  sa3_visited[u] = true;
  int cnt = 1; // 自分自身を含めるので1からスタート
  rep(i, sa3_graph_count[u])
  {
    if (!sa3_visited[sa3_graph[u][i]]) {
      cnt += sa3_dfs(sa3_graph[u][i], u);
    }
  }
  sa3_subtree_size[u] = cnt;
  return cnt;
}

void sa3_dfs2(int u, int p, int g_num)
{
  ans[u] = g_num;
  ans_nums[g_num].push_back(u);
  rep(i, sa3_graph_count[u])
  {
    if (sa3_graph[u][i] != p) {
      ans_edges[g_num].emplace_back(u, sa3_graph[u][i]);
      sa3_dfs2(sa3_graph[u][i], u, g_num);
    }
  }
}

void SimulatedAnnealing3(Hypers hypers)
{
  double timeLimit = TL;
  double startTime = GetNowTime();

  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0 || sort_g[m - 1] >= 100) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    ll tmpScore = ansScore;

    // 近傍解作成
    int raMode = RandXor() % hypers.Partition[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;
    if (raMode < hypers.Partition[0]) {
      if (m == 1)continue;
      ra1 = RandXor() % m;
      ra2 = RandXor() % m;
      while (ra1 == ra2) {
        ra1 = RandXor() % m;
        ra2 = RandXor() % m;
      }

      int sz1 = sort_g[ra1];
      int sz2 = sort_g[ra2];

      auto vec = ans_nums[ra1];
      vec.insert(vec.end(), ans_nums[ra2].begin(), ans_nums[ra2].end());
      //auto zen = BuildMST(vec);
      auto zen = BuildMST(ra1, ra2);

      for (auto num : vec) {
        sa3_graph_count[num] = 0;
        sa3_visited[num] = false;
      }
      for (auto p : zen) {
        sa3_graph[p.first][sa3_graph_count[p.first]] = p.second;
        sa3_graph_count[p.first]++;
        sa3_graph[p.second][sa3_graph_count[p.second]] = p.first;
        sa3_graph_count[p.second]++;
      }

      sa3_dfs(vec[0], -1);

      int cutNum = -1;
      int ra1Root = -1;
      int ra2Root = -1;
      double cutLen = -1;
      rep(i, zen.size())
      {
        int u = zen[i].first;
        int v = zen[i].second;
        if (sa3_parent[v] == u) {
          if (sa3_subtree_size[v] == sz1 || sa3_subtree_size[v] == sz2) {
            if (Distance(GetPoint(u), GetPoint(v)) > cutLen) {
              cutLen = Distance(GetPoint(u), GetPoint(v));
              cutNum = i;
              if (sa3_subtree_size[v] == sz1) {
                ra1Root = v;
                ra2Root = u;
              }
              else {
                ra1Root = u;
                ra2Root = v;
              }
            }
          }
        }
        else {
          if (sa3_subtree_size[u] == sz1 || sa3_subtree_size[u] == sz2) {
            if (Distance(GetPoint(u), GetPoint(v)) > cutLen) {
              cutLen = Distance(GetPoint(u), GetPoint(v));
              cutNum = i;
              if (sa3_subtree_size[u] == sz1) {
                ra1Root = u;
                ra2Root = v;
              }
              else {
                ra1Root = v;
                ra2Root = u;
              }
            }
          }
        }
      }

      if (cutNum == -1) continue;

      ans_nums[ra1].clear();
      ans_edges[ra1].clear();
      ans_nums[ra2].clear();
      ans_edges[ra2].clear();
      sa3_dfs2(ra1Root, ra2Root, ra1);
      sa3_dfs2(ra2Root, ra1Root, ra2);
    }
    else if (raMode < hypers.Partition[1]) {

    }

    // 焼きなまし
    double diffScore = (ansScore - tmpScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    //if (prob > Rand01()) {
    if (true) {
      // 採用
      //ansScore = tmpScore;

      // Best解よりもいいか
      //if (ansScore < best_ansScore) {
      //  CopyToBest();
      //}
    }
    else {
      // 元に戻す
      if (raMode < hypers.Partition[0]) {

      }
      else if (raMode < hypers.Partition[1]) {
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << "sa3 : " << loop << endl;
  }

  //CopyToAns();

  rep(i, m)
  {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }
  ansScore = CalcScoreAll();

  CopyToBest();
}

// 問題を解く関数
ll Solve(int problem_num, Hypers hypers)
{
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  //if (l >= 5)return 0;

  // 出力ファイルストリームオープン
  OpenOfs(problem_num, ofs);

  // 初期解生成
  Method1();

  // 焼きなまし
  int aespa = 2;
  if (l <= 4)aespa = 4;
  rep(i, aespa)
  {
    Method1_Query(0, q / aespa * (i + 1));
    SimulatedAnnealing0(hypers, TL / 2 / aespa * (i + 1));
  }
  SimulatedAnnealing1(hypers);
  //SimulatedAnnealing2(hypers);
  SimulatedAnnealing3(hypers);

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    cout << CalcScoreAll() << endl;
    score = CalcScoreLocal();

    //rep(i, 100)
    //{
    //  cout << "(" << true_x[i] << ", " << true_y[i] << ") (" << GetMIddlePoint(i).first << ", " << GetMIddlePoint(i).second << ") (" << pred_x[i] << ", " << pred_y[i] << ")" << endl;
    //}
  }
  return score;
}

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;
  isSimulateTruePoint = false;

  Hypers HYPERS;
  HYPERS.StartTemp = 248.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 12345.0;
  HYPERS.Partition[0] = 100;
  HYPERS.Partition[1] = 200;
  HYPERS.Partition[2] = 300;
  HYPERS.Partition[3] = 400;
  HYPERS.Partition[4] = 500;
  HYPERS.Partition[5] = 600;
  HYPERS.Partition[6] = 700;
  HYPERS.Partition[7] = 800;
  HYPERS.Partition[8] = 900;
  HYPERS.Partition[9] = 1000;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 150)
    {
      ll score = Solve(i, HYPERS);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
  }
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition[0] = RandXor() % 101;

      ll sum = 0;
      srep(i, 0, 15)
      {
        ll score = Solve(i, hypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition[0]
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
