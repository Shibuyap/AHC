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
typedef pair<P, P> PP;

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

static double RandRange(double l, double r)
{
  return l + (r - l) * Rand01();
}

// [l, r]
static uint32_t RandRange(uint32_t l, uint32_t r)
{
  return l + Rand() % (r - l + 1);
}


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタ
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

// 2次元キュー
int queueArr2[10000][2];
int queueHead2 = 0;
int queueTail2 = 0;
void ClearQueue2()
{
  queueHead2 = 0;
  queueTail2 = 0;
}
int Front2X()
{
  return queueArr2[queueHead2][0];
}
int Front2Y()
{
  return queueArr2[queueHead2][1];
}
void Push2(int x, int y)
{
  queueArr2[queueTail2][0] = x;
  queueArr2[queueTail2][1] = y;
  queueTail2++;
}
void Pop2()
{
  queueHead2++;
}
int Size2()
{
  return queueTail2 - queueHead2;
}

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
    if (x == y) { return; }

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


const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;


const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;

std::chrono::steady_clock::time_point startTimeClock;

void ResetTime()
{
  startTimeClock = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTimeClock;
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
vector<P> queryPredMST[q];
int  queryPredMSTSum[q];
vector<int> num_queries[n];
set<vector<int>> querySet;

set<ll> edgeCompareSet;
set<ll> shorterEdges[n * n];
set<ll> longerEdges[n * n];

int true_x[n], true_y[n];

ll ansScore;
double dScore;
int ans[n];
vector<int> ans_nums[MAX_M];
vector<P> ans_edges[MAX_M];
int ans_MSTSums[MAX_M];

ll best_ansScore;
double best_dScore;
int best_ans[n];
vector<int> best_ans_nums[MAX_M];
vector<P> best_ans_edges[MAX_M];
int best_ans_MSTSums[MAX_M];

int free_use_graph[n][n];
int free_use_graph_count[n];
bool free_use_visited[n];
int free_use_parent[n];

void CopyToBest()
{
  best_ansScore = ansScore;
  best_dScore = dScore;
  for (int i = 0; i < n; ++i) {
    best_ans[i] = ans[i];
  }
  for (int i = 0; i < m; ++i) {
    best_ans_nums[i] = ans_nums[i];
    best_ans_edges[i] = ans_edges[i];
    best_ans_MSTSums[i] = ans_MSTSums[i];
  }
}

void CopyToAns()
{
  ansScore = best_ansScore;
  dScore = best_dScore;
  for (int i = 0; i < n; ++i) {
    ans[i] = best_ans[i];
  }
  for (int i = 0; i < m; ++i) {
    ans_nums[i] = best_ans_nums[i];
    ans_edges[i] = best_ans_edges[i];
    ans_MSTSums[i] = best_ans_MSTSums[i];
  }
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = INT_INF;
  for (int i = 0; i < MAX_M; ++i) {
    ans_nums[i].clear();
    ans_edges[i].clear();
  }
  CopyToBest();

  queryCount = 0;
  for (int i = 0; i < q; ++i) {
    queries[i].clear();
    queryAnswers[i].clear();
    queryPredMST[i].clear();
  }
  for (int i = 0; i < n; ++i) {
    num_queries[i].clear();
  }
  querySet.clear();
  edgeCompareSet.clear();
  for (int i = 0; i < n * n; ++i) {
    shorterEdges[i].clear();
    longerEdges[i].clear();
  }
}

constexpr ll pown1 = 800;
constexpr ll pown2 = 640000;
constexpr ll pown3 = 512000000;
ll MakeLL(int a, int b, int c, int d)
{
  return a * pown3 + b * pown2 + c * pown1 + d;
}

void InitShorterEdges()
{
  for (auto num : edgeCompareSet) {
    ll a = num / pown3;
    ll b = num % pown3 / pown2;
    ll ab = a * pown1 + b;
    ll cd = num % pown2;
    shorterEdges[cd].insert(ab);
    longerEdges[ab].insert(cd);
  }
}

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

void BuildMST_InitPoints(const vector<int>& nums, int startIndex, bool isTrue)
{
  for (int i = 0; i < nums.size(); ++i) {
    if (isTrue) {
      buildMST_points[startIndex + i] = GetTruePoint(nums[i]);
    }
    else {
      buildMST_points[startIndex + i] = GetPoint(nums[i]);
    }
    buildMST_idxDict[nums[i]] = startIndex + i;
  }
}

vector<P> BuildMST(const vector<int>& nums, bool isTrue = false)
{
  buildMST_sum = 0;

  UF_Init(nums.size());

  BuildMST_InitPoints(nums, 0, isTrue);

  int edgeCount = 0;
  for (int i = 0; i < nums.size(); ++i) {
    for (int j = i + 1; j < nums.size(); ++j) {
      buildMST_edges[edgeCount].dist = Distance(buildMST_points[i], buildMST_points[j]);
      buildMST_edges[edgeCount].u = i;
      buildMST_edges[edgeCount].v = j;
      edgeCount++;
    }
  }

  sort(buildMST_edges, buildMST_edges + edgeCount);

  vector<P> res(nums.size() - 1);
  int resCount = 0;

  for (int i = 0; i < edgeCount; ++i) {
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

int buildMSTWithOnePoint_used[n];
vector<P> BuildMSTWithOnePoint(const vector<int>& _nums, vector<int>& newNums, bool isTrue = false)
{
  int vCount = _nums.size();
  auto nums = _nums;
  for (int i = 0; i < m; ++i) {
    if (sort_g[i] == 1) {
      nums.push_back(ans_nums[i][0]);
    }
    else {
      break;
    }
  }

  UF_Init(nums.size());

  BuildMST_InitPoints(nums, 0, isTrue);

  int edgeCount = 0;
  for (int i = 0; i < nums.size(); ++i) {
    for (int j = i + 1; j < nums.size(); ++j) {
      buildMST_edges[edgeCount].dist = Distance(buildMST_points[i], buildMST_points[j]);
      buildMST_edges[edgeCount].u = i;
      buildMST_edges[edgeCount].v = j;
      edgeCount++;
    }
  }

  sort(buildMST_edges, buildMST_edges + edgeCount);

  vector<P> res;
  vector<P> idxs;
  int root = -1;
  int rootNum = -1;
  int sz = -1;

  for (int i = 0; i < edgeCount; ++i) {
    if (!UF_Same(buildMST_edges[i].u, buildMST_edges[i].v)) {
      UF_Unite(buildMST_edges[i].u, buildMST_edges[i].v);
      res.push_back(P(nums[buildMST_edges[i].u], nums[buildMST_edges[i].v]));
      idxs.push_back(P(buildMST_edges[i].u, buildMST_edges[i].v));
      if (UF_Count(buildMST_edges[i].u) >= vCount) {
        root = UF_Find(buildMST_edges[i].u);
        rootNum = nums[root];
        sz = UF_Count(buildMST_edges[i].u);
        break;
      }
    }
  }

  if (sz != vCount) {
    return {};
  }

  buildMST_sum = 0;
  newNums.clear();

  vector<P> es;
  for (int i = 0; i < res.size(); ++i) {
    if (UF_Find(idxs[i].first) == root) {
      es.push_back(res[i]);
    }
  }

  for (auto p : es) {
    free_use_graph_count[p.first] = 0;
    free_use_graph_count[p.second] = 0;
  }

  for (auto p : es) {
    free_use_graph[p.first][free_use_graph_count[p.first]] = p.second;;
    free_use_graph[p.second][free_use_graph_count[p.second]] = p.first;
    free_use_graph_count[p.first]++;
    free_use_graph_count[p.second]++;
  }

  priority_queue<pair<int, P>, vector<pair<int, P>>, greater<pair<int, P>>> pque;
  pair<int, P> pp;
  for (int i = 0; i < free_use_graph_count[rootNum]; ++i) {
    int y = free_use_graph[rootNum][i];
    pp.first = Distance(rootNum, y);
    pp.second.first = rootNum;
    pp.second.second = y;
    pque.emplace(pp);
  }
  newNums.push_back(rootNum);

  vector<P> es2;
  while (pque.size()) {
    pp = pque.top();
    pque.pop();
    es2.emplace_back(pp.second);
    buildMST_sum += sqrt(pp.first);
    newNums.push_back(pp.second.second);

    if (es2.size() == vCount - 1) {
      break;
    }

    int par = pp.second.first;
    int x = pp.second.second;
    for (int i = 0; i < free_use_graph_count[x]; ++i) {
      int y = free_use_graph[x][i];
      if (y == par)continue;
      pp.first = Distance(x, y);
      pp.second.first = x;
      pp.second.second = y;
      pque.emplace(pp);
    }
  }

  int now = 0;
  for (auto num : nums) {
    buildMSTWithOnePoint_used[num] = true;
  }
  for (auto num : newNums) {
    buildMSTWithOnePoint_used[num] = false;
  }
  for (auto num : nums) {
    if (buildMSTWithOnePoint_used[num]) {
      ans_nums[now][0] = num;
      ans[num] = now;
      now++;
    }
  }

  return es2;
}


vector<P> BuildMSTWithEdgeCompare(const vector<int>& nums, bool isTrue = false)
{
  buildMST_sum = 0;

  UF_Init(nums.size());

  BuildMST_InitPoints(nums, 0, isTrue);

  int edgeCount = 0;
  for (int i = 0; i < nums.size(); ++i) {
    for (int j = i + 1; j < nums.size(); ++j) {
      buildMST_edges[edgeCount].dist = Distance(buildMST_points[i], buildMST_points[j]);
      buildMST_edges[edgeCount].u = i;
      buildMST_edges[edgeCount].v = j;
      edgeCount++;
    }
  }

  sort(buildMST_edges, buildMST_edges + edgeCount);

  vector<P> res(nums.size() - 1);
  int resCount = 0;

  vector<int> numsa, numsb;
  for (int i = 0; i < edgeCount; ++i) {
    int ia = buildMST_edges[i].u;
    int ib = buildMST_edges[i].v;
    if (!UF_Same(ia, ib)) {
      for (int _ = 0; _ < 5; ++_) {
        int a = nums[ia];
        int b = nums[ib];
        if (a > b) {
          swap(a, b);
          swap(ia, ib);
        }

        numsa.clear();
        numsb.clear();
        numsa.push_back(ia);
        numsb.push_back(ib);
        if (UF_Count(ia) > 1 || UF_Count(ib) > 1) {
          int ra = UF_Find(ia);
          int rb = UF_Find(ib);
          for (int j = 0; j < nums.size(); ++j) {
            if (j == ia || j == ib)continue;
            int rj = UF_Find(j);
            if (rj == ra) {
              numsa.push_back(j);
            }
            else if (rj == rb) {
              numsb.push_back(j);
            }
          }
        }

        int ab = a * pown1 + b;
        int nia = -1;
        int nib = -1;
        for (auto ic : numsa) {
          for (auto id : numsb) {
            if (ic == ia && id == ib)continue;
            int c = nums[ic];
            int d = nums[id];
            if (c > d) {
              swap(c, d);
            }
            int cd = c * pown1 + d;
            if (longerEdges[cd].size() < shorterEdges[ab].size()) {
              if (longerEdges[cd].find(ab) != longerEdges[cd].end()) {
                nia = ic;
                nib = id;
                break;
              }
            }
            else {
              if (shorterEdges[ab].find(cd) != shorterEdges[ab].end()) {
                nia = ic;
                nib = id;
                break;
              }
            }
          }
          if (nia != -1)break;
        }

        if (nia != -1) {
          if (!UF_Same(ia, nia) || !UF_Same(ib, nib)) {
            cerr << "NG" << endl;
          }
          ia = nia;
          ib = nib;
        }
        else {
          break;
        }
      }

      buildMST_sum += sqrt(buildMST_edges[i].dist);
      UF_Unite(ia, ib);
      res[resCount] = P(nums[ia], nums[ib]);
      resCount++;
      if (UF_Count(ia) == nums.size()) {
        break;
      }
    }
  }

  return res;
}

vector<P> BuildMST(const vector<int>& nums1, const vector<P>& edges1, const vector<int>& nums2, const vector<P>& edges2, bool isTrue = false)
{
  int n1 = nums1.size();
  int n2 = nums2.size();

  buildMST_sum = 0;

  UF_Init(n1 + n2);

  BuildMST_InitPoints(nums1, 0, isTrue);
  BuildMST_InitPoints(nums2, n1, isTrue);

  int maxDist = -1;
  int edgeCount = 0;
  for (auto p : edges1) {
    buildMST_edges[edgeCount].dist = Distance(buildMST_points[buildMST_idxDict[p.first]], buildMST_points[buildMST_idxDict[p.second]]);
    maxDist = max(maxDist, buildMST_edges[edgeCount].dist);
    buildMST_edges[edgeCount].u = buildMST_idxDict[p.first];
    buildMST_edges[edgeCount].v = buildMST_idxDict[p.second];
    edgeCount++;
  }
  for (auto p : edges2) {
    buildMST_edges[edgeCount].dist = Distance(buildMST_points[buildMST_idxDict[p.first]], buildMST_points[buildMST_idxDict[p.second]]);
    maxDist = max(maxDist, buildMST_edges[edgeCount].dist);
    buildMST_edges[edgeCount].u = buildMST_idxDict[p.first];
    buildMST_edges[edgeCount].v = buildMST_idxDict[p.second];
    edgeCount++;
  }

  Edge minE;
  minE.dist = INT_INF;
  bool contains = false;
  for (int ii = 0; ii < nums1.size(); ++ii) {
    for (int jj = 0; jj < nums2.size(); ++jj) {
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

  vector<P> res((ll)n1 + n2 - 1);
  int resCount = 0;

  for (int i = 0; i < edgeCount; ++i) {
    if (!UF_Same(buildMST_edges[i].u, buildMST_edges[i].v)) {
      buildMST_sum += sqrt(buildMST_edges[i].dist);
      UF_Unite(buildMST_edges[i].u, buildMST_edges[i].v);
      int uu, vv;
      if (buildMST_edges[i].u < n1) {
        uu = nums1[buildMST_edges[i].u];
      }
      else {
        uu = nums2[(ll)buildMST_edges[i].u - n1];
      }
      if (buildMST_edges[i].v < n1) {
        vv = nums1[buildMST_edges[i].v];
      }
      else {
        vv = nums2[(ll)buildMST_edges[i].v - n1];
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

vector<P> BuildMST(const int g_num1, const int g_num2, bool isTrue = false)
{
  return BuildMST(ans_nums[g_num1], ans_edges[g_num1], ans_nums[g_num2], ans_edges[g_num2]);
}

bool isSimulateTruePoint = false;
// 入力を受け取る関数
void Input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int _n, _q;
    cin >> _n >> m >> _q >> l >> w;
    for (int i = 0; i < m; ++i)cin >> g[i];
    for (int i = 0; i < n; ++i) {
      cin >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
  }
  else {
    // ファイル入力
    int _n, _q;
    ifs >> _n >> m >> _q >> l >> w;
    for (int i = 0; i < m; ++i)ifs >> g[i];
    for (int i = 0; i < n; ++i) {
      ifs >> lx[i] >> rx[i] >> ly[i] >> ry[i];
    }
    for (int i = 0; i < n; ++i)ifs >> true_x[i] >> true_y[i];
  }

  if (isSimulateTruePoint) {
    for (int i = 0; i < n; ++i) {
      lx[i] = true_x[i];
      rx[i] = true_x[i];
      ly[i] = true_y[i];
      ry[i] = true_y[i];
    }
  }

  for (int i = 0; i < n; ++i) {
    pred_x[i] = GetMIddlePoint(i).first;
    pred_y[i] = GetMIddlePoint(i).second;
  }

  vector<P> vp;
  for (int i = 0; i < m; ++i) {
    vp.emplace_back(g[i], i);
  }
  sort(vp.begin(), vp.end());
  for (int i = 0; i < m; ++i) {
    sort_g[i] = vp[i].first;
    arg_g[i] = vp[i].second;
    arg_rev_g[vp[i].second] = i;
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int case_num, ofstream& ofs)
{
  if (ofs.is_open()) {
    ofs.close();
  }

  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore(int g_num, bool isTrue = false)
{
  int res = 0;

  for (auto e : ans_edges[g_num]) {
    if (isTrue) {
      res += sqrt(Distance(GetTruePoint(e.first), GetTruePoint(e.second)));
    }
    else {
      res += sqrt(Distance(GetPoint(e.first), GetPoint(e.second)));
    }
  }

  return res;
}

int CalcScoreAll()
{
  int res = 0;
  for (int i = 0; i < m; ++i) {
    res += CalcScore(i);
  }
  return res;
}

int CalcScoreLocal()
{
  int res = 0;
  for (int i = 0; i < m; ++i) {
    res += CalcScore(i, true);
  }
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    cout << '!' << endl;
    for (int i = 0; i < m; ++i) {
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
    for (int i = 0; i < m; ++i) {
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

vector<int> divideTreeByEdge_graph[n];
void DivideTreeByEdge(const vector<int>& nums, const vector<P>& edges, int rootEdgeNum, vector<int>& outNums1, vector<int>& outNums2)
{
  for (auto num : nums) {
    divideTreeByEdge_graph[num].clear();
  }
  for (auto e : edges) {
    divideTreeByEdge_graph[e.first].push_back(e.second);
    divideTreeByEdge_graph[e.second].push_back(e.first);
  }

  outNums1.clear();
  outNums2.clear();

  outNums1.push_back(edges[rootEdgeNum].first);
  ClearQueue2();
  Push2(edges[rootEdgeNum].first, edges[rootEdgeNum].second);
  while (Size2()) {
    int x = Front2X();
    int par = Front2Y();
    Pop2();
    for (auto y : divideTreeByEdge_graph[x]) {
      if (y == par)continue;
      outNums1.push_back(y);
      Push2(y, x);
    }
  }

  outNums2.push_back(edges[rootEdgeNum].second);
  ClearQueue2();
  Push2(edges[rootEdgeNum].second, edges[rootEdgeNum].first);
  while (Size2()) {
    int x = Front2X();
    int par = Front2Y();
    Pop2();
    for (auto y : divideTreeByEdge_graph[x]) {
      if (y == par)continue;
      outNums2.push_back(y);
      Push2(y, x);
    }
  }
}

void Query()
{
  if (mode == 0) {
    cout << "? " << l;
    for (int j = 0; j < l; ++j) {
      cout << ' ' << queries[queryCount][j];
    }
    cout << endl;
    fflush(stdout);

    for (int j = 0; j < l - 1; ++j) {
      int a, b;
      cin >> a >> b;
      queryAnswers[queryCount].push_back(P(a, b));
    }
  }
  else {
    ofs << "? " << l;
    for (int j = 0; j < l; ++j) {
      ofs << ' ' << queries[queryCount][j];
    }
    ofs << endl;

    queryAnswers[queryCount] = BuildMST(queries[queryCount], true);
  }

  queryPredMST[queryCount] = BuildMST(queries[queryCount], false);
  queryPredMSTSum[queryCount] = buildMST_sum;

  for (int i = 0; i < queryAnswers[queryCount].size(); ++i) {
    P edge = queryAnswers[queryCount][i];
    int a = edge.first;
    int b = edge.second;
    if (a > b)swap(a, b);
    vector<int> nums1, nums2;
    DivideTreeByEdge(queries[queryCount], queryAnswers[queryCount], i, nums1, nums2);
    for (auto c : nums1) {
      for (auto d : nums2) {
        if (c > d)swap(c, d);
        if (c == a && d == b)continue;
        edgeCompareSet.insert(MakeLL(a, b, c, d));
      }
    }
  }

  queryCount++;
}

// 面積大きい順にQ個クエリを投げる
// Lは最大まで使う
// 残りL-1頂点は近い順
void Method1_Query(int start, int queryEnd)
{
  vector<P> vp;
  for (int i = 0; i < n; ++i) {
    vp.emplace_back((rx[i] - lx[i] + 1) * (ry[i] - ly[i] + 1), i);
  }
  sort(vp.begin(), vp.end(), greater<P>());

  for (int i = start; i < n; ++i) {
    int num = vp[i].second;
    auto po = GetPoint(num);
    vector<P> dists;
    for (int j = 0; j < n; ++j) {
      dists.emplace_back(Distance(po, GetPoint(j)), j);
    }
    sort(dists.begin(), dists.end());
    queries[queryCount].clear();
    for (int j = 0; j < l; ++j) {
      queries[queryCount].push_back(dists[j].second);
    }
    sort(queries[queryCount].begin(), queries[queryCount].end());
    if (querySet.find(queries[queryCount]) != querySet.end()) {
      continue;
    }
    for (int j = 0; j < l; ++j) {
      num_queries[dists[j].second].push_back(queryCount);
    }
    querySet.insert(queries[queryCount]);
    Query();
    if (queryCount == queryEnd) {
      if (mode == 2) {
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
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < sort_g[i]; ++j) {
      ans[now] = i;
      ans_nums[i].push_back(now);
      now++;
    }
  }

  for (int i = 0; i < m; ++i) {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }

  ansScore = CalcScoreAll();
  CopyToBest();
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp[10];
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

vector<int> sa0_graph[n];
vector<vector<int>> sa0_DivideTreeNums;
vector<vector<P>> sa0_DivideTreeEdges;
void sa0_DivideTree(const vector<int>& nums, const vector<P>& edges, int root)
{
  sa0_DivideTreeNums.clear();
  sa0_DivideTreeEdges.clear();

  for (auto num : nums) {
    sa0_graph[num].clear();
  }
  for (auto e : edges) {
    sa0_graph[e.first].push_back(e.second);
    sa0_graph[e.second].push_back(e.first);
  }

  vector<int> dTreeNums;
  vector<P> dTreeEdges;

  dTreeNums.push_back(root);
  sa0_DivideTreeNums.push_back(dTreeNums);
  sa0_DivideTreeEdges.push_back(dTreeEdges);

  for (auto num : sa0_graph[root]) {
    dTreeNums.clear();
    dTreeEdges.clear();
    ClearQueue2();

    Push2(num, root);
    dTreeNums.push_back(num);

    while (Size2()) {
      int x = Front2X();
      int par = Front2Y();
      Pop2();
      for (auto y : sa0_graph[x]) {
        if (y == par)continue;
        dTreeEdges.emplace_back(x, y);
        dTreeNums.push_back(y);
        Push2(y, x);
      }
    }

    sa0_DivideTreeNums.push_back(dTreeNums);
    sa0_DivideTreeEdges.push_back(dTreeEdges);
  }
}

vector<P> sa0_tmpMST[q];
int sa0_tmpMSTSum[q];
void SimulatedAnnealing0(Hypers hypers, double timeLimit)
{
  double startTime = GetNowTime();

  dScore = 0;
  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp[0];
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) { break; }
    }

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    double tmpScore = dScore;

    // 近傍解作成
    int mode = Rand() % hypers.Partition[0];
    int idx1, idx2, nx, ny, tmp;
    int prev_x, prev_y, tmp3, tmp4, tmp5;

    if (mode < hypers.Partition[0]) {
      idx1 = Rand() % n;
      while (num_queries[idx1].empty()) {
        idx1 = Rand() % n;
      }
      if (Rand() % 100 < 50) {
        nx = Rand() % (rx[idx1] - lx[idx1] + 1) + lx[idx1];
        ny = Rand() % (ry[idx1] - ly[idx1] + 1) + ly[idx1];
      }
      else {
        int len_x = (rx[idx1] - lx[idx1] + 1) / 5 + 1;
        int len_y = (ry[idx1] - ly[idx1] + 1) / 5 + 1;
        nx = Rand() % (len_x * 2 + 1) - len_x + pred_x[idx1];
        ny = Rand() % (len_y * 2 + 1) - len_y + pred_y[idx1];
      }
      nx = max(nx, lx[idx1]);
      nx = min(nx, rx[idx1]);
      ny = max(ny, ly[idx1]);
      ny = min(ny, ry[idx1]);

      for (int i = 0; i < num_queries[idx1].size(); ++i) {
        int q_num = num_queries[idx1][i];
        sa0_tmpMST[i] = queryPredMST[q_num];
        sa0_tmpMSTSum[i] = queryPredMSTSum[q_num];
        double querySum = 0;
        for (auto p : queryAnswers[q_num]) {
          querySum += sqrt(Distance(GetPoint(p.first), GetPoint(p.second)));
        }
        tmpScore -= querySum / sa0_tmpMSTSum[i];
      }

      prev_x = pred_x[idx1];
      prev_y = pred_y[idx1];
      pred_x[idx1] = nx;
      pred_y[idx1] = ny;

      for (int i = 0; i < num_queries[idx1].size(); ++i) {
        int q_num = num_queries[idx1][i];
        queryPredMST[q_num] = BuildMST(queries[q_num]);
        queryPredMSTSum[q_num] = buildMST_sum;
        double querySum = 0;
        for (auto p : queryAnswers[q_num]) {
          querySum += sqrt(Distance(GetPoint(p.first), GetPoint(p.second)));
        }
        tmpScore += querySum / queryPredMSTSum[q_num];
      }
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
      if (mode < hypers.Partition[0]) {
        pred_x[idx1] = prev_x;
        pred_y[idx1] = prev_y;
        for (int i = 0; i < num_queries[idx1].size(); ++i) {
          int q_num = num_queries[idx1][i];
          queryPredMST[q_num] = sa0_tmpMST[i];
          queryPredMSTSum[q_num] = sa0_tmpMSTSum[i];
        }
      }
    }
  }

  if (mode == 2) {
    cout << "sa0 : " << loop << endl;
  }

  CopyToAns();

  for (int i = 0; i < m; ++i) {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }
  ansScore = CalcScoreAll();

  CopyToBest();
}

void SimulatedAnnealing1(Hypers hypers, double timeLimit)
{
  double startTime = GetNowTime();

  CopyToBest();

  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp[1];
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) { break; }
    }

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    ll tmpScore = ansScore;

    // 近傍解作成
    int mode = Rand() % hypers.Partition[0];
    int idx1, idx2, nx, ny, tmp;
    int prev_x, prev_y, tmp3, tmp4, tmp5;
    if (m == 1)continue;
    idx1 = Rand() % n;
    idx2 = Rand() % n;
    while (ans[idx1] == ans[idx2]) {
      idx1 = Rand() % n;
      idx2 = Rand() % n;
    }

    int g_num1 = ans[idx1];
    int g_num2 = ans[idx2];

    tmpScore -= BoundingScore(g_num1);
    tmpScore -= BoundingScore(g_num2);

    ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), idx1));
    ans_nums[g_num1].push_back(idx2);
    ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), idx2));
    ans_nums[g_num2].push_back(idx1);

    tmpScore += BoundingScore(g_num1);
    tmpScore += BoundingScore(g_num2);

    swap(ans[idx1], ans[idx2]);

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
      int g_num1 = ans[idx1];
      int g_num2 = ans[idx2];

      ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), idx1));
      ans_nums[g_num1].push_back(idx2);
      ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), idx2));
      ans_nums[g_num2].push_back(idx1);

      swap(ans[idx1], ans[idx2]);
    }
  }

  if (mode == 2) {
    cout << "sa1 : " << loop << endl;
  }

  CopyToAns();

  for (int i = 0; i < m; ++i) {
    ans_edges[i] = BuildMST(ans_nums[i]);
    ans_MSTSums[i] = buildMST_sum;
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
  const double START_TEMP = hypers.StartTemp[2];
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0 || sort_g[m - 1] >= 100) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) { break; }
    }

    double progressRatio = (nowTime - startTime) / (timeLimit - startTime);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    double tmpScore = INT_INF;

    // 近傍解作成
    int idx1, idx2, tmp3, tmp4, tmp5;
    int keep1, keep2, keep3, keep4, keep5;

    if (m == 1)continue;
    idx1 = Rand() % n;
    idx2 = Rand() % n;
    while (ans[idx1] == ans[idx2]) {
      idx1 = Rand() % n;
      idx2 = Rand() % n;
    }

    int g_num1 = ans[idx1];
    int g_num2 = ans[idx2];

    tmpScore = ansScore;
    tmpScore -= CalcScore(g_num1);
    tmpScore -= CalcScore(g_num2);

    ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), idx1));
    ans_nums[g_num1].push_back(idx2);
    ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), idx2));
    ans_nums[g_num2].push_back(idx1);

    ans_edges[g_num1] = BuildMST(ans_nums[g_num1]);
    ans_edges[g_num2] = BuildMST(ans_nums[g_num2]);

    tmpScore += CalcScore(g_num1);
    tmpScore += CalcScore(g_num2);

    swap(ans[idx1], ans[idx2]);

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
      int g_num1 = ans[idx1];
      int g_num2 = ans[idx2];

      ans_nums[g_num1].erase(find(ans_nums[g_num1].begin(), ans_nums[g_num1].end(), idx1));
      ans_nums[g_num1].push_back(idx2);
      ans_nums[g_num2].erase(find(ans_nums[g_num2].begin(), ans_nums[g_num2].end(), idx2));
      ans_nums[g_num2].push_back(idx1);

      ans_edges[g_num1] = BuildMST(ans_nums[g_num1]);
      ans_edges[g_num2] = BuildMST(ans_nums[g_num2]);

      swap(ans[idx1], ans[idx2]);
    }
  }

  if (mode == 2) {
    cout << "sa2 : " << loop << endl;
  }

  CopyToAns();
  ansScore = CalcScoreAll();
  CopyToBest();
}

int sa3_subtree_size[n];
int sa3_dfs(int u, int p)
{
  free_use_parent[u] = p;
  free_use_visited[u] = true;
  int cnt = 1; // 自分自身を含めるので1からスタート
  for (int i = 0; i < free_use_graph_count[u]; ++i) {
    if (!free_use_visited[free_use_graph[u][i]]) {
      cnt += sa3_dfs(free_use_graph[u][i], u);
    }
  }
  sa3_subtree_size[u] = cnt;
  return cnt;
}

void sa3_dfs2(int u, int p, int g_num)
{
  ans[u] = g_num;
  ans_nums[g_num].push_back(u);
  for (int i = 0; i < free_use_graph_count[u]; ++i) {
    if (free_use_graph[u][i] != p) {
      ans_edges[g_num].emplace_back(u, free_use_graph[u][i]);
      sa3_dfs2(free_use_graph[u][i], u, g_num);
    }
  }
}

void SimulatedAnnealing3(Hypers hypers, double timeLimit)
{
  double startTime = GetNowTime();

  CopyToBest();

  double nowTime = GetNowTime();

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0 || sort_g[m - 1] >= 100) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) { break; }
    }

    // 近傍解作成
    int mode = Rand() % hypers.Partition[1];
    int g1, g2, tmp3, tmp4, tmp5;
    int keep1, keep2, keep3, keep4, keep5;
    if (true || mode < hypers.Partition[0]) {
      if (m == 1)continue;
      g1 = Rand() % m;
      g2 = Rand() % m;
      while (g1 == g2) {
        g1 = Rand() % m;
        g2 = Rand() % m;
      }

      int sz1 = sort_g[g1];
      int sz2 = sort_g[g2];

      auto vec = ans_nums[g1];
      vec.insert(vec.end(), ans_nums[g2].begin(), ans_nums[g2].end());

      vector<P> mst_edges;
      if (Rand() % 100 < 100) {
        if (Rand() % 100 < 90) {
          mst_edges = BuildMST(g1, g2);
        }
        else {
          mst_edges = BuildMST(vec);
        }
      }
      else {
        mst_edges = BuildMSTWithEdgeCompare(vec);
      }

      for (auto num : vec) {
        free_use_graph_count[num] = 0;
        free_use_visited[num] = false;
      }
      for (auto p : mst_edges) {
        free_use_graph[p.first][free_use_graph_count[p.first]] = p.second;
        free_use_graph_count[p.first]++;
        free_use_graph[p.second][free_use_graph_count[p.second]] = p.first;
        free_use_graph_count[p.second]++;
      }

      sa3_dfs(vec[0], -1);

      int cutNum = -1;
      int root1 = -1;
      int root2 = -1;
      double cutLen = -1;
      for (int i = 0; i < mst_edges.size(); ++i) {
        int u = mst_edges[i].first;
        int v = mst_edges[i].second;
        if (free_use_parent[v] == u) {
          if (sa3_subtree_size[v] == sz1 || sa3_subtree_size[v] == sz2) {
            if (Distance(GetPoint(u), GetPoint(v)) > cutLen) {
              cutLen = Distance(GetPoint(u), GetPoint(v));
              cutNum = i;
              if (sa3_subtree_size[v] == sz1) {
                root1 = v;
                root2 = u;
              }
              else {
                root1 = u;
                root2 = v;
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
                root1 = u;
                root2 = v;
              }
              else {
                root1 = v;
                root2 = u;
              }
            }
          }
        }
      }

      if (cutNum == -1) { continue; }

      ans_nums[g1].clear();
      ans_edges[g1].clear();
      ans_nums[g2].clear();
      ans_edges[g2].clear();
      sa3_dfs2(root1, root2, g1);
      sa3_dfs2(root2, root1, g2);
    }
    else if (mode < hypers.Partition[1]) {
      g1 = Rand() % m;
      while (sort_g[g1] == 1) {
        g1 = Rand() % m;
      }

      auto res = BuildMSTWithOnePoint(ans_nums[g1], ans_nums[g1]);
      if (!res.empty()) {
        ans_edges[g1] = res;
        ans_MSTSums[g1] = buildMST_sum;
        for (auto num : ans_nums[g1]) {
          ans[num] = g1;
        }
      }
    }


  }

  if (mode == 2) {
    cout << "sa3 : " << loop << endl;
  }

  for (int i = 0; i < m; ++i) {
    ans_edges[i] = BuildMST(ans_nums[i]);
  }
  ansScore = CalcScoreAll();

  CopyToBest();
}

void SimulatedAnnealing4(Hypers hypers, double timeLimit)
{
  double startTime = GetNowTime();
  double nowTime = GetNowTime();

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0 || sort_g[m - 1] >= 100) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) { break; }
    }

    // 近傍解作成
    int g_idx = Rand() % m;
    ans_edges[g_idx] = BuildMSTWithEdgeCompare(ans_nums[g_idx]);
    ans_MSTSums[g_idx] = buildMST_sum;
  }

  if (mode == 2) {
    cout << "sa4 : " << loop << endl;
  }

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

  // 出力ファイルストリームオープン
  OpenOfs(problem_num, ofs);

  // 初期解生成
  Method1();

  // 焼きなまし
  int phase_cnt = 2;
  if (l <= 4)phase_cnt = 4;
  for (int i = 0; i < phase_cnt; ++i) {
    Method1_Query(0, q / phase_cnt * (i + 1));
    SimulatedAnnealing0(hypers, (TL * 0.6) / phase_cnt * (i + 1));
  }
  InitShorterEdges();
  SimulatedAnnealing1(hypers, TL * 0.7);
  //SimulatedAnnealing2(hypers);
  SimulatedAnnealing3(hypers, TL * 0.8);
  SimulatedAnnealing4(hypers, TL * 1.0);

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScoreLocal();
  }
  return score;
}

int main()
{
  mode = 2;
  isSimulateTruePoint = true;

  Hypers HYPERS;
  HYPERS.StartTemp[0] = 248.0;
  HYPERS.StartTemp[1] = 20048.0;
  HYPERS.StartTemp[2] = 248.0;
  HYPERS.StartTemp[3] = 248.0;
  HYPERS.StartTemp[4] = 248.0;
  HYPERS.StartTemp[5] = 248.0;
  HYPERS.StartTemp[6] = 248.0;
  HYPERS.StartTemp[7] = 248.0;
  HYPERS.StartTemp[8] = 248.0;
  HYPERS.StartTemp[9] = 248.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 12345.0;
  HYPERS.Partition[0] = 100;
  HYPERS.Partition[1] = 110;
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
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(5, HYPERS);
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
      hypers.StartTemp[0] = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition[0] = Rand() % 101;

      ll sum = 0;
      for (int i = 0; i < 10; ++i) {
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
        << ", StartTemp = " << hypers.StartTemp[0]
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
