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
const int INF = 1000000000;
const ll LLINF = 1001001001001001001;

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

struct edge
{
  int to;
  int cost;
  int id;
};

//struct RadixHeap {
//  using uint = unsigned;
//  vector<pair<uint, int>> v[33];
//  uint size, last;
//
//  RadixHeap() : size(0), last(0) {}
//
//  bool empty() const { return size == 0; }
//
//  inline int getbit(int a) { return a ? 32 - __builtin_clz(a) : 0; }
//
//  void push(uint key, const int &value) {
//    ++size;
//    v[getbit(key ^ last)].emplace_back(key, value);
//  }
//
//  pair<uint, int> pop() {
//    if (v[0].empty()) {
//      int idx = 1;
//      while (v[idx].empty()) ++idx;
//      last = min_element(begin(v[idx]), end(v[idx]))->first;
//      for (auto &p : v[idx]) v[getbit(p.first ^ last)].emplace_back(p);
//      v[idx].clear();
//    }
//    --size;
//    auto ret = v[0].back();
//    v[0].pop_back();
//    return ret;
//  }
//};

namespace
{
  int N, M, D, K;
  int u[3100], v[3100], w[3100];
  int X[1100], Y[1100];
  int distRank[1100][1100];
  int ans[3100];
  int dayCount[32];
  vector<edge> G[1100];
  ll idea;
  ll minScore;
  ll minScoreDays[32];
}  // namespace

int dist[1100];
ll Dijkstra(int start, int day)
{
  rep(i, N) { dist[i] = INF; }
  dist[start] = 0;
  priority_queue<P, vector<P>, greater<P>> pque;
  pque.push(P(0, start));
  while (!pque.empty()) {
    P p = pque.top();
    pque.pop();
    if (p.first != dist[p.second]) continue;
    int x = p.second;
    rep(i, G[x].size())
    {
      int y = G[x][i].to;
      if (ans[G[x][i].id] == day) {
        continue;
      }
      int cost = G[x][i].cost;
      if (dist[y] > dist[x] + cost) {
        dist[y] = dist[x] + cost;
        pque.push(P(dist[y], y));
      }
    }
  }
  ll sum = 0;
  rep(i, N) sum += dist[i];
  return sum;
}

int FINISH_COUNT = 50;
double FIRST_HALF = 0.0;
double nowTime;
//ll Dijkstra2(int start, int day) {
//  rep(i, N) { dist[i] = INF; }
//  dist[start] = 0;
//  // priority_queue<P, vector<P>, greater<P>> pque;
//  RadixHeap heap;
//  heap.push(0, start);
//  // pque.push(P(0, start));
//
//  ll sum = 0;
//  int cnt = 0;
//  // while (!pque.empty()) {
//  while (!heap.empty()) {
//    // P p = pque.top();
//    P p = heap.pop();
//    // pque.pop();
//    if (p.first != dist[p.second]) continue;
//    int x = p.second;
//    if (distRank[start][x] < FINISH_COUNT) {
//      sum += dist[x];
//      cnt++;
//    }
//    if (cnt == FINISH_COUNT) {
//      break;
//    }
//    rep(i, G[x].size()) {
//      int y = G[x][i].to;
//      if (ans[G[x][i].id] == day) {
//        continue;
//      }
//      int cost = G[x][i].cost;
//      if (dist[y] > dist[x] + cost) {
//        dist[y] = dist[x] + cost;
//        // pque.push(P(dist[y], y));
//        heap.push(dist[y], y);
//      }
//    }
//  }
//
//  if (cnt < min(FINISH_COUNT, N)) {
//    sum += (ll)(FINISH_COUNT - cnt) * INF;
//  }
//
//  return sum;
//}

int NG;
ll Dijkstra22(int start, int day)
{
  rep(i, N) { dist[i] = INF; }
  dist[start] = 0;
  priority_queue<P, vector<P>, greater<P>> pque;
  pque.push(P(0, start));

  ll sum = 0;
  int cnt = 0;
  while (!pque.empty()) {
    P p = pque.top();
    pque.pop();
    if (p.first != dist[p.second]) continue;
    int x = p.second;
    if (distRank[start][x] < FINISH_COUNT) {
      sum += dist[x];
      cnt++;
    }
    if (cnt == FINISH_COUNT) {
      break;
    }
    rep(i, G[x].size())
    {
      int y = G[x][i].to;
      if (ans[G[x][i].id] == day) {
        continue;
      }
      int cost = G[x][i].cost;
      if (dist[y] > dist[x] + cost) {
        dist[y] = dist[x] + cost;
        pque.push(P(dist[y], y));
      }
    }
  }

  if (cnt < min(FINISH_COUNT, N)) {
    NG = 1;
    sum += (ll)(FINISH_COUNT - cnt) * INF;
  }

  return sum;
}

int queArr[10000];
ll Dijkstra01(int start, int day)
{
  bitset<1100> visited(0);
  dist[start] = 0;
  int ite = 0;
  int sz = 0;
  queArr[sz] = start;
  sz++;
  visited[start] = true;

  ll sum = 0;
  int cnt = 0;
  while (ite < sz) {
    int x = queArr[ite];
    ite++;
    if (distRank[start][x] < FINISH_COUNT) {
      sum += dist[x];
      cnt++;
    }
    if (cnt == FINISH_COUNT) {
      break;
    }
    rep(i, G[x].size())
    {
      int y = G[x][i].to;
      if (ans[G[x][i].id] == day) {
        continue;
      }
      if (!visited[y]) {
        visited[y] = true;
        dist[y] = dist[x] + 1;
        queArr[sz] = y;
        sz++;
      }
    }
  }

  if (cnt < min(FINISH_COUNT, N)) {
    NG = 1;
    sum += (ll)(FINISH_COUNT - cnt) * INF;
  }

  return sum;
}


ll CalcScoreOneDayMini(int day, int id)
{
  ll sum = Dijkstra(u[id], day) + Dijkstra(v[id], day);
  return sum;
}

ll CalcScoreOneDayMini2(int day, int id)
{
  if (nowTime < FIRST_HALF) {
    return Dijkstra01(u[id], day) + Dijkstra01(v[id], day);
  }
  ll sum = Dijkstra22(u[id], day) + Dijkstra22(v[id], day);
  return sum;
}

ll CalcScoreOneDayMiniVertex(int day, int id)
{
  if (nowTime < FIRST_HALF) {
    return Dijkstra01(id, day);
  }
  ll sum = Dijkstra22(id, day);
  return sum;
}


ll CalcScoreOneDay(int day)
{
  ll sum = 0;
  rep(i, N) { sum += Dijkstra(i, day); }
  return sum;
}

ll CalcScoreReal()
{
  ll sum = 0;

  rep(day, D)
  {
    sum += CalcScoreOneDay(day);
    sum -= idea;
  }

  ll roundSum = round((double)sum * 1000 / D / (N * (N - 1)));
  return roundSum;
}

ll CalcScoreRealSubstitute()
{
  ll sum = minScore;
  sum -= idea * D;
  ll roundSum = round((double)sum * 1000 / D / (N * (N - 1)));
  return roundSum;
}

ll CalcScore()
{
  ll sum = 0;
  rep(day, D) { sum += CalcScoreOneDay(day); }
  return sum;
}

// 初期状態作成（これを呼べばスタート位置に戻れることを想定、real_maxScore等は戻さない）
void Init()
{
  // distRank
  rep(i, N)
  {
    vector<P> vec;
    rep(j, N)
    {
      vec.push_back(
        P((X[j] - X[i]) * (X[j] - X[i]) + (Y[j] - Y[i]) * (Y[j] - Y[i]), j));
    }
    sort(vec.begin(), vec.end());
    rep(j, N) { distRank[i][vec[j].second] = j; }
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
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
    cin >> N >> M >> D >> K;
    rep(i, M)
    {
      cin >> u[i] >> v[i] >> w[i];
      u[i]--;
      v[i]--;
    }
    rep(i, N) { cin >> X[i] >> Y[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N >> M >> D >> K;
    rep(i, M)
    {
      ifs >> u[i] >> v[i] >> w[i];
      u[i]--;
      v[i]--;
    }
    rep(i, N) { ifs >> X[i] >> Y[i]; }
  }

  Init();
}

void InputAns(int problemNum)
{
  string fileNameIfs = "./out/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // ファイル入力する
  if (ifs.is_open()) {
    rep(i, M)
    {
      ifs >> ans[i];
      ans[i]--;
    }
  }
}

// 解答出力
void Output(int mode, int problemNum)
{
  if (mode == 0) {
    rep(i, M) { cout << ans[i] + 1 << ' '; }
    cout << endl;
  }

  // ファイル出力
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

    ofstream ofs(fileNameOfs);

    rep(i, M) { ofs << ans[i] + 1 << ' '; }
    ofs << endl;
    ofs.close();
  }
}

void Method1()
{
  rep(i, M) { ans[i] = i / K; }
}

void Method2()
{
  rep(i, M) { ans[i] = i % D; }
}

// ある辺をある日にしたらその周辺4辺はその日にしない
void Method3()
{
  // シャッフル
  mt19937 mt(Rand());

  rep(i, M) ans[i] = -1;

  vector<int> vec[32];
  rep(i, D)
  {
    rep(j, M) vec[i].push_back(j);
    std::shuffle(vec[i].begin(), vec[i].end(), mt);
  }

  int flag[32][3100];
  rep(i, D)
  {
    rep(j, M) { flag[i][j] = 0; }
  }

  int cnt[32] = {};
  rep(ii, M)
  {
    rep(j, D)
    {
      if (cnt[j] == K) continue;
      int i = vec[j][ii];
      if (ans[i] != -1 || flag[j][i]) continue;
      ans[i] = j;
      cnt[j]++;
      // 周囲4辺をNGに
      int uu = u[i];
      int vv = v[i];
      int uSize = G[uu].size();
      int vSize = G[vv].size();
      rep(k, uSize)
      {
        if (G[uu][k].id == i) {
          int right = G[uu][(k + 1) % uSize].id;
          int left = G[uu][(k + uSize - 1) % uSize].id;
          flag[j][right] = 1;
          flag[j][left] = 1;
        }
      }
      rep(k, vSize)
      {
        if (G[vv][k].id == i) {
          int right = G[vv][(k + 1) % vSize].id;
          int left = G[vv][(k + vSize - 1) % vSize].id;
          flag[j][right] = 1;
          flag[j][left] = 1;
        }
      }
    }
  }
}

// 各頂点でばらけさせる
void Method4()
{
  rep(i, M) ans[i] = -1;
  int now = 0;
  rep(i, N)
  {
    rep(j, G[i].size())
    {
      int id = G[i][j].id;
      if (ans[id] != -1) continue;
      ans[id] = now;
      now = (now + 1) % D;
    }
  }
}

// ある辺の日にちを変える
void InnerMethod1()
{
  int id = Rand() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = Rand() % D;
    if (dayCount[newDay] == K) {
      newDay = day;
    }
  }

  ans[id] = newDay;

  ll newOldDayScore = CalcScoreOneDay(day);
  ll newNewDayScore = CalcScoreOneDay(newDay);
  ll diff = (newOldDayScore + newNewDayScore) -
    (minScoreDays[day] + minScoreDays[newDay]);

  // 小さくなれば更新
  if (diff <= 0) {
    minScoreDays[day] = newOldDayScore;
    minScoreDays[newDay] = newNewDayScore;
    minScore += diff;
    dayCount[day]--;
    dayCount[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// 隣接する2つの辺の日にちをスワップする
void InnerMethod2()
{
  int id = Rand() % M;
  int vertex = Rand() % 2 ? u[id] : v[id];
  int id2 = id;
  while (id2 == id) {
    int sz = G[vertex].size();
    id2 = G[vertex][Rand() % sz].id;
  }
  if (ans[id] == ans[id2]) return;

  swap(ans[id], ans[id2]);

  ll Day1Score = CalcScoreOneDay(ans[id]);
  ll Day2Score = CalcScoreOneDay(ans[id2]);
  ll diff = (Day1Score + Day2Score) -
    (minScoreDays[ans[id]] + minScoreDays[ans[id2]]);

  // 小さくなれば更新
  if (diff <= 0) {
    minScoreDays[ans[id]] = Day1Score;
    minScoreDays[ans[id2]] = Day2Score;
    minScore += diff;
  }
  else {
    swap(ans[id], ans[id2]);
  }
}

// ある辺の日にちを変える
void InnerMethod3()
{
  int id = Rand() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = Rand() % D;
    if (dayCount[newDay] == K) {
      newDay = day;
    }
  }

  ll oldOldDayScore = CalcScoreOneDayMini(day, id);
  ll oldNewDayScore = CalcScoreOneDayMini(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = CalcScoreOneDayMini(day, id);
  ll newNewDayScore = CalcScoreOneDayMini(newDay, id);
  ll diff =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  // 小さくなれば更新
  if (diff <= 0) {
    dayCount[day]--;
    dayCount[newDay]++;
  }
  else {
    ans[id] = day;
  }
}


void InnerMethod4(double temperature)
{
  int id = Rand() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = Rand() % D;
    if (dayCount[newDay] == K) {
      newDay = day;
    }
  }

  ll oldOldDayScore = CalcScoreOneDayMini2(day, id);
  ll oldNewDayScore = CalcScoreOneDayMini2(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = CalcScoreOneDayMini2(day, id);
  ll newNewDayScore = CalcScoreOneDayMini2(newDay, id);
  ll diffScore =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  double prob = exp((double)-diffScore / temperature);
  if (prob > Rand01()) {
    dayCount[day]--;
    dayCount[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// ある辺の日にちを隣接する辺と同じにする
void InnerMethod5(double temperature)
{
  int id = Rand() % M;
  int day = ans[id];
  set<int> se;
  for (auto e : G[u[id]]) {
    se.insert(ans[e.id]);
  }
  for (auto e : G[v[id]]) {
    se.insert(ans[e.id]);
  }
  vector<int> vec;
  for (auto x : se) {
    if (x != day) vec.push_back(x);
  }
  if (vec.size() == 0) return;
  int newDay = vec[Rand() % vec.size()];
  if (dayCount[newDay] == K) {
    return;
  }

  ll oldOldDayScore = CalcScoreOneDayMini2(day, id);
  ll oldNewDayScore = CalcScoreOneDayMini2(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = CalcScoreOneDayMini2(day, id);
  ll newNewDayScore = CalcScoreOneDayMini2(newDay, id);
  ll diffScore =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  double prob = exp((double)-diffScore / temperature);
  if (prob > Rand01()) {
    dayCount[day]--;
    dayCount[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// 隣接する2つの辺の日にちをスワップする
void InnerMethod6(double temperature)
{
  int id1 = Rand() % M;
  int vertex = Rand() % 2 ? u[id1] : v[id1];
  int id2 = id1;
  while (id2 == id1) {
    int sz = G[vertex].size();
    id2 = G[vertex][Rand() % sz].id;
  }
  if (ans[id1] == ans[id2]) return;

  int day1 = ans[id1];
  int day2 = ans[id2];



  ll oldOldDay1Score = CalcScoreOneDayMini2(day1, id1);
  ll oldNewDay1Score = CalcScoreOneDayMini2(day2, id1);
  ans[id1] = day2;
  ll newOldDay1Score = CalcScoreOneDayMini2(day1, id1);
  ll newNewDay1Score = CalcScoreOneDayMini2(day2, id1);

  ll oldOldDay2Score = CalcScoreOneDayMini2(day1, id2);
  ll oldNewDay2Score = CalcScoreOneDayMini2(day2, id2);
  ans[id2] = day1;
  ll newOldDay2Score = CalcScoreOneDayMini2(day1, id2);
  ll newNewDay2Score = CalcScoreOneDayMini2(day2, id2);

  ll diffScore = (newOldDay1Score + newNewDay1Score) - (oldOldDay1Score + oldNewDay1Score)
    + (newOldDay2Score + newNewDay2Score) - (oldOldDay2Score + oldNewDay2Score);

  double prob = exp((double)-diffScore / temperature);
  if (prob > Rand01()) {
    ;
  }
  else {
    swap(ans[id1], ans[id2]);
  }
}

// まとめて引っ越し
void InnerMethod7(double temperature)
{
  int id = Rand() % M;
  int day = ans[id];
  queue<int> que;
  que.push(id);
  set<int> edges;
  edges.insert(id);
  while (que.size()) {
    int e_id = que.front();
    que.pop();
    for (auto e : G[u[e_id]]) {
      if (ans[e.id] != day) continue;
      if (edges.find(e.id) == edges.end()) {
        que.push(e.id);
        edges.insert(e.id);
      }
    }
    for (auto e : G[v[e_id]]) {
      if (ans[e.id] != day) continue;
      if (edges.find(e.id) == edges.end()) {
        que.push(e.id);
        edges.insert(e.id);
      }
    }
  }

  int newDay = day;
  int sz = edges.size();
  rep(_, 30)
  {
    newDay = Rand() % D;
    if (dayCount[newDay] + sz > K) newDay = day;
    if (newDay != day) break;
  }
  if (newDay == day) return;

  set<int> vertices;
  for (auto e_id : edges) {
    vertices.insert(u[e_id]);
    vertices.insert(v[e_id]);
  }

  ll diffScore = 0;
  for (auto x : vertices) {
    diffScore -= (CalcScoreOneDayMiniVertex(day, x) + CalcScoreOneDayMiniVertex(newDay, x));
  }
  for (auto e_id : edges) {
    ans[e_id] = newDay;
  }
  for (auto x : vertices) {
    NG = 0;
    diffScore += (CalcScoreOneDayMiniVertex(day, x) + CalcScoreOneDayMiniVertex(newDay, x));
    if (NG) {
      NG = 0;
      for (auto e_id : edges) {
        ans[e_id] = day;
      }
      return;
    }
  }

  double prob = exp((double)-diffScore / temperature);
  if (prob > Rand01()) {
    dayCount[day] -= sz;
    dayCount[newDay] += sz;
  }
  else {
    for (auto e_id : edges) {
      ans[e_id] = day;
    }
  }
}

void CalcIdea()
{
  idea = 0;
  rep(i, N) { idea += Dijkstra(i, -1); }
}

int Solve(int mode, int problemNum)
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  // 入力
  Input(problemNum);

  // G作成
  {
    rep(i, N)
    {
      G[i].clear();
    }
    rep(i, M)
    {
      edge e;
      e.cost = w[i];
      e.id = i;
      e.to = v[i];
      G[u[i]].push_back(e);
      e.to = u[i];
      G[v[i]].push_back(e);
    }
    // 偏角ソートしておく
    rep(i, N)
    {
      vector<pair<double, int>> GG;
      vector<edge> keepG;
      rep(j, G[i].size())
      {
        pair<double, int> p;
        p.first = atan2(Y[G[i][j].to] - Y[i], X[G[i][j].to] - X[i]);
        p.second = j;
        keepG.push_back(G[i][j]);
        GG.push_back(p);
      }
      G[i].clear();
      sort(GG.begin(), GG.end());
      rep(j, GG.size()) { G[i].push_back(keepG[GG[j].second]); }
    }
  }

#if 0
  // 理論値計算
  CalcIdea();
#endif

  // Method1();
  // Method2();
  Method3();
  // Method4();

#if 0
  // 過去のスコアをインプット
  if (mode != 0) {
    InputAns(problemNum);
  }
#endif

#if 0
  minScore = 0;
  rep(i, D)
  {
    minScoreDays[i] = CalcScoreOneDay(i);
    minScore += minScoreDays[i];
  }
#endif

  rep(i, D) { dayCount[i] = 0; }
  rep(i, M) { dayCount[ans[i]]++; }

  // 焼きなまし
  endTime = clock();
  nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  double TL = 5.8;
  double nowProgress = nowTime / TL;
  int loop = 0;
  double startTemperature = 5000;
  double endTemperature = 0;
  while (true) {
    loop++;
    if (loop % 1 == 0) {
      endTime = clock();
      nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
      nowProgress = nowTime / TL;
    }
    if (nowProgress > 1.0) break;

#if 0
    if (Rand() % 2 == 0) {
      InnerMethod1();
    }
    else {
      InnerMethod2();
    }
#endif
    FINISH_COUNT = nowProgress * 100 + 50;
    double temperature = startTemperature + (endTemperature - startTemperature) * nowProgress;
    if (nowTime < FIRST_HALF) {
      temperature = 200 + (0 - 200) * (nowTime / FIRST_HALF);
    }
    int ra = Rand() % 100;
    if (ra < 60) {
      InnerMethod5(temperature);
    }
    else if (ra < 90) {
      InnerMethod7(temperature);
    }
    else if (ra < 92) {
      InnerMethod6(temperature);
    }
    else {
      InnerMethod4(temperature);
    }

#if 0
    if (loop % 10 == 0) {
      Output(mode, problemNum);
      cout << loop << "  " << CalcScoreRealSubstitute();
      rep(i, D) cout << " " << dayCount[i];
      cout << endl;
    }
#endif
  }  // while文ここまで（メインループ）

  // 出力
  Output(mode, problemNum);

  if (mode != 0) {
    cout << "loop = " << loop << endl;
    CalcIdea();
    minScore = 0;
    rep(i, D)
    {
      minScoreDays[i] = CalcScoreOneDay(i);
      minScore += minScoreDays[i];
    }
    cout << "Score = " << CalcScoreRealSubstitute() << endl;
    return CalcScoreRealSubstitute();
  }

  return 0;
}

int main()
{
  int mode = 0;
  if (mode == 0) {
    Solve(0, 2);
  }
  else if (mode == 1) {
    Solve(1, 2);
  }
  else {
    ll sum = 0;
    rep(_, 10)
    {
      sum += Solve(1, _);
    }
    cout << "sum = " << sum << endl;
  }

  return 0;
}
