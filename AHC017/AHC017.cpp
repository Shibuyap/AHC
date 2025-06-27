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
const int INF = 1000000000;
const ll LLINF = 1001001001001001001;

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

struct edge
{
  int to;
  int cost;
  int id;
};

namespace
{
  int N, M, D, K;
  int u[3100], v[3100], w[3100];
  int X[1100], Y[1100];
  int dist_rank[1100][1100];
  int ans[3100];
  int day_count[32];
  vector<edge> G[1100];
  ll idea;
  ll min_score;
  ll min_score_days[32];
}  // namespace

int dist[1100];
ll dijkstra(int start, int day)
{
  for (int i = 0; i < N; ++i) { dist[i] = INF; }
  dist[start] = 0;
  priority_queue<P, vector<P>, greater<P>> pque;
  pque.push(P(0, start));
  while (!pque.empty()) {
    P p = pque.top();
    pque.pop();
    if (p.first != dist[p.second]) { continue; }
    int x = p.second;
    for (int i = 0; i < G[x].size(); ++i) {
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
  for (int i = 0; i < N; ++i) sum += dist[i];
  return sum;
}

int FINISH_COUNT = 50;
double FIRST_HALF = 0.0;
double current_time;

int NG;
ll dijkstra_limited(int start, int day)
{
  for (int i = 0; i < N; ++i) { dist[i] = INF; }
  dist[start] = 0;
  priority_queue<P, vector<P>, greater<P>> pque;
  pque.push(P(0, start));

  ll sum = 0;
  int cnt = 0;
  while (!pque.empty()) {
    P p = pque.top();
    pque.pop();
    if (p.first != dist[p.second]) { continue; }
    int x = p.second;
    if (dist_rank[start][x] < FINISH_COUNT) {
      sum += dist[x];
      cnt++;
    }
    if (cnt == FINISH_COUNT) {
      break;
    }
    for (int i = 0; i < G[x].size(); ++i) {
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
ll dijkstra_bfs(int start, int day)
{
  bitset<1100> visited(0);
  dist[start] = 0;
  int idx = 0;
  int sz = 0;
  queArr[sz] = start;
  sz++;
  visited[start] = true;

  ll sum = 0;
  int cnt = 0;
  while (idx < sz) {
    int x = queArr[idx];
    idx++;
    if (dist_rank[start][x] < FINISH_COUNT) {
      sum += dist[x];
      cnt++;
    }
    if (cnt == FINISH_COUNT) {
      break;
    }
    for (int i = 0; i < G[x].size(); ++i) {
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


ll calc_score_one_day_mini(int day, int id)
{
  ll sum = dijkstra(u[id], day) + dijkstra(v[id], day);
  return sum;
}

ll calc_score_one_day_mini2(int day, int id)
{
  if (current_time < FIRST_HALF) {
    return dijkstra_bfs(u[id], day) + dijkstra_bfs(v[id], day);
  }
  ll sum = dijkstra_limited(u[id], day) + dijkstra_limited(v[id], day);
  return sum;
}

ll calc_score_one_day_mini_vertex(int day, int id)
{
  if (current_time < FIRST_HALF) {
    return dijkstra_bfs(id, day);
  }
  ll sum = dijkstra_limited(id, day);
  return sum;
}


ll calc_score_one_day(int day)
{
  ll sum = 0;
  for (int i = 0; i < N; ++i) { sum += dijkstra(i, day); }
  return sum;
}

ll calc_score_real()
{
  ll sum = 0;

  for (int day = 0; day < D; ++day) {
    sum += calc_score_one_day(day);
    sum -= idea;
  }

  ll roundSum = round((double)sum * 1000 / D / (N * (N - 1)));
  return roundSum;
}

ll calc_score_real_substitute()
{
  ll sum = min_score;
  sum -= idea * D;
  ll roundSum = round((double)sum * 1000 / D / (N * (N - 1)));
  return roundSum;
}

ll calc_score()
{
  ll sum = 0;
  for (int day = 0; day < D; ++day) { sum += calc_score_one_day(day); }
  return sum;
}

// 初期状態作成（これを呼べばスタート位置に戻れることを想定、real_maxScore等は戻さない）
void init()
{
  // dist_rank
  for (int i = 0; i < N; ++i) {
    vector<P> vec;
    for (int j = 0; j < N; ++j) {
      vec.push_back(
        P((X[j] - X[i]) * (X[j] - X[i]) + (Y[j] - Y[i]) * (Y[j] - Y[i]), j));
    }
    sort(vec.begin(), vec.end());
    for (int j = 0; j < N; ++j) { dist_rank[i][vec[j].second] = j; }
  }
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    cin >> N >> M >> D >> K;
    for (int i = 0; i < M; ++i) {
      cin >> u[i] >> v[i] >> w[i];
      u[i]--;
      v[i]--;
    }
    for (int i = 0; i < N; ++i) { cin >> X[i] >> Y[i]; }
  }
  // ファイル入力する
  else {
    ifs >> N >> M >> D >> K;
    for (int i = 0; i < M; ++i) {
      ifs >> u[i] >> v[i] >> w[i];
      u[i]--;
      v[i]--;
    }
    for (int i = 0; i < N; ++i) { ifs >> X[i] >> Y[i]; }
  }

  init();
}

void input_ans(int problemNum)
{
  string fileNameIfs = "./out/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  // ファイル入力する
  if (ifs.is_open()) {
    for (int i = 0; i < M; ++i) {
      ifs >> ans[i];
      ans[i]--;
    }
  }
}

// 解答出力
void output_data(int mode, int case_num)
{
  if (mode == 0) {
    for (int i = 0; i < M; ++i) { cout << ans[i] + 1 << ' '; }
    cout << endl;
  }

  // ファイル出力
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < M; ++i) { ofs << ans[i] + 1 << ' '; }
    ofs << endl;
    ofs.close();
  }
}

void method1()
{
  for (int i = 0; i < M; ++i) { ans[i] = i / K; }
}

void method2()
{
  for (int i = 0; i < M; ++i) { ans[i] = i % D; }
}

// ある辺をある日にしたらその周辺4辺はその日にしない
void method3()
{
  // シャッフル
  mt19937 mt(rand32());

  for (int i = 0; i < M; ++i) ans[i] = -1;

  vector<int> vec[32];
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < M; ++j) vec[i].push_back(j);
    std::shuffle(vec[i].begin(), vec[i].end(), mt);
  }

  int flag[32][3100];
  for (int i = 0; i < D; ++i) {
    for (int j = 0; j < M; ++j) { flag[i][j] = 0; }
  }

  int cnt[32] = {};
  for (int ii = 0; ii < M; ++ii) {
    for (int j = 0; j < D; ++j) {
      if (cnt[j] == K) { continue; }
      int i = vec[j][ii];
      if (ans[i] != -1 || flag[j][i]) { continue; }
      ans[i] = j;
      cnt[j]++;
      // 周囲4辺をNGに
      int uu = u[i];
      int vv = v[i];
      int uSize = G[uu].size();
      int vSize = G[vv].size();
      for (int k = 0; k < uSize; ++k) {
        if (G[uu][k].id == i) {
          int right = G[uu][(k + 1) % uSize].id;
          int left = G[uu][(k + uSize - 1) % uSize].id;
          flag[j][right] = 1;
          flag[j][left] = 1;
        }
      }
      for (int k = 0; k < vSize; ++k) {
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
void method4()
{
  for (int i = 0; i < M; ++i) ans[i] = -1;
  int now = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < G[i].size(); ++j) {
      int id = G[i][j].id;
      if (ans[id] != -1) { continue; }
      ans[id] = now;
      now = (now + 1) % D;
    }
  }
}

// ある辺の日にちを変える
void inner_method1()
{
  int id = rand32() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = rand32() % D;
    if (day_count[newDay] == K) {
      newDay = day;
    }
  }

  ans[id] = newDay;

  ll newOldDayScore = calc_score_one_day(day);
  ll newNewDayScore = calc_score_one_day(newDay);
  ll diff = (newOldDayScore + newNewDayScore) -
    (min_score_days[day] + min_score_days[newDay]);

  // 小さくなれば更新
  if (diff <= 0) {
    min_score_days[day] = newOldDayScore;
    min_score_days[newDay] = newNewDayScore;
    min_score += diff;
    day_count[day]--;
    day_count[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// 隣接する2つの辺の日にちをスワップする
void inner_method2()
{
  int id = rand32() % M;
  int vertex = rand32() % 2 ? u[id] : v[id];
  int id2 = id;
  while (id2 == id) {
    int sz = G[vertex].size();
    id2 = G[vertex][rand32() % sz].id;
  }
  if (ans[id] == ans[id2]) { return; }

  swap(ans[id], ans[id2]);

  ll Day1Score = calc_score_one_day(ans[id]);
  ll Day2Score = calc_score_one_day(ans[id2]);
  ll diff = (Day1Score + Day2Score) -
    (min_score_days[ans[id]] + min_score_days[ans[id2]]);

  // 小さくなれば更新
  if (diff <= 0) {
    min_score_days[ans[id]] = Day1Score;
    min_score_days[ans[id2]] = Day2Score;
    min_score += diff;
  }
  else {
    swap(ans[id], ans[id2]);
  }
}

// ある辺の日にちを変える
void inner_method3()
{
  int id = rand32() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = rand32() % D;
    if (day_count[newDay] == K) {
      newDay = day;
    }
  }

  ll oldOldDayScore = calc_score_one_day_mini(day, id);
  ll oldNewDayScore = calc_score_one_day_mini(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = calc_score_one_day_mini(day, id);
  ll newNewDayScore = calc_score_one_day_mini(newDay, id);
  ll diff =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  // 小さくなれば更新
  if (diff <= 0) {
    day_count[day]--;
    day_count[newDay]++;
  }
  else {
    ans[id] = day;
  }
}


void inner_method4(double temperature)
{
  int id = rand32() % M;
  int day = ans[id];
  int newDay = day;
  while (newDay == day) {
    newDay = rand32() % D;
    if (day_count[newDay] == K) {
      newDay = day;
    }
  }

  ll oldOldDayScore = calc_score_one_day_mini2(day, id);
  ll oldNewDayScore = calc_score_one_day_mini2(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = calc_score_one_day_mini2(day, id);
  ll newNewDayScore = calc_score_one_day_mini2(newDay, id);
  ll diffScore =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  double prob = exp((double)-diffScore / temperature);
  if (prob > rand_01()) {
    day_count[day]--;
    day_count[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// ある辺の日にちを隣接する辺と同じにする
void inner_method5(double temperature)
{
  int id = rand32() % M;
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
  if (vec.size() == 0) { return; }
  int newDay = vec[rand32() % vec.size()];
  if (day_count[newDay] == K) {
    return;
  }

  ll oldOldDayScore = calc_score_one_day_mini2(day, id);
  ll oldNewDayScore = calc_score_one_day_mini2(newDay, id);

  ans[id] = newDay;

  ll newOldDayScore = calc_score_one_day_mini2(day, id);
  ll newNewDayScore = calc_score_one_day_mini2(newDay, id);
  ll diffScore =
    (newOldDayScore + newNewDayScore) - (oldOldDayScore + oldNewDayScore);

  double prob = exp((double)-diffScore / temperature);
  if (prob > rand_01()) {
    day_count[day]--;
    day_count[newDay]++;
  }
  else {
    ans[id] = day;
  }
}

// 隣接する2つの辺の日にちをスワップする
void inner_method6(double temperature)
{
  int id1 = rand32() % M;
  int vertex = rand32() % 2 ? u[id1] : v[id1];
  int id2 = id1;
  while (id2 == id1) {
    int sz = G[vertex].size();
    id2 = G[vertex][rand32() % sz].id;
  }
  if (ans[id1] == ans[id2]) { return; }

  int day1 = ans[id1];
  int day2 = ans[id2];



  ll oldOldDay1Score = calc_score_one_day_mini2(day1, id1);
  ll oldNewDay1Score = calc_score_one_day_mini2(day2, id1);
  ans[id1] = day2;
  ll newOldDay1Score = calc_score_one_day_mini2(day1, id1);
  ll newNewDay1Score = calc_score_one_day_mini2(day2, id1);

  ll oldOldDay2Score = calc_score_one_day_mini2(day1, id2);
  ll oldNewDay2Score = calc_score_one_day_mini2(day2, id2);
  ans[id2] = day1;
  ll newOldDay2Score = calc_score_one_day_mini2(day1, id2);
  ll newNewDay2Score = calc_score_one_day_mini2(day2, id2);

  ll diffScore = (newOldDay1Score + newNewDay1Score) - (oldOldDay1Score + oldNewDay1Score)
    + (newOldDay2Score + newNewDay2Score) - (oldOldDay2Score + oldNewDay2Score);

  double prob = exp((double)-diffScore / temperature);
  if (prob > rand_01()) {
    ;
  }
  else {
    swap(ans[id1], ans[id2]);
  }
}

// まとめて引っ越し
void inner_method7(double temperature)
{
  int id = rand32() % M;
  int day = ans[id];
  queue<int> que;
  que.push(id);
  set<int> edges;
  edges.insert(id);
  while (que.size()) {
    int e_id = que.front();
    que.pop();
    for (auto e : G[u[e_id]]) {
      if (ans[e.id] != day) { continue; }
      if (edges.find(e.id) == edges.end()) {
        que.push(e.id);
        edges.insert(e.id);
      }
    }
    for (auto e : G[v[e_id]]) {
      if (ans[e.id] != day) { continue; }
      if (edges.find(e.id) == edges.end()) {
        que.push(e.id);
        edges.insert(e.id);
      }
    }
  }

  int newDay = day;
  int sz = edges.size();
  for (int _ = 0; _ < 30; ++_) {
    newDay = rand32() % D;
    if (day_count[newDay] + sz > K) newDay = day;
    if (newDay != day) { break; }
  }
  if (newDay == day) { return; }

  set<int> vertices;
  for (auto e_id : edges) {
    vertices.insert(u[e_id]);
    vertices.insert(v[e_id]);
  }

  ll diffScore = 0;
  for (auto x : vertices) {
    diffScore -= (calc_score_one_day_mini_vertex(day, x) + calc_score_one_day_mini_vertex(newDay, x));
  }
  for (auto e_id : edges) {
    ans[e_id] = newDay;
  }
  for (auto x : vertices) {
    NG = 0;
    diffScore += (calc_score_one_day_mini_vertex(day, x) + calc_score_one_day_mini_vertex(newDay, x));
    if (NG) {
      NG = 0;
      for (auto e_id : edges) {
        ans[e_id] = day;
      }
      return;
    }
  }

  double prob = exp((double)-diffScore / temperature);
  if (prob > rand_01()) {
    day_count[day] -= sz;
    day_count[newDay] += sz;
  }
  else {
    for (auto e_id : edges) {
      ans[e_id] = day;
    }
  }
}

void calc_idea()
{
  idea = 0;
  for (int i = 0; i < N; ++i) { idea += dijkstra(i, -1); }
}

int solve(int mode, int problemNum)
{
  Timer timer;
  timer.start();

  // 入力
  input_data(problemNum);

  // G作成
  {
    for (int i = 0; i < N; ++i) {
      G[i].clear();
    }
    for (int i = 0; i < M; ++i) {
      edge e;
      e.cost = w[i];
      e.id = i;
      e.to = v[i];
      G[u[i]].push_back(e);
      e.to = u[i];
      G[v[i]].push_back(e);
    }
    // 偏角ソートしておく
    for (int i = 0; i < N; ++i) {
      vector<pair<double, int>> GG;
      vector<edge> keepG;
      for (int j = 0; j < G[i].size(); ++j) {
        pair<double, int> p;
        p.first = atan2(Y[G[i][j].to] - Y[i], X[G[i][j].to] - X[i]);
        p.second = j;
        keepG.push_back(G[i][j]);
        GG.push_back(p);
      }
      G[i].clear();
      sort(GG.begin(), GG.end());
      for (int j = 0; j < GG.size(); ++j) { G[i].push_back(keepG[GG[j].second]); }
    }
  }

#if 0
  // 理論値計算
  calc_idea();
#endif

  // Method1();
  // Method2();
  method3();
  // Method4();

#if 0
  // 過去のスコアをインプット
  if (mode != 0) {
    input_ans(problemNum);
  }
#endif

#if 0
  min_score = 0;
  for (int i = 0; i < D; ++i) {
    min_score_days[i] = calc_score_one_day(i);
    min_score += min_score_days[i];
  }
#endif

  for (int i = 0; i < D; ++i) { day_count[i] = 0; }
  for (int i = 0; i < M; ++i) { day_count[ans[i]]++; }

  // 焼きなまし
  current_time = timer.get_elapsed_time();
  double TL = 5.8;
  double nowProgress = current_time / TL;
  int loop = 0;
  double startTemperature = 5000;
  double endTemperature = 0;
  while (true) {
    loop++;
    if (loop % 1 == 0) {
      current_time = timer.get_elapsed_time();
      nowProgress = current_time / TL;
    }
    if (nowProgress > 1.0) {
      break;
    }

#if 0
    if (rand32() % 2 == 0) {
      inner_method1();
    }
    else {
      inner_method2();
    }
#endif
    FINISH_COUNT = nowProgress * 100 + 50;
    double temperature = startTemperature + (endTemperature - startTemperature) * nowProgress;
    if (current_time < FIRST_HALF) {
      temperature = 200 + (0 - 200) * (current_time / FIRST_HALF);
    }
    int ra = rand32() % 100;
    if (ra < 60) {
      inner_method5(temperature);
    }
    else if (ra < 90) {
      inner_method7(temperature);
    }
    else if (ra < 92) {
      inner_method6(temperature);
    }
    else {
      inner_method4(temperature);
    }

#if 0
    if (loop % 10 == 0) {
      output_data(mode, problemNum);
      cout << loop << "  " << calc_score_real_substitute();
      for (int i = 0; i < D; ++i) cout << " " << day_count[i];
      cout << endl;
    }
#endif
  }  // while文ここまで（メインループ）

  // 出力
  output_data(mode, problemNum);

  if (mode != 0) {
    cout << "loop = " << loop << endl;
    calc_idea();
    min_score = 0;
    for (int i = 0; i < D; ++i) {
      min_score_days[i] = calc_score_one_day(i);
      min_score += min_score_days[i];
    }
    cout << "Score = " << calc_score_real_substitute() << endl;
    return calc_score_real_substitute();
  }

  return 0;
}

int main()
{
  int mode = 0;
  if (mode == 0) {
    solve(0, 2);
  }
  else if (mode == 1) {
    solve(1, 2);
  }
  else {
    ll sum = 0;
    for (int _ = 0; _ < 10; ++_) {
      sum += solve(1, _);
    }
    cout << "sum = " << sum << endl;
  }

  return 0;
}
