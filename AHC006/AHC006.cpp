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
#define MAX_N 200005

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

namespace /* 乱数ライブラリ */
{
  static uint32_t Rand() {
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

struct Point
{
  int x;
  int y;
};

namespace /* 変数 */
{
  // 入力用変数
  int n = 1000;
  int m = 60;
  vector<int> start_x(n), start_y(n), goal_x(n), goal_y(n);
  vector<int> argA;

  // 解答用変数
  ll maxScore;
  ll minTime;
  vector<P> current_path;
  vector<int> argAns(140);
  vector<int> argAns2(140);
  vector<int> use(1100);
  vector<int> pair_(2100);

  // 焼きなまし用変数
  ll real_maxScore;
  ll real_minTime;
  vector<P> best_path;
  vector<int> real_argAns(140);
  vector<int> real_argAns2(140);
  vector<int> real_use(1100);
  vector<int> real_pair_(2100);

  void ResetParam() {
    argA.clear();
    maxScore = 0;
    minTime = 0;
    current_path.clear();
    argAns.resize(140);
    argAns2.resize(140);
    use.resize(1100);
    pair_.resize(2100);
    real_maxScore = 0;
    real_minTime = 0;
    best_path.clear();
  }

}  // namespace

// スコア計算
ll CalcScore(ll time) {
  ll res = round(100000000.0 / (1000.0 + time));
  return res;
}

ll compute_path_time() {
  ll timeSum = 0;
  srep(i, 1, current_path.size()) {
    timeSum += abs(current_path[i].first - current_path[i - 1].first) + abs(current_path[i].second - current_path[i - 1].second);
  }
  return timeSum;
}

inline bool is_inside_rect(int ite, int L, int R, int U, int D) {
  if (start_x[ite] < L || R < start_x[ite]) return false;
  if (start_y[ite] < U || D < start_y[ite]) return false;
  if (goal_x[ite] < L || R < goal_x[ite]) return false;
  if (goal_y[ite] < U || D < goal_y[ite]) return false;
  return true;
}

void filter_input_points() {
  vector<int> aa, bb, cc, dd;
  int mi = 1001001;
  int loop = 0;
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      if (get_elapsed_time() > 0.2) {
        break;
      }
    }

    vector<int> aaa, bbb, ccc, ddd, argAA;
    int L1 = Rand() % 801;
    int R1 = Rand() % 801;
    int U1 = Rand() % 801;
    int D1 = Rand() % 801;
    if (L1 > R1) swap(L1, R1);
    if (U1 > D1) swap(U1, D1);

    int sz = (R1 - L1) + (D1 - U1);
    if (sz >= mi) {
      continue;
    }

    rep(i, n) {
      if (is_inside_rect(i, L1, R1, U1, D1)) {
        aaa.push_back(start_x[i]);
        bbb.push_back(start_y[i]);
        ccc.push_back(goal_x[i]);
        ddd.push_back(goal_y[i]);
        argAA.push_back(i);
      }
    }
    if (aaa.size() >= m + 1) {
      aa = aaa;
      bb = bbb;
      cc = ccc;
      dd = ddd;
      argA = argAA;
      mi = sz;
    }
  }

  if (aa.size() >= m) {
    start_x = aa;
    start_y = bb;
    goal_x = cc;
    goal_y = dd;
    n = start_x.size();
  }
}

///////////////////////////////////////////////////////////////////////////////
// 始まったらやること
// 1. 入力部
// 2. maxScoreとans
// 3. 出力部
// 4. 愚直解
// 5. スコア計算関数
// 6. 貪欲解
// 7. 山登り
// 8. 焼きなまし
///////////////////////////////////////////////////////////////////////////////
int Solve(int mode) {
  start_timer();

  // 入力部
  string fileNameIfs = "./in/0000.txt";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    rep(i, n) {
      cin >> start_x[i] >> start_y[i] >> goal_x[i] >> goal_y[i];
    }
  }
  else {  // ファイル入力する
    rep(i, n) {
      ifs >> start_x[i] >> start_y[i] >> goal_x[i] >> goal_y[i];
    }
  }
  ifs.close();

  filter_input_points();

  // 愚直解
  current_path.push_back(P(400, 400));
  rep(i, m) {
    current_path.push_back(P(start_x[i], start_y[i]));
    argAns[i] = i;
    use[i] = 1;
  }
  rep(i, m) {
    current_path.push_back(P(goal_x[i], goal_y[i]));
  }
  current_path.push_back(P(400, 400));
  minTime = compute_path_time();
  maxScore = CalcScore(minTime);

  best_path = current_path;
  real_minTime = minTime;
  real_maxScore = maxScore;

  // 山登り解、焼きなまし解
  double TL = 1.8;
  double start_temp = 48;
  double end_temp = 0.0001;
  int loop = 0;
  double now_time = get_elapsed_time();
  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
      if (now_time > 0.4) {
        break;
      }
    }

    int x = Rand() % n;

    while (use[x]) {
      x = Rand() % n;
    }

    int ite = Rand() % m;

    int keepNum = argAns[ite];
    int diff = 0;
    diff += abs(start_x[x] - current_path[ite].first) + abs(start_y[x] - current_path[ite].second);
    diff += abs(start_x[x] - current_path[ite + 2].first) + abs(start_y[x] - current_path[ite + 2].second);
    diff += abs(goal_x[x] - current_path[ite + m].first) + abs(goal_y[x] - current_path[ite + m].second);
    diff += abs(goal_x[x] - current_path[ite + m + 2].first) + abs(goal_y[x] - current_path[ite + m + 2].second);

    diff -= abs(current_path[ite + 1].first - current_path[ite].first) + abs(current_path[ite + 1].second - current_path[ite].second);
    diff -= abs(current_path[ite + 1].first - current_path[ite + 2].first) + abs(current_path[ite + 1].second - current_path[ite + 2].second);
    diff -= abs(current_path[ite + m + 1].first - current_path[ite + m].first) + abs(current_path[ite + m + 1].second - current_path[ite + m].second);
    diff -= abs(current_path[ite + m + 1].first - current_path[ite + m + 2].first) + abs(current_path[ite + m + 1].second - current_path[ite + m + 2].second);

    int tmpTime = minTime + diff;
    int diffScore = -diff;

    double temp = start_temp + (end_temp - start_temp) * now_time / TL;
    double prob = exp((double)diffScore / temp);
    if (prob > Rand01()) {
      minTime = tmpTime;
      use[argAns[ite]] = 0;
      use[x] = 1;
      argAns[ite] = x;
      current_path[ite + 1] = P(start_x[x], start_y[x]);
      current_path[ite + m + 1] = P(goal_x[x], goal_y[x]);
      maxScore = CalcScore(minTime);
      if (maxScore > real_maxScore) {
        real_minTime = minTime;
        real_maxScore = maxScore;
        best_path = current_path;
        real_argAns = argAns;
        real_use = use;
      }
    }
    else {
      // 元に戻す
      ;
    }
  }

  // 最高スコアを戻す
  current_path = best_path;
  minTime = real_minTime;
  maxScore = real_maxScore;
  argAns = real_argAns;
  use = real_use;

  rep(i, m) {
    pair_[argAns[i]] = i + m;
    pair_[argAns[i] + 1000] = i;
  }
  rep(i, m) {
    argAns[i + m] = argAns[i] + 1000;
  }
  real_argAns = argAns;
  real_pair_ = pair_;

  argAns2 = argAns;
  real_argAns2 = argAns2;

  // それぞれをTSP
  rep(ui_tei, 10) {
    while (true) {
      loop++;
      if (loop % 100 == 1) {
        now_time = get_elapsed_time();
      }
      if (now_time > (TL / 15) * (6 + ui_tei)) break;

      int ite1 = Rand() % (m * 2);
      int ite2 = Rand() % (m * 2);
      while (ite1 == ite2) {
        ite2 = Rand() % (m * 2);
      }
      if (ite1 > ite2) swap(ite1, ite2);
      if (ite2 - ite1 == 1) continue;

      if (argAns[ite1] + 1000 == argAns[ite2]) continue;
      if (argAns[ite1] < 1000) {
        int pa = pair_[argAns[ite1]];
        if (pa < ite2) {
          continue;
        }
      }
      if (argAns[ite2] >= 1000) {
        int pa = pair_[argAns[ite2]];
        if (ite1 < pa) {
          continue;
        }
      }

      int diff = 0;

      if (true || loop % 2 == 0) {
        diff += abs(current_path[ite1 + 1].first - current_path[ite2].first) + abs(current_path[ite1 + 1].second - current_path[ite2].second);
        diff += abs(current_path[ite1 + 1].first - current_path[ite2 + 2].first) + abs(current_path[ite1 + 1].second - current_path[ite2 + 2].second);
        diff -= abs(current_path[ite1 + 1].first - current_path[ite1].first) + abs(current_path[ite1 + 1].second - current_path[ite1].second);
        diff -= abs(current_path[ite1 + 1].first - current_path[ite1 + 2].first) + abs(current_path[ite1 + 1].second - current_path[ite1 + 2].second);

        diff += abs(current_path[ite2 + 1].first - current_path[ite1].first) + abs(current_path[ite2 + 1].second - current_path[ite1].second);
        diff += abs(current_path[ite2 + 1].first - current_path[ite1 + 2].first) + abs(current_path[ite2 + 1].second - current_path[ite1 + 2].second);
        diff -= abs(current_path[ite2 + 1].first - current_path[ite2].first) + abs(current_path[ite2 + 1].second - current_path[ite2].second);
        diff -= abs(current_path[ite2 + 1].first - current_path[ite2 + 2].first) + abs(current_path[ite2 + 1].second - current_path[ite2 + 2].second);

        int tmpTime = minTime + diff;
        int diffScore = -diff;

        double temp = start_temp + (end_temp - start_temp) * now_time / TL;
        double prob = exp((double)diffScore / temp);
        if (prob > Rand01()) {
          minTime = tmpTime;
          if (argAns[ite1] < 1000) {
            pair_[argAns[ite1] + 1000] = ite2;
          }
          else {
            pair_[argAns[ite1] - 1000] = ite2;
          }
          if (argAns[ite2] < 1000) {
            pair_[argAns[ite2] + 1000] = ite1;
          }
          else {
            pair_[argAns[ite2] - 1000] = ite1;
          }
          swap(argAns[ite1], argAns[ite2]);
          swap(current_path[ite1 + 1], current_path[ite2 + 1]);
          maxScore = CalcScore(minTime);
          if (maxScore > real_maxScore) {
            real_minTime = minTime;
            real_maxScore = maxScore;
            best_path = current_path;
            real_argAns = argAns;
            real_use = use;
            real_pair_ = pair_;
          }
        }
        else {
          // 元に戻す
          ;
        }
      }
    }

    // 最高スコアを戻す
    current_path = best_path;
    minTime = real_minTime;
    maxScore = real_maxScore;
    argAns = real_argAns;
    argAns2 = real_argAns2;
    use = real_use;
    pair_ = real_pair_;
  }

  TL = 1.95;

  const int INF = 1001001001;
  int mi = INF;
  argAns2 = argAns;
  vector<P> ans2 = current_path;

  random_device seed_gen;
  mt19937 engine(seed_gen());

  vector<int> v;
  rep(i, m) {
    v.push_back(i);
  }

  while (true) {
    loop++;
    if (loop % 100 == 1) {
      now_time = get_elapsed_time();
    }
    if (now_time > TL) break;

    std::shuffle(v.begin(), v.end(), engine);
    vector<int> vv;
    rep(i, 50) vv.push_back(v[i]);
    sort(vv.begin(), vv.end());

    set<int> useIte;
    int cnt = 0;
    int now = 0;
    rep(i, m * 2) {
      if (argAns[i] < 1000) {
        if (cnt == vv[now]) {
          useIte.insert(argAns[i]);
          now++;
        }
        cnt++;
      }
    }

    int tmpTime = 0;
    int nowx = 400;
    int nowy = 400;
    vector<P> ans3;
    ans3.push_back(P(400, 400));
    rep(i, m * 2) {
      int ite = argAns[i];
      if (ite >= 1000) ite -= 1000;
      if (useIte.find(ite) != useIte.end()) {
        ans3.push_back(current_path[i + 1]);
        tmpTime += abs(nowx - current_path[i + 1].first) + abs(nowy - current_path[i + 1].second);
        nowx = current_path[i + 1].first;
        nowy = current_path[i + 1].second;
      }
    }
    ans3.push_back(P(400, 400));

    if (tmpTime < mi) {
      mi = tmpTime;
      ans2 = ans3;
      argAns2.clear();

      for (auto ite : useIte) {
        argAns2.push_back(ite);
      }
    }
  }

  current_path = ans2;
  argAns = argAns2;

  minTime = compute_path_time();
  maxScore = CalcScore(minTime);

  // 解の出力
  if (mode == 0) {
    cout << 50;
    rep(i, 50) {
      if (argAns[i] < 1000) {
        cout << " " << argA[argAns[i]] + 1;
      }
    }
    cout << endl;
    cout << current_path.size();
    rep(i, current_path.size()) {
      cout << " " << current_path[i].first << " " << current_path[i].second;
    }
    cout << endl;
  }

  // デバッグ用
  if (mode != 0) {
    cout << "loop = " << loop << endl;
    cout << maxScore << endl;
    cout << get_elapsed_time() << "sec." << endl;
  }

  // ファイル出力
  if (true) {
    string fileNameOfs = "0000_out.txt";
    ofstream ofs(fileNameOfs);

    ofs << 50;
    rep(i, 50) {
      if (argAns[i] < 1000) {
        ofs << " " << argA[argAns[i]] + 1;
      }
    }
    ofs << endl;
    ofs << current_path.size();
    rep(i, current_path.size()) {
      ofs << " " << current_path[i].first << " " << current_path[i].second;
    }
    ofs << endl;

    ofs.close();
  }

  return maxScore;
}

int main() {
  int mode = 1;
  m = 55;
  if (mode == 0) {
    Solve(mode);
  }
  else if (mode == 1) {
    Solve(mode);
  }

  return 0;
}
