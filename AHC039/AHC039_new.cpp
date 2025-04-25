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
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

typedef long long int ll;

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;
  void start_timer() { start_time_clock = std::chrono::steady_clock::now(); }
  double get_elapsed_time() { std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock; return elapsed.count(); }
}

// 乱数
namespace
{
  static uint32_t rand_xorshift() {
    static uint32_t x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    uint32_t t = x ^ (x << 11);
    x = y; y = z; z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
  }
  static double rand_01() { return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX); }
  static double rand_range(double l, double r) { return l + (r - l) * rand_01(); }
  static uint32_t rand_range(uint32_t l, uint32_t r) { return l + rand_xorshift() % (r - l + 1); }
  void shuffle_array(int* arr, int n) { for (int i = n - 1; i >= 0; i--) { int j = rand_xorshift() % (i + 1); swap(arr[i], arr[j]); } }
}

struct Point { int x, y; };

enum class FishType : uint8_t { MACKEREL = 0, SARDINE = 1 };
struct Fish { Point p; FishType type; };

// 入力データ
int N_fish_half;                   // 「サバ」と「イワシ」それぞれの匹数 (最大 1000)
vector<Fish> fishes;               // サイズ 2*N_fish_half

// 出力データ（ポリゴン）
vector<Point> polygon;             // 時計回り（反時計でも OK）頂点数 m (m ≤ 1000)

// ─────────────────────────────────────────────
// Field : 矩形内のサバ/イワシ数を数える簡易実装
//   ※ N ≤ 2000 なので毎回 O(N) で十分高速
// ─────────────────────────────────────────────
class Field
{
  const vector<Fish>& fs;
public:
  explicit Field(const vector<Fish>& f) : fs(f) {}

  /**
   * @brief [x1,x2]×[y1,y2] に含まれる (サバ,イワシ) を数える
   */
  pair<int, int> countRect(int x1, int y1, int x2, int y2) const {
    int a = 0, b = 0;
    for (const auto& f : fs) {
      if (f.p.x >= x1 && f.p.x <= x2 && f.p.y >= y1 && f.p.y <= y2) {
        (f.type == FishType::MACKEREL ? a : b)++;
      }
    }
    return { a,b };
  }

  vector<Point> getMackerelCoords() const {
    vector<Point> res; res.reserve(fs.size() / 2);
    for (const auto& f : fs) if (f.type == FishType::MACKEREL) res.push_back(f.p);
    return res;
  }
};

// ─────────────────────────────────────────────
// Solver : とりあえず全サバを覆う最小長方形を返す簡易版
//           必要に応じて SA などに差し替え可能
// ─────────────────────────────────────────────
class Solver
{
  const Field& fld;
  vector<Point> bestPoly;
  ll bestScore = -1e18;
public:
  explicit Solver(const Field& f) : fld(f) {}

  /** sc = max(0, a - b + 1) */
  static ll score(int a, int b) { return max(0, a - b + 1); }

  /** 初期長方形を作成 */
  void buildInitial() {
    auto mac = fld.getMackerelCoords();
    if (mac.empty()) {
      bestPoly = { {0,0},{0,1},{1,1},{1,0} };
      bestScore = 1; return;
    }
    int minx = INT_MAX, miny = INT_MAX, maxx = INT_MIN, maxy = INT_MIN;
    for (auto& p : mac) {
      minx = min(minx, p.x); maxx = max(maxx, p.x);
      miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    if (minx == maxx) maxx++; // 幅ゼロ回避
    if (miny == maxy) maxy++;
    auto cnt_ab = fld.countRect(minx, miny, maxx, maxy);
    int a_m = cnt_ab.first;
    int b_m = cnt_ab.second;
    bestScore = score(a_m, b_m);
    bestPoly = { {minx,miny},{minx,maxy},{maxx,maxy},{maxx,miny} };
  }

  /**
   * @return 得られたポリゴン
   */
  vector<Point> solve(double timelimit_sec) {
    buildInitial();

    // --- Greedy 1‑step edge shrink to try to remove sardines without losing mackerel
    auto mac = fld.getMackerelCoords();
    int left = bestPoly[0].x, bottom = bestPoly[0].y;
    int right = bestPoly[2].x, top = bestPoly[2].y; // 座標系: y 上向きでも下向きでも長方形なら OK

    auto allMacInside = [&](int l, int b, int r, int t) {
      for (auto& p : mac) { if (p.x<l || p.x>r || p.y<b || p.y>t) return false; } return true; };

    while (get_elapsed_time() < timelimit_sec) {
      bool improved = false;
      // 4 方向それぞれ 1 ステップずつ縮めてみる
      array<tuple<int, int, int, int>, 4> cand = {
        make_tuple(left + 1,bottom,right,top),
        make_tuple(left,bottom + 1,right,top),
        make_tuple(left,bottom,right - 1,top),
        make_tuple(left,bottom,right,top - 1)
      };
      for (const auto& tpl : cand) {
        int cl, cb, cr, ct;
        std::tie(cl, cb, cr, ct) = tpl;   // tuple → 個別変数へ展開
        if (cl >= cr || cb >= ct) continue;
        if (!allMacInside(cl, cb, cr, ct)) continue;
        auto cnt = fld.countRect(cl, cb, cr, ct);
        int a_m = cnt.first;
        int sard_m = cnt.second;
        ll sc = score(a_m, sard_m);
        if (sc > bestScore) {
          bestScore = sc;
          left = cl; bottom = cb; right = cr; top = ct;
          improved = true;
          break;
        }
      }
      if (!improved) break;
    }
    bestPoly = { {left,bottom},{left,top},{right,top},{right,bottom} };
    return bestPoly;
  }
};

const double TIME_LIMIT = 1.8;
int exec_mode;

/**
 * @brief 標準入力 or ファイルから問題入力を読み込む。
 * @details 読み込んだ結果はグローバル変数 fishes に格納。
 */
void input_data(int case_num) {
  fishes.clear();

  // どこから読むか決定
  istream* pin = &cin;
  ifstream ifs;
  if (exec_mode != 0) {
    ostringstream oss;
    oss << "./in/" << setw(4) << setfill('0') << case_num << ".txt";
    ifs.open(oss.str());
    if (ifs.is_open()) pin = &ifs; // ファイルが存在すればそちらを優先
  }

  // ここから読み取り開始
  if (!(*pin)) return; // 失敗
  (*pin) >> N_fish_half;
  fishes.reserve(2 * N_fish_half);

  for (int i = 0; i < 2 * N_fish_half; ++i) {
    int x, y; (*pin) >> x >> y;
    fishes.push_back({ {x, y}, (i < N_fish_half ? FishType::MACKEREL : FishType::SARDINE) });
  }
}

/**
 * @brief ポリゴンを (標準 or ファイル) に出力。
 * @details polygon が空の場合はダミーの 1×1 正方形を出力。
 */
void output_data(int case_num) {
  if (polygon.empty()) {
    // TODO: アルゴリズム実装後は消去。暫定ダミーポリゴン
    polygon = { {0,0},{0,1},{1,1},{1,0} };
  }

  ostream* pout = &cout;
  ofstream ofs;
  if (exec_mode != 0) {
    ostringstream oss;
    oss << "./out/" << setw(4) << setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
    if (ofs.is_open()) pout = &ofs;
  }

  // 出力フォーマット: m & 各頂点
  (*pout) << polygon.size() << '\n';
  for (auto& pt : polygon) (*pout) << pt.x << ' ' << pt.y << '\n';
}

// ─────────────────────────────────────────────
// スコア計算（オフラインデバッグ用）
// ※ 本番では呼ばれないので未実装のままでも OK
// ─────────────────────────────────────────────
ll calculate_score() {
  // TODO: 必要なら実装。ここでは 0 を返す。
  return 0;
}

// ─────────────────────────────────────────────
// 1 ケース解くラッパ関数
// ─────────────────────────────────────────────
ll solve_case(int case_num) {
  start_timer();
  polygon.clear();

  input_data(case_num);

  Field fld(fishes);
  Solver solver(fld);
  polygon = solver.solve(TIME_LIMIT);

  output_data(case_num);

  ll score = 0;
  if (exec_mode != 0) score = calculate_score();
  return score;
}

// ─────────────────────────────────────────────
// main – テンプレートに合わせて main_new → main に変更
// ─────────────────────────────────────────────
int main() {
  exec_mode = 2; // 0 / 1 / 2 を状況に応じて変更

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
      ll score = solve_case(i);
      sum_score += score;
      cerr << "case = " << setw(2) << i << ", "
        << "score = " << setw(6) << score << ", "
        << "sum = " << setw(7) << sum_score << ", "
        << "time = " << fixed << setprecision(2) << get_elapsed_time() << "\n";
    }
  }
  return 0;
}
