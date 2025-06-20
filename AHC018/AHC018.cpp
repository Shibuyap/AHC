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

//------------------------------------------------------------------------------
// (1) Union-Find関連：unnamed namespace
//------------------------------------------------------------------------------
namespace
{
  // 定数 (UPPER_SNAKE_CASE)
  const int UNION_FIND_MAX = 50000;

  // Union-Findで使用する配列 (snake_case)
  int uf_parent[UNION_FIND_MAX];
  int uf_rank[UNION_FIND_MAX];
  int uf_count[UNION_FIND_MAX];

  // n要素で初期化
  void uf_init()
  {
    for (int i = 0; i < UNION_FIND_MAX; i++) {
      uf_parent[i] = i;
      uf_rank[i] = 0;
      uf_count[i] = 1;
    }
  }

  // 木の根を求める
  int uf_find(int x)
  {
    if (uf_parent[x] == x) {
      return x;
    }
    else {
      return uf_parent[x] = uf_find(uf_parent[x]);
    }
  }

  // xとyの属する集合を併合
  void uf_unite(int x, int y)
  {
    x = uf_find(x);
    y = uf_find(y);
    if (x == y) { return; }

    if (uf_rank[x] < uf_rank[y]) {
      uf_parent[x] = y;
      uf_count[y] += uf_count[x];
    }
    else {
      uf_parent[y] = x;
      uf_count[x] += uf_count[y];
      if (uf_rank[x] == uf_rank[y]) uf_rank[x]++;
    }
  }

  // xとyが同じ集合に属するか否か
  bool uf_same(int x, int y)
  {
    return uf_find(x) == uf_find(y);
  }
}

//------------------------------------------------------------------------------
// (2) 便利系・定数関連：unnamed namespace
//------------------------------------------------------------------------------
namespace
{
  const int INF = 1001001001;

  // 移動用の差分 (UPPER_SNAKE_CASE)
  const int DELTA_X[4] = { -1, 0, 1, 0 };
  const int DELTA_Y[4] = { 0, -1, 0, 1 };

  // (例) 移動方向を示す文字列 (使っていない可能性あり)
  const char MOVE_CHARS[4] = { 'U', 'L', 'D', 'R' };
}

//------------------------------------------------------------------------------
// (3) 乱数ライブラリ：unnamed namespace
//------------------------------------------------------------------------------
namespace
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
}

//------------------------------------------------------------------------------
// (4) 主要グローバル変数：unnamed namespace
//------------------------------------------------------------------------------
namespace
{
  // もとのコードの MIN_ROCK, MAX_ROCK は使っていないため置いておきます
  const int MIN_ROCK = 10;
  const int MAX_ROCK = 5000;

  // 旧: CC (攻撃コストの候補)
  const int C_OPTIONS[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  // 実行モード (0: 提出用 / それ以外: ローカル)
  int mode = 0;
  ofstream ofs;

  // 入力関連 (旧N, W, K, C)
  int n;          // N
  int w;         // W
  int k;         // K
  int attack_cost;     // C

  // 水源の座標 (旧: a[i], b[i])
  int water_x[20];
  int water_y[20];

  // 家の座標 (旧: c[i], d[i])
  int house_x[20];
  int house_y[20];

  // 岩の強度マップ (旧: S)
  int rock_strength[210][210];

  // 各マスが壊れているかどうかなど (旧: f, minS, maxS)
  int is_broken[210][210];
  int min_strength[210][210];
  int max_strength[210][210];

  // 各家が最も近い水源を記録 (旧: nearX, nearY)
  int nearest_water_x[20];
  int nearest_water_y[20];

  // スコア計算用: HP (旧: hp)
  int health_points;

  // 攻撃力の候補 (旧: attack_power[8])
  const int ATTACK_POWER_VALUES[8] = { 50, 50, 50, 50, 50, 100, 100, 100 };

  // 実際に使う攻撃力 (旧: ATTACK_POWER)
  int attack_power_global = 100;
}

//------------------------------------------------------------------------------
// (5) ヘルパー関数
//------------------------------------------------------------------------------
inline bool is_out_of_bounds(int x, int y)
{
  return (x < 0 || x >= n || y < 0 || y >= n);
}

// x,yが水とつながっているか (旧: IsUniteWater)
inline bool is_united_with_water(int x, int y)
{
  // UF_MAX - 1 = UNION_FIND_MAX - 1
  // x * n + y が UNION_FIND_MAX - 1 と同じ連結成分か
  return uf_same(x * n + y, UNION_FIND_MAX - 1);
}

// マンハッタン距離 (旧: Manhattan)
inline int manhattan_distance(int x1, int y1, int x2, int y2)
{
  return abs(x1 - x2) + abs(y1 - y2);
}

//------------------------------------------------------------------------------
// (6) スコア計算 (旧: CalcScore)
//------------------------------------------------------------------------------
int calc_score()
{
  return health_points;
}

//------------------------------------------------------------------------------
// (7) 攻撃処理 (旧: Attack)
//------------------------------------------------------------------------------
int attack_cell(int x, int y, int power)
{
  // min_strengthを今のmax_strengthに更新
  min_strength[x][y] = max_strength[x][y];
  // max_strengthに今回の追加攻撃量を加算
  max_strength[x][y] += power;
  // HP(スコア)を加算
  health_points += attack_cost + power;

  int res = 0;
  if (mode == 0) {
    // 提出用: 標準出力に攻撃内容を出し、結果(壊れたかどうか)を受け取る
    cout << x << ' ' << y << ' ' << power << endl;
    fflush(stdout);
    cin >> res;
  }
  else {
    // ローカル実行: rock_strengthを減らして壊れたか判定
    ofs << x << ' ' << y << ' ' << power << endl;
    rock_strength[x][y] -= power;
    if (rock_strength[x][y] <= 0) {
      res = 1; // 壊れた
    }
  }

  // 壊れた場合
  if (res != 0) {
    is_broken[x][y] = 1;

    // 周囲マスがすでに壊れていれば Union-Find でつなぐ
    for (int i = 0; i < 4; ++i) {
      int nx = x + DELTA_X[i];
      int ny = y + DELTA_Y[i];
      if (is_out_of_bounds(nx, ny)) { continue; }
      if (is_broken[nx][ny] != 0) {
        uf_unite(x * n + y, nx * n + ny);
      }
    }

    // もし水源マスだったら、UF_MAX-1 とつなぐ
    for (int i = 0; i < w; ++i) {
      if (x == water_x[i] && y == water_y[i]) {
        uf_unite(x * n + y, UNION_FIND_MAX - 1);
      }
    }

    // ここで再度「全ての家が繋がったかどうか」を判定しているらしい
    res = 2;
    for (int i = 0; i < k; ++i) {
      if (!is_united_with_water(house_x[i], house_y[i])) {
        res = 1; // まだ未接続の家がある
      }
    }
  }

  return res;
}

//------------------------------------------------------------------------------
// (8) 連続攻撃で壊れるまで掘り続ける (旧: Challenge)
//------------------------------------------------------------------------------
int challenge_cell(int x, int y, int power)
{
  int res = is_broken[x][y];
  while (res == 0) {
    res = attack_cell(x, y, power);
  }
  return res;
}

//------------------------------------------------------------------------------
// (9) 複数ケース用にすべてクリア (旧: AllClear_MultiCase)
//------------------------------------------------------------------------------
void clear_all_multicase()
{
  // 実装無し(サンプル)
}

//------------------------------------------------------------------------------
// (10) 初期状態作成 (旧: Init)
//------------------------------------------------------------------------------
void init(int problem_num)
{
  uf_init();
  health_points = 0;

  // 状態クリア
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      is_broken[i][j] = 0;
      min_strength[i][j] = 0;
      max_strength[i][j] = 0;
    }
  }

  // 攻撃コストCに応じて攻撃力を選択
  for (int i = 0; i < 8; ++i) {
    if (attack_cost == C_OPTIONS[i]) {
      attack_power_global = ATTACK_POWER_VALUES[i];
    }
  }

  // ローカル実行時はファイル出力をオープン
  if (mode != 0) {
    string file_name_ofs = "./out/";
    {
      // problem_numを4桁化してファイル名を作る
      int tmp = problem_num;
      string str_num;
      for (int i = 0; i < 4; ++i) {
        str_num += char((tmp % 10) + '0');
        tmp /= 10;
      }
      reverse(str_num.begin(), str_num.end());
      file_name_ofs += str_num + ".txt";
    }

    ofs.open(file_name_ofs);
  }
}

//------------------------------------------------------------------------------
// (11) 入力受け取り (旧: Input)
//------------------------------------------------------------------------------
void read_input(int problem_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problem_num << ".txt";
  ifstream ifs(oss.str());

  // 実行モード=0またはファイルが開けなかったら標準入力から読む
  if (mode == 0 || !ifs.is_open()) {
    cin >> n >> w >> k >> attack_cost;
    for (int i = 0; i < w; ++i) {
      cin >> water_x[i] >> water_y[i];
    }
    for (int i = 0; i < k; ++i) {
      cin >> house_x[i] >> house_y[i];
    }
  }
  else {
    // ローカル時はファイル入力
    ifs >> n >> w >> k >> attack_cost;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        ifs >> rock_strength[i][j];
      }
    }
    for (int i = 0; i < w; ++i) {
      ifs >> water_x[i] >> water_y[i];
    }
    for (int i = 0; i < k; ++i) {
      ifs >> house_x[i] >> house_y[i];
    }
  }
}

//------------------------------------------------------------------------------
// (12) 解答出力 (旧: Output)
//------------------------------------------------------------------------------
void write_output(int /*problem_num*/)
{
  // ローカル実行時はファイル出力クローズ
  if (mode != 0) {
    ofs.close();
  }
}

//------------------------------------------------------------------------------
// (13) 実際の処理 (旧: Solve)
//------------------------------------------------------------------------------
int solve_problem(int problem_num = 0)
{
  clock_t start_time = clock();
  clock_t end_time = clock();
  (void)start_time; // 使わないならvoidキャストで消す
  (void)end_time;

  // 初期化
  init(problem_num);

  // 各家の一番近い水源を探す
  for (int i = 0; i < k; ++i) {
    int dist = INF;
    for (int j = 0; j < w; ++j) {
      int mdist = manhattan_distance(house_x[i], house_y[i], water_x[j], water_y[j]);
      if (mdist < dist) {
        dist = mdist;
        nearest_water_x[i] = water_x[j];
        nearest_water_y[i] = water_y[j];
      }
    }
  }

  // 家を水源につなげる
  for (int i = 0; i < k; ++i) {
    int phase = 0;
    int now_x = house_x[i], now_y = house_y[i];
    int next_x = -1, next_y = -1;
    int dir = -1;

    // その家がまだ水と繋がっていない間は掘り続ける
    while (!is_united_with_water(house_x[i], house_y[i])) {
      if (phase == 0) {
        // 家のマスを掘る
        challenge_cell(now_x, now_y, attack_power_global);
        phase = 1;
      }
      else if (phase == 1) {
        // 次の「ターゲット地点(next_x, next_y)」を決める
        int diff_x = nearest_water_x[i] - now_x;
        if (abs(diff_x) > 20) {
          diff_x = (diff_x > 0) ? 20 : -20;
        }
        int diff_y = nearest_water_y[i] - now_y;
        if (abs(diff_y) > 20) {
          diff_y = (diff_y > 0) ? 20 : -20;
        }

        // 既に壊れて水路が通っているマスが近くにあれば、そちらを優先
        int diff_sum = abs(diff_x) + abs(diff_y);
        for (int k = 1; k < diff_sum; ++k) {
          bool found = false;
          for (int j = 0; j < 4; ++j) {
            int nx = now_x + DELTA_X[j] * k;
            int ny = now_y + DELTA_Y[j] * k;
            if (!is_out_of_bounds(nx, ny) && is_united_with_water(nx, ny)) {
              diff_x = nx - now_x;
              diff_y = ny - now_y;
              found = true;
              break;
            }
          }
          if (found) { break; }
        }

        // x方向 or y方向だけの場合
        if (diff_x == 0) {
          next_x = now_x;
          next_y = now_y + diff_y;
          challenge_cell(next_x, next_y, attack_power_global);
        }
        else if (diff_y == 0) {
          next_x = now_x + diff_x;
          next_y = now_y;
          challenge_cell(next_x, next_y, attack_power_global);
        }
        else {
          // 斜めに進むときはまず x方向に掘って、それから y方向に掘る or 逆順
          int nx1 = now_x + diff_x;
          int ny1 = now_y;
          challenge_cell(nx1, ny1, attack_power_global);

          int nx2 = now_x;
          int ny2 = now_y + diff_y;
          challenge_cell(nx2, ny2, attack_power_global);

          // どちらのマスのほうが max_strength が小さいかで次の現在地を決める
          if (max_strength[nx1][ny1] <= max_strength[nx2][ny2]) {
            next_x = nx1;
            next_y = ny1;
          }
          else {
            next_x = nx2;
            next_y = ny2;
          }
        }

        // 移動方向dirを決める (DELTA_X/DELTA_Yの添字を対応)
        if (next_x - now_x < 0) {
          dir = 0; // 上方向(DELTA_X[0],DELTA_Y[0])
        }
        else if (next_x - now_x > 0) {
          dir = 2; // 下方向(DELTA_X[2],DELTA_Y[2])
        }
        else if (next_y - now_y < 0) {
          dir = 1; // 左方向(DELTA_X[1],DELTA_Y[1])
        }
        else if (next_y - now_y > 0) {
          dir = 3; // 右方向(DELTA_X[3],DELTA_Y[3])
        }

        phase = 2;
      }
      else if (phase == 2) {
        // next_x, next_y に到達済なら再度ターゲットを再計算
        if (now_x == next_x && now_y == next_y) {
          phase = 1;
        }
        else {
          // next_x, next_y に向かって一歩進む
          int nx = now_x + DELTA_X[dir];
          int ny = now_y + DELTA_Y[dir];
          if (max_strength[nx][ny] == 0) {
            // まだ攻撃されていないマスなら，周囲の状況に応じた攻撃量計算
            int d1 = manhattan_distance(now_x, now_y, nx, ny);
            int d2 = manhattan_distance(nx, ny, next_x, next_y);
            int p1 = min_strength[now_x][now_y];
            int p2 = min_strength[next_x][next_y];
            int power = (d2 * p1 + d1 * p2) / (d1 + d2);
            power = max(power, 10);
            attack_cell(nx, ny, power);
          }
          else {
            // すでに攻撃されたことのあるマスなら、壊れるまで掘る
            challenge_cell(nx, ny, attack_power_global);
            now_x = nx;
            now_y = ny;
          }
        }
      }
    }
  }

  // ローカル実行時はデバッグ出力
  if (mode != 0) {
    cerr << "problem_num = " << problem_num
      << ", health_points = " << health_points << endl;
  }

  return calc_score();
}

//------------------------------------------------------------------------------
// (14) ラッパ関数 (旧: SolveOuter)
//------------------------------------------------------------------------------
int solve_outer(int problem_num = 0)
{
  // 入力受け取り
  read_input(problem_num);

  // 実行
  int score = solve_problem(problem_num);

  // 出力
  write_output(problem_num);

  return score;
}

//------------------------------------------------------------------------------
// (15) main
//------------------------------------------------------------------------------
int main()
{
  // mode = 0 => 提出用
  mode = 2;

  if (mode == 0) {
    // 提出用: 単一ケース
    solve_outer();
  }
  else if (mode == 1) {
    // 1ケースのみファイル読み・書き
    solve_outer(0);
  }
  else if (mode == 2) {
    // 複数ケース
    ll score_sum = 0;
    for (int i = 0; i < 10; ++i) {
      score_sum += solve_outer(i);
      clear_all_multicase();
    }
    cout << "scoreSum = " << score_sum << endl;
  }

  return 0;
}
