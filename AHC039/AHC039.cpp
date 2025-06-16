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
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

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

// 二次元座標を表す構造体
struct Point
{
public:
  int x;
  int y;

  Point() { x = 0; y = 0; }                  // デフォルトコンストラクタ
  Point(int _x, int _y) { x = _x; y = _y; }  // 座標を指定するコンストラクタ
};

const int MAX_N = 30;  // 未使用の定数（今後の拡張用？）

const int n = 5000;    // 魚の数（サバとイワシそれぞれの数）

vector<Point> saba, iwashi;  // サバとイワシの座標を格納するベクター

vector<Point> ans;  // 出力するポリゴンの頂点座標を格納するベクター

int ansScore;       // 現在のスコア
int best_ansScore;  // 最良のスコア（未使用？）

int f[510][510];
int best_f[510][510];

// 現在の解答を最良の解答として保存する関数（未使用）
void CopyToBest(int blockSize)
{
  best_ansScore = ansScore;
  rep(i, blockSize + 2)
  {
    rep(j, blockSize + 2)
    {
      best_f[i][j] = f[i][j];
    }
  }
}

// 最良の解答を現在の解答として設定する関数（未使用）
void CopyToAns(int blockSize)
{
  ansScore = best_ansScore;
  rep(i, blockSize + 2)
  {
    rep(j, blockSize + 2)
    {
      f[i][j] = best_f[i][j];
    }
  }
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;   // スコアの初期化
  rep(i, 510)rep(j, 510)best_f[i][j] = 0;

  ans.clear();    // 解答の初期化
  saba.clear();   // サバの座標の初期化
  iwashi.clear(); // イワシの座標の初期化
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());  // ファイルストリームを開く

  saba.resize(n);    // サバのベクターをリサイズ
  iwashi.resize(n);  // イワシのベクターをリサイズ

  if (!ifs.is_open()) {
    // 標準入力からの読み込み
    int _n;
    cin >> _n;  // 魚の数（未使用）
    rep(i, n) cin >> saba[i].x >> saba[i].y;     // サバの座標を読み込む
    rep(i, n) cin >> iwashi[i].x >> iwashi[i].y; // イワシの座標を読み込む
  }
  else {
    // ファイルからの読み込み
    int _n;
    ifs >> _n;  // 魚の数（未使用）
    rep(i, n) ifs >> saba[i].x >> saba[i].y;     // サバの座標を読み込む
    rep(i, n) ifs >> iwashi[i].x >> iwashi[i].y; // イワシの座標を読み込む
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    // モードが0以外の場合、出力ファイルを作成
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
ll CalcScore()
{
  ll res = ansScore + 1;  // 問題文の得点計算式に基づく（max(0, a - b + 1)）
  return res;             // スコアを返す
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力に出力
    cout << ans.size() << endl;  // ポリゴンの頂点数を出力
    for (auto p : ans) cout << p.x << ' ' << p.y << endl;  // 各頂点の座標を出力
  }
  else {
    // ファイルに出力
    ofs << ans.size() << endl;  // ポリゴンの頂点数を出力
    for (auto p : ans) ofs << p.x << ' ' << p.y << endl;  // 各頂点の座標を出力
  }
}

// ポリゴンの辺の総和が制約を満たしているか確認する関数
bool IsLengthOK(vector<Point> vp)
{
  int len = 0;
  // 各辺の長さを計算
  rep(i, vp.size() - 1)
  {
    len += abs(vp[i + 1].x - vp[i].x);  // x座標の差の絶対値を加算
    len += abs(vp[i + 1].y - vp[i].y);  // y座標の差の絶対値を加算
  }

  // 最後の頂点と最初の頂点を結ぶ辺の長さを加算
  len += abs(vp[0].x - vp.back().x);
  len += abs(vp[0].y - vp.back().y);

  // 総和が400,000以下であればtrueを返す
  return len <= 400000;
}

// 座標がグリッドの範囲外かを確認する関数
bool IsNG(int x, int y, int blockSize)
{
  if (x < 0 || blockSize <= x || y < 0 || blockSize <= y) return true;  // 範囲外ならtrue
  return false;  // 範囲内ならfalse
}

int block[510][510];  // 各グリッドセルのスコアを格納する配列
void InitBlock(const int blockSize)
{
  rep(i, blockSize) rep(j, blockSize) block[i][j] = 0;  // 初期化

  // サバとイワシの位置から、各セルのスコアを計算
  rep(i, n)
  {
    {
      // サバの座標をセルに割り当て、スコアを加算
      int xx = saba[i].x / (100000 / blockSize);
      xx = min(xx, blockSize - 1);  // 最大値を超えないように調整
      int yy = saba[i].y / (100000 / blockSize);
      yy = min(yy, blockSize - 1);

      block[xx][yy]++;  // サバがいるセルのスコアを+1
    }

    {
      // イワシの座標をセルに割り当て、スコアを減算
      int xx = iwashi[i].x / (100000 / blockSize);
      xx = min(xx, blockSize - 1);
      int yy = iwashi[i].y / (100000 / blockSize);
      yy = min(yy, blockSize - 1);

      block[xx][yy]--;  // イワシがいるセルのスコアを-1
    }
  }
}

int haba[510][510];  // 幅優先探索用の配列
void Method3_SA(const int xx1, const int xx2, const int yy1, const int yy2, const int blockSize, int& loop2, double timeLimit)
{
  queue<P> que;  // 幅優先探索のためのキュー

  double nowTime = GetNowTime();  // 現在の経過時間
  const double START_TEMP = 200.0;  // 焼きなまし法の開始温度
  const double END_TEMP = 0.1;      // 焼きなまし法の終了温度
  double temp = START_TEMP + (END_TEMP - START_TEMP) * nowTime / timeLimit;  // 現在の温度

  while (true) {
    if (loop2 % 1 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;  // 時間制限を過ぎたらループを抜ける
    }
    loop2++;

    int rax = Rand() % blockSize;  // ランダムなセルのxインデックス
    int ray = Rand() % blockSize;  // ランダムなセルのyインデックス

    int ng = 1;  // 変更が可能かどうかのフラグ
    rep(i, 4)
    {
      int nx = rax + dx[i];
      int ny = ray + dy[i];
      if (IsNG(nx, ny, blockSize)) continue;  // 範囲外は無視
      if (f[nx + 1][ny + 1] != f[rax + 1][ray + 1]) ng = 0;  // 隣接セルが異なる状態なら変更可能
    }

    if (ng) continue;  // 変更不可なら次のループへ

    int tmpScore = ansScore;
    if (f[rax + 1][ray + 1] == 0) {
      tmpScore += block[rax][ray];  // セルを追加した場合のスコア
    }
    else {
      tmpScore += -block[rax][ray];  // セルを削除した場合のスコア
    }

    const double progressRatio = nowTime / TL;  // 進捗率（0.0〜1.0）
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;  // 温度の更新

    double diff = tmpScore - ansScore;  // スコアの差分
    double prob = exp(diff / temp);     // 焼きなまし法の採用確率
    f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];  // 状態を反転
    int upd = 0;  // 解答を更新するかどうかのフラグ
    if (prob > Rand01()) {
      // 解答を更新する場合の処理
      // 幅優先探索で連結成分の数を確認
      rep(i, blockSize + 2) rep(j, blockSize + 2) haba[i][j] = 0;  // 初期化
      int now = 1;
      upd = 1;
      srep(i, 1, blockSize + 1)
      {
        srep(j, 1, blockSize + 1)
        {
          if (haba[i][j] != 0) continue;  // 既に探索済みならスキップ
          if (now == 3) {
            upd = 0;  // 連結成分が2つを超える場合は更新不可
            break;
          }

          haba[i][j] = now;
          que.push(P(i, j));
          while (que.size()) {
            int x = que.front().first;
            int y = que.front().second;
            que.pop();
            rep(k, 4)
            {
              int nx = x + dx[k];
              int ny = y + dy[k];
              if (IsNG(nx - 1, ny - 1, blockSize)) continue;
              if (haba[nx][ny] == 0 && f[nx][ny] == f[i][j]) {
                haba[nx][ny] = now;
                que.push(P(nx, ny));
              }
            }
          }

          now++;
        }
        if (upd == 0) break;
      }
      if (now > 3) upd = 0;  // 連結成分が3つ以上なら更新不可
    }

    if (upd) {
      auto ans2 = ans;  // 現在の解答を一時保存

      ans.clear();

      int sx = -1, sy = -1;
      int befx = -1, befy = -1;
      // 境界の始点を探す
      srep(i, 1, blockSize + 1)
      {
        srep(j, 1, blockSize + 1)
        {
          if (f[i][j] == 1 && f[i][j - 1] == 0) {
            sx = i;
            sy = j;
            befx = i + 1;
            befy = j;
          }
        }
      }

      if (sx == -1) {
        assert(false);  // 始点が見つからない場合はエラー
      }

      vector<P> vp;  // ポリゴンの頂点を格納
      vp.emplace_back(sx, sy);
      while (true) {
        int x = -1, y = -1;

        // 境界をたどる（周囲4方向を確認）
        if (x == -1) {
          int nx = sx - 1;
          int ny = sy;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy - 1] != f[sx - 1][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx + 1;
          int ny = sy;
          if (!(nx == befx && ny == befy)) {
            if (f[sx][sy - 1] != f[sx][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx;
          int ny = sy - 1;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy - 1] != f[sx][sy - 1]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == -1) {
          int nx = sx;
          int ny = sy + 1;
          if (!(nx == befx && ny == befy)) {
            if (f[sx - 1][sy] != f[sx][sy]) {
              x = nx;
              y = ny;
            }
          }
        }

        if (x == vp[0].first && y == vp[0].second) break;  // 始点に戻ったらループ終了

        if (x == -1) {
          assert(false);  // 次の点が見つからない場合はエラー
          for (auto p : vp) cout << p.first << ' ' << p.second << endl;
          cout << sx << ' ' << sy << ' ' << befx << ' ' << befy << endl;
          srep(i, 1, blockSize + 1)
          {
            srep(j, 1, blockSize + 1)
            {
              cout << f[i][j];
            }
            cout << endl;
          }
        }

        befx = sx;
        befy = sy;
        sx = x;
        sy = y;
        vp.emplace_back(sx, sy);  // 頂点を追加
      }

      for (auto p : vp) {
        // グリッドのインデックスを実座標に変換
        int x = (p.first - 1) * (100000 / blockSize);
        int y = (p.second - 1) * (100000 / blockSize);
        ans.emplace_back(x, y);  // ポリゴンの頂点として追加
      }

      if (IsLengthOK(ans)) {
        upd = 1;  // 制約を満たしていれば更新
      }
      else {
        upd = 0;  // 制約を満たしていなければ更新しない
        ans = ans2;  // 元の解答に戻す
      }
    }

    if (upd) {
      ansScore = tmpScore;  // スコアを更新
      if (ansScore > best_ansScore) {
        CopyToBest(blockSize);
      }
    }
    else {
      f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];  // 状態を元に戻す
    }
  }

  CopyToAns(blockSize);
}

// 解法のメイン部分を実装する関数
int ff[510][510];
void Method3()
{
  const int blockSize = 20;  // グリッドの分割数（20×20のグリッド）

  InitBlock(blockSize);

  int xx1, yy1, xx2, yy2;  // 最良の矩形領域の座標を格納する変数

  ansScore = 0;  // スコアの初期化
  int loop1 = 0;  // ループ回数のカウンタ
  while (true) {
    if (loop1 % 100 == 0) {
      if (GetNowTime() > TL * 0.1) break;  // 時間制限の半分を過ぎたらループを抜ける
    }
    loop1++;
    // ランダムに矩形領域を選択
    int x1 = Rand() % blockSize;
    int x2 = Rand() % blockSize;
    int y1 = Rand() % blockSize;
    int y2 = Rand() % blockSize;
    if (x1 > x2) swap(x1, x2);  // x1とx2を小さい順に並べ替え
    if (y1 > y2) swap(y1, y2);  // y1とy2を小さい順に並べ替え

    int cnt = 0;  // 選択した領域のスコア
    srep(i, x1, x2 + 1)
    {
      srep(j, y1, y2 + 1)
      {
        cnt += block[i][j];  // 選択したセルのスコアを合計
      }
    }

    if (cnt > ansScore) {
      // スコアが改善された場合、解答を更新
      ansScore = cnt;
      ans.clear();
      // 矩形の四隅の座標を計算し、解答に追加
      ans.emplace_back(x1 * (100000 / blockSize), y1 * (100000 / blockSize));
      ans.emplace_back((x2 + 1) * (100000 / blockSize), y1 * (100000 / blockSize));
      ans.emplace_back((x2 + 1) * (100000 / blockSize), (y2 + 1) * (100000 / blockSize));
      ans.emplace_back(x1 * (100000 / blockSize), (y2 + 1) * (100000 / blockSize));
      // 最良の矩形領域のインデックスを保存
      xx1 = x1;
      yy1 = y1;
      xx2 = x2;
      yy2 = y2;
    }
  }


  rep(i, blockSize + 2)
  {
    rep(j, blockSize + 2)
    {
      f[i][j] = 0;  // 初期化
    }
  }

  // 最良の矩形領域をf配列に設定
  srep(i, xx1, xx2 + 1)
  {
    srep(j, yy1, yy2 + 1)
    {
      f[i + 1][j + 1] = 1;  // 選択した領域を1とする
    }
  }

  CopyToBest(blockSize);

  int loop2 = 0;
  Method3_SA(xx1, xx2, yy1, yy2, blockSize, loop2, TL * 0.75);

  int loop3 = 0;
  int blockSize40 = 40;
  rep(i, blockSize + 2)
  {
    rep(j, blockSize + 2)
    {
      ff[i][j] = f[i][j];
    }
  }
  rep(i, blockSize40 + 2)
  {
    rep(j, blockSize40 + 2)
    {
      f[i][j] = 0;
    }
  }
  rep(i, blockSize + 2)
  {
    rep(j, blockSize + 2)
    {
      if (ff[i][j] == 1) {
        f[i * 2 - 1][j * 2 - 1] = 1;
        f[i * 2 - 1][j * 2] = 1;
        f[i * 2][j * 2 - 1] = 1;
        f[i * 2][j * 2] = 1;
      }
    }
  }
  CopyToBest(blockSize40);

  InitBlock(blockSize40);
  Method3_SA(xx1, xx2, yy1, yy2, blockSize40, loop3, TL * 1.0);

  int loop4 = 0;
  int blockSize80 = 80;
  rep(i, blockSize40 + 2)
  {
    rep(j, blockSize40 + 2)
    {
      ff[i][j] = f[i][j];
    }
  }
  rep(i, blockSize80 + 2)
  {
    rep(j, blockSize80 + 2)
    {
      f[i][j] = 0;
    }
  }
  rep(i, blockSize40 + 2)
  {
    rep(j, blockSize40 + 2)
    {
      if (ff[i][j] == 1) {
        f[i * 2 - 1][j * 2 - 1] = 1;
        f[i * 2 - 1][j * 2] = 1;
        f[i * 2][j * 2 - 1] = 1;
        f[i * 2][j * 2] = 1;
      }
    }
  }
  CopyToBest(blockSize80);

  InitBlock(blockSize80);
  Method3_SA(xx1, xx2, yy1, yy2, blockSize80, loop4, TL);


  if (mode != 0) {
    // デバッグ用の出力
    cout << "loop1 = " << loop1 << ", ";
    cout << "loop2 = " << loop2 << ", ";
    cout << "loop3 = " << loop3 << ", ";
    cout << "loop4 = " << loop4 << ", ";
    cout << endl;
    srep(i, 1, blockSize80 + 1)
    {
      srep(j, 1, blockSize80 + 1)
      {
        cout << f[i][j];
      }
      cout << endl;
    }
  }
}

// 問題を解く関数
ll Solve(int problem_num)
{
  ResetTime();  // 時間計測のリセット

  SetUp();  // 内部状態の初期化

  Input(problem_num);  // 入力の受け取り

  ofstream ofs;
  OpenOfs(problem_num, ofs);  // 出力ファイルストリームのオープン

  Method3();  // 解法の実行

  Output(ofs);  // 解答の出力

  if (ofs.is_open()) {
    ofs.close();  // 出力ファイルストリームのクローズ
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();  // スコアの計算
  }
  return score;  // スコアを返す
}

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main_old()
{
  mode = 2;

  if (mode == 0) {
    Solve(0);  // 問題番号0を解く
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);  // 問題番号iを解く
      sum += score;         // スコアの合計を更新
      if (mode == 1) {
        cout << score << endl;  // スコアを出力
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

  return 0;
}
