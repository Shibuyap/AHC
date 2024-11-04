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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

// 型定義のエイリアス
typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// 乱数生成（XorShift法による擬似乱数生成器）
static uint32_t randxor()
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

// 0以上1未満の小数を返す乱数関数
static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = randxor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化（使用されていない）
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;   // 非常に大きな値（オーバーフローに注意）
const int INT_INF = 1001001001;       // int型の非常に大きな値

// 移動方向の配列（使用されていない）
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;  // 時間制限（Time Limit）
int mode;         // 実行モード
std::chrono::steady_clock::time_point startTime, endTime;  // 時間計測用

// 時間計測をリセットする関数
void ResetTime()
{
  startTime = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime()
{
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTime - startTime;
  return elapsed.count();
}

// 箱の位置を表す構造体
struct Point {
  int x;  // 山の番号
  int y;  // 山の中での高さ（下からの位置）
};

const int MAX_N = 30;  // 使用されていない

const int n = 200;  // 箱の総数
const int m = 10;   // 山の総数

vector<int> init_b[m];   // 初期状態の各山に積まれた箱の番号リスト
vector<Point> init_c;    // 初期状態の各箱の位置情報

// 問題の状態を表す構造体
struct Problem {
  vector<int> b[m];      // 現在の各山の状態
  vector<Point> c;       // 現在の各箱の位置
  vector<P> ans;         // 操作列の記録（解答）
};

// 不使用の関数（おそらく将来的な最良解の保存用）
void CopyToBest()
{
}

// 不使用の関数（おそらく現在の解答をコピーするためのもの）
void CopyToAns()
{
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  rep(i, m) {
    init_b[i].clear();
  }
  init_c.clear();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  rep(i, m) {
    init_b[i].resize(n / m);  // 各山のサイズを設定
  }
  init_c.resize(n);  // 箱の位置情報のサイズを設定

  // 標準入力から受け取る場合
  if (!ifs.is_open()) {
    int _n, _m;
    cin >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        cin >> init_b[i][j];
        init_b[i][j]--;  // 0-indexedに変換
      }
    }
  }
  // ファイルから受け取る場合
  else {
    int _n, _m;
    ifs >> _n >> _m;
    rep(i, m)
    {
      rep(j, n / m)
      {
        ifs >> init_b[i][j];
        init_b[i][j]--;  // 0-indexedに変換
      }
    }
  }

  // 各箱の位置情報を設定
  rep(i, m) {
    rep(j, n / m) {
      init_c[init_b[i][j]].x = i;  // 山の番号
      init_c[init_b[i][j]].y = j;  // 高さ
    }
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
int CalcScore(const vector<P>& ans)
{
  vector<int> tmp_b[m];          // 一時的な山の状態
  vector<Point> tmp_c = init_c;  // 一時的な箱の位置情報

  int cnt = 0;  // 運び出した箱の数
  rep(i, m)
  {
    tmp_b[i] = init_b[i];  // 初期状態をコピー
  }

  int res = 10000;  // 初期スコア
  rep(i, ans.size())
  {
    int num = ans[i].first;  // 操作する箱の番号
    int x = tmp_c[num].x;    // 箱の現在の山の番号
    int y = tmp_c[num].y;    // 箱の現在の高さ
    int nx = ans[i].second;  // 移動先の山の番号
    if (nx == -1) {
      // 操作2：箱を運び出す
      tmp_b[x].pop_back();
      cnt++;
    }
    else {
      // 操作1：箱とその上の箱を移動
      res--;
      rep(j, tmp_b[x].size() - y)
      {
        int num2 = tmp_b[x][y + j];
        tmp_c[num2].x = nx;           // 箱の新しい山の番号を設定
        tmp_c[num2].y = tmp_b[nx].size();  // 新しい高さを設定
        tmp_b[nx].push_back(num2);    // 移動先の山に箱を追加
        res--;
      }
      tmp_b[x].resize(y);  // 元の山から移動した分を削除
    }
  }
  if (cnt != n) return -1;  // 全ての箱を運び出せなかった場合
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs, const vector<P>& ans)
{
  if (mode == 0) {
    // 標準出力に出力
    rep(i, ans.size()) { cout << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
  else {
    // ファイルに出力
    rep(i, ans.size()) { ofs << ans[i].first + 1 << ' ' << ans[i].second + 1 << endl; }
  }
}

// 1ターンでの操作を実行する関数
void MoveOneTurn(Problem& prob, int turn, int from, const vector<int>& idx) {
  int turnY = prob.c[turn].y;  // 現在のターンで運び出す箱の高さ

  // 運び出す箱の上にある箱を移動
  for (int k = prob.b[from].size() - 1; k > turnY;) {
    int to = idx[k];  // 移動先の山
    while (idx[k - 1] == to) k--;  // 同じ移動先の箱をまとめる

    // 箱を移動
    srep(l, k, prob.b[from].size()) {
      int num = prob.b[from][l];
      prob.c[num].x = to;           // 新しい山の番号
      prob.c[num].y = prob.b[to].size();  // 新しい高さ
      prob.b[to].push_back(num);    // 移動先の山に追加
    }

    // 操作を記録
    prob.ans.emplace_back(prob.b[from][k], to);
    prob.b[from].resize(k);  // 元の山から削除

    k--;
  }

  // 現在のターンの箱を運び出す
  prob.b[from].pop_back();
  prob.ans.emplace_back(turn, -1);  // 操作を記録
}

// 箱の移動先を決定する関数
void DecideIdx(const Problem& prob, vector<int>& idx, const vector<P>& minIs, int from, int k) {
  int bb = prob.b[from][k];  // 移動する箱の番号
  int id = 1;
  while (id + 1 < m && minIs[id].first < bb) id++;
  idx[k] = minIs[id].second;  // 移動先の山の番号を設定
}

// 残りの操作をシミュレーションする関数（プレイアウト）
int PlayOut(Problem prob, int _turn, int _k, vector<int> _idx, vector<P> _minIs, int _from) {
  int turnY = prob.c[_turn].y;
  // 移動先を決定
  srep(k, _k + 1, prob.b[_from].size()) {
    DecideIdx(prob, _idx, _minIs, _from, k);
  }
  // 1ターン分の操作を実行
  MoveOneTurn(prob, _turn, _from, _idx);

  // 残りのターンを順次実行
  srep(turn, _turn + 1, n) {
    vector<P> minIs(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, prob.b[i].size()) {
        minI = min(minI, prob.b[i][j]);  // 各山の最小の箱の番号を取得
      }
      minIs[i] = P(minI, i);
    }

    sort(minIs.begin(), minIs.end());  // 最小の箱の番号でソート

    int from = minIs[0].second;  // 次に操作する山
    int turnY = prob.c[turn].y;

    vector<int> idx(prob.b[from].size(), -1);
    // 箱の移動先を決定
    srep(k, turnY + 1, prob.b[from].size()) {
      DecideIdx(prob, idx, minIs, from, k);
    }

    // 1ターン分の操作を実行
    MoveOneTurn(prob, turn, from, idx);
  }

  return CalcScore(prob.ans);  // スコアを計算して返す
}

// kotatsugameさんの貪欲法を実装した関数
Problem Method2() {
  Problem prob;
  rep(i, m) prob.b[i] = init_b[i];  // 初期状態をコピー
  prob.c = init_c;

  rep(turn, n) {  // 各ターン（各箱）について
    vector<P> minIs(m);
    rep(i, m) {
      int minI = INT_INF;
      rep(j, prob.b[i].size()) {
        minI = min(minI, prob.b[i][j]);  // 各山の最小の箱の番号を取得
      }
      minIs[i] = P(minI, i);
    }

    sort(minIs.begin(), minIs.end());  // 最小の箱の番号でソート

    int from = minIs[0].second;  // 現在のターンで操作する山
    int turnY = prob.c[turn].y;

    vector<int> keepIdx(prob.b[from].size(), -1);
    // 箱の移動先を決定
    srep(k, turnY + 1, prob.b[from].size()) {
      DecideIdx(prob, keepIdx, minIs, from, k);
    }

    vector<int> idx(prob.b[from].size(), -1);
    // 各箱について最適な移動先を探索
    srep(k, turnY + 1, prob.b[from].size()) {
      int maxScore = -1;
      int maxId = -1;

      // 各可能な移動先についてプレイアウト
      srep(l, 1, m) {
        idx[k] = minIs[l].second;
        int tmpScore = PlayOut(prob, turn, k, idx, minIs, from);
        if (tmpScore > maxScore) {
          maxScore = tmpScore;
          maxId = minIs[l].second;
        }
      }

      idx[k] = maxId;  // 最良の移動先を設定
    }

    // 現在のスコアを計算
    int score = PlayOut(prob, turn, prob.b[from].size() - 1, idx, minIs, from);

    // ランダムに移動先を変更して探索（焼きなまし的な手法）
    if (prob.b[from].size() - (turnY + 1) >= 2 && GetNowTime() < TL) {
      rep(aespa, 500) {
        int randomK = randxor() % (prob.b[from].size() - (turnY + 1)) + turnY + 1;
        int keep = idx[randomK];

        if (randxor() % 2 == 0) {
          // ランダムに移動先を変更
          int randomIdx = randxor() % (m - 1) + 1;
          idx[randomK] = minIs[randomIdx].second;
        }
        else {
          // 上下の箱の移動先に合わせる
          int randomDir = -1;
          if (randomK == turnY + 1) {
            randomDir = 1;
          }
          else if (randomK != turnY + 1 && randomK != prob.b[from].size() - 1) {
            if (randxor() % 2 == 0) {
              randomDir = 1;
            }
          }
          idx[randomK] = idx[randomK + randomDir];
        }

        int tmpScore = PlayOut(prob, turn, prob.b[from].size() - 1, idx, minIs, from);
        if (tmpScore >= score) {
          score = tmpScore;
        }
        else {
          idx[randomK] = keep;  // 改善しなければ元に戻す
        }
      }
    }

    // 決定した移動先で操作を実行
    MoveOneTurn(prob, turn, from, idx);
  }

  return prob;
}

// 問題を解く関数
ll Solve(int probNum)
{
  ResetTime();  // 時間計測をリセット

  // 内部状態を初期化
  SetUp();

  // 入力を受け取る
  Input(probNum);

  // 出力ファイルストリームを開く
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解を生成
  auto problem = Method2();

  // 解答を出力
  Output(ofs, problem.ans);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore(problem.ans);  // スコアを計算
  }
  return score;
}

/////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////
int main()
{
  srand((unsigned)time(NULL));  // 乱数の種を設定
  while (rand() % 100) {
    randxor();  // 乱数を進めておく
  }

  mode = 2;  // 実行モードの設定

  if (mode == 0) {
    Solve(0);  // 単一のケースを解く
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);  // 複数のケースを解く
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
