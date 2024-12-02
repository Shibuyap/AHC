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
static uint32_t RandXor() {
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
static double Rand01() {
  return (RandXor() + 0.5) * (1.0 / UINT_MAX);
}

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double rot) {
  return l + (rot - l) * Rand01();
}

// 配列をシャッフルする関数（Fisher-Yatesアルゴリズム）
void FisherYates(int* data, int n) {
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

// 非常に大きな値
const int INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 2.8;                                       // 時間制限（Time Limit）
int mode;                                              // 実行モード
std::chrono::steady_clock::time_point startTimeClock;  // 時間計測用

// 時間計測をリセットする関数
void ResetTime() {
  startTimeClock = std::chrono::steady_clock::now();
}

// 現在の経過時間を取得する関数
double GetNowTime() {
  auto endTimeClock = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endTimeClock - startTimeClock;
  return elapsed.count();
}

// 区間の交差判定
bool isCrossing(double l1, double r1, double l2, double r2) {
  return (std::max(l1, l2) < std::min(r1, r2));
}

const int MAX_N = 100;
const int MAX_T = 400;
const int MAX_WIDTH = 100000;
const int MIN_WIDTH = 10000;
const int MAX_HEIGHT = 100000;
const int MIN_HEIGHT = 10000;

int n, t, sigma;
int w[MAX_N], h[MAX_N];

int W[MAX_N], H[MAX_N];
int dW[MAX_T], dH[MAX_T];

class Piece {
public:
  int num;
  int rot;
  int dir;
  int base;

  int width() const {
    if (rot == 1) {
      return h[num];
    }
    return w[num];
  }

  int height() const {
    if (rot == 1) {
      return w[num];
    }
    return h[num];
  }

  // for sorting
  bool operator<(const Piece& other) const {
    return num < other.num;
  }
};

class Block {
public:
  Piece piece1;
  Piece piece2;

  Block() {
    piece1.num = -1;
    piece2.num = -1;
  }

  int count() const {
    if (piece1.num == -1) {
      return 0;
    }
    else if (piece2.num == -1) {
      return 1;
    }
    return 2;
  }

  int width() const {
    int cnt = count();
    if (cnt == 0) {
      return 0;
    }
    else if (cnt == 1) {
      return piece1.width();
    }
    return max(piece1.width(), piece2.width());
  }

  int height() const {
    int cnt = count();
    if (cnt == 0) {
      return 0;
    }
    else if (cnt == 1) {
      return piece1.height();
    }
    return piece1.height() + piece2.height();
  }

  void clear() {
    piece1.num = -1;
    piece2.num = -1;
  }

  void SetBase(int base) {
    piece1.base = base;
    piece2.base = base;
  }
};

class Row {
public:
  Block blocks[MAX_N];
  int sz;
  int sumWidth;
  int maxHeight;

public:
  Row() {
    sz = 0;
    sumWidth = 0;
    maxHeight = 0;
  }

  int count() const {
    return sz;
  }

  void clear() {
    sz = 0;
    sumWidth = 0;
    maxHeight = 0;
  }

  void add(const Block& block) {
    blocks[sz] = block;
    sz++;
    sumWidth += block.width();
    maxHeight = max(maxHeight, block.height());

    int base;
    if (sz == 1) {
      base = -1;
    }
    else {
      base = blocks[sz - 2].piece1.num;
      if (blocks[sz - 2].count() == 2 && blocks[sz - 2].piece2.width() > blocks[sz - 2].piece1.width()) {
        base = blocks[sz - 2].piece2.num;
      }
    }
    blocks[sz - 1].SetBase(base);
  }

  void addPiece(int index, const Piece& piece) {
    int beforeWidth = blocks[index].width();
    blocks[index].piece2 = piece;
    blocks[index].piece2.base = blocks[index].piece1.base;
    sumWidth += blocks[index].width() - beforeWidth;
    maxHeight = max(maxHeight, blocks[index].height());
  }

  int GetSumWidth() const {
    return sumWidth;
  }

  int GetMaxHeight() const {
    return maxHeight;
  }
};

class Board {
private:
  Row rows[MAX_N];
  int sz;
  int maxWidth;
  int sumHeight;

public:
  Board() {
    sz = 0;
    maxWidth = 0;
    sumHeight = 0;
  }

  void clear() {
    sz = 0;
    maxWidth = 0;
    sumHeight = 0;
  }

  int count() const {
    return sz;
  }

  // ※厳密なスコアではない
  int score() const {
    return maxWidth + sumHeight;
  }

  void Add(const Row& row) {
    rows[sz] = row;
    sz++;
    maxWidth = max(maxWidth, row.GetSumWidth());
    sumHeight += row.GetMaxHeight();
  }

  void Add(int index, const Block& block) {
    int beforeMaxHeight = rows[index].GetMaxHeight();
    rows[index].add(block);
    maxWidth = max(maxWidth, rows[index].GetSumWidth());
    sumHeight += rows[index].GetMaxHeight() - beforeMaxHeight;
  }

  int GetMaxWidth() const {
    return maxWidth;
  }

  int GetSumHeight() const {
    return sumHeight;
  }

  vector<Piece> CreateQuery() const {
    vector<Piece> pieces;
    rep(i, sz) {
      rep(j, rows[i].sz) {
        if (rows[i].blocks[j].count() >= 1) {
          pieces.push_back(rows[i].blocks[j].piece1);
        }
        if (rows[i].blocks[j].count() >= 2) {
          pieces.push_back(rows[i].blocks[j].piece2);
        }
      }
    }
    sort(pieces.begin(), pieces.end());
    return pieces;
  }

  Row& back() {
    return rows[sz - 1];
  }
};

struct Score {
  int score;
  int ww;
  int hh;
};

int queryCount;
Score tScores[MAX_T];

void CopyToBest() {
}

void CopyToAns() {
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp() {
  queryCount = 0;
}

// 入力を受け取る関数
void Input(int problemNum) {
  if (mode == 0) {
    // 標準入力
    cin >> n >> t >> sigma;
    rep(i, n) {
      cin >> w[i] >> h[i];
    }
  }
  else {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
    ifstream ifs(oss.str());
    // ファイル入力
    ifs >> n >> t >> sigma;
    rep(i, n) {
      ifs >> w[i] >> h[i];
    }
    rep(i, n) {
      ifs >> W[i] >> H[i];
    }
  }

  rep(i, n) {
    w[i] = max(MIN_WIDTH, w[i]);
    w[i] = min(MAX_WIDTH, w[i]);
    h[i] = max(MIN_HEIGHT, h[i]);
    h[i] = min(MAX_HEIGHT, h[i]);
  }
}

// 出力ファイルストリームを開く関数
void OpenOfs(int probNum, ofstream& ofs) {
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

Score CalcTScore() {
  Score score;
  score.score = INF;
  rep(i, queryCount) {
    if (tScores[i].score < score.score) {
      score = tScores[i];
    }
  }
  return score;
}

int cs_use[MAX_N] = {};
int cs_up[MAX_N], cs_down[MAX_N], cs_left[MAX_N], cs_right[MAX_N];
// スコアを計算する関数
Score CalcScore(const vector<Piece>& pieces, bool cheat) {
  rep(i, MAX_N) {
    cs_up[i] = -1;
    cs_down[i] = -1;
    cs_left[i] = -1;
    cs_right[i] = -1;
  }
  for (auto col : pieces) {
    int num = col.num;
    int wid = w[num];
    int hei = h[num];
    if (mode != 0 && cheat == true) {
      wid = W[num];
      hei = H[num];
    }
    if (col.rot == 1) swap(wid, hei);
    if (col.dir == 0) {
      cs_left[num] = 0;
      if (col.base != -1) {
        cs_left[num] = cs_right[col.base];
      }
      cs_right[num] = cs_left[num] + wid;

      cs_up[num] = 0;
      rep(i, num) {
        if (cs_use[i]) {
          if (isCrossing(cs_left[num], cs_right[num], cs_left[i], cs_right[i])) {
            cs_up[num] = max(cs_up[num], cs_down[i]);
          }
        }
      }
      cs_down[num] = cs_up[num] + hei;
    }
    else {
      cs_up[num] = 0;
      if (col.base != -1) {
        cs_up[num] = cs_down[col.base];
      }
      cs_down[num] = cs_up[num] + hei;

      cs_left[num] = 0;
      rep(i, num) {
        if (cs_use[i]) {
          if (isCrossing(cs_up[num], cs_down[num], cs_up[i], cs_down[i])) {
            cs_left[num] = max(cs_left[num], cs_right[i]);
          }
        }
      }
      cs_right[num] = cs_left[num] + wid;
    }

    cs_use[num] = 1;
  }

  int maxDown = 0, maxRight = 0;
  rep(i, n) {
    if (cs_use[i]) {
      maxDown = max(maxDown, cs_down[i]);
      maxRight = max(maxRight, cs_right[i]);
    }
  }

  Score score;

  score.ww = maxRight;
  score.hh = maxDown;
  score.score = maxDown + maxRight;
  rep(i, n) {
    if (cs_use[i] == 0) {
      score.score += w[i] + h[i];
    }
  }

  return score;
}

Score Print(const vector<Piece>& pieces, ofstream& ofs) {
  Score score;

  if (mode == 0) {
    cout << pieces.size() << endl;
    rep(i, pieces.size()) {
      cout << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }
    fflush(stdout);

    cin >> score.ww >> score.hh;
    score.score = score.ww + score.hh;
  }
  else {
    ofs << pieces.size() << endl;
    rep(i, pieces.size()) {
      ofs << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }

    score = CalcScore(pieces, true);
    score.ww += dW[queryCount];
    score.hh += dH[queryCount];
    score.score = score.ww + score.hh;
  }

  tScores[queryCount] = score;
  queryCount++;

  return score;
}

int bestsCount;
vector<Piece> bests[MAX_T];
int bestScores[MAX_T];

void Method1_Shoki1_Internal1(vector<Piece>& tmp) {
  int widSum = 0;
  int heiSum = 0;

  int widLimit = RandXor() % 1000000 + 200000;
  int now = 0;
  int maxHeight = 0;
  rep(i, n) {
    int wid = w[i];
    int hei = h[i];
    if (tmp[i].rot == 1) {
      wid = h[i];
      hei = w[i];
    }

    if (now + wid <= widLimit) {
      tmp[i].base = i - 1;
      now += wid;
      maxHeight = max(maxHeight, hei);
    }
    else {
      tmp[i].base = -1;
      widSum = max(widSum, now);
      now = wid;
      heiSum += maxHeight;
      maxHeight = hei;
    }
  }

  widSum = max(widSum, now);
  heiSum += maxHeight;

  {
    int tmpScore = widSum + heiSum;
    if (bestsCount < t / 2) {
      bests[bestsCount] = tmp;
      bestScores[bestsCount] = tmpScore;
      bestsCount++;
    }
    else if (tmpScore < bestScores[bestsCount - 1]) {
      bests[bestsCount - 1] = tmp;
      bestScores[bestsCount - 1] = tmpScore;
    }
    int now = bestsCount - 1;
    while (now >= 1) {
      if (bestScores[now] < bestScores[now - 1]) {
        swap(bests[now], bests[now - 1]);
        swap(bestScores[now], bestScores[now - 1]);
        now--;
      }
      else {
        break;
      }
    }
  }
}

int keepRot[MAX_N];
int keepRotCount = 0;
void Method1_Shoki1_Internal2(vector<Piece>& tmp) {
  int widSum = 0;
  int heiSum = 0;

  int widLimit = RandXor() % 1000000 + 200000;

  int now = 0;
  int last = -1;
  int maxHeight = 0;

  int beforeNow = INF;
  int beforeLast = -1;
  int beforeMaxHeight = 0;

  int isRandomRot = RandXor() % 10;

  keepRotCount = 0;
  rep(i, n) {
    if (isRandomRot >= 1 && RandXor() % n <= isRandomRot) {
      tmp[i].rot = 1 - tmp[i].rot;
      keepRot[keepRotCount] = i;
      keepRotCount++;
    }

    int wid = w[i];
    int hei = h[i];
    if (tmp[i].rot == 1) {
      wid = h[i];
      hei = w[i];
    }

    if (now < widLimit * 0.9 && beforeNow + wid <= widLimit) {
      tmp[i].base = beforeLast;
      beforeLast = i;
      beforeNow += wid;
      widSum = max(widSum, beforeNow);
      if (hei > beforeMaxHeight) {
        heiSum += hei - beforeMaxHeight;
        beforeMaxHeight = hei;
      }
    }
    else if (now + wid <= widLimit) {
      tmp[i].base = last;
      last = i;
      now += wid;
      maxHeight = max(maxHeight, hei);
    }
    else {
      beforeLast = last;
      beforeMaxHeight = maxHeight;
      beforeNow = now;

      tmp[i].base = -1;
      last = i;
      widSum = max(widSum, now);
      now = wid;
      heiSum += maxHeight;
      maxHeight = hei;
    }
  }

  widSum = max(widSum, now);
  heiSum += maxHeight;

  {
    int tmpScore = widSum + heiSum;
    if (bestsCount < t / 2) {
      bests[bestsCount] = tmp;
      bestScores[bestsCount] = tmpScore;
      bestsCount++;
    }
    else if (tmpScore < bestScores[bestsCount - 1]) {
      bests[bestsCount - 1] = tmp;
      bestScores[bestsCount - 1] = tmpScore;
    }
    int now = bestsCount - 1;
    while (now >= 1) {
      if (bestScores[now] < bestScores[now - 1]) {
        swap(bests[now], bests[now - 1]);
        swap(bestScores[now], bestScores[now - 1]);
        now--;
      }
      else {
        break;
      }
    }
  }

  rep(i, keepRotCount) {
    int ii = keepRot[i];
    tmp[ii].rot = 1 - tmp[ii].rot;
  }
}

void Method1_Shoki1() {
  bestsCount = 0;

  vector<Piece> tmp(n);
  rep(i, n) {
    tmp[i].num = i;
    tmp[i].rot = 0;
    if (w[i] > h[i]) tmp[i].rot = 1;
    tmp[i].dir = 0;
    tmp[i].base = i - 1;
  }

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30) {
        break;
      }
    }

    loop++;

    Method1_Shoki1_Internal1(tmp);
  }

  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30 * 25) {
        break;
      }
    }

    loop++;

    Method1_Shoki1_Internal2(tmp);
  }
}

void Method1(ofstream& ofs) {
  vector<Piece> best(n);
  int bestScore = INF;

  vector<Piece> tmp(n);

  Method1_Shoki1();

  rep(aespa, bestsCount) {
    tmp = bests[aespa];
    Score score = Print(tmp, ofs);
    if (score.score < bestScore) {
      bestScore = score.score;
      best = tmp;
    }
  }

  srep(aespa, bestsCount, t) {
    int raMode = RandXor() % 2;

    tmp = best;
    if (raMode == 0) {
      int ra = RandXor() % n;
      tmp[ra].rot = 1 - tmp[ra].rot;
    }
    else {
      vector<vector<Piece>> vvc;
      vector<Piece> vc;
      for (auto col : best) {
        if (col.base == -1 && !vc.empty()) {
          vvc.push_back(vc);
          vc.clear();
        }
        vc.push_back(col);
      }
      vvc.push_back(vc);
      if (vvc.size() >= 2) {
        while (true) {
          int ra1 = RandXor() % (vvc.size() - 1);
          int ra2 = RandXor() % 2;
          if (ra2 == 0) {
            if (vvc[ra1].size() >= 2) {
              vvc[ra1 + 1].insert(vvc[ra1 + 1].begin(), vvc[ra1].back());
              vvc[ra1].pop_back();
              break;
            }
          }
          else {
            if (vvc[ra1 + 1].size() >= 2) {
              vvc[ra1].push_back(vvc[ra1 + 1][0]);
              vvc[ra1 + 1].erase(vvc[ra1 + 1].begin());
              break;
            }
          }
        }
      }
      tmp.clear();
      rep(i, vvc.size()) {
        rep(j, vvc[i].size()) {
          Piece col = vvc[i][j];
          if (j == 0) {
            col.base = -1;
          }
          else {
            col.base = vvc[i][j - 1].num;
          }
          tmp.push_back(col);
        }
      }
    }

    Score score = Print(tmp, ofs);
    if (score.score < bestScore) {
      bestScore = score.score;
      best = tmp;
    }
  }
}

vector<Piece> initialPieces;
void InitializePieces() {
  initialPieces.resize(n);
  rep(i, n) {
    initialPieces[i].num = i;
    initialPieces[i].rot = 0;
    if (w[i] > h[i]) initialPieces[i].rot = 1;
    initialPieces[i].dir = 0;
    initialPieces[i].base = i - 1;
  }
}

int sizeRank[MAX_N];
void InitializeSizeRank() {
  vector<P> vp;
  rep(i, n) {
    vp.emplace_back(max(w[i], h[i]), i);
  }
  sort(vp.begin(), vp.end());
  rep(i, n) {
    sizeRank[vp[i].second] = i;
  }
}

Block block;
Row row;
Board board;

int bestsCount2;
Board bests2[MAX_T];
void Method2_Shoki1_Internal1() {
  block.clear();
  row.clear();
  board.clear();

  int widLimit = RandXor() % 1000000 + 200000;

  rep(i, n) {
    block.piece1 = initialPieces[i];
    if (row.GetSumWidth() + block.piece1.width() <= widLimit) {
      row.add(block);
    }
    else {
      board.Add(row);
      row.clear();
      row.add(block);
    }
  }

  board.Add(row);

  // bests更新
  if (bestsCount2 < t / 2) {
    bests2[bestsCount2] = board;
    bestsCount2++;
  }
  else if (board.score() < bests2[bestsCount2 - 1].score()) {
    bests2[bestsCount2 - 1] = board;
  }
  int now = bestsCount2 - 1;
  while (now >= 1) {
    if (bests2[now].score() < bests2[now - 1].score()) {
      swap(bests2[now], bests2[now - 1]);
      now--;
    }
    else {
      break;
    }
  }
}

void Method2_Shoki1_Internal2() {
  block.clear();
  row.clear();
  board.clear();

  int widLimit = RandXor() % 1000000 + 200000;

  int isRandomRot = RandXor() % 10;

  rep(i, n) {
    block.piece1 = initialPieces[i];

    if (isRandomRot >= 1 && RandXor() % n <= isRandomRot) {
      block.piece1.rot = 1 - block.piece1.rot;
    }

    if (row.GetSumWidth() < widLimit * 0.9 && board.count() >= 1 && board.back().GetSumWidth() + block.piece1.width() <= widLimit) {
      board.Add(board.count() - 1, block);
    }
    else if (row.GetSumWidth() + block.piece1.width() <= widLimit) {
      row.add(block);
    }
    else {
      board.Add(row);
      row.clear();
      row.add(block);
    }
  }

  board.Add(row);

  // bests更新
  if (bestsCount2 < t / 2) {
    bests2[bestsCount2] = board;
    bestsCount2++;
  }
  else if (board.score() < bests2[bestsCount2 - 1].score()) {
    bests2[bestsCount2 - 1] = board;
  }
  int now = bestsCount2 - 1;
  while (now >= 1) {
    if (bests2[now].score() < bests2[now - 1].score()) {
      swap(bests2[now], bests2[now - 1]);
      now--;
    }
    else {
      break;
    }
  }
}

void Method2_Shoki1_Internal3() {
  block.clear();
  row.clear();
  board.clear();
  Piece piece;

  int widLimit = RandXor() % 1000000 + 200000;

  int isRandomRot = RandXor() % 10;

  int isUsePiece2_1 = RandXor() % 200;

  int bestsLimit = t / 2;
  int ng = 0;

  rep(i, n) {
    if (bestsCount2 == bestsLimit && board.score() + row.maxHeight >= bests2[bestsCount2 - 1].score()) {
      ng = 1;
      break;
    }

    piece = initialPieces[i];

    if (isRandomRot >= 1 && RandXor() % n <= isRandomRot) {
      piece.rot = 1 - piece.rot;
    }

    if (isUsePiece2_1 < 50 && sizeRank[i] < 10) {
      if (row.sz > 0 && row.blocks[row.sz - 1].count() == 1 && row.blocks[row.sz - 1].height() + piece.height() < MAX_HEIGHT * 1.1) {
        int ok = 1;
        if (RandXor() % 25 > isUsePiece2_1) {
          ok = 0;
        }
        if (ok) {
          row.addPiece(row.sz - 1, piece);
          continue;
        }
      }
    }
    else if (isUsePiece2_1 < 100 && sizeRank[i] < 10) {
      if (row.sz >= 1 && row.blocks[row.sz - 1].count() == 1 && row.blocks[row.sz - 1].height() + piece.height() < MAX_HEIGHT * 1.0) {
        int ok = 1;
        if (RandXor() % 2 == 0) {
          ok = 0;
        }
        if (ok) {
          row.addPiece(row.sz - 1, piece);
          continue;
        }
      }

      if (row.sz >= 2) {
        int isAdd = 0;
        rep(j, row.sz - 1) {
          if (row.blocks[j].count() == 1 && row.blocks[j].height() + piece.height() < MAX_HEIGHT * 1.0 && piece.width() < row.blocks[j].width()) {
            int ok = 1;
            if (RandXor() % 2 == 0) {
              ok = 0;
            }
            if (ok) {
              row.addPiece(j, piece);
              isAdd = 1;
              break;
            }
          }
        }
        if (isAdd) {
          continue;
        }
      }
    }

    block.piece1 = piece;
    if (row.GetSumWidth() < widLimit * 0.9 && board.count() >= 1 && board.back().GetSumWidth() + block.piece1.width() <= widLimit) {
      board.Add(board.count() - 1, block);
    }
    else if (row.GetSumWidth() + block.piece1.width() <= widLimit) {
      row.add(block);
    }
    else {
      board.Add(row);
      row.clear();
      row.add(block);
    }
  }

  if (ng) {
    return;
  }

  board.Add(row);

  // bests更新
  if (bestsCount2 < t / 2) {
    bests2[bestsCount2] = board;
    bestsCount2++;
  }
  else if (board.score() < bests2[bestsCount2 - 1].score()) {
    bests2[bestsCount2 - 1] = board;
  }
  int now = bestsCount2 - 1;
  while (now >= 1) {
    if (bests2[now].score() < bests2[now - 1].score()) {
      swap(bests2[now], bests2[now - 1]);
      now--;
    }
    else {
      break;
    }
  }
}


void Method2_Shoki1() {
  bestsCount2 = 0;

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30) {
        break;
      }
    }

    loop++;

    Method2_Shoki1_Internal1();
  }

  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > TL / 30 * 25) {
        break;
      }
    }

    loop++;

    //Method2_Shoki1_Internal2();
    Method2_Shoki1_Internal3();
  }

  if (mode != 0) {
    cout << "loop = " << loop << endl;
  }
}

void Method2(ofstream& ofs) {
  Board best;
  int bestScore = INF;

  Method2_Shoki1();

  int aespa = 0;
  while (aespa < bestsCount2) {
    Score score = Print(bests2[aespa].CreateQuery(), ofs);
    if (score.score < bestScore) {
      bestScore = score.score;
      best = bests2[aespa];
    }
    aespa++;
  }

  vector<Piece> bestPieces = best.CreateQuery();
  vector<Piece> tmp = bestPieces;

  int loop = 0;
  double startTime = GetNowTime();
  double nowTime = GetNowTime();  // 現在の経過時間
  const double START_TEMP = 0.0;  // 焼きなまし法の開始温度
  const double END_TEMP = 200000.0;      // 焼きなまし法の終了温度
  const double timeLimit = TL - startTime;
  double temp = START_TEMP;  // 現在の温度
  while (aespa < t) {
    loop++;
    if (loop % 100 == 0) {
      nowTime = GetNowTime();
    }

    int raMode = RandXor() % 3;

    if (RandXor() % 2 == 0) {
      tmp = bestPieces;
    }

    if (raMode == 0) {
      int ra = RandXor() % n;
      tmp[ra].rot = 1 - tmp[ra].rot;
    }
    else if (raMode == 1) {
      int ra = RandXor() % n;
      int ra2 = RandXor() % (ra + 1) - 1;
      if (tmp[ra].base == ra2)continue;
      tmp[ra].base = ra2;
    }
    else {
      int ra = RandXor() % n;
      tmp[ra].dir = 1 - tmp[ra].dir;
      if (RandXor() % 2 == 0) {
        int ra2 = RandXor() % (ra + 1) - 1;
        tmp[ra].base = ra2;
      }
    }

    auto preScore = CalcScore(tmp, false);

    // だんだん温度が上がる焼きなまし
    const double progressRatio = (nowTime - startTime) / timeLimit;  // 進捗率（0.0～1.0）
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio * progressRatio * progressRatio;  // 温度の更新
    double diff = bestScore - preScore.score - 100000;  // スコアの差分
    double prob = exp(diff / temp);     // 焼きなまし法の採用確率
    if (preScore.score < bestScore || prob > Rand01() || nowTime > TL) {
      Score score = Print(tmp, ofs);
      if (score.score < bestScore) {
        if (mode != 0) {
          cout << raMode << ' ' << diff << ' ' << temp << ' ' << prob << ' ' << nowTime << endl;
        }
        bestScore = score.score;
        bestPieces = tmp;
      }
      aespa++;
    }
  }

  if (mode != 0) {
    cout << "aespa loop = " << loop << endl;
  }
}

// 問題を解く関数
ll Solve(int problem_num) {
  ResetTime();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(problem_num, ofs);

  // 初期解生成
  InitializePieces();
  InitializeSizeRank();
  Method2(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  Score score;
  score.score = 0;
  if (mode != 0) {
    score = CalcTScore();
    cout << "ww = " << score.ww << ", hh = " << score.hh << ", score = " << score.score << endl;
  }
  return score.score;
}

/////////////////////////////////////////////////////////////////////////
/*
メモ

*/
/////////////////////////////////////////////////////////////////////////
int main() {
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  mode = 2;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 3) {
    rep(_, 10) {
      ll sum = 0;
      srep(i, 0, 100) {
        ll score = Solve(i);
        sum += score;
      }
      cout << sum << endl;
    }
  }
  else {
    ll sum = 0;
    srep(i, 0, 100) {
      ll score = Solve(i);
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

  return 0;
}
