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
static double Rand01()
{
  return (RandXor() + 0.5) * (1.0 / UINT_MAX);
}

// l以上r未満の実数をとる乱数
static double RandUniform(double l, double rot)
{
  return l + (rot - l) * Rand01();
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

// 非常に大きな値
const int INF = 1001001001;

// 移動方向の配列
const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 2.8;                                       // 時間制限（Time Limit）
int mode;                                              // 実行モード
std::chrono::steady_clock::time_point startTimeClock;  // 時間計測用

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

// 区間の交差判定
bool isCrossing(double l1, double r1, double l2, double r2)
{
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

class Piece
{
public:
  int num;
  int rot;
  int dir;
  int base;

  int width() const
  {
    if (rot == 1) {
      return h[num];
    }
    return w[num];
  }

  int height() const
  {
    if (rot == 1) {
      return w[num];
    }
    return h[num];
  }

  // for sorting
  bool operator<(const Piece& other) const
  {
    return num < other.num;
  }
};

class Block
{
public:
  Piece piece1;
  Piece piece2;

  Block()
  {
    piece1.num = -1;
    piece2.num = -1;
  }

  int count() const
  {
    if (piece1.num == -1) {
      return 0;
    }
    else if (piece2.num == -1) {
      return 1;
    }
    return 2;
  }

  int width() const
  {
    int cnt = count();
    if (cnt == 0) {
      return 0;
    }
    else if (cnt == 1) {
      return piece1.width();
    }
    return max(piece1.width(), piece2.width());
  }

  int height() const
  {
    int cnt = count();
    if (cnt == 0) {
      return 0;
    }
    else if (cnt == 1) {
      return piece1.height();
    }
    return piece1.height() + piece2.height();
  }

  void clear()
  {
    piece1.num = -1;
    piece2.num = -1;
  }

  void SetBase(int base)
  {
    piece1.base = base;
    piece2.base = base;
  }
};

class Row
{
public:
  Block blocks[MAX_N];
  int sz;
  int sumWidth;
  int maxHeight;

public:
  Row()
  {
    sz = 0;
    sumWidth = 0;
    maxHeight = 0;
  }

  int count() const
  {
    return sz;
  }

  void clear()
  {
    sz = 0;
    sumWidth = 0;
    maxHeight = 0;
  }

  void add(const Block& block)
  {
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

  void rotateBack()
  {
    sumWidth -= blocks[sz - 1].width();
    sumWidth += blocks[sz - 1].height();
    blocks[sz - 1].piece1.rot = 1 - blocks[sz - 1].piece1.rot;
  }

  void addPiece(int index, const Piece& piece)
  {
    int beforeWidth = blocks[index].width();
    blocks[index].piece2 = piece;
    blocks[index].piece2.base = blocks[index].piece1.base;
    sumWidth += blocks[index].width() - beforeWidth;
    maxHeight = max(maxHeight, blocks[index].height());
  }

  int GetSumWidth() const
  {
    return sumWidth;
  }

  int GetMaxHeight() const
  {
    return maxHeight;
  }
};

class Board
{
public:
  Row rows[MAX_N];
  int sz;
  int maxWidth;
  int sumHeight;

public:
  Board()
  {
    sz = 0;
    maxWidth = 0;
    sumHeight = 0;
  }

  void clear()
  {
    sz = 0;
    maxWidth = 0;
    sumHeight = 0;
  }

  int count() const
  {
    return sz;
  }

  // ※厳密なスコアではない
  int score() const
  {
    return maxWidth + sumHeight;
  }

  void Add(const Row& row)
  {
    rows[sz] = row;
    sz++;
    maxWidth = max(maxWidth, row.GetSumWidth());
    sumHeight += row.GetMaxHeight();
  }

  void Add(int index, const Block& block)
  {
    int beforeMaxHeight = rows[index].GetMaxHeight();
    rows[index].add(block);
    maxWidth = max(maxWidth, rows[index].GetSumWidth());
    sumHeight += rows[index].GetMaxHeight() - beforeMaxHeight;
  }

  int GetMaxWidth() const
  {
    return maxWidth;
  }

  int GetSumHeight() const
  {
    return sumHeight;
  }

  vector<Piece> CreateQuery() const
  {
    vector<Piece> pieces;
    rep(i, sz)
    {
      rep(j, rows[i].sz)
      {
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

  Row& back()
  {
    return rows[sz - 1];
  }
};

struct Score
{
  int score;
  int ww;
  int hh;
};

int queryCount;
Score tScores[MAX_T];

void CopyToBest()
{
}

void CopyToAns()
{
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  queryCount = 0;
}

// 入力を受け取る関数
void Input(int problemNum)
{
  if (mode == 0) {
    // 標準入力
    cin >> n >> t >> sigma;
    rep(i, n)
    {
      cin >> w[i] >> h[i];
    }
  }
  else {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
    ifstream ifs(oss.str());
    // ファイル入力
    ifs >> n >> t >> sigma;
    rep(i, n)
    {
      ifs >> w[i] >> h[i];
    }
    rep(i, n)
    {
      ifs >> W[i] >> H[i];
    }
    rep(i, t)
    {
      ifs >> dW[i] >> dH[i];
    }
  }

  rep(i, n)
  {
    w[i] = max(MIN_WIDTH, w[i]);
    w[i] = min(MAX_WIDTH, w[i]);
    h[i] = max(MIN_HEIGHT, h[i]);
    h[i] = min(MAX_HEIGHT, h[i]);
  }

  if (mode == 4) {
    rep(i, n)
    {
      w[i] = W[i];
      h[i] = H[i];
    }
    rep(i, t)
    {
      dW[i] = 0;
      dH[i] = 0;
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

Score CalcTScore()
{
  Score score;
  score.score = INF;
  rep(i, queryCount)
  {
    if (tScores[i].score < score.score) {
      score = tScores[i];
    }
  }
  return score;
}

int cs_use[MAX_N] = {};
int cs_up[MAX_N], cs_down[MAX_N], cs_left[MAX_N], cs_right[MAX_N];
int cs_max_down[MAX_N], cs_max_right[MAX_N];
// スコアを計算する関数
Score CalcScore(const vector<Piece>& pieces, bool cheat)
{
  int sz = pieces.size();
  rep(i, n)
  {
    cs_use[i] = 0;
    cs_up[i] = -1;
    cs_down[i] = -1;
    cs_left[i] = -1;
    cs_right[i] = -1;
    cs_max_down[i] = -1;
    cs_max_right[i] = -1;
  }
  rep(i, sz)
  {
    Piece col = pieces[i];
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
      drep(j, num)
      {
        if (cs_use[j]) {
          if (isCrossing(cs_left[num], cs_right[num], cs_left[j], cs_right[j])) {
            cs_up[num] = max(cs_up[num], cs_down[j]);
          }
        }
        if (cs_max_down[j] <= cs_up[num])break;
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
      drep(j, num)
      {
        if (cs_use[j]) {
          if (isCrossing(cs_up[num], cs_down[num], cs_up[j], cs_down[j])) {
            cs_left[num] = max(cs_left[num], cs_right[j]);
          }
        }
        if (cs_max_right[j] <= cs_left[num])break;
      }
      cs_right[num] = cs_left[num] + wid;
    }

    cs_max_down[num] = cs_down[num];
    cs_max_right[num] = cs_right[num];
    if (i > 0) {
      int befNum = pieces[i - 1].num;
      srep(j, befNum + 1, num + 1)
      {
        cs_max_down[j] = max(cs_max_down[j], cs_max_down[befNum]);
        cs_max_right[j] = max(cs_max_right[j], cs_max_right[befNum]);
      }
    }

    cs_use[num] = 1;
  }

  int maxDown = 0, maxRight = 0;
  rep(i, n)
  {
    if (cs_use[i]) {
      maxDown = max(maxDown, cs_down[i]);
      maxRight = max(maxRight, cs_right[i]);
    }
  }

  Score score;

  score.ww = maxRight;
  score.hh = maxDown;
  score.score = maxDown + maxRight;
  rep(i, n)
  {
    if (cs_use[i] == 0) {
      score.score += w[i] + h[i];
    }
  }

  return score;
}

Score Print(const vector<Piece>& pieces, ofstream& ofs)
{
  Score score;

  if (mode == 0) {
    cout << pieces.size() << endl;
    rep(i, pieces.size())
    {
      cout << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }
    fflush(stdout);

    cin >> score.ww >> score.hh;
    score.score = score.ww + score.hh;
    tScores[queryCount] = score;
  }
  else {
    ofs << pieces.size() << endl;
    rep(i, pieces.size())
    {
      ofs << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }

    score = CalcScore(pieces, true);
    tScores[queryCount] = score;
    score.ww += dW[queryCount];
    score.hh += dH[queryCount];
    score.score = score.ww + score.hh;
  }

  queryCount++;

  return score;
}

int bestsCount;
vector<Piece> bests[MAX_T];
int bestScores[MAX_T];

void Method1_Shoki1_Internal1(vector<Piece>& tmp)
{
  int widSum = 0;
  int heiSum = 0;

  int widLimit = RandXor() % 1000000 + 200000;
  int now = 0;
  int maxHeight = 0;
  rep(i, n)
  {
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
void Method1_Shoki1_Internal2(vector<Piece>& tmp)
{
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
  rep(i, n)
  {
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

  rep(i, keepRotCount)
  {
    int ii = keepRot[i];
    tmp[ii].rot = 1 - tmp[ii].rot;
  }
}

void Method1_Shoki1()
{
  bestsCount = 0;

  vector<Piece> tmp(n);
  rep(i, n)
  {
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

void Method1(ofstream& ofs)
{
  vector<Piece> best(n);
  int bestScore = INF;

  vector<Piece> tmp(n);

  Method1_Shoki1();

  rep(aespa, bestsCount)
  {
    tmp = bests[aespa];
    Score score = Print(tmp, ofs);
    if (score.score < bestScore) {
      bestScore = score.score;
      best = tmp;
    }
  }

  srep(aespa, bestsCount, t)
  {
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
      rep(i, vvc.size())
      {
        rep(j, vvc[i].size())
        {
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
void InitializePieces()
{
  initialPieces.resize(n);
  rep(i, n)
  {
    initialPieces[i].num = i;
    initialPieces[i].rot = 0;
    if (w[i] > h[i]) initialPieces[i].rot = 1;
    initialPieces[i].dir = 0;
    initialPieces[i].base = i - 1;
  }
}

int sizeRank[MAX_N];
int sizeSize[MAX_N];
void InitializeSizeRank()
{
  vector<P> vp;
  rep(i, n)
  {
    vp.emplace_back(max(w[i], h[i]), i);
  }
  sort(vp.begin(), vp.end());
  rep(i, n)
  {
    sizeRank[vp[i].second] = i;
    sizeSize[vp[i].second] = vp[i].first;
  }
}

Block block;
Row row;
Board board;

int bestsCount2;
Board bests2[MAX_T];
int FIRST_ROUND = MAX_T;
void Method2_Shoki1_Internal1()
{
  FIRST_ROUND = t / 2;

  block.clear();
  row.clear();
  board.clear();

  int widLimit = RandXor() % 1000000 + 200000;

  rep(i, n)
  {
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
  if (bestsCount2 < FIRST_ROUND) {
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

void Method2_Shoki1_Internal3(double progressRatio)
{
  FIRST_ROUND = t / 2;

  block.clear();
  row.clear();
  board.clear();
  Piece piece;

  int widLimit = RandXor() % 1000000 + 200000;
  if (progressRatio > 0.5) {
    widLimit = RandXor() % 200000 - 100000 + bests2[0].GetMaxWidth();
  }

  int isRandomRot = RandXor() % 10;

  int isUsePiece2_1 = RandXor() % 200;

  double maxRatio = 0.9 + Rand01() * 0.2;

  int bestsLimit = t / 2;
  int ng = 0;

  int useSizeRank = RandXor() % 10 + 5;

  int isLie = RandXor() % 100;

  rep(i, n)
  {
    if (bestsCount2 == bestsLimit && board.score() + row.maxHeight >= bests2[bestsCount2 - 1].score()) {
      ng = 1;
      break;
    }

    piece = initialPieces[i];

    if (isRandomRot >= 1 && RandXor() % n <= isRandomRot) {
      piece.rot = 1 - piece.rot;
    }

    if (isUsePiece2_1 < 50 && (sizeRank[i] < useSizeRank || sizeSize[i] < MAX_HEIGHT * maxRatio / 2)) {
      if (row.sz > 0 && row.blocks[row.sz - 1].count() == 1 && row.blocks[row.sz - 1].height() + piece.height() < MAX_HEIGHT * maxRatio) {
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
    else if (isUsePiece2_1 < 100 && (sizeRank[i] < useSizeRank || sizeSize[i] < MAX_HEIGHT * maxRatio / 2)) {
      if (row.sz >= 1 && row.blocks[row.sz - 1].count() == 1 && row.blocks[row.sz - 1].height() + piece.height() < MAX_HEIGHT * maxRatio) {
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
        rep(j, row.sz - 1)
        {
          if (row.blocks[j].count() == 1 && row.blocks[j].height() + piece.height() < MAX_HEIGHT * maxRatio && piece.width() < row.blocks[j].width()) {
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
    else if (isUsePiece2_1 < 150 && (sizeRank[i] < useSizeRank || sizeSize[i] < MAX_HEIGHT * maxRatio / 2)) {
      if (row.sz >= 1 && row.blocks[row.sz - 1].count() == 1 && row.blocks[row.sz - 1].height() + piece.height() < MAX_HEIGHT * maxRatio) {
        int ok = 1;
        if (RandXor() % 2 == 0) {
          ok = 0;
        }
        if (ok) {
          row.addPiece(row.sz - 1, piece);
          continue;
        }
      }

      if (board.sz >= 1) {
        int isAdd = 0;
        int nowWidth = 0;
        rep(j, board.rows[board.sz - 1].sz - 1)
        {
          if (nowWidth < row.sumWidth) {
            nowWidth += board.rows[board.sz - 1].blocks[j].width();
            continue;
          }
          if (board.rows[board.sz - 1].blocks[j].count() == 1
            && board.rows[board.sz - 1].blocks[j].height() + piece.height() < MAX_HEIGHT * maxRatio
            && piece.width() < board.rows[board.sz - 1].blocks[j].width()) {
            int ok = 1;
            if (RandXor() % 2 == 0) {
              ok = 0;
            }
            if (ok) {
              board.rows[board.sz - 1].addPiece(j, piece);
              isAdd = 1;
              break;
            }
          }
        }
        if (isAdd) {
          continue;
        }
      }

      if (row.sz >= 2) {
        int isAdd = 0;
        rep(j, row.sz - 1)
        {
          if (row.blocks[j].count() == 1 && row.blocks[j].height() + piece.height() < MAX_HEIGHT * maxRatio && piece.width() < row.blocks[j].width()) {
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

    if (isLie < 100) {
      if (row.sz > 0 && row.blocks[row.sz - 1].count() == 1) {
        int height = max(row.blocks[row.sz - 1].height(), piece.height());
        int newHeight = max(row.blocks[row.sz - 1].width(), piece.width());
        if (abs(row.blocks[row.sz - 1].height() - piece.height()) < 20000 && 70000 < newHeight && newHeight < MAX_HEIGHT * 1.5) {
          if (RandXor() % 50 < isLie) {
            row.rotateBack();
            piece.rot = 1 - piece.rot;
            row.addPiece(row.sz - 1, piece);
            continue;
          }
        }
      }
    }

    block.piece1 = piece;
    if (row.GetSumWidth() < widLimit * 0.9 && board.count() >= 1 && board.back().GetSumWidth() + block.piece1.width() <= widLimit && RandXor() % 4 != 0) {
      board.Add(board.count() - 1, block);
    }
    else if (row.GetSumWidth() + block.piece1.width() <= widLimit && RandXor() % (n) != 0) {
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
  if (bestsCount2 < FIRST_ROUND) {
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


void Method2_Shoki1()
{
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

  double proglessRatio = 0.0;
  double timeLimit =  TL / 30 * 25;
  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      proglessRatio = nowTime / timeLimit;
      if (nowTime > timeLimit) {
        break;
      }
    }

    loop++;

    Method2_Shoki1_Internal3(proglessRatio);
  }

  if (mode >= 2) {
    cout << "loop = " << loop << endl;
  }
}

struct Ans
{
  vector<Piece> pieces;
  int score;

  // for sorting
  bool operator<(const Ans& other) const
  {
    return score < other.score;
  }
};

void Method2(ofstream& ofs)
{
  Method2_Shoki1();

  vector<Ans> answers;
  rep(i, bestsCount2)
  {
    Ans ans;
    ans.pieces = bests2[i].CreateQuery();
    ans.score = CalcScore(ans.pieces, false).score;
    answers.push_back(ans);
  }

  int karina = 0;
  double timeLimit = TL * 25 / 30;
  while (false) {
    karina++;
    if (karina % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > timeLimit) {
        break;
      }
    }

    int raQ = RandXor() % bestsCount2;

    int raMode = RandXor() % 3;
    int ra = RandXor() % n;
    Piece keep = answers[raQ].pieces[ra];

    if (raMode == 0) {
      answers[raQ].pieces[ra].rot = 1 - answers[raQ].pieces[ra].rot;
    }
    else if (raMode == 1) {
      int ra2 = RandXor() % (ra + 1) - 1;
      if (answers[raQ].pieces[ra].base == ra2)continue;
      answers[raQ].pieces[ra].base = ra2;
    }
    else {
      answers[raQ].pieces[ra].dir = 1 - answers[raQ].pieces[ra].dir;
      if (RandXor() % 2 == 0) {
        int ra2 = RandXor() % (ra + 1) - 1;
        answers[raQ].pieces[ra].base = ra2;
      }
    }

    auto preScore = CalcScore(answers[raQ].pieces, false);

    if (preScore.score <= answers[raQ].score) {
      answers[raQ].score = preScore.score;
    }
    else {
      answers[raQ].pieces[ra] = keep;
    }
  }

  if (mode >= 2) {
    cout << "karina = " << karina << endl;
  }

  sort(answers.begin(), answers.end());

  int aespa = 0;
  int SECOND_ROUND = t / 2;
  while (aespa < SECOND_ROUND) {
    answers[aespa].score = Print(answers[aespa].pieces, ofs).score;
    aespa++;
  }

  sort(answers.begin(), answers.end());

  int loop = 0;
  double startTime = GetNowTime();
  double nowTime = GetNowTime();  // 現在の経過時間
  const double START_TEMP = 0.0;  // 焼きなまし法の開始温度
  const double END_TEMP = 200000.0;      // 焼きなまし法の終了温度
  timeLimit = TL - startTime;
  double temp = START_TEMP;  // 現在の温度

  int raQ = 0;
  vector<Piece> keep = answers[raQ].pieces;
  int raQCount = min(bestsCount2, 3);
  while (aespa < t) {
    loop++;
    if (loop % 100 == 0) {
      nowTime = GetNowTime();
    }

    int raMode = RandXor() % 3;

    if (RandXor() % 2 == 0) {
      answers[raQ].pieces = keep;
      raQ = RandXor() % raQCount;
      keep = answers[raQ].pieces;
    }

    if (raMode == 0) {
      int ra = RandXor() % n;
      answers[raQ].pieces[ra].rot = 1 - answers[raQ].pieces[ra].rot;
    }
    else if (raMode == 1) {
      int ra = RandXor() % n;
      int ra2 = RandXor() % (ra + 1) - 1;
      if (answers[raQ].pieces[ra].base == ra2)continue;
      answers[raQ].pieces[ra].base = ra2;
    }
    else {
      int ra = RandXor() % n;
      answers[raQ].pieces[ra].dir = 1 - answers[raQ].pieces[ra].dir;
      if (RandXor() % 2 == 0) {
        int ra2 = RandXor() % (ra + 1) - 1;
        answers[raQ].pieces[ra].base = ra2;
      }
    }

    auto preScore = CalcScore(answers[raQ].pieces, false);

    // だんだん温度が上がる焼きなまし
    const double progressRatio = (nowTime - startTime) / timeLimit;  // 進捗率（0.0～1.0）
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio * progressRatio * progressRatio;  // 温度の更新
    double diff = answers[raQ].score - preScore.score - 100000;  // スコアの差分
    double prob = exp(diff / temp);     // 焼きなまし法の採用確率
    if (preScore.score < answers[raQ].score || prob > Rand01() || nowTime > TL) {
      Score score = Print(answers[raQ].pieces, ofs);
      if (score.score < answers[raQ].score) {
        if (mode >= 2) {
          cout << raQ << ' ' << raMode << ' ' << diff << ' ' << temp << ' ' << prob << ' ' << nowTime << endl;
        }
        answers[raQ].score = score.score;
        keep = answers[raQ].pieces;
      }
      aespa++;
    }
  }

  if (mode >= 2) {
    cout << "aespa loop = " << loop << endl;
  }
}

bool use[MAX_T * 2][MAX_N * 2];
int sum[MAX_T * 2];
int measure[MAX_T * 2];
void Yamanobori(ofstream& ofs)
{
  int TT = t / 3;

  int originH[MAX_N], originW[MAX_N];
  rep(i, n)
  {
    originW[i] = w[i];
    originH[i] = h[i];
  }

  rep(i, TT * 2)
  {
    rep(j, n * 2)
    {
      use[i][j] = false;
    }
  }

  rep(i, TT * 2)
  {
    sum[i] = 0;
    measure[i] = 0;
  }

  int useCount[MAX_N * 2];
  rep(i, n * 2)useCount[i] = 0;

  rep(aespa, TT)
  {
    vector<Piece> cols;
    rep(i, n)
    {
      if (RandXor() % 1 == 0) {
        Piece col;
        col.num = i;
        col.rot = RandXor() % 2;
        col.dir = RandXor() % 2;
        col.base = -1;
        cols.push_back(col);
      }
    }

    Score score1 = Print(cols, ofs);
    Score score2 = Print(cols, ofs);
    Score score3 = Print(cols, ofs);

    measure[aespa * 2] = (score1.hh + score2.hh + score3.hh) / 3;
    measure[aespa * 2 + 1] = (score1.ww + score2.ww + score3.ww) / 3;
    bool head = true;
    for (auto col : cols) {
      if (head) {
        if (col.rot == 0) {
          sum[aespa * 2] += h[col.num];
          use[aespa * 2][col.num * 2] = 1;
          sum[aespa * 2 + 1] += w[col.num];
          use[aespa * 2 + 1][col.num * 2 + 1] = 1;

          useCount[col.num * 2]++;
          useCount[col.num * 2 + 1]++;
        }
        else {
          sum[aespa * 2] += w[col.num];
          use[aespa * 2][col.num * 2 + 1] = 1;
          sum[aespa * 2 + 1] += h[col.num];
          use[aespa * 2 + 1][col.num * 2] = 1;

          useCount[col.num * 2]++;
          useCount[col.num * 2 + 1]++;
        }
        head = false;
      }
      else {
        if (col.dir == 0) {
          if (col.rot == 0) {
            sum[aespa * 2] += h[col.num];
            use[aespa * 2][col.num * 2] = 1;

            useCount[col.num * 2]++;
          }
          else {
            sum[aespa * 2] += w[col.num];
            use[aespa * 2][col.num * 2 + 1] = 1;

            useCount[col.num * 2 + 1]++;
          }
        }
        else {
          if (col.rot == 0) {
            sum[aespa * 2 + 1] += w[col.num];
            use[aespa * 2 + 1][col.num * 2 + 1] = 1;

            useCount[col.num * 2 + 1]++;
          }
          else {
            sum[aespa * 2 + 1] += h[col.num];
            use[aespa * 2 + 1][col.num * 2] = 1;

            useCount[col.num * 2]++;
          }
        }
      }
    }
  }

  double timeLimit = TL / 3;
  int loop = 0;
  double nowTime = GetNowTime();  // 現在の経過時間

  double diffSum[MAX_N * 2];

  while (true) {
    if (loop % 100 == 0) {
      auto nowTime = GetNowTime();
      if (nowTime > timeLimit) {
        break;
      }
    }
    loop++;

    rep(i, n * 2)
    {
      diffSum[i] = 0;
    }

    int diffSize = RandXor() % 3;
    int raDiff = RandXor() % 2001 - 1000;
    if (diffSize == 1) {
      raDiff = RandXor() % 201 - 100;
    }
    else if (diffSize == 2) {
      raDiff = RandXor() % 21 - 10;
    }

    double maxDiffSum = -INF;
    int argMax = -1;

    rep(i, n * 2)
    {
      if (useCount[i] == 0)continue;
      if (i % 2 == 0) {
        if (h[i / 2] + raDiff > MAX_HEIGHT)continue;
        if (h[i / 2] + raDiff < MIN_HEIGHT)continue;
      }
      else {
        if (w[i / 2] + raDiff > MAX_WIDTH)continue;
        if (w[i / 2] + raDiff < MIN_WIDTH)continue;
      }

      rep(j, TT * 2)
      {
        if (use[j][i]) {
          double beforeDiff = abs(sum[j] - measure[j]);
          double afterDiff = abs((sum[j] + raDiff) - measure[j]);
          diffSum[i] += (beforeDiff - afterDiff) / useCount[i];
        }
      }

      if (maxDiffSum < diffSum[i]) {
        maxDiffSum = diffSum[i];
        argMax = i;
      }
    }

    if (maxDiffSum >= 0) {
      if (argMax % 2 == 0) {
        h[argMax / 2] += raDiff;
      }
      else {
        w[argMax / 2] += raDiff;
      }

      rep(j, TT * 2)
      {
        if (use[j][argMax]) {
          sum[j] += raDiff;
        }
      }
    }
  }

  if (mode >= 2) {
    int diffSum1 = 0, diffSum2 = 0;
    //cout << "元 最終 正解" << endl;
    rep(i, n)
    {
      cout << originH[i] << ' ' << h[i] << ' ' << H[i] << endl;
      diffSum1 += abs(originH[i] - H[i]);
      diffSum2 += abs(h[i] - H[i]);
    }
    cout << "loop = " << loop << ", t = " << TT << ' ';
    cout << diffSum1 << "  " << diffSum2 << endl;
  }
}

// 問題を解く関数
ll Solve(int problem_num)
{
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
  //Yamanobori(ofs);
  Method2(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  Score score = CalcTScore();
  if (mode >= 2) {
    cout << "ww = " << score.ww << ", hh = " << score.hh << ", score = " << score.score << endl;
  }
  return score.score;
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

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 3) {
    rep(_, 10)
    {
      ll sum = 0;
      srep(i, 0, 100)
      {
        ll score = Solve(i);
        sum += score;
      }
      cout << sum << endl;
    }
  }
  else {
    ll sum = 0;
    srep(i, 0, 100)
    {
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
    cout << "sum = " << sum << endl;
  }

  return 0;
}
