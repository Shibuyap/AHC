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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

static uint32_t RandXor()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;
  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

static double Rand01()
{
  return (RandXor() + 0.5) * (1.0 / UINT_MAX);
}

static double RandUniform(double l, double r)
{
  return l + (r - l) * Rand01();
}

void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = RandXor() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

const int INF = 1001001001;
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

double TL = 2.8;
int executionMode;
std::chrono::steady_clock::time_point globalStartTimePoint;

void ResetTime()
{
  globalStartTimePoint = std::chrono::steady_clock::now();
}

double GetNowTime()
{
  auto endPoint = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = endPoint - globalStartTimePoint;
  return elapsed.count();
}

bool isCrossing(double l1, double r1, double l2, double r2)
{
  return (std::max(l1, l2) < std::min(r1, r2));
}

class RectanglePiece
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

  bool operator<(const RectanglePiece& other) const
  {
    return num < other.num;
  }
};

class StackUnit
{
public:
  RectanglePiece piece1;
  RectanglePiece piece2;

  StackUnit()
  {
    piece1.num = -1;
    piece2.num = -1;
  }

  int count() const
  {
    if (piece1.num == -1) return 0;
    else if (piece2.num == -1) return 1;
    return 2;
  }

  int width() const
  {
    int cnt = count();
    if (cnt == 0) return 0;
    else if (cnt == 1) return piece1.width();
    return max(piece1.width(), piece2.width());
  }

  int height() const
  {
    int cnt = count();
    if (cnt == 0) return 0;
    else if (cnt == 1) return piece1.height();
    return piece1.height() + piece2.height();
  }

  void clear()
  {
    piece1.num = -1;
    piece2.num = -1;
  }

  void SetBase(int baseIndex)
  {
    piece1.base = baseIndex;
    piece2.base = baseIndex;
  }
};

class Shelf
{
public:
  StackUnit blocks[MAX_N];
  int sz;
  int sumWidth;
  int maxHeight;

  Shelf()
  {
    sz = 0; sumWidth = 0; maxHeight = 0;
  }

  int count() const { return sz; }

  void clear()
  {
    sz = 0; sumWidth = 0; maxHeight = 0;
  }

  void add(const StackUnit& block)
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

  void addPiece(int index, const RectanglePiece& piece)
  {
    int beforeWidth = blocks[index].width();
    blocks[index].piece2 = piece;
    blocks[index].piece2.base = blocks[index].piece1.base;
    sumWidth += blocks[index].width() - beforeWidth;
    maxHeight = max(maxHeight, blocks[index].height());
  }

  int GetSumWidth() const { return sumWidth; }
  int GetMaxHeight() const { return maxHeight; }
};

class Layout
{
public:
  Shelf shelves[MAX_N];
  int sz;
  int maxWidth;
  int sumHeight;

  Layout()
  {
    sz = 0; maxWidth = 0; sumHeight = 0;
  }

  void clear()
  {
    sz = 0; maxWidth = 0; sumHeight = 0;
  }

  int count() const { return sz; }

  int score() const
  {
    return maxWidth + sumHeight;
  }

  void Add(const Shelf& shelf)
  {
    shelves[sz] = shelf;
    sz++;
    maxWidth = max(maxWidth, shelf.GetSumWidth());
    sumHeight += shelf.GetMaxHeight();
  }

  void Add(int index, const StackUnit& block)
  {
    int beforeMaxHeight = shelves[index].GetMaxHeight();
    shelves[index].add(block);
    maxWidth = max(maxWidth, shelves[index].GetSumWidth());
    sumHeight += shelves[index].GetMaxHeight() - beforeMaxHeight;
  }

  int GetMaxWidth() const { return maxWidth; }
  int GetSumHeight() const { return sumHeight; }

  vector<RectanglePiece> CreateQuery() const
  {
    vector<RectanglePiece> pieces;
    rep(i, sz)
    {
      rep(j, shelves[i].sz)
      {
        if (shelves[i].blocks[j].count() >= 1) {
          pieces.push_back(shelves[i].blocks[j].piece1);
        }
        if (shelves[i].blocks[j].count() >= 2) {
          pieces.push_back(shelves[i].blocks[j].piece2);
        }
      }
    }
    sort(pieces.begin(), pieces.end());
    return pieces;
  }

  Shelf& back()
  {
    return shelves[sz - 1];
  }
};

struct ScoreStruct
{
  int score;
  int ww;
  int hh;
};

int queryCounter;
ScoreStruct tScores[MAX_T];

void InitializeGlobalState()
{
  queryCounter = 0;
}

void LoadInputData(int problemNum)
{
  if (executionMode == 0) {
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

  if (executionMode == 4) {
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

void OpenOutputStream(int probNum, ofstream& ofs)
{
  if (executionMode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

ScoreStruct FindBestQueryScore()
{
  ScoreStruct score;
  score.score = INF;
  rep(i, queryCounter)
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

ScoreStruct EvaluateScore(const vector<RectanglePiece>& pieces, bool cheat)
{
  int sz = (int)pieces.size();
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

  RectanglePiece currentRectPiece;
  rep(i, sz)
  {
    currentRectPiece = pieces[i];
    int num = currentRectPiece.num;
    int currentWidth = w[num];
    int currentHeight = h[num];
    if (executionMode != 0 && cheat == true) {
      currentWidth = W[num];
      currentHeight = H[num];
    }
    if (currentRectPiece.rot == 1) swap(currentWidth, currentHeight);
    if (currentRectPiece.dir == 0) {
      cs_left[num] = 0;
      if (currentRectPiece.base != -1) {
        cs_left[num] = cs_right[currentRectPiece.base];
      }
      cs_right[num] = cs_left[num] + currentWidth;

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
      cs_down[num] = cs_up[num] + currentHeight;
    }
    else {
      cs_up[num] = 0;
      if (currentRectPiece.base != -1) {
        cs_up[num] = cs_down[currentRectPiece.base];
      }
      cs_down[num] = cs_up[num] + currentHeight;

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
      cs_right[num] = cs_left[num] + currentWidth;
    }

    cs_max_down[num] = cs_down[num];
    cs_max_right[num] = cs_right[num];
    if (i > 0) {
      int previousNum = pieces[i - 1].num;
      srep(j, previousNum + 1, num + 1)
      {
        cs_max_down[j] = max(cs_max_down[j], cs_max_down[previousNum]);
        cs_max_right[j] = max(cs_max_right[j], cs_max_right[previousNum]);
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

  ScoreStruct scoreResult;
  scoreResult.ww = maxRight;
  scoreResult.hh = maxDown;
  scoreResult.score = maxDown + maxRight;
  rep(i, n)
  {
    if (cs_use[i] == 0) {
      scoreResult.score += w[i] + h[i];
    }
  }

  return scoreResult;
}

ScoreStruct PrintFinalArrangement(const vector<RectanglePiece>& pieces, ofstream& ofs)
{
  ScoreStruct resultScore;

  if (executionMode == 0) {
    cout << pieces.size() << endl;
    rep(i, pieces.size())
    {
      cout << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }
    fflush(stdout);

    cin >> resultScore.ww >> resultScore.hh;
    resultScore.score = resultScore.ww + resultScore.hh;
    tScores[queryCounter] = resultScore;
  }
  else {
    ofs << "# " << pieces.size() << endl;
    rep(i, pieces.size())
    {
      ofs << "# " << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }

    ofs << pieces.size() << endl;
    rep(i, pieces.size())
    {
      ofs << pieces[i].num << ' ' << pieces[i].rot << ' ' << (pieces[i].dir == 0 ? 'U' : 'L') << ' ' << pieces[i].base << endl;
    }

    resultScore = EvaluateScore(pieces, true);
    tScores[queryCounter] = resultScore;
    resultScore.ww += dW[queryCounter];
    resultScore.hh += dH[queryCounter];
    resultScore.score = resultScore.ww + resultScore.hh;
  }

  queryCounter++;

  return resultScore;
}

vector<RectanglePiece> basePieces;
void SetupBasePieces()
{
  basePieces.resize(n);
  rep(i, n)
  {
    basePieces[i].num = i;
    basePieces[i].rot = 0;
    if (w[i] > h[i]) basePieces[i].rot = 1;
    basePieces[i].dir = 0;
    basePieces[i].base = i - 1;
  }
}

int pieceSizeOrder[MAX_N];
int pieceLargestDimension[MAX_N];
void SetupPieceSizeOrder()
{
  vector<P> pieceDimensionIndexPairs;
  rep(i, n)
  {
    pieceDimensionIndexPairs.emplace_back(max(w[i], h[i]), i);
  }
  sort(pieceDimensionIndexPairs.begin(), pieceDimensionIndexPairs.end());
  rep(i, n)
  {
    pieceSizeOrder[pieceDimensionIndexPairs[i].second] = i;
    pieceLargestDimension[pieceDimensionIndexPairs[i].second] = pieceDimensionIndexPairs[i].first;
  }
}

StackUnit currentStackUnit;
Shelf currentShelf;
Layout currentLayout;

int candidateLayoutCount;
Layout candidateLayouts[MAX_T];
int INITIAL_ROUND_LIMIT = MAX_T;

void BuildInitialMethod2LayoutsPhase1()
{
  INITIAL_ROUND_LIMIT = t / 2;

  currentStackUnit.clear();
  currentShelf.clear();
  currentLayout.clear();

  int widthLimit = RandXor() % 1000000 + 200000;

  rep(i, n)
  {
    currentStackUnit.piece1 = basePieces[i];
    if (currentShelf.GetSumWidth() + currentStackUnit.piece1.width() <= widthLimit) {
      currentShelf.add(currentStackUnit);
    }
    else {
      currentLayout.Add(currentShelf);
      currentShelf.clear();
      currentShelf.add(currentStackUnit);
    }
  }

  currentLayout.Add(currentShelf);

  if (candidateLayoutCount < INITIAL_ROUND_LIMIT) {
    candidateLayouts[candidateLayoutCount] = currentLayout;
    candidateLayoutCount++;
  }
  else if (currentLayout.score() < candidateLayouts[candidateLayoutCount - 1].score()) {
    candidateLayouts[candidateLayoutCount - 1] = currentLayout;
  }
  int pos = candidateLayoutCount - 1;
  while (pos >= 1) {
    if (candidateLayouts[pos].score() < candidateLayouts[pos - 1].score()) {
      swap(candidateLayouts[pos], candidateLayouts[pos - 1]);
      pos--;
    }
    else {
      break;
    }
  }
}

void UpdateMethod2LayoutsPhase(double progressRatio)
{
  INITIAL_ROUND_LIMIT = t / 2;

  currentStackUnit.clear();
  currentShelf.clear();
  currentLayout.clear();
  RectanglePiece piece;

  int widthLimit = RandXor() % 1000000 + 200000;
  if (progressRatio > 0.5) {
    widthLimit = RandXor() % 200000 - 100000 + candidateLayouts[0].GetMaxWidth();
  }

  int rotationProbabilityThreshold = RandXor() % 10;
  int secondStackingThreshold = RandXor() % 200;
  double heightRatioLimit = 0.9 + Rand01() * 0.2;
  int bestsLimit = t / 2;
  int terminateFlag = 0;

  int rankThreshold = RandXor() % 10 + 5;
  int lyingConditionThreshold = RandXor() % 100;

  rep(i, n)
  {
    if (candidateLayoutCount == bestsLimit && currentLayout.score() + currentShelf.maxHeight >= candidateLayouts[candidateLayoutCount - 1].score()) {
      terminateFlag = 1;
      break;
    }

    piece = basePieces[i];

    if (rotationProbabilityThreshold >= 1 && RandXor() % n <= rotationProbabilityThreshold) {
      piece.rot = 1 - piece.rot;
    }

    // ������isPlacementAllowed��p���ĉǐ�����
    // �S�Ă�ok�ϐ���isPlacementAllowed�ɕύX
    int isPlacementAllowed = 1;

    if (secondStackingThreshold < 50 && (pieceSizeOrder[i] < rankThreshold || pieceLargestDimension[i] < MAX_HEIGHT * heightRatioLimit / 2)) {
      if (currentShelf.sz > 0 && currentShelf.blocks[currentShelf.sz - 1].count() == 1 && currentShelf.blocks[currentShelf.sz - 1].height() + piece.height() < MAX_HEIGHT * heightRatioLimit) {
        isPlacementAllowed = 1;
        if (RandXor() % 25 > secondStackingThreshold) {
          isPlacementAllowed = 0;
        }
        if (isPlacementAllowed) {
          currentShelf.addPiece(currentShelf.sz - 1, piece);
          continue;
        }
      }
    }
    else if (secondStackingThreshold < 100 && (pieceSizeOrder[i] < rankThreshold || pieceLargestDimension[i] < MAX_HEIGHT * heightRatioLimit / 2)) {
      if (currentShelf.sz >= 1 && currentShelf.blocks[currentShelf.sz - 1].count() == 1 && currentShelf.blocks[currentShelf.sz - 1].height() + piece.height() < MAX_HEIGHT * heightRatioLimit) {
        isPlacementAllowed = 1;
        if (RandXor() % 2 == 0) {
          isPlacementAllowed = 0;
        }
        if (isPlacementAllowed) {
          currentShelf.addPiece(currentShelf.sz - 1, piece);
          continue;
        }
      }

      {
        int isAdd = 0;
        if (currentShelf.sz >= 2) {
          rep(j, currentShelf.sz - 1)
          {
            if (currentShelf.blocks[j].count() == 1 && currentShelf.blocks[j].height() + piece.height() < MAX_HEIGHT * heightRatioLimit && piece.width() < currentShelf.blocks[j].width()) {
              isPlacementAllowed = 1;
              if (RandXor() % 2 == 0) {
                isPlacementAllowed = 0;
              }
              if (isPlacementAllowed) {
                currentShelf.addPiece(j, piece);
                isAdd = 1;
                break;
              }
            }
          }
        }
        if (isAdd) {
          continue;
        }
      }
    }
    else if (secondStackingThreshold < 150 && (pieceSizeOrder[i] < rankThreshold || pieceLargestDimension[i] < MAX_HEIGHT * heightRatioLimit / 2)) {
      if (currentShelf.sz >= 1 && currentShelf.blocks[currentShelf.sz - 1].count() == 1 && currentShelf.blocks[currentShelf.sz - 1].height() + piece.height() < MAX_HEIGHT * heightRatioLimit) {
        isPlacementAllowed = 1;
        if (RandXor() % 2 == 0) {
          isPlacementAllowed = 0;
        }
        if (isPlacementAllowed) {
          currentShelf.addPiece(currentShelf.sz - 1, piece);
          continue;
        }
      }

      {
        int isAdd = 0;
        if (currentLayout.sz >= 1) {
          int accumulatedWidth = 0;
          rep(j, currentLayout.shelves[currentLayout.sz - 1].sz - 1)
          {
            if (accumulatedWidth < currentShelf.sumWidth) {
              accumulatedWidth += currentLayout.shelves[currentLayout.sz - 1].blocks[j].width();
              continue;
            }
            if (currentLayout.shelves[currentLayout.sz - 1].blocks[j].count() == 1
              && currentLayout.shelves[currentLayout.sz - 1].blocks[j].height() + piece.height() < MAX_HEIGHT * heightRatioLimit
              && piece.width() < currentLayout.shelves[currentLayout.sz - 1].blocks[j].width()) {
              isPlacementAllowed = 1;
              if (RandXor() % 2 == 0) {
                isPlacementAllowed = 0;
              }
              if (isPlacementAllowed) {
                currentLayout.shelves[currentLayout.sz - 1].addPiece(j, piece);
                isAdd = 1;
                break;
              }
            }
          }
          if (isAdd) {
            continue;
          }
        }

        if (currentShelf.sz >= 2) {
          isAdd = 0;
          rep(j, currentShelf.sz - 1)
          {
            if (currentShelf.blocks[j].count() == 1 && currentShelf.blocks[j].height() + piece.height() < MAX_HEIGHT * heightRatioLimit && piece.width() < currentShelf.blocks[j].width()) {
              isPlacementAllowed = 1;
              if (RandXor() % 2 == 0) {
                isPlacementAllowed = 0;
              }
              if (isPlacementAllowed) {
                currentShelf.addPiece(j, piece);
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
    }

    if (lyingConditionThreshold < 100) {
      if (currentShelf.sz > 0 && currentShelf.blocks[currentShelf.sz - 1].count() == 1) {
        int heightVal = max(currentShelf.blocks[currentShelf.sz - 1].height(), piece.height());
        int newHeight = max(currentShelf.blocks[currentShelf.sz - 1].width(), piece.width());
        if (abs(currentShelf.blocks[currentShelf.sz - 1].height() - piece.height()) < 20000 && 70000 < newHeight && newHeight < MAX_HEIGHT * 1.5) {
          isPlacementAllowed = 1;
          if (RandXor() % 50 < lyingConditionThreshold) {
            currentShelf.rotateBack();
            piece.rot = 1 - piece.rot;
            currentShelf.addPiece(currentShelf.sz - 1, piece);
            continue;
          }
        }
      }
    }

    currentStackUnit.piece1 = piece;
    if (currentShelf.GetSumWidth() < widthLimit * 0.9 && currentLayout.count() >= 1 && currentLayout.back().GetSumWidth() + currentStackUnit.piece1.width() <= widthLimit && RandXor() % 4 != 0) {
      currentLayout.Add(currentLayout.count() - 1, currentStackUnit);
    }
    else if (currentShelf.GetSumWidth() + currentStackUnit.piece1.width() <= widthLimit && RandXor() % (n) != 0) {
      currentShelf.add(currentStackUnit);
    }
    else {
      currentLayout.Add(currentShelf);
      currentShelf.clear();
      currentShelf.add(currentStackUnit);
    }
  }

  if (terminateFlag) {
    return;
  }

  currentLayout.Add(currentShelf);

  if (candidateLayoutCount < INITIAL_ROUND_LIMIT) {
    candidateLayouts[candidateLayoutCount] = currentLayout;
    candidateLayoutCount++;
  }
  else if (currentLayout.score() < candidateLayouts[candidateLayoutCount - 1].score()) {
    candidateLayouts[candidateLayoutCount - 1] = currentLayout;
  }
  int pos = candidateLayoutCount - 1;
  while (pos >= 1) {
    if (candidateLayouts[pos].score() < candidateLayouts[pos - 1].score()) {
      swap(candidateLayouts[pos], candidateLayouts[pos - 1]);
      pos--;
    }
    else {
      break;
    }
  }
}

void BuildInitialMethod2Layouts()
{
  candidateLayoutCount = 0;

  int processingLoopCount = 0;
  while (true) {
    if (processingLoopCount % 100 == 0) {
      auto currentElapsedTime = GetNowTime();
      if (currentElapsedTime > TL / 30) {
        break;
      }
    }

    processingLoopCount++;

    BuildInitialMethod2LayoutsPhase1();
  }

  double progressRatio = 0.0;
  double timeLimit = TL / 30 * 25;
  while (true) {
    if (processingLoopCount % 100 == 0) {
      auto currentElapsedTime = GetNowTime();
      progressRatio = currentElapsedTime / timeLimit;
      if (currentElapsedTime > timeLimit) {
        break;
      }
    }

    processingLoopCount++;

    UpdateMethod2LayoutsPhase(progressRatio);
  }

  if (executionMode >= 2) {
    cout << "loop = " << processingLoopCount << endl;
  }
}

struct Ans
{
  vector<RectanglePiece> pieces;
  int score;

  bool operator<(const Ans& other) const
  {
    return score < other.score;
  }
};

void RefineAndPrintSolutions(ofstream& ofs)
{
  BuildInitialMethod2Layouts();

  vector<Ans> solutionCandidates;
  rep(i, candidateLayoutCount)
  {
    Ans ans;
    ans.pieces = candidateLayouts[i].CreateQuery();
    ans.score = EvaluateScore(ans.pieces, false).score;
    solutionCandidates.push_back(ans);
  }

  int improvementStepCount = 0;
  double timeLimit = TL * 27 / 30;

  Ans keepAns;
  while (false) {
    improvementStepCount++;
    if (improvementStepCount % 100 == 0) {
      auto currentElapsedTime = GetNowTime();
      if (currentElapsedTime > timeLimit) {
        break;
      }
    }

    int randomCandidateIndex = RandXor() % candidateLayoutCount;
    int randomModificationMode = RandXor() % 4;
    int randomPieceIndex = RandXor() % n;
    RectanglePiece keep = solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex];

    if (randomModificationMode == 0) {
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].rot = 1 - solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].rot;
    }
    else if (randomModificationMode == 1) {
      int newBaseIndex = RandXor() % (randomPieceIndex + 1) - 1;
      if (solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base == newBaseIndex)continue;
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base = newBaseIndex;
    }
    else if (randomModificationMode == 2) {
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].dir = 1 - solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].dir;
      if (RandXor() % 2 == 0) {
        int newBaseIndex = RandXor() % (randomPieceIndex + 1) - 1;
        solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base = newBaseIndex;
      }
    }
    else if (randomModificationMode == 3) {
      keepAns = solutionCandidates[randomCandidateIndex];

      int randomPieceIndex2 = RandXor() % n;
      while (randomPieceIndex2 == randomPieceIndex) {
        randomPieceIndex2 = RandXor() % n;
      }
      if (randomPieceIndex2 < randomPieceIndex) {
        swap(randomPieceIndex, randomPieceIndex2);
      }

      if (solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex2].base == randomPieceIndex) {
        continue;
      }

      // 1��2�̍s��1�ԑO
      int base1 = solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex2].base;
      while (base1 > randomCandidateIndex) {
        base1 = solutionCandidates[randomCandidateIndex].pieces[base1].base;
      }
      if (base1 == randomCandidateIndex) {
        continue;
      }

      // 2��1�̍s��1�Ԍ��
      int base2_1 =  randomCandidateIndex;
      int base2_21 =  -2;
      int base2_22 =  -2;
      srep(i, randomCandidateIndex + 1, randomPieceIndex2)
      {
        if (solutionCandidates[randomCandidateIndex].pieces[i].base == base2_1) {
          if (base2_21 == -2) {
            base2_21 = i;
          }
          else {
            base2_22 = i;
          }
        }
        else if (solutionCandidates[randomCandidateIndex].pieces[i].base == base2_21) {
          base2_1 = base2_21;
          base2_21 = i;
          base2_22 = -2;
        }
        else if (solutionCandidates[randomCandidateIndex].pieces[i].base == base2_22) {
          base2_1 = base2_22;
          base2_21 = i;
          base2_22 = -2;
        }
      }

      int base2 = solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base;
      if (base2_21 != -2) {
        base2 = base2_21;
      }

      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base = base1;
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex2].base = base2;
    }

    auto preScore = EvaluateScore(solutionCandidates[randomCandidateIndex].pieces, false);

    if (preScore.score <= solutionCandidates[randomCandidateIndex].score) {
      solutionCandidates[randomCandidateIndex].score = preScore.score;
    }
    else {
      if (randomModificationMode == 3) {
        solutionCandidates[randomCandidateIndex] = keepAns;
      }
      else {
        solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex] = keep;
      }
    }
  }

  if (executionMode >= 2) {
    cout << "improvementStepCount = " << improvementStepCount << endl;
  }

  sort(solutionCandidates.begin(), solutionCandidates.end());

  int printedSolutionCount = 0;
  int SECOND_ROUND = t / 2;
  while (printedSolutionCount < SECOND_ROUND) {
    solutionCandidates[printedSolutionCount].score = PrintFinalArrangement(solutionCandidates[printedSolutionCount].pieces, ofs).score;
    printedSolutionCount++;
  }

  sort(solutionCandidates.begin(), solutionCandidates.end());

  int mainLoopCounter = 0;
  double startTime = GetNowTime();
  double currentElapsedTime = GetNowTime();
  const double START_TEMP = 0.0;
  const double END_TEMP = 200000.0;
  timeLimit = TL - startTime;
  double currentTemperature = START_TEMP;

  int randomCandidateIndex = 0;
  vector<RectanglePiece> keep = solutionCandidates[randomCandidateIndex].pieces;
  int raQCount = min(candidateLayoutCount, 3);
  int outputIteration = printedSolutionCount;
  while (outputIteration < t) {
    mainLoopCounter++;
    if (mainLoopCounter % 100 == 0) {
      currentElapsedTime = GetNowTime();
    }

    int randomModificationMode = RandXor() % 3;

    if (RandXor() % 2 == 0) {
      solutionCandidates[randomCandidateIndex].pieces = keep;
      randomCandidateIndex = RandXor() % raQCount;
      keep = solutionCandidates[randomCandidateIndex].pieces;
    }

    int randomPieceIndex = RandXor() % n;
    if (randomModificationMode == 0) {
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].rot = 1 - solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].rot;
    }
    else if (randomModificationMode == 1) {
      int newBaseIndex = RandXor() % (randomPieceIndex + 1) - 1;
      if (solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base == newBaseIndex)continue;
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base = newBaseIndex;
    }
    else {
      solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].dir = 1 - solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].dir;
      if (RandXor() % 2 == 0) {
        int newBaseIndex = RandXor() % (randomPieceIndex + 1) - 1;
        solutionCandidates[randomCandidateIndex].pieces[randomPieceIndex].base = newBaseIndex;
      }
    }

    auto preScore = EvaluateScore(solutionCandidates[randomCandidateIndex].pieces, false);

    double progressRatio = (currentElapsedTime - startTime) / timeLimit;
    currentTemperature = START_TEMP + (END_TEMP - START_TEMP) * progressRatio * progressRatio * progressRatio;
    double scoreDifference = solutionCandidates[randomCandidateIndex].score - preScore.score - 100000;
    double acceptanceProbability = exp(scoreDifference / currentTemperature);

    if (preScore.score < solutionCandidates[randomCandidateIndex].score || acceptanceProbability > Rand01() || currentElapsedTime > TL) {
      ScoreStruct score = PrintFinalArrangement(solutionCandidates[randomCandidateIndex].pieces, ofs);
      if (score.score < solutionCandidates[randomCandidateIndex].score) {
        if (executionMode >= 2) {
          cout << randomCandidateIndex << ' ' << randomModificationMode << ' ' << scoreDifference << ' ' << currentTemperature << ' ' << acceptanceProbability << ' ' << currentElapsedTime << endl;
        }
        solutionCandidates[randomCandidateIndex].score = score.score;
        keep = solutionCandidates[randomCandidateIndex].pieces;
      }
      outputIteration++;
    }
  }

  if (executionMode >= 2) {
    cout << "finalOutputLoopCount = " << mainLoopCounter << endl;
  }
}

ll ExecuteSolution(int problem_num)
{
  ResetTime();
  InitializeGlobalState();
  LoadInputData(problem_num);

  ofstream ofs;
  OpenOutputStream(problem_num, ofs);

  SetupBasePieces();
  SetupPieceSizeOrder();

  RefineAndPrintSolutions(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ScoreStruct score = FindBestQueryScore();
  if (executionMode >= 2) {
    cout << "ww = " << score.ww << ", hh = " << score.hh << ", score = " << score.score << endl;
  }
  return score.score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    RandXor();
  }

  executionMode = 2;

  if (executionMode == 0) {
    ExecuteSolution(0);
  }
  else if (executionMode == 3) {
    rep(_, 10)
    {
      ll totalScoreSum = 0;
      srep(i, 0, 100)
      {
        ll score = ExecuteSolution(i);
        totalScoreSum += score;
      }
      cout << totalScoreSum << endl;
    }
  }
  else {
    ll totalScoreSum = 0;
    srep(i, 0, 100)
    {
      ll score = ExecuteSolution(i);
      totalScoreSum += score;
      if (executionMode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << totalScoreSum << ", ";
        cout << "time = " << setw(5) << GetNowTime() << ", ";
        cout << endl;
      }
    }
    cout << "sum = " << totalScoreSum << endl;
  }

  return 0;
}
