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


void FisherYates(int* data, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// ランダムデバイスとメルセンヌ・ツイスタの初期化
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);


const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;


const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char dc[4] = { 'U','L','D','R' };

double TL = 1.8;
int mode;
int mode2;

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

//const int MAX_N = 30;

const int n = 20;
int m;

int a[n + 2][n + 2];
int a2[n + 2][n + 2];
vector<int> posX, posY;

int init_a[n + 2][n + 2];
int init_a2[n + 2][n + 2];
vector<int> init_posX, init_posY;

int holeX[20], holeY[20];

int ansScore;
int ans[11000][2];
int ansCount;
int fallCount;
int crystalCount;
vector<P> order;

int best_ansScore;
int best_ans[11000][2];
int best_ansCount;
int best_fallCount;
int best_crystalCount;
vector<P> best_order;

void CopyToBest()
{
  best_ansScore = ansScore;
  best_ansCount = ansCount;
  best_fallCount = fallCount;
  best_crystalCount = crystalCount;
  for (int i = 0; i < ansCount; ++i)
  {
    best_ans[i][0] = ans[i][0];
    best_ans[i][1] = ans[i][1];
  }
  best_order = order;
}

void CopyToAns()
{
  ansScore = best_ansScore;
  ansCount = best_ansCount;
  fallCount = best_fallCount;
  crystalCount = best_crystalCount;
  for (int i = 0; i < ansCount; ++i)
  {
    ans[i][0] = best_ans[i][0];
    ans[i][1] = best_ans[i][1];
  }
  order = best_order;
}

bool IsNG(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void SetUp()
{
  ansScore = 0;
  for (int i = 0; i < n + 2; ++i)
  {
    for (int j = 0; j < n + 2; ++j)
    {
      a[i][j] = -2;
      a2[i][j] = -1;
    }
  }
  posX.clear();
  posY.clear();
  order.clear();
}

// 入力を受け取る関数
void Input(int problemNum)
{
  std::ostringstream oss;
  if (mode2 == 0) {
    oss << "./inA/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  }
  else if (mode2 == 1) {
    oss << "./inB/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  }
  else {
    oss << "./inC/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  }

  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nn;
    cin >> nn >> m;
    for (int i = 0; i < n; ++i)
    {
      string s;
      cin >> s;
      for (int j = 0; j < n; ++j)
      {
        a[i + 1][j + 1] = -1;
        if (s[j] == '@') {
          a[i + 1][j + 1] = 0;
        }
        else if (s[j] == 'a') {
          a[i + 1][j + 1] = 1;
        }
        else if (s[j] == 'b') {
          a[i + 1][j + 1] = 2;
        }
        else if (s[j] == 'c') {
          a[i + 1][j + 1] = 3;
        }
        else if (s[j] == 'A') {
          a[i + 1][j + 1] = 11;
        }
        else if (s[j] == 'B') {
          a[i + 1][j + 1] = 12;
        }
        else if (s[j] == 'C') {
          a[i + 1][j + 1] = 13;
        }
      }
    }
  }
  else {
    // ファイル入力
    int nn;
    ifs >> nn >> m;
    for (int i = 0; i < n; ++i)
    {
      string s;
      ifs >> s;
      for (int j = 0; j < n; ++j)
      {
        a[i + 1][j + 1] = -1;
        if (s[j] == '@') {
          a[i + 1][j + 1] = 0;
        }
        else if (s[j] == 'a') {
          a[i + 1][j + 1] = 1;
        }
        else if (s[j] == 'b') {
          a[i + 1][j + 1] = 2;
        }
        else if (s[j] == 'c') {
          a[i + 1][j + 1] = 3;
        }
        else if (s[j] == 'A') {
          a[i + 1][j + 1] = 11;
        }
        else if (s[j] == 'B') {
          a[i + 1][j + 1] = 12;
        }
        else if (s[j] == 'C') {
          a[i + 1][j + 1] = 13;
        }
      }
    }
  }

  for (int i = 0; i < n + 2; ++i)
  {
    for (int j = 0; j < n + 2; ++j)
    {
      if (1 <= a[i][j] && a[i][j] <= 3) {
        posX.push_back(i);
        posY.push_back(j);
        a2[i][j] = posX.size() - 1;
      }
    }
  }


  for (int i = 0; i < n + 2; ++i)
  {
    for (int j = 0; j < n + 2; ++j)
    {
      init_a[i][j] = a[i][j];
      init_a2[i][j] = a2[i][j];
    }
  }
  init_posX = posX;
  init_posY = posY;

  for (int i = 0; i < n + 2; ++i)
  {
    for (int j = 0; j < n + 2; ++j)
    {
      if (a[i][j] >= 11) {
        holeX[a[i][j] - 10] = i;
        holeY[a[i][j] - 10] = j;
      }
    }
  }
}

void InitA()
{
  for (int i = 0; i < n + 2; ++i)
  {
    for (int j = 0; j < n + 2; ++j)
    {
      a[i][j] = init_a[i][j];
      a2[i][j] = init_a2[i][j];
    }
  }
  posX = init_posX;
  posY = init_posY;
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
ll CalcScore()
{
  ll res = 0;
  if (fallCount < crystalCount) {
    res = round(1000000.0 * fallCount / crystalCount);
  }
  else {
    res = round(1000000.0 * (1.0 + log2(10000.0 / ansCount)));
  }
  return res;
}

// 解答を出力する関数
void Output(ofstream& ofs)
{
  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < ansCount; ++i)
    {
      cout << ans[i][0] << ' ' << dc[ans[i][1]] << endl;
    }
  }
  else {
    // ファイル出力
    ofs << ansCount << endl;
    for (int i = 0; i < ansCount; ++i)
    {
      ofs << ans[i][0] << ' ' << dc[ans[i][1]] << endl;
    }
  }
}

int queueArr[1000][2];
int queueHead = 0;
int queueTail = 0;
void ClearQueue()
{
  queueHead = 0;
  queueTail = 0;
}
int FrontX()
{
  return queueArr[queueHead][0];
}
int FrontY()
{
  return queueArr[queueHead][1];
}
void Push(int x, int y)
{
  queueArr[queueTail][0] = x;
  queueArr[queueTail][1] = y;
  queueTail++;
}
void Pop()
{
  queueHead++;
}
int Size()
{
  return queueTail - queueHead;
}

// ハイパーパラメータ
struct Hypers
{
  double StartTemp;
  double EndTemp;
  double MultipleValue;
  int Partition;
};

int sousa[1000];
int sousaCount;

int visited[n + 2][n + 2];
int hukugen[n + 2][n + 2];
int bfsCount = 0;

void Method1_3_Fall(int hx, int hy, int& nowX, int& nowY, int sx, int sy, int id, int opDir)
{
  int U = hx;
  int D = hx;
  int L = hy;
  int R = hy;
  while (a[U - 1][hy] == -1 || a2[U - 1][hy] == id)U--;
  while (a[D + 1][hy] == -1 || a2[D + 1][hy] == id)D++;
  while (a[hx][L - 1] == -1 || a2[hx][L - 1] == id)L--;
  while (a[hx][R + 1] == -1 || a2[hx][R + 1] == id)R++;

  // まず移動
  while (nowX > sx) {
    ans[ansCount][0] = 1;
    ans[ansCount][1] = 0;
    ansCount++;
    nowX--;
  }
  while (sx > nowX) {
    ans[ansCount][0] = 1;
    ans[ansCount][1] = 2;
    ansCount++;
    nowX++;
  }
  while (nowY > sy) {
    ans[ansCount][0] = 1;
    ans[ansCount][1] = 1;
    ansCount++;
    nowY--;
  }
  while (sy > nowY) {
    ans[ansCount][0] = 1;
    ans[ansCount][1] = 3;
    ansCount++;
    nowY++;
  }

  if (opDir < 4) {
    ans[ansCount][0] = 3;
    ans[ansCount][1] = opDir;
    ansCount++;
    int x = nowX;
    int y = nowY;
    while (true) {
      int nx = x + dx[opDir];
      int ny = y + dy[opDir];
      if (a[nx][ny] >= 10) {
        if (a[nx][ny] - 10 == a[nowX][nowY]) {
          fallCount++;
        }
        a[nowX][nowY] = -1;
        a2[nowX][nowY] = -1;
        posX[id] = -1;
        posY[id] = -1;
        break;
      }
      else if (a[nx][ny] == -2 || (0 <= a[nx][ny] && a[nx][ny] <= 3)) {
        swap(a[nowX][nowY], a[x][y]);
        swap(a2[nowX][nowY], a2[x][y]);
        posX[id] = x;
        posY[id] = y;
        break;
      }
      else {
        x = nx;
        y = ny;
      }
    }
    return;
  }

  // 鉱石を運ぶbfs
  int saitan = 1;
  int tx = -1;
  int ty = -1;

  if (abs(hx - nowX) <= abs(hy - nowY)) {
    tx = hx;
    ty = nowY;
    if (nowX < hx) {
      for (int i = nowX + 1; i < hx + 1; ++i)
      {
        if (a[i][nowY] != -1) {
          saitan = 0;
          break;
        }
        hukugen[i][nowY] = 0;
      }
    }
    else {
      for (int i = hx; i < nowX; ++i)
      {
        if (a[i][nowY] != -1) {
          saitan = 0;
          break;
        }
        hukugen[i][nowY] = 2;
      }
    }
  }
  else {
    tx = nowX;
    ty = hy;
    if (nowY < hy) {
      for (int i = nowY + 1; i < hy + 1; ++i)
      {
        if (a[nowX][i] != -1) {
          saitan = 0;
          break;
        }
        hukugen[nowX][i] = 1;
      }
    }
    else {
      for (int i = hy; i < nowY; ++i)
      {
        if (a[nowX][i] != -1) {
          saitan = 0;
          break;
        }
        hukugen[nowX][i] = 3;
      }
    }
  }

  if (saitan == 0) {
    tx = -1;
    ty = -1;

    bfsCount++;
    visited[nowX][nowY] = bfsCount;
    ClearQueue();
    Push(nowX, nowY);
    while (Size() > 0) {
      int x = FrontX();
      int y = FrontY();
      Pop();

      if (x == hx && L <= y && y <= R) {
        tx = x;
        ty = y;
        break;
      }
      else if (y == hy && U <= x && x <= D) {
        tx = x;
        ty = y;
        break;
      }

      for (int j = 0; j < 4; ++j)
      {
        int nx = x + dx[j];
        int ny = y + dy[j];
        if (a[nx][ny] != -1)continue;
        if (visited[nx][ny] != bfsCount) {
          visited[nx][ny] = bfsCount;
          Push(nx, ny);
          hukugen[nx][ny] = (j + 2) % 4;
        }
      }
    }
  }

  if (tx == -1) {
    //cerr << "return" << endl;
    return;
  }

  sousaCount = 0;
  int xx = tx;
  int yy = ty;
  while (xx != nowX || yy != nowY) {
    sousa[sousaCount] = (hukugen[xx][yy] + 2) % 4;
    sousaCount++;
    int nx = xx + dx[hukugen[xx][yy]];
    int ny = yy + dy[hukugen[xx][yy]];
    xx = nx;
    yy = ny;
  }
  for (int j = sousaCount - 1; j >= 0; --j)
  {
    ans[ansCount][0] = 2;
    ans[ansCount][1] = sousa[j];
    ansCount++;
    int nx = nowX + dx[sousa[j]];
    int ny = nowY + dy[sousa[j]];
    swap(a[nowX][nowY], a[nx][ny]);
    swap(a2[nowX][nowY], a2[nx][ny]);
    if (id >= 0) {
      posX[id] = nx;
      posY[id] = ny;
    }
    nowX = nx;
    nowY = ny;
  }

  // 穴に転がして入れる
  if (nowX == hx && nowY < hy) {
    ans[ansCount][0] = 3;
    ans[ansCount][1] = 3;
    ansCount++;
    a[nowX][nowY] = -1;
    a2[nowX][nowY] = -1;
    //cerr << 3 << endl;
    if (id >= 0) {
      posX[id] = -1;
      posY[id] = -1;
      fallCount++;
    }
  }
  else  if (nowX == hx && nowY > hy) {
    ans[ansCount][0] = 3;
    ans[ansCount][1] = 1;
    ansCount++;
    a[nowX][nowY] = -1;
    a2[nowX][nowY] = -1;
    //cerr << 1 << endl;
    if (id >= 0) {
      posX[id] = -1;
      posY[id] = -1;
      fallCount++;
    }
  }
  else if (nowY == hy && nowX < hx) {
    ans[ansCount][0] = 3;
    ans[ansCount][1] = 2;
    ansCount++;
    a[nowX][nowY] = -1;
    a2[nowX][nowY] = -1;
    //cerr << 2 << endl;
    if (id >= 0) {
      posX[id] = -1;
      posY[id] = -1;
      fallCount++;
    }
  }
  else  if (nowY == hy && nowX > hx) {
    ans[ansCount][0] = 3;
    ans[ansCount][1] = 0;
    ansCount++;
    a[nowX][nowY] = -1;
    a2[nowX][nowY] = -1;
    //cerr << 0 << endl;
    if (id >= 0) {
      posX[id] = -1;
      posY[id] = -1;
      fallCount++;
    }
  }
  else {
    cerr << "NG" << endl;
  }
}

void Method1_3()
{
  int nowX = holeX[1];
  int nowY = holeY[1];

  int isAllConnected = false;

  for (int i = 0; i < order.size(); ++i)
  {
    // 十字に落とせる岩があれば落とす
    while (true) {
      int tx = -1;
      int ty = -1;
      int dir = -1;
      int isCrystal = 0;
      for (int j = 0; j < m; ++j)
      {
        int hx = holeX[j + 1];
        int hy = holeY[j + 1];
        if (dir == -1) {
          for (int k = hx - 1; k >= 1; --k)
          {
            if (a[k][hy] == 0) {
              tx = k;
              ty = hy;
              dir = 2;
              break;
            }
            else if (a[k][hy] == j + 1) {
              tx = k;
              ty = hy;
              dir = 2;
              isCrystal = 1;
              break;
            }
            else if (a[k][hy] != -1) {
              break;
            }
          }
        }
        if (dir == -1) {
          for (int k = hx + 1; k < n + 1; ++k)
          {
            if (a[k][hy] == 0) {
              tx = k;
              ty = hy;
              dir = 0;
              break;
            }
            else if (a[k][hy] == j + 1) {
              tx = k;
              ty = hy;
              dir = 0;
              isCrystal = 1;
              break;
            }
            else if (a[k][hy] != -1) {
              break;
            }
          }
        }
        if (dir == -1) {
          for (int k = hy - 1; k >= 1; --k)
          {
            if (a[hx][k] == 0) {
              tx = hx;
              ty = k;
              dir = 3;
              break;
            }
            else if (a[hx][k] == j + 1) {
              tx = hx;
              ty = k;
              dir = 3;
              isCrystal = 1;
              break;
            }
            else if (a[hx][k] != -1) {
              break;
            }
          }
        }
        if (dir == -1) {
          for (int k = hy + 1; k < n + 1; ++k)
          {
            if (a[hx][k] == 0) {
              tx = hx;
              ty = k;
              dir = 1;
              break;
            }
            else if (a[hx][k] == j + 1) {
              tx = hx;
              ty = k;
              dir = 1;
              isCrystal = 1;
              break;
            }
            else if (a[hx][k] != -1) {
              break;
            }
          }
        }
      }

      if (dir == -1)break;

      // まず移動
      while (nowX > tx) {
        ans[ansCount][0] = 1;
        ans[ansCount][1] = 0;
        ansCount++;
        nowX--;
      }
      while (tx > nowX) {
        ans[ansCount][0] = 1;
        ans[ansCount][1] = 2;
        ansCount++;
        nowX++;
      }
      while (nowY > ty) {
        ans[ansCount][0] = 1;
        ans[ansCount][1] = 1;
        ansCount++;
        nowY--;
      }
      while (ty > nowY) {
        ans[ansCount][0] = 1;
        ans[ansCount][1] = 3;
        ansCount++;
        nowY++;
      }

      // 転がして落とす
      if (isCrystal == 1) {
        int id = a2[nowX][nowY];
        posX[id] = -1;
        posY[id] = -1;
        fallCount++;
      }
      ans[ansCount][0] = 3;
      ans[ansCount][1] = dir;
      ansCount++;
      a[nowX][nowY] = -1;
      a2[nowX][nowY] = -1;
    }

    // C対策：全部の部屋を繋げる
    if (!isAllConnected && mode2 == 2) {
      while (true) {
        int hx = holeX[1];
        int hy = holeY[1];
        // 今の数
        int nowCount = 0;
        {
          bfsCount++;
          ClearQueue();
          Push(hx, hy);
          visited[hx][hy] = bfsCount;
          while (Size() > 0) {
            int x = FrontX();
            int y = FrontY();
            Pop();
            for (int j = 0; j < 4; ++j)
            {
              int nx = x + dx[j];
              int ny = y + dy[j];
              if (a[nx][ny] == -2 || a[nx][ny] == 0)continue;
              if (visited[nx][ny] != bfsCount) {
                visited[nx][ny] = bfsCount;
                Push(nx, ny);
                if (a[nx][ny] == 1)nowCount++;
              }
            }
          }
        }
        if (nowCount == crystalCount - fallCount) {
          isAllConnected = true;
          break;
        }

        // どこか一列空ける
        int mi = 1001001;
        int dir = -1;
        int idx = -1;
        int maxCount = 0;
        // 行
        for (int k = 1; k < n + 1; ++k)
        {
          int tmpCount = 0;
          int rockCount = 0;
          bfsCount++;
          ClearQueue();
          Push(hx, hy);
          visited[hx][hy] = bfsCount;
          while (Size() > 0) {
            int x = FrontX();
            int y = FrontY();
            Pop();
            for (int j = 0; j < 4; ++j)
            {
              int nx = x + dx[j];
              int ny = y + dy[j];
              if (a[nx][ny] == -2)continue;
              if (a[nx][ny] == 0 && nx != k)continue;
              if (visited[nx][ny] != bfsCount) {
                visited[nx][ny] = bfsCount;
                Push(nx, ny);
                if (a[nx][ny] == 1)tmpCount++;
                if (a[nx][ny] == 0)rockCount++;
              }
            }
          }
          if (tmpCount > nowCount && rockCount < mi) {
            mi = rockCount;
            dir = 0;
            idx = k;
            maxCount = tmpCount;
          }
        }
        // 列
        for (int k = 1; k < n + 1; ++k)
        {
          int tmpCount = 0;
          int rockCount = 0;
          bfsCount++;
          ClearQueue();
          Push(hx, hy);
          visited[hx][hy] = bfsCount;
          while (Size() > 0) {
            int x = FrontX();
            int y = FrontY();
            Pop();
            for (int j = 0; j < 4; ++j)
            {
              int nx = x + dx[j];
              int ny = y + dy[j];
              if (a[nx][ny] == -2)continue;
              if (a[nx][ny] == 0 && ny != k)continue;
              if (visited[nx][ny] != bfsCount) {
                visited[nx][ny] = bfsCount;
                Push(nx, ny);
                if (a[nx][ny] == 1)tmpCount++;
                if (a[nx][ny] == 0)rockCount++;
              }
            }
          }
          if (tmpCount > nowCount && rockCount < mi) {
            mi = rockCount;
            dir = 1;
            idx = k;
            maxCount = tmpCount;
          }
        }

        if (dir == -1) {
          cerr << "NGROCK" << endl;
          break;
        }

        //cout << nowCount << ' ' << maxCount << ' ' << mi << ' ' << dir << ' ' << idx << endl;

        // その行・列を全部落とす
        if (dir == 0) {
          for (int k = hy - 1; k >= 1; --k)
          {
            int huyou = 1;
            if (k == 1)huyou = 0;
            for (int l = k - 1; l >= 1; --l)
            {
              for (int z = 0; z < 4; ++z)
              {
                int nx = idx + dx[z];
                int ny = l + dy[z];
                if (a[nx][ny] != 0 && a[nx][ny] != -2) {
                  huyou = 0;
                }
              }
              if (huyou == 0)break;
            }
            if (huyou == 1) {
              break;
            }
            if (a[idx][k] == 0) {
              Method1_3_Fall(hx, hy, nowX, nowY, idx, k, -999, 5);
            }
            else if (a[idx][k] == 1) {
              Method1_3_Fall(hx, hy, nowX, nowY, idx, k, a2[idx][k], 5);
            }
          }
          for (int k = hy + 1; k < n + 1; ++k)
          {
            int huyou = 1;
            if (k == n)huyou = 0;
            for (int l = k + 1; l < n + 1; ++l)
            {
              for (int z = 0; z < 4; ++z)
              {
                int nx = idx + dx[z];
                int ny = l + dy[z];
                if (a[nx][ny] != 0 && a[nx][ny] != -2) {
                  huyou = 0;
                }
              }
              if (huyou == 0)break;
            }
            if (huyou == 1) {
              break;
            }
            if (a[idx][k] == 0) {
              Method1_3_Fall(hx, hy, nowX, nowY, idx, k, -999, 5);
            }
            else if (a[idx][k] == 1) {
              Method1_3_Fall(hx, hy, nowX, nowY, idx, k, a2[idx][k], 5);
            }
          }
        }
        else {
          for (int k = hx - 1; k >= 1; --k)
          {
            int huyou = 1;
            if (k == 1)huyou = 0;
            for (int l = k - 1; l >= 1; --l)
            {
              for (int z = 0; z < 4; ++z)
              {
                int nx = l + dx[z];
                int ny = idx + dy[z];
                if (a[nx][ny] != 0 && a[nx][ny] != -2) {
                  huyou = 0;
                }
              }
              if (huyou == 0)break;
            }
            if (huyou == 1) {
              break;
            }
            if (a[k][idx] == 0) {
              Method1_3_Fall(hx, hy, nowX, nowY, k, idx, -999, 5);
            }
            else if (a[k][idx] == 1) {
              Method1_3_Fall(hx, hy, nowX, nowY, k, idx, a2[k][idx], 5);
            }
          }
          for (int k = hx + 1; k < n + 1; ++k)
          {
            int huyou = 1;
            if (k == n)huyou = 0;
            for (int l = k + 1; l < n + 1; ++l)
            {
              for (int z = 0; z < 4; ++z)
              {
                int nx = l + dx[z];
                int ny = idx + dy[z];
                if (a[nx][ny] != 0 && a[nx][ny] != -2) {
                  huyou = 0;
                }
              }
              if (huyou == 0)break;
            }
            if (huyou == 1) {
              break;
            }
            if (a[k][idx] == 0) {
              Method1_3_Fall(hx, hy, nowX, nowY, k, idx, -999, 5);
            }
            else if (a[k][idx] == 1) {
              Method1_3_Fall(hx, hy, nowX, nowY, k, idx, a2[k][idx], 5);
            }
          }
        }
      }
    }

    int id = order[i].first;
    if (posX[id] == -1)continue;
    int num = a[posX[id]][posY[id]];
    int hx = holeX[num];
    int hy = holeY[num];
    Method1_3_Fall(hx, hy, nowX, nowY, posX[id], posY[id], id, order[i].second);
  }
}

int ErasePos[1000000];

// ナイーブな解法
void Method1(Hypers hypers)
{
  ansCount = 0;
  crystalCount = posX.size();
  for (int i = 0; i < crystalCount; ++i)order.push_back(P(i, 5));
  fallCount = 0;
  ansScore = 0;
  CopyToBest();


  double nowTime = GetNowTime();
  const double START_TEMP = hypers.StartTemp;
  const double END_TEMP = hypers.EndTemp;

  int loop = 0;
  vector<P> keepOrder;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL) { break; }
    }

    // 戻す
    //if (ansScore * 1.2 < best_ansScore) {
    //  CopyToAns();
    //}

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

    keepOrder = order;
    int raMode = Rand() % 100;
    if (raMode < hypers.Partition) {
      // 近傍解作成
      int ra1 = Rand() % crystalCount;
      int ra2 = Rand() % crystalCount;

      swap(order[ra1], order[ra2]);
    }
    else if (raMode < 75) {
      int ra1 = Rand() % crystalCount;
      int raPos = Rand() % order.size();
      int raDir = Rand() % 4;
      order.insert(order.begin() + raPos, P(ra1, raDir));
    }
    else if (raMode < 100) {
      int eCount = 0;
      for (int i = 0; i < order.size(); ++i)
      {
        if (order[i].second < 4) {
          ErasePos[eCount] = i;
          eCount++;
        }
      }
      if (eCount == 0)continue;
      int raIdx = Rand() % eCount;
      int pos = ErasePos[raIdx];
      order.erase(order.begin() + pos);
    }

    ansCount = 0;
    fallCount = 0;
    InitA();

    Method1_3();

    // スコア計算
    double tmpScore = CalcScore();

    // 焼きなまし
    double diffScore = (tmpScore - ansScore) * hypers.MultipleValue;
    double prob = exp(diffScore / temp);
    if (prob > Rand01()) {
      // 採用
      ansScore = tmpScore;

      // Best解よりもいいか
      if (ansScore > best_ansScore) {
        CopyToBest();
      }
    }
    else {
      // 元に戻す
      //swap(order[ra1], order[ra2]);
      order = keepOrder;
    }
  }


  CopyToAns();

  if (mode != 0 && mode != 3) {
    cout << loop << endl;
    for (int i = 0; i < order.size(); ++i)
    {
      if (order[i].second < 4) {
        cout << order[i].first << ':' << order[i].second << endl;
      }
    }
  }
}

// 問題を解く関数
ll Solve(int problem_num, Hypers hypers)
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
  Method1(hypers);

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

int main()
{
  mode = 2;
  mode2 = 1;

  Hypers HYPERS;
  HYPERS.StartTemp = 2000.0;
  HYPERS.EndTemp = 0.0;
  HYPERS.MultipleValue = 1.0;
  HYPERS.Partition = 50;

  if (mode == 0) {
    Solve(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i)
    {
      ll score = Solve(i, HYPERS);
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
  else if (mode == 3) {
    int loop = 0;
    Hypers bestHypers;
    ll bestSumScore = 0;

    while (true) {
      Hypers hypers;
      hypers.StartTemp = pow(2.0, Rand01() * 20);
      hypers.EndTemp = 0.0;
      hypers.MultipleValue = pow(2.0, Rand01() * 20);
      hypers.Partition = Rand() % 101;

      ll sum = 0;
      for (int i = 8; i < 15; ++i)
      {
        ll score = Solve(i, hypers);
        sum += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout
        << "Loop = " << loop
        << ", Sum = " << sum
        << ", StartTemp = " << hypers.StartTemp
        << ", EndTemp = " << hypers.EndTemp
        << ", MultipleValue = " << hypers.MultipleValue
        << ", Partition1 = " << hypers.Partition
        << endl;

      if (sum > bestSumScore) {
        bestSumScore = sum;
        bestHypers = hypers;
      }

      loop++;
    }
  }

  return 0;
}
