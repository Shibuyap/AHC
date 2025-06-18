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

// Removed macro definitions - using standard for loops instead

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

std::chrono::steady_clock::time_point start_time_clock;

void start_timer()
{
  start_time_clock = std::chrono::steady_clock::now();
}

double get_elapsed_time()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
  return elapsed.count();
}

static uint32_t rand_xorshift()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  return w;
}

static double rand_01()
{
  return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
}

static double rand_range(double l, double r)
{
  return l + (r - l) * rand_01();
}

static uint32_t rand_range(uint32_t l, uint32_t r)
{
  return l + rand_xorshift() % (r - l + 1); // [l, r]
}

void shuffle_array(int* arr, int n)
{
  for (int i = n - 1; i >= 0; i--) {
    int j = rand_xorshift() % (i + 1);
    int swa = arr[i];
    arr[i] = arr[j];
    arr[j] = swa;
  }
}

std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

struct Point { int x, y; };

double get_distance(Point p1, Point p2)
{
  return sqrt((ll)(p1.x - p2.x) * (p1.x - p2.x) + (ll)(p1.y - p2.y) * (p1.y - p2.y));
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

const double TIME_LIMIT = 1.9;
int exec_mode;

const int H = 1000000;
const int W = 1000000;
const int MAX_N = 100;
const int MAX_T = 10000;
int A, B, C;
Point a[MAX_N], b[MAX_N], c[MAX_N];
vector<Point> ua, ub, uc;
vector<Point> da, db, dc;

int a_flag[MAX_N];
int b_flag[MAX_N];

void reset_flag()
{
  for (int i = 0; i < A; ++i) a_flag[i] = 0;
  for (int i = 0; i < B; ++i) b_flag[i] = 0;
}

Point box_a[MAX_N][2];
int box_a_count[MAX_N];
int box_a_flag[MAX_N][MAX_N];
Point box_b[MAX_N][2];
int box_b_count[MAX_N];
int box_b_flag[MAX_N][MAX_N];

double current_score;
int moves_a_count;
int moves_b_count;
Point moves_a[MAX_T][2];
Point moves_b[MAX_T][2];
int nums_a[MAX_N];
int nums_b[MAX_N];

double keep_score;
int keep_moves_a_count;
int keep_moves_b_count;
Point keep_moves_a[MAX_T][2];
Point keep_moves_b[MAX_T][2];
int keep_nums_a[MAX_N];
int keep_nums_b[MAX_N];

double best_score;
int best_moves_a_count;
int best_moves_b_count;
Point best_moves_a[MAX_T][2];
Point best_moves_b[MAX_T][2];
int best_nums_a[MAX_N];
int best_nums_b[MAX_N];

void store_keep_score()
{
  keep_score = current_score;
  keep_moves_a_count = moves_a_count;
  keep_moves_b_count = moves_b_count;
  for (int i = 0; i < moves_a_count; ++i) {
    keep_moves_a[i][0] = moves_a[i][0];
    keep_moves_a[i][1] = moves_a[i][1];
  }
  for (int i = 0; i < moves_b_count; ++i) {
    keep_moves_b[i][0] = moves_b[i][0];
    keep_moves_b[i][1] = moves_b[i][1];
  }
  for (int i = 0; i < A; ++i) {
    keep_nums_a[i] = nums_a[i];
  }
  for (int i = 0; i < B; ++i) {
    keep_nums_b[i] = nums_b[i];
  }

}

void restore_keep_score()
{
  current_score = keep_score;
  moves_a_count = keep_moves_a_count;
  moves_b_count = keep_moves_b_count;
  for (int i = 0; i < moves_a_count; ++i) {
    moves_a[i][0] = keep_moves_a[i][0];
    moves_a[i][1] = keep_moves_a[i][1];
  }
  for (int i = 0; i < moves_b_count; ++i) {
    moves_b[i][0] = keep_moves_b[i][0];
    moves_b[i][1] = keep_moves_b[i][1];
  }
  for (int i = 0; i < A; ++i) {
    nums_a[i] = keep_nums_a[i];
  }
  for (int i = 0; i < B; ++i) {
    nums_b[i] = keep_nums_b[i];
  }

}

void store_best_score()
{
  best_score = current_score;
  best_moves_a_count = moves_a_count;
  best_moves_b_count = moves_b_count;
  for (int i = 0; i < moves_a_count; ++i) {
    best_moves_a[i][0] = moves_a[i][0];
    best_moves_a[i][1] = moves_a[i][1];
  }
  for (int i = 0; i < moves_b_count; ++i) {
    best_moves_b[i][0] = moves_b[i][0];
    best_moves_b[i][1] = moves_b[i][1];
  }
  for (int i = 0; i < A; ++i) {
    best_nums_a[i] = nums_a[i];
  }
  for (int i = 0; i < B; ++i) {
    best_nums_b[i] = nums_b[i];
  }

}

void restore_best_score()
{
  current_score = best_score;
  moves_a_count = best_moves_a_count;
  moves_b_count = best_moves_b_count;
  for (int i = 0; i < moves_a_count; ++i) {
    moves_a[i][0] = best_moves_a[i][0];
    moves_a[i][1] = best_moves_a[i][1];
  }
  for (int i = 0; i < moves_b_count; ++i) {
    moves_b[i][0] = best_moves_b[i][0];
    moves_b[i][1] = best_moves_b[i][1];
  }
  for (int i = 0; i < A; ++i) {
    nums_a[i] = best_nums_a[i];
  }
  for (int i = 0; i < B; ++i) {
    nums_b[i] = best_nums_b[i];
  }

}

bool is_out_of_range(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

void initialize_state()
{
  current_score = 0;
  moves_a_count = 0;
  moves_b_count = 0;
  store_best_score();
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./inC/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> A >> B >> C;
    for (int i = 0; i < A; ++i) {
      cin >> a[i].x >> a[i].y;
    }
    for (int i = 0; i < B; ++i) {
      cin >> b[i].x >> b[i].y;
    }
    for (int i = 0; i < C; ++i) {
      cin >> c[i].x >> c[i].y;
    }
  }
  else {
    // ファイル入力
    ifs >> A >> B >> C;
    for (int i = 0; i < A; ++i) {
      ifs >> a[i].x >> a[i].y;
    }
    for (int i = 0; i < B; ++i) {
      ifs >> b[i].x >> b[i].y;
    }
    for (int i = 0; i < C; ++i) {
      ifs >> c[i].x >> c[i].y;
    }
  }
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

double calculate_score()
{
  double res = 0;
  for (int i = 1; i < max(moves_a_count, moves_b_count); ++i) {
    double d = 0;
    if (i < moves_a_count) {
      d = max(d, get_distance(moves_a[i][0], moves_a[i - 1][0]) + get_distance(moves_a[i][1], moves_a[i - 1][1]));
    }
    if (i < moves_b_count) {
      d = max(d, get_distance(moves_b[i][0], moves_b[i - 1][0]) + get_distance(moves_b[i][1], moves_b[i - 1][1]));
    }
    res += d;
  }
  return res;
}



double simulate_score()
{
  reset_flag();

  moves_a_count = 0;

  int x = a[nums_a[0]].x;
  int y = a[nums_a[0]].y;
  for (int i = 0; i < A; ++i) {
    int idx = nums_a[i];
    if (a_flag[idx] == 1) {
      //cout << "OK" << endl;
      continue;
    }
    int ux = box_a[idx][0].x;
    int dx = box_a[idx][1].x;
    int uy = box_a[idx][0].y;
    int dy = box_a[idx][1].y;
    moves_a[moves_a_count][0].x = ux;
    moves_a[moves_a_count][0].y = uy;
    moves_a[moves_a_count][1].x = ux;
    moves_a[moves_a_count][1].y = uy;
    moves_a_count++;
    moves_a[moves_a_count][0].x = ux;
    moves_a[moves_a_count][0].y = dy;
    moves_a[moves_a_count][1].x = dx;
    moves_a[moves_a_count][1].y = uy;
    moves_a_count++;
    moves_a[moves_a_count][0].x = dx;
    moves_a[moves_a_count][0].y = dy;
    moves_a[moves_a_count][1].x = dx;
    moves_a[moves_a_count][1].y = dy;
    moves_a_count++;
    for (int j = 0; j < A; ++j) {
      if (box_a_flag[idx][j] == 1) {
        a_flag[j] = 1;
      }
    }
  }
  //int xxx;
  //cin >> xxx;
  moves_b_count = 0;

  x = b[nums_b[0]].x;
  y = b[nums_b[0]].y;
  for (int i = 0; i < B; ++i) {
    int idx = nums_b[i];
    if (b_flag[idx] == 1) {
      continue;
    }
    int ux = box_b[idx][0].x;
    int dx = box_b[idx][1].x;
    int uy = box_b[idx][0].y;
    int dy = box_b[idx][1].y;
    moves_b[moves_b_count][0].x = ux;
    moves_b[moves_b_count][0].y = uy;
    moves_b[moves_b_count][1].x = ux;
    moves_b[moves_b_count][1].y = uy;
    moves_b_count++;
    moves_b[moves_b_count][0].x = ux;
    moves_b[moves_b_count][0].y = dy;
    moves_b[moves_b_count][1].x = dx;
    moves_b[moves_b_count][1].y = uy;
    moves_b_count++;
    moves_b[moves_b_count][0].x = dx;
    moves_b[moves_b_count][0].y = dy;
    moves_b[moves_b_count][1].x = dx;
    moves_b[moves_b_count][1].y = dy;
    moves_b_count++;
    for (int j = 0; j < B; ++j) {
      if (box_b_flag[idx][j] == 1) {
        b_flag[j] = 1;
      }
    }
  }

  double score = calculate_score();
  score = moves_a_count + moves_b_count;
  return score;
}

double calculate_score_one_point(int idx)
{
  double res = 0;
  for (int i = idx; i < idx + 2; ++i) {
    double d = 0;
    if (i < moves_a_count) {
      d = max(d, get_distance(moves_a[i][0], moves_a[i - 1][0]) + get_distance(moves_a[i][1], moves_a[i - 1][1]));
    }
    if (i < moves_b_count) {
      d = max(d, get_distance(moves_b[i][0], moves_b[i - 1][0]) + get_distance(moves_b[i][1], moves_b[i - 1][1]));
    }
    res += d;
  }
  return res;
}

void output_data(ofstream& ofs)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < max(moves_a_count, moves_b_count); ++i) {
      if (i < moves_a_count) {
        cout << moves_a[i][0].x << " " << moves_a[i][0].y << " " << moves_a[i][1].x << " " << moves_a[i][1].y << ' ';
      }
      else {
        cout << moves_a[moves_a_count - 1][0].x << " " << moves_a[moves_a_count - 1][0].y << " " << moves_a[moves_a_count - 1][1].x << " " << moves_a[moves_a_count - 1][1].y << ' ';
      }
      if (i < moves_b_count) {
        cout << moves_b[i][0].x << " " << moves_b[i][0].y << " " << moves_b[i][1].x << " " << moves_b[i][1].y << endl;
      }
      else {
        cout << moves_b[moves_b_count - 1][0].x << " " << moves_b[moves_b_count - 1][0].y << " " << moves_b[moves_b_count - 1][1].x << " " << moves_b[moves_b_count - 1][1].y << endl;
      }
    }
  }
  else {
    // ファイル出力
    for (int i = 0; i < max(moves_a_count, moves_b_count); ++i) {
      if (i < moves_a_count) {
        ofs << moves_a[i][0].x << " " << moves_a[i][0].y << " " << moves_a[i][1].x << " " << moves_a[i][1].y << ' ';
      }
      else {
        ofs << moves_a[moves_a_count - 1][0].x << " " << moves_a[moves_a_count - 1][0].y << " " << moves_a[moves_a_count - 1][1].x << " " << moves_a[moves_a_count - 1][1].y << ' ';
      }
      if (i < moves_b_count) {
        ofs << moves_b[i][0].x << " " << moves_b[i][0].y << " " << moves_b[i][1].x << " " << moves_b[i][1].y << endl;
      }
      else {
        ofs << moves_b[moves_b_count - 1][0].x << " " << moves_b[moves_b_count - 1][0].y << " " << moves_b[moves_b_count - 1][1].x << " " << moves_b[moves_b_count - 1][1].y << endl;
      }
    }
  }
}

void init_box()
{
  for (int i = 0; i < A; ++i) {
    box_a[i][0] = a[i];
    box_a[i][1] = a[i];
    box_a_count[i] = 1;
    for (int j = 0; j < A; ++j) {
      box_a_flag[i][j] = 0;
    }
    box_a_flag[i][i] = 1;
  }
  int loop = 0;

  int flag[MAX_N];
  while (true) {
    loop++;
    if (loop % 100 == 0) {
      if (get_elapsed_time() > 0.5) break;
    }

    int up = a[rand_xorshift() % A].x;
    int down = a[rand_xorshift() % A].x;
    int left = a[rand_xorshift() % A].y;
    int right = a[rand_xorshift() % A].y;
    if (up > down) std::swap(up, down);
    if (left > right) std::swap(left, right);

    int ok = 1;
    for (int i = 0; i < B; ++i) {
      if (up <= b[i].x && b[i].x <= down && left <= b[i].y && b[i].y <= right) {
        ok = 0;
        break;
      }
    }
    if (ok == 0) { continue; }
    for (int i = 0; i < C; ++i) {
      if (up <= c[i].x && c[i].x <= down && left <= c[i].y && c[i].y <= right) {
        ok = 0;
        break;
      }
    }
    if (ok == 0) { continue; }
    int cnt = 0;
    for (int i = 0; i < A; ++i) {
      if (up <= a[i].x && a[i].x <= down && left <= a[i].y && a[i].y <= right) {
        flag[cnt] = i;
        cnt++;
      }
    }
    if (cnt == 0) { continue; }
    for (int i = 0; i < cnt; ++i) {
      int idx = flag[i];
      if (cnt > box_a_count[idx]) {
        box_a_count[idx] = cnt;
        box_a[idx][0].x = up;
        box_a[idx][0].y = left;
        box_a[idx][1].x = down;
        box_a[idx][1].y = right;
      }
      for (int j = 0; j < A; ++j) {
        box_a_flag[idx][j] = 0;
      }
      for (int j = 0; j < cnt; ++j) {
        int idx2 = flag[j];
        box_a_flag[idx][idx2] = 1;
      }
    }
  }

  //for (int i = 0; i < (A); ++i) {
  //  cout << box_a_count[i] << " ";
  //}
  //cout << endl;

  if (B == 0) {
    return;
  }

  for (int i = 0; i < B; ++i) {
    box_b[i][0] = b[i];
    box_b[i][1] = b[i];
    box_b_count[i] = 1;
    for (int j = 0; j < B; ++j) {
      box_b_flag[i][j] = 0;
    }
    box_b_flag[i][i] = 1;
  }

  loop = 0;
  while (true) {
    loop++;
    if (loop % 100 == 0) {
      if (get_elapsed_time() > 1.0) break;
    }

    int up = a[rand_xorshift() % B].x;
    int down = a[rand_xorshift() % B].x;
    int left = a[rand_xorshift() % B].y;
    int right = a[rand_xorshift() % B].y;
    if (up > down) std::swap(up, down);
    if (left > right) std::swap(left, right);

    int ok = 1;
    for (int i = 0; i < A; ++i) {
      if (up <= a[i].x && a[i].x <= down && left <= a[i].y && a[i].y <= right) {
        ok = 0;
        break;
      }
    }
    if (ok == 0) { continue; }
    for (int i = 0; i < C; ++i) {
      if (up <= c[i].x && c[i].x <= down && left <= c[i].y && c[i].y <= right) {
        ok = 0;
        break;
      }
    }
    if (ok == 0) { continue; }
    int cnt = 0;
    for (int i = 0; i < B; ++i) {
      if (up <= b[i].x && b[i].x <= down && left <= b[i].y && b[i].y <= right) {
        cnt++;
        flag[cnt] = i;
      }
    }
    if (cnt == 0) { continue; }
    for (int i = 0; i < cnt; ++i) {
      int idx = flag[i];
      if (cnt > box_b_count[idx]) {
        box_b_count[idx] = cnt;
        box_b[idx][0].x = up;
        box_b[idx][0].y = left;
        box_b[idx][1].x = down;
        box_b[idx][1].y = right;
      }
      for (int j = 0; j < B; ++j) {
        box_b_flag[idx][j] = 0;
      }
      for (int j = 0; j < cnt; ++j) {
        int idx2 = flag[j];
        box_b_flag[idx][idx2] = 1;
      }
    }
  }
}

void merge_gomi()
{
  ua.clear();
  da.clear();
  for (int i = 0; i < A; ++i) {
    ua.push_back(a[i]);
    da.push_back(a[i]);
  }


  while (true) {
    ll mi = INF;
    int idx1 = -1;
    int idx2 = -1;
    Point mi_nua, mi_nda;
    for (int i = 0; i < ua.size(); ++i) {
      for (int j = i + 1; j < ua.size(); ++j) {
        if (i == j) { continue; }
        Point nua = ua[i];
        nua.x = min(nua.x, ua[j].x);
        nua.y = min(nua.y, ua[j].y);
        Point nda = da[i];
        nda.x = max(nda.x, da[j].x);
        nda.y = max(nda.y, da[j].y);
        ll tmp = (nda.x - nua.x) * (nda.y - nua.y);
        if (tmp >= mi) { continue; }
        int ok = 1;
        for (int k = 0; k < B; ++k) {
          if (nua.x <= b[k].x && b[k].x <= nda.x && nua.y <= b[k].y && b[k].y <= nda.y) {
            ok = 0;
            break;
          }
        }
        if (ok == 0) { continue; }
        for (int k = 0; k < C; ++k) {
          if (nua.x <= c[k].x && c[k].x <= nda.x && nua.y <= c[k].y && c[k].y <= nda.y) {
            ok = 0;
            break;
          }
        }
        if (ok == 0) { continue; }

        mi = tmp;
        idx1 = i;
        idx2 = j;
        mi_nua = nua;
        mi_nda = nda;
      }
    }

    if (idx1 == -1) {
      break;
    }

    ua.erase(ua.begin() + idx2);
    ua.erase(ua.begin() + idx1);
    da.erase(da.begin() + idx2);
    da.erase(da.begin() + idx1);
    ua.push_back(mi_nua);
    da.push_back(mi_nda);
  }

  ub.clear();
  db.clear();
  for (int i = 0; i < B; ++i) {
    ub.push_back(b[i]);
    db.push_back(b[i]);
  }


  while (true) {
    ll mi = INF;
    int idx1 = -1;
    int idx2 = -1;
    Point mi_nub, mi_ndb;
    for (int i = 0; i < ub.size(); ++i) {
      for (int j = i + 1; j < ub.size(); ++j) {
        if (i == j) { continue; }
        Point nub = ub[i];
        nub.x = min(nub.x, ub[j].x);
        nub.y = min(nub.y, ub[j].y);
        Point ndb = db[i];
        ndb.x = max(ndb.x, db[j].x);
        ndb.y = max(ndb.y, db[j].y);
        ll tmp = (ndb.x - nub.x) * (ndb.y - nub.y);
        if (tmp >= mi) { continue; }
        int ok = 1;
        for (int k = 0; k < A; ++k) {
          if (nub.x <= a[k].x && a[k].x <= ndb.x && nub.y <= a[k].y && a[k].y <= ndb.y) {
            ok = 0;
            break;
          }
        }
        if (ok == 0) { continue; }
        for (int k = 0; k < C; ++k) {
          if (nub.x <= c[k].x && c[k].x <= ndb.x && nub.y <= c[k].y && c[k].y <= ndb.y) {
            ok = 0;
            break;
          }
        }
        if (ok == 0) { continue; }

        mi = tmp;
        idx1 = i;
        idx2 = j;
        mi_nub = nub;
        mi_ndb = ndb;
      }
    }

    if (idx1 == -1) {
      break;
    }

    ub.erase(ub.begin() + idx2);
    ub.erase(ub.begin() + idx1);
    db.erase(db.begin() + idx2);
    db.erase(db.begin() + idx1);
    ub.push_back(mi_nub);
    db.push_back(mi_ndb);
  }
}

void build_initial_solution(double ttll)
{
  moves_a_count = 0;
  for (int i = 0; i < A; ++i) {
    moves_a[moves_a_count][0] = a[i];
    moves_a[moves_a_count][1] = a[i];
    nums_a[i] = i;
    moves_a_count++;
  }
  for (int i = 0; i < B; ++i) {
    moves_b[moves_b_count][0] = b[i];
    moves_b[moves_b_count][1] = b[i];
    nums_b[i] = i;
    moves_b_count++;
  }
  if (moves_b_count == 0) {
    nums_b[0] = 0;
    moves_b[moves_b_count][0].x = 0;
    moves_b[moves_b_count][0].y = 0;
    moves_b[moves_b_count][1].x = 0;
    moves_b[moves_b_count][1].y = 0;
  }

  current_score = calculate_score();
  store_best_score();
}

void build_initial_solution2()
{
  moves_a_count = 0;
  for (int i = 0; i < ua.size(); ++i) {
    nums_a[i] = i;
    moves_a[moves_a_count][0] = ua[i];
    moves_a[moves_a_count][1] = ua[i];
    moves_a_count++;
    moves_a[moves_a_count][0].x = ua[i].x;
    moves_a[moves_a_count][0].y = da[i].y;
    moves_a[moves_a_count][1].x = da[i].x;
    moves_a[moves_a_count][1].y = ua[i].y;
    moves_a_count++;
    moves_a[moves_a_count][0].x = da[i].x;
    moves_a[moves_a_count][0].y = da[i].y;
    moves_a[moves_a_count][1].x = da[i].x;
    moves_a[moves_a_count][1].y = da[i].y;
    moves_a_count++;
  }
  moves_b_count = 0;
  for (int i = 0; i < ub.size(); ++i) {
    nums_b[i] = i;
    moves_b[moves_b_count][0] = ub[i];
    moves_b[moves_b_count][1] = ub[i];
    moves_b_count++;
    moves_b[moves_b_count][0].x = ub[i].x;
    moves_b[moves_b_count][0].y = db[i].y;
    moves_b[moves_b_count][1].x = db[i].x;
    moves_b[moves_b_count][1].y = ub[i].y;
    moves_b_count++;
    moves_b[moves_b_count][0].x = db[i].x;
    moves_b[moves_b_count][0].y = db[i].y;
    moves_b[moves_b_count][1].x = db[i].x;
    moves_b[moves_b_count][1].y = db[i].y;
    moves_b_count++;
  }
  if (moves_b_count == 0) {
    nums_b[0] = 0;
    moves_b[moves_b_count][0].x = 0;
    moves_b[moves_b_count][0].y = 0;
    moves_b[moves_b_count][1].x = 0;
    moves_b[moves_b_count][1].y = 0;
  }

  current_score = calculate_score();
  store_best_score();
}


struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams, double tl)
{
  current_score = calculate_score();
  store_best_score();
  //cout << current_score << endl;

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > tl) break;
    }

    if (rand_xorshift() % 321212 == 0) {
      restore_best_score();
      continue;
    }

    double progress_ratio = now_time / tl;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    double tmp_score = current_score;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % 2;
      if (ra1 == 0) {
        ra2 = rand_xorshift() % (A - 2) + 1;
        ra3 = rand_xorshift() % (A - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        tmp_score -= calculate_score_one_point(ra2);
        tmp_score -= calculate_score_one_point(ra3);
        swap(moves_a[ra2][0], moves_a[ra3][0]);
        swap(moves_a[ra2][1], moves_a[ra3][1]);
        tmp_score += calculate_score_one_point(ra2);
        tmp_score += calculate_score_one_point(ra3);
        //tmp_score = calculate_score();
      }
      else {
        if (B == 0) {
          continue;
        }
        ra2 = rand_xorshift() % (B - 2) + 1;
        ra3 = rand_xorshift() % (B - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        tmp_score -= calculate_score_one_point(ra2);
        tmp_score -= calculate_score_one_point(ra3);
        swap(moves_b[ra2][0], moves_b[ra3][0]);
        swap(moves_b[ra2][1], moves_b[ra3][1]);
        tmp_score += calculate_score_one_point(ra2);
        tmp_score += calculate_score_one_point(ra3);
        //tmp_score = calculate_score();
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算


    // 焼きなましで採用判定
    double diff_score = (current_score - tmp_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score < best_score) {
        store_best_score();
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        if (ra1 == 0) {
          swap(moves_a[ra2][0], moves_a[ra3][0]);
          swap(moves_a[ra2][1], moves_a[ra3][1]);
        }
        else {
          swap(moves_b[ra2][0], moves_b[ra3][0]);
          swap(moves_b[ra2][1], moves_b[ra3][1]);
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();

  for (int i = 0; i < A; ++i) {
    for (int j = 0; j < A; ++j) {
      if (moves_a[i][0].x == a[j].x) {
        nums_a[i] = j;
      }
    }
  }
  for (int i = 0; i < B; ++i) {
    for (int j = 0; j < B; ++j) {
      if (moves_b[i][0].x == b[j].x) {
        nums_b[i] = j;
      }
    }
  }

  store_best_score();
}

void run_simulated_annealing_merge(AnnealingParams annealingParams)
{
  A = ua.size();
  B = ub.size();

  current_score = calculate_score();
  store_best_score();
  //cout << current_score << endl;

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > 1.9) break;
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    double tmp_score = current_score;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % 2;
      if (ra1 == 0) {
        ra2 = rand_xorshift() % (A - 2) + 1;
        ra3 = rand_xorshift() % (A - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        for (int j = 0; j < 3; ++j) {
          swap(moves_a[ra2 * 3 + j][0], moves_a[ra3 * 3 + j][0]);
          swap(moves_a[ra2 * 3 + j][1], moves_a[ra3 * 3 + j][1]);
        }
        tmp_score = calculate_score();
      }
      else {
        if (B == 0) {
          continue;
        }
        ra2 = rand_xorshift() % (B - 2) + 1;
        ra3 = rand_xorshift() % (B - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        for (int j = 0; j < 3; ++j) {
          swap(moves_b[ra2 * 3 + j][0], moves_b[ra3 * 3 + j][0]);
          swap(moves_b[ra2 * 3 + j][1], moves_b[ra3 * 3 + j][1]);
        }
        tmp_score = calculate_score();
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算


    // 焼きなましで採用判定
    double diff_score = (current_score - tmp_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score < best_score) {
        store_best_score();
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        if (ra1 == 0) {
          for (int j = 0; j < 3; ++j) {
            swap(moves_a[ra2 * 3 + j][0], moves_a[ra3 * 3 + j][0]);
            swap(moves_a[ra2 * 3 + j][1], moves_a[ra3 * 3 + j][1]);
          }
        }
        else {
          for (int j = 0; j < 3; ++j) {
            swap(moves_b[ra2 * 3 + j][0], moves_b[ra3 * 3 + j][0]);
            swap(moves_b[ra2 * 3 + j][1], moves_b[ra3 * 3 + j][1]);
          }
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();

  for (int i = 0; i < A; ++i) {
    for (int j = 0; j < A; ++j) {
      if (moves_a[i][0].x == a[j].x) {
        nums_a[i] = j;
      }
    }
  }
  for (int i = 0; i < B; ++i) {
    for (int j = 0; j < B; ++j) {
      if (moves_b[i][0].x == b[j].x) {
        nums_b[i] = j;
      }
    }
  }

  store_best_score();
}

void run_simulated_annealing_box(AnnealingParams annealingParams)
{
  current_score = calculate_score();
  store_best_score();
  //cout << current_score << endl;

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  int isFirst = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > 1.9) break;
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    double tmp_score = current_score;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % 2;
      if (ra1 == 0) {
        ra2 = rand_xorshift() % (A - 2) + 1;
        ra3 = rand_xorshift() % (A - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        store_keep_score();
        swap(nums_a[ra2], nums_a[ra3]);
        tmp_score = simulate_score();
      }
      else {
        if (B == 0) {
          continue;
        }
        ra2 = rand_xorshift() % (B - 2) + 1;
        ra3 = rand_xorshift() % (B - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }


        store_keep_score();
        swap(nums_b[ra2], nums_b[ra3]);
        tmp_score = simulate_score();
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算


    // 焼きなましで採用判定
    double diff_score = (current_score - tmp_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01() || isFirst) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score < best_score || isFirst) {
        cerr << current_score << endl;
        store_best_score();
        isFirst = 0;
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        if (ra1 == 0) {
          restore_keep_score();
        }
        else {
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << "box = " << loop << endl;
  }

  restore_best_score();
  //for (int i = 0; i < (A); ++i) {
  //  cout << nums_a[i] << " ";
  //}
  //cout << endl;
}


void run_simulated_annealing_B(AnnealingParams annealingParams)
{
  for (int i = 0; i < B; ++i) {
    moves_b[i][0].x = 0;
    moves_b[i][0].y = 0;
    moves_b[i][1].x = 0;
    moves_b[i][1].y = 0;
  }
  moves_b_count = B;

  current_score = calculate_score();
  store_best_score();

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) break;
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[0];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    double tmp_score = current_score;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % 1;
      if (ra1 == 0) {
        ra2 = rand_xorshift() % (A - 2) + 1;
        ra3 = rand_xorshift() % (A - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        tmp_score -= calculate_score_one_point(ra2);
        tmp_score -= calculate_score_one_point(ra3);
        swap(moves_a[ra2][0], moves_a[ra3][0]);
        swap(moves_a[ra2][1], moves_a[ra3][1]);
        tmp_score += calculate_score_one_point(ra2);
        tmp_score += calculate_score_one_point(ra3);
      }
      else {
        if (B == 0) {
          continue;
        }
        ra2 = rand_xorshift() % (B - 2) + 1;
        ra3 = rand_xorshift() % (B - 2) + 1;
        if (ra2 > ra3) {
          swap(ra2, ra3);
        }
        if (ra3 <= ra2 + 1) {
          continue;
        }

        tmp_score -= calculate_score_one_point(ra2);
        tmp_score -= calculate_score_one_point(ra3);
        swap(moves_b[ra2][0], moves_b[ra3][0]);
        swap(moves_b[ra2][1], moves_b[ra3][1]);
        tmp_score += calculate_score_one_point(ra2);
        tmp_score += calculate_score_one_point(ra3);

      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算


    // 焼きなましで採用判定
    double diff_score = (current_score - tmp_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score < best_score) {
        store_best_score();
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        if (ra1 == 0) {
          swap(moves_a[ra2][0], moves_a[ra3][0]);
          swap(moves_a[ra2][1], moves_a[ra3][1]);
        }
        else {
          swap(moves_b[ra2][0], moves_b[ra3][0]);
          swap(moves_b[ra2][1], moves_b[ra3][1]);
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  //init_box();

  if (true) {
    build_initial_solution(1.9);
    run_simulated_annealing(annealingParams, 1.5);

    int loop = 0;
    while (moves_a_count < MAX_T && moves_b_count < MAX_T) {
      loop++;
      if (loop % 100 == 0) {
        if (get_elapsed_time() > 1.90) break;
      }

      int ra = rand_xorshift() % 2;
      store_keep_score();
      if (ra == 0) {
        int ra2 = rand_xorshift() % (moves_a_count - 1);
        for (int i = moves_a_count - 1; i >= 0; --i) {
          moves_a[i + 1][0] = moves_a[i][0];
          moves_a[i + 1][1] = moves_a[i][1];
          if (i == ra2) {
            break;
          }
        }
      }
      else {
        int ra2 = rand_xorshift() % (moves_b_count - 1);
        for (int i = moves_b_count - 1; i >= 0; --i) {
          moves_b[i + 1][0] = moves_b[i][0];
          moves_b[i + 1][1] = moves_b[i][1];
          if (i == ra2) {
            break;
          }
        }
      }

      double score = calculate_score();
      if (score < current_score) {
        current_score = score;
        //cout << current_score << endl;
      }
      else {
        restore_keep_score();
      }
    }

  }
  else {
    merge_gomi();
    build_initial_solution2();
    run_simulated_annealing_merge(annealingParams);
  }

  // 焼きなまし実行
  //if (B == 0 && C == 0) {
  //  run_simulated_annealing_B(annealingParams);
  //  moves_b_count = 0;
  //  for (int i = 0; i < A; ++i) {
  //    moves_b[moves_b_count][0].x = 0;
  //    moves_b[moves_b_count][0].y = 0;
  //    moves_b[moves_b_count][1].x = 0;
  //    moves_b[moves_b_count][1].y = 0;
  //    moves_b_count++;
  //  }
  //  moves_b[moves_b_count][0].x = 1000000;
  //  moves_b[moves_b_count][0].y = 0;
  //  moves_b[moves_b_count][1].x = 0;
  //  moves_b[moves_b_count][1].y = 1000000;
  //  moves_b_count++;
  //  moves_b[moves_b_count][0].x = 1000000;
  //  moves_b[moves_b_count][0].y = 1000000;
  //  moves_b[moves_b_count][1].x = 1000000;
  //  moves_b[moves_b_count][1].y = 1000000;
  //  moves_b_count++;
  //}
  //else {
  //  run_simulated_annealing(annealingParams);
  //  //run_simulated_annealing_box(annealingParams);
  //}

  // 解答を出力
  output_data(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main()
{
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 100000048.0;
  annealingParams.start_temperature[1] = 2048.0;
  annealingParams.start_temperature[2] = 2048.0;
  annealingParams.start_temperature[3] = 2048.0;
  annealingParams.start_temperature[4] = 2048.0;
  annealingParams.start_temperature[5] = 2048.0;
  annealingParams.start_temperature[6] = 2048.0;
  annealingParams.start_temperature[7] = 2048.0;
  annealingParams.start_temperature[8] = 2048.0;
  annealingParams.start_temperature[9] = 2048.0;
  annealingParams.end_temperature = 0.0;
  annealingParams.score_scale = 1234.0;
  annealingParams.operation_thresholds[0] = 100;
  annealingParams.operation_thresholds[1] = 200;
  annealingParams.operation_thresholds[2] = 300;
  annealingParams.operation_thresholds[3] = 400;
  annealingParams.operation_thresholds[4] = 500;
  annealingParams.operation_thresholds[5] = 600;
  annealingParams.operation_thresholds[6] = 700;
  annealingParams.operation_thresholds[7] = 800;
  annealingParams.operation_thresholds[8] = 900;
  annealingParams.operation_thresholds[9] = 1000;

  if (exec_mode == 0) {
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = solve_case(i, annealingParams);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(10) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << get_elapsed_time() << ", "
          << endl;
      }
    }
  }
  else if (exec_mode == 3) {
    int loop_count = 0;
    AnnealingParams best_annealingParams;
    ll best_sum_score = 0;

    while (true) {
      AnnealingParams new_annealingParams;
      new_annealingParams.start_temperature[0] = pow(2.0, rand_01() * 20);
      new_annealingParams.end_temperature = 0.0;
      new_annealingParams.score_scale = pow(2.0, rand_01() * 20);
      new_annealingParams.operation_thresholds[0] = rand() % 101;

      ll sum_score = 0;
      for (int i = 0; i < 15; ++i) {
        ll score = solve_case(i, new_annealingParams);
        sum_score += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cerr << "loop_count = " << loop_count
        << ", sum_score = " << sum_score
        << ", start_temperature = " << new_annealingParams.start_temperature[0]
        << ", end_temperature = " << new_annealingParams.end_temperature
        << ", score_scale = " << new_annealingParams.score_scale
        << ", operation_thresholds = " << new_annealingParams.operation_thresholds[0]
        << endl;

      if (sum_score > best_sum_score) {
        best_sum_score = sum_score;
        best_annealingParams = new_annealingParams;
      }

      loop_count++;
    }
  }

  return 0;
}
