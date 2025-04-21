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
#include <utility>
#include <vector>
// clang-format off
#ifdef ONLINE_JUDGE
const int LOCAL = 0;
const double CLOCK_RATIO = 1000;
#else
const int LOCAL = 1;
const double CLOCK_RATIO = 1000;
#endif

using namespace std;
using ll = long long;
using ull = unsigned long long;
using pii = pair<int, int>;
const int inf = 1e9;
struct point { int x; int y; };
bool operator==(const point& lhs, const point& rhs) { return (lhs.x == rhs.x && lhs.y == rhs.y); }
bool operator!=(const point& lhs, const point& rhs) { return !(lhs == rhs); }
bool operator<(const point& lhs, const point& rhs) {
  if (lhs.x != rhs.x) { return lhs.x < rhs.x; }
  return lhs.y < rhs.y;
}
point operator-(const point& self) {
  return { -self.x, -self.y };
}
point operator+(const point& lhs, const point& rhs) {
  return { lhs.x + rhs.x, lhs.y + rhs.y };
}
point operator-(const point& lhs, const point& rhs) {
  return lhs + (-rhs);
}
std::ostream& operator<<(std::ostream& os, point& pt) {
  string s;
  s = "(" + to_string(int(pt.x)) + ", " + to_string(int(pt.y)) + ")";
  return os << s;
};

#define debug1(a) {if(LOCAL){cerr<<#a<<":"<<a<<endl; }}
#define debug2(a,b) {if(LOCAL){cerr<<#a<<":"<<a<<" "<<#b<<":"<<b<<endl; }}
#define debug3(a,b,c) {if(LOCAL){cerr<<#a<<":"<<a<<" "<<#b<<":"<<b<<" "<<#c<<":"<<c<<endl; }}
#define debug4(a,b,c,d) {if(LOCAL){cerr<<#a<<":"<<a<<" "<<#b<<":"<<b<<" "<<#c<<":"<<c<<" "<<#d<<":"<<d<<endl; }}
#define debug5(a,b,c,d,e) {if(LOCAL){cerr<<#a<<":"<<a<<" "<<#b<<":"<<b<<" "<<#c<<":"<<c<<" "<<#d<<":"<<d<<" "<<#e<<":"<<e<<endl; }}
// clang-format on

int B;
int C;
const int A = 100;
const int BLIM = 100;
const int CLIM = 100;
const int POSLIM = 1000000;
vector<point> INIT_APOS;
vector<point> INIT_BPOS;
vector<point> INIT_CPOS;

namespace Solver {
  struct hand {
    point l;
    point r;
  };
  struct operation_t {
    hand alc;
    hand bob;
  };
  bool operator<(const operation_t& lhs, const operation_t& rhs) {
    return false;
  }
  struct state_t {
    hand alc;
    hand bob;
    vector<point> apos;
    vector<point> bpos;
    // vector<point> cpos;
  };
  int calc_distance(point a, point b) {
    return sqrt(ll(a.x - b.x) * ll(a.x - b.x) + ll(a.y - b.y) * ll(a.y - b.y));
  }
  int orient(point a, point b, point p) {
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
  }
  bool contians(point p, point a, point b, point c) {
    if (orient(a, b, c) == 0) {
      if (orient(a, b, p) != 0 || orient(a, c, p) != 0) {
        return false;
      }
      vector<int> xs = { a.x, b.x, c.x };
      vector<int> ys = { a.y, b.y, c.y };
      // return min(xs) <= p[0] <= max(xs) and min(ys) <= p[1] <= max(ys)
      int min_xs = min(a.x, min(b.x, c.x));
      int max_xs = max(a.x, max(b.x, c.x));
      int min_ys = min(a.y, min(b.y, c.y));
      int max_ys = max(a.y, max(b.y, c.y));
      return min_xs <= p.x && p.x <= max_xs && min_ys <= p.y && p.y <= max_ys;
    }

    int c1 = orient(a, b, p);
    int c2 = orient(b, c, p);
    int c3 = orient(c, a, p);
    return (c1 >= 0 && c2 >= 0 && c3 >= 0) || (c1 <= 0 && c2 <= 0 && c3 <= 0);
  }
  /*

  def orient(a, b, p):
      return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])

  def contains(p, a, b, c):
      # degenerate (colinear) case
      if orient(a, b, c) == 0:
          if orient(a, b, p) != 0 or orient(a, c, p) != 0:
              return False
          xs = (a[0], b[0], c[0])
          ys = (a[1], b[1], c[1])
          return min(xs) <= p[0] <= max(xs) and min(ys) <= p[1] <= max(ys)

      # non]degenerate: check same side of all edges
      c1 = orient(a, b, p)
      c2 = orient(b, c, p)
      c3 = orient(c, a, p)
      return (c1 >= 0 and c2 >= 0 and c3 >= 0) or (c1 <= 0 and c2 <= 0 and c3 <= 0)
  */
  class Solver {
  public:
    point rot_point(int rot, point p) {
      if (rot == 0) {
        return p;  // ‚»‚Ì‚Ü‚Ü
      }
      if (rot == 1) {
        return { POSLIM - p.x, POSLIM - p.y };  // x”½“]
      }
      if (rot == 2) {
        return { p.y, p.x };  // “]’u
      }
      if (rot == 3) {
        return { POSLIM - p.y, POSLIM - p.x };  //  p.x};  // “]’u‚µ‚Äx”½“]
      }
      assert(0);
      return p;
    }
    void solve() {
      pair<double, pair<int, vector<operation_t>>> bestops;
      bestops.first = 1e18;

      int bymin = POSLIM;
      int bymax = 0;

      for (int rot = 0; rot < 4; rot++) {
        vector<operation_t> ops;
        auto st = state_t();
        double totalcost = 0;
        vector<pair<int, point>> allpos;
        {
          for (auto a : INIT_APOS) {
            allpos.push_back({ 0, rot_point(rot, a) });
          }
          for (auto b : INIT_BPOS) {
            allpos.push_back({ 1, rot_point(rot, b) });
          }

          int axmin = inf;
          for (auto a : INIT_APOS) {
            axmin = min(rot_point(rot, a).x, axmin);
          }
          int bxmin = inf;
          for (auto b : INIT_BPOS) {
            bxmin = min(rot_point(rot, b).x, bxmin);
          }
          for (auto b : INIT_BPOS) {
            bymin = min(rot_point(rot, b).y, bymin);
            bymax = max(rot_point(rot, b).y, bymax);
          }

          int mx = min(axmin, bxmin);
          st.alc = { {mx - 1, 0}, {mx - 1, POSLIM} };
          st.bob = { {mx - 1, bymin}, {mx - 1, bymax} };
        }
        ops.push_back({ st.alc, st.bob });

        sort(allpos.begin(), allpos.end(), [](pair<int, point> p, pair<int, point> q) {
          return p.second.x < q.second.x;
          });

        map<point, int> owner;
        while (allpos.size()) {
          {
            auto old_alc = st.alc;
            auto old_bob = st.bob;

            int nxta = st.alc.l.x;
            int nxtb = st.bob.l.x;
            array<int, 2> minpos = { inf, inf };
            array<int, 2> maxpos = { -1, -1 };
            for (auto p : allpos) {
              minpos[p.first] = min(minpos[p.first], p.second.x);
              maxpos[p.first] = max(maxpos[p.first], p.second.x);
            }
            bool e_a = minpos[0] >= 0;
            bool e_b = minpos[1] >= 0;
            if (e_a && e_b) {
              if (minpos[0] == minpos[1]) {
                int mx = minpos[0];
                vector<pair<int, point>> next_allpos;
                for (auto p : allpos) {
                  if (p.second.x == minpos[0]) {
                    if (p.first == 0) {
                      st.alc = { p.second, p.second };
                    }
                    else {
                      st.bob = { p.second, p.second };
                    }
                  }
                  else {
                    next_allpos.push_back(p);
                  }
                }
                allpos = next_allpos;
                totalcost += max(calc_distance(old_alc.l, st.alc.l) + calc_distance(old_alc.r, st.alc.r),  //
                  calc_distance(old_bob.l, st.bob.l) + calc_distance(old_bob.r, st.bob.r));
                ops.push_back({ st.alc, st.bob });
                old_alc = st.alc;
                old_bob = st.bob;

                st.alc.l.y = 0;
                st.alc.r.y = POSLIM;
                st.bob.l.y = bymin;
                st.bob.r.y = bymax;

                totalcost += max(calc_distance(old_alc.l, st.alc.l) + calc_distance(old_alc.r, st.alc.r),  //
                  calc_distance(old_bob.l, st.bob.l) + calc_distance(old_bob.r, st.bob.r));
                ops.push_back({ st.alc, st.bob });
                continue;
              }
              else {
                nxta = min(maxpos[0], minpos[1] - 1);
                for (auto p : allpos) {
                  if (p.second.x >= minpos[1] && p.first == 0) {
                    if (bymin <= p.second.y && p.second.y <= bymax) {
                      break;
                    }
      
                  }
                  nxtb = p.second.x;
                }
                // debug2(nxta, nxtb);
              }
            }
            else if (e_a) {
              nxta = maxpos[0];
            }
            else if (e_b) {
              nxtb = maxpos[1];
            }
            else {
              assert(0);
            }

            // hand old_alc = st.alc;
            // hand old_bob = st.bob;
            st.alc = { {nxta, 0}, {nxta, POSLIM} };
            st.bob = { {nxtb, bymin}, {nxtb, bymax} };
            hand new_alc = st.alc;
            hand new_bob = st.bob;

            ops.push_back({ st.alc, st.bob });

            vector<pair<int, point>> next_allpos;
            for (auto p : allpos) {
              if (old_alc.l.x <= p.second.x && p.second.x <= new_alc.l.x) {
                if (p.first == 1) {
                  debug1("failure!");
                }
                continue;
              }
              else if (old_bob.l.x <= p.second.x && p.second.x <= new_bob.l.x) {
                if (p.first == 0) {
                  debug1("failure!");
                }
                continue;
              }
              next_allpos.push_back(p);
            }
            allpos = next_allpos;
            totalcost += max(calc_distance(old_alc.l, st.alc.l) + calc_distance(old_alc.r, st.alc.r),  //
              calc_distance(old_bob.l, st.bob.l) + calc_distance(old_bob.r, st.bob.r));
          }
        }
        debug3(bestops.first, totalcost, rot);
        if (bestops.first > totalcost) {
          bestops = { totalcost, {rot, ops} };
        }
      }
      {
        auto ops = bestops.second.second;
        int rot = bestops.second.first;

        for (auto op : ops) {
          op.alc.l = rot_point(rot, op.alc.l);
          op.alc.r = rot_point(rot, op.alc.r);
          op.bob.l = rot_point(rot, op.bob.l);
          op.bob.r = rot_point(rot, op.bob.r);
          cout << op.alc.l.x << " " << op.alc.l.y << " " << op.alc.r.x << " " << op.alc.r.y << " " <<  //
            op.bob.l.x << " " << op.bob.l.y << " " << op.bob.r.x << " " << op.bob.r.y << endl;
        }
      }
    }
  };
}  // namespace Solver
int main() {
  int a;
  cin >> a >> B >> C;
  for (int i = 0; i < A; i++) {
    int x, y;
    cin >> x >> y;
    INIT_APOS.push_back({ x, y });
  }
  for (int i = 0; i < B; i++) {
    int x, y;
    cin >> x >> y;
    INIT_BPOS.push_back({ x, y });
  }
  for (int i = 0; i < C; i++) {
    int x, y;
    cin >> x >> y;
    INIT_CPOS.push_back({ x, y });
  }
  Solver::Solver().solve();

  return 0;
}
