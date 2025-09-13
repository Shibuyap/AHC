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

#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// タイマー
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

Timer timer;

// 乱数
namespace
{
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

  static uint64_t rand_xorshift64star()
  {
    static uint64_t x = 88172645463325252ULL;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * 2685821657736338717ULL;
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
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

const double TIME_LIMIT = 1.8;
int exec_mode;

const ll n = 500;
const ll m = 50;
const ll L = 998000000000000;
const ll U = 1002000000000000;
const ll DIFF = U - L;
const ll UNIT = 1000000000000;

ll a[n] = {
100000,
100000,
100000,
100000,
105855,
119903,
124287,
124411,
137492,
149959,
197789,
257971,
285076,
305435,
313930,
390921,
464925,
489809,
514063,
521340,
524841,
526629,
531781,
534073,
543632,
561574,
597385,
699379,
710630,
766063,
780154,
783157,
788456,
855506,
889247,
900170,
914851,
918450,
937027,
951458,
1166757,
1320817,
1339880,
1356478,
1456512,
1517641,
1564141,
1579488,
1618269,
1726656,
1772931,
1786328,
1822696,
2101573,
2207248,
2285293,
2398233,
2491637,
2512831,
2549226,
2589372,
2610774,
2915572,
2982112,
3128600,
3168420,
3168575,
3245550,
3267399,
3391751,
3484444,
3566922,
3779328,
3906501,
3972507,
4197597,
4210809,
4277123,
4814583,
5072247,
5085387,
5110837,
5121776,
5230687,
5233542,
5853322,
6035757,
6313160,
6328578,
6512756,
6539735,
6745251,
6822815,
6873882,
6905242,
7074878,
7076919,
7321605,
8084771,
8533625,
9244674,
10581539,
10796839,
10825622,
12300616,
12455032,
12747777,
13126139,
13218926,
13772112,
14417688,
14497066,
14623491,
14786050,
15249946,
16009858,
16806188,
17938477,
18073278,
18868424,
21462754,
22006928,
22375435,
23283377,
23965320,
24047151,
26056106,
26933696,
27984027,
28845988,
29598216,
30082942,
30151715,
31001412,
31418394,
34212884,
34288021,
36708621,
37752519,
37901216,
39090785,
39307749,
43459598,
44464344,
45118874,
48774934,
49736568,
50946407,
55611570,
56470212,
57866639,
57997295,
61336892,
62244598,
65468165,
69064504,
69071423,
72998412,
75019115,
76575542,
82918744,
85788423,
91856056,
95906781,
96942490,
112593826,
119812122,
129882225,
133636568,
136303067,
169441692,
171590532,
174199450,
178433178,
181552988,
181932487,
184553348,
191093899,
191928521,
195530356,
204178053,
207635792,
216171516,
217532322,
222241771,
222558466,
230333611,
231931978,
258839052,
270333747,
297002046,
304821545,
350830177,
375885406,
391496468,
403721179,
410055140,
423213380,
430103546,
473367870,
475099615,
488434220,
512349154,
519049506,
527852633,
552418287,
585235975,
717920238,
731373774,
750569608,
804277050,
820354823,
871854730,
902364867,
911875064,
916614350,
925008962,
963378073,
970419094,
1020128699,
1081893708,
1106509195,
1177166079,
1241456185,
1247328759,
1281766974,
1343777049,
1515455363,
1538649300,
1577245711,
1589471603,
1757079738,
1780511440,
1835200349,
1859332987,
1962387643,
1980721803,
2157888494,
2269413808,
2296122708,
2432604399,
2466808980,
2510864838,
2760960410,
2772803900,
2877539527,
3027119546,
3039494797,
3086824086,
3113358140,
3214301429,
3295234985,
3296592311,
3402548712,
3508255712,
3560836659,
3595084860,
3638398148,
3709020906,
3744070552,
3753226000,
4370887109,
4538070610,
4729204497,
5067577907,
5120357717,
5553355890,
5836940792,
6111073300,
6514215734,
6891500193,
7000368974,
7147419044,
7398553684,
7518262882,
8585075875,
9122302127,
9496593995,
9803715130,
10170948434,
10330546888,
10605305313,
10803872786,
11190424758,
11850864203,
11910098524,
11946537852,
12167079039,
12540752518,
12611139463,
12616227253,
13069470783,
13871741404,
14184776143,
14604798394,
15194620980,
15471268042,
17415398511,
17935906175,
17980214214,
18321629211,
19999773449,
20028757088,
20348473260,
20970743387,
20994853874,
21359241038,
21487983005,
21776691620,
22521770362,
23490612331,
24131594957,
24493123243,
24838478442,
26766760193,
26782443085,
28396718555,
28899444290,
28958356939,
30716581915,
31226672539,
33424439804,
35894612591,
36214583652,
36458966977,
38464277415,
38538572587,
38989837457,
39143914811,
44521390327,
46825594074,
46959446254,
47693927940,
47954094521,
50283978234,
50474694035,
50764718412,
54050994372,
54546476924,
56801736785,
58466480830,
59793674595,
60595606715,
61006091974,
61940637944,
68798246728,
70970096287,
74029757907,
79047287036,
79919942381,
82250769747,
84578985875,
87805780721,
92172975849,
92664333764,
93593800276,
99232430171,
99526041228,
100266742074,
100371670377,
102796338063,
104630351963,
105566380155,
108693814083,
111168185063,
119950020824,
124685918237,
128943914146,
131026223334,
136962652942,
138516352502,
141072968696,
142345497365,
150476854689,
150976440923,
152845814301,
173665498540,
173823879002,
182814520429,
188242419431,
199954772557,
202480257010,
228083856246,
230519409805,
230938180452,
235420955952,
236033862174,
244874654788,
248338561827,
250781650204,
256624973827,
272572244184,
275369632084,
277620306456,
283755154902,
287641429084,
288190359261,
292401180194,
292964885899,
317418673767,
330965884120,
359157511305,
360999982088,
370568259636,
379929298459,
383113131933,
383581935141,
396507510907,
420352471804,
423559500947,
445394210321,
458312491426,
460222839990,
465127980273,
468518524113,
468861848647,
487621613912,
508984993677,
513933003344,
546169757325,
547275919318,
564869692262,
578720120718,
595678638027,
598360720506,
606914346758,
611692561769,
631578741326,
639876618806,
659135142952,
662346100385,
719718548214,
738398296848,
739669048807,
752021981839,
759914447805,
778585513109,
820341180869,
851599442663,
880156747718,
890762058768,
896465998506,
907424826486,
921088036701,
922747251211,
923451551490,
929246564464,
971516328554,
981380427226,
981647169442,
998000034044677,
998002236767672,
998011403790918,
998031226254839,
998051800164007,
998096000278709,
998121206024899,
998146201152865,
998180773923707,
998200718745267,
998293693020828,
998310143647080,
998364786965019,
998371737100426,
998437352608622,
998478256108139,
998579222942559,
998661634101574,
998673458028637,
998761267751780,
998805541136735,
998830327469815,
998890773682802,
998891158849364,
998991127842631,
999170378608032,
999266758902520,
999357229535022,
999412199695845,
999469501786374,
999499611523256,
999667139041269,
999675767499375,
999733603389245,
999852297167028,
999913326758096,
999998169593035,
1000118743946853,
1000124272540824,
1000134376739151,
1000401167679874,
1000597689498456,
1000630005980648,
1000653867293697,
1000676919596635,
1000743017628855,
1001048666362414,
1001187843989340,
1001286916709463,
1001420809639081,
};

ll b[m + 10];
ll x[n];
ll sums[m + 10];
ll e_sum;

int current_score;

int best_score;

void store_best_score()
{
  best_score = current_score;
}

void restore_best_score()
{
  current_score = best_score;
}

void initialize_state()
{
  current_score = 0;
  rep(i, m)
  {
    sums[i] = 0;
  }
  e_sum = 0;
  rep(i, n)
  {
    x[i] = m;
  }
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    ll _n, _m, _l, _u;
    cin >> _n >> _m >> _l >> _u;
  }
  else {
    // ファイル入力
    ll _n, _m, _l, _u;
    ifs >> _n >> _m >> _l >> _u;
    rep(i, m)
    {
      ifs >> b[i];
    }
    ifs.close();
  }
}

void input_data2(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    rep(i, m)
    {
      cin >> b[i];
    }
  }
  else {
    // ファイル入力
    // 何もしない
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

void output_data(ofstream& ofs)
{
  if (exec_mode == 0) {
    // 標準出力
    rep(i, n)
    {
      cout << a[i] << ' ';
    }
    cout << endl;
    fflush(stdout);
  }
  else {
    // ファイル出力
    rep(i, n)
    {
      ofs << a[i] << ' ';
    }
    ofs << endl;
  }
}

void output_data2(ofstream& ofs)
{
  if (exec_mode == 0) {
    // 標準出力
    rep(i, n)
    {
      if (x[i] == m) {
        cout << 0 << ' ';
      }
      else {
        cout << x[i] + 1 << ' ';
      }
    }
    cout << endl;
    fflush(stdout);
  }
  else {
    // ファイル出力
    rep(i, n)
    {
      if (x[i] == m) {
        ofs << 0 << ' ';
      }
      else {
        ofs << x[i] + 1 << ' ';
      }
    }
    ofs << endl;
  }
}

ll calculate_score()
{
  ll res = round((20.0 - log10(1 + e_sum)) * 5 * 1e7);
  return res;
}

void calc_sum_e()
{
  e_sum = 0;
  rep(i, m)
  {
    e_sum += abs(sums[i] - b[i]);
  }
}

void build_initial_solution()
{
  const ll all_sum = 500 * 1000000000000000;

  ll aa[m];
  rep(j, m)
  {
    aa[j] = U;
  }
  ll aaa[m];
  rep(i, 1500)
  {
    rep(j, m)
    {
      aaa[j] = L + rand_xorshift64star() % DIFF;
    }
    sort(aaa, aaa + m);
    rep(j, m)
    {
      aa[j] = min(aa[j], aaa[j]);
    }
  }

  rep(i, n)
  {
    if (i < m) {
      a[i] = aa[i];
    }
    else {
      a[i] = rand_xorshift64star() % (UNIT);
      int ra = rand_xorshift() % 24;
      ra = max(0, ra - 6);
      rep(j, ra)
      {
        a[i] /= 2;
      }
      a[i] += 1;
      a[i] = max(a[i], 100000LL);
    }
  }
  sort(a, a + n);
}

void build_x()
{
  drep(i, n)
  {
    ll ma = -1;
    ll ma_j = -1;
    if (n - i <= m) {
      ma = DIFF + 1;
      ma_j = i - (n - m);
    }
    else {
      drep(j, m)
      {
        if (sums[j] + a[i] <= b[j]) {
          if (b[j] - sums[j] > ma) {
            ma = b[j] - sums[j];
            ma_j = j;
          }
        }
      }
    }
    if (ma_j != -1) {
      x[i] = ma_j;
      sums[ma_j] += a[i];
    }
  }

  drep(i, n)
  {
    if (x[i] != m) continue;
    ll ma = -1;
    ll ma_j = -1;
    drep(j, m)
    {
      ll diff1 = abs(sums[j] - b[j]);
      ll diff2 = abs(sums[j] + a[i] - b[j]);
      if (diff2 < diff1) {
        if (diff1 - diff2 > ma) {
          ma = diff1 - diff2;
          ma_j = j;
        }
      }
    }
    if (ma_j != -1) {
      x[i] = ma_j;
      sums[ma_j] += a[i];
    }
  }

  calc_sum_e();
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams)
{
  store_best_score();

  double now_time = timer.get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = timer.get_elapsed_time();
      if (now_time > TIME_LIMIT) { break; }
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[2];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5, keep6;

    ll now_diff;
    ll next_diff;
    vector<ll> use1;
    vector<ll> use2;

    int ramode = 0;
    if (now_time < 0.5 || ra_exec_mode < annealingParams.operation_thresholds[0]) {
      ramode = 0;
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      ramode = 1;
    }
    else {
      ramode = 2;
    }

    if (ramode == 0) {
      // 近傍操作1
      ra1 = rand_xorshift() % m;
      vector<int> candidates;
      vector<int> use;
      rep(i, n)
      {
        if (x[i] == m) {
          candidates.push_back(i);
        }
        if (x[i] == ra1) {
          use.push_back(i);
        }
      }
      if (use.size() < 3) {
        continue;
      }
      ra2 = rand_xorshift() % (use.size() - 1);
      ra3 = rand_xorshift() % (use.size() - 1);
      while (ra2 == ra3) {
        ra3 = rand_xorshift() % (use.size() - 1);
      }
      keep1 = use[ra2];
      keep2 = use[ra3];

      now_diff = abs(sums[ra1] - b[ra1]);
      sums[ra1] -= a[keep1];
      sums[ra1] -= a[keep2];
      x[keep1] = m;
      x[keep2] = m;
      next_diff = abs(sums[ra1] - b[ra1]);
      keep3 = -1;
      keep4 = -1;
      rep(i, candidates.size())
      {
        ll new_diff = abs(sums[ra1] + a[candidates[i]] - b[ra1]);
        if (new_diff < next_diff) {
          next_diff = new_diff;
          keep3 = candidates[i];
        }
      }
      rep(i, candidates.size())
      {
        srep(j, i + 1, candidates.size())
        {
          ll new_diff = abs(sums[ra1] + a[candidates[i]] + a[candidates[j]] - b[ra1]);
          if (new_diff < next_diff) {
            next_diff = new_diff;
            keep3 = candidates[i];
            keep4 = candidates[j];
          }
        }
      }
    }
    else if (ramode == 1) {
      // 近傍操作2
      ra1 = rand_xorshift() % m;
      vector<int> candidates;
      vector<int> use;
      rep(i, n)
      {
        if (x[i] == m) {
          candidates.push_back(i);
        }
        if (x[i] == ra1) {
          use.push_back(i);
        }
      }
      int remov = rand_xorshift() % 3 + 1;
      if (use.size() < remov + 1) {
        continue;
      }
      ra2 = rand_xorshift() % (use.size() - 1);
      ra3 = -1;
      if (remov >= 2) {
        ra3 = rand_xorshift() % (use.size() - 1);
        while (ra2 == ra3) {
          ra3 = rand_xorshift() % (use.size() - 1);
        }
      }
      ra4 = -1;
      if (remov == 3) {
        ra4 = rand_xorshift() % (use.size() - 1);
        while (ra2 == ra4 || ra3 == ra4) {
          ra4 = rand_xorshift() % (use.size() - 1);
        }
      }
      keep1 = use[ra2];
      keep2 = -1;
      if (ra3 != -1) {
        keep2 = use[ra3];
      }
      keep5 = -1;
      if (ra4 != -1) {
        keep5 = use[ra4];
      }

      now_diff = abs(sums[ra1] - b[ra1]);
      sums[ra1] -= a[keep1];
      if (keep2 != -1) {
        sums[ra1] -= a[keep2];
      }
      if (keep5 != -1) {
        sums[ra1] -= a[keep5];
      }
      x[keep1] = m;
      if (keep2 != -1) {
        x[keep2] = m;
      }
      if (keep5 != -1) {
        x[keep5] = m;
      }
      next_diff = abs(sums[ra1] - b[ra1]);
      keep3 = -1;
      keep4 = -1;
      keep6 = -1;
      rep(i, candidates.size())
      {
        ll new_diff = abs(sums[ra1] + a[candidates[i]] - b[ra1]);
        if (new_diff < next_diff) {
          next_diff = new_diff;
          keep3 = candidates[i];
        }
      }
      rep(i, candidates.size())
      {
        srep(j, i + 1, candidates.size())
        {
          ll new_diff = abs(sums[ra1] + a[candidates[i]] + a[candidates[j]] - b[ra1]);
          if (new_diff < next_diff) {
            next_diff = new_diff;
            keep3 = candidates[i];
            keep4 = candidates[j];
          }
        }
      }
      rep(i, candidates.size())
      {
        srep(j, i + 1, candidates.size())
        {
          srep(k, j + 1, candidates.size())
          {
            ll new_diff = abs(sums[ra1] + a[candidates[i]] + a[candidates[j]] + a[candidates[k]] - b[ra1]);
            if (new_diff < next_diff) {
              next_diff = new_diff;
              keep3 = candidates[i];
              keep4 = candidates[j];
              keep6 = candidates[k];
            }
          }
        }
      }
    }
    else if (ramode == 2) {
      ra1 = rand_xorshift() % m;
      ra2 = rand_xorshift() % m;
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % m;
      }

      now_diff = abs(sums[ra1] - b[ra1]) + abs(sums[ra2] - b[ra2]);
      next_diff = now_diff;

      use1.clear();
      use2.clear();


      rep(i, n)
      {
        if (x[i] == ra1) {
          use1.push_back(i);
        }
        if (x[i] == ra2) {
          use2.push_back(i);
        }
      }
      if (use1.size() == 0 || use2.size() == 0) {
        continue;
      }

      vector<ll> candidates1;
      rep(i, (1 << use1.size()))
      {
        bitset<20> bs(i);
        ll sum1 = 0;
        rep(j, use1.size())
        {
          if (bs.test(j)) {
            sum1 += a[use1[j]];
          }
        }
        candidates1.push_back(sum1);
      }
      vector<ll> candidates2;
      rep(i, (1 << use2.size()))
      {
        bitset<20> bs(i);
        ll sum2 = 0;
        rep(j, use2.size())
        {
          if (bs.test(j)) {
            sum2 += a[use2[j]];
          }
        }
        candidates2.push_back(sum2);
      }

      keep1 = -1;
      keep2 = -1;
      rep(i, candidates1.size())
      {
        rep(j, candidates2.size())
        {
          ll new_diff = abs((sums[ra1] - candidates1[i] + candidates2[j]) - b[ra1]) + abs((sums[ra2] - candidates2[j] + candidates1[i]) - b[ra2]);
          if (new_diff < next_diff) {
            next_diff = new_diff;
            keep1 = i;
            keep2 = j;
          }
        }
      }
    }

    // スコア計算
    ll tmp_score = now_diff - next_diff;

    // 焼きなましで採用判定
    double diff_score = tmp_score * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      //if (tmp_score >= 0) {
        // 採用
      if (ramode == 0) {
        // 近傍操作1 の適用
        if (keep3 != -1) {
          sums[ra1] += a[keep3];
          x[keep3] = ra1;
        }
        if (keep4 != -1) {
          sums[ra1] += a[keep4];
          x[keep4] = ra1;
        }
      }
      else if (ramode == 1) {
        // 近傍操作2 の適用
        if (keep3 != -1) {
          sums[ra1] += a[keep3];
          x[keep3] = ra1;
        }
        if (keep4 != -1) {
          sums[ra1] += a[keep4];
          x[keep4] = ra1;
        }
        if (keep6 != -1) {
          sums[ra1] += a[keep6];
          x[keep6] = ra1;
        }
      }
      else if (ramode == 2) {
        // 近傍操作3 の適用
        if (keep1 != -1) {
          //cerr << "keep1 = " << keep1 << ", keep2 = " << keep2 << endl;
          ll sum1 = 0;
          bitset<20> bs(keep1);
          rep(j, use1.size())
          {
            if (bs.test(j)) {
              sum1 += a[use1[j]];
              sums[ra1] -= a[use1[j]];
              x[use1[j]] = ra2;
            }
          }
          sums[ra2] += sum1;
        }
        if (keep2 != -1) {
          ll sum2 = 0;
          bitset<20> bs(keep2);
          rep(j, use2.size())
          {
            if (bs.test(j)) {
              sum2 += a[use2[j]];
              sums[ra2] -= a[use2[j]];
              x[use2[j]] = ra1;
            }
          }
          sums[ra1] += sum2;
        }
      }
    }
    else {
      // 元に戻す
      if (ramode == 0) {

        // 近傍操作1 の巻き戻し
        sums[ra1] += a[keep1];
        sums[ra1] += a[keep2];
        x[keep1] = ra1;
        x[keep2] = ra1;
      }
      else if (ramode == 1) {
        // 近傍操作2 の巻き戻し
        sums[ra1] += a[keep1];
        if (keep2 != -1) {
          sums[ra1] += a[keep2];
        }
        if (keep5 != -1) {
          sums[ra1] += a[keep5];
        }
        x[keep1] = ra1;
        if (keep2 != -1) {
          x[keep2] = ra1;
        }
        if (keep5 != -1) {
          x[keep5] = ra1;
        }
      }
      else if (ramode == 2) {
        // 近傍操作3 の巻き戻し
      }
    }
    calc_sum_e();
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  timer.start();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  //build_initial_solution();

  // 解答を出力
  output_data(ofs);

  input_data2(case_num);

  build_x();

  run_simulated_annealing(annealingParams);

  output_data2(ofs);

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

  ll best_a[n];
  ll best_score = 0;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 2e8;
  annealingParams.start_temperature[1] = 20.0;
  annealingParams.start_temperature[2] = 2048.0;
  annealingParams.start_temperature[3] = 2048.0;
  annealingParams.start_temperature[4] = 2048.0;
  annealingParams.start_temperature[5] = 2048.0;
  annealingParams.start_temperature[6] = 2048.0;
  annealingParams.start_temperature[7] = 2048.0;
  annealingParams.start_temperature[8] = 2048.0;
  annealingParams.start_temperature[9] = 2048.0;
  annealingParams.end_temperature = 0.0;
  annealingParams.score_scale = 12345.0;
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

  build_initial_solution();
  if (exec_mode == 0) {
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    rep(loop, 1)
    {
      ll sum_score = 0;
      for (int i = 62; i < 100; ++i) {
        ll score = solve_case(i, annealingParams);
        sum_score += score;
        if (exec_mode == 1) {
          cerr << score << endl;
        }
        else {
          cerr << "case = " << setw(2) << i << ", "
            << "score = " << setw(4) << score << ", "
            << "sum = " << setw(5) << sum_score << ", "
            << "time = " << setw(5) << timer.get_elapsed_time() << ", "
            << endl;
        }
      }
      //cout << "loop = " << loop << ", sum_score = " << sum_score << endl;
      //if (sum_score > best_score) {
      //  best_score = sum_score;
      //  rep(i, n)
      //  {
      //    best_a[i] = a[i];
      //  }

      //  cout << "loop = " << loop << ", best_score = " << best_score << endl;
      //  cout << "ll a[n] = {" << endl;
      //  rep(i, n)
      //  {
      //    cout << best_a[i] << ", " << endl;
      //  }
      //  cout << "};" << endl;
      //}
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
      for (int i = 5; i < 6; ++i) {
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
