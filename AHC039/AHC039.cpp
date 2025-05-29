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

static uint32_t Rand() {
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

static double Rand01() {
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}

static double RandRange(double l, double r) {
  return l + (r - l) * Rand01();
}

// [l, r]
static uint32_t RandRange(uint32_t l, uint32_t r) {
  return l + Rand() % (r - l + 1);
}


void FisherYates(int* data, int n) {
  for (int i = n - 1; i >= 0; i--) {
    int j = Rand() % (i + 1);
    int swa = data[i];
    data[i] = data[j];
    data[j] = swa;
  }
}

// �����_���f�o�C�X�ƃ����Z���k�E�c�C�X�^
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

void ResetTime() {
  startTimeClock = std::chrono::steady_clock::now();
}

double GetNowTime() {
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTimeClock;
  return elapsed.count();
}

// �񎟌����W��\���\����
struct Point {
public:
  int x;
  int y;

  Point() { x = 0; y = 0; }                  // �f�t�H���g�R���X�g���N�^
  Point(int _x, int _y) { x = _x; y = _y; }  // ���W���w�肷��R���X�g���N�^
};

const int MAX_N = 30;  // ���g�p�̒萔�i����̊g���p�H�j

const int n = 5000;    // ���̐��i�T�o�ƃC���V���ꂼ��̐��j

vector<Point> saba, iwashi;  // �T�o�ƃC���V�̍��W���i�[����x�N�^�[

vector<Point> ans;  // �o�͂���|���S���̒��_���W���i�[����x�N�^�[

int ansScore;       // ���݂̃X�R�A
int best_ansScore;  // �ŗǂ̃X�R�A�i���g�p�H�j

int f[510][510];
int best_f[510][510];

// ���݂̉𓚂��ŗǂ̉𓚂Ƃ��ĕۑ�����֐��i���g�p�j
void CopyToBest(int blockSize) {
  best_ansScore = ansScore;
  rep(i, blockSize + 2) {
    rep(j, blockSize + 2) {
      best_f[i][j] = f[i][j];
    }
  }
}

// �ŗǂ̉𓚂����݂̉𓚂Ƃ��Đݒ肷��֐��i���g�p�j
void CopyToAns(int blockSize) {
  ansScore = best_ansScore;
  rep(i, blockSize + 2) {
    rep(j, blockSize + 2) {
      f[i][j] = best_f[i][j];
    }
  }
}

// �����̃P�[�X����������ۂɁA������Ԃ�����������֐�
void SetUp() {
  ansScore = 0;   // �X�R�A�̏�����
  rep(i, 510)rep(j, 510)best_f[i][j] = 0;

  ans.clear();    // �𓚂̏�����
  saba.clear();   // �T�o�̍��W�̏�����
  iwashi.clear(); // �C���V�̍��W�̏�����
}

// ���͂��󂯎��֐�
void Input(int problemNum) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";  // ���̓t�@�C���̃p�X���쐬
  ifstream ifs(oss.str());  // �t�@�C���X�g���[�����J��

  saba.resize(n);    // �T�o�̃x�N�^�[�����T�C�Y
  iwashi.resize(n);  // �C���V�̃x�N�^�[�����T�C�Y

  if (!ifs.is_open()) {
    // �W�����͂���̓ǂݍ���
    int _n;
    cin >> _n;  // ���̐��i���g�p�j
    rep(i, n) cin >> saba[i].x >> saba[i].y;     // �T�o�̍��W��ǂݍ���
    rep(i, n) cin >> iwashi[i].x >> iwashi[i].y; // �C���V�̍��W��ǂݍ���
  }
  else {
    // �t�@�C������̓ǂݍ���
    int _n;
    ifs >> _n;  // ���̐��i���g�p�j
    rep(i, n) ifs >> saba[i].x >> saba[i].y;     // �T�o�̍��W��ǂݍ���
    rep(i, n) ifs >> iwashi[i].x >> iwashi[i].y; // �C���V�̍��W��ǂݍ���
  }
}

// �o�̓t�@�C���X�g���[�����J���֐�
void OpenOfs(int probNum, ofstream& ofs) {
  if (mode != 0) {
    // ���[�h��0�ȊO�̏ꍇ�A�o�̓t�@�C�����쐬
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// �X�R�A���v�Z����֐�
ll CalcScore() {
  ll res = ansScore + 1;  // ��蕶�̓��_�v�Z���Ɋ�Â��imax(0, a - b + 1)�j
  return res;             // �X�R�A��Ԃ�
}

// �𓚂��o�͂���֐�
void Output(ofstream& ofs) {
  if (mode == 0) {
    // �W���o�͂ɏo��
    cout << ans.size() << endl;  // �|���S���̒��_�����o��
    for (auto p : ans) cout << p.x << ' ' << p.y << endl;  // �e���_�̍��W���o��
  }
  else {
    // �t�@�C���ɏo��
    ofs << ans.size() << endl;  // �|���S���̒��_�����o��
    for (auto p : ans) ofs << p.x << ' ' << p.y << endl;  // �e���_�̍��W���o��
  }
}

// �|���S���̕ӂ̑��a������𖞂����Ă��邩�m�F����֐�
bool IsLengthOK(vector<Point> vp) {
  int len = 0;
  // �e�ӂ̒������v�Z
  rep(i, vp.size() - 1) {
    len += abs(vp[i + 1].x - vp[i].x);  // x���W�̍��̐�Βl�����Z
    len += abs(vp[i + 1].y - vp[i].y);  // y���W�̍��̐�Βl�����Z
  }

  // �Ō�̒��_�ƍŏ��̒��_�����ԕӂ̒��������Z
  len += abs(vp[0].x - vp.back().x);
  len += abs(vp[0].y - vp.back().y);

  // ���a��400,000�ȉ��ł����true��Ԃ�
  return len <= 400000;
}

// ���W���O���b�h�͈̔͊O�����m�F����֐�
bool IsNG(int x, int y, int blockSize) {
  if (x < 0 || blockSize <= x || y < 0 || blockSize <= y) return true;  // �͈͊O�Ȃ�true
  return false;  // �͈͓��Ȃ�false
}

int block[510][510];  // �e�O���b�h�Z���̃X�R�A���i�[����z��
void InitBlock(const int blockSize) {
  rep(i, blockSize) rep(j, blockSize) block[i][j] = 0;  // ������

  // �T�o�ƃC���V�̈ʒu����A�e�Z���̃X�R�A���v�Z
  rep(i, n) {
    {
      // �T�o�̍��W���Z���Ɋ��蓖�āA�X�R�A�����Z
      int xx = saba[i].x / (100000 / blockSize);
      xx = min(xx, blockSize - 1);  // �ő�l�𒴂��Ȃ��悤�ɒ���
      int yy = saba[i].y / (100000 / blockSize);
      yy = min(yy, blockSize - 1);

      block[xx][yy]++;  // �T�o������Z���̃X�R�A��+1
    }

    {
      // �C���V�̍��W���Z���Ɋ��蓖�āA�X�R�A�����Z
      int xx = iwashi[i].x / (100000 / blockSize);
      xx = min(xx, blockSize - 1);
      int yy = iwashi[i].y / (100000 / blockSize);
      yy = min(yy, blockSize - 1);

      block[xx][yy]--;  // �C���V������Z���̃X�R�A��-1
    }
  }
}

int haba[510][510];  // ���D��T���p�̔z��
void Method3_SA(const int xx1, const int xx2, const int yy1, const int yy2, const int blockSize, int& loop2, double timeLimit) {
  queue<P> que;  // ���D��T���̂��߂̃L���[

  double nowTime = GetNowTime();  // ���݂̌o�ߎ���
  const double START_TEMP = 200.0;  // �Ă��Ȃ܂��@�̊J�n���x
  const double END_TEMP = 0.1;      // �Ă��Ȃ܂��@�̏I�����x
  double temp = START_TEMP + (END_TEMP - START_TEMP) * nowTime / timeLimit;  // ���݂̉��x

  while (true) {
    if (loop2 % 1 == 0) {
      nowTime = GetNowTime();
      if (nowTime > timeLimit) break;  // ���Ԑ������߂����烋�[�v�𔲂���
    }
    loop2++;

    int rax = Rand() % blockSize;  // �����_���ȃZ����x�C���f�b�N�X
    int ray = Rand() % blockSize;  // �����_���ȃZ����y�C���f�b�N�X

    int ng = 1;  // �ύX���\���ǂ����̃t���O
    rep(i, 4) {
      int nx = rax + dx[i];
      int ny = ray + dy[i];
      if (IsNG(nx, ny, blockSize)) continue;  // �͈͊O�͖���
      if (f[nx + 1][ny + 1] != f[rax + 1][ray + 1]) ng = 0;  // �אڃZ�����قȂ��ԂȂ�ύX�\
    }

    if (ng) continue;  // �ύX�s�Ȃ玟�̃��[�v��

    int tmpScore = ansScore;
    if (f[rax + 1][ray + 1] == 0) {
      tmpScore += block[rax][ray];  // �Z����ǉ������ꍇ�̃X�R�A
    }
    else {
      tmpScore += -block[rax][ray];  // �Z�����폜�����ꍇ�̃X�R�A
    }

    const double progressRatio = nowTime / TL;  // �i�����i0.0�`1.0�j
    temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;  // ���x�̍X�V

    double diff = tmpScore - ansScore;  // �X�R�A�̍���
    double prob = exp(diff / temp);     // �Ă��Ȃ܂��@�̗̍p�m��
    f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];  // ��Ԃ𔽓]
    int upd = 0;  // �𓚂��X�V���邩�ǂ����̃t���O
    if (prob > Rand01()) {
      // �𓚂��X�V����ꍇ�̏���
      // ���D��T���ŘA�������̐����m�F
      rep(i, blockSize + 2) rep(j, blockSize + 2) haba[i][j] = 0;  // ������
      int now = 1;
      upd = 1;
      srep(i, 1, blockSize + 1) {
        srep(j, 1, blockSize + 1) {
          if (haba[i][j] != 0) continue;  // ���ɒT���ς݂Ȃ�X�L�b�v
          if (now == 3) {
            upd = 0;  // �A��������2�𒴂���ꍇ�͍X�V�s��
            break;
          }

          haba[i][j] = now;
          que.push(P(i, j));
          while (que.size()) {
            int x = que.front().first;
            int y = que.front().second;
            que.pop();
            rep(k, 4) {
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
      if (now > 3) upd = 0;  // �A��������3�ȏ�Ȃ�X�V�s��
    }

    if (upd) {
      auto ans2 = ans;  // ���݂̉𓚂��ꎞ�ۑ�

      ans.clear();

      int sx = -1, sy = -1;
      int befx = -1, befy = -1;
      // ���E�̎n�_��T��
      srep(i, 1, blockSize + 1) {
        srep(j, 1, blockSize + 1) {
          if (f[i][j] == 1 && f[i][j - 1] == 0) {
            sx = i;
            sy = j;
            befx = i + 1;
            befy = j;
          }
        }
      }

      if (sx == -1) {
        assert(false);  // �n�_��������Ȃ��ꍇ�̓G���[
      }

      vector<P> vp;  // �|���S���̒��_���i�[
      vp.emplace_back(sx, sy);
      while (true) {
        int x = -1, y = -1;

        // ���E�����ǂ�i����4�������m�F�j
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

        if (x == vp[0].first && y == vp[0].second) break;  // �n�_�ɖ߂����烋�[�v�I��

        if (x == -1) {
          assert(false);  // ���̓_��������Ȃ��ꍇ�̓G���[
          for (auto p : vp) cout << p.first << ' ' << p.second << endl;
          cout << sx << ' ' << sy << ' ' << befx << ' ' << befy << endl;
          srep(i, 1, blockSize + 1) {
            srep(j, 1, blockSize + 1) {
              cout << f[i][j];
            }
            cout << endl;
          }
        }

        befx = sx;
        befy = sy;
        sx = x;
        sy = y;
        vp.emplace_back(sx, sy);  // ���_��ǉ�
      }

      for (auto p : vp) {
        // �O���b�h�̃C���f�b�N�X�������W�ɕϊ�
        int x = (p.first - 1) * (100000 / blockSize);
        int y = (p.second - 1) * (100000 / blockSize);
        ans.emplace_back(x, y);  // �|���S���̒��_�Ƃ��Ēǉ�
      }

      if (IsLengthOK(ans)) {
        upd = 1;  // ����𖞂����Ă���΍X�V
      }
      else {
        upd = 0;  // ����𖞂����Ă��Ȃ���΍X�V���Ȃ�
        ans = ans2;  // ���̉𓚂ɖ߂�
      }
    }

    if (upd) {
      ansScore = tmpScore;  // �X�R�A���X�V
      if (ansScore > best_ansScore) {
        CopyToBest(blockSize);
      }
    }
    else {
      f[rax + 1][ray + 1] = 1 - f[rax + 1][ray + 1];  // ��Ԃ����ɖ߂�
    }
  }

  CopyToAns(blockSize);
}

// ��@�̃��C����������������֐�
int ff[510][510];
void Method3() {
  const int blockSize = 20;  // �O���b�h�̕������i20�~20�̃O���b�h�j

  InitBlock(blockSize);

  int xx1, yy1, xx2, yy2;  // �ŗǂ̋�`�̈�̍��W���i�[����ϐ�

  ansScore = 0;  // �X�R�A�̏�����
  int loop1 = 0;  // ���[�v�񐔂̃J�E���^
  while (true) {
    if (loop1 % 100 == 0) {
      if (GetNowTime() > TL * 0.1) break;  // ���Ԑ����̔������߂����烋�[�v�𔲂���
    }
    loop1++;
    // �����_���ɋ�`�̈��I��
    int x1 = Rand() % blockSize;
    int x2 = Rand() % blockSize;
    int y1 = Rand() % blockSize;
    int y2 = Rand() % blockSize;
    if (x1 > x2) swap(x1, x2);  // x1��x2�����������ɕ��בւ�
    if (y1 > y2) swap(y1, y2);  // y1��y2�����������ɕ��בւ�

    int cnt = 0;  // �I�������̈�̃X�R�A
    srep(i, x1, x2 + 1) {
      srep(j, y1, y2 + 1) {
        cnt += block[i][j];  // �I�������Z���̃X�R�A�����v
      }
    }

    if (cnt > ansScore) {
      // �X�R�A�����P���ꂽ�ꍇ�A�𓚂��X�V
      ansScore = cnt;
      ans.clear();
      // ��`�̎l���̍��W���v�Z���A�𓚂ɒǉ�
      ans.emplace_back(x1 * (100000 / blockSize), y1 * (100000 / blockSize));
      ans.emplace_back((x2 + 1) * (100000 / blockSize), y1 * (100000 / blockSize));
      ans.emplace_back((x2 + 1) * (100000 / blockSize), (y2 + 1) * (100000 / blockSize));
      ans.emplace_back(x1 * (100000 / blockSize), (y2 + 1) * (100000 / blockSize));
      // �ŗǂ̋�`�̈�̃C���f�b�N�X��ۑ�
      xx1 = x1;
      yy1 = y1;
      xx2 = x2;
      yy2 = y2;
    }
  }


  rep(i, blockSize + 2) {
    rep(j, blockSize + 2) {
      f[i][j] = 0;  // ������
    }
  }

  // �ŗǂ̋�`�̈��f�z��ɐݒ�
  srep(i, xx1, xx2 + 1) {
    srep(j, yy1, yy2 + 1) {
      f[i + 1][j + 1] = 1;  // �I�������̈��1�Ƃ���
    }
  }

  CopyToBest(blockSize);

  int loop2 = 0;
  Method3_SA(xx1, xx2, yy1, yy2, blockSize, loop2, TL * 0.75);

  int loop3 = 0;
  int blockSize40 = 40;
  rep(i, blockSize + 2) {
    rep(j, blockSize + 2) {
      ff[i][j] = f[i][j];
    }
  }
  rep(i, blockSize40 + 2) {
    rep(j, blockSize40 + 2) {
      f[i][j] = 0;
    }
  }
  rep(i, blockSize + 2) {
    rep(j, blockSize + 2) {
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
  rep(i, blockSize40 + 2) {
    rep(j, blockSize40 + 2) {
      ff[i][j] = f[i][j];
    }
  }
  rep(i, blockSize80 + 2) {
    rep(j, blockSize80 + 2) {
      f[i][j] = 0;
    }
  }
  rep(i, blockSize40 + 2) {
    rep(j, blockSize40 + 2) {
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
    // �f�o�b�O�p�̏o��
    cout << "loop1 = " << loop1 << ", ";
    cout << "loop2 = " << loop2 << ", ";
    cout << "loop3 = " << loop3 << ", ";
    cout << "loop4 = " << loop4 << ", ";
    cout << endl;
    srep(i, 1, blockSize80 + 1) {
      srep(j, 1, blockSize80 + 1) {
        cout << f[i][j];
      }
      cout << endl;
    }
  }
}

// ���������֐�
ll Solve(int problem_num) {
  ResetTime();  // ���Ԍv���̃��Z�b�g

  SetUp();  // ������Ԃ̏�����

  Input(problem_num);  // ���͂̎󂯎��

  ofstream ofs;
  OpenOfs(problem_num, ofs);  // �o�̓t�@�C���X�g���[���̃I�[�v��

  Method3();  // ��@�̎��s

  Output(ofs);  // �𓚂̏o��

  if (ofs.is_open()) {
    ofs.close();  // �o�̓t�@�C���X�g���[���̃N���[�Y
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();  // �X�R�A�̌v�Z
  }
  return score;  // �X�R�A��Ԃ�
}

/////////////////////////////////////////////////////////////////////////
/*
����

*/
/////////////////////////////////////////////////////////////////////////
int main_old() {
  mode = 2;

  if (mode == 0) {
    Solve(0);  // ���ԍ�0������
  }
  else {
    ll sum = 0;
    srep(i, 0, 100) {
      ll score = Solve(i);  // ���ԍ�i������
      sum += score;         // �X�R�A�̍��v���X�V
      if (mode == 1) {
        cout << score << endl;  // �X�R�A���o��
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
