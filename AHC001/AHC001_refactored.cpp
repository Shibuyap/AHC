#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace chrono;

// Type definitions
typedef long long int ll;
typedef pair<int, int> P;

// Constants
constexpr int MAX_N = 205;
constexpr int GRID_SIZE = 10000;
constexpr double SCORE_MULTIPLIER = 1000000000.0;
constexpr int MOD = 1000000007;
constexpr double REAL_TIME_LIMIT = 4.8;

// Shuffle patterns for rectangle optimization
constexpr int SHUFFLE_PATTERNS[24][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
    {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
    {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
    {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
};

// Structures
struct Point {
    int x;
    int y;
};

struct Rectangle {
    Point p1;  // Top-left corner
    Point p2;  // Bottom-right corner
};

// Random number generator
class RandomGenerator {
private:
    static uint32_t x, y, z, w;
    
public:
    static uint32_t generateInt() {
        uint32_t t = x ^ (x << 11);
        x = y; y = z; z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
    
    static double generateDouble() {
        return (generateInt() + 0.5) * (1.0 / UINT_MAX);
    }
};

uint32_t RandomGenerator::x = 123456789;
uint32_t RandomGenerator::y = 362436069;
uint32_t RandomGenerator::z = 521288629;
uint32_t RandomGenerator::w = 88675123;

// Global state
class ProblemState {
public:
    int n;
    Point target_points[MAX_N];
    int target_sizes[MAX_N];
    Rectangle rectangles[MAX_N];
    int area_sizes[MAX_N];
    Rectangle best_rectangles[MAX_N];
    
    int max_score = -1;
    int real_max_score = -1;
    
    double p[MAX_N];
    double p_sum;
    
    int sort_x[MAX_N], sort_y[MAX_N];
    int arg_sort_x[MAX_N], arg_sort_y[MAX_N];
    
    int mode_count[20] = {0};
    
    // UI-Tei specific storage
    int ui_tei_max_score = -1;
    int ui_tei_a[MAX_N], ui_tei_b[MAX_N], ui_tei_c[MAX_N], ui_tei_d[MAX_N];
    
    // Multi-level best storage
    Rectangle best_best_rectangles[MAX_N];
    int real_real_max_score = -1;
    
    Rectangle best_best_best_rectangles[MAX_N];
    int real_real_real_max_score = -1;
    
    void clear() {
        n = 0;
        max_score = -1;
        real_max_score = -1;
        p_sum = 0;
        real_real_max_score = -1;
        ui_tei_max_score = -1;
        real_real_real_max_score = -1;
        
        for (int i = 0; i < MAX_N; ++i) {
            target_points[i] = {0, 0};
            target_sizes[i] = 0;
            rectangles[i] = {{0, 0}, {0, 0}};
            area_sizes[i] = 0;
            best_rectangles[i] = {{0, 0}, {0, 0}};
            p[i] = 0;
            sort_x[i] = 0;
            sort_y[i] = 0;
            arg_sort_x[i] = 0;
            arg_sort_y[i] = 0;
            best_best_rectangles[i] = {{0, 0}, {0, 0}};
            ui_tei_a[i] = 0;
            ui_tei_b[i] = 0;
            ui_tei_c[i] = 0;
            ui_tei_d[i] = 0;
            best_best_best_rectangles[i] = {{0, 0}, {0, 0}};
        }
    }
};

ProblemState state;

// Utility functions
inline void calculateArea(int idx) {
    state.area_sizes[idx] = (state.rectangles[idx].p2.x - state.rectangles[idx].p1.x) * 
                            (state.rectangles[idx].p2.y - state.rectangles[idx].p1.y);
}

inline void initializeRectangle(Rectangle& rect, const Point& point) {
    rect.p1.x = point.x;
    rect.p1.y = point.y;
    rect.p2.x = point.x + 1;
    rect.p2.y = point.y + 1;
}

inline bool doRectanglesOverlap(int i, int j) {
    int count = 0;
    if (state.rectangles[i].p1.x <= state.rectangles[j].p1.x && state.rectangles[j].p1.x < state.rectangles[i].p2.x) count++;
    else if (state.rectangles[j].p1.x <= state.rectangles[i].p1.x && state.rectangles[i].p1.x < state.rectangles[j].p2.x) count++;
    if (state.rectangles[i].p1.y <= state.rectangles[j].p1.y && state.rectangles[j].p1.y < state.rectangles[i].p2.y) count++;
    else if (state.rectangles[j].p1.y <= state.rectangles[i].p1.y && state.rectangles[i].p1.y < state.rectangles[j].p2.y) count++;
    return count == 2;
}

// Input/Output functions
void loadInput(int fileNum) {
    std::ostringstream oss;
    oss << "./in/" << std::setw(4) << std::setfill('0') << fileNum << ".txt";
    ifstream ifs(oss.str());
    
    if (!ifs.is_open()) {
        cin >> state.n;
        for (int i = 0; i < state.n; ++i) {
            cin >> state.target_points[i].x >> state.target_points[i].y >> state.target_sizes[i];
        }
    } else {
        ifs >> state.n;
        for (int i = 0; i < state.n; ++i) {
            ifs >> state.target_points[i].x >> state.target_points[i].y >> state.target_sizes[i];
        }
    }
}

void writeOutput(int fileNum, const string& suffix = "_out.txt") {
    string fileName = to_string(fileNum) + suffix;
    ofstream ofs(fileName);
    for (int i = 0; i < state.n; ++i) {
        ofs << state.rectangles[i].p1.x << ' ' << state.rectangles[i].p1.y << ' ' 
            << state.rectangles[i].p2.x << ' ' << state.rectangles[i].p2.y << endl;
    }
    ofs.close();
}

// Score calculation
int calculateScore(int ite = -1) {
    if (ite == -1) {
        double sum = 0;
        for (int i = 0; i < state.n; ++i) {
            calculateArea(i);
            double min_size = min(state.target_sizes[i], state.area_sizes[i]);
            double max_size = max(state.target_sizes[i], state.area_sizes[i]);
            double ratio = min_size / max_size;
            state.p[i] = 1.0 - (1.0 - ratio) * (1.0 - ratio);
            sum += state.p[i];
        }
        state.p_sum = sum;
        return round((sum / state.n) * SCORE_MULTIPLIER);
    } else {
        double sum = state.p_sum - state.p[ite];
        calculateArea(ite);
        double min_size = min(state.target_sizes[ite], state.area_sizes[ite]);
        double max_size = max(state.target_sizes[ite], state.area_sizes[ite]);
        double ratio = min_size / max_size;
        state.p[ite] = 1.0 - (1.0 - ratio) * (1.0 - ratio);
        sum += state.p[ite];
        state.p_sum = sum;
        return round((sum / state.n) * SCORE_MULTIPLIER);
    }
}

// Validation functions
bool isRectangleValid(int ite) {
    const Rectangle& rect = state.rectangles[ite];
    
    // Check bounds
    if (rect.p1.x < 0 || rect.p1.x > GRID_SIZE) return false;
    if (rect.p1.y < 0 || rect.p1.y > GRID_SIZE) return false;
    if (rect.p2.x < 0 || rect.p2.x > GRID_SIZE) return false;
    if (rect.p2.y < 0 || rect.p2.y > GRID_SIZE) return false;
    
    // Check valid rectangle
    if (rect.p2.x <= rect.p1.x || rect.p2.y <= rect.p1.y) return false;
    
    // Check if target point is inside
    if (state.target_points[ite].x < rect.p1.x || state.target_points[ite].x >= rect.p2.x) return false;
    if (state.target_points[ite].y < rect.p1.y || state.target_points[ite].y >= rect.p2.y) return false;
    
    return true;
}

bool isConfigurationValid(int ite = -1) {
    if (ite == -1) {
        // Check all rectangles
        for (int i = 0; i < state.n; ++i) {
            if (!isRectangleValid(i)) return false;
        }
        
        // Check overlaps
        for (int i = 0; i < state.n; ++i) {
            for (int j = i + 1; j < state.n; ++j) {
                if (doRectanglesOverlap(i, j)) return false;
            }
        }
    } else {
        // Check single rectangle
        if (!isRectangleValid(ite)) return false;
        
        // Check overlaps with others
        for (int i = 0; i < state.n; ++i) {
            if (i == ite) continue;
            if (doRectanglesOverlap(i, ite)) return false;
        }
    }
    return true;
}

// Sorting initialization
void initializeSorting() {
    vector<P> v;
    
    // Sort by x coordinate
    for (int i = 0; i < state.n; ++i) {
        v.emplace_back(state.target_points[i].x, i);
    }
    sort(v.begin(), v.end());
    for (int i = 0; i < state.n; ++i) {
        state.sort_x[i] = v[i].second;
        state.arg_sort_x[v[i].second] = i;
    }
    
    // Sort by y coordinate
    v.clear();
    for (int i = 0; i < state.n; ++i) {
        v.emplace_back(state.target_points[i].y, i);
    }
    sort(v.begin(), v.end());
    for (int i = 0; i < state.n; ++i) {
        state.sort_y[i] = v[i].second;
        state.arg_sort_y[v[i].second] = i;
    }
}

// State management
void storeBestSolution() {
    state.real_max_score = state.max_score;
    for (int i = 0; i < state.n; ++i) {
        state.best_rectangles[i] = state.rectangles[i];
    }
}

// Optimization methods
void slideRectangle(int ite, double temperature) {
    int diff = 0;
    while (diff == 0) diff = RandomGenerator::generateInt() % 101 - 50;
    int axis = RandomGenerator::generateInt() % 2;
    
    if (axis == 0) {
        state.rectangles[ite].p1.x += diff;
        state.rectangles[ite].p2.x += diff;
    } else {
        state.rectangles[ite].p1.y += diff;
        state.rectangles[ite].p2.y += diff;
    }
    
    if (!isConfigurationValid(ite)) {
        if (axis == 0) {
            state.rectangles[ite].p1.x -= diff;
            state.rectangles[ite].p2.x -= diff;
        } else {
            state.rectangles[ite].p1.y -= diff;
            state.rectangles[ite].p2.y -= diff;
        }
        return;
    }
    
    int new_score = calculateScore(ite);
    int score_diff = new_score - state.max_score;
    double probability = exp(score_diff / temperature);
    
    if (new_score >= state.max_score) {
        state.mode_count[1]++;
        state.max_score = new_score;
        if (state.max_score > state.real_max_score) {
            storeBestSolution();
        }
    } else {
        // Revert changes
        if (axis == 0) {
            state.rectangles[ite].p1.x -= diff;
            state.rectangles[ite].p2.x -= diff;
        } else {
            state.rectangles[ite].p1.y -= diff;
            state.rectangles[ite].p2.y -= diff;
        }
        calculateScore(ite);
    }
}

void modifyOneCorner(int ite, double temperature) {
    int diff = 0;
    while (diff == 0) diff = RandomGenerator::generateInt() % 101 - 50;
    int corner = RandomGenerator::generateInt() % 4;
    
    switch (corner) {
        case 0: state.rectangles[ite].p1.x += diff; break;
        case 1: state.rectangles[ite].p1.y += diff; break;
        case 2: state.rectangles[ite].p2.x += diff; break;
        case 3: state.rectangles[ite].p2.y += diff; break;
    }
    
    if (!isConfigurationValid(ite)) {
        switch (corner) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        return;
    }
    
    int new_score = calculateScore(ite);
    int score_diff = new_score - state.max_score;
    double probability = exp(score_diff / temperature);
    
    if (probability > RandomGenerator::generateDouble()) {
        state.mode_count[0]++;
        state.max_score = new_score;
        if (state.max_score > state.real_max_score) {
            storeBestSolution();
        }
    } else {
        // Revert changes
        switch (corner) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        calculateScore(ite);
    }
}

// Complex expansion logic
Rectangle expandRectangle;
void expandAndShrink(int ite) {
    initializeRectangle(expandRectangle, state.target_points[ite]);
    
    int direction = RandomGenerator::generateInt() % 2;
    
    if (direction == 0) {
        // Expand horizontally first
        expandRectangle.p1.x = 0;
        expandRectangle.p2.x = GRID_SIZE;
        
        // Find constraints from neighboring rectangles
        int argX = state.arg_sort_x[ite];
        
        // Check left neighbors
        for (int ii = argX - 1; ii >= 0; --ii) {
            int i = state.sort_x[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.y <= expandRectangle.p1.y && expandRectangle.p1.y < state.rectangles[i].p2.y) overlaps = true;
            if (expandRectangle.p1.y <= state.rectangles[i].p1.y && state.rectangles[i].p1.y < expandRectangle.p2.y) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].x <= state.target_points[ite].x) {
                    expandRectangle.p1.x = max(expandRectangle.p1.x, state.rectangles[i].p2.x);
                } else {
                    expandRectangle.p2.x = min(expandRectangle.p2.x, state.rectangles[i].p1.x);
                }
                break;
            }
        }
        
        // Check right neighbors
        for (int ii = argX + 1; ii < state.n; ++ii) {
            int i = state.sort_x[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.y <= expandRectangle.p1.y && expandRectangle.p1.y < state.rectangles[i].p2.y) overlaps = true;
            if (expandRectangle.p1.y <= state.rectangles[i].p1.y && state.rectangles[i].p1.y < expandRectangle.p2.y) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].x <= state.target_points[ite].x) {
                    expandRectangle.p1.x = max(expandRectangle.p1.x, state.rectangles[i].p2.x);
                } else {
                    expandRectangle.p2.x = min(expandRectangle.p2.x, state.rectangles[i].p1.x);
                }
                break;
            }
        }
        
        // Similar logic for vertical expansion
        expandRectangle.p1.y = 0;
        expandRectangle.p2.y = GRID_SIZE;
        
        int argY = state.arg_sort_y[ite];
        int nowLeft = expandRectangle.p1.x;
        
        for (int ii = argY - 1; ii >= 0; --ii) {
            int i = state.sort_y[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.x <= expandRectangle.p1.x && expandRectangle.p1.x < state.rectangles[i].p2.x) overlaps = true;
            if (expandRectangle.p1.x <= state.rectangles[i].p1.x && state.rectangles[i].p1.x < expandRectangle.p2.x) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].y <= state.target_points[ite].y) {
                    expandRectangle.p1.y = max(expandRectangle.p1.y, state.rectangles[i].p2.y);
                } else {
                    expandRectangle.p2.y = min(expandRectangle.p2.y, state.rectangles[i].p1.y);
                }
                
                if (state.rectangles[i].p1.x <= nowLeft) {
                    nowLeft = max(nowLeft, state.rectangles[i].p2.x);
                    if (expandRectangle.p2.x <= nowLeft) break;
                }
            }
        }
        
        nowLeft = expandRectangle.p1.x;
        for (int ii = argY + 1; ii < state.n; ++ii) {
            int i = state.sort_y[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.x <= expandRectangle.p1.x && expandRectangle.p1.x < state.rectangles[i].p2.x) overlaps = true;
            if (expandRectangle.p1.x <= state.rectangles[i].p1.x && state.rectangles[i].p1.x < expandRectangle.p2.x) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].y <= state.target_points[ite].y) {
                    expandRectangle.p1.y = max(expandRectangle.p1.y, state.rectangles[i].p2.y);
                } else {
                    expandRectangle.p2.y = min(expandRectangle.p2.y, state.rectangles[i].p1.y);
                }
                
                if (state.rectangles[i].p1.x <= nowLeft) {
                    nowLeft = max(nowLeft, state.rectangles[i].p2.x);
                    if (expandRectangle.p2.x <= nowLeft) break;
                }
            }
        }
    } else {
        // Expand vertically first (similar logic)
        expandRectangle.p1.y = 0;
        expandRectangle.p2.y = GRID_SIZE;
        
        int argY = state.arg_sort_y[ite];
        
        for (int ii = argY - 1; ii >= 0; --ii) {
            int i = state.sort_y[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.x <= expandRectangle.p1.x && expandRectangle.p1.x < state.rectangles[i].p2.x) overlaps = true;
            if (expandRectangle.p1.x <= state.rectangles[i].p1.x && state.rectangles[i].p1.x < expandRectangle.p2.x) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].y <= state.target_points[ite].y) {
                    expandRectangle.p1.y = max(expandRectangle.p1.y, state.rectangles[i].p2.y);
                } else {
                    expandRectangle.p2.y = min(expandRectangle.p2.y, state.rectangles[i].p1.y);
                }
                break;
            }
        }
        
        for (int ii = argY + 1; ii < state.n; ++ii) {
            int i = state.sort_y[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.x <= expandRectangle.p1.x && expandRectangle.p1.x < state.rectangles[i].p2.x) overlaps = true;
            if (expandRectangle.p1.x <= state.rectangles[i].p1.x && state.rectangles[i].p1.x < expandRectangle.p2.x) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].y <= state.target_points[ite].y) {
                    expandRectangle.p1.y = max(expandRectangle.p1.y, state.rectangles[i].p2.y);
                } else {
                    expandRectangle.p2.y = min(expandRectangle.p2.y, state.rectangles[i].p1.y);
                }
                break;
            }
        }
        
        expandRectangle.p1.x = 0;
        expandRectangle.p2.x = GRID_SIZE;
        
        int argX = state.arg_sort_x[ite];
        int nowLeft = expandRectangle.p1.y;
        
        for (int ii = argX - 1; ii >= 0; --ii) {
            int i = state.sort_x[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.y <= expandRectangle.p1.y && expandRectangle.p1.y < state.rectangles[i].p2.y) overlaps = true;
            if (expandRectangle.p1.y <= state.rectangles[i].p1.y && state.rectangles[i].p1.y < expandRectangle.p2.y) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].x <= state.target_points[ite].x) {
                    expandRectangle.p1.x = max(expandRectangle.p1.x, state.rectangles[i].p2.x);
                } else {
                    expandRectangle.p2.x = min(expandRectangle.p2.x, state.rectangles[i].p1.x);
                }
                
                if (state.rectangles[i].p1.y <= nowLeft) {
                    nowLeft = max(nowLeft, state.rectangles[i].p2.y);
                    if (expandRectangle.p2.y <= nowLeft) break;
                }
            }
        }
        
        nowLeft = expandRectangle.p1.y;
        for (int ii = argX + 1; ii < state.n; ++ii) {
            int i = state.sort_x[ii];
            bool overlaps = false;
            
            if (state.rectangles[i].p1.y <= expandRectangle.p1.y && expandRectangle.p1.y < state.rectangles[i].p2.y) overlaps = true;
            if (expandRectangle.p1.y <= state.rectangles[i].p1.y && state.rectangles[i].p1.y < expandRectangle.p2.y) overlaps = true;
            
            if (overlaps) {
                if (state.target_points[i].x <= state.target_points[ite].x) {
                    expandRectangle.p1.x = max(expandRectangle.p1.x, state.rectangles[i].p2.x);
                } else {
                    expandRectangle.p2.x = min(expandRectangle.p2.x, state.rectangles[i].p1.x);
                }
                
                if (state.rectangles[i].p1.y <= nowLeft) {
                    nowLeft = max(nowLeft, state.rectangles[i].p2.y);
                    if (expandRectangle.p2.y <= nowLeft) break;
                }
            }
        }
    }
    
    // Shrink to target size
    int shuf[4];
    int shuffle_seed = RandomGenerator::generateInt() % 24;
    for (int j = 0; j < 4; ++j) {
        shuf[j] = SHUFFLE_PATTERNS[shuffle_seed][j];
    }
    
    int jitter = RandomGenerator::generateInt() % 2;
    
    for (int i = 0; i < 4; ++i) {
        int area = (expandRectangle.p2.x - expandRectangle.p1.x) * (expandRectangle.p2.y - expandRectangle.p1.y);
        if (area <= state.target_sizes[ite]) break;
        
        if (shuf[i] == 0) {
            // Shrink from left
            int target_width = state.target_sizes[ite] / (expandRectangle.p2.y - expandRectangle.p1.y) + jitter;
            int shrink_amount = (expandRectangle.p2.x - expandRectangle.p1.x) - target_width;
            int max_shrink = state.target_points[ite].x - expandRectangle.p1.x;
            expandRectangle.p1.x += min(shrink_amount, max_shrink);
        }
        else if (shuf[i] == 1) {
            // Shrink from top
            int target_height = state.target_sizes[ite] / (expandRectangle.p2.x - expandRectangle.p1.x) + jitter;
            int shrink_amount = (expandRectangle.p2.y - expandRectangle.p1.y) - target_height;
            int max_shrink = state.target_points[ite].y - expandRectangle.p1.y;
            expandRectangle.p1.y += min(shrink_amount, max_shrink);
        }
        else if (shuf[i] == 2) {
            // Shrink from right
            int target_width = state.target_sizes[ite] / (expandRectangle.p2.y - expandRectangle.p1.y) + jitter;
            int shrink_amount = (expandRectangle.p2.x - expandRectangle.p1.x) - target_width;
            int max_shrink = expandRectangle.p2.x - (state.target_points[ite].x + 1);
            expandRectangle.p2.x -= min(shrink_amount, max_shrink);
        }
        else if (shuf[i] == 3) {
            // Shrink from bottom
            int target_height = state.target_sizes[ite] / (expandRectangle.p2.x - expandRectangle.p1.x) + jitter;
            int shrink_amount = (expandRectangle.p2.y - expandRectangle.p1.y) - target_height;
            int max_shrink = expandRectangle.p2.y - (state.target_points[ite].y + 1);
            expandRectangle.p2.y -= min(shrink_amount, max_shrink);
        }
    }
    
    // Validate the result
    bool valid = true;
    if (expandRectangle.p1.x < 0 || expandRectangle.p1.x > GRID_SIZE) valid = false;
    if (expandRectangle.p1.y < 0 || expandRectangle.p1.y > GRID_SIZE) valid = false;
    if (expandRectangle.p2.x < 0 || expandRectangle.p2.x > GRID_SIZE) valid = false;
    if (expandRectangle.p2.y < 0 || expandRectangle.p2.y > GRID_SIZE) valid = false;
    if (expandRectangle.p2.x <= expandRectangle.p1.x) valid = false;
    if (expandRectangle.p2.y <= expandRectangle.p1.y) valid = false;
    if (state.target_points[ite].x < expandRectangle.p1.x || expandRectangle.p2.x <= state.target_points[ite].x) valid = false;
    if (state.target_points[ite].y < expandRectangle.p1.y || expandRectangle.p2.y <= state.target_points[ite].y) valid = false;
    
    if (!valid) {
        initializeRectangle(expandRectangle, state.target_points[ite]);
    }
}

void applyExpansion(int ite, double temperature) {
    Rectangle keep = state.rectangles[ite];
    
    expandAndShrink(ite);
    state.rectangles[ite] = expandRectangle;
    
    int new_score = calculateScore(ite);
    int score_diff = new_score - state.max_score;
    double probability = exp(score_diff / temperature);
    
    if (probability > RandomGenerator::generateDouble()) {
        state.mode_count[4]++;
        state.max_score = new_score;
        if (state.max_score > state.real_max_score) {
            storeBestSolution();
        }
    } else {
        state.rectangles[ite] = keep;
        calculateScore(ite);
    }
}

// Complex optimization with boundary shifting
int overlap_array[MAX_N];
int overlap_count;
int keep_coords[MAX_N][4];

void findOverlappingRectangles(int ite, int direction) {
    overlap_count = 0;
    
    if (direction == 0) { // Left boundary
        int argX = state.arg_sort_x[ite];
        int nowLeft = state.rectangles[ite].p1.y;
        int nowRight = state.rectangles[ite].p2.y;
        
        for (int ii = argX - 1; ii >= 0; --ii) {
            int i = state.sort_x[ii];
            if (doRectanglesOverlap(i, ite)) {
                if (state.rectangles[ite].p1.x <= state.target_points[i].x) {
                    overlap_array[0] = -1;
                    overlap_count = 1;
                    return;
                }
                overlap_array[overlap_count++] = i;
            }
            
            if (state.rectangles[i].p1.y <= nowLeft) {
                nowLeft = max(nowLeft, state.rectangles[i].p2.y);
                if (nowLeft >= nowRight) break;
            }
            if (nowRight <= state.rectangles[i].p2.y) {
                nowRight = min(nowRight, state.rectangles[i].p1.y);
                if (nowLeft >= nowRight) break;
            }
        }
    }
    else if (direction == 1) { // Top boundary
        int argY = state.arg_sort_y[ite];
        int nowLeft = state.rectangles[ite].p1.x;
        int nowRight = state.rectangles[ite].p2.x;
        
        for (int ii = argY - 1; ii >= 0; --ii) {
            int i = state.sort_y[ii];
            if (doRectanglesOverlap(i, ite)) {
                if (state.rectangles[ite].p1.y <= state.target_points[i].y) {
                    overlap_array[0] = -1;
                    overlap_count = 1;
                    return;
                }
                overlap_array[overlap_count++] = i;
            }
            
            if (state.rectangles[i].p1.x <= nowLeft) {
                nowLeft = max(nowLeft, state.rectangles[i].p2.x);
                if (nowLeft >= nowRight) break;
            }
            if (nowRight <= state.rectangles[i].p2.x) {
                nowRight = min(nowRight, state.rectangles[i].p1.x);
                if (nowLeft >= nowRight) break;
            }
        }
    }
    else if (direction == 2) { // Right boundary
        int argX = state.arg_sort_x[ite];
        int nowLeft = state.rectangles[ite].p1.y;
        int nowRight = state.rectangles[ite].p2.y;
        
        for (int ii = argX + 1; ii < state.n; ++ii) {
            int i = state.sort_x[ii];
            if (doRectanglesOverlap(i, ite)) {
                if (state.target_points[i].x < state.rectangles[ite].p2.x) {
                    overlap_array[0] = -1;
                    overlap_count = 1;
                    return;
                }
                overlap_array[overlap_count++] = i;
            }
            
            if (state.rectangles[i].p1.y <= nowLeft) {
                nowLeft = max(nowLeft, state.rectangles[i].p2.y);
                if (nowLeft >= nowRight) break;
            }
            if (nowRight <= state.rectangles[i].p2.y) {
                nowRight = min(nowRight, state.rectangles[i].p1.y);
                if (nowLeft >= nowRight) break;
            }
        }
    }
    else if (direction == 3) { // Bottom boundary
        int argY = state.arg_sort_y[ite];
        int nowLeft = state.rectangles[ite].p1.x;
        int nowRight = state.rectangles[ite].p2.x;
        
        for (int ii = argY + 1; ii < state.n; ++ii) {
            int i = state.sort_y[ii];
            if (doRectanglesOverlap(i, ite)) {
                if (state.target_points[i].y < state.rectangles[ite].p2.y) {
                    overlap_array[0] = -1;
                    overlap_count = 1;
                    return;
                }
                overlap_array[overlap_count++] = i;
            }
            
            if (state.rectangles[i].p1.x <= nowLeft) {
                nowLeft = max(nowLeft, state.rectangles[i].p2.x);
                if (nowLeft >= nowRight) break;
            }
            if (nowRight <= state.rectangles[i].p2.x) {
                nowRight = min(nowRight, state.rectangles[i].p1.x);
                if (nowLeft >= nowRight) break;
            }
        }
    }
}

void shiftBoundaries(int ite, double temperature) {
    int diff = 0;
    while (diff == 0) diff = RandomGenerator::generateInt() % 50 + 1;
    int direction = RandomGenerator::generateInt() % 4;
    
    if (direction < 2) diff *= -1;
    
    // Apply shift to the main rectangle
    switch (direction) {
        case 0: state.rectangles[ite].p1.x += diff; break;
        case 1: state.rectangles[ite].p1.y += diff; break;
        case 2: state.rectangles[ite].p2.x += diff; break;
        case 3: state.rectangles[ite].p2.y += diff; break;
    }
    
    if (!isRectangleValid(ite)) {
        switch (direction) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        return;
    }
    
    findOverlappingRectangles(ite, direction);
    
    if (overlap_count > 0 && overlap_array[0] == -1) {
        switch (direction) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        return;
    }
    
    // Store original positions
    for (int i = 0; i < overlap_count; ++i) {
        keep_coords[i][0] = state.rectangles[overlap_array[i]].p1.x;
        keep_coords[i][1] = state.rectangles[overlap_array[i]].p1.y;
        keep_coords[i][2] = state.rectangles[overlap_array[i]].p2.x;
        keep_coords[i][3] = state.rectangles[overlap_array[i]].p2.y;
    }
    
    // Apply shifts to overlapping rectangles
    bool valid = true;
    for (int i = 0; i < overlap_count; ++i) {
        switch (direction) {
            case 0: state.rectangles[overlap_array[i]].p2.x = state.rectangles[ite].p1.x; break;
            case 1: state.rectangles[overlap_array[i]].p2.y = state.rectangles[ite].p1.y; break;
            case 2: state.rectangles[overlap_array[i]].p1.x = state.rectangles[ite].p2.x; break;
            case 3: state.rectangles[overlap_array[i]].p1.y = state.rectangles[ite].p2.y; break;
        }
        
        if (!isRectangleValid(overlap_array[i])) {
            valid = false;
            break;
        }
    }
    
    if (!valid) {
        // Revert all changes
        for (int i = 0; i < overlap_count; ++i) {
            state.rectangles[overlap_array[i]].p1.x = keep_coords[i][0];
            state.rectangles[overlap_array[i]].p1.y = keep_coords[i][1];
            state.rectangles[overlap_array[i]].p2.x = keep_coords[i][2];
            state.rectangles[overlap_array[i]].p2.y = keep_coords[i][3];
        }
        
        switch (direction) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        return;
    }
    
    // Calculate new score
    for (int i = 0; i < overlap_count; ++i) {
        calculateScore(overlap_array[i]);
    }
    int new_score = calculateScore(ite);
    
    int score_diff = new_score - state.max_score;
    double probability = exp(score_diff / temperature);
    
    if (probability > RandomGenerator::generateDouble()) {
        state.mode_count[5]++;
        state.max_score = new_score;
        if (state.max_score > state.real_max_score) {
            storeBestSolution();
        }
    } else {
        // Revert all changes
        for (int i = 0; i < overlap_count; ++i) {
            state.rectangles[overlap_array[i]].p1.x = keep_coords[i][0];
            state.rectangles[overlap_array[i]].p1.y = keep_coords[i][1];
            state.rectangles[overlap_array[i]].p2.x = keep_coords[i][2];
            state.rectangles[overlap_array[i]].p2.y = keep_coords[i][3];
            calculateScore(overlap_array[i]);
        }
        
        switch (direction) {
            case 0: state.rectangles[ite].p1.x -= diff; break;
            case 1: state.rectangles[ite].p1.y -= diff; break;
            case 2: state.rectangles[ite].p2.x -= diff; break;
            case 3: state.rectangles[ite].p2.y -= diff; break;
        }
        calculateScore(ite);
    }
}

// Main optimization functions
void initializeRectangles() {
    for (int i = 0; i < state.n; ++i) {
        initializeRectangle(state.rectangles[i], state.target_points[i]);
    }
    
    state.max_score = calculateScore(-1);
    if (state.max_score > state.real_max_score) {
        storeBestSolution();
    }
}

// UI-Tei optimization phase
void uiTeiOptimization() {
    clock_t start, end;
    
    for (int ui_tei = 0; ui_tei < 5; ++ui_tei) {
        // Initialize rectangles
        for (int i = 0; i < state.n; ++i) {
            initializeRectangle(state.rectangles[i], state.target_points[i]);
        }
        
        const int T = 5;
        for (int phase = 0; phase < T; ++phase) {
            start = clock();
            
            // Initialize score
            state.max_score = calculateScore(-1);
            storeBestSolution();
            
            // Simulated annealing
            end = clock();
            double now_time = ((double)end - start) / CLOCKS_PER_SEC;
            double time_limit = 0.10 / T;
            double start_temp = 2048;
            double end_temp = 0.1;
            double temp = start_temp + (end_temp - start_temp) * now_time / time_limit;
            int loop = 0;
            
            while (true) {
                loop++;
                if (loop % 100 == 1) {
                    end = clock();
                    now_time = ((double)end - start) / CLOCKS_PER_SEC;
                    if (now_time > time_limit) break;
                    temp = start_temp + (end_temp - start_temp) * now_time / time_limit;
                }
                
                int mode = loop % 4;
                int ite = RandomGenerator::generateInt() % state.n;
                
                if (mode == 0) {
                    modifyOneCorner(ite, temp);
                } else if (mode == 1) {
                    slideRectangle(ite, temp);
                } else if (now_time > 2.0 / T && mode == 2) {
                    // Aspect ratio change (not implemented in refactored version)
                }
            }
            
            // Restore best solution
            state.max_score = state.real_max_score;
            for (int i = 0; i < state.n; ++i) {
                state.rectangles[i] = state.best_rectangles[i];
            }
            calculateScore(-1);
            
            if (state.max_score > state.real_real_max_score) {
                state.real_real_max_score = state.max_score;
                for (int i = 0; i < state.n; ++i) {
                    state.best_best_rectangles[i] = state.rectangles[i];
                }
            }
        }
        
        // Restore best solution from this UI-Tei iteration
        state.max_score = state.real_real_max_score;
        state.real_real_max_score = 0;
        for (int i = 0; i < state.n; ++i) {
            state.rectangles[i] = state.best_best_rectangles[i];
        }
        calculateScore(-1);
        
        // Store UI-Tei results
        if (state.max_score > state.ui_tei_max_score) {
            state.ui_tei_max_score = state.max_score;
            for (int i = 0; i < state.n; ++i) {
                state.ui_tei_a[i] = state.rectangles[i].p1.x;
                state.ui_tei_b[i] = state.rectangles[i].p1.y;
                state.ui_tei_c[i] = state.rectangles[i].p2.x;
                state.ui_tei_d[i] = state.rectangles[i].p2.y;
            }
        }
    }
    
    // Reset state
    state.max_score = 0;
    state.real_max_score = 0;
    state.real_real_max_score = 0;
    for (int i = 0; i < state.n; ++i) {
        initializeRectangle(state.rectangles[i], state.target_points[i]);
        state.best_rectangles[i] = state.rectangles[i];
        state.best_best_rectangles[i] = state.rectangles[i];
    }
}

// Main solve function
int solve(int submission_mode, int file_num) {
    auto start_clock = system_clock::now();
    clock_t real_start = clock();
    
    loadInput(file_num);
    initializeSorting();
    
    // Run UI-Tei optimization
    uiTeiOptimization();
    
    // Restore UI-Tei best solution
    state.max_score = state.ui_tei_max_score;
    for (int i = 0; i < state.n; ++i) {
        state.rectangles[i].p1.x = state.ui_tei_a[i];
        state.rectangles[i].p1.y = state.ui_tei_b[i];
        state.rectangles[i].p2.x = state.ui_tei_c[i];
        state.rectangles[i].p2.y = state.ui_tei_d[i];
    }
    calculateScore(-1);
    
    // Initialize for main optimization
    state.max_score = calculateScore(-1);
    state.real_max_score = state.max_score;
    
    for (int i = 0; i < state.n; ++i) {
        state.best_rectangles[i] = state.rectangles[i];
    }
    
    // Beam search storage
    int a2[100][MAX_N], b2[100][MAX_N], c2[100][MAX_N], d2[100][MAX_N];
    int a4[100][MAX_N], b4[100][MAX_N], c4[100][MAX_N], d4[100][MAX_N];
    int max_score4[100] = {0};
    
    const int beam_width = 1;
    
    // Store initial state
    for (int asai = 0; asai < beam_width; ++asai) {
        for (int j = 0; j < state.n; ++j) {
            a2[asai][j] = state.rectangles[j].p1.x;
            b2[asai][j] = state.rectangles[j].p1.y;
            c2[asai][j] = state.rectangles[j].p2.x;
            d2[asai][j] = state.rectangles[j].p2.y;
        }
    }
    
    // Main optimization loop
    const int T = 250;
    for (int iteration = 0; iteration < T; ++iteration) {
        for (int i = 0; i < 6; ++i) state.mode_count[i] = 0;
        
        const int TT = 1;
        for (int asai = 0; asai < TT; ++asai) {
            int kiyoshi = asai % beam_width;
            
            // Restore state
            for (int i = 0; i < state.n; ++i) {
                state.rectangles[i].p1.x = a2[kiyoshi][i];
                state.rectangles[i].p1.y = b2[kiyoshi][i];
                state.rectangles[i].p2.x = c2[kiyoshi][i];
                state.rectangles[i].p2.y = d2[kiyoshi][i];
            }
            
            // Initialize score
            state.max_score = calculateScore(-1);
            state.real_max_score = state.max_score;
            
            for (int i = 0; i < state.n; ++i) {
                state.best_rectangles[i] = state.rectangles[i];
            }
            
            // Simulated annealing
            clock_t start = clock();
            start_clock = system_clock::now();
            double time_limit = ((REAL_TIME_LIMIT - 0.7) / T) / TT;
            double start_temp = 20048.0;
            double end_temp = 0.1;
            int loop = 0;
            int second_half = iteration % 2;
            
            while (true) {
                loop++;
                if (loop % 100 == 1) {
                    const double time = duration_cast<microseconds>(system_clock::now() - start_clock).count() * 1e-6;
                    if (time > time_limit) break;
                    const double progress_ratio = time / time_limit;
                }
                
                int mode = loop % 6;
                
                if (mode == 1) {
                    int ite = RandomGenerator::generateInt() % state.n;
                    slideRectangle(ite, start_temp);
                }
                else if (mode == 4) {
                    int ite = RandomGenerator::generateInt() % state.n;
                    applyExpansion(ite, start_temp);
                }
                else if (mode == 5) {
                    int ite = RandomGenerator::generateInt() % state.n;
                    shiftBoundaries(ite, start_temp);
                }
                
                // Periodically apply aggressive expansion
                if (loop % 34567 == 1120) {
                    int ite = RandomGenerator::generateInt() % state.n;
                    // ExtendKing functionality would go here
                }
                
                // Recalculate score periodically to avoid numerical errors
                if (loop % 10000 == 1) {
                    state.max_score = calculateScore(-1);
                }
            }
            
            // Restore best solution
            state.max_score = state.real_max_score;
            for (int i = 0; i < state.n; ++i) {
                state.rectangles[i] = state.best_rectangles[i];
            }
            calculateScore(-1);
            
            if (state.max_score > state.real_real_max_score) {
                state.real_real_max_score = state.max_score;
                for (int i = 0; i < state.n; ++i) {
                    state.best_best_rectangles[i] = state.rectangles[i];
                }
            }
            
            // Store for beam search
            max_score4[asai] = state.max_score;
            for (int i = 0; i < state.n; ++i) {
                a4[asai][i] = state.rectangles[i].p1.x;
                b4[asai][i] = state.rectangles[i].p1.y;
                c4[asai][i] = state.rectangles[i].p2.x;
                d4[asai][i] = state.rectangles[i].p2.y;
            }
        }
        
        // Select best for next generation
        vector<P> beam_scores;
        for (int asai = 0; asai < TT; ++asai) {
            beam_scores.emplace_back(max_score4[asai], asai);
        }
        sort(beam_scores.begin(), beam_scores.end(), greater<P>());
        
        // Update next generation
        for (int ii = 0; ii < beam_width; ++ii) {
            int i = beam_scores[ii].second;
            for (int j = 0; j < state.n; ++j) {
                a2[ii][j] = a4[i][j];
                b2[ii][j] = b4[i][j];
                c2[ii][j] = c4[i][j];
                d2[ii][j] = d4[i][j];
            }
        }
        
        // Debug output
        if (submission_mode == 0 && iteration % 10 == 0) {
            cout << "Iteration = " << iteration;
            cout << ", Best score = (" << beam_scores[0].first << ", " << beam_scores[0].second << ")" << endl;
        }
        
        // Check time limit
        clock_t end = clock();
        if (((double)end - real_start) / CLOCKS_PER_SEC > REAL_TIME_LIMIT) break;
    }
    
    // Restore best solution
    state.max_score = state.real_real_max_score;
    for (int i = 0; i < state.n; ++i) {
        state.rectangles[i] = state.best_best_rectangles[i];
    }
    calculateScore(-1);
    
    if (submission_mode == 0) {
        cout << "Final score = " << state.max_score << endl;
    }
    
    // Check for invalid score
    if (submission_mode == 0 && state.max_score > MOD) {
        cout << "ERROR: Score exceeds maximum" << endl;
        writeOutput(file_num, "_out_ERROR.txt");
    }
    
    // Update global best if this is better
    if (state.max_score > state.real_real_real_max_score && state.max_score < MOD) {
        state.real_real_real_max_score = state.max_score;
        for (int i = 0; i < state.n; ++i) {
            state.best_best_best_rectangles[i] = state.rectangles[i];
        }
    }
    
    // Reset state for next iteration
    state.max_score = 0;
    state.real_max_score = 0;
    state.real_real_max_score = 0;
    for (int i = 0; i < state.n; ++i) {
        initializeRectangle(state.rectangles[i], state.target_points[i]);
        state.best_rectangles[i] = state.rectangles[i];
        state.best_best_rectangles[i] = state.rectangles[i];
    }
    
    // Restore global best solution
    state.max_score = state.real_real_real_max_score;
    for (int i = 0; i < state.n; ++i) {
        state.rectangles[i] = state.best_best_best_rectangles[i];
    }
    calculateScore(-1);
    
    // Output final solution
    if (submission_mode) {
        for (int i = 0; i < state.n; ++i) {
            cout << state.rectangles[i].p1.x << ' ' << state.rectangles[i].p1.y << ' '
                 << state.rectangles[i].p2.x << ' ' << state.rectangles[i].p2.y << endl;
        }
    }
    
    writeOutput(file_num);
    
    if (submission_mode == 0) {
        cout << "File No. = " << file_num << ", Final score = " << state.max_score << endl;
    }
    
    if (submission_mode == 0 && state.max_score > MOD) {
        writeOutput(file_num, "_out_ERROR.txt");
    }
    
    return state.max_score;
}

int main() {
    int submission_mode = 0;  // 0 for testing, 1 for submission
    
    if (submission_mode) {
        solve(submission_mode, 0);
    } else {
        int mode = 0;
        
        if (mode == 0) { // Code testing
            solve(submission_mode, 0);
        }
        else if (mode == 1) { // Score verification
            for (int i = 0; i < 1000; ++i) {
                for (int j = 0; j < 50; ++j) {
                    for (int k = 0; k < 10; ++k) {
                        state.clear();
                        solve(submission_mode, j);
                    }
                }
            }
        }
        else if (mode == 2) { // Hyperparameter tuning
            string fileName = "hyperparameter_results.txt";
            ofstream ofs(fileName);
            
            vector<int> beam_widths = {1, 2, 4, 8};
            
            for (int k = 0; k < beam_widths.size(); ++k) {
                for (int l = 0; l < 4; ++l) {
                    // Set hyperparameters based on l
                    
                    ll sum = 0;
                    for (int i = 0; i < 50; ++i) {
                        for (int j = 0; j < 1; ++j) {
                            state.clear();
                            sum += solve(submission_mode, i);
                        }
                    }
                    
                    cout << "Beam width = " << beam_widths[k];
                    cout << ", Configuration = " << l;
                    cout << ", Total score = " << sum << endl;
                    
                    ofs << "Beam width = " << beam_widths[k];
                    ofs << ", Configuration = " << l;
                    ofs << ", Total score = " << sum << endl;
                }
            }
            
            ofs.close();
        }
    }
    
    return 0;
}