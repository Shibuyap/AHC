#pragma once

// 二次元座標
struct Point
{
  int x;
  int y;

  Point() : x(0), y(0) {}
  Point(int x_, int y_) : x(x_), y(y_) {}
};

// 長方形
struct Rectangle
{
  // 最小点
  Point min_point;
  // 最大点
  Point max_point;
};

// inner が outer に内包されているかを判定する関数
bool is_rectangle_inside(const Rectangle& inner, const Rectangle& outer)
{
  return (outer.min_point.x <= inner.min_point.x &&
    inner.max_point.x <= outer.max_point.x &&
    outer.min_point.y <= inner.min_point.y &&
    inner.max_point.y <= outer.max_point.y);
}

// 二つの長方形が交差しているかを判定する関数
bool are_rectangles_intersecting(const Rectangle& rect1, const Rectangle& rect2)
{
  if (rect1.max_point.x < rect2.min_point.x || rect2.max_point.x < rect1.min_point.x) {
    return false;
  }
  if (rect1.max_point.y < rect2.min_point.y || rect2.max_point.y < rect1.min_point.y) {
    return false;
  }
  return true;
}
