#pragma once

// �񎟌����W
struct Point
{
  int x;
  int y;

  Point() : x(0), y(0) {}
  Point(int x_, int y_) : x(x_), y(y_) {}
};

// �����`
struct Rectangle
{
  // �ŏ��_
  Point min_point;
  // �ő�_
  Point max_point;
};

// inner �� outer �ɓ����Ă��邩�𔻒肷��֐�
bool is_rectangle_inside(const Rectangle& inner, const Rectangle& outer)
{
  return (outer.min_point.x <= inner.min_point.x &&
    inner.max_point.x <= outer.max_point.x &&
    outer.min_point.y <= inner.min_point.y &&
    inner.max_point.y <= outer.max_point.y);
}

// ��̒����`���������Ă��邩�𔻒肷��֐�
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
