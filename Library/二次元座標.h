#pragma once

// �񎟌����W
struct Point
{
public:
  int x;
  int y;

  Point() { x = 0; y = 0; }
  Point(int _x, int _y) { x = _x; y = _y; }
};

// �����`
struct Rectangle
{
  Point p1; // �ŏ��_
  Point p2; // �ő�_
};

// inner��outer�ɓ����Ă��邩�𔻒肷��֐�
bool IsRectangleInside(const Rectangle& inner, const Rectangle& outer)
{
  return outer.p1.x <= inner.p1.x && inner.p2.x <= outer.p2.x && outer.p1.y <= inner.p1.y && inner.p2.y <= outer.p2.y;
}

// ��̒����`���������Ă��邩�𔻒肷��֐�
bool AreRectanglesIntersecting(const Rectangle& rect1, const Rectangle& rect2)
{
  if (rect1.p2.x < rect2.p1.x || rect2.p2.x < rect1.p1.x) return false;
  if (rect1.p2.y < rect2.p1.y || rect2.p2.y < rect1.p1.y) return false;
  return true;
}
