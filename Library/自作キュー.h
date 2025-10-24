#pragma once

// 1次元キューのクラス
class Queue1D
{
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue1D() : head(0), tail(0) {}

  void clear_queue()
  {
    head = 0;
    tail = 0;
  }

  int front() const
  {
    return arr[head];
  }

  void push(int val)
  {
    arr[tail] = val;
    tail++;
  }

  void pop()
  {
    head++;
  }

  int size() const
  {
    return tail - head;
  }

  bool empty() const
  {
    return head == tail;
  }
};

Queue1D que;

// 2次元キューのクラス
class Queue2D
{
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE][2];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue2D() : head(0), tail(0) {}

  void clear_queue()
  {
    head = 0;
    tail = 0;
  }

  int front_x() const
  {
    return arr[head][0];
  }

  int front_y() const
  {
    return arr[head][1];
  }

  void push(int x, int y)
  {
    arr[tail][0] = x;
    arr[tail][1] = y;
    tail++;
  }

  void pop()
  {
    head++;
  }

  int size() const
  {
    return tail - head;
  }

  bool empty() const
  {
    return head == tail;
  }
};

Queue2D que;
