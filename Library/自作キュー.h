#pragma once

// 1次元キュー
int queue_arr[10000];
int queue_head = 0;
int queue_tail = 0;

void clear_queue()
{
  queue_head = 0;
  queue_tail = 0;
}

int front()
{
  return queue_arr[queue_head];
}

void push(int val)
{
  queue_arr[queue_tail] = val;
  queue_tail++;
}

void pop()
{
  queue_head++;
}

int size()
{
  return queue_tail - queue_head;
}

// 2次元キュー
int queue_arr_2d[10000][2];
int queue_head_2d = 0;
int queue_tail_2d = 0;

void clear_queue_2d()
{
  queue_head_2d = 0;
  queue_tail_2d = 0;
}

int front_2d_x()
{
  return queue_arr_2d[queue_head_2d][0];
}

int front_2d_y()
{
  return queue_arr_2d[queue_head_2d][1];
}

void push_2d(int x, int y)
{
  queue_arr_2d[queue_tail_2d][0] = x;
  queue_arr_2d[queue_tail_2d][1] = y;
  queue_tail_2d++;
}

void pop_2d()
{
  queue_head_2d++;
}

int size_2d()
{
  return queue_tail_2d - queue_head_2d;
}
