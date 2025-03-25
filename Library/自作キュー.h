#pragma once

// 1ŽŸŒ³
int queueArr[10000];
int queueHead = 0;
int queueTail = 0;
void ClearQueue()
{
  queueHead = 0;
  queueTail = 0;
}
int Front()
{
  return queueArr[queueHead];
}
void Push(int val)
{
  queueArr[queueTail] = val;
  queueTail++;
}
void Pop()
{
  queueHead++;
}
int Size()
{
  return queueTail - queueHead;
}

// 2ŽŸŒ³
int queueArr2[10000][2];
int queueHead2 = 0;
int queueTail2 = 0;
void ClearQueue()
{
  queueHead2 = 0;
  queueTail2 = 0;
}
int FrontX()
{
  return queueArr2[queueHead2][0];
}
int FrontY()
{
  return queueArr2[queueHead2][1];
}
void Push(int x, int y)
{
  queueArr2[queueTail2][0] = x;
  queueArr2[queueTail2][1] = y;
  queueTail2++;
}
void Pop()
{
  queueHead2++;
}
int Size()
{
  return queueTail2 - queueHead2;
}
