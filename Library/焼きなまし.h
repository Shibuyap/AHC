#pragma once
#include <functional>

const double TL = 1.8;
double GetNowTime();
double Rand01();

// Ä‚«‚È‚Ü‚µ
void SimulatedAnnealing()
{
  double nowTime = GetNowTime();
  const double START_TEMP = 2048.0;
  const double END_TEMP = 0.1;

  int loop = 0;
  while (true) {
    if (loop % 100 == 0) {
      nowTime = GetNowTime();
      if (nowTime > TL / 2) break;
    }
    loop++;

    double diffScore = 0;

    double progressRatio = nowTime / TL;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
    double prob = exp(diffScore / temp);

    if (prob > Rand01()) {
      // Ì—p
    }
    else {
      // Œ³‚É–ß‚·
    }
  }
}
