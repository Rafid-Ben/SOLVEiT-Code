/* 
 * File:   timer.h
 * Author: Marshall
 *
 * Created on November 26, 2021, 9:18 PM
 */

#ifndef TIMER_H
#define TIMER_H

#include <cstdlib>
#include <cstdio>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#ifndef __USE_BSD
#define __USE_BSD
#endif
#include <sys/time.h>
#endif

class Timer{
    const char* eventName;
    const char* detailName;
    int detailCount;
    int timerPaused = 1;
public:
    void startTimer();
    double getTimer();
    void pauseTimer();
    void printTimeRecorded();
    void resetTimer();
    void changeCount(int n);
    void setEventName(const char* name);
    void setDetailName(const char* name);
    Timer(const char* name = "", const char* detail = "", int n = -1);
private:
    double tempTimeMilliseconds = 0.0;
    double tempTimeSeconds = 0.0;
#ifdef WIN32
    double PCFreq = 0.0;
    __int64 timerStart = 0;
#else
    struct timeval timerStart;
#endif
};

#endif /* TIMER_H */

