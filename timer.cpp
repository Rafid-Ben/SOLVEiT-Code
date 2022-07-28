#include "timer.h"

/**
 * Constructor for Timer
 * @param name - name of timer
 * @param detail - additional details about timer
 * @param n - count
 */
Timer::Timer(const char* name, const char* detail, int n){
    eventName = name;
    detailName = detail;
    detailCount = n;
}

/**
 * Set the event name for the timer
 * @param name - event name
 */
void Timer::setEventName(const char* name){
    eventName = name;
}

/**
 * Set the detail name for the timer
 * @param name - detail name
 */
void Timer::setDetailName(const char* name){
    detailName = name;
}

/**
 * Set the count for the timer
 * @param n - count
 */
void Timer::changeCount(int n){
    detailCount = n;
}

/**
 * Reset the timer to zero
 */
void Timer::resetTimer(){
    timerPaused = 1;
    tempTimeSeconds = 0.0;
    tempTimeMilliseconds = 0.0;
}

/**
 * Start the timer
 */
void Timer::startTimer()
{
    timerPaused = 0;
#ifdef WIN32
    LARGE_INTEGER li;
    if (!QueryPerformanceFrequency(&li))
        printf("QueryPerformanceFrequency failed!\n");

    PCFreq = (double)li.QuadPart / 1000.0;

    QueryPerformanceCounter(&li);
    timerStart = li.QuadPart;
#else
    gettimeofday(&timerStart, NULL);
#endif
};

/**
 * Get the time elapsed in ms
 * @return - time elapsed in ms
 */
double Timer::getTimer()
{
#ifdef WIN32
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return (double)(li.QuadPart - timerStart) / PCFreq + tempTimeMilliseconds;
#else
    struct timeval timerStop, timerElapsed;
    gettimeofday(&timerStop, NULL);
    timersub(&timerStop, &timerStart, &timerElapsed);
    return timerElapsed.tv_sec * 1000.0 + timerElapsed.tv_usec / 1000.0 + tempTimeMilliseconds;
#endif
};

/**
 * Pause the timer
 */
void Timer::pauseTimer(){
    if (timerPaused == 0){
#ifdef WIN32
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    tempTimeMilliseconds += (double)(li.QuadPart - timerStart) / PCFreq;
#else
    struct timeval timerStop, timerElapsed;
    gettimeofday(&timerStop, NULL);
    timersub(&timerStop, &timerStart, &timerElapsed);
    tempTimeMilliseconds += timerElapsed.tv_sec * 1000.0 + timerElapsed.tv_usec / 1000.0;
#endif
    }
    else{
        printf("Timer already paused.\n");
    }
    timerPaused = 1;
}

/**
 * Print the time recorded by the timer
 */
void Timer::printTimeRecorded(){
//    if (timerPaused == 0){
//        pauseTimer();
//    }
    if (timerPaused == 0){
        // If the timer isn't paused, 
        tempTimeMilliseconds = getTimer();
        startTimer();
    }
    
    tempTimeSeconds = tempTimeMilliseconds / 1000.0;
    if (detailCount != -1){
        if (detailName[0] != '\0'){
            if (tempTimeSeconds < 60){
                printf("Time recorded for %s %s %d: %.3f seconds.\n",eventName, detailName, detailCount, tempTimeSeconds);
            }
            else if (tempTimeSeconds < 3600){
                int minutes = (int) tempTimeSeconds / 60;
                float seconds = tempTimeSeconds - 60 * minutes;
                printf("Time recorded for %s %s %d: %d minutes, %.3f seconds.\n", eventName, detailName, detailCount, minutes, seconds);
            }
            else{
                int hours = (int) tempTimeSeconds / 3600;
                int minutes = (int) (tempTimeSeconds - 3600 * hours) / 60;
                float seconds = tempTimeSeconds - 3600 * hours - 60 * minutes;
                printf("Time recorded for %s %s %d: %d hours, %d minutes, %.3f seconds.\n",eventName, detailName, detailCount,hours,minutes,seconds);
            }
        }
        else {
            if (tempTimeSeconds < 60){
                printf("Time recorded for %s %d: %.3f seconds.\n",eventName, detailCount, tempTimeSeconds);
            }
            else if (tempTimeSeconds < 3600){
                int minutes = (int) tempTimeSeconds / 60;
                float seconds = tempTimeSeconds - 60 * minutes;
                printf("Time recorded for %s %d: %d minutes, %.3f seconds.\n", eventName, detailCount, minutes, seconds);
            }
            else{
                int hours = (int) tempTimeSeconds / 3600;
                int minutes = (int) (tempTimeSeconds - 3600 * hours) / 60;
                float seconds = tempTimeSeconds - 3600 * hours - 60 * minutes;
                printf("Time recorded for %s %d: %d hours, %d minutes, %.3f seconds.\n",eventName, detailCount,hours,minutes,seconds);
            }
        }
    }
    else if (detailName[0] != '\0'){
        if (tempTimeSeconds < 60){
            printf("Time recorded for %s %s: %.3f seconds.\n",eventName, detailName, tempTimeSeconds);
        }
        else if (tempTimeSeconds < 3600){
            int minutes = (int) tempTimeSeconds / 60;
            float seconds = tempTimeSeconds - 60 * minutes;
            printf("Time recorded for %s %s: %d minutes, %.3f seconds.\n", eventName, detailName, minutes, seconds);
        }
        else{
            int hours = (int) tempTimeSeconds / 3600;
            int minutes = (int) (tempTimeSeconds - 3600 * hours) / 60;
            float seconds = tempTimeSeconds - 3600 * hours - 60 * minutes;
            printf("Time recorded for %s %s: %d hours, %d minutes, %.3f seconds.\n",eventName, detailName,hours,minutes,seconds);
        }
    }
    else if (eventName[0] != '\0'){
        if (tempTimeSeconds < 60){
            printf("Time recorded for %s: %.3f seconds.\n",eventName, tempTimeSeconds);
        }
        else if (tempTimeSeconds < 3600){
            int minutes = (int) tempTimeSeconds / 60;
            float seconds = tempTimeSeconds - 60 * minutes;
            printf("Time recorded for %s: %d minutes, %.3f seconds.\n", eventName, minutes, seconds);
        }
        else{
            int hours = (int) tempTimeSeconds / 3600;
            int minutes = (int) (tempTimeSeconds - 3600 * hours) / 60;
            float seconds = tempTimeSeconds - 3600 * hours - 60 * minutes;
            printf("Time recorded for %s: %d hours, %d minutes, %.3f seconds.\n",eventName,hours,minutes,seconds);
        }
    }
    else{
        if (tempTimeSeconds < 60){
            printf("Time recorded: %.3f seconds.\n", tempTimeSeconds);
        }
        else if (tempTimeSeconds < 3600){
            int minutes = (int) tempTimeSeconds / 60;
            float seconds = tempTimeSeconds - 60 * minutes;
            printf("Time recorded: %d minutes, %.3f seconds.\n", minutes, seconds);
        }
        else{
            int hours = (int) tempTimeSeconds / 3600;
            int minutes = (int) (tempTimeSeconds - 3600 * hours) / 60;
            float seconds = tempTimeSeconds - 3600 * hours - 60 * minutes;
            printf("Time recorded: %d hours, %d minutes, %.3f seconds.\n",hours,minutes,seconds);
        }
    }
}