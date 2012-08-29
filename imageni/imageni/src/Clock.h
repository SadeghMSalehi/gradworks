/*
 * CClock.h
 *
 *  Created on: Jul 20, 2012
 *      Author: joohwile
 */

#ifndef CCLOCK_H_
#define CCLOCK_H_

#include "config.h"
#include "time.h"

#ifdef __MACH__
#include <sys/time.h>
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void get_current_time(struct timespec* ts);

/**
 * measure time between tick() and tock()
 */
#ifdef ANDROID
class CClock {
public:
	void tick() {}
	int tock() {return 0;}
};
#endif
#ifdef _MACH_
void get_current_time(struct timespec* ts) {
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	ts->tv_sec = mts.tv_sec;
	ts->tv_nsec = mts.tv_nsec;
}
#endif
#ifdef WIN32

typedef union _LARGE_INTEGER {
  struct {
    unsigned long LowPart;
    long  HighPart;
  };
  struct {
	unsigned long LowPart;
    long  HighPart;
  } u;
  __int64 QuadPart;
} LARGE_INTEGER, *PLARGE_INTEGER;

extern "C" {
__declspec(dllimport)
int
__stdcall
QueryPerformanceCounter(LARGE_INTEGER *lpPerformanceCount);

__declspec(dllimport)
int
__stdcall
QueryPerformanceFrequency(LARGE_INTEGER *lpFrequency);
}

class CClock {
private:
	LARGE_INTEGER frequency;
	LARGE_INTEGER start, end;

public:
	CClock() {
		if (::QueryPerformanceFrequency(&frequency) == false) {
			throw "Food";
		}
	}

	void tick() {
		if (::QueryPerformanceCounter(&start) == false) {
			throw "Food";
		}
	}
// Calculation.

	unsigned int tock() {
		if (::QueryPerformanceCounter(&end) == false) {
			throw "foo";
		}

		return static_cast<unsigned int>((end.QuadPart - start.QuadPart) * 1000
				/ frequency.QuadPart);
	}
};
#endif
#ifdef CYGWIN
class CClock {
private:
	struct timespec ts1;
	struct timespec ts2;
	struct timespec ts;
public:
	static void elapsed_time(struct timespec &ts, struct timespec &ts1,
			struct timespec &ts2) {
		const int NANO_SECONDS_IN_SEC = 1000000000;
		ts.tv_sec = ts1.tv_sec - ts2.tv_sec;
		ts.tv_nsec = ts1.tv_nsec - ts2.tv_nsec;
		if (ts.tv_nsec < 0) {
			ts.tv_sec--;
			ts.tv_nsec += NANO_SECONDS_IN_SEC;
		}
	}

	static int timespec_milliseconds(struct timespec &a) {
		return a.tv_sec * 1000 + a.tv_nsec / 1000000;
	}

	void tick() {
		get_current_time(&ts2);
	}

	int tock() {
		clock_gettime(CLOCK_MONOTONIC, &ts1);
		elapsed_time(ts, ts1, ts2);
		return timespec_milliseconds(ts);
	}
};
#endif /* ANDROID */
#endif /* CCLOCK_H_ */
