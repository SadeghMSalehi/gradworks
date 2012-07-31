/*
 * CClock.h
 *
 *  Created on: Jul 20, 2012
 *      Author: joohwile
 */

#ifndef CCLOCK_H_
#define CCLOCK_H_

#include "time.h"
#include <sys/time.h>

void get_current_time(struct timespec* ts);

/**
 * measure time between tick() and tock()
 */
class CClock {
#ifdef ANDROID
public:
    void tick() {}
    int tock() { return 0; }
#else
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
		get_current_time(&ts1);
		elapsed_time(ts, ts1, ts2);
		return timespec_milliseconds(ts);
	}
#endif /* ANDROID */
};
#endif /* CCLOCK_H_ */
