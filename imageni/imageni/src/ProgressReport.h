/*
 * ProgressReport.h
 *
 *  Created on: Aug 11, 2012
 *      Author: joohwile
 */

#ifndef PROGRESSREPORT_H_
#define PROGRESSREPORT_H_

class ProgressReport {
public:
	ProgressReport() {
	}
	virtual ~ProgressReport() {
	}
	virtual void report(int iter, double value) = 0;
};


#endif /* PROGRESSREPORT_H_ */
