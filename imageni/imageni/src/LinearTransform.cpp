/*
 * LinearTransform.cpp
 *
 *  Created on: Jul 29, 2012
 *      Author: joohwile
 */

#include "LinearTransform.h"

namespace cmath {

std::ostream& operator<<(std::ostream& os, const TransformParamType& param) {
	std::cout << "{ ";
	std::cout << param._x[0];
	for (int i = 1; i < 8; i++) {
		if (false && i == 2) {
			std::cout << ", " << param._x[i] * 180 / PI;
		} else {
			std::cout << ", " << param._x[i];
		}
	}
	std::cout << " }";
	return os;
}

}
