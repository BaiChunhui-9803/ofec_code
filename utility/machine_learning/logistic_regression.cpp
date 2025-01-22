#include "logistic_regression.h"

namespace ofec {
	void check_it(bool flag) {
		if (!flag) fprintf(stderr, "logistic regression error! \n");
	}
}