#pragma once
#include <functional>
#include "Types.hpp"

namespace flo {
	void setThreadCount(uint thread_count);

	void runInParallel(const std::function<void(int, int)>& function, int problem_size);
}