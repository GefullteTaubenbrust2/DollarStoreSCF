#include "Parallel.hpp"
#include <thread>

namespace flo {
	uint thread_count = 8;

	void setThreadCount(uint _thread_count) {
		thread_count = _thread_count;
	}

	void runInParallel(const std::function<void(int, int)>& function, int problem_size) {
		std::vector<std::thread*> threads;
		threads.resize(thread_count - 1);
		for (int i = 0; i < thread_count - 1; ++i) {
			threads[i] = new std::thread(function, (i + 1) * problem_size / thread_count, (i + 2) * problem_size / thread_count);
		}
		function(0, problem_size / thread_count);
		for (int i = 0; i < thread_count - 1; ++i) {
			threads[i]->join();
			delete threads[i];
		}
	}
}