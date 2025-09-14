#pragma once
#include <chrono>
#include "Types.hpp"

namespace flo {
	/// <summary> Units of time measurement </summary>
	enum TimeUnit {
		nanosecond = 0,
		microsecond = 1,
		millisecond = 2,
		second = 3
	};

	/// <summary> A basic struct for measuring time
	struct Time {
	protected:
		u64 nanoseconds = 0;

	public:
		Time() = default;

		/// <summary> Convert an amount of time to a Time object </summary>
		/// <param name="amount"> The amount of units to convert to a Time object </param>
		/// <param name="unit"> The time unit to be used </param>
		Time(double amount, TimeUnit unit) {
			switch (unit) {
			case nanosecond :
				nanoseconds = amount;
				break;
			case microsecond :
				nanoseconds = amount * 1000;
				break;
			case millisecond :
				nanoseconds = amount * 1000000;
				break;
			case second :
				nanoseconds = amount * 1000000000;
				break;
			}
		}

		/// <summary> Convert an amount of time to a Time object </summary>
		/// <param name="amount"> The amount of units to convert to a Time object </param>
		Time(u64 amount) {
			nanoseconds = amount;
		}

		/// <summary> Convert seconds to a Time object </summary>
		/// <param name="seconds The amount of seconds
		explicit Time(double seconds) : nanoseconds(seconds * 1000000000) {}

		double asSeconds() {
			return (double)nanoseconds / 1000000000.;
		}

		double asMilliseconds() {
			return (double)nanoseconds / 1000000.;
		}

		double asMicroseconds() {
			return (double)nanoseconds / 1000.;
		}
		
		double asNanoSeconds() {
			return (double)nanoseconds;
		}

		Time operator+(const Time time) {
			return Time(static_cast<u64>(nanoseconds + time.nanoseconds));
		}

		Time operator-(const Time time) {
			return Time(static_cast<u64>(nanoseconds - time.nanoseconds));
		}

		Time operator*(const double right) {
			return Time(static_cast<u64>(nanoseconds * right));
		}

		Time operator/(const double right) {
			return Time(static_cast<u64>(nanoseconds / right));
		}

		bool operator==(const Time right) {
			return nanoseconds == right.nanoseconds;
		}

		bool operator>=(const Time right) {
			return nanoseconds >= right.nanoseconds;
		}

		bool operator<=(const Time right) {
			return nanoseconds <= right.nanoseconds;
		}

		bool operator>(const Time right) {
			return nanoseconds > right.nanoseconds;
		}

		bool operator<(const Time right) {
			return nanoseconds < right.nanoseconds;
		}

		/// <summary> Convert a Time object to seconds by casting </summary>
		operator double() const {
			return (double)nanoseconds / 1000000000.;
		}
	};

	/// <summary> A simple struct for measuring time intervals </summary>
	/// \see Time
	struct Stopclock {
	protected:
		std::chrono::time_point<std::chrono::high_resolution_clock> start;

	public:
		/// <summary> Construct and start the stopclock </summary>
		Stopclock() {
			start = std::chrono::high_resolution_clock::now();
		}

		/// <summary> If this function is called, the stopclock measures its time from that point on </summary>
		void reset() {
			start = std::chrono::high_resolution_clock::now();
		}

		/// <summary> This function does not actually stop the stopclock, but measures the time since it was started or last reset </summary>
		/// <returns> The time since it was started or last reset </returns>
		/// \see reset
		Time stop() {
			return Time(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count(), nanosecond);
		}
	};
}