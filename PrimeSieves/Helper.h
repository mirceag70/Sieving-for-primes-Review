#pragma once

#pragma once

#include <iostream>
#include <iomanip>
#include <queue>
#include <cassert>
#include <array>
#include <thread>
#include <numeric>
#include <fstream>
#include <chrono>
#include <limits>
#include <functional>
#include <typeinfo>

class cTimer
{
private:

	enum class InternalStates { Initial, Iddle, Counting, Pause, Stop };

	InternalStates timer_state = InternalStates::Initial;

	long time_accumulated = 0;
	clock_t mark_start = 0;
	mutable clock_t mark_lap = 0;

	constexpr long TimeDiff2Miliseconds(clock_t const tmDiff) const
	{
		return tmDiff * (1'000 / CLOCKS_PER_SEC);
	}

	void OutputTime(long const tm, std::string_view const msg = "") const
	{
		std::cout << " ( ";
		if (not msg.empty()) std::cout << "[" << msg << "] ";
		std::cout << std::fixed << std::setprecision(0) << tm;
		std::cout << " ms ) ";
	}

public:

	cTimer(void) : timer_state(InternalStates::Iddle) {};

	void Start(void)
	{
		assert(timer_state == InternalStates::Iddle);
		timer_state = InternalStates::Counting;

		// mark down current clock
		mark_lap = mark_start = clock();	// use portable clock()
	}

	void Stop(bool const print_time = false)
	{
		switch (timer_state)
		{
		case InternalStates::Counting:
			time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
			mark_start = 0;
			[[fallthrough]];
		case InternalStates::Pause:
			timer_state = InternalStates::Stop;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(time_accumulated, "StopTime");
	}

	long LapTime(bool print_time = false) const
	{
		assert(timer_state == InternalStates::Counting);

		// get current clock
		clock_t tm = clock();
		long laptime = TimeDiff2Miliseconds(tm - mark_lap);
		//prepare for next lap
		mark_lap = tm;

		if (print_time)
			OutputTime(laptime, "LapTime");

		return laptime;
	}

	long GetTime(bool print_time = false) const
	{
		clock_t tm = 0;

		switch (timer_state)
		{
		case InternalStates::Counting:
			tm = TimeDiff2Miliseconds(tm - mark_start);
			[[fallthrough]];
		case InternalStates::Pause:
		case InternalStates::Stop:
			tm += time_accumulated;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(tm, "Time");

		return tm;
	}

	long Pause(bool print_time = false)
	{
		assert(timer_state == InternalStates::Counting);
		timer_state = InternalStates::Pause;

		time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
		mark_start = 0;

		if (print_time)
			OutputTime(time_accumulated, "PauseTime");

		return time_accumulated;
	}

	void Resume(void)
	{
		assert(timer_state == InternalStates::Pause);
		timer_state = InternalStates::Counting;

		mark_lap = mark_start = clock();
	}

	void Reset(void)
	{
		assert(timer_state == InternalStates::Stop);
		timer_state = InternalStates::Iddle;

		time_accumulated = mark_lap = mark_start = 0;
	}
};

constexpr std::size_t default_iterations = 5;

inline const void nln(bool delimiter = false)
{
	std::cout << std::endl;
	if (delimiter)
		std::cout << "\t--------------------------" << std::endl;
}

template<std::size_t range_length, typename INTEGRAL = int>
consteval auto Range(const int offset = 0)
{
	static_assert(std::is_integral<INTEGRAL>::value,
		"Integral type required for Range.");

	std::array<INTEGRAL, range_length> arr{};
	for (int i = 0; i < range_length; i++)
		arr[i] = i + offset;
	return arr;
}

template<class NUMBER>
constexpr double Average(std::vector<NUMBER> vec)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for Average.");

	return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

inline uint64_t PI_Nmax(uint64_t Nmax)
{
	double logNmax = log(Nmax);
	return (uint64_t)(Nmax / logNmax * (1 + 1.2762 / logNmax));
}

// scalar data type used for prime numbers
//go with 32b for small experiments
#define USE_64_BITS_PRIMES
//we need 64b for larger numbers		
#ifndef USE_64_BITS_PRIMES
typedef		uint32_t		tpPrime;
#else
typedef		uint64_t		tpPrime;
#endif // !USE_64_BITS

class cChecker
{
	std::ifstream file_with_primes;
	tpPrime max_value = 0;

public:
	tpPrime check_next_prime(tpPrime nGenerated)
	{
		if (nGenerated == 2)
			BackToZero();

		tpPrime nRead;
		file_with_primes >> nRead;
		if (nRead <= max_value)
			return nRead == nGenerated;
		else
			return true;
	}

	void BackToZero(void)
	{
		file_with_primes.seekg(0, file_with_primes.beg);
	}

	cChecker(std::string file_name = "E:/list/50MilPrimes.txt", tpPrime maxPrime = 982'451'653) :
		max_value(maxPrime)
	{
		file_with_primes.open(file_name, std::ifstream::in);
		assert(file_with_primes.is_open());
	}

	~cChecker() { file_with_primes.close(); }
};

inline uint8_t idx_last_primes = 0;
inline std::array<uint64_t, 256> last_primes;
inline void AddPrime(uint64_t prime)
{
#ifdef _DEBUG
	static cChecker ckr;
	assert(ckr.check_next_prime(prime));
#endif // _DEBUG

	last_primes[idx_last_primes++] = prime;
}

template<typename INTEGRAL1, typename INTEGRAL2, typename INTEGRAL3,
	std::size_t size1, std::size_t size2 = 0, std::size_t size3 = 0>
	void Try_Sieve(const uint64_t LIMIT, std::string message,
		std::function<uint64_t(uint64_t, INTEGRAL1[], INTEGRAL2[], INTEGRAL3[])> Sieve, bool show_all = true)
{
	static_assert(std::is_integral<INTEGRAL1>::value, "Integral type required for INTEGRAL1.");
	static_assert(std::is_integral<INTEGRAL2>::value, "Integral type required for INTEGRAL2.");
	static_assert(std::is_integral<INTEGRAL3>::value, "Integral type required for INTEGRAL3.");

	INTEGRAL1* v1 = NULL; INTEGRAL2* v2 = NULL; INTEGRAL3* v3 = NULL;

	if (size1 > 0) v1 = new INTEGRAL1[size1];
	if (size2 > 0) v2 = new INTEGRAL2[size2];
	if (size3 > 0) v3 = new INTEGRAL3[size3];

	nln(true);
	cTimer tmr;
	std::vector<decltype(tmr.LapTime())> times;

	tmr.Start();
	std::cout << " - " << message << " - ";
	for (auto i : Range<5>())
	{

		auto numPrimes = Sieve(LIMIT, v1, v2, v3);

		if (i == 0 or show_all)
		{
			nln();
			std::cout << numPrimes << " primes up to " << LIMIT;
			times.push_back(tmr.LapTime(true));
		}
		else
		{
			times.push_back(tmr.LapTime(false));
		}
	}
	nln(true);
	tmr.Stop();
	std::cout << "Average compute time: " << Average(times);

	if (v1) delete[] v1;
	if (v2) delete[] v2;
	if (v3) delete[] v3;
}
