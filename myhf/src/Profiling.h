#pragma once


#if PROFILING_ENABLED

#include <chrono>
#include <fstream>
#include <thread>

#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

#define PROFILE_BEGIN_SESSION(name, filepath) ::myhf::Profiler::Get().BeginSession(name, filepath)
#define PROFILE_END_SESSION() ::myhf::Profiler::Get().EndSession()
#define PROFILE_SCOPE(name) ::myhf::ProfilerTimer TOKENPASTE2(timer, __LINE__)(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCTION__)

#else
#define PROFILE_BEGIN_SESSION(name, filepath)
#define PROFILE_END_SESSION()
#define PROFILE_SCOPE(name)
#define PROFILE_FUNCTION()
#endif

#if PROFILING_ENABLED
namespace myhf
{
struct ProfileResult
{
	std::string name;
	long long start;
	long long end;
	uint32_t threadID;
};

class Profiler
{
public:
	void BeginSession(std::string_view name, const std::string& outputFilename = "results.json")
	{
		m_outputStream.open(outputFilename);
		WriteHeader();
		m_sessionName = name;
	}

	void EndSession()
	{
		WriteFooter();
		m_outputStream.close();
		m_profileCount = 0;
	}

	void WriteProfile(const ProfileResult& result)
	{
		if (m_profileCount++ > 0)
			m_outputStream << ',';

		std::string name = result.name;
		std::replace(name.begin(), name.end(), '"', '\'');

		m_outputStream << "{";
		m_outputStream << "\"cat\":\"function\",";
		m_outputStream << "\"dur\":" << (result.end - result.start) << ',';
		m_outputStream << "\"name\":\"" << name << "\",";
		m_outputStream << "\"ph\":\"X\",";
		m_outputStream << "\"pid\":0,";
		m_outputStream << "\"tid\":" << result.threadID << ',';
		m_outputStream << "\"ts\":" << result.start;
		m_outputStream << "}";

		m_outputStream.flush();
	}

	void WriteHeader()
	{
		m_outputStream << "{\"otherData\": {},\"traceEvents\":[";
		m_outputStream.flush();
	}

	void WriteFooter()
	{
		m_outputStream << "]}";
		m_outputStream.flush();
	}

	static Profiler& Get()
	{
		static Profiler profiler;
		return profiler;
	}

private:
	std::string m_sessionName{};
	std::ofstream m_outputStream{};
	int m_profileCount{ 0 };
};

class ProfilerTimer
{
public:
	ProfilerTimer(std::string_view name) :
		m_name(name), 
		m_stopped(false), 
		m_startTimePoint(std::chrono::high_resolution_clock::now())
	{}
	~ProfilerTimer()
	{
		if (!m_stopped)
			Stop();
	}

	void Stop()
	{
		auto endTimePoint = std::chrono::high_resolution_clock::now();
		long long start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTimePoint).time_since_epoch().count();
		long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimePoint).time_since_epoch().count();

		uint32_t threadID = std::hash<std::thread::id>{}(std::this_thread::get_id());
		Profiler::Get().WriteProfile({ m_name, start, end, threadID });

		m_stopped = true;
	}

private:
	std::string m_name;
	bool m_stopped;
	std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimePoint;
};


}
#endif