#include "timer.h"
#include <chrono>
#include <print>
#include <unordered_map>

class TimerSummary {
public:
  static TimerSummary &instance() {
    static TimerSummary instance;
    return instance;
  }

  void add_time(const std::string &name, double milliseconds) {
    m_times[name] += milliseconds;
  }

  void print_summary() const {
    std::println("=== Timer Summary ===");
    for (const auto &[name, time] : m_times) {
      std::println("  '{}' : {:.3f} ms", name, time);
    }
  }

  TimerSummary(const TimerSummary &) = delete;
  TimerSummary &operator=(const TimerSummary &) = delete;
  TimerSummary(TimerSummary &&) = delete;
  TimerSummary &operator=(TimerSummary &&) = delete;

private:
  TimerSummary() = default;
  ~TimerSummary() { print_summary(); }

  std::unordered_map<std::string, double> m_times;
};

ScopedTimer::ScopedTimer(std::string_view name)
    : m_name(name), m_start(std::chrono::steady_clock::now()) {}

ScopedTimer::~ScopedTimer() {
  const auto end = std::chrono::steady_clock::now();
  const std::chrono::duration<double, std::milli> elapsed = end - m_start;
  TimerSummary::instance().add_time(m_name, elapsed.count());
  std::println("[Timer] '{}' took {:.3f} ms", m_name, elapsed.count());
}
