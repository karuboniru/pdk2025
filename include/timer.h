#pragma once

#include <chrono>
#include <string>
#include <string_view>

class ScopedTimer {
public:
  explicit ScopedTimer(std::string_view name);
  ~ScopedTimer();

  ScopedTimer(const ScopedTimer &) = delete;
  ScopedTimer &operator=(const ScopedTimer &) = delete;
  ScopedTimer(ScopedTimer &&) = delete;
  ScopedTimer &operator=(ScopedTimer &&) = delete;

private:
  std::string m_name;
  std::chrono::steady_clock::time_point m_start;
};