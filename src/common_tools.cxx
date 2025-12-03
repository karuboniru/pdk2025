#include "common_tools.hxx"
#include <optional>
#include <print>
#include <string>
#include <string_view>
#include <unordered_map>

extern char *const *const environ;

class env_handler {
public:
  static env_handler &get_instance();
  std::optional<std::string> get_env(const std::string &key);

private:
  env_handler();
  std::unordered_map<std::string, std::string> envs{};
  env_handler(const env_handler &other) = delete;
  env_handler &operator=(const env_handler &other) = delete;
};

env_handler::env_handler() {
  for (auto s = environ; *s; s++) {
    std::string str{*s};
    auto pos = str.find('=');
    if (pos != std::string_view::npos) {
      envs[str.substr(0, pos)] = str.substr(pos + 1);
    }
  }
}

std::optional<std::string> env_handler::get_env(const std::string &key) {
  auto it = envs.find(std::string{key});
  if (it != envs.end()) {
    return it->second;
  }
  return std::nullopt;
}

env_handler &env_handler::get_instance() {
  static env_handler instance;
  return instance;
}

size_t guess_nproc_from_env() {
  constexpr auto keys = std::to_array(
      {"NPROC", "OMP_NUM_THREADS", "GSL_NUM_THREADS", "MKL_NUM_THREADS",
       "JULIA_NUM_THREADS", "TF_NUM_THREADS", "GOMAXPROCS",
       "SLURM_CPUS_ON_NODE", "SLURM_NTASKS", "SLURM_NPROCS"});
  for (auto &&key : keys) {
    if (auto val = env_handler::get_instance().get_env(key); val.has_value()) {
      try {
        size_t nproc = std::stoul(val.value());
        if (nproc > 0) {
          return nproc;
        }
      } catch (const std::exception &) {
        // ignore invalid values
        std::println("Warning: invalid value for {}: {}", key, val.value());
        continue;
      }
    }
  }
  // 0 means let ROOT decide
  return 0;
}

double rayleigh_log_likelihood(double value, double sigma) {
  if (value < 0 || sigma <= 0) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::log(value) - ((value * value) / (2 * sigma * sigma));
}

double rayleigh_log_likelihood_normalized(double value, double sigma) {
  return -2. * (rayleigh_log_likelihood(value, sigma) -
                rayleigh_log_likelihood(sigma, sigma));
}

double angle_between(const momentum_t &p4_1, const momentum_t &p4_2) {
  auto dot_product = p4_1.Vect().Dot(p4_2.Vect());
  auto magnitude_product = p4_1.P() * p4_2.P();
  auto cos_angle = dot_product / magnitude_product;
  // Clamp the value to the valid range for acos
  cos_angle = std::clamp(cos_angle, -1.0, 1.0);
  return std::acos(cos_angle);
}
