#include "local_rand.h"
#include <random>

TRandom3 &get_thread_local_random() {
  static thread_local TRandom3 rand(std::random_device{}());
  return rand;
}
