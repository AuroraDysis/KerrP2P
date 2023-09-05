#pragma once

#include <oneapi/tbb/concurrent_queue.h>

template<typename T>
class ObjectPool {
  std::shared_ptr<tbb::concurrent_bounded_queue<std::shared_ptr<T>>>
      pool;
public:
  explicit ObjectPool()
      : pool(new tbb::concurrent_bounded_queue<std::shared_ptr<T>>()) {
  }

  // Create overloads with different amount of templated parameters.
  std::shared_ptr<T> create() {
    std::shared_ptr<T> obj;
    if (!pool->try_pop(obj))
      obj = std::make_shared<T>();

    // Automatically collects obj.
    return std::shared_ptr<T>(obj.get(), [=](T *) { pool->push(obj); });
  }

  void clear() {
    pool->clear();
  }
};
