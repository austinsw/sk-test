#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

#include <iostream>
#include <memory>
using namespace std;

struct AllocationMetrics {
  uint32_t TotalAllocated = 0;
  uint32_t TotalFreed = 0;
  uint32_t CurrentUsage() { 
    return TotalAllocated - TotalFreed;
  }
};

static AllocationMetrics s_AllocationMetrics;

void* operator new(size_t size) {
  s_AllocationMetrics.TotalAllocated += size;
  return malloc(size);
}

static void PrintMemoryUsage() {
  std::cout << "Memory Usage" << s_AllocationMetrics.CurrentUsage() << "bytes\n";
}

void operator delete(void* memory, size_t size) {
  s_AllocationMetrics.TotalFreed += size;;
  free(memory);
}

class ref_func {
private:
  int i;
  float f;
  char c;
public:
  ref_func(int i, float f, char c) {
    this->i = i;
    this->f = f;
    this->c = c;
  }
  ref_func &seti(int p) {
    i = p;
    return *this;
  }
  ref_func &setf(float q) {
    this->i = i++;
    f = q;
    return *this;
  }
  ref_func &setc(char r) {
    c = r;
    return *this; }
  void disp_val() {
    cout << "The integer value is = " << i << endl;
    cout << "The float value is = " << f << endl;
    cout << "The character value is = " << c << endl;
  }
};

#endif