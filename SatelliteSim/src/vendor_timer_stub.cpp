// stubs/vendor_timer_stub.cpp

#include "vendor_timer.h"
#include <chrono>
#include <thread>

/// Returns a monotonic “tick” count in milliseconds.
/// In real firmware this would read a hardware timer register.
uint64_t vendor_get_time_ms() {
    using namespace std::chrono;
    auto now = steady_clock::now().time_since_epoch();
    return static_cast<uint64_t>(duration_cast<milliseconds>(now).count());
}

/// Busy‐wait or sleep for the specified number of milliseconds.
/// In real firmware this might spin on a timer or invoke an RTOS delay.
void vendor_delay_ms(uint64_t ms) {
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}
