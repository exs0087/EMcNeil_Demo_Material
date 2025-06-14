// include/vendor_timer.h
#pragma once

#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif

    /// @returns     A monotonic tick count in milliseconds.
    uint64_t vendor_get_time_ms();

    /// Delay for the given number of milliseconds.
    void vendor_delay_ms(uint64_t ms);

#ifdef __cplusplus
}
#endif
