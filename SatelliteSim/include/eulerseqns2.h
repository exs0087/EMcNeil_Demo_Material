#pragma once

#include <vector>
#include "adcs_hal.h"   // for ActuatorCommands
#include "Vec7.h"

#ifdef __cplusplus
extern "C" {
#endif

    /// Dynamics + control call (RK4 integrates this)
    Vec7 eulerseqns2(double t, const Vec7& y);

    /// Retrieve the one‐per‐step actuator history (after RK4)
    const std::vector<ActuatorCommands>& getActHistory();

    /// **NEW**: Reset the history & call count before a fresh run
    void clearActHistory();

#ifdef __cplusplus
}
#endif
