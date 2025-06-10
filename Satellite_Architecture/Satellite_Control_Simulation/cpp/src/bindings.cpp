#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "control_sim.hpp"

namespace py = pybind11;

PYBIND11_MODULE(satellite_kernels, m) {
    m.doc() = "Satellite attitude control simulation";

    m.def("simulate_control",
        [](double t0, double tf, double dt,
           std::array<double,3> w0,
           std::array<double,4> q0,
           bool k_bdot)
        {
            State x0{w0,q0};
            History H = simulateControl(t0, tf, dt, x0, k_bdot);

            py::list py_states;
            for(auto &S: H.states) {
                py::dict sd;
                sd["omega"] = S.omega;
                sd["quat"]  = S.quat;
                py_states.append(sd);
            }

            py::list py_torques;
            for(auto &L: H.control_torques)
                py_torques.append(L);

            py::dict result;
            result["states"]          = py_states;
            result["control_torques"] = py_torques;
            return result;
        },
        py::arg("t0"), py::arg("tf"), py::arg("dt"),
        py::arg("w0"), py::arg("q0"), py::arg("k_bdot")
    );
}
