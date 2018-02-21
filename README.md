# Bloch

This is a package of various sequences I have modeled using Bloch simulations.  This simulation is a bit more sophisticated than a typical Bloch simulation (e.g., the Hargreaves MATLAB implementation) as it includes the actual impacts of the RF pulse and slice profile.  As a result, to modify the sequence you will need to modify the actual `c++` code.

To build this package you will need to install both `pybind11` and `armadillo`.  The `pybind11` package is found through Anaconda (`conda install pybind11`).  Armadillo is available through `apt-get` for Ubuntu users, though building it from source is very easy.