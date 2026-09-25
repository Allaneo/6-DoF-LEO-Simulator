# 6-DoF LEO Satellite Simulator

During my supervised engineering internship at the Aerospace Technology Centre of the National University of La Plata, I developed a MATLAB/Simulink simulator for the coupled orbital and attitude motion of a satellite in low Earth orbit. Its purpose is to provide a configurable environment for evaluating attitude determination and control designs.

The model includes gravity, atmospheric drag, solar radiation pressure, and geomagnetic effects. Several perturbations offer alternative model choices, so the level of detail can be matched to the analysis. The main Simulink model is [`sim_orbit_min.slx`](Simulink%20test/sim_orbit_min.slx); its scenario and model settings are defined in [`INIT_parametros.m`](Simulink%20test/Inicializaci%C3%B3n/INIT_parametros.m).

## Start here

The project documentation is in Spanish:

- [Quick user guide](Gui%CC%81a%20Ra%CC%81pida%20de%20Usuario%20para%20el%20Simulador%20Orbital.docx): configuration, model options, and representative simulation campaigns.
- [Technical report](Informe%20Simulador%206DoF.docx): the physics, modelling decisions, and simulator architecture.
- [Verification and validation report](Informe%20de%20Verificacio%CC%81n%20y%20Validacio%CC%81n%20del%20Simulador%20Orbital.docx): comparisons with MATLAB Aerospace Toolbox and checks of the simulated perturbations and orbital behaviour.

The source files and analysis scripts are under [`Simulink test`](Simulink%20test/). The reports explain the configuration and the scope of the checks in more detail than this overview.

This is a design and analysis simulator. Its verification checks the implementation and physical consistency of selected cases; it has not been validated as an operational orbit-prediction tool against satellite tracking data.

For a short account of the project and what I learned while building it, see [the companion article](https://olivitoallan.space/posts/6dof-leo-satellite-simulator/).
