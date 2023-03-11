from fluids.scenarios import FluidBlock

scenario = FluidBlock(
      density=1e3,
      kinematic_viscosity=1e-6,
      position=(0., 0., 0.),
      dimensions=(1.0, 0.3, 0.3),
      inital_velocity=(0.,0.,0.))

simulation_output = scenario.simulate(device="cpu",
                                      simulation_time=2,
                                      output_time_step=0.03,
                                      cfl_method="no",
                                      particle_radius=0.02)
