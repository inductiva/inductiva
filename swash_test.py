import inductiva

inductiva.api_key = "b7177b3febd6784807c0eed2f38df9ac327eec729e91232087963fa7cf6a7baa"

scenario = inductiva.fluids.scenarios.WindTunnel("motorbike.obj")

scenario.simulate(output_dir="wind_tunnel_output")
