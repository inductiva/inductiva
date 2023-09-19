import inductiva

inductiva.api_url = "http://localhost:7000"
inductiva.api_key = "1234"

scenario = inductiva.fluids.FluidTank()

task = scenario.simulate(
)

output=task.get_output()
print(output)