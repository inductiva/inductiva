import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)


# Initialize the SWASH simulator
swash = inductiva.simulators.SWASH()

task = swash.run(input_dir="swash-resources-example",
                    sim_config_filename="input.sws",
                   )

task.wait()
