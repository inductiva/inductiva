import inductiva

inductiva.api_url = "http://localhost:7000"
inductiva.api_key = "1234"

insulin = inductiva.molecules.utils.download_pdb_from_rcsb(pdb_id="1ZNI")

scenario = inductiva.molecules.ProteinSolvation(protein_pdb="1ZNI.pdb", temperature=300)

task = scenario.simulate(simulation_time_ns = 0.01)

output = task.get_output(all_files=True)

print(output)