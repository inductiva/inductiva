from typing import Any, Dict, List

import inductiva
from inductiva import types

def run (
        name: str,
        input_files: types.Path,
        machines: List[Dict[str, Any]],
        simulator: inductiva.simulators.Simulator,
        input_args: Dict[str, Any] = None,
        machine_args: Dict[str, Any] = None,
        simulator_args: Dict[str, Any] = None
        ):
    
    #TODO
    # Use Projects
    # Add docstring
    # Add quotas check
    # write to csv

    base_input_files=input_files
    
    for machine in machines:

        #Iterate over machine_args and check if any value is callable,
        #if it is call it with machine["machine_type"] as parameter
        for key, value in machine_args.items():
            if callable(value):
                machine_args[key] = value(machine["machine_type"])

        machine_group = inductiva.resources.MachineGroup(machine_type=machine["machine_type"], **machine_args)

        #Iterate over input_args and simulator_args check if any value is callable,
        #if it is call it with machine_group as parameter
        for key, value in input_args.items():
            if callable(value):
                input_args[key] = value(machine_group)
        for key, value in simulator_args.items():
            if callable(value):
                simulator_args[key] = value(machine_group)

        if input_args is not None:
            template_manager = inductiva.TemplateManager(
                template_dir=base_input_files)
            
            template_manager.set_root_dir(name)
            input_files = template_manager.render_dir(**input_args)
            # Deals with templating
            print()
        
        simulator.run(machine_group, input_files, **simulator_args)

    pass