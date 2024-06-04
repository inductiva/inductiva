
from typing import Any, Dict, List

import inductiva

def query(provider:str ="GCP", query_filter = None) -> List[Dict[str,Any]]:
    
    queried_machines=[]
    
    all_machines = inductiva.resources.machine_types.get_available_machine_types(
        provider)
    
    if query_filter is None:
        return all_machines
    for machine in all_machines:
        if callable(query_filter):
            if query_filter(machine["machine_type"]):
                queried_machines.append(machine)
        else:
            raise ValueError("query_filter must be a callable")
    return queried_machines