"""Available machine types and their number of cores."""


def list_available_machines():
    """Lists the available machines from the descriptive dict.
    
    Each machine type is a key in the dict with respective possible values are:
    - vcpus: list of integers with the available vcpu possibilities.
    - memory: list of string with the available memory types. Default is
        ["highmem", "standard", "highcpu"].
    - extra_configs: tuple with two lists of strings and integers. The first
        list is extra memory types that are the only ones containing the
        possible vcpus in the second list.
    - lssd: boolean indicating if the machine type has the local SSD
        option. Default is False.
    """

    available_machines = []

    for machine_type, configs in AVAILABLE_MACHINES.items():
        vcpu_list = configs["vcpus"]
        memory_list = configs.get("memory", ["highmem", "standard", "highcpu"])
        lssd = configs.get("lssd", False)
        extra_configs = configs.get("extra_configs", None)

        for memory in memory_list:
            for vcpu in vcpu_list:
                available_machines.append(machine_type + "-" + memory + "-" +
                                          str(vcpu))

        if extra_configs:
            for memory in extra_configs[0]:
                for vcpu in extra_configs[1]:
                    available_machines.append(machine_type + "-" + memory +
                                              "-" + str(vcpu))
        if lssd:
            for vcpu in vcpu_list[1:]:
                available_machines.append(machine_type + "-standard-" +
                                          str(vcpu) + "-lssd")

    return tuple(available_machines)


AVAILABLE_MACHINES = {
    "c2": {
        "vcpus": [4, 8, 16, 30, 60],
        "memory": ["standard"]
    },
    "c3": {
        "vcpus": [4, 8, 22, 44, 88, 176],
        "lssd": True
    },
    "c2d": {
        "vcpus": [2, 4, 8, 16, 32, 56, 112]
    },
    "c3d": {
        "vcpus": [4, 8, 16, 30, 60, 90, 180, 360],
        "lssd": True
    },
    "e2": {
        "vcpus": [2, 4, 8, 16],
        "extra_configs": (["standard", "highcpu"], [32])
    },
    "n2": {
        "vcpus": [2, 4, 8, 16, 32, 48, 64, 80, 96],
        "extra_configs": (["standard", "highmem"], [128])
    },
    "n2d": {
        "vcpus": [2, 4, 8, 16, 32, 48, 64, 80, 96],
        "extra_configs": (["standard", "highcpu"], [128, 224])
    },
    "n1": {
        "vcpus": [1, 2, 4, 8, 16, 32, 64, 96]
    }
}
