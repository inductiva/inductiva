import enum

AVAILABLE_MACHINES = [
    ("c2-standard-", [4, 8, 16, 30, 60]),
    ("c3-standard-", [4, 8, 22, 44, 88, 176]),
    ("c2d-standard-", [2, 4, 8, 16, 32, 56, 112]),
    ("c2d-highcpu-", [2, 4, 8, 16, 32, 56, 112]),
    ("e2-standard-", [2, 4, 8, 16, 32]),
    ("n2-standard-", [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]),
    ("n2d-standard-", [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]),
    ("n1-standard-", [1, 2, 4, 8, 16, 32, 64, 96]),
]

MachineType = enum.Enum("MachineType", AVAILABLE_MACHINES)
