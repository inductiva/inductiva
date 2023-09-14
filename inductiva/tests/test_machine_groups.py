"""Test for MachineGroup class."""
import inductiva


def test_machine_group_start(mg: inductiva.resources.MachineGroup):

    mg.start()

    # When the machine group is successfully started, it returns
    # a name of that group. So if the name exists the machine group
    # was created without any errors.
    assert mg.name is not None


def test_get_machine_group(mg: inductiva.resources.MachineGroup):
    machine_group = inductiva.resources.machine_groups.get(mg.name)

    # Check if the parameters of the machine group retrieved from the API are
    # the same as the one created in the test.
    assert machine_group.name == mg.name
    assert machine_group.machine_type == mg.machine_type
    assert machine_group.num_machines == mg.num_machines
    assert machine_group.spot == mg.spot
    assert machine_group.disk_size_gb == mg.disk_size_gb
    assert machine_group.zone == mg.zone


def test_get_all_machine_groups():
    machine_groups = inductiva.resources.machine_groups.get_all()

    assert len(machine_groups) > 0


def test_machine_groups_cost(mg: inductiva.resources.MachineGroup):
    cost = mg.estimate_cloud_cost()

    # Since we don't have a proper way to test if the printed cost is right,
    # here, only if the cost is returned checked.
    assert cost > 0


def test_machine_group_termination(mg: inductiva.resources.MachineGroup):
    mg.terminate()
    machine_groups = inductiva.resources.machine_groups.get_all()

    # If the machine group was successfully terminated, the list of machine
    # groups should be empty.
    assert len(machine_groups) == 0


def test_machine_group():
    """Test the MachienGroup obect methods.

    All the methods to be tested are grouped in this function, so that the
    machine group is only created once and all the tests are ran with it."""
    mg = inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                          num_machines=2,
                                          spot=True,
                                          disk_size_gb=50)

    try:
        test_machine_group_start(mg)
        test_get_machine_group(mg)
        test_get_all_machine_groups()
        test_machine_groups_cost(mg)
        test_machine_group_termination(mg)

    except AssertionError as e:
        # If any error occured, the machine group is terminated in order not to
        # leave any machine group running while debugging.
        mg.terminate()
        raise e
