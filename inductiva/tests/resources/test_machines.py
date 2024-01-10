"""Unit tests for the Computational Resources classes."""
import inductiva


def test_machines__mpicluster__register(mocker):
    """Check the registering of a MPICluster.
    
    Goal: Verify that the MPICluster is initializating and the registration is
    processed, by mocking it instead of calling the API.
    """
    mock_register = mocker.patch.object(inductiva.resources.MPICluster,
                                        "_register_machine_group",
                                        return_value=("id-resource",
                                                      "name-resource"))
    cluster = inductiva.resources.MPICluster(machine_type="c2-standard-16",
                                             num_machines=2)

    # Check that the mock was called once
    mock_register.assert_called_once()

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"
    assert cluster.id == "id-resource"
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.type == "mpi"


def test_machines__machine_group__register(mocker):
    """Check the registering of a MachineGroup.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration.
    """
    mock_register = mocker.patch.object(inductiva.resources.MachineGroup,
                                        "_register_machine_group",
                                        return_value=("id-resource",
                                                      "name-resource"))
    cluster = inductiva.resources.MachineGroup(machine_type="c2-standard-16",
                                               num_machines=2)

    # Check that the mock was called once
    mock_register.assert_called_once()

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"
    assert cluster.id == "id-resource"
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.type == "standard"
