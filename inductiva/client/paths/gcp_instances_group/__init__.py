# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from inductiva.client.paths.gcp_instances_group import Api

from inductiva.client.paths import PathValues

path = PathValues.GCP_INSTANCES_GROUP