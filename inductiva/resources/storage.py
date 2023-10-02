"""Methods to interact with the user storage resources."""
import inductiva
from inductiva.client.apis.tags import instance_api
from inductiva.utils import format_utils


def get_space_used():
    """Returns the occupied storage size in GB."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.get_storage_size()
        storage_used = round(float(response.body), 3)
        print(f"Total storage used: {storage_used} GB")
        return storage_used
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def list_contents(max_results: int = 3,
                  sort_by: str = "size",
                  ascending: bool = False):
    """Lists the contents of the user storage."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        contents = api.list_storage_contents({
            "max_results": max_results,
            "sort_by": sort_by,
            "ascending": ascending
        }).body
        all_contents = []
        for content_name, info in contents.items():
            size = info['size']
            creation_time = info['creation_time']
            all_contents.append({
                'content_name': content_name,
                'size': round(float(size), 3),
                'creation_time': creation_time
            })
        print(_print_contents_table(all_contents))
        return all_contents
    except inductiva.client.ApiException as api_exception:
        raise api_exception


def _print_contents_table(contents):
    columns = ["Name", "Size", "Creation Time"]
    rows = []

    for content in contents:
        row = [
            content['content_name'], content['size'], content['creation_time']
        ]
        rows.append(row)
    override_col_space = {
        "Name": 22,
        "Size": 20,
        "Creation Time": 20,
    }
    formatters = {
        "Creation Time": format_utils.datetime_formatter,
        "Size": format_utils.bytes_formatter
    }

    return format_utils.get_tabular_str(
        rows,
        columns,
        default_col_space=15,
        override_col_space=override_col_space,
        formatters=formatters,
    )
