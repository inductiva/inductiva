# coding: utf-8

"""
    InductivaWebAPI

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501

from __future__ import annotations
import pprint
import re  # noqa: F401
import json

from pydantic import BaseModel, ConfigDict, StrictBool, StrictFloat, StrictInt, StrictStr
from typing import Any, ClassVar, Dict, List, Optional, Union
from inductiva.client.models.compression_method import CompressionMethod
from typing import Optional, Set
from typing_extensions import Self


class TaskRequest(BaseModel):
    """
    TaskRequest
    """

  # noqa: E501
    simulator: StrictStr
    simulator_name_alias: Optional[StrictStr] = None
    resource_pool: StrictStr
    storage_path_prefix: Optional[StrictStr] = ''
    project: Optional[StrictStr] = None
    container_image: StrictStr
    time_to_live_seconds: Optional[Union[StrictFloat, StrictInt]] = None
    resubmit_on_preemption: Optional[StrictBool] = False
    input_resources: Optional[List[StrictStr]] = None
    stream_zip: Optional[StrictBool] = True
    compress_with: Optional[CompressionMethod] = None
    extra_params: Optional[Dict[str, Any]] = None
    __properties: ClassVar[List[str]] = [
        "simulator", "simulator_name_alias", "resource_pool",
        "storage_path_prefix", "project", "container_image",
        "time_to_live_seconds", "resubmit_on_preemption", "input_resources",
        "stream_zip", "compress_with", "extra_params"
    ]

    model_config = ConfigDict(
        populate_by_name=True,
        validate_assignment=True,
        protected_namespaces=(),
    )

    def to_str(self) -> str:
        """Returns the string representation of the model using alias"""
        return pprint.pformat(self.model_dump(by_alias=True))

    def to_json(self) -> str:
        """Returns the JSON representation of the model using alias"""
        # TODO: pydantic v2: use .model_dump_json(by_alias=True, exclude_unset=True) instead
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_str: str) -> Optional[Self]:
        """Create an instance of TaskRequest from a JSON string"""
        return cls.from_dict(json.loads(json_str))

    def to_dict(self) -> Dict[str, Any]:
        """Return the dictionary representation of the model using alias.

        This has the following differences from calling pydantic's
        `self.model_dump(by_alias=True)`:

        * `None` is only added to the output dict for nullable fields that
          were set at model initialization. Other fields with value `None`
          are ignored.
        """
        excluded_fields: Set[str] = set([])

        _dict = self.model_dump(
            by_alias=True,
            exclude=excluded_fields,
            exclude_none=True,
        )
        # set to None if simulator_name_alias (nullable) is None
        # and model_fields_set contains the field
        if self.simulator_name_alias is None and "simulator_name_alias" in self.model_fields_set:
            _dict['simulator_name_alias'] = None

        # set to None if project (nullable) is None
        # and model_fields_set contains the field
        if self.project is None and "project" in self.model_fields_set:
            _dict['project'] = None

        # set to None if time_to_live_seconds (nullable) is None
        # and model_fields_set contains the field
        if self.time_to_live_seconds is None and "time_to_live_seconds" in self.model_fields_set:
            _dict['time_to_live_seconds'] = None

        # set to None if extra_params (nullable) is None
        # and model_fields_set contains the field
        if self.extra_params is None and "extra_params" in self.model_fields_set:
            _dict['extra_params'] = None

        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of TaskRequest from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "simulator":
                obj.get("simulator"),
            "simulator_name_alias":
                obj.get("simulator_name_alias"),
            "resource_pool":
                obj.get("resource_pool"),
            "storage_path_prefix":
                obj.get("storage_path_prefix")
                if obj.get("storage_path_prefix") is not None else '',
            "project":
                obj.get("project"),
            "container_image":
                obj.get("container_image"),
            "time_to_live_seconds":
                obj.get("time_to_live_seconds"),
            "resubmit_on_preemption":
                obj.get("resubmit_on_preemption")
                if obj.get("resubmit_on_preemption") is not None else False,
            "input_resources":
                obj.get("input_resources"),
            "stream_zip":
                obj.get("stream_zip")
                if obj.get("stream_zip") is not None else True,
            "compress_with":
                obj.get("compress_with"),
            "extra_params":
                obj.get("extra_params")
        })
        return _obj
