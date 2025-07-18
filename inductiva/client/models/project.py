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

from datetime import datetime
from pydantic import BaseModel, ConfigDict, Field, StrictFloat, StrictInt, StrictStr
from typing import Any, ClassVar, Dict, List, Optional, Union
from typing_extensions import Annotated
from inductiva.client.models.project_statistics import ProjectStatistics
from inductiva.client.models.project_type import ProjectType
from typing import Optional, Set
from typing_extensions import Self


class Project(BaseModel):
    """
    Project
    """

  # noqa: E501
    name: Annotated[str, Field(min_length=1, strict=True, max_length=128)]
    project_type: ProjectType
    id: StrictStr
    created_at: datetime
    num_tasks: StrictInt
    estimated_computation_cost: Optional[Union[StrictFloat, StrictInt]] = 0.0
    task_status_overview: Optional[Dict[str, Optional[StrictInt]]] = None
    project_metadata: Optional[Dict[str, StrictStr]] = None
    statistics: Optional[ProjectStatistics] = None
    __properties: ClassVar[List[str]] = [
        "name", "project_type", "id", "created_at", "num_tasks",
        "estimated_computation_cost", "task_status_overview",
        "project_metadata", "statistics"
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
        """Create an instance of Project from a JSON string"""
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
        # override the default output from pydantic by calling `to_dict()` of statistics
        if self.statistics:
            _dict['statistics'] = self.statistics.to_dict()
        # set to None if project_metadata (nullable) is None
        # and model_fields_set contains the field
        if self.project_metadata is None and "project_metadata" in self.model_fields_set:
            _dict['project_metadata'] = None

        # set to None if statistics (nullable) is None
        # and model_fields_set contains the field
        if self.statistics is None and "statistics" in self.model_fields_set:
            _dict['statistics'] = None

        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of Project from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "name":
                obj.get("name"),
            "project_type":
                obj.get("project_type"),
            "id":
                obj.get("id"),
            "created_at":
                obj.get("created_at"),
            "num_tasks":
                obj.get("num_tasks"),
            "estimated_computation_cost":
                obj.get("estimated_computation_cost")
                if obj.get("estimated_computation_cost") is not None else 0.0,
            "task_status_overview":
                obj.get("task_status_overview"),
            "project_metadata":
                obj.get("project_metadata"),
            "statistics":
                ProjectStatistics.from_dict(obj["statistics"])
                if obj.get("statistics") is not None else None
        })
        return _obj
