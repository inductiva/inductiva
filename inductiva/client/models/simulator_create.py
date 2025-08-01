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

from pydantic import BaseModel, ConfigDict, StrictStr
from typing import Any, ClassVar, Dict, List, Optional
from inductiva.client.models.processor_type import ProcessorType
from inductiva.client.models.simulator_version import SimulatorVersion
from typing import Optional, Set
from typing_extensions import Self


class SimulatorCreate(BaseModel):
    """
    Schema for creating a new simulator.
    """

  # noqa: E501
    name: StrictStr
    display_name: Optional[StrictStr] = None
    recommended_machine_type: Optional[List[StrictStr]] = None
    processor_type: ProcessorType
    description: Optional[StrictStr] = None
    versions: Optional[List[SimulatorVersion]] = None
    __properties: ClassVar[List[str]] = [
        "name", "display_name", "recommended_machine_type", "processor_type",
        "description", "versions"
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
        """Create an instance of SimulatorCreate from a JSON string"""
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
        # override the default output from pydantic by calling `to_dict()` of each item in versions (list)
        _items = []
        if self.versions:
            for _item_versions in self.versions:
                if _item_versions:
                    _items.append(_item_versions.to_dict())
            _dict['versions'] = _items
        # set to None if display_name (nullable) is None
        # and model_fields_set contains the field
        if self.display_name is None and "display_name" in self.model_fields_set:
            _dict['display_name'] = None

        # set to None if recommended_machine_type (nullable) is None
        # and model_fields_set contains the field
        if self.recommended_machine_type is None and "recommended_machine_type" in self.model_fields_set:
            _dict['recommended_machine_type'] = None

        # set to None if description (nullable) is None
        # and model_fields_set contains the field
        if self.description is None and "description" in self.model_fields_set:
            _dict['description'] = None

        # set to None if versions (nullable) is None
        # and model_fields_set contains the field
        if self.versions is None and "versions" in self.model_fields_set:
            _dict['versions'] = None

        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of SimulatorCreate from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "name":
                obj.get("name"),
            "display_name":
                obj.get("display_name"),
            "recommended_machine_type":
                obj.get("recommended_machine_type"),
            "processor_type":
                obj.get("processor_type"),
            "description":
                obj.get("description"),
            "versions": [
                SimulatorVersion.from_dict(_item) for _item in obj["versions"]
            ] if obj.get("versions") is not None else None
        })
        return _obj
