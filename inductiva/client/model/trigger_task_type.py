# coding: utf-8
"""
    InductivaWebAPI

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)  # noqa: E501

    The version of the OpenAPI document: 0.1.0
    Generated by: https://openapi-generator.tech
"""

from datetime import date, datetime  # noqa: F401
import decimal  # noqa: F401
import functools  # noqa: F401
import io  # noqa: F401
import re  # noqa: F401
import typing  # noqa: F401
import typing_extensions  # noqa: F401
import uuid  # noqa: F401

import frozendict  # noqa: F401

from inductiva.client import schemas  # noqa: F401


class TriggerTaskType(schemas.EnumBase, schemas.StrSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.

    Task triggers.
    """

    class MetaOapg:
        enum_value_to_name = {
            "task_output_uploaded": "OUTPUT_UPLOADED",
            "task_failed": "FAILED",
            "task_started": "STARTED",
        }

    @schemas.classproperty
    def OUTPUT_UPLOADED(cls):
        return cls("task_output_uploaded")

    @schemas.classproperty
    def FAILED(cls):
        return cls("task_failed")

    @schemas.classproperty
    def STARTED(cls):
        return cls("task_started")
