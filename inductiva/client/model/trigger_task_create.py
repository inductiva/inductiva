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


class TriggerTaskCreate(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "trigger_type",
            "task_id",
            "trigger",
        }

        class properties:

            class trigger_type(schemas.EnumBase, schemas.StrSchema):

                class MetaOapg:
                    enum_value_to_name = {
                        "task": "TASK",
                    }

                @schemas.classproperty
                def TASK(cls):
                    return cls("task")

            task_id = schemas.StrSchema

            @staticmethod
            def trigger() -> typing.Type['TriggerTaskType']:
                return TriggerTaskType

            __annotations__ = {
                "trigger_type": trigger_type,
                "task_id": task_id,
                "trigger": trigger,
            }

    trigger_type: MetaOapg.properties.trigger_type
    task_id: MetaOapg.properties.task_id
    trigger: 'TriggerTaskType'

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["trigger_type"]
    ) -> MetaOapg.properties.trigger_type:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["task_id"]
    ) -> MetaOapg.properties.task_id:
        ...

    @typing.overload
    def __getitem__(
            self,
            name: typing_extensions.Literal["trigger"]) -> 'TriggerTaskType':
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "trigger_type",
        "task_id",
        "trigger",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["trigger_type"]
    ) -> MetaOapg.properties.trigger_type:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["task_id"]
    ) -> MetaOapg.properties.task_id:
        ...

    @typing.overload
    def get_item_oapg(
            self,
            name: typing_extensions.Literal["trigger"]) -> 'TriggerTaskType':
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "trigger_type",
        "task_id",
        "trigger",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        trigger_type: typing.Union[
            MetaOapg.properties.trigger_type,
            str,
        ],
        task_id: typing.Union[
            MetaOapg.properties.task_id,
            str,
        ],
        trigger: 'TriggerTaskType',
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'TriggerTaskCreate':
        return super().__new__(
            cls,
            *_args,
            trigger_type=trigger_type,
            task_id=task_id,
            trigger=trigger,
            _configuration=_configuration,
            **kwargs,
        )


from inductiva.client.model.trigger_task_type import TriggerTaskType
