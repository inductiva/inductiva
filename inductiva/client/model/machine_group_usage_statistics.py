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


class MachineGroupUsageStatistics(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "idle_time",
            "total_time",
            "task_time",
        }

        class properties:

            @staticmethod
            def total_time() -> typing.Type['MachineGroupTimeStat']:
                return MachineGroupTimeStat

            @staticmethod
            def task_time() -> typing.Type['MachineGroupTimeStat']:
                return MachineGroupTimeStat

            @staticmethod
            def idle_time() -> typing.Type['MachineGroupTimeStat']:
                return MachineGroupTimeStat

            __annotations__ = {
                "total_time": total_time,
                "task_time": task_time,
                "idle_time": idle_time,
            }

    idle_time: 'MachineGroupTimeStat'
    total_time: 'MachineGroupTimeStat'
    task_time: 'MachineGroupTimeStat'

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def __getitem__(
            self, name: typing_extensions.Literal["task_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def __getitem__(
            self, name: typing_extensions.Literal["idle_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "total_time",
        "task_time",
        "idle_time",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def get_item_oapg(
            self, name: typing_extensions.Literal["task_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def get_item_oapg(
            self, name: typing_extensions.Literal["idle_time"]
    ) -> 'MachineGroupTimeStat':
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "total_time",
        "task_time",
        "idle_time",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        idle_time: 'MachineGroupTimeStat',
        total_time: 'MachineGroupTimeStat',
        task_time: 'MachineGroupTimeStat',
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'MachineGroupUsageStatistics':
        return super().__new__(
            cls,
            *_args,
            idle_time=idle_time,
            total_time=total_time,
            task_time=task_time,
            _configuration=_configuration,
            **kwargs,
        )


from inductiva.client.model.machine_group_time_stat import MachineGroupTimeStat