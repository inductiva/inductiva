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


class TaskStatusInfo(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "operations",
            "alias",
            "description",
            "status",
            "timestamp",
        }

        class properties:

            @staticmethod
            def status() -> typing.Type['TaskStatusCode']:
                return TaskStatusCode

            timestamp = schemas.DateTimeSchema
            description = schemas.StrSchema
            alias = schemas.StrSchema

            class operations(schemas.ListSchema):

                class MetaOapg:

                    @staticmethod
                    def items() -> typing.Type['TaskOperation']:
                        return TaskOperation

                def __new__(
                    cls,
                    _arg: typing.Union[typing.Tuple['TaskOperation'],
                                       typing.List['TaskOperation']],
                    _configuration: typing.Optional[
                        schemas.Configuration] = None,
                ) -> 'operations':
                    return super().__new__(
                        cls,
                        _arg,
                        _configuration=_configuration,
                    )

                def __getitem__(self, i: int) -> 'TaskOperation':
                    return super().__getitem__(i)

            class machine_id(
                    schemas.UUIDBase,
                    schemas.ComposedSchema,
            ):

                class MetaOapg:
                    format = 'uuid'
                    any_of_0 = schemas.StrSchema
                    any_of_1 = schemas.NoneSchema

                    @classmethod
                    @functools.lru_cache()
                    def any_of(cls):
                        # we need this here to make our import statements work
                        # we must store _composed_schemas in here so the code is only run
                        # when we invoke this method. If we kept this at the class
                        # level we would get an error because the class level
                        # code would be run when this module is imported, and these composed
                        # classes don't exist yet because their module has not finished
                        # loading
                        return [
                            cls.any_of_0,
                            cls.any_of_1,
                        ]

                def __new__(
                    cls,
                    *_args: typing.Union[
                        dict,
                        frozendict.frozendict,
                        str,
                        date,
                        datetime,
                        uuid.UUID,
                        int,
                        float,
                        decimal.Decimal,
                        bool,
                        None,
                        list,
                        tuple,
                        bytes,
                        io.FileIO,
                        io.BufferedReader,
                    ],
                    _configuration: typing.Optional[
                        schemas.Configuration] = None,
                    **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                                           frozendict.frozendict, str, date,
                                           datetime, uuid.UUID, int, float,
                                           decimal.Decimal, None, list, tuple,
                                           bytes],
                ) -> 'machine_id':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class end_timestamp(
                    schemas.DateTimeBase,
                    schemas.ComposedSchema,
            ):

                class MetaOapg:
                    format = 'date-time'
                    any_of_0 = schemas.StrSchema
                    any_of_1 = schemas.NoneSchema

                    @classmethod
                    @functools.lru_cache()
                    def any_of(cls):
                        # we need this here to make our import statements work
                        # we must store _composed_schemas in here so the code is only run
                        # when we invoke this method. If we kept this at the class
                        # level we would get an error because the class level
                        # code would be run when this module is imported, and these composed
                        # classes don't exist yet because their module has not finished
                        # loading
                        return [
                            cls.any_of_0,
                            cls.any_of_1,
                        ]

                def __new__(
                    cls,
                    *_args: typing.Union[
                        dict,
                        frozendict.frozendict,
                        str,
                        date,
                        datetime,
                        uuid.UUID,
                        int,
                        float,
                        decimal.Decimal,
                        bool,
                        None,
                        list,
                        tuple,
                        bytes,
                        io.FileIO,
                        io.BufferedReader,
                    ],
                    _configuration: typing.Optional[
                        schemas.Configuration] = None,
                    **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                                           frozendict.frozendict, str, date,
                                           datetime, uuid.UUID, int, float,
                                           decimal.Decimal, None, list, tuple,
                                           bytes],
                ) -> 'end_timestamp':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            __annotations__ = {
                "status": status,
                "timestamp": timestamp,
                "description": description,
                "alias": alias,
                "operations": operations,
                "machine_id": machine_id,
                "end_timestamp": end_timestamp,
            }

    operations: MetaOapg.properties.operations
    alias: MetaOapg.properties.alias
    description: MetaOapg.properties.description
    status: 'TaskStatusCode'
    timestamp: MetaOapg.properties.timestamp

    @typing.overload
    def __getitem__(
            self,
            name: typing_extensions.Literal["status"]) -> 'TaskStatusCode':
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["timestamp"]
    ) -> MetaOapg.properties.timestamp:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["description"]
    ) -> MetaOapg.properties.description:
        ...

    @typing.overload
    def __getitem__(
            self, name: typing_extensions.Literal["alias"]
    ) -> MetaOapg.properties.alias:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["operations"]
    ) -> MetaOapg.properties.operations:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["machine_id"]
    ) -> MetaOapg.properties.machine_id:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["end_timestamp"]
    ) -> MetaOapg.properties.end_timestamp:
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "status",
        "timestamp",
        "description",
        "alias",
        "operations",
        "machine_id",
        "end_timestamp",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
            self,
            name: typing_extensions.Literal["status"]) -> 'TaskStatusCode':
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["timestamp"]
    ) -> MetaOapg.properties.timestamp:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["description"]
    ) -> MetaOapg.properties.description:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: typing_extensions.Literal["alias"]
    ) -> MetaOapg.properties.alias:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["operations"]
    ) -> MetaOapg.properties.operations:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["machine_id"]
    ) -> typing.Union[MetaOapg.properties.machine_id, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["end_timestamp"]
    ) -> typing.Union[MetaOapg.properties.end_timestamp, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "status",
        "timestamp",
        "description",
        "alias",
        "operations",
        "machine_id",
        "end_timestamp",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        operations: typing.Union[
            MetaOapg.properties.operations,
            list,
            tuple,
        ],
        alias: typing.Union[
            MetaOapg.properties.alias,
            str,
        ],
        description: typing.Union[
            MetaOapg.properties.description,
            str,
        ],
        status: 'TaskStatusCode',
        timestamp: typing.Union[
            MetaOapg.properties.timestamp,
            str,
            datetime,
        ],
        machine_id: typing.Union[MetaOapg.properties.machine_id, dict,
                                 frozendict.frozendict, str, date, datetime,
                                 uuid.UUID, int, float, decimal.Decimal, bool,
                                 None, list, tuple, bytes, io.FileIO,
                                 io.BufferedReader,
                                 schemas.Unset] = schemas.unset,
        end_timestamp: typing.Union[MetaOapg.properties.end_timestamp, dict,
                                    frozendict.frozendict, str, date, datetime,
                                    uuid.UUID, int, float, decimal.Decimal,
                                    bool, None, list, tuple, bytes, io.FileIO,
                                    io.BufferedReader,
                                    schemas.Unset] = schemas.unset,
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'TaskStatusInfo':
        return super().__new__(
            cls,
            *_args,
            operations=operations,
            alias=alias,
            description=description,
            status=status,
            timestamp=timestamp,
            machine_id=machine_id,
            end_timestamp=end_timestamp,
            _configuration=_configuration,
            **kwargs,
        )


from inductiva.client.model.task_operation import TaskOperation
from inductiva.client.model.task_status_code import TaskStatusCode