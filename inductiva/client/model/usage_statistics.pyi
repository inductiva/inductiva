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


class UsageStatistics(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:

        class properties:

            class period(
                    schemas.ComposedSchema,):

                class MetaOapg:
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
                ) -> 'period':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_core_hours_failed(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'total_core_hours_failed':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_core_hours(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'total_core_hours':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_tasks(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.IntSchema
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
                ) -> 'total_tasks':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_failed_tasks(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.IntSchema
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
                ) -> 'total_failed_tasks':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class avg_total_tasks_duration(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'avg_total_tasks_duration':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class avg_computation_seconds_task_duration(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'avg_computation_seconds_task_duration':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class success_rate(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'success_rate':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class avg_estimated_cost(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'avg_estimated_cost':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_estimated_cost(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'total_estimated_cost':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class total_failed_cost(
                    schemas.ComposedSchema,):

                class MetaOapg:
                    any_of_0 = schemas.NumberSchema
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
                ) -> 'total_failed_cost':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            __annotations__ = {
                "period":
                    period,
                "total_core_hours_failed":
                    total_core_hours_failed,
                "total_core_hours":
                    total_core_hours,
                "total_tasks":
                    total_tasks,
                "total_failed_tasks":
                    total_failed_tasks,
                "avg_total_tasks_duration":
                    avg_total_tasks_duration,
                "avg_computation_seconds_task_duration":
                    avg_computation_seconds_task_duration,
                "success_rate":
                    success_rate,
                "avg_estimated_cost":
                    avg_estimated_cost,
                "total_estimated_cost":
                    total_estimated_cost,
                "total_failed_cost":
                    total_failed_cost,
            }

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["period"]
    ) -> MetaOapg.properties.period:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_core_hours_failed"]
    ) -> MetaOapg.properties.total_core_hours_failed:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_core_hours"]
    ) -> MetaOapg.properties.total_core_hours:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_tasks"]
    ) -> MetaOapg.properties.total_tasks:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_failed_tasks"]
    ) -> MetaOapg.properties.total_failed_tasks:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["avg_total_tasks_duration"]
    ) -> MetaOapg.properties.avg_total_tasks_duration:
        ...

    @typing.overload
    def __getitem__(
        self,
        name: typing_extensions.Literal["avg_computation_seconds_task_duration"]
    ) -> MetaOapg.properties.avg_computation_seconds_task_duration:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["success_rate"]
    ) -> MetaOapg.properties.success_rate:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["avg_estimated_cost"]
    ) -> MetaOapg.properties.avg_estimated_cost:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_estimated_cost"]
    ) -> MetaOapg.properties.total_estimated_cost:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_failed_cost"]
    ) -> MetaOapg.properties.total_failed_cost:
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "period",
        "total_core_hours_failed",
        "total_core_hours",
        "total_tasks",
        "total_failed_tasks",
        "avg_total_tasks_duration",
        "avg_computation_seconds_task_duration",
        "success_rate",
        "avg_estimated_cost",
        "total_estimated_cost",
        "total_failed_cost",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["period"]
    ) -> typing.Union[MetaOapg.properties.period, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_core_hours_failed"]
    ) -> typing.Union[MetaOapg.properties.total_core_hours_failed,
                      schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_core_hours"]
    ) -> typing.Union[MetaOapg.properties.total_core_hours, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_tasks"]
    ) -> typing.Union[MetaOapg.properties.total_tasks, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_failed_tasks"]
    ) -> typing.Union[MetaOapg.properties.total_failed_tasks, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["avg_total_tasks_duration"]
    ) -> typing.Union[MetaOapg.properties.avg_total_tasks_duration,
                      schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self,
        name: typing_extensions.Literal["avg_computation_seconds_task_duration"]
    ) -> typing.Union[MetaOapg.properties.avg_computation_seconds_task_duration,
                      schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["success_rate"]
    ) -> typing.Union[MetaOapg.properties.success_rate, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["avg_estimated_cost"]
    ) -> typing.Union[MetaOapg.properties.avg_estimated_cost, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_estimated_cost"]
    ) -> typing.Union[MetaOapg.properties.total_estimated_cost, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_failed_cost"]
    ) -> typing.Union[MetaOapg.properties.total_failed_cost, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "period",
        "total_core_hours_failed",
        "total_core_hours",
        "total_tasks",
        "total_failed_tasks",
        "avg_total_tasks_duration",
        "avg_computation_seconds_task_duration",
        "success_rate",
        "avg_estimated_cost",
        "total_estimated_cost",
        "total_failed_cost",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        period: typing.Union[MetaOapg.properties.period, dict,
                             frozendict.frozendict, str, date, datetime,
                             uuid.UUID, int, float, decimal.Decimal, bool, None,
                             list, tuple, bytes, io.FileIO, io.BufferedReader,
                             schemas.Unset] = schemas.unset,
        total_core_hours_failed: typing.Union[
            MetaOapg.properties.total_core_hours_failed, dict,
            frozendict.frozendict, str, date, datetime, uuid.UUID, int, float,
            decimal.Decimal, bool, None, list, tuple, bytes, io.FileIO,
            io.BufferedReader, schemas.Unset] = schemas.unset,
        total_core_hours: typing.Union[MetaOapg.properties.total_core_hours,
                                       dict, frozendict.frozendict, str, date,
                                       datetime, uuid.UUID, int, float,
                                       decimal.Decimal, bool, None, list, tuple,
                                       bytes, io.FileIO, io.BufferedReader,
                                       schemas.Unset] = schemas.unset,
        total_tasks: typing.Union[MetaOapg.properties.total_tasks, dict,
                                  frozendict.frozendict, str, date, datetime,
                                  uuid.UUID, int, float, decimal.Decimal, bool,
                                  None, list, tuple, bytes, io.FileIO,
                                  io.BufferedReader,
                                  schemas.Unset] = schemas.unset,
        total_failed_tasks: typing.Union[MetaOapg.properties.total_failed_tasks,
                                         dict, frozendict.frozendict, str, date,
                                         datetime, uuid.UUID, int, float,
                                         decimal.Decimal, bool, None, list,
                                         tuple, bytes, io.FileIO,
                                         io.BufferedReader,
                                         schemas.Unset] = schemas.unset,
        avg_total_tasks_duration: typing.Union[
            MetaOapg.properties.avg_total_tasks_duration, dict,
            frozendict.frozendict, str, date, datetime, uuid.UUID, int, float,
            decimal.Decimal, bool, None, list, tuple, bytes, io.FileIO,
            io.BufferedReader, schemas.Unset] = schemas.unset,
        avg_computation_seconds_task_duration: typing.Union[
            MetaOapg.properties.avg_computation_seconds_task_duration, dict,
            frozendict.frozendict, str, date, datetime, uuid.UUID, int, float,
            decimal.Decimal, bool, None, list, tuple, bytes, io.FileIO,
            io.BufferedReader, schemas.Unset] = schemas.unset,
        success_rate: typing.Union[MetaOapg.properties.success_rate, dict,
                                   frozendict.frozendict, str, date, datetime,
                                   uuid.UUID, int, float, decimal.Decimal, bool,
                                   None, list, tuple, bytes, io.FileIO,
                                   io.BufferedReader,
                                   schemas.Unset] = schemas.unset,
        avg_estimated_cost: typing.Union[MetaOapg.properties.avg_estimated_cost,
                                         dict, frozendict.frozendict, str, date,
                                         datetime, uuid.UUID, int, float,
                                         decimal.Decimal, bool, None, list,
                                         tuple, bytes, io.FileIO,
                                         io.BufferedReader,
                                         schemas.Unset] = schemas.unset,
        total_estimated_cost: typing.Union[
            MetaOapg.properties.total_estimated_cost, dict,
            frozendict.frozendict, str, date, datetime, uuid.UUID, int, float,
            decimal.Decimal, bool, None, list, tuple, bytes, io.FileIO,
            io.BufferedReader, schemas.Unset] = schemas.unset,
        total_failed_cost: typing.Union[MetaOapg.properties.total_failed_cost,
                                        dict, frozendict.frozendict, str, date,
                                        datetime, uuid.UUID, int, float,
                                        decimal.Decimal, bool, None, list,
                                        tuple, bytes, io.FileIO,
                                        io.BufferedReader,
                                        schemas.Unset] = schemas.unset,
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'UsageStatistics':
        return super().__new__(
            cls,
            *_args,
            period=period,
            total_core_hours_failed=total_core_hours_failed,
            total_core_hours=total_core_hours,
            total_tasks=total_tasks,
            total_failed_tasks=total_failed_tasks,
            avg_total_tasks_duration=avg_total_tasks_duration,
            avg_computation_seconds_task_duration=
            avg_computation_seconds_task_duration,
            success_rate=success_rate,
            avg_estimated_cost=avg_estimated_cost,
            total_estimated_cost=total_estimated_cost,
            total_failed_cost=total_failed_cost,
            _configuration=_configuration,
            **kwargs,
        )
