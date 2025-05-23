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


class MachineGroupCostPerHour(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "min",
            "max",
            "max_reason",
            "min_reason",
        }

        class properties:
            min = schemas.NumberSchema
            max = schemas.NumberSchema
            min_reason = schemas.StrSchema
            max_reason = schemas.StrSchema

            class current(
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
                ) -> 'current':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            class currency(
                    schemas.ComposedSchema,):

                class MetaOapg:

                    @classmethod
                    @functools.lru_cache()
                    def all_of(cls):
                        # we need this here to make our import statements work
                        # we must store _composed_schemas in here so the code is only run
                        # when we invoke this method. If we kept this at the class
                        # level we would get an error because the class level
                        # code would be run when this module is imported, and these composed
                        # classes don't exist yet because their module has not finished
                        # loading
                        return [
                            CurrencyCode,
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
                ) -> 'currency':
                    return super().__new__(
                        cls,
                        *_args,
                        _configuration=_configuration,
                        **kwargs,
                    )

            __annotations__ = {
                "min": min,
                "max": max,
                "min_reason": min_reason,
                "max_reason": max_reason,
                "current": current,
                "currency": currency,
            }

    min: MetaOapg.properties.min
    max: MetaOapg.properties.max
    max_reason: MetaOapg.properties.max_reason
    min_reason: MetaOapg.properties.min_reason

    @typing.overload
    def __getitem__(
            self,
            name: typing_extensions.Literal["min"]) -> MetaOapg.properties.min:
        ...

    @typing.overload
    def __getitem__(
            self,
            name: typing_extensions.Literal["max"]) -> MetaOapg.properties.max:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["min_reason"]
    ) -> MetaOapg.properties.min_reason:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["max_reason"]
    ) -> MetaOapg.properties.max_reason:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["current"]
    ) -> MetaOapg.properties.current:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["currency"]
    ) -> MetaOapg.properties.currency:
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "min",
        "max",
        "min_reason",
        "max_reason",
        "current",
        "currency",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
            self,
            name: typing_extensions.Literal["min"]) -> MetaOapg.properties.min:
        ...

    @typing.overload
    def get_item_oapg(
            self,
            name: typing_extensions.Literal["max"]) -> MetaOapg.properties.max:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["min_reason"]
    ) -> MetaOapg.properties.min_reason:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["max_reason"]
    ) -> MetaOapg.properties.max_reason:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["current"]
    ) -> typing.Union[MetaOapg.properties.current, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["currency"]
    ) -> typing.Union[MetaOapg.properties.currency, schemas.Unset]:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "min",
        "max",
        "min_reason",
        "max_reason",
        "current",
        "currency",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        min: typing.Union[
            MetaOapg.properties.min,
            decimal.Decimal,
            int,
            float,
        ],
        max: typing.Union[
            MetaOapg.properties.max,
            decimal.Decimal,
            int,
            float,
        ],
        max_reason: typing.Union[
            MetaOapg.properties.max_reason,
            str,
        ],
        min_reason: typing.Union[
            MetaOapg.properties.min_reason,
            str,
        ],
        current: typing.Union[MetaOapg.properties.current, dict,
                              frozendict.frozendict, str, date, datetime,
                              uuid.UUID, int, float, decimal.Decimal, bool,
                              None, list, tuple, bytes, io.FileIO,
                              io.BufferedReader, schemas.Unset] = schemas.unset,
        currency: typing.Union[MetaOapg.properties.currency, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, bool,
                               None, list, tuple, bytes, io.FileIO,
                               io.BufferedReader,
                               schemas.Unset] = schemas.unset,
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'MachineGroupCostPerHour':
        return super().__new__(
            cls,
            *_args,
            min=min,
            max=max,
            max_reason=max_reason,
            min_reason=min_reason,
            current=current,
            currency=currency,
            _configuration=_configuration,
            **kwargs,
        )


from inductiva.client.model.currency_code import CurrencyCode
