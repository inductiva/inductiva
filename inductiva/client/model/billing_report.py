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


class BillingReport(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "start_ts",
            "total_cost",
            "end_ts",
            "report",
        }

        class properties:
            start_ts = schemas.DateTimeSchema
            end_ts = schemas.DateTimeSchema
            total_cost = schemas.NumberSchema

            class report(schemas.ListSchema):

                class MetaOapg:

                    @staticmethod
                    def items() -> typing.Type['BillingReportEntity']:
                        return BillingReportEntity

                def __new__(
                    cls,
                    _arg: typing.Union[typing.Tuple['BillingReportEntity'],
                                       typing.List['BillingReportEntity']],
                    _configuration: typing.Optional[
                        schemas.Configuration] = None,
                ) -> 'report':
                    return super().__new__(
                        cls,
                        _arg,
                        _configuration=_configuration,
                    )

                def __getitem__(self, i: int) -> 'BillingReportEntity':
                    return super().__getitem__(i)

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
                "start_ts": start_ts,
                "end_ts": end_ts,
                "total_cost": total_cost,
                "report": report,
                "currency": currency,
            }

    start_ts: MetaOapg.properties.start_ts
    total_cost: MetaOapg.properties.total_cost
    end_ts: MetaOapg.properties.end_ts
    report: MetaOapg.properties.report

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["start_ts"]
    ) -> MetaOapg.properties.start_ts:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["end_ts"]
    ) -> MetaOapg.properties.end_ts:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["total_cost"]
    ) -> MetaOapg.properties.total_cost:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["report"]
    ) -> MetaOapg.properties.report:
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
        "start_ts",
        "end_ts",
        "total_cost",
        "report",
        "currency",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["start_ts"]
    ) -> MetaOapg.properties.start_ts:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["end_ts"]
    ) -> MetaOapg.properties.end_ts:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["total_cost"]
    ) -> MetaOapg.properties.total_cost:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["report"]
    ) -> MetaOapg.properties.report:
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
        "start_ts",
        "end_ts",
        "total_cost",
        "report",
        "currency",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        start_ts: typing.Union[
            MetaOapg.properties.start_ts,
            str,
            datetime,
        ],
        total_cost: typing.Union[
            MetaOapg.properties.total_cost,
            decimal.Decimal,
            int,
            float,
        ],
        end_ts: typing.Union[
            MetaOapg.properties.end_ts,
            str,
            datetime,
        ],
        report: typing.Union[
            MetaOapg.properties.report,
            list,
            tuple,
        ],
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
    ) -> 'BillingReport':
        return super().__new__(
            cls,
            *_args,
            start_ts=start_ts,
            total_cost=total_cost,
            end_ts=end_ts,
            report=report,
            currency=currency,
            _configuration=_configuration,
            **kwargs,
        )


from inductiva.client.model.billing_report_entity import BillingReportEntity
from inductiva.client.model.currency_code import CurrencyCode