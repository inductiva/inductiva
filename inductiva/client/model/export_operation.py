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


class ExportOperation(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.

    Body of the request to the export files endpoint.
    """

    class MetaOapg:
        required = {
            "path",
            "dest_url",
        }

        class properties:
            path = schemas.StrSchema

            class dest_url(schemas.StrSchema):

                class MetaOapg:
                    format = 'uri'
                    max_length = 65536
                    min_length = 1

            __annotations__ = {
                "path": path,
                "dest_url": dest_url,
            }

    path: MetaOapg.properties.path
    dest_url: MetaOapg.properties.dest_url

    @typing.overload
    def __getitem__(
            self, name: typing_extensions.Literal["path"]
    ) -> MetaOapg.properties.path:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["dest_url"]
    ) -> MetaOapg.properties.dest_url:
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "path",
        "dest_url",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
            self, name: typing_extensions.Literal["path"]
    ) -> MetaOapg.properties.path:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["dest_url"]
    ) -> MetaOapg.properties.dest_url:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "path",
        "dest_url",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        path: typing.Union[
            MetaOapg.properties.path,
            str,
        ],
        dest_url: typing.Union[
            MetaOapg.properties.dest_url,
            str,
        ],
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'ExportOperation':
        return super().__new__(
            cls,
            *_args,
            path=path,
            dest_url=dest_url,
            _configuration=_configuration,
            **kwargs,
        )