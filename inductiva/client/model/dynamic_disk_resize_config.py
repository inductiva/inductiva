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


class DynamicDiskResizeConfig(schemas.DictSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """

    class MetaOapg:
        required = {
            "free_space_threshold_gb",
            "size_increment_gb",
            "max_disk_size_gb",
        }

        class properties:
            free_space_threshold_gb = schemas.IntSchema
            size_increment_gb = schemas.IntSchema
            max_disk_size_gb = schemas.IntSchema
            __annotations__ = {
                "free_space_threshold_gb": free_space_threshold_gb,
                "size_increment_gb": size_increment_gb,
                "max_disk_size_gb": max_disk_size_gb,
            }

    free_space_threshold_gb: MetaOapg.properties.free_space_threshold_gb
    size_increment_gb: MetaOapg.properties.size_increment_gb
    max_disk_size_gb: MetaOapg.properties.max_disk_size_gb

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["free_space_threshold_gb"]
    ) -> MetaOapg.properties.free_space_threshold_gb:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["size_increment_gb"]
    ) -> MetaOapg.properties.size_increment_gb:
        ...

    @typing.overload
    def __getitem__(
        self, name: typing_extensions.Literal["max_disk_size_gb"]
    ) -> MetaOapg.properties.max_disk_size_gb:
        ...

    @typing.overload
    def __getitem__(self, name: str) -> schemas.UnsetAnyTypeSchema:
        ...

    def __getitem__(self, name: typing.Union[typing_extensions.Literal[
        "free_space_threshold_gb",
        "size_increment_gb",
        "max_disk_size_gb",
    ], str]):
        # dict_instance[name] accessor
        return super().__getitem__(name)

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["free_space_threshold_gb"]
    ) -> MetaOapg.properties.free_space_threshold_gb:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["size_increment_gb"]
    ) -> MetaOapg.properties.size_increment_gb:
        ...

    @typing.overload
    def get_item_oapg(
        self, name: typing_extensions.Literal["max_disk_size_gb"]
    ) -> MetaOapg.properties.max_disk_size_gb:
        ...

    @typing.overload
    def get_item_oapg(
            self, name: str
    ) -> typing.Union[schemas.UnsetAnyTypeSchema, schemas.Unset]:
        ...

    def get_item_oapg(self, name: typing.Union[typing_extensions.Literal[
        "free_space_threshold_gb",
        "size_increment_gb",
        "max_disk_size_gb",
    ], str]):
        return super().get_item_oapg(name)

    def __new__(
        cls,
        *_args: typing.Union[
            dict,
            frozendict.frozendict,
        ],
        free_space_threshold_gb: typing.Union[
            MetaOapg.properties.free_space_threshold_gb,
            decimal.Decimal,
            int,
        ],
        size_increment_gb: typing.Union[
            MetaOapg.properties.size_increment_gb,
            decimal.Decimal,
            int,
        ],
        max_disk_size_gb: typing.Union[
            MetaOapg.properties.max_disk_size_gb,
            decimal.Decimal,
            int,
        ],
        _configuration: typing.Optional[schemas.Configuration] = None,
        **kwargs: typing.Union[schemas.AnyTypeSchema, dict,
                               frozendict.frozendict, str, date, datetime,
                               uuid.UUID, int, float, decimal.Decimal, None,
                               list, tuple, bytes],
    ) -> 'DynamicDiskResizeConfig':
        return super().__new__(
            cls,
            *_args,
            free_space_threshold_gb=free_space_threshold_gb,
            size_increment_gb=size_increment_gb,
            max_disk_size_gb=max_disk_size_gb,
            _configuration=_configuration,
            **kwargs,
        )