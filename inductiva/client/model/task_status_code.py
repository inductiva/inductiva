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


class TaskStatusCode(schemas.EnumBase, schemas.StrSchema):
    """NOTE: This class is auto generated by OpenAPI Generator.
    Ref: https://openapi-generator.tech

    Do not edit the class manually.

    Possible task status codes.
    """

    class MetaOapg:
        enum_value_to_name = {
            "pending-input":
                "PENDINGINPUT",
            "submitted":
                "SUBMITTED",
            "started":
                "STARTED",
            "success":
                "SUCCESS",
            "failed":
                "FAILED",
            "pending-kill":
                "PENDINGKILL",
            "killed":
                "KILLED",
            "spot-instance-preempted":
                "SPOTINSTANCEPREEMPTED",
            "executer-terminated":
                "EXECUTERTERMINATED",
            "executer-terminated-by-user":
                "EXECUTERTERMINATEDBYUSER",
            "executer-terminated-ttl-exceeded":
                "EXECUTERTERMINATEDTTLEXCEEDED",
            "executer-terminated-credits-exhausted":
                "EXECUTERTERMINATEDCREDITSEXHAUSTED",
            "executer-failed":
                "EXECUTERFAILED",
            "zombie":
                "ZOMBIE",
            "computation-started":
                "COMPUTATIONSTARTED",
            "computation-ended":
                "COMPUTATIONENDED",
            "ttl-exceeded":
                "TTLEXCEEDED",
        }

    @schemas.classproperty
    def PENDINGINPUT(cls):
        return cls("pending-input")

    @schemas.classproperty
    def SUBMITTED(cls):
        return cls("submitted")

    @schemas.classproperty
    def STARTED(cls):
        return cls("started")

    @schemas.classproperty
    def SUCCESS(cls):
        return cls("success")

    @schemas.classproperty
    def FAILED(cls):
        return cls("failed")

    @schemas.classproperty
    def PENDINGKILL(cls):
        return cls("pending-kill")

    @schemas.classproperty
    def KILLED(cls):
        return cls("killed")

    @schemas.classproperty
    def SPOTINSTANCEPREEMPTED(cls):
        return cls("spot-instance-preempted")

    @schemas.classproperty
    def EXECUTERTERMINATED(cls):
        return cls("executer-terminated")

    @schemas.classproperty
    def EXECUTERTERMINATEDBYUSER(cls):
        return cls("executer-terminated-by-user")

    @schemas.classproperty
    def EXECUTERTERMINATEDTTLEXCEEDED(cls):
        return cls("executer-terminated-ttl-exceeded")

    @schemas.classproperty
    def EXECUTERTERMINATEDCREDITSEXHAUSTED(cls):
        return cls("executer-terminated-credits-exhausted")

    @schemas.classproperty
    def EXECUTERFAILED(cls):
        return cls("executer-failed")

    @schemas.classproperty
    def ZOMBIE(cls):
        return cls("zombie")

    @schemas.classproperty
    def COMPUTATIONSTARTED(cls):
        return cls("computation-started")

    @schemas.classproperty
    def COMPUTATIONENDED(cls):
        return cls("computation-ended")

    @schemas.classproperty
    def TTLEXCEEDED(cls):
        return cls("ttl-exceeded")
