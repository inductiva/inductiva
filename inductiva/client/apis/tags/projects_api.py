# coding: utf-8
"""
    InductivaWebAPI

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)  # noqa: E501

    The version of the OpenAPI document: 0.1.0
    Generated by: https://openapi-generator.tech
"""

from inductiva.client.paths.projects.post import CreateProject
from inductiva.client.paths.projects_name.get import GetProject
from inductiva.client.paths.projects.get import GetUserProjects


class ProjectsApi(
        CreateProject,
        GetProject,
        GetUserProjects,
):
    """NOTE: This class is auto generated by OpenAPI Generator
    Ref: https://openapi-generator.tech

    Do not edit the class manually.
    """
    pass