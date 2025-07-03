"""Tests for the projects class"""
import inductiva
from inductiva.client.models.event_create import EventCreate
from inductiva.events import triggers, actions
import pytest


def test_create_event():
    """Tests if can create new project."""

    task = inductiva.tasks.Task("123")

    trigger = triggers.TaskOutputUploaded(task_id=task.id)
    action = actions.WebhookNotification(webhook_url="https://example.com")
    event = EventCreate.from_dict({
        "trigger": trigger.get_trigger(),
        "action": action.get_action()
    })

    assert event is not None

    assert trigger.model_dump().get("task_id") == trigger.get_trigger().get(
        "task_id")
    assert action.model_dump().get("webhook_url") == action.get_action().get(
        "webhook_url")

    action = actions.EmailNotification(email_address="email@email.com")
    event = EventCreate.from_dict({
        "trigger": trigger.get_trigger(),
        "action": action.get_action()
    })

    assert action.model_dump().get("email_address") == action.get_action().get(
        "email_address")
