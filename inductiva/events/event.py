"""Event module"""
import json
import logging

import inductiva.client
from inductiva import api as inductiva_api
from inductiva.client.models.event_create import EventCreate
from inductiva.client import ApiException
from inductiva.events.triggers import Trigger
from inductiva.events.actions import Action
from inductiva.events.parser import parse_event_info

_logger = logging.getLogger(__name__)


def list_events():
    """Gets all the user's events."""
    try:
        _logger.debug("Trying to get events")
        api = inductiva.client.EventsApi(inductiva_api.get_client())
        resp = api.get_events_without_preload_content()
        body = json.loads(resp.data)
    except ApiException as ex:
        _logger.error("Failed to get events", exc_info=ex)
        raise ex

    for event in body:
        print(parse_event_info(event))


def remove(event_id: str):
    """Removes an event with the given ID."""
    try:
        _logger.debug("Trying to remove event %s", event_id)
        api = inductiva.client.EventsApi(inductiva_api.get_client())
        api.delete_event(event_id=event_id)
    except ApiException as ex:
        _logger.error("Failed to remove event %s", event_id, exc_info=ex)
        raise ex
    print(f"Event {event_id} removed.")


def get(event_id: str):
    """Gets the event with the given ID."""
    try:
        _logger.debug("Trying to get event %s", event_id)
        api = inductiva.client.EventsApi(inductiva_api.get_client())
        response = api.get_event_without_preload_content(event_id=event_id)
        body = json.loads(response.data)
    except ApiException as ex:
        _logger.error("Failed to get event %s", event_id, exc_info=ex)
        raise ex
    print(parse_event_info(json.loads(body)))


def register(trigger: Trigger, action: Action):
    """Registers a trigger and action to an event."""
    try:
        _logger.debug("Trying to register trigger %s and action %s", trigger,
                      action)
        api = inductiva.client.EventsApi(inductiva_api.get_client())
        event = EventCreate.from_dict({
            "trigger": trigger.get_trigger(),
            "action": action.get_action()
        })
        response = api.create_event_without_preload_content(event)
    except ApiException as ex:
        _logger.error("Failed to register trigger %s and action %s",
                      trigger,
                      action,
                      exc_info=ex)
        raise ex
    return json.loads(response.data)
