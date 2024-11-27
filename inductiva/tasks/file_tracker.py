"""File Tracker module for connecting to a running task via WebRTC."""
import asyncio
import json
import uuid
import aiohttp
import enum
import os
from aiortc import RTCPeerConnection, RTCSessionDescription

from inductiva.client.apis.tags import tasks_api

SIGNALING_SERVER = os.environ.get("INDUCTIVA_API_URL",
                                  "https://api.inductiva.ai")
API_KEY = os.environ.get("INDUCTIVA_API_KEY", None)

# STUN/TURN server configuration
ICE_SERVERS = [{
    "urls": ["stun:34.79.246.4:3478"]
}, {
    "urls": ["turn:34.79.246.4:3478"]
}]


class Operations(enum.Enum):
    LIST = "ls"
    TAIL = "tail"


class FileTracker:
    """File Tracker class for connecting to a running task via WebRTC."""

    def __init__(self, api):
        self.pc = RTCPeerConnection()
        self.pc.configuration = {"iceServers": ICE_SERVERS}
        self._message = None
        self._api = tasks_api.TasksApi(api)

    async def setup_channel(self, operation, **kwargs):
        channel = self.pc.createDataChannel("file_transfer")
        fut = asyncio.Future()

        @channel.on("open")
        def on_open():
            request = operation.value
            if kwargs:
                request += ":" + json.dumps(",".join(
                    kwargs.values())).strip("\"")
            channel.send(request)

        @channel.on("message")
        def on_message(message):
            self._message = json.loads(message)
            fut.set_result(self._message)

        return fut

    async def connect_to_task(self, task_id):
        connection_id = str(uuid.uuid4())
        path_params={"task_id": task_id}
        self._api.register_task(body={"sender_id": connection_id}, path_params=path_params)

        offer = await self.pc.createOffer()
        await self.pc.setLocalDescription(offer)

        self._api.offer_task(body={
                                "sender_id": connection_id,
                                "receiver_id": task_id,
                                "type": "offer",
                                "sdp": self.pc.localDescription.sdp
                                }, path_params=path_params)

        resp = self._api.get_message(query_params={"client": connection_id}, path_params=path_params)
        data = resp.body
        if data["type"] == "answer":
            await self.pc.setRemoteDescription(
                RTCSessionDescription(sdp=data["sdp"],
                                        type=data["type"]))

    async def cleanup(self):
        await self.pc.close()
