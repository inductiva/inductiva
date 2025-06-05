"""File Tracker module for connecting to a running task via WebRTC."""
import asyncio
import json
import uuid
import enum
import logging
import warnings
import inductiva.client
from inductiva import constants
import inductiva.client.models
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import aiortc

# STUN/TURN server configuration
ICE_SERVERS = [
    aiortc.RTCIceServer("stun:" + constants.TURN_SERVER_URL),
    aiortc.RTCIceServer("turn:" + constants.TURN_SERVER_URL)
]

aiortc_logger = logging.getLogger("aioice")
aiortc_logger.setLevel(logging.WARNING)


class Operations(enum.Enum):
    LIST = "ls"
    TAIL = "tail"
    TOP = "top"
    LAST_MODIFIED_FILE = "last_modified_file"


class FileTracker:
    """File Tracker class for connecting to a running task via WebRTC."""

    def __init__(self):
        self.peer_connections = []

    async def setup_channel(self, pc, operation, follow=False, **kwargs):
        channel = pc.createDataChannel("file_transfer")
        queue = asyncio.Queue()
        end_event = asyncio.Event()

        @channel.on("open")
        def on_open():
            request = {"type": operation.value, "follow": follow}
            request["args"] = kwargs
            channel.send(json.dumps(request))

        @channel.on("message")
        async def on_message(message):
            await queue.put(json.loads(message))
            if not follow:
                end_event.set()

        @channel.on("close")
        async def on_close():
            await queue.put(None)

        return queue, end_event

    def create_peer_connection(self):
        pc = aiortc.RTCPeerConnection(
            aiortc.RTCConfiguration(iceServers=ICE_SERVERS))
        self.peer_connections.append(pc)
        return pc

    async def connect_to_task(
        self,
        api: inductiva.client.TasksApi,
        pc,
        task_id,
    ):
        connection_id = str(uuid.uuid4())
        api.register_task(
            task_connection=inductiva.client.models.TaskConnection(
                sender_id=connection_id),
            task_id=task_id,
        )

        offer = await pc.createOffer()
        await pc.setLocalDescription(offer)

        api.offer_task(
            task_connection=inductiva.client.models.TaskConnection(
                sender_id=connection_id,
                receiver_id=task_id,
                type="offer",
                sdp=pc.localDescription.sdp,
            ),
            task_id=task_id,
        )

        resp = api.get_message_without_preload_content(client=connection_id,
                                                       task_id=task_id)

        if resp.status == 204:
            return False

        data = json.loads(resp.data)
        if data["type"] == "answer":
            await pc.setRemoteDescription(
                aiortc.RTCSessionDescription(sdp=data["sdp"],
                                             type=data["type"]))

        return True

    async def cleanup(self):
        for pc in self.peer_connections:
            await pc.close()
