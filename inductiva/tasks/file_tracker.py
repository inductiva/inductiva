"""File Tracker module for connecting to a running task via WebRTC."""
import asyncio
import json
import uuid
import enum
import logging

aiortc_imported = True

try:
    import aiortc
except ImportError:
    aiortc_imported = False

# STUN/TURN server configuration
if aiortc_imported:
    ICE_SERVERS = [
        aiortc.RTCIceServer("stun:webrtc.inductiva.ai:3478"),
        aiortc.RTCIceServer("turn:webrtc.inductiva.ai:3478")
    ]

aiortc_logger = logging.getLogger("aioice")
aiortc_logger.setLevel(logging.WARNING)


class Operations(enum.Enum):
    LIST = "ls"
    TAIL = "tail"


class FileTracker:
    """File Tracker class for connecting to a running task via WebRTC."""

    def __init__(self):
        if not aiortc_imported:
            raise NotImplementedError("Feature not available for this version.")

        self.pc = aiortc.RTCPeerConnection(
            aiortc.RTCConfiguration(iceServers=ICE_SERVERS))

    async def setup_channel(self, operation, follow=False, **kwargs):
        channel = self.pc.createDataChannel("file_transfer")
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

    async def connect_to_task(self, api, task_id):
        connection_id = str(uuid.uuid4())
        path_params = {"task_id": task_id}
        api.register_task(body={"sender_id": connection_id},
                          path_params=path_params)

        offer = await self.pc.createOffer()
        await self.pc.setLocalDescription(offer)

        api.offer_task(body={
            "sender_id": connection_id,
            "receiver_id": task_id,
            "type": "offer",
            "sdp": self.pc.localDescription.sdp
        },
                       path_params=path_params)

        resp = api.get_message(query_params={"client": connection_id},
                               path_params=path_params)

        if resp.response.status == 204:
            return False

        data = resp.body
        if data["type"] == "answer":
            await self.pc.setRemoteDescription(
                aiortc.RTCSessionDescription(sdp=data["sdp"],
                                             type=data["type"]))

        return True

    async def cleanup(self):
        await self.pc.close()