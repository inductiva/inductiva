import asyncio
import json
import uuid
from aiortc import RTCPeerConnection, RTCSessionDescription
import aiohttp
import enum
import logging

SIGNALING_SERVER = "http://34.79.246.4:6000"

# STUN/TURN server configuration
ICE_SERVERS = [
    {"urls": ["stun:34.79.246.4:3478"]},
    {"urls": ["turn:34.79.246.4:3478"]}
]

class Operations(enum.Enum):
    LIST = "ls"
    TAIL = "tail"

class FileTracker:
    def __init__(self):
        self.pc = RTCPeerConnection()
        self.pc.configuration = {"iceServers": ICE_SERVERS}
        self._message = None


    async def setup_channel(self, operation, **kwargs):
        channel = self.pc.createDataChannel("file_transfer")
        fut = asyncio.Future()

        @channel.on("open")
        def on_open():
            request = operation.value
            if kwargs:
                request += ":" + json.dumps(",".join(kwargs.values())).strip('"')
            channel.send(request)

        @channel.on("message")
        def on_message(message):
            self._message = json.loads(message)
            channel.close()

        @channel.on("close")
        def on_close():
            fut.set_result(self._message)

        return fut
    
    async def connect_to_task(self, task_id):
        client_id = str(uuid.uuid4())
        async with aiohttp.ClientSession() as session:
            await session.post(f"{SIGNALING_SERVER}/register", json={"clientId": client_id})
            
            offer = await self.pc.createOffer()
            await self.pc.setLocalDescription(offer)
            
            await session.post(f"{SIGNALING_SERVER}/offer", json={
                "receiverId": task_id,
                "senderId": client_id,
                "type": "offer",
                "sdp": self.pc.localDescription.sdp
            })

            async with session.get(f"{SIGNALING_SERVER}/message?clientId={client_id}") as resp:
                if resp.status == 200:
                    data = await resp.json()
                    if data['type'] == 'answer':
                        await self.pc.setRemoteDescription(RTCSessionDescription(sdp=data['sdp'], type=data['type']))