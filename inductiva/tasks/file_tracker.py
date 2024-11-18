import asyncio
import json
import uuid
from aiortc import RTCPeerConnection, RTCSessionDescription
import aiohttp
from asyncio import Future
import sys
import enum

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
    def init(self):
        self.peer_id = str(uuid.uuid4())
        self.pc = RTCPeerConnection()
        self._message = None


    async def create_peer_connection(self,operation, **kwargs):
        self.pc.configuration = {"iceServers": ICE_SERVERS}
        channel = self.pc.createDataChannel("file_transfer")
        fut = asyncio.Future()
        @channel.on("open")
        def on_open():
            request = operation.value
            if kwargs:
                request += ":" + json.dumps(",".join(kwargs.values())).strip('"')
            channel.send(request)

        @channel.on("message")
        async def on_message(message):
            if operation == Operations.LIST:
                self.message = json.loads(message)
            elif operation == Operations.TAIL:
                self.message = message.decode()
            channel.close()

        @channel.on("close")
        async def on_close():
            print("Channel closed")
            fut.set_result(self.message)

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