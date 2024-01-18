import websocket


def stream_task_logs(task_id):
    """Stream the logs of a task through a websocket.
    
    Args:
        task_id (int): ID of the task to print the logs of."""

    logs_url = "ws://127.0.0.1:5000/logs"

    websocket_url = f"{logs_url}/{task_id}"

    def on_message(ws, message):
        print(f"Received message: {message}")

    def on_error(ws, error):
        print(f"WebSocket Error: {error}")

    def on_close(ws, close_status_code, close_msg):
        print("Closed WebSocket")

    def on_open(ws):
        print("WebSocket connection opened")

    ws = websocket.WebSocketApp(websocket_url,
                                on_open=on_open,
                                on_message=on_message,
                                on_error=on_error,
                                on_close=on_close)

    ws.run_forever()
