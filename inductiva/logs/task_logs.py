"""Module for streaming the logs of a task through a websocket."""
import websocket
import inductiva


def on_message(ws, message):
        print(f"Received message: {message}.")


def on_error(ws, error):
    print(f"WebSocket Error: {error}.")


def on_close(ws, status_code, message):
    print(f"Closed WebSocket with {status_code}: {message}.")


def on_open(ws):
    print("WebSocket connection opened.")


def stream_task_logs(task_id):
    """Stream the logs of a task through a websocket.
    
    Args:
        task_id (int): ID of the task to print the logs of."""

    websocket_url = f"{inductiva.logs_websocket}/{task_id}"

    ws = websocket.WebSocketApp(websocket_url,
                                on_open=on_open,
                                on_message=on_message,
                                on_error=on_error,
                                on_close=on_close)

    ws.run_forever()
