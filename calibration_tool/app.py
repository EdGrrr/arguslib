import os
import json
import datetime as dt
import numpy as np
from PIL import Image
import io
from flask import Flask, render_template, request, jsonify, send_file

from arguslib import UndistortedCamera
from arguslib.misc import plotting as pl

app = Flask(__name__)

CAMERA_ID = 'COBALT_5-3'

def dataset_file(camera_id: str) -> str:
    return f'calibration_dataset_{camera_id}.json'

def load_dataset(camera_id: str):
    path = dataset_file(camera_id)
    if os.path.exists(path):
        with open(path, 'r') as f:
            try:
                return json.load(f)
            except json.JSONDecodeError:
                return []
    return []

def save_dataset_all(camera_id: str, data):
    path = dataset_file(camera_id)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)

def toggle_dataset_entry(camera_id: str, new_entry):
    data = load_dataset(camera_id)
    existing_index = -1
    for i, item in enumerate(data):
        if (item['timestamp'] == new_entry['timestamp'] and
            item['star_index'] == new_entry['star_index']):
            existing_index = i
            break

    if existing_index != -1:
        del data[existing_index]
        action = 'removed'
    else:
        data.append(new_entry)
        action = 'added'

    save_dataset_all(camera_id, data)
    return action

def get_cam(camera_id: str):
    return UndistortedCamera.from_config(*camera_id.split('_'), calibration_images=True)

@app.route('/')
def index():
    start_time = dt.datetime(2025, 3, 2, 0).isoformat()
    return render_template('index.html', start_time=start_time, start_camera=CAMERA_ID)

@app.route('/api/image')
def get_image():
    ts_str = request.args.get('time')
    camera_id = request.args.get('camera', CAMERA_ID)
    ts = dt.datetime.fromisoformat(ts_str)

    cam = get_cam(camera_id)
    im_array = cam.get_data_time(ts)

    img = Image.fromarray(im_array)
    file_object = io.BytesIO()
    img.save(file_object, 'JPEG', quality=85)
    file_object.seek(0)
    return send_file(file_object, mimetype='image/jpeg')

@app.route('/api/stars')
def get_stars():
    ts_str = request.args.get('time')
    camera_id = request.args.get('camera', CAMERA_ID)
    ts = dt.datetime.fromisoformat(ts_str)

    try:
        cam = get_cam(camera_id)
        star_px = pl.get_star_pixel_coords(cam, ts)
        stars_list = star_px.tolist()
    except Exception as e:
        return jsonify({'error': str(e), 'stars': [], 'calibrations': {}})

    dataset = load_dataset(camera_id)
    calibrations = {}
    for entry in dataset:
        if entry.get('timestamp') == ts_str:
            calibrations[entry['star_index']] = entry['actual_px']

    return jsonify({'stars': stars_list, 'calibrations': calibrations})

@app.route('/api/submit', methods=['POST'])
def submit_point():
    payload = request.json
    camera_id = payload.get('camera', CAMERA_ID)
    entry = {
        'timestamp': payload['timestamp'],
        'star_index': payload['star_index'],
        'expected_px': payload['expected'],
        'actual_px': payload['actual'],
    }
    action = toggle_dataset_entry(camera_id, entry)
    return jsonify({'status': action, 'entry': entry})

@app.route('/api/clear_all', methods=['POST'])
def clear_all():
    payload = request.json or {}
    camera_id = payload.get('camera', CAMERA_ID)
    save_dataset_all(camera_id, [])
    return jsonify({'status': f'Cleared all detections for camera {camera_id}.'})

if __name__ == '__main__':
    app.run(debug=True, port=8875)