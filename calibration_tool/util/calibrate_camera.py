# %%
import json
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from arguslib import UndistortedCamera
from arguslib.instruments.instruments import rotation_matrix_i_to_g

CAMERA_ID = 'COBALT_5-4'
DATASET_FILE = f'calibration_dataset_{CAMERA_ID}.json'

def load_data_and_vectors():
    print(f"Initializing Camera {CAMERA_ID} to recover vectors...")
    cam_orig = UndistortedCamera.from_config(*CAMERA_ID.split('_'), calibration_images=True)
    R_orig = rotation_matrix_i_to_g(*cam_orig.rotation)

    with open(DATASET_FILE, 'r') as f:
        raw_data = json.load(f)

    observations = []
    print(f"Recovering global star vectors for {len(raw_data)} points...")

    for entry in raw_data:
        expected_px = np.array(entry['expected_px'])
        expected_px_intrinsic = expected_px * cam_orig.scale_factor
        v_inst = cam_orig.intrinsic.image_to_view(expected_px_intrinsic, norm=True)
        v_global = R_orig @ v_inst
        observations.append({
            'timestamp': entry['timestamp'],
            'v_global': v_global,
            'actual_px': np.array(entry['actual_px'])
        })
    return observations, cam_orig

def objective_function(params, cam, observations):
    el, az, roll = params
    try:
        R_new = rotation_matrix_i_to_g(el, az, roll)
    except Exception:
        return 1e6

    R_global_to_inst = R_new.T

    errors = []
    for obs in observations:
        v_global = obs['v_global']
        v_inst_new = R_global_to_inst @ v_global
        if v_inst_new[2] <= 0:
            errors.append(1000.0)
            continue
        proj_px_intrinsic = cam.intrinsic.view_to_image(v_inst_new)
        proj_px = proj_px_intrinsic / cam.scale_factor
        dist = np.linalg.norm(proj_px - obs['actual_px'])
        errors.append(dist)
    return np.mean(errors)

def run_parameter_sweep(observations, cam):
    print("\n--- Starting Grid Search ---")
    azimuths = np.linspace(0, 360, 72)
    elevations = [80, 85, 90]
    rolls = [-180, -90, 0, 90, 180]

    best_error = float('inf')
    best_params = cam.rotation

    plot_data_az = []
    plot_data_err = []

    for el in elevations:
        for az in azimuths:
            local_best_err = float('inf')
            for r in rolls:
                err = objective_function([el, az, r], cam, observations)
                if err < local_best_err:
                    local_best_err = err
                    if err < best_error:
                        best_error = err
                        best_params = [el, az, r]
                        print(f"New Best: El={el}, Az={az:.1f}, Roll={r} -> Err={err:.2f} px")
            plot_data_az.append(az)
            plot_data_err.append(local_best_err)

    plt.figure(figsize=(10, 5))
    plt.scatter(plot_data_az, plot_data_err, s=10, alpha=0.7)
    plt.xlabel('Azimuth (deg)')
    plt.ylabel('Error (px)')
    plt.title('Coarse Grid Search Results')
    plt.grid(True)
    plt.show()

    return best_params

def fine_tune(start_params, observations, cam):
    print(f"\n--- Fine Tuning via Nelder-Mead from {start_params} ---")
    res = minimize(
        objective_function,
        start_params,
        args=(cam, observations),
        method='Nelder-Mead',
        tol=1e-4,
        options={'maxiter': 500}
    )
    print("Optimization Complete.")
    print(f"Final Params: {res.x}")
    print(f"Final Mean Error: {res.fun:.4f} pixels")
    return res.x

if __name__ == "__main__":
    obs, cam = load_data_and_vectors()
    if not obs:
        print("No calibration data found. Run the flask app first!")
        exit()

    best_coarse = run_parameter_sweep(obs, cam)
    final_params = fine_tune(best_coarse, obs, cam)
    print("\n" + "="*40)
    print("RECOMMENDED CONFIGURATION UPDATE")
    print("="*40)
    print(f"elevation: {final_params[0]:.4f}")
    print(f"azimuth:   {final_params[1]:.4f}")
    print(f"roll:      {final_params[2]:.4f}")
    print("="*40)
# %%
