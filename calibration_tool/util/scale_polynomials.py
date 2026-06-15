# %%



def scale_camera_coefficients(angle_to_r, r_to_z, scale_factor):
    # 1. Scale angle_to_radius: Simply multiply all coefficients by S
    # r_new = S * r_old
    new_angle_to_r = [c * scale_factor for c in angle_to_r]

    # 2. Scale radius_to_z: Apply the power rule
    # Coeff_new[i] = Coeff_old[i] / S^(i-1)
    new_r_to_z = []
    for i, c in enumerate(r_to_z):
        exponent = i - 1
        divisor = scale_factor ** exponent
        new_r_to_z.append(c / divisor)

    return new_angle_to_r, new_r_to_z

# --- Configuration ---
SCALE_FACTOR = 1.02


input = """
poly_incident_angle_to_radius:
- 0.0
- 980.8407917545795
- -18.09642824474127
- 45.84020924445242
- -36.58603464202598
poly_radius_to_z:
- 982.0389053824945
- 0.0
- -0.00042127044198000447
- 1.8698855138170602e-07
- -1.4664704293655303e-10
"""

poly_incident_angle_to_radius = input.split("poly_incident_angle_to_radius:")[1].split("poly_radius_to_z:")[0].strip().split("\n")
poly_incident_angle_to_radius = [float(c[2:]) for c in poly_incident_angle_to_radius]
poly_radius_to_z = input.split("poly_radius_to_z:")[1].strip().split("\n")
poly_radius_to_z = [float(c[2:]) for c in poly_radius_to_z]


# --- Execution ---
new_ar, new_rz = scale_camera_coefficients(
    poly_incident_angle_to_radius, 
    poly_radius_to_z, 
    SCALE_FACTOR
)

# --- Print Output ---
print("poly_incident_angle_to_radius:")
for c in new_ar:
    print(f"- {c}")
print("poly_radius_to_z:")
for c in new_rz:
    print(f"- {c}")
# %%
