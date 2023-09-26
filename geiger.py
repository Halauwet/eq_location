import numpy as np

# Simulated seismic data (Replace with real data)
seismometer_data = {
    "seismometer_1": {"latitude": 35.0, "longitude": -120.0, "arrival_time": 10.0},
    "seismometer_2": {"latitude": 36.0, "longitude": -119.0, "arrival_time": 9.8},
    "seismometer_3": {"latitude": 34.5, "longitude": -121.0, "arrival_time": 10.5},
}

# Speed of seismic waves (km/s)
wave_speed = 8.0

def estimate_epicenter(seismometer_data):
    # Convert data to arrays for easy calculations
    latitudes = np.array([data["latitude"] for data in seismometer_data.values()])
    longitudes = np.array([data["longitude"] for data in seismometer_data.values()])
    arrival_times = np.array([data["arrival_time"] for data in seismometer_data.values()])

    # Calculate time differences
    time_diffs = arrival_times - arrival_times.min()

    # Calculate distances from the epicenter (distance = time * speed)
    distances = time_diffs * wave_speed

    # Calculate the weighted average of latitudes and longitudes
    weighted_latitude = np.sum(latitudes / distances) / np.sum(1 / distances)
    weighted_longitude = np.sum(longitudes / distances) / np.sum(1 / distances)

    return weighted_latitude, weighted_longitude


estimated_latitude, estimated_longitude = estimate_epicenter(seismometer_data)
print("Estimated Epicenter Latitude:", estimated_latitude)
print("Estimated Epicenter Longitude:", estimated_longitude)