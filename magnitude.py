"""Local Magnitude
 began life as the Richter magnitude, however, this is a magnitude scale with specific scaling parameters for Southern California. The same equation form used in the initial Richter scale is used elsewhere, but with appropriate attenuation parameters: if the original Southern California Richter scale is used in New Zealand it doesn't give the correct result (as well will see when we try to use it!).

 is based on the measurement of the peak amplitude in a waveform recorded by a Wood-Anderson seismometer. There are very few Wood-Anderson seismometers still in operation today, but we know what the instrument response of a Wood-Anderson instrument is, and can therefore simulate Wood-Anderson waveforms from other instruments by correcting for the original instruments response (as above) and convolving the resulting ground motion with the response of a Wood-Anderson instrument.

The original Richter local magnitude scale is:

where
 is the peak amplitude (in micro-meters) on a Wood-Anderson seismometer,
 is the epicentral distance. This only works well for shallow earthquakes. The second term effectively corrects for geometrical spreading, while the final factor is included to scale the magnitude.

The adapted local magnitude scale for New Zealand is:
where S is a station correction term. See Ristau, 2009 for a discussion of some of the magnitudes in use in New Zealand.

Recent work by Michailos et al., 2019 derived a robust local magnitude scale for the Southern Alps tied to moment magnitude:

where
 is the hypocentral distance,
 is an attenuation parameter, which, for
 is
 and for
 is
.
 is a site dependent correction which incorporates local site effects. Geometrical spreading is assumed to be logarithmically related to
, with a gradient of 1.

The reason that we end up with multiple magnitude scales is similar to the reason we need to include station correction terms: we don't do a good job of correcting for path effects. When the paths change, the geology the waves encounter is different and so the attenuation is different. Although we correct for some of the attenuation in the attenuation correction terms, these don't capture the 3D range of variability we see in the geology. This is also part of the reason why we get different magnitude estimates from different stations (another major factor is that amplitudes vary with azimuth from the source, as well as directivity affects).

Lets look at how we would pick an amplitude for local magnitude for a nearby earthquake. We will use the seismic_picker.py applet again which has been extended to allow amplitude and duration picks to be made."""

from obspy.clients.fdsn import Client
from gphs445_utilities.plot_event import get_geonet_waveforms

client = Client("GEONET")
event = client.get_events(eventid="2019p304574")[0]
# Lets just use the five closest weak motion stations
clean_event = event.copy()
clean_event.picks = [
    p for p in sorted(event.picks, key=lambda p: p.time)
    if p.waveform_id.channel_code[0] != "B"][0:5]

# We want to remove the amplitudes already picked, and magnitudes so that we can overwrite them with our own.
clean_event.amplitudes = []
clean_event.station_magnitudes = []

st = get_geonet_waveforms(clean_event, all_components=True)
st = st.detrend().taper(0.05)
#st.plot()  # For some reason on Windows I need to plot in this first cell to get later plots to work.
# We need to response information for the stations so that we can correct
# the repsonse and simulate a wood anderson instrument.
bulk = [
    (p.waveform_id.network_code, p.waveform_id.station_code,
     p.waveform_id.location_code, p.waveform_id.channel_code[0:-1] + "?",
     p.time - 60, p.time + 200) for p in clean_event.picks]
inv = client.get_stations_bulk(bulk, level="response")

paz_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
          'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
st_wa = st.copy()
st_wa = st_wa.remove_response(inv, "DISP")
st_wa = st_wa.simulate(paz_remove=None, paz_simulate=paz_WA)

import matplotlib.pyplot as plt
from gphs445_utilities.seismic_picker import SeismicPicker

picker = SeismicPicker(st_wa, event_in=clean_event)
event_out = SeismicPicker(st_wa, event_in=clean_event).pick()
print(event_out)

event_picked = event_out.copy()
# Quick plot to check that the picks are recorded in the right place
SeismicPicker(st_wa, event_picked).show()


for amplitude in event_picked.amplitudes:
    if amplitude.type != "END":
        print("Amplitude: {0:.2g} m".format(amplitude.generic_amplitude))

from gphs445_utilities.location import Geographic
from math import log10


def _distance(point_1, point_2):
    """
    Calcuate hypocentral distance from Geographic points

    :type point_1: `coordinates.Geographic`
    :type point_2: `coordinates.Geographic`

    :returns: float
    """
    point_2_xyz = point_2.to_xyz(origin=point_1, strike=0, dip=90)
    return (point_2_xyz.x ** 2 + point_2_xyz.y ** 2 + point_2_xyz.z ** 2) ** 0.5


origin = clean_event.preferred_origin()
origin = Geographic(
    latitude=origin.latitude, longitude=origin.longitude,
    depth=origin.depth / 1000.)
magnitude = 0
used_station_count = 0
for amplitude in event_picked.amplitudes:
    if amplitude.type == 'END':
        continue
    pick = amplitude.pick_id.get_referred_object()
    station_loc = inv.get_coordinates(pick.waveform_id.get_seed_string(),
                                      pick.time)
    station_loc = Geographic(
        latitude=station_loc["latitude"], longitude=station_loc["longitude"],
        depth=(station_loc["local_depth"] - station_loc["elevation"]) / 1000.)
    distance = _distance(origin, station_loc)
    print("Amplitude {0:.2g} m at {1:.2f} km".format(
        amplitude.generic_amplitude, distance))
    station_magnitude = (
            log10(abs(amplitude.generic_amplitude) * 1e6) + log10(distance) -
            0.0029 * distance + 0)
    print("Using the Richter scale gives: {0:.2f}".format(station_magnitude))
    magnitude += station_magnitude
    used_station_count += 1

if used_station_count == 0:
    print(f"No amplitude picks found")
else:
    magnitude /= used_station_count
    print("Average magnitude: {0:.2f}".format(magnitude))
    print("GeoNet magnitude: {0:.2f}".format(clean_event.preferred_magnitude().mag))

