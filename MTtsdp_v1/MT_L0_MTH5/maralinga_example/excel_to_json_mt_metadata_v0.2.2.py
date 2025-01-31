import pandas as pd
import numpy as np
import json
import glob
import os

excel_spreadsheet = '/scratch/abc/abc123/.../SA_Maralinga_demo_mt_metadata.xlsx'
station_json_from_spreadsheet_directory = '/scratch/abc/abc123/.../station_json_outdir/import/'
json_export_directory = '/scratch/abc/abc123/.../station_json_outdir/export/'
base_json_file = '/scratch/abc/abc123/.../base_v0.2.1.json'


if not os.path.exists(station_json_from_spreadsheet_directory):
    os.makedirs(station_json_from_spreadsheet_directory)
else:
    print(f'Directory "{station_json_from_spreadsheet_directory}" already exists. Continuing...')

if not os.path.exists(json_export_directory):
    os.makedirs(json_export_directory)
else:
    print(f'Directory "{json_export_directory}" already exists. Continuing...')


df = pd.read_excel(excel_spreadsheet)

for dataframe in df.iloc:
    filename = dataframe["station.id"]+".json"
    dataframe.to_json(station_json_from_spreadsheet_directory+filename)#,orient='records')


station_import_json_files = sorted(glob.glob(station_json_from_spreadsheet_directory+'*.json'))


for station in station_import_json_files:
    outfilename = os.path.basename(station)
    outfile_fullpath = json_export_directory+outfilename
    with open(station, 'r') as fh:
        json_export = json.load(fh)
    with open(base_json_file, 'r') as gh:
        json_data = json.load(gh)
    for key, value in json_data.items():
        if key == "survey":
            value["acquired_by.author"] = json_export["survey.acquired_by.author"]
            value["acquired_by.comments"] = json_export["survey.acquired_by.comments"]
            value["acquired_by.organization"] = json_export["survey.acquired_by.organization"]
            value["citation_dataset.doi"] = json_export["survey.citation_dataset.doi"]
            value["citation_journal.doi"] = json_export["survey.citation_journal.doi"]
            value["comments"] = json_export["survey.comments"]
            value["country"] = json_export["survey.country"]
            value["datum"] = json_export["survey.datum"]
            value["funding_source.organization"] = json_export["survey.funding_source.organization"]
            value["name"] = json_export["survey.name"]
            value["project"] = json_export["survey.project"]
            value["project_lead.author"] = json_export["survey.project_lead.author"]
            value["project_lead.email"] = json_export["survey.project_lead.email"]
            value["project_lead.organization"] = json_export["survey.project_lead.organization"]
            value["release_license"] = json_export["survey.release_license"]
            value["state"] = json_export["survey.state"]
            value["summary"] = json_export["survey.summary"]
            value["time_period.end_date"] = json_export["survey.time_period.end_date"]
            value["time_period.start_date"] = json_export["survey.time_period.start_date"]

        if key == "station":
            value["acquired_by.name"] = json_export["station.acquired_by.name"]
            value["acquired_by.comments"] = json_export["station.acquired_by.comments"]
            value["channel_layout"] = json_export["station.channel_layout"]
            #value["channels_recorded"] = json_export["station.channels_recorded"]
            value["channels_recorded"] = ["Ex, Ey, Bx, By, Bz"]
            value["comments"] = json_export["station.comments"]
            value["data_type"] = json_export["station.data_type"]
            value["geographic_name"] = json_export["station.geographic_name"]
            value["id"] = json_export["station.id"]
            value["location.declination.model"] = json_export["station.location.declination.model"]
            value["location.declination.value"] = json_export["station.location.declination.value"]
            value["location.declination.comments"] = json_export["station.location.declination.comments"]
            value["location.elevation"] = json_export["station.location.elevation"]
            value["location.latitude"] = json_export["station.location.latitude"]
            value["location.longitude"] = json_export["station.location.longitude"]
            value["orientation.reference_frame"] = json_export["station.orientation.reference_frame"]
            value["release_license"] = json_export["station.release_license"]
            value["run_list"] = json_export["station.run_list"]
            value["time_period.end"] = json_export["station.time_period.end"]
            value["time_period.start"] = json_export["station.time_period.start"]

        elif key == "run":
            value["acquired_by.author"] = json_export["run.acquired_by.author"]
            #value["channels_recorded_electric"] = json_export["run.channels_recorded_electric"]
            value["channels_recorded_electric"] = ["Ex, Ey"]
            #value["channels_recorded_magnetic"] = json_export["run.channels_recorded_magnetic"]
            value["channels_recorded_magnetic"] = ["Bx, By, Bz"]
            value["data_logger.id"] = json_export["run.data_logger.id"]
            value["data_logger.manufacturer"] = json_export["run.data_logger.manufacturer"]
            value["data_logger.model"] = json_export["run.data_logger.model"]
            value["data_logger.power_source.type"] = json_export["run.data_logger.power_source.type"]
            value["data_logger.timing_system.type"] = json_export["run.data_logger.timing_system.type"]
            value["data_logger.type"] = json_export["run.data_logger.type"]
            value["data_type"] = json_export["run.data_type"]
            value["id"] = json_export["run.id"]
            value["metadata_by.author"] = json_export["run.metadata_by.author"]
            value["sample_rate"] = json_export["run.sample_rate"]
            value["time_period.end"] = json_export["run.time_period.end"]
            value["time_period.start"] = json_export["run.time_period.start"]

        elif key == "electric_ex":
            value["dipole_length"] = json_export["electric.dipole_length_ex"]
            value["measurement_azimuth"] = json_export["electric.measurement_azimuth_ex"]
            value["negative.manufacturer"] = json_export["electric.negative.manufacturer"]
            value["negative.model"] = json_export["electric.negative.model"]
            value["negative.type"] = json_export["electric.negative.type"]
            value["positive.manufacturer"] = json_export["electric.positive.manufacturer"]
            value["positive.model"] = json_export["electric.positive.model"]
            value["positive.type"] = json_export["electric.positive.type"]
            value["sample_rate"] = json_export["electric.sample_rate_ex"]
            value["time_period.end"] = json_export["electric.time_period.end"]
            value["time_period.start"] = json_export["electric.time_period.start"]
            value["units"] = json_export["electric.units"]

        elif key == "electric_ey":
            value["dipole_length"] = json_export["electric.dipole_length_ey"]
            value["measurement_azimuth"] = json_export["electric.measurement_azimuth_ey"]
            value["negative.manufacturer"] = json_export["electric.negative.manufacturer"]
            value["negative.model"] = json_export["electric.negative.model"]
            value["negative.type"] = json_export["electric.negative.type"]
            value["positive.manufacturer"] = json_export["electric.positive.manufacturer"]
            value["positive.model"] = json_export["electric.positive.model"]
            value["positive.type"] = json_export["electric.positive.type"]
            value["sample_rate"] = json_export["electric.sample_rate_ey"]
            value["time_period.end"] = json_export["electric.time_period.end"]
            value["time_period.start"] = json_export["electric.time_period.start"]
            value["units"] = json_export["electric.units"]

        elif key == "magnetic_bx":
            value["measurement_azimuth"] = json_export["magnetic.measurement_azimuth_bx"]
            value["sample_rate"] = json_export["magnetic.sample_rate_bx"]
            value["sensor.id"] = json_export["magnetic.sensor.id_bx"]
            value["sensor.manufacturer"] = json_export["magnetic.sensor.manufacturer"]
            value["sensor.model"] = json_export["magnetic.sensor.model"]
            value["sensor.type"] = json_export["magnetic.sensor.type"]
            value["time_period.end"] = json_export["magnetic.time_period.end"]
            value["time_period.start"] = json_export["magnetic.time_period.start"]
            value["translated_azimuth"] = json_export["magnetic.translated_azimuth"]
            value["units"] = json_export["magnetic.units"]

        elif key == "magnetic_by":
            value["measurement_azimuth"] = json_export["magnetic.measurement_azimuth_by"]
            value["sample_rate"] = json_export["magnetic.sample_rate_by"]
            value["sensor.id"] = json_export["magnetic.sensor.id_by"]
            value["sensor.manufacturer"] = json_export["magnetic.sensor.manufacturer"]
            value["sensor.model"] = json_export["magnetic.sensor.model"]
            value["sensor.type"] = json_export["magnetic.sensor.type"]
            value["time_period.end"] = json_export["magnetic.time_period.end"]
            value["time_period.start"] = json_export["magnetic.time_period.start"]
            value["translated_azimuth"] = json_export["magnetic.translated_azimuth"]
            value["units"] = json_export["magnetic.units"]

        elif key == "magnetic_bz":
            value["measurement_azimuth"] = json_export["magnetic.measurement_azimuth_bz"]
            value["sample_rate"] = json_export["magnetic.sample_rate_bz"]
            value["sensor.id"] = json_export["magnetic.sensor.id_bz"]
            value["sensor.manufacturer"] = json_export["magnetic.sensor.manufacturer"]
            value["sensor.model"] = json_export["magnetic.sensor.model"]
            value["sensor.type"] = json_export["magnetic.sensor.type"]
            value["time_period.end"] = json_export["magnetic.time_period.end"]
            value["time_period.start"] = json_export["magnetic.time_period.start"]
            value["translated_azimuth"] = json_export["magnetic.translated_azimuth"]
            value["units"] = json_export["magnetic.units"]


    with open(outfile_fullpath,'w') as fp:
        json.dump(json_data, fp, indent=4, sort_keys=True)
