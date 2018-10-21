import requests
import json

# The base endpoint
endpt= "https://api.gdc.cancer.gov/"

"""
filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "files.access",
                "value": ["open"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_category",
                "value": ["Simple Nucleotide Variation"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_type",
                "value": ["Masked Somatic Mutation"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "program.name",
                "value": ["TCGA"]
            }
        }
    ]
}
"""

def download_files(filters, num_files):
    """
        params: filters
                    The filter for which to find files to download
                num_files
                    the number of files in which to download
        output: Generates all the downloaded files
    """
    import re
    files_endpt = endpt + "files"
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id",
        "format": "JSON",
        "size": num_files
    }
    response = requests.get(files_endpt, params = params)
    file_uuid_list = []
    for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
        file_uuid_list.append(file_entry["file_id"])

        data_endpt = "https://api.gdc.cancer.gov/data"

        params = {"ids": file_uuid_list}

        response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

        response_head_cd = response.headers["Content-Disposition"]

        file_name = re.findall("filename=(.+)", response_head_cd)[0]

        with open("datasets/" + file_name, "wb") as output_file:
            output_file.write(response.content)

"""
Sample fields for patient data

fields = [
        'case_id',
        'project.name',
        'project.project_id',
        'primary_site',
        'diagnoses.age_at_diagnosis',
        'diagnoses.tumor_stage',
        'demographic.gender',
        'demographic.race',
        'exposures.alcohol_history',
        'exposures.cigarettes_per_day',
        'exposures.years_smoked'
    ]

Sample filter for patient data

filters = {
    "op": "in",
    "content": {
        "field": "primary_site",
        "value": ["Lung"]
    }
}
"""

def fetch_patient_data(fields, size, data_format = "TSV", filters = None):
    """
        params: fields
                    The fields for which you want to grab certain patient data
                filters
                    The particular filter for which type of data you want to fetch
        returns: response data
                    The response data wwith all the values
    """
    cases_endpt = endpt + 'cases'
    fields = ','.join(fields)
    params = None
    if filters is None:
        params = {
            "fields": fields,
            "format": data_format,
            "size": size
        }
    else:
        params = {
            "filters": json.dumps(filters),
            "fields": fields,
            "format": data_format,
            "size": size
        }
    response = requests.get(cases_endpt, params=params)
    return response.text
