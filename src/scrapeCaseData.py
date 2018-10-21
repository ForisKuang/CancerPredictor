import requests
import json

cases_endpt = 'https://api.gdc.cancer.gov/cases'

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

fields = ','.join(fields)
filters = {
    "op": "in",
    "content": {
        "field": "primary_site",
        "value": ["Lung"]
    }
}

params= {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "TSV",
    "size": 580
}

response = requests.get(cases_endpt, params=params)

f = open('LungData.tsv', "w+")
f.write(response.content)
