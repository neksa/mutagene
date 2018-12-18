import requests
import json


TOKEN_FILE_PATH = "/Users/agoncear/data/gdc-user-token.txt"


def gdc_read_file(file_id="11443f3c-9b8b-4e47-b5b7-529468fec098"):
    data_endpt = "https://api.gdc.cancer.gov/slicing/view/{}".format(file_id)

    with open(TOKEN_FILE_PATH, "r") as token:
        token_string = str(token.read().strip())

    params = {"gencode": ["BRCA1", "BRCA2"]}

    response = requests.post(
        data_endpt,
        data=json.dumps(params),
        headers={
            "Content-Type": "application/json",
            "X-Auth-Token": token_string
        })
    return response.content
