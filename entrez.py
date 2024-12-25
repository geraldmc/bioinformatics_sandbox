import os
import time
import warnings
import json
import requests
from urllib.parse import urlencode
from dotenv import load_dotenv, dotenv_values 

# ------------------------------------------------------------------------------
# NCBI/Entrez - Entrez functions based on Biopython but using the Requests lib.
# NCBI defaults to XML. Try to use json wherever possible.
# ------------------------------------------------------------------------------

load_dotenv() # get env vars

email = 'gerald.mccollam@gmail.com'
max_tries = 3
sleep_between_tries = 15
tool = "gene_genie"
api_key = os.getenv("NCBI_key")
local_cache = None

def esearch(db, term, **keywords):
  '''This function searches and retrieves primary IDs (for use in EFetch,
     ELink and ESummary).
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  variables = {"db": db, "term": term}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))
  return response

def efetch(db, **keywords):
  '''This function retrieves records in the requested format from a list or
     set of one or more UIs or from user's environment.
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  variables = {"db": db}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))

def esummary(**keywords):
  '''This function retrieves document summaries from a list of primary IDs or
     from the user's environment.
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
  variables = {}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))

def einfo(**keywords):
  '''This function returns a summary of the Entrez databases as a results handle.
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
  variables = {}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))

def elink(**keywords):
  '''This function checks for the existence of an external or Related Articles
     link from a list of one or more primary IDs;  retrieves IDs and relevancy
     scores for links to Entrez databases or Related Articles.
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
  variables = {}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))

def epost(db, **keywords):
  '''Posts a file containing a list of primary IDs for future use.
  '''
  base_cgi = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi"
  variables = {"db": db}
  variables.update(keywords)
  request = _construct(base_cgi, variables)
  try:
    response = requests.get(base_cgi, request)  
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))
def _get_params(params, join_ids=True):
  pass

def _format_ids(ids):
    '''
    Input: a single ID (int or str), or iterable of strings/ints,
    or a string of IDs separated by commas.
    '''
    if isinstance(ids, int):
        return str(ids) # Single int, convert it to str.

    if isinstance(ids, str):
        # Multiple IDs, remove white space.
        return ",".join(id.strip() for id in ids.split(","))

    # Not a string or integer, assume iterable
    return ",".join(map(str, ids))

def _has_api_key(request):
  pass

def _construct(base_cgi, params=None, join_ids=True):
  '''
  '''
  if params is None:
    params = {}
  params.setdefault("tool", tool)
  params.setdefault("email", email)
  params.setdefault("api_key", api_key)
  params.setdefault("usehistory", 'y')
  params.setdefault("retmode", 'json')

  if join_ids and "id" in params:
      params["id"] = _format_ids(params["id"])

  return params

def _process(request):
  return request

def test_handle_exception(resp):
 try:
    resp.raise_for_status()
 except requests.exceptions.RequestException as e:
    print("Request failed:", str(e))

if __name__ == "__main__":
  resp = esearch(db="pubmed", retmax=5, term="cancer[mesh] epigenomics[mesh] 2024[pdat]")
  resp.json()