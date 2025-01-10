#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# Class to retrieve, parse and format entries from UNIPROT
#------------------------------------------------------------------------------

# API Documentation: https://www.uniprot.org/help/api
# API Examples: https://www.uniprot.org/help/api_retrieve_entries

import json, re
from typing import Union, Any

from requests import get
from requests.exceptions import RequestException
import xmltodict
from contextlib import closing

class UniProt(object):
  """USAGE (from app root):
  uni = UniProt_API(args.term)
  resp = uni.make_request_txt()
  terms = uni.extract_txt(resp)
  tags = uni.filter_tags(terms)
  print (json.dumps(tags))
  """

  def __init__(self, term, ext='.txt'):
    self.base = 'https://www.uniprot.org/uniprot/'
    self.extension = ext
    self.term = term
    self.url = self.base + term + self.extension

  @property
  def make_request_xml(self):
    """ Returns a dictionary of the retrieved xml.
    """
    #print (self.url)
    try:
      with closing(get(self.url, stream=True)) as resp: #returns b`xml`
        if self.is_good_enough_xml(resp):
          return xmltodict.parse(resp.content.decode("utf-8"))
        else:
          return None
    except RequestException as e:
      print(f'Error during requests to {self.url} : {str(e)}')
      return None

  @property
  def make_request_txt(self):
    """ Handling the returned content as text is easier. 
    """
    #print (self.url)
    try:
      with closing(get(self.url, stream=True)) as resp: 
        if self.is_txt(resp):
          return resp.content.decode("utf-8")
        else:
          return None
    except RequestException as e:
      print('Error during requests to {0} : {1}'.format(self.url, str(e)))
      return None
        
  def extract_xml(self, xml_dict):
    """
    """
    return_dict = {}
    for ref in xml_dict['uniprot']['entry']['dbReference']:
      try:
        return_dict[ref['@type']] = ref['@id']
      except KeyError:  # there may be mutiple entries per dbReference
        pass  # grab just the first one and go
    return_dict['ID'] = self.term
    return return_dict

  @staticmethod
  def is_txt(resp):
    return True # TBD!

  @staticmethod
  def is_good_enough_xml(resp):
    """
    Returns True if response looks like XML, otherwise False.
    """
    content_type = resp.headers['Content-Type'].lower()
    
    return (resp.status_code == 200 
            and content_type is not None 
            and content_type.find('xml') > -1)

  @staticmethod
  def extract_txt(resp):
    """ I love text!
    """
    txt_list = list()
    tags_dict = dict()
    li = resp.split('\n') # create a list
    for txt in li:
      if re.search(r'^ID', txt):
        tags_dict['ID'] = txt[5:].split()[0]
      elif re.search(r'^DE', txt):
        tags_dict['SubName'] = txt[5:].split(',')[0]
        #print (txt[5:].split(',')[0])
      elif re.search(r'^DR', txt):
        txt_list.append(txt[5:].split(';')[:2])
        #print (txt[5:].split(';')[:2])
    for item in txt_list:
      try:
        tags_dict[item[0]] = item[1].strip()
      except KeyError: # repeated elements, i.e. GO terms, grab just one
        pass
    return tags_dict

  @staticmethod
  def filter_tags(tags_dict):
    filtered_tags = dict()
    white_list = ['ID','EMBL','RefSeq','MaxQB', 'PeptideAtlas','ChEMBL','DrugBank', 'PIRSF',
                'Ensembl','GeneID', 'KEGG','HGNC','MIM','KO','Reactome','GO','Pfam',
                'OrthoDB','SMART','PROSITE','PRIDE', 'ProteomicsDB', 'HOGENOM', 'PANTHER']
    for tag in white_list:
      try:
        filtered_tags[tag] = tags_dict[tag]
      except KeyError:
        pass
    return filtered_tags


if __name__ == '__main__':
  import argparse
  import random
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--term', help='Search term')
  args = parser.parse_args()

  uniprot_id_list = ['B4DNQ5_HUMAN', 'B4DVM1_HUMAN', 'B5MC21_HUMAN', 
                      'IL6_HUMAN', 'Q75MH2_HUMAN']
  
  if args.term:
    print(f"Extracting text (as json) from a UniProt query, w/ command line args:\n")
    uni = UniProt(args.term)
    response = uni.make_request_txt
    terms = uni.extract_txt(response)
    tags = uni.filter_tags(terms)
    print (json.dumps(tags))
  else:
    print(f"Extracting text (as json) from a UniProt query, w/o command line args:\n")
    arg = random.choice(uniprot_id_list)
    uni = UniProt(arg)
    response= uni.make_request_txt
    terms = uni.extract_txt(response)
    tags = uni.filter_tags(terms)
    print (json.dumps(tags))

  if args.term:
    print(f"\nExtracting xml (as a dictionary) from a UniProt query, w/ command line args:\n")
    uni = UniProt(args.term, ext='.xml')
    xml_dict = uni.make_request_xml
    extracted = uni.extract_xml(xml_dict)
    print(extracted)
  else:
    print(f"\nExtracting xml (as a dictionary) from a UniProt query, w/o command line args:\n")
    arg = random.choice(uniprot_id_list)
    uni = UniProt(arg, ext='.xml')
    xml_dict = uni.make_request_xml
    extracted = uni.extract_xml(xml_dict)
    print(extracted)