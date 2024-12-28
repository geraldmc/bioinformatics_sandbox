import os

from Bio import Entrez
from dotenv import load_dotenv 
import pandas as pd
import numpy as np

load_dotenv() # get env vars

# Use Biopython's Esearch with PubMed.
def search(query, retmax=10):
    Entrez.email = os.getenv("email")
    handle = Entrez.esearch(db='pubmed',
                            tool = "gene_genie",
                            api_key = os.getenv("NCBI_key"),
                            sort='relevance',
                            retmax=retmax,
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

# Use Biopython's Efetch with PubMed, for IDs from a previous Esearch.
def fetch(id_list):
    ids = ','.join(id_list)
    Entrez.email = os.getenv("email")
    handle = Entrez.efetch(db='pubmed',
                           tool = "gene_genie",
                           api_key = os.getenv("NCBI_key"),
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

# Use Biopython's ESummary with PubMed, for IDs from a previous Esearch.
def summary(id_list):
    ids = ','.join(id_list)
    Entrez.email = os.getenv("email")
    handle = Entrez.esummary(db='pubmed',
                           tool = "gene_genie",
                           api_key = os.getenv("NCBI_key"),
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

if __name__ == "__main__":
  ''' Run Esearch, retrieve from Efetch, get PubMed urls.
      For each url, parse and ... 
  '''
  title_list= []
  abstract_list=[]
  journal_list = []
  pubdate_list = []

  results = search('epigenome', retmax=3)
  id_list = results['IdList']

  articles = fetch(id_list)

for idx, article in enumerate(articles['PubmedArticle']):
  title_list.append(article['MedlineCitation']['Article']['ArticleTitle'].replace(u'\xa0', ' '))
  try:
    abstract_list.append(article['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
  except:
    abstract_list.append('No Abstract')
  journal_list.append(article['MedlineCitation']['Article']['Journal']['Title'])
  try:
    year = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
  except:
    year = 'No Data'
  try:
    month = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
  except:
    month = 'No Data'
  pubdate_list.append(year + '/' + month)

df = pd.DataFrame(list(zip(title_list, abstract_list, journal_list, pubdate_list)), 
                  columns=['Title', 'Abstract', 'Journal', 'PubDate'])