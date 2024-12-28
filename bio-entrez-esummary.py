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
  ''' Run Esummary, get PubMed article data.
      For each url, parse and ... 
  '''
  PMID_list= []
  doi_list = []
  abstract_list=[]
  author_list = []
  journal_list = []
  pubdate_list = []

  results = search('epigenome', retmax=5) # throttle search here
  id_list = results['IdList']

  articles = summary(id_list)

  for article in articles:
    print(f"PubMed ID: {article['ArticleIds']['pubmed'][0]}")
    print(f"DOI: {article['ArticleIds']['doi']}")
    print(f"Title: {article['Title']}")
    print("Authors: ")
    print(*article['AuthorList'], sep="; ")
    print(f"Journal: {article['Source']}")
    print(f"PubDate: {article['EPubDate']}\n")

    PMID_list.append(article['ArticleIds']['pubmed'][0])
    doi_list.append(article['ArticleIds']['doi'])
    if int(articles[0]['HasAbstract']) == 1:
       abstract_list.append("True")
    else: abstract_list.append("False") 
    author_list.append(article['AuthorList'])
    journal_list.append(article['Source'])
    pubdate_list.append(article['EPubDate'])


df = pd.DataFrame(list(zip(PMID_list, doi_list, abstract_list, author_list, journal_list, pubdate_list)), 
                  columns=['PMID', 'DOI', 'Has_Abstract', 'Authors', 'Journal', 'PubDate'])