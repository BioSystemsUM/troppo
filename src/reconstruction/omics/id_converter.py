"""
 Created by Jorge Gomes on 06/06/2018
 id_converter

"""
import pandas as pd
import urllib.request
import urllib.error
import shutil
from pathlib import Path
from datetime import date


"""
This python file contains two functions that rely on the HGNC complete set file.

idConverter:
    This function converts the ids from a given omics dataset into the desired ones to better match a metabolic model.
    Conversion is done based on the HGNC database.
    
searchNomenclature:
    This function searches which gene identification nomenclature is used on the provided ids.
"""


def _get_HGNC():
    # Download the file from `url` and save it locally under `hgnc_complete_set_[Date]`:
    now = date.today()
    url = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
    path = "hgnc_complete_set_" + str(now) + ".tsv"
    file = Path(path)

    if file.is_file():  # if file has already been downloaded today skip this step
        return path
    else:
        try:
            file = "hgnc_complete_set_" + str(now) + ".tsv"
            with urllib.request.urlopen(url) as response, open(file, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

        except urllib.error.URLError as e:
            file = "hgnc_complete_set_2018-09-13.tsv"

            print('Please Check Internet Connection, using locally available HGNC file')
        return file


def idConverter(ids, old, new):
    """
    This function converts the ids from a given omics dataset into the desired ones to better match a metabolic model.
    Conversion is done based on the HGNC database.

    Args:
        ids: list or set, containing the ids to be converted
        old: string - exact match, the nomenclature designation of the input IDS.
                      Must be different from new and contained in NOMENCLATURES
        new: string - exact match, the nomenclature designation of the output IDs.
                      Must be different from old and contained in NOMENCLATURES

    NOMENCLATURES:["hgnc_id","symbol","name","entrez_id","ensembl_gene_id","vega_id","ucsc_id","ccds_id",
                   "uniprot_ids","pubmed_id","omim_id","locus_group","locus_type","alias_symbol","alias_name",
                   "prev_symbol","prev_name","ena","refseq_accession","rna_central_ids"]

    Returns:
        List of converted ids.
    """
    file = _get_HGNC()
    ds = pd.read_csv(file, sep='\t', low_memory=False)
    try:
        d = dict(zip(ds[old.lower()], ds[new.lower()]))
    except KeyError:
        print('The new ID designation is incorrect, must be one of the following:\n',
              ds.columns.values.tolist())
        return

    res = {x: str(d[x]).split('.')[0] for x in ids if x in d}
    return res


def searchNomenclature(ids):
    """
    This function searches which gene identification nomenclature is used on the provided ids.
    When ids from different nomenclatures are input, the result will be the nomenclature with the most matches.
    Also handles cases where some ids do not match but others do.

    Args:
        ids: list, list of ids (all using the same nomenclature)

    Returns: string, the nomenclature designation according to HGNC complete set table.
    """

    file = _get_HGNC()
    found = False  # some ids may not be contained in HGNC
    nomenclature_col = None
    matches = {}  # workaround for mixed ids

    while not found:  # cross each id with every line till a match comes up
        test = ids.pop()
        with open(file, 'r', encoding='utf8') as f:
            for line in f.readlines():
                if test in line.split('\t'):  # ensures an exact match
                    nomenclature_col = line.split('\t').index(test)
                    if nomenclature_col not in matches:
                        matches[nomenclature_col] = 1
                        break
                    else:
                        matches[nomenclature_col] += 1
                        break
            if matches != {}:
                if len(ids) < 10:
                    threshold = 1
                else:
                    threshold = 10
                if matches[max(matches, key=matches.get)] >= threshold:
                    found = True
                    nomenclature_col = [x for x, y in matches.items() if y == matches[max(matches, key=matches.get)]][0]

            if len(ids) == 0 and not found:
                print('No match was found for the provided ids')
                break

    ds = pd.read_csv(file, sep='\t', low_memory=False)
    if found:
        return ds.columns.values[nomenclature_col]
    else:
        return None


if __name__ == "__main__":
    a = ['A1BG','HGNC:1']
    print(idConverter(a, 'symbol', 'entrez_id'))