# -*- coding: utf-8 -*-
import json
import logging
import os
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from typing import Dict, List, Optional, Tuple, Union

import requests

from gseapy.utils import DEFAULT_CACHE_PATH, unique


def gsea_cls_parser(cls: str) -> Tuple[str]:
    """Extract class(phenotype) name from .cls file.

    :param cls: the a class list instance or .cls file which is identical to GSEA input .
    :return: phenotype name and a list of class vector.
    """

    if not isinstance(cls, str) and isinstance(cls, Iterable):
        classes = list(cls)
        sample_name = unique(classes)
    elif isinstance(cls, str):
        with open(cls) as c:
            file = c.readlines()
        classes = file[2].strip().split()
        sample_name = file[1].strip("#").strip().split()

        tmp = set(sample_name) & set(classes)
        if len(tmp) < 2:  # classes and sample_name are different
            s1 = classes[0]
            for i, c in enumerate(classes):
                if c == s1:
                    classes[i] = sample_name[0]
                else:
                    classes[i] = sample_name[1]
    else:
        raise Exception("Error parsing sample name!")

    if len(sample_name) != 2:
        raise Exception("Input groups have to be 2!")

    return sample_name[0], sample_name[1], classes


def gsea_edb_parser(results_path: str) -> Dict[str, List[str]]:
    """Parse results.edb file stored under **edb** file folder.

    :param results_path: the path of results.edb file.
    :return:
        a dict contains { enrichment_term: [es, nes, pval, fdr, fwer, hit_ind]}
    """

    xtree = ET.parse(results_path)
    xroot = xtree.getroot()
    res = {}
    # dict_keys(['RANKED_LIST', 'GENESET', 'FWER', 'ES_PROFILE',
    # 'HIT_INDICES', 'ES', 'NES', 'TEMPLATE', 'RND_ES', 'RANK_SCORE_AT_ES',
    # 'NP', 'RANK_AT_ES', 'FDR'])
    for node in xroot.findall("DTG"):
        enrich_term = node.attrib.get("GENESET").split("#")[1]
        es_profile = node.attrib.get("ES_PROFILE").split(" ")
        # esnull = term.get('RND_ES').split(" ")
        hit_ind = node.attrib.get("HIT_INDICES").split(" ")
        es_profile = [float(i) for i in es_profile]
        hit_ind = [int(i) for i in hit_ind]
        # esnull = [float(i) for i in esnull ]
        es = float(node.attrib.get("ES"))
        nes = float(node.attrib.get("NES"))
        pval = float(node.attrib.get("NP"))
        fdr = float(node.attrib.get("FDR"))
        fwer = float(node.attrib.get("FWER"))
        logging.debug("Enriched Gene set is: " + enrich_term)
        res[enrich_term] = [es, nes, pval, fdr, fwer, hit_ind]
    return res


def read_gmt(path: str) -> Dict[str, List[str]]:
    """Read GMT file

    :param str path: the path to a gmt file.
    :return: a dict object
    """
    if path.lower().endswith("gmt"):
        return get_library(name=path, min_size=0, max_size=100000, gene_list=None)
    else:
        raise ValueError("Please input a gmt file")
    return


def get_library(
    name: str,
    organism: str = "Human",
    min_size: int = 0,
    max_size: int = 2000,
    gene_list: Optional[List[str]] = None,
) -> Dict[str, List[str]]:
    """Parse gene_sets.gmt(gene set database) file or download from enrichr server.

    :param str name: the gene_sets.gmt file or an enrichr library name.
                    checkout full enrichr library name here: https://maayanlab.cloud/Enrichr/#libraries

    :param str organism: choose one from { 'Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm' }.
                         This arugment has not effect if input is a `.gmt` file.

    :param min_size: Minimum allowed number of genes for each gene set. Default: 0.
    :param max_size: Maximum allowed number of genes for each gene set. Default: 2000.

    :param gene_list: if input a gene list, min and max overlapped genes between gene set and gene_list are kept.

    :return dict: Return a filtered gene set database dictionary.

    Note: **DO NOT** filter gene sets, when use :func:`replot`. Because ``GSEA`` Desktop have already
    done this for you.

    """
    genesets_dict = {}
    if name.lower().endswith(".gmt"):
        logging.info("User Defined gene sets is given.......continue..........")
        with open(name) as genesets:
            for line in genesets:
                entries = line.strip().split("\t")
                key = entries[0]
                genesets_dict[key] = entries[2:]
    else:
        # get gene sets from enrichr libary
        names = get_library_name(organism=organism)
        if name in names:
            logging.info("Downloading and generating Enrichr library gene sets...")
            genesets_dict = download_library(name, organism=organism)
        else:
            raise ValueError(
                "Sorry. The input: %s could be be found given organism: %s"
                % (name, organism)
            )

    # filtering gene_sets
    total = len(genesets_dict)
    keys = list(genesets_dict.keys())
    if gene_list is None:
        for k in keys:
            if min_size <= len(genesets_dict[k]) <= max_size:
                continue
            del genesets_dict[k]
    else:
        # given a gene_list, filter gene sets by gene_overlap numbers
        gene_dict = {g: i for i, g in enumerate(gene_list)}
        for subset in keys:
            subset_list = set(genesets_dict[subset])  # remove duplicates
            # drop genes not found in the gene_dict
            gene_overlap = [g for g in subset_list if g in gene_dict]
            tag_len = len(gene_overlap)
            if (min_size <= tag_len <= max_size) and tag_len < len(gene_list):
                # tag_len should < gene_list
                genesets_dict[subset] = gene_overlap
                continue
            del genesets_dict[subset]

    filsets_num = total - len(genesets_dict)
    if filsets_num > 0:
        logging.info(
            "%04d gene_sets have been filtered out when max_size=%s and min_size=%s"
            % (filsets_num, max_size, min_size)
        )

    if filsets_num == len(genesets_dict):
        raise Exception(
            "No gene sets passed throught filtering condition!!!, try new paramters again!\n"
            + "Note: Gene names for gseapy is case sensitive."
        )

    return genesets_dict


def get_library_name(organism: str = "Human") -> List[str]:
    """return enrichr active enrichr library name.
    see also: https://maayanlab.cloud/modEnrichr/

    :param str organism: Select one from { 'Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm' }
    :return: a list of enrichr libraries from selected database

    """
    default = [
        "human",
        "mouse",
        "hs",
        "mm",
        "homo sapiens",
        "mus musculus",
        "h. sapiens",
        "m. musculus",
    ]
    _organisms = {
        "Fly": ["fly", "d. melanogaster", "drosophila melanogaster"],
        "Yeast": ["yeast", "s. cerevisiae", "saccharomyces cerevisiae"],
        "Worm": ["worm", "c. elegans", "caenorhabditis elegans", "nematode"],
        "Fish": ["fish", "d. rerio", "danio rerio", "zebrafish"],
    }
    ENRICHR_URL = "http://maayanlab.cloud"
    database = ""
    if organism.lower() in default:
        database = "Enrichr"
    else:
        for k, v in _organisms.items():
            if organism.lower() in v:
                database = k + "Enrichr"
                break

    if not database.endswith("Enrichr"):
        raise LookupError(
            """No supported database. Please input one of these:
                            ('Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm') """
        )
    # make a get request to get the gmt names and meta data from Enrichr
    # old code
    # response = requests.get('http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=meta')
    # gmt_data = response.json()
    # # generate list of lib names
    # libs = []
    # # get library names
    # for inst_gmt in gmt_data['libraries']:
    #     # only include active gmts
    #     if inst_gmt['isActive'] == True:
    #         libs.append(inst_gmt['libraryName'])
    lib_url = "%s/%s/datasetStatistics" % (ENRICHR_URL, database)
    response = requests.get(lib_url, verify=True)
    if not response.ok:
        raise Exception("Error getting the Enrichr libraries")
    libs_json = json.loads(response.text)
    libs = [lib["libraryName"] for lib in libs_json["statistics"]]

    return sorted(libs)


def download_library(name: str, organism: str = "human") -> Dict[str, List[str]]:
    """download enrichr libraries.

    :param str name: the enrichr library name. see `gseapy.get_library_name()`.
    :param str organism: Select one from { 'Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm' }
    :return dict: gene_sets of the enrichr library from selected organism


    """

    default = [
        "human",
        "mouse",
        "hs",
        "mm",
        "homo sapiens",
        "mus musculus",
        "h. sapiens",
        "m. musculus",
    ]
    _organisms = {
        "Fly": ["fly", "d. melanogaster", "drosophila melanogaster"],
        "Yeast": ["yeast", "s. cerevisiae", "saccharomyces cerevisiae"],
        "Worm": ["worm", "c. elegans", "caenorhabditis elegans", "nematode"],
        "Fish": ["fish", "d. rerio", "danio rerio", "zebrafish"],
    }
    ENRICHR_URL = "http://maayanlab.cloud"
    database = ""
    if organism.lower() in default:
        database = "Enrichr"
    else:
        for k, v in _organisms.items():
            if organism.lower() in v:
                database = k + "Enrichr"
                break

    if not database.endswith("Enrichr"):
        raise LookupError(
            """No supported database. Please input one of these:
                            ('Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm') """
        )

    tmpname = "%s.%s.gmt" % (database, name)
    tempath = os.path.join(DEFAULT_CACHE_PATH, tmpname)
    if os.path.isfile(tempath):
        logging.info("Library is already downloaded in: %s, use local file" % tempath)
        genesets_dict = {}
        with open(tempath) as genesets:
            for line in genesets:
                entries = line.strip().split("\t")
                key = entries[0]
                genesets_dict[key] = entries[2:]
        return genesets_dict
    # queery string
    ENRICHR_URL = ENRICHR_URL + "/%s/geneSetLibrary" % database
    query_string = "?mode=text&libraryName=%s"
    # get
    response = requests.get(
        ENRICHR_URL + query_string % name, timeout=None, stream=True
    )
    if not response.ok:
        raise Exception(
            "Error fetching gene set library, input name is correct for the organism you've set?."
        )
    # reformat to dict
    genesets_dict = {}
    # outname = os.path.join(DEFAULT_CACHE_PATH,  "enrichr.%s.gmt" % libname)
    for line in response.iter_lines(chunk_size=1024, decode_unicode="utf-8"):
        line = line.strip().split("\t")
        k = line[0]
        v = map(lambda x: x.split(",")[0], line[2:])
        v = list(filter(lambda x: True if len(x) else False, v))
        genesets_dict[k] = v

    return genesets_dict
