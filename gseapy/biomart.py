from collections.abc import Iterable
from io import StringIO

import pandas as pd
import requests
from bioservices import BioMart


class Biomart(BioMart):
    """query from BioMart"""

    def __init__(self, host="www.ensembl.org", verbose=False):
        """A wrapper of BioMart() from bioseverices.

        How to query validated dataset, attributes, filters.
        Example::
        >>> from gseapy.parser import Biomart
        >>> bm = Biomart(verbose=False, host="ensembl.org")
        >>> ## view validated marts
        >>> marts = bm.get_marts()
        >>> ## view validated dataset
        >>> datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
        >>> ## view validated attributes
        >>> attrs = bm.get_attributes(dataset='hsapiens_gene_ensembl')
        >>> ## view validated filters
        >>> filters = bm.get_filters(dataset='hsapiens_gene_ensembl')
        >>> ## query results
        >>> queries = ['ENSG00000125285','ENSG00000182968'] # a python list
        >>> results = bm.query(dataset='hsapiens_gene_ensembl',
                            attributes=['entrezgene_id', 'go_id'],
                            filters={'ensembl_gene_id': queries}
                            )
        """
        super(Biomart, self).__init__(host=host, verbose=verbose)
        hosts = ["www.ensembl.org", "asia.ensembl.org", "useast.ensembl.org"]
        # if host not work, select next
        i = 0
        while (self.host is None) and (i < 3):
            self.host = hosts[i]
            i += 1

        self.attributes_xml = []
        self.filters_xml = []
        self.dataset_xml = ""

        params = {
            "version": "1.0",
            "virtualSchemaName": "default",
            "formatter": "TSV",
            "header": 0,
            "uniqueRows": 0,
            "configVersion": "0.6",
            "completionStamp": "1",
        }

        self.header = (
            """http://ensembl.org/biomart/martservice?query="""
            """<?xml version="%(version)s" encoding="UTF-8"?>"""
            """<!DOCTYPE Query>"""
            """<Query virtualSchemaName="%(virtualSchemaName)s" formatter="%(formatter)s" """
            """header="%(header)s" uniqueRows="%(uniqueRows)s" count="" """
            """datasetConfigVersion="%(configVersion)s" completionStamp="%(completionStamp)s">"""
        )

        self.header = self.header % params
        self.footer = "</Dataset></Query>"
        self.reset()

    def add_filter(self, name, value):

        _filter = ""
        if "=" in value:
            _filter = """<Filter name="%s" %s/>""" % (name, value)
        else:
            _filter = """<Filter name="%s" value="%s"/>""" % (name, value)
        self.filters_xml.append(_filter)

    def add_attribute(self, attribute):
        _attr = """<Attribute name="%s"/>""" % attribute

        self.attributes_xml.append(_attr)

    def add_dataset(self, dataset):
        self.dataset_xml = """<Dataset name="%s" interface="default" >""" % dataset

    def reset(self):
        self.attributes_xml = []
        self.filters_xml = []
        self.dataset_xml = ""

    def get_xml(self):
        xml = self.header
        xml += self.dataset_xml
        for line in self.filters_xml:
            xml += line
        for line in self.attributes_xml:
            xml += line
        xml += self.footer
        return xml

    def get_marts(self):
        """Get available marts and their names.

        URL:
        /martservice/marts  -> return xml format
        /martservice/marts.json -> return json format
        """

        mart_names = pd.Series(self.names, name="Name")
        mart_descriptions = pd.Series(self.displayNames, name="Description")

        return pd.concat([mart_names, mart_descriptions], axis=1)

    def get_datasets(self, mart="ENSEMBL_MART_ENSEMBL"):
        """Get available datasets from mart you've selected
        URL Example:
        /martservice/datasets?config=snp_config
        """
        datasets = super(Biomart, self).datasets(mart, raw=True)
        return pd.read_csv(
            StringIO(datasets),
            header=None,
            usecols=[1, 2],
            names=["Name", "Description"],
            sep="\t",
        )

    def get_attributes(self, dataset):
        """Get available attritbutes from dataset you've selected
        URL Example:
        /martservice/attributes?datasets=btaurus_snp&config=snp_config
        """
        attributes = super(Biomart, self).attributes(dataset)
        attr_ = [(k, v[0]) for k, v in attributes.items()]
        return pd.DataFrame(attr_, columns=["Attribute", "Description"])

    def get_filters(self, dataset):
        """Get available filters from dataset you've selected

         URL Example:
        /martservice/filters?datasets=btaurus_snp&config=snp_config

        """
        filters = super(Biomart, self).filters(dataset)
        filt_ = [(k, v[0]) for k, v in filters.items()]
        return pd.DataFrame(filt_, columns=["Filter", "Description"])

    def query(
        self, dataset="hsapiens_gene_ensembl", attributes=[], filters={}, filename=None
    ):
        """mapping ids using BioMart.

        :param dataset: str, default: 'hsapiens_gene_ensembl'
        :param attributes: str, list, tuple
        :param filters: dict, {'filter name': list(filter value)}
        :param host: www.ensembl.org, asia.ensembl.org, useast.ensembl.org
        :return: a dataframe contains all attributes you selected.

        """
        if not attributes:
            attributes = [
                "ensembl_gene_id",
                "external_gene_name",
                "entrezgene_id",
                "go_id",
            ]

        self.new_query()
        # 'mmusculus_gene_ensembl'
        self.add_dataset_to_xml(dataset)
        for at in attributes:
            self.add_attribute_to_xml(at)
        # add filters
        if filters:
            for k, v in filters.items():
                if isinstance(v, str) or not isinstance(v, Iterable):
                    continue
                v = ",".join(list(v))
                self.add_filter_to_xml(k, v)

        xml_query = super(Biomart, self).get_xml()
        results = super(Biomart, self).query(xml_query)
        if str(results).startswith("Query ERROR"):
            print(results)
            return results

        df = pd.read_csv(
            StringIO(results), header=None, sep="\t", index_col=None, names=attributes
        )
        if "entrezgene_id" in df.columns:
            df["entrezgene_id"] = df["entrezgene_id"].astype(pd.Int32Dtype())

        self.results = df
        # save file to cache path.
        if filename is not None:
            # mkdirs(DEFAULT_CACHE_PATH)
            # filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
            df.to_csv(filename, sep="\t", index=False)

        return df

    def query_simple(
        self, dataset="hsapiens_gene_ensembl", attributes=[], filters={}, filename=None
    ):
        """
        This function is a simple version of BioMart REST API.
        sample parameter to query().

        However, you could get cross page of mapping. such as Mouse 2 human gene names

        **Note**: it will take a couple of minutes to get the results.
        A xml template for querying biomart. (see https://gist.github.com/keithshep/7776579)

        Example::
            >>> from gseapy import Biomart
            >>> bm = Biomart()
            >>> results = bm.query_simple(dataset='mmusculus_gene_ensembl',
                                          attributes=['ensembl_gene_id',
                                                      'external_gene_name',
                                                      'hsapiens_homolog_associated_gene_name',
                                                      'hsapiens_homolog_ensembl_gene'])
        """
        self.reset()
        self.add_dataset(dataset)
        for at in attributes:
            self.add_attribute(at)
        for n, v in filters.items():
            self.add_filter(n, v)
        _xml = self.get_xml()
        response = requests.get(_xml, stream=True)
        if response.ok:
            df = pd.read_table(StringIO(response.text), header=None, names=attributes)
            if filename is not None:
                df.to_csv(filename, sep="\t", index=False)
            self.results = df
            return df
        return response.text
