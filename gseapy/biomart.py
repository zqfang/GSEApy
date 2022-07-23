import pandas as pd
from io import StringIO
from collections.abc import Iterable
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

    def get_marts(self):
        """Get available marts and their names."""

        mart_names = pd.Series(self.names, name="Name")
        mart_descriptions = pd.Series(self.displayNames, name="Description")

        return pd.concat([mart_names, mart_descriptions], axis=1)

    def get_datasets(self, mart="ENSEMBL_MART_ENSEMBL"):
        """Get available datasets from mart you've selected"""
        datasets = self.datasets(mart, raw=True)
        return pd.read_csv(
            StringIO(datasets),
            header=None,
            usecols=[1, 2],
            names=["Name", "Description"],
            sep="\t",
        )

    def get_attributes(self, dataset):
        """Get available attritbutes from dataset you've selected"""
        attributes = self.attributes(dataset)
        attr_ = [(k, v[0]) for k, v in attributes.items()]
        return pd.DataFrame(attr_, columns=["Attribute", "Description"])

    def get_filters(self, dataset):
        """Get available filters from dataset you've selected"""
        filters = self.filters(dataset)
        filt_ = [(k, v[0]) for k, v in filters.items()]
        return pd.DataFrame(filt_, columns=["Filter", "Description"])

    def query(
        self, dataset="hsapiens_gene_ensembl", 
        attributes=[], 
        filters={}, 
        filename=None
    ):
        """mapping ids using BioMart.  

        :param dataset: str, default: 'hsapiens_gene_ensembl'
        :param attributes: str, list, tuple
        :param filters: dict, {'filter name': list(filter value)}
        :param host: www.ensembl.org, asia.ensembl.org, useast.ensembl.org
        :return: a dataframe contains all attributes you selected.

        **Note**: it will take a couple of minutes to get the results.
        A xml template for querying biomart. (see https://gist.github.com/keithshep/7776579)
        Example::
        >>> import requests
        >>> exampleTaxonomy = "mmusculus_gene_ensembl"
        >>> exampleGene = "ENSMUSG00000086981,ENSMUSG00000086982,ENSMUSG00000086983"
        >>> urlTemplate = \
            '''http://ensembl.org/biomart/martservice?query=''' \
            '''<?xml version="1.0" encoding="UTF-8"?>''' \
            '''<!DOCTYPE Query>''' \
            '''<Query virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">''' \
            '''<Dataset name="%s" interface="default"><Filter name="ensembl_gene_id" value="%s"/>''' \
            '''<Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/>''' \
            '''<Attribute name="transcript_start"/><Attribute name="transcript_end"/>''' \
            '''<Attribute name="exon_chrom_start"/><Attribute name="exon_chrom_end"/>''' \
            '''</Dataset>''' \
            '''</Query>'''    
        >>> exampleURL = urlTemplate % (exampleTaxonomy, exampleGene)
        >>> req = requests.get(exampleURL, stream=True)

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

        xml_query = self.get_xml()
        results = super(Biomart, self).query(xml_query)
        if str(results).startswith("Query ERROR"):
            print(results)
            return results

        df = pd.read_csv( StringIO(results), 
                          header=None, 
                          sep="\t", 
                          index_col=None, 
                          names=attributes)
        if "entrezgene_id" in df.columns:
            df["entrezgene_id"] = df["entrezgene_id"].astype(pd.Int32Dtype())

        self.results = df
        # save file to cache path.
        if filename is not None:
            # mkdirs(DEFAULT_CACHE_PATH)
            # filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
            df.to_csv(filename, sep="\t", index=False)

        return df
