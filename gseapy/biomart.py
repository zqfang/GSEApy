import logging
import os
from collections.abc import Iterable
from io import StringIO
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from xml.etree import cElementTree as ET

import pandas as pd
import requests

from gseapy.utils import log_init, mkdirs, retry


class Biomart:
    """query from BioMart"""

    CHUNK_SIZE = 300  # Default chunk size for batching queries
    MAX_WORKERS: Optional[int] = None  # Default to auto; set an int to override

    def __init__(self, host: str = "www.ensembl.org", verbose: bool = False):
        """simple API to BioMart services.

        How to query validated dataset, attributes, filters.
        Example::
        >>> from gseapy import Biomart
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
        self._id = str(id(self))
        self._logger = log_init(
            name="Biomart" + self._id,
            log_level=logging.INFO if verbose else logging.WARNING,
            filename=None,
        )
        self._set_host(host)

        self.attributes_xml = []
        self.filters_xml = []
        self.dataset_xml = ""

        params = {
            "host": self.host,
            "version": "1.0",
            "virtualSchemaName": "default",
            "formatter": "TSV",
            "header": 0,
            "uniqueRows": 1,
            "configVersion": "0.6",
            "completionStamp": 0,
        }

        self.header = (
            """https://%(host)s/biomart/martservice?query="""
            """<?xml version="%(version)s" encoding="UTF-8"?>"""
            """<!DOCTYPE Query>"""
            """<Query virtualSchemaName="%(virtualSchemaName)s" formatter="%(formatter)s" """
            """header="%(header)s" uniqueRows="%(uniqueRows)s" count="" """
            """datasetConfigVersion="%(configVersion)s" completionStamp="%(completionStamp)s">"""
        )

        self.header = self.header % params
        self.footer = "</Dataset></Query>"
        self.reset()

        # get supported marts
        self._marts = None

    def __del__(self):
        handlers = self._logger.handlers[:]
        for handler in handlers:
            handler.close()  # close file
            self._logger.removeHandler(handler)

    def _set_host(self, host: str):
        """set host"""

        hosts = ["useast.ensembl.org", "asia.ensembl.org"]
        hosts.insert(0, host)
        secure = "s"
        # if host not work, select next
        i = 0
        while i < len(hosts):
            url = "http{}://{}/biomart/martservice?type=registry".format(
                secure, hosts[i]
            )
            request = requests.get(url)
            # '<html>\n\n<head>\n  <title>Service unavailable</title>\n
            # "\n<MartRegistry>\n"
            if request.ok and request.text.startswith("\n<MartRegistry>\n"):
                self.host = hosts[i]
                self._marts = self._get_mart(request.text)
                break
            self._logger.warning(
                "host {} is not reachable, try {} ".format(
                    hosts[i], hosts[(i + 1) % len(hosts)]
                )
            )
            i += 1
        if i == len(hosts):
            self._logger.warning("hosts is not reachable. Please try again later.")

    def add_filter(self, name: str, value: Iterable[str]):
        """
        key: filter names
        value: Iterable[str]
        """
        if isinstance(value, str):
            pass  # keep string as is
        elif isinstance(value, Iterable):
            value = ",".join([str(v) for v in list(value)])
        _filter = ""
        if name.lower().startswith("with"):
            _filter = """<Filter name="%s" excluded="%s"/>""" % (name, value)
        else:
            _filter = """<Filter name="%s" value="%s"/>""" % (name, value)
        self.filters_xml.append(_filter)

    def add_attribute(self, attribute: str):
        _attr = """<Attribute name="%s"/>""" % attribute
        self.attributes_xml.append(_attr)

    def add_dataset(self, dataset: str):
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

    def get_xml_body(self):
        """Return only the XML body without the URL prefix.

        This is suitable for POST requests where the XML is sent in the body
        as the 'query' form field to the biomart endpoint.
        """
        # self.header starts with 'https://{host}/biomart/martservice?query=' then the XML header
        # Split once on 'query=' to strip the URL and keep the XML portion intact
        try:
            xml_header = self.header.split("query=", 1)[1]
        except Exception:
            # Fallback: if format changes in the future, best effort to keep original behavior
            xml_header = self.header
        xml = xml_header
        xml += self.dataset_xml
        for line in self.filters_xml:
            xml += line
        for line in self.attributes_xml:
            xml += line
        xml += self.footer
        return xml

    def _build_xml_strings(
        self,
        dataset: str,
        attributes: List[str],
        filters: Dict[str, Iterable[str]],
    ) -> Dict[str, str]:
        """Build XML strings without mutating instance state.

        Returns a dict with keys:
          - body: XML suitable for POST body (no URL prefix)
          - url: Full GET URL with query parameter
        """
        # Build dataset and tags locally to keep thread-safety
        dataset_xml = f'<Dataset name="{dataset}" interface="default" >'

        filter_tags: List[str] = []
        for name, value in filters.items():
            if isinstance(value, str):
                v = value
            elif isinstance(value, Iterable):
                v = ",".join([str(vv) for vv in list(value)])
            else:
                v = str(value)
            if name.lower().startswith("with"):
                filter_tags.append(f'<Filter name="{name}" excluded="{v}"/>')
            else:
                filter_tags.append(f'<Filter name="{name}" value="{v}"/>')

        attribute_tags = [f'<Attribute name="{a}"/>' for a in attributes]

        # Header body portion (post-'query=') derived from self.header
        try:
            xml_header = self.header.split("query=", 1)[1]
        except Exception:
            xml_header = self.header
        body = (
            xml_header
            + dataset_xml
            + "".join(filter_tags)
            + "".join(attribute_tags)
            + self.footer
        )
        url = f"https://{self.host}/biomart/martservice?query=" + body
        return {"body": body, "url": url}

    def _get_mart(self, text: str):
        """
        Parse the xml text and return a dataframe of supported marts.

        Parameters
        ----------
        text : str
            a xml text

        Returns
        -------
        marts : pd.DataFrame
            a dataframe of supported marts with columns:
                - Mart: the name of mart
                - Version: the version of mart
        """
        marts = [e.attrib for e in ET.XML(text)]
        marts = pd.DataFrame(marts)
        required_columns = ["database", "displayName", "name"]
        missing = [col for col in required_columns if col not in marts.columns]
        if missing:
            raise ValueError(
                f"BioMart registry XML missing columns: {missing}. Schema may have changed."
            )
        marts = marts.loc[:, required_columns]
        marts.columns = ["Version", "DisplayName", "Mart"]
        # get supported marts
        return marts.loc[:, ["Mart", "Version"]]

    def get_marts(self):
        """Get available marts and their names."""
        url = "https://{host}/biomart/martservice?type=registry&requestid=gseapy{i}".format(
            host=self.host, i=self._id
        )
        if self._marts is not None:
            return self._marts
        resp = requests.get(url)
        if resp.ok and resp.text.startswith("\n<MartRegistry>\n"):
            self._marts = self._get_mart(resp.text)
            return self._marts

        return resp.text

    def get_datasets(self, mart: str = "ENSEMBL_MART_ENSEMBL"):
        """Get available datasets from mart you've selected"""

        marts = self.get_marts()
        if mart not in marts["Mart"].values:
            raise ValueError(
                "Provided mart name (%s) is not valid. see 'names' attribute" % mart
            )

        url = "https://{host}/biomart/martservice?type=datasets&mart={mart}".format(
            host=self.host, mart=mart
        )
        resp = requests.get(url)
        if resp.ok:
            if resp.text.startswith("Problem"):
                return resp.text
            datasets = [
                record.split("\t")
                for record in resp.text.split("\n")
                if len(record) > 1
            ]
            datasets = pd.DataFrame(datasets).iloc[:, 1:3]
            datasets.columns = ["Dataset", "Description"]
            return datasets
        return resp.text

    def get_attributes(self, dataset: str = "hsapiens_gene_ensembl"):
        """Get available attritbutes from dataset you've selected"""
        # assert dataset in

        url = "https://{host}/biomart/martservice?type=attributes&dataset={dataset}".format(
            host=self.host, dataset=dataset
        )
        resp = requests.get(url)
        if resp.ok:
            attributes = [text.split("\t") for text in resp.text.strip().split("\n")]
            attributes = pd.DataFrame(attributes).iloc[:, :3]
            attributes.columns = ["Attribute", "Description", "Additional"]
            return attributes
        return resp.text

    def get_filters(self, dataset: str = "hsapiens_gene_ensembl"):
        """Get available filters from dataset you've selected"""
        # filters = super(Biomart, self).filters(dataset)
        # if dataset not in [x for k in self.valid_attributes.keys() for x in self.valid_attributes[k]]:
        #     raise ValueError("provided dataset (%s) is not found. see valid_attributes" % dataset)
        url = (
            "https://{host}/biomart/martservice?type=filters&dataset={dataset}".format(
                host=self.host, dataset=dataset
            )
        )

        resp = requests.get(url)
        if resp.ok:
            if str(resp.text).startswith("Query ERROR"):
                return resp.text
            filters = [text.split("\t") for text in resp.text.strip().split("\n")]
            df_filters = pd.DataFrame(filters)
            # Check if there are enough columns before selecting
            expected_indices = [0, 1, 3, 5]
            if df_filters.shape[1] >= max(expected_indices) + 1:
                df_filters = df_filters.iloc[:, expected_indices]
                df_filters.columns = [
                    "Filter",
                    "Description",
                    "Additional",
                    "InputType",
                ]
                return df_filters
            else:
                self._logger.warning(
                    f"Filter response has {df_filters.shape[1]} columns, expected at least {max(expected_indices)+1}. Returning raw DataFrame."
                )
                return df_filters
        return resp.text

    def query(
        self,
        dataset: str = "hsapiens_gene_ensembl",
        attributes: Optional[List[str]] = None,
        filters: Optional[Dict[str, Iterable[str]]] = None,
        filename: Optional[str] = None,
    ):
        """mapping ids using BioMart.

        :param dataset: str, default: 'hsapiens_gene_ensembl'
        :param attributes: str, list, tuple
        :param filters: dict, {'filter name': list(filter value)}
        :param host: www.ensembl.org, asia.ensembl.org, useast.ensembl.org
        :return: a dataframe contains all attributes you selected.

        Example:

            >>> queries = {'ensembl_gene_id': ['ENSG00000125285','ENSG00000182968'] } # need to be a python dict
            >>> results = bm.query(dataset='hsapiens_gene_ensembl',
                                   attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
                                   filters=queries)
        """
        if attributes is None or not attributes:
            attributes = [
                "ensembl_gene_id",
                "external_gene_name",
                "entrezgene_id",
                "go_id",
            ]
        if isinstance(attributes, str):
            attributes = attributes.split(",")

        if filters is None:
            filters = {}
        if not isinstance(filters, dict):
            raise ValueError("filters only accept a dict object")

        df = self.query_simple(
            dataset=dataset, filters=filters, attributes=attributes, filename=None
        )
        if df is None:
            return
        elif isinstance(df, str):
            print(df)
        if "entrezgene_id" in df.columns:
            if not pd.api.types.is_integer_dtype(df["entrezgene_id"]):
                df["entrezgene_id"] = df["entrezgene_id"].astype(pd.Int32Dtype())
        if "entrezgene_id" in df.columns:
            df["entrezgene_id"] = df["entrezgene_id"].astype(pd.Int32Dtype())

        self.results = df
        # save file to cache path.
        if filename is not None:
            mkdirs(os.path.dirname(filename))
            df.to_csv(filename, sep="\t", index=False)

        return df

    def query_simple(
        self,
        dataset: str = "hsapiens_gene_ensembl",
        attributes: Optional[List[str]] = None,
        filters: Optional[Dict[str, Iterable[str]]] = None,
        filename: Optional[str] = None,
    ):
        """
        This function is a simple version of BioMart REST API.
        same parameter to query().

        However, you could get cross page of mapping. such as Mouse 2 human gene names

        **Note**: it will take a couple of minutes to get the results.
        A xml template for querying biomart. (see https://gist.github.com/keithshep/7776579)

        Example::
            >>> from gseapy import Biomart
            >>> bm = Biomart()
            >>> results = bm.query_simple(dataset='mmusculus_gene_ensembl',
                                          attributes=['ensembl_gene_id',
                                                      'external_gene_name',
                                                      'hsapiens_homolog_associated_gene_name'])
        """
        if attributes is None:
            attributes = []
        if filters is None:
            filters = {}

        self.reset()
        self.add_dataset(dataset)
        self._logger.debug("Add attributes")
        for at in attributes:
            self.add_attribute(at)
        self._logger.debug("Add filters")
        for n, v in filters.items():
            self.add_filter(n, v)
        self._logger.debug("Build xml")
        # Build XML for both GET (URL) and POST (body)
        self._xml = self.get_xml()
        xml_body = self.get_xml_body()

        endpoint = f"https://{self.host}/biomart/martservice"
        s = retry(num=5)

        def _parse_response_text(text: str):
            if str(text).startswith("Query ERROR"):
                self._logger.error(text)
                return None
            df_local = pd.read_table(StringIO(text), header=None)
            # Validate column count before assigning names
            if len(df_local.columns) == len(attributes):
                df_local.columns = attributes
            else:
                self._logger.warning(
                    f"Response column count ({len(df_local.columns)}) does not match attributes ({len(attributes)})."
                )
            return df_local

        def _try_post(xml_body_local: str):
            try:
                resp = s.post(endpoint, data={"query": xml_body_local})
            except Exception as e:
                self._logger.warning(f"POST request failed: {e}")
                return None, None
            if resp.ok:
                df_local = _parse_response_text(resp.text)
                if df_local is not None:
                    return df_local, None
                return None, resp.text
            return None, resp.text

        def _try_get(xml_url_local: str):
            try:
                resp = s.get(xml_url_local)
            except Exception as e:
                self._logger.warning(f"GET request failed: {e}")
                return None, None
            if resp.ok:
                df_local = _parse_response_text(resp.text)
                if df_local is not None:
                    return df_local, None
                return None, resp.text
            return None, resp.text

        # Strategy:
        # 1) Prefer POST to avoid 414 (URL too long)
        # 2) Fallback to GET for small queries
        # 3) If still failing (e.g., payload too large), batch the largest list-like filter

        # First attempt: POST once with the whole query
        df = None
        err_text = None

        # Heuristic: if URL is short, it's safe to try GET directly too, but POST is fine for both
        df, err_text = _try_post(xml_body)
        if df is None and (err_text is None or "414" in str(err_text)):
            # Try GET if POST didn't return a valid table (network issue) or explicit 414 indicates GET may also fail
            df, err_text = _try_get(self._xml)

        if df is not None:
            if filename is not None:
                df.to_csv(filename, sep="\t", index=False)
            self.results = df
            return df

        # If we reached here, try batching on the largest list-like filter value
        # Identify list-like filters
        list_like_keys = []
        for k, v in filters.items():
            if isinstance(v, (list, tuple, set, pd.Series)):
                list_like_keys.append((k, list(v)))
        if not list_like_keys:
            # No list-like filters to batch; return the last error text
            return err_text if err_text is not None else "BioMart request failed."

        # Choose the largest filter to batch
        key_to_batch, values_to_batch = max(list_like_keys, key=lambda kv: len(kv[1]))
        if len(values_to_batch) == 0:
            return err_text if err_text is not None else "No values to query."

        # Reasonable chunk size that works reliably with BioMart
        dfs: List[pd.DataFrame] = []

        # Determine chunk size; allow runtime override via instance/class attribute if present
        chunk_size = getattr(self, "CHUNK_SIZE", 300)
        # Build chunk specs first (immutable)
        chunks = []
        for i in range(0, len(values_to_batch), chunk_size):
            sub_values = values_to_batch[i : i + chunk_size]
            chunk_filters = {**filters, key_to_batch: sub_values}
            xmls = self._build_xml_strings(dataset, attributes, chunk_filters)
            chunks.append((i // chunk_size + 1, xmls["body"], xmls["url"]))

        def _fetch_chunk(idx: int, body_xml: str, url_xml: str):
            df_chunk, err_text_chunk = _try_post(body_xml)
            if df_chunk is None:
                df_chunk, err_text_chunk = _try_get(url_xml)
            if df_chunk is None:
                self._logger.warning(f"BioMart chunk {idx} failed: {err_text_chunk}")
            return df_chunk

        # Parallelize chunk fetching with a modest concurrency to improve throughput
        # Keep concurrency bounded to be nice to the BioMart service
        max_workers = getattr(self, "MAX_WORKERS", None)
        if max_workers is None:
            # Enable moderate parallelism for very large inputs; otherwise keep sequential
            max_workers = 1 if len(values_to_batch) < 2000 else 8

        if max_workers > 1 and len(chunks) > 1:
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                futures = [
                    ex.submit(_fetch_chunk, idx, body, url) for idx, body, url in chunks
                ]
                for fut in as_completed(futures):
                    df_chunk = fut.result()
                    if df_chunk is not None:
                        dfs.append(df_chunk)
        else:
            for idx, body, url in chunks:
                df_chunk = _fetch_chunk(idx, body, url)
                if df_chunk is not None:
                    dfs.append(df_chunk)

        if not dfs:
            return (
                err_text
                if err_text is not None
                else "BioMart request failed in all chunks."
            )

        df_all = pd.concat(dfs, ignore_index=True).drop_duplicates()
        if filename is not None:
            df_all.to_csv(filename, sep="\t", index=False)
        self.results = df_all
        return df_all
