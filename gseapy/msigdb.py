import re
from io import StringIO

import pandas as pd
import requests


class Msigdb:
    url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    _pattern = re.compile(r"(\w.+)\.(v\d.+)\.(entrez|symbols)\.gmt")

    def __init__(self, dbver: str = "2023.1.Hs"):
        """
        dbver: MSIGDB version number. default: 2023.1.Hs
        """
        self._db_version = self._get_db_version()
        if self._db_version is None:
            raise Exception("Failed to fetch available MSIGDB versions")
        self.categoires = self.list_category(dbver)

    @classmethod
    def _get_db_version(cls):
        """
        Get all available MSIGDB versions

        Return:
            A pd.DataFrame of all available MSIGDB versions.
            If failed to fetch, return None.
        """
        resp = requests.get(cls.url)
        if resp.ok:
            d = pd.read_html(StringIO(resp.text))[0]
            # remove item : parent dictory and NA columns
            d = d.dropna(how="all").iloc[1:, 1:3].reset_index(drop=True)
            d.iloc[:, 0] = d.iloc[:, 0].str.rstrip("/")
            return d
        return None

    @classmethod
    def get_gmt(
        cls, category: str = "h.all", dbver: str = "2023.1.Hs", entrez: bool = False
    ):
        """
        :params category: choose one from .list_category()
        :params dbver: choose one from .list_dbver()

        An example of query url: "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/c2.cp.kegg.v2023.1.Hs.entrez.gmt"
        """
        identifier = "symbols"
        if entrez:
            identifier = "entrez"
        url = f"{cls.url}/{dbver}/{category}.v{dbver}.{identifier}.gmt"
        resp = requests.get(url)
        if resp.ok:
            d = {}
            for line in resp.text.strip().split("\n"):
                row = line.split("\t")
                d[row[0]] = row[2:]
            return d
        return None

    def list_dbver(self):
        # self._db_version.columns = ["dbver", "date"]
        """
        Return a pd.DataFrame of all available MSIGDB versions

        Return:
            A pd.DataFrame of all available MSIGDB versions.
            If failed to fetch, return None.
        """
        return self._db_version

    @classmethod
    def list_category(cls, dbver: str = "2023.1.Hs"):
        """
        dbver: MSIGDB version number. default: 2023.1.Hs
        see a list of dbver, call .list_dbver()
        """
        d = cls.list_gmt(dbver)
        if d is not None:
            categories = (
                d.iloc[:, 0]
                .apply(lambda s: cls._pattern.match(s).groups()[0])
                .drop_duplicates()
            )
            return categories.to_list()
        return None

    @classmethod
    def list_gmt(cls, db: str):
        """
        list all gmt files in MSIGDB database.

        :param db: MSIGDB version number. default: 2023.1.Hs
        :return: a pandas DataFrame object
        """

        url = cls.url + db
        resp = requests.get(url)
        if resp.ok:
            d = pd.read_html(StringIO(resp.text))[0]
            # remove item : parent dictory and NA columns
            d = d.dropna(how="all").iloc[1:, 1:4]
            d = d[d.iloc[:, 0].str.match(cls._pattern)]
            return d.reset_index(drop=True)
        return None
