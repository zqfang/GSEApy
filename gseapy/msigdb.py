import re

import pandas as pd
import requests


class Msigdb:
    def __init__(self, dbver: str = "2023.1.Hs"):
        """
        dbver: MSIGDB version number. default: 2023.1.Hs
        """
        self.url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
        self._pattern = re.compile("(\w.+)\.(v\d.+)\.(entrez|symbols)\.gmt")
        self._db_version = self._get_db_version()
        self.categoires = self.list_category(dbver)

    def _get_db_version(self):
        resp = requests.get(self.url)
        if resp.ok:
            d = pd.read_html(resp.text)[0]
            # remove item : parent dictory and NA columns
            d = d.dropna(how="all").iloc[1:, 1:3].reset_index(drop=True)
            d.iloc[:, 0] = d.iloc[:, 0].str.rstrip("/")
            return d
        return None

    def get_gmt(
        self, category: str = "h.all", dbver: str = "2023.1.Hs", entrez: bool = False
    ):
        """
        :params category: choose one from .list_category()
        :params dbver: choose one from .list_dbver()

        An example of query url: "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/c2.cp.kegg.v2023.1.Hs.entrez.gmt"
        """
        identifier = "symbols"
        if entrez:
            identifier = "entrez"
        url = f"{self.url}/{dbver}/{category}.v{dbver}.{identifier}.gmt"
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
        return self._db_version

    def list_category(self, dbver: str = "2023.1.Hs"):
        """
        dbver: MSIGDB version number. default: 2023.1.Hs
        see a list of dbver, call .list_dbver()
        """
        d = self.list_gmt(dbver)
        if d is not None:
            categories = (
                d.iloc[:, 0]
                .apply(lambda s: self._pattern.match(s).groups()[0])
                .drop_duplicates()
            )
            return categories.to_list()
        return None

    def list_gmt(self, db: str):
        url = self.url + db
        resp = requests.get(url)
        if resp.ok:
            d = pd.read_html(resp.text)[0]
            # remove item : parent dictory and NA columns
            d = d.dropna(how="all").iloc[1:, 1:4]
            d = d[d.iloc[:, 0].str.match(self._pattern)]
            return d.reset_index(drop=True)
        return None
