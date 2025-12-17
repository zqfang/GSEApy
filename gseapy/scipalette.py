import importlib
import json
from typing import List

from matplotlib.colors import LinearSegmentedColormap


class SciPalette:
    def __init__(self):
        """
        Cateogorical color palette collection of popular-sci journals.
        """
        handle = importlib.resources.open_binary("gseapy.data", "palette.json")
        self._db = json.load(handle)
        handle.close()

    def __repr__(self):
        # Name that Color. http://chir.ag/projects/name-that-color/
        return str(f"{self.__class__.__name__} color palette collection")

    def name_color(self, hex=None):
        """naming  a color

        :param hex: hex code of color
        """
        return "go to: http://chir.ag/projects/name-that-color/"

    @staticmethod
    def create_colormap(
        colors: List[str] = ["#000080", "#ffffff", "#8b0000"],
        positions: List[float] = [0.0, 0.5, 1.0],
        name: str = "navyblue2darkred",
    ):
        """create colormap given 3 color and break points

        :param colors: a list of color names. Default: ["#000080", "#ffffff", "#8b0000"]
        :param position: position of each color in range [0,1]. Default: [0.0, 0.5, 1.0]
        :param name: name of the return cmap
        :return: matplotlib cmap object
        """
        # cmap = ListedColormap(["#000080", "#ffffff", "#8b0000"])
        # positions = [0.0, 0.5, 1.0]
        if name is None:
            name = "mycmap"
        if colors is None:
            colors = ["#000080", "#ffffff", "#8b0000"]  # navyblue, white, darkred
        if positions is None:
            return LinearSegmentedColormap.from_list(name, colors)
        return LinearSegmentedColormap.from_list(name, list(zip(positions, colors)))

    @property
    def npg(self) -> List[str]:
        """Discrete Color Palettes inspired by plots in Nature Reviews Cancer"""
        return list(self._db["npg"].values())

    @property
    def aaas(self) -> List[str]:
        """Color palette inspired by plots in Science from AAAS"""
        return list(self._db["aaas"].values())

    @property
    def nejm(self) -> List[str]:
        """Color palette inspired by plots in The New England Journal of Medicine"""
        return list(self._db["nejm"].values())

    @property
    def lancet(self) -> List[str]:
        """Color palette inspired by plots in Lancet Oncology"""
        return list(self._db["lancet"].values())

    @property
    def jama(self) -> List[str]:
        """Color palette inspired by plots in The Journal of the American Medical Association"""
        return list(self._db["jama"].values())

    @property
    def jco(self) -> List[str]:
        """Color palette inspired by plots in Journal of Clinical Oncology"""
        return list(self._db["jco"].values())

    @property
    def ucscgb(self) -> List[str]:
        """Color palette inspired by UCSC Genome Browser Chromosome Colors"""
        return list(self._db["ucscgb"].values())

    def d3js(self, category: str = "c20a") -> List[str]:
        """
        choose category from (c10, c20a, c20b, c20c)
        """
        if category in self._db["d3js"]:
            return list(self._db["d3js"][category].values())
        return []

    @property
    def igv(self) -> List[str]:
        """Color palette inspired by IGV"""
        return list(self._db["igv"].values())

    @property
    def igv_alternating(self) -> List[str]:
        """Color palette inspired by IGV"""
        return list(self._db["igv_alternating"].values())

    @property
    def locuszoom(self) -> List[str]:
        """Color palette inspired by LocusZoom"""
        return list(self._db["locuszoom"].values())

    def uchicago(self, category: str = "default") -> List[str]:
        """
        Color palette inspired by University of Chicago Color Palette

        choose category from (light, dark, default)
        """
        if category in self._db["uchicago"]:
            return list(self._db["uchicago"][category].values())
        return []

    def hallmark(self, category: str = "dark") -> List[str]:
        """
        choose category from (dark, light)
        """
        if category in self._db["hallmark"]:
            return list(self._db["hallmark"][category].values())
        return []

    @property
    def cosmic(self) -> List[str]:
        """Color palette inspired by COSMIC Hallmarks of Cancer"""
        return list(self._db["cosmic"].values())

    @property
    def simpsons(self) -> List[str]:
        """Color palette inspired by The Simpsons"""
        return list(self._db["simpsons"].values())

    @property
    def futurama(self) -> List[str]:
        """Color palette inspired by Futurama"""
        return list(self._db["futurama"].values())

    @property
    def rickandmorty(self) -> List[str]:
        """Color palette inspired by Rick and Morty"""
        return list(self._db["rickandmorty"].values())

    @property
    def startrek(self) -> List[str]:
        """Color palette inspired by Star Trek"""
        return list(self._db["startrek"].values())

    @property
    def tron(self) -> List[str]:
        """Color palette inspired by Tron Legacy"""
        return list(self._db["tron"].values())

    @property
    def gsea(self) -> List[str]:
        """Color palette inspired by heatmaps generated by GSEA GenePattern"""
        return list(self._db["gsea"].values())

    def material(self, category: str = "indigo") -> List[str]:
        """
        choose category from
        (red, pink, purple, deeppurple, blue, lightblue, indigo,
        cyan, teal, lime, green, yellow, amber, organge, deeporange,
        brown, gray, bluegray)
        """
        if category in self._db["material"]:
            return list(self._db["material"][category].values())
        return []

    @property
    def zeileis(self) -> List[str]:
        """
        https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
        orig reference http://epub.wu.ac.at/1692/1/document.pdf
        """
        return list(self._db["zeileis"].values())

    @property
    def godsnot(self) -> List[str]:
        """
        take from http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
        or
        http://godsnotwheregodsnot.blogspot.com/2013/11/kmeans-color-quantization-seeding.html
        """
        return list(self._db["godsnot"].values())

    @property
    def boynton(self) -> List[str]:
        return list(self._db["boynton"].values())

    @property
    def kelly(self) -> List[str]:
        """
        The set of 22 colours of maximum contrast proposed by Kenneth Kelly in the work:
        http://www.iscc-archive.org/pdf/PC54_1724_001.pdf
        """
        return list(self._db["kelly"].values())

    @property
    def watlington(self) -> List[str]:
        return list(self._db["watlington"].values())

    @property
    def glasbey(self) -> List[str]:
        return list(self._db["glasbey"].values())
