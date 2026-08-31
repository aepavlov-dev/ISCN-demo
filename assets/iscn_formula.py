"""ISCN 2024 formula generator, reverse parser and ambiguity analyser.

Purpose
-------
Turn a copy-number / structural result from sequencing (PGT-A, PGT-SR) into
every ISCN-conformant written form, then verify each form by parsing it back
with an independent grammar and by measuring what the form cannot express.

Design
------
1.  One typed input (``Case``) is the single source of every form, so the
    forms cannot disagree with each other.
2.  ``render(case, form)`` never writes free text: every token is produced by
    a rule with a citation into ISCN 2024 (see ``RULES``).
3.  ``project(case, form)`` states, mechanically, which facts a form is able
    to carry.  ``parse(text)`` recovers facts from a string with a grammar
    that does not call the renderer.  Reverse validation is then an exact
    identity: ``parse(render(case, f)) == project(case, f)``.
4.  Ambiguity is measured, not asserted: a form is ambiguous when ``project``
    maps two different cases onto the same value (collision), and
    ``LOSS_MATRIX`` records which attribute each form drops.

Rule citations refer to ISCN 2024 (Karger, 2024) section numbers.

Coordinates
-----------
All coordinates in ``Event`` are 1-based inclusive nucleotide positions, the
convention ISCN shares with HGVS (11.2.3.1).  Callers holding 0-based
half-open (BED) intervals must use :func:`from_bed`.
"""

from __future__ import annotations

import argparse
import bisect
import dataclasses
import gzip
import json
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "Case", "Event", "Form", "FORMS", "RULES", "LOSS_MATRIX",
    "CytobandTable", "render", "render_all", "parse", "project",
    "roundtrip_check", "ambiguity_report", "validate_case", "from_bed",
]

TIMES = "\u00d7"          # multiplication sign, ISCN 4.4.3
ARROW = "\u2192"          # from-to, detailed system, 5.4.2.2
NDASH = "\u2013"          # 1–22 chromosome range in the abbreviated system

# ---------------------------------------------------------------------------
# Rule register: every rule the renderer or parser enforces, with its citation.
# ---------------------------------------------------------------------------

RULES: Dict[str, str] = {
    '8.2.2': 'Abnormal copy number results: the abnormality is described by the band designation, the nucleotides and the copy number.',
    '8.2.5': 'Mixed cell populations and uncertain copy number.',
    '4.4.1': 'Spaces are not given in the ISCN description except for three scenarios (see 4.4.1a-c).',
    '4.7.1': 'The karyotype format describes the abnormality using chromosome bands (short and detailed systems).',
    '5.3': 'Numerical abnormalities.',
    '5.4.2.2': 'Structurally altered chromosomes are defined by their band composition in the detailed system (karyotype format).',
    '11.2.2b': 'Where several abnormalities are described in the ISCN, a mixture of the short and detailed systems with the karyotype format may be used.',
    '5.1r': 'In the interest of clarity, complex rearrangements necessitating descriptions using the detailed system (karyotype format) should be written out in full the first time they are used in the report.',
    '4.4.5h': 'The karyotype and microarray formats may be combined in one description.',
    '4.5.1': 'A clone is a cell population of common origin, described as its own line.',
    '4.5.3d': 'Clones are listed in decreasing size; the normal clone is last.',
    '4.5.3h': 'mos precedes the first clone in mosaicism, chi in chimerism.',
    '4.6c': 'A change below the banding resolution is not included in the banded karyotype.',
    '5.3.1.1': 'Sex chromosome aneuploidy is expressed by the complement itself, not by plus and minus signs.',
    '5.4.1': 'The resolution of the record follows from the number of bands per haploid set.',
    '5.6': 'Several copies of a rearranged chromosome are denoted by the multiplication sign and the number of copies.',
    '5.7': 'The ploidy level is given in angle brackets when it differs from diploid.',
    '8.1.1b': 'A normal microarray result is arr, a space, the sex chromosomes, then the autosomes, the multiplication sign and the copy number.',
    '8.2.6': 'Zygosity is denoted hmz and htz; for a whole complement it follows the copy number.',
    '8.3b': 'For shallow sequencing the resolution of the method is about 5-10 Mb.',
    '11.2.2': 'The sequencing line is designated seq, or sseq for shallow sequencing.',
    '11.3': 'A normal sequencing result is written in the abbreviated system.',
    '11.4.1': 'Whole-chromosome aneuploidy from sequencing is written in the abbreviated system.',
    '11.4.2.1': 'A segmental change from sequencing is written in the short or extended system.',
    '11.4.2.3': 'If multiple copy number alterations are detected without any structural information, the nomenclature can be written similar to the microarray format nomenclature; examples describe copy number alterations associated with structural rearrangements.',
    "4.2.1a": "Chromosome count first, comma, sex chromosome complement, comma, "
              "then the abnormalities; no spaces around commas.",
    "4.2.1b": "If a sex chromosome is structurally abnormal, the normal sex "
              "chromosome (if present) is listed first.",
    "4.2.1c": "Metaphase counts are not given for constitutional samples unless "
              "there is clinically significant mosaicism.",
    "4.3b":   "Sex chromosome aberrations first (X before Y), then autosomes in "
              "numerical order, irrespective of aberration type.",
    "4.3c":   "For each chromosome, numerical abnormalities precede structural "
              "ones; constitutional precede the same acquired abnormality.",
    "4.3d":   "Multiple structural changes of homologous chromosomes are listed "
              "alphabetically by the abbreviated term.",
    "4.3e":   "Multiple changes within the same chromosome are listed pter to "
              "qter, overriding alphabetical order.",
    "4.4.1a": "A space separates the technique token (with genome build) from "
              "the result.",
    "4.4.1b": "A space separates two abbreviations when no period or comma "
              "intervenes.",
    "4.4.1c": "A space precedes and follows chi, con, mos, or, sep.",
    "4.4.3a": "The multiplication sign gives the number of copies and follows "
              "the abnormality.",
    "4.4.4d": "For a deletion or duplication with both breakpoints in one band, "
              "the band is repeated in karyotype format and given once in "
              "microarray format.",
    "4.4.5a": "Genomic coordinates follow the NCBI/UCSC translation table "
              "cytoBand.txt for the stated build.",
    "4.4.5b": "When coordinates are used, the genome build is given in square "
              "brackets immediately after the technique, then a space.",
    "4.4.5c": "The build is omitted when no coordinates are reported.",
    "4.4.5d": "Thousands separators are recommended, but are not used in the "
              "extended system because they confound the copy number.",
    "4.4.5e": "A nucleotide span is separated by an underscore.",
    "4.4.5f": "Only abnormal regions are named; sex chromosome regions first.",
    "4.5.2c": "mos is not used in microarray nomenclature.",
    "4.5.3e": "The normal cell line is listed last even when it is the largest.",
    "4.5.3n": "For molecular cytogenomic results the proportion of the sample "
              "is given in square brackets after the abnormality; a question "
              "mark stands for an undeterminable proportion; a copy-number "
              "range with a tilde is the alternative; mos is not used.",
    "4.7.2b": "Abbreviated system: no breakpoints and no nucleotides.",
    "4.7.2c": "Short system: bands, nucleotides and copy number.",
    "4.7.2d": "Extended system: the abnormal span plus the flanking normal "
              "nucleotides.",
    "5.2":    "46,U designates a normal complement with undisclosed sex "
              "chromosomes.",
    "5.5.11d": "i may be used for classical isochromosomes, e.g. i(X)(q10); "
               "otherwise der is used.",
    "5.5.15c": "Without an identified parental rearrangement the chromosome is "
               "a der, not a rec.",
    "5.7a":   "All changes are expressed relative to the ploidy level.",
    "8.1.1e": "A mixed cell population is given either as the estimated "
              "proportion in square brackets or as a copy-number range with a "
              "tilde.",
    "8.1.1": "Three systems exist in microarray format: abbreviated, short, "
              "extended.",
    "8.1.1m": "sseq replaces arr when shallow genome sequencing is used.",
    "8.2.3":  "Parental origin follows the copy number with no space "
              "(dn, mat, pat, inh, dmat, dpat, dinh, umat, upat).",
    "8.2.6a": "Zygosity of a region is hmz or htz.",
    "8.2.6c": "Uniparental disomy with known origin is umat or upat.",
    "8.3e":   "Polar body results use the chromatid term cht.",
    "11.2.2a": "The sequencing line starts with seq; the build follows in "
               "square brackets when coordinates are given; sseq is used for "
               "shallow sequencing.",
    "11.2.2c": "When copy number is detected without exact breakpoints or "
               "structural information (read-depth calling), the microarray "
               "format is used inside the seq line.",
    "11.2.2d": "A karyotype result is given first, then a period, then seq.",
    "11.2.2e": "Sex chromosomes are omitted unless a sex chromosome variant is "
               "reported or the complement is needed for interpretation.",
}

# ---------------------------------------------------------------------------
# Cytobands
# ---------------------------------------------------------------------------


class CytobandTable:
    """Cytoband coordinates for one genome build (ISCN 4.4.5a).

    The UCSC file is 0-based half-open; positions handed to and returned by
    this class are 1-based inclusive.
    """

    def __init__(self, path: str, build: str = "GRCh38") -> None:
        self.build = build
        self.path = path
        self.bands: Dict[str, List[Tuple[int, int, str]]] = {}
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if not line.strip():
                    continue
                chrom, start, end, name, _stain = line.rstrip("\n").split("\t")
                if "_" in chrom:                    # alt/random contigs
                    continue
                c = chrom[3:] if chrom.startswith("chr") else chrom
                self.bands.setdefault(c, []).append((int(start) + 1, int(end), name))
        for c in self.bands:
            self.bands[c].sort()
        self._starts = {c: [b[0] for b in v] for c, v in self.bands.items()}
        self.lengths = {c: v[-1][1] for c, v in self.bands.items()}

    # -- queries ----------------------------------------------------------
    def chromosomes(self) -> List[str]:
        return sorted(self.bands, key=chrom_sort_key)

    def band_at(self, chrom: str, pos: int) -> str:
        """Band name containing a 1-based position."""
        v = self.bands[chrom]
        i = bisect.bisect_right(self._starts[chrom], pos) - 1
        if i < 0 or pos > v[i][1]:
            raise ValueError(f"position {pos} outside {chrom} (len {self.lengths[chrom]})")
        return v[i][2]

    def band_span(self, chrom: str, start: int, end: int) -> Tuple[str, str]:
        return self.band_at(chrom, start), self.band_at(chrom, end)

    def band_bounds(self, chrom: str, band: str) -> Tuple[int, int]:
        for s, e, name in self.bands[chrom]:
            if name == band:
                return s, e
        raise ValueError(f"unknown band {chrom}{band}")

    def arm_bounds(self, chrom: str, arm: str) -> Tuple[int, int]:
        sel = [(s, e) for s, e, n in self.bands[chrom] if n.startswith(arm)]
        if not sel:
            raise ValueError(f"{chrom} has no {arm} arm in the table")
        return sel[0][0], sel[-1][1]

    def arm_of(self, chrom: str, pos: int) -> str:
        return self.band_at(chrom, pos)[0]

    def whole_chromosome(self, chrom: str) -> Tuple[int, int]:
        return 1, self.lengths[chrom]

    def covers_whole(self, chrom: str, start: int, end: int, tol: int = 3_000_000) -> bool:
        return start <= tol and end >= self.lengths[chrom] - tol

    def covers_arm(self, chrom: str, start: int, end: int,
                   tol: int = 3_000_000) -> Optional[str]:
        for arm in ("p", "q"):
            try:
                a, b = self.arm_bounds(chrom, arm)
            except ValueError:
                continue
            if abs(start - a) <= tol and abs(end - b) <= tol:
                return arm
        return None


def chrom_sort_key(chrom: str) -> Tuple[int, int]:
    """ISCN 4.3b order: X, then Y, then autosomes ascending."""
    if chrom == "X":
        return (0, 0)
    if chrom == "Y":
        return (0, 1)
    return (1, int(chrom))


# ---------------------------------------------------------------------------
# Event model
# ---------------------------------------------------------------------------

CN_TYPES = {"gain", "loss"}
STRUCT_SYMBOLS = {
    "del", "dup", "trp", "qdp", "inv", "ins", "i", "idic", "ider", "dic",
    "der", "rec", "rob", "r", "t", "add", "fis", "trc", "mar", "ace",
}
INHERITANCE = {"dn", "mat", "pat", "inh", "dmat", "dpat", "dinh", "umat", "upat",
               "mat pat"}
ZYGOSITY = {"hmz", "htz", "hmz htz"}


@dataclass
class Event:
    """One reported abnormality.

    chrom / start / end / copy_number describe the dosage fact (what read
    depth measures).  ``structure`` describes the mechanism when an
    independent channel established it; it drives the karyotype-format forms
    and is refused when absent (5.5.15c).
    """

    chrom: str
    copy_number: int
    start: Optional[int] = None            # 1-based inclusive
    end: Optional[int] = None
    scale: str = "segment"                 # "chromosome" | "arm" | "segment"
    arm: Optional[str] = None
    flank_left: Optional[int] = None       # extended system (4.7.2d)
    flank_right: Optional[int] = None
    flank_copy_number: Optional[int] = None
    copy_number_range: Optional[Tuple[int, int]] = None   # 8.1.1e tilde form
    mosaic_fraction: Optional[float] = None               # 4.5.3n bracket form
    mosaic_fraction_range: Optional[Tuple[float, float]] = None
    mosaic_fraction_unknown: bool = False                 # renders [?]
    inheritance: Optional[str] = None
    zygosity: Optional[str] = None
    structure: Optional[Dict[str, Any]] = None
    supernumerary: bool = False            # extra chromosome, changes the count
    explained_by: Optional[str] = None     # dosage accounted for by another
                                           # event's structure: not listed twice
    label: str = ""                        # free text for tables, never rendered

    # -- derived ----------------------------------------------------------
    @property
    def direction(self) -> str:
        cn = max(self.copy_number_range) if self.copy_number_range else self.copy_number
        return "gain" if cn > self.baseline else "loss"

    baseline: int = 2                      # copies expected for this locus

    @property
    def is_numerical(self) -> bool:
        return self.scale == "chromosome"

    @property
    def count_delta(self) -> int:
        """Effect on the karyotype chromosome count (4.2.1a, 5.7a)."""
        if self.scale == "chromosome":
            return self.copy_number - self.baseline
        return 1 if self.supernumerary else 0


@dataclass
class Clone:
    """One cell line of a mosaic or chimeric result (4.5.1, 4.5.3).

    Clone-level input is needed whenever more than one line is abnormal:
    45,X/46,X,i(X)(q10) cannot be derived from a single dosage profile,
    because the dosage profile is the weighted sum of the two lines.
    """

    events: List[Event] = field(default_factory=list)
    sex_chromosomes: Optional[str] = None
    cells: Optional[int] = None          # metaphases, karyotype format only
    fraction: Optional[float] = None     # estimated share of the sample
    ploidy: int = 2
    normal: bool = False


@dataclass
class Case:
    """A complete result for one sample."""

    sample: str
    events: List[Event] = field(default_factory=list)
    clones: List[Clone] = field(default_factory=list)
    assembly: str = "GRCh38"
    technique: str = "seq"                 # seq | sseq | arr | ogm
    ploidy: int = 2                        # 5.7a
    sex_chromosomes: Optional[str] = None  # "XX", "XY", "X", "XXY", "U", None
    sex_disclosed: bool = True
    karyotype_available: bool = False      # was banded analysis performed?
    cell_counts: Optional[Dict[str, int]] = None   # clone -> metaphases
    chromatid: bool = False                # polar body, 8.3e
    genome_zygosity: Optional[str] = None  # 8.2.6: hmz / "hmz htz" (moles)
    notes: str = ""

    def sex_event_present(self) -> bool:
        return any(e.chrom in ("X", "Y") for e in self.events)

    def normal(self) -> bool:
        return not self.events and not any(c.events for c in self.clones)


def from_bed(chrom: str, start0: int, end: int, **kw) -> Event:
    """Build an Event from a 0-based half-open interval."""
    return Event(chrom=chrom, start=start0 + 1, end=end, **kw)


# ---------------------------------------------------------------------------
# Forms
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Form:
    key: str
    name_ru: str
    fmt: str          # "microarray" | "karyotype"
    system: str       # "abbreviated" | "short" | "extended" | "detailed"
    coords: bool
    sex_policy: str   # "auto" | "always" | "never" | "required"
    needs_structure: bool = False
    citation: str = ""


FORMS: Dict[str, Form] = {
    "abbrev": Form("abbrev", "сокращённая система, микрочиповый формат",
                   "microarray", "abbreviated", False, "auto",
                   citation="4.7.2b, 11.4.1"),
    "abbrev_sex": Form("abbrev_sex", "сокращённая система с явным половым набором",
                       "microarray", "abbreviated", False, "always",
                       citation="4.7.2b, 8.2.2, 11.2.2e"),
    "short": Form("short", "короткая система с координатами",
                  "microarray", "short", True, "auto",
                  citation="4.7.2c, 11.4.2.1"),
    "extended": Form("extended", "расширенная система с фланками",
                     "microarray", "extended", True, "auto",
                     citation="4.7.2d, 11.4.2.1"),
    "abbrev_range": Form("abbrev_range", "сокращённая система, мозаицизм диапазоном копийности",
                         "microarray", "abbreviated", False, "auto",
                         citation="4.7.2b, 8.1.1e, 8.2.5"),
    "short_range": Form("short_range", "короткая система, мозаицизм диапазоном копийности",
                        "microarray", "short", True, "auto",
                        citation="4.7.2c, 8.1.1e, 8.2.5"),
    "karyo_short": Form("karyo_short", "короткая система, кариотипный формат (цитобенды без координат)",
                        "karyotype", "short", False, "auto", True,
                        citation="4.7.1, 11.4.2.1"),
    "karyo_detailed": Form("karyo_detailed", "подробная система, кариотипный формат",
                           "karyotype", "detailed", False, "auto", True,
                           citation="4.7.1, 5.4.2.2, 11.4.2.3"),
    "mixed": Form("mixed", "смешанная запись: механизм кариотипным форматом, "
                           "остальные события микрочиповым",
                  "mixed", "short", True, "auto", False,
                  citation="4.4.5h, 11.2.2b"),
    "karyotype": Form("karyotype", "классическая кариотипная строка",
                      "karyotype", "short", False, "required", True,
                      citation="4.2.1a, 5.2, 5.3"),
    "karyotype_mos": Form("karyotype_mos", "кариотипная строка с клонами (mos)",
                          "karyotype", "short", False, "required", True,
                          citation="4.5.3e, 4.5.3h"),
}

# Which attribute of the input each form can carry.  This is the executable
# form of the loss matrix: project() below is generated from it.
ATTRIBUTES = ["chromosome", "direction", "copy_number", "breakpoint_bands",
              "coordinates", "flanking_normal", "mosaic_fraction",
              "mechanism", "parental_origin", "zygosity", "sex_complement",
              "segment_orientation"]

LOSS_MATRIX: Dict[str, Dict[str, bool]] = {
    "abbrev":        dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=False, coordinates=False,
                          flanking_normal=False, mosaic_fraction=True,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "abbrev_sex":    dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=False, coordinates=False,
                          flanking_normal=False, mosaic_fraction=True,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=True, segment_orientation=False),
    "short":         dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=True, coordinates=True,
                          flanking_normal=False, mosaic_fraction=True,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "extended":      dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=True, coordinates=True,
                          flanking_normal=True, mosaic_fraction=True,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "abbrev_range":  dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=False, coordinates=False,
                          flanking_normal=False, mosaic_fraction=False,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "short_range":   dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=True, coordinates=True,
                          flanking_normal=False, mosaic_fraction=False,
                          mechanism=False, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "karyo_short":   dict(chromosome=True, direction=True, copy_number=False,
                          breakpoint_bands=True, coordinates=False,
                          flanking_normal=False, mosaic_fraction=False,
                          mechanism=True, parental_origin=True, zygosity=False,
                          sex_complement=False, segment_orientation=False),
    "karyo_detailed": dict(chromosome=True, direction=True, copy_number=False,
                           breakpoint_bands=True, coordinates=False,
                           flanking_normal=False, mosaic_fraction=False,
                           mechanism=True, parental_origin=True, zygosity=False,
                           sex_complement=False, segment_orientation=True),
    "mixed":         dict(chromosome=True, direction=True, copy_number=True,
                          breakpoint_bands=True, coordinates=True,
                          flanking_normal=False, mosaic_fraction=True,
                          mechanism=True, parental_origin=True, zygosity=True,
                          sex_complement=False, segment_orientation=False),
    "karyotype":     dict(chromosome=True, direction=True, copy_number=False,
                          breakpoint_bands=True, coordinates=False,
                          flanking_normal=False, mosaic_fraction=False,
                          mechanism=True, parental_origin=True, zygosity=False,
                          sex_complement=True, segment_orientation=False),
    "karyotype_mos": dict(chromosome=True, direction=True, copy_number=False,
                          breakpoint_bands=True, coordinates=False,
                          flanking_normal=False, mosaic_fraction=True,
                          mechanism=True, parental_origin=True, zygosity=False,
                          sex_complement=True, segment_orientation=False),
}


# ---------------------------------------------------------------------------
# Input validation (before any rendering)
# ---------------------------------------------------------------------------


def validate_case(case: Case, bands: CytobandTable) -> List[str]:
    """Structural problems that make a formula impossible or wrong."""
    problems: List[str] = []
    for i, e in enumerate(case.events):
        tag = f"событие {i + 1} ({e.chrom})"
        if e.chrom not in bands.bands:
            problems.append(f"{tag}: хромосома отсутствует в таблице цитобендов")
            continue
        if e.scale != "chromosome" or e.start is not None:
            if e.start is None or e.end is None:
                problems.append(f"{tag}: масштаб {e.scale} требует координат")
            elif e.start < 1 or e.end > bands.lengths[e.chrom] or e.start > e.end:
                problems.append(
                    f"{tag}: координаты {e.start}_{e.end} вне 1..{bands.lengths[e.chrom]}")
        if e.copy_number < 0:
            problems.append(f"{tag}: отрицательная копийность")
        if (e.copy_number == e.baseline and not e.zygosity
                and not e.copy_number_range and not e.structure):
            problems.append(f"{tag}: копийность равна ожидаемой и зиготность не задана — "
                            "заявлять нечего")
        if e.mosaic_fraction is not None and not 0 < e.mosaic_fraction <= 1:
            problems.append(f"{tag}: доля мозаицизма {e.mosaic_fraction} вне (0;1]")
        if e.inheritance and e.inheritance not in INHERITANCE:
            problems.append(f"{tag}: неизвестное обозначение происхождения {e.inheritance!r}")
        if e.zygosity and e.zygosity not in ZYGOSITY:
            problems.append(f"{tag}: неизвестное обозначение зиготности {e.zygosity!r}")
        if e.structure:
            chain = e.structure.get("chain") or [e.structure]
            for term in chain:
                sym = term.get("symbol")
                if sym not in STRUCT_SYMBOLS:
                    problems.append(f"{tag}: неизвестный символ перестройки {sym!r}")
    if case.technique not in ("seq", "sseq", "arr", "ogm"):
        problems.append(f"неизвестная методика {case.technique!r}")
    if case.sex_chromosomes and not re.fullmatch(r"X+Y*|Y|U", case.sex_chromosomes):
        problems.append(f"недопустимый половой набор {case.sex_chromosomes!r}")
    return problems


# ---------------------------------------------------------------------------
# Ordering
# ---------------------------------------------------------------------------


def order_key_microarray(e: Event) -> Tuple:
    """4.3b + 4.4.5f + 8.2.2xxiii: sex chromosomes first, autosomes ascending,
    then pter to qter within a chromosome irrespective of gain or loss."""
    return chrom_sort_key(e.chrom) + (e.start or 0, e.end or 0)


def order_key_karyotype(e: Event) -> Tuple:
    """4.3b/c/d/e: chromosome order; within a chromosome numerical before
    structural, gains before losses, then pter to qter, then alphabetical."""
    numerical = 0 if e.is_numerical else 1
    gain = 0 if e.direction == "gain" else 1
    sym = (e.structure or {}).get("symbol", "")
    return chrom_sort_key(e.chrom) + (numerical, gain, e.start or 0, sym)


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------


def _thousands(n: int) -> str:
    return f"{n:,}"


def _copy_spec(e: Event, ascii_times: bool = False) -> str:
    t = "x" if ascii_times else TIMES
    if e.copy_number_range:
        lo, hi = e.copy_number_range
        return f"{t}{lo}~{hi}"
    return f"{t}{e.copy_number}"


def _suffixes(e: Event) -> str:
    """Zygosity, then parental origin, then sample proportion (8.2.3, 4.5.3n)."""
    out = ""
    if e.zygosity:
        out += e.zygosity                     # no space: ×2hmz (8.2.6)
    if e.inheritance:
        out += (" " if out else "") + e.inheritance   # 4.4.1b
    if e.mosaic_fraction_unknown:
        out += "[?]"
    elif e.mosaic_fraction_range:
        lo, hi = e.mosaic_fraction_range
        out += f"[{lo:g}~{hi:g}]"
    elif e.mosaic_fraction is not None:
        out += f"[{e.mosaic_fraction:g}]"
    return out


def _band_designation(e: Event, bands: CytobandTable, microarray: bool) -> str:
    """Band or band range for a segment (4.4.4d)."""
    b1, b2 = bands.band_span(e.chrom, e.start, e.end)
    if b1 == b2:
        return f"{e.chrom}{b1}" if microarray else f"{b1}{b2}"
    return f"{e.chrom}{b1}{b2}" if microarray else f"{b1}{b2}"


def _tech_token(case: Case, coords: bool) -> str:
    """Technique, build and the separating space (4.4.1a, 4.4.5b, 4.4.5c)."""
    tok = case.technique
    if coords:
        tok += f"[{case.assembly}]"
    return tok + " "


def _sex_prefix_needed(case: Case, form: Form) -> bool:
    if form.sex_policy == "never":
        return False
    if form.sex_policy == "always":
        return True
    # "auto" — 11.2.2e / 8.1.1c: only when a sex chromosome is involved or the
    # complement is needed for interpretation.
    return case.sex_event_present() or not case.sex_disclosed is True and False


def _sex_group(case: Case) -> str:
    """Sex complement as an abbreviated-system group, e.g. (X)x2,(Y)x1."""
    sx = case.sex_chromosomes or ""
    if sx in ("", "U"):
        return ""
    nx, ny = sx.count("X"), sx.count("Y")
    parts = []
    if nx:
        parts.append(f"(X){TIMES}{nx}")
    if ny:
        parts.append(f"(Y){TIMES}{ny}")
    return ",".join(parts)


# ---------------------------------------------------------------------------
# Microarray-format renderers
# ---------------------------------------------------------------------------


def _render_abbrev(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.7.2b: chromosome or arm level, no nucleotides, no build."""
    items: List[str] = []
    if form.sex_policy == "always" and case.sex_chromosomes not in (None, "U"):
        sg = _sex_group(case)
        if sg:
            items.append(sg)
    groups: Dict[Tuple, List[str]] = {}
    for e in sorted(case.events, key=order_key_microarray):
        if e.scale == "chromosome":
            token = e.chrom
        elif e.scale == "arm" or (e.start and bands.covers_arm(e.chrom, e.start, e.end)):
            arm = e.arm or bands.covers_arm(e.chrom, e.start, e.end)
            token = f"{e.chrom}{arm}"
        else:
            raise Unrenderable("abbrev", "сегментное событие не выражается "
                                         "сокращённой системой (4.7.2b)")
        key = (_copy_spec(e), _suffixes(e))
        groups.setdefault(key, []).append(token)
    for (cs, sfx), toks in groups.items():
        items.append(f"({','.join(toks)}){cs}{sfx}")
    cht = "cht" if case.chromatid else ""
    return _tech_token(case, coords=False) + cht + ",".join(items)


def _render_short(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.7.2c: bands, nucleotides with thousands separators, copy number."""
    items: List[str] = []
    if form.sex_policy == "always":
        sg = _sex_group(case)
        if sg:
            items.append(sg)
    for e in sorted(case.events, key=order_key_microarray):
        if e.start is None:
            # 8.2.2xxi: a whole-chromosome or whole-arm change keeps the
            # abbreviated group form inside a short-system line
            if e.scale == "chromosome":
                tokname = e.chrom
            elif e.scale == "arm" and e.arm:
                tokname = f"{e.chrom}{e.arm}"
            else:
                raise Unrenderable(
                    "short", "сегментное событие без координат: короткая система "
                             "требует фактических крайних опрошенных позиций "
                             "(8.2.2, 11.4.2.1)")
            items.append(f"({tokname}){_copy_spec(e)}{_suffixes(e)}")
            continue
        band = _band_designation(e, bands, microarray=True)
        items.append(f"{band}({_thousands(e.start)}_{_thousands(e.end)})"
                     f"{_copy_spec(e)}{_suffixes(e)}")
    coords = any(e.start is not None for e in case.events)
    return _tech_token(case, coords=coords) + ",".join(items)


def _render_extended(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.7.2d: abnormal span plus flanking normal positions, no thousands
    separators (4.4.5d)."""
    items: List[str] = []
    for e in sorted(case.events, key=order_key_microarray):
        if e.start is None or e.end is None:
            raise Unrenderable("extended", "расширенная система требует координат")
        if e.flank_left is None or e.flank_right is None:
            raise Unrenderable("extended", "расширенная система требует фланкирующих "
                                           "нормальных позиций (4.7.2d)")
        band = _band_designation(e, bands, microarray=True)
        fcn = e.flank_copy_number if e.flank_copy_number is not None else e.baseline
        inner = (f"{e.flank_left}{TIMES}{fcn},"
                 f"{e.start}_{e.end}{_copy_spec(e)},"
                 f"{e.flank_right}{TIMES}{fcn}")
        items.append(f"{band}({inner}){_suffixes(e)}")
    return _tech_token(case, coords=True) + ",".join(items)


# ---------------------------------------------------------------------------
# Karyotype-format renderers
# ---------------------------------------------------------------------------


class Unrenderable(Exception):
    """The form cannot express this case; the reason is part of the result."""

    def __init__(self, form: str, reason: str) -> None:
        super().__init__(f"{form}: {reason}")
        self.form = form
        self.reason = reason


def _struct_chain(e: Event, bands: CytobandTable) -> List[Dict[str, Any]]:
    """Normalize a structure record to a chain of rearrangement terms.

    A chain is how ISCN writes several rearrangements belonging to ONE
    derivative chromosome: the terms are juxtaposed with no comma between
    them (4.3e and the der() examples in 5.5.3).
    """
    st = e.structure
    if not st:
        raise Unrenderable("karyotype", "механизм не установлен: по правилу 5.5.15c "
                                        "запись структурной перестройки недопустима")
    if "chain" in st:
        chain = [dict(x) for x in st["chain"]]
    else:
        chain = [{k: v for k, v in st.items() if k in ("symbol", "chroms",
                                                       "breakpoints", "detailed")}]
    # Breakpoints are derived from coordinates only for symbols whose
    # breakpoints ARE the segment boundaries; der/t/ins/rec/rob describe a
    # junction and must be given explicitly.
    DERIVABLE = {"del", "dup", "inv", "i", "idic", "add", "r", "trp", "qdp"}
    NO_BREAKPOINTS = {"der", "mar", "rec", "ider"}
    for term in chain:
        term.setdefault("chroms", [e.chrom])
        if term.get("detailed") or term.get("breakpoints") is not None:
            continue
        if term["symbol"] in NO_BREAKPOINTS:
            term["breakpoints"] = []
            continue
        if term["symbol"] not in DERIVABLE:
            raise Unrenderable(
                "karyotype",
                f"для символа {term['symbol']} точки разрыва должны быть заданы "
                "явно: из границ сегмента они не выводятся")
        if e.start is None:
            raise Unrenderable("karyotype", "не заданы ни точки разрыва, ни "
                                            "координаты для их вывода")
        b1, b2 = bands.band_span(e.chrom, e.start, e.end)
        term["breakpoints"] = ([b1] if (b1 == b2 and term["symbol"] in
                                        ("i", "idic", "add")) else [b1, b2])
    return chain


def _bp_groups(term: Dict[str, Any]) -> Tuple[Tuple[str, ...], ...]:
    """Breakpoints grouped by chromosome.

    ISCN writes breakpoints of one chromosome juxtaposed and separates
    chromosomes with a semicolon: dup(1)(q25q42) is one chromosome with two
    breakpoints, r(1;3)(p36.1q23;q21q27) is two chromosomes with two each.
    A flat list is taken as one group per chromosome when the counts match.
    """
    bps = term.get("breakpoints")
    if not bps:
        return ()
    if isinstance(bps[0], (list, tuple)):
        return tuple(tuple(g) for g in bps)
    chroms = list(term.get("chroms") or [])
    if len(chroms) > 1 and len(bps) == len(chroms):
        return tuple((b,) for b in bps)
    return (tuple(bps),)


def _term_text(term: Dict[str, Any], detailed: bool = False) -> str:
    sym = term["symbol"]
    if detailed and term.get("detailed"):
        return f"{sym}({term['detailed']})"
    inner = ";".join(list(term.get("chroms") or []))
    groups = _bp_groups(term)
    if not groups:
        return f"{sym}({inner})"
    bptxt = ";".join("".join(g) for g in groups)
    return f"{sym}({inner})({bptxt})"


def _struct_token(e: Event, bands: CytobandTable, detailed: bool = False) -> str:
    """Structural abnormality in karyotype format (5.4.2.1, 5.4.2.2).

    Terms of one derivative chromosome are juxtaposed without a comma.
    """
    chain = _struct_chain(e, bands)
    txt = "".join(_term_text(t, detailed) for t in chain)
    return txt + _struct_tail(e)


def _struct_tail(e: Event) -> str:
    """Constitutional marker and parental origin written after the term
    (4.2.1e, 4.2.1g, 8.2.3)."""
    return (e.structure or {}).get("tail") or e.inheritance or ""


KARYOTYPE_RESOLUTION = 5_000_000
"""Smallest change a banded karyotype can carry.  ISCN 4.6c: a cryptic
abnormality that cannot be visualized is not listed in the banded karyotype;
5.4.1 ties resolution to bands per haploid set."""


def _karyotype_events(case: Case, banded: bool = True
                      ) -> Tuple[List[Event], List[Tuple[Event, str]]]:
    """Split events into those the line may carry and those it must leave out.

    ``banded=True`` is a microscope result and therefore also drops changes
    below the banding resolution (4.6c).  A molecular line (seq in karyotype
    format, or the mixed form) has no such limit and drops only dosage that
    another event's structure already accounts for.
    """
    keep: List[Event] = []
    omitted: List[Tuple[Event, str]] = []
    for e in case.events:
        if e.explained_by:
            omitted.append((e, f"доза объясняется перестройкой на chr{e.explained_by} "
                               "и отдельной записи не получает"))
        elif (banded and e.scale != "chromosome" and not e.structure
              and e.start is not None
              and (e.end - e.start + 1) < KARYOTYPE_RESOLUTION):
            omitted.append((e, f"{(e.end - e.start + 1) / 1e6:.2f} Мб — ниже предела "
                               "разрешения бэндирования, в кариотипную строку не "
                               "включается (4.6c)"))
        else:
            keep.append(e)
    return keep, omitted


def karyotype_omissions(case: Case) -> List[str]:
    return [f"chr{e.chrom}: {why}" for e, why in _karyotype_events(case)[1]]


def _require_mechanism(events: Sequence[Event]) -> None:
    missing = [e.chrom for e in events if e.scale != "chromosome" and not e.structure]
    if missing:
        raise Unrenderable(
            "karyotype",
            "механизм не установлен для событий на " +
            ", ".join(f"chr{c}" for c in sorted(set(missing), key=chrom_sort_key)) +
            ": по правилу 5.5.15c запись структурной перестройки недопустима")


def _render_mixed(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.4.5h: karyotype and microarray systems may be combined in one
    description, provided a single abnormality is not described by both."""
    keep, _ = _karyotype_events(case, banded=False)
    if not any(e.structure for e in keep):
        raise Unrenderable("mixed", "ни одно событие не имеет установленного "
                                    "механизма: смешанная запись не даёт выигрыша")
    items: List[str] = []
    for e in sorted(keep, key=order_key_microarray):
        if e.structure:
            tok = _struct_token(e, bands)
            if e.supernumerary:
                tok = "+" + tok
            items.append(tok)
        else:
            if e.start is None:
                raise Unrenderable("mixed", "событие без координат и без механизма")
            band = _band_designation(e, bands, microarray=True)
            items.append(f"{band}({_thousands(e.start)}_{_thousands(e.end)})"
                         f"{_copy_spec(e)}{_suffixes(e)}")
    return _tech_token(case, coords=True) + ",".join(items)


def _render_karyo_bands(case: Case, bands: CytobandTable, form: Form) -> str:
    """seq[build] del(10)(p15.3) — karyotype format inside a molecular line."""
    keep, _ = _karyotype_events(case, banded=False)
    _require_mechanism(keep)
    items = []
    for e in sorted(keep, key=order_key_karyotype):
        tok = _struct_token(e, bands, detailed=(form.system == "detailed"))
        if e.supernumerary:
            tok = "+" + tok                       # 11.4.2.3
        items.append(tok)
    prefix = _tech_token(case, coords=True)
    sexpart = ""
    if any(e.chrom in ("X", "Y") for e in keep) and case.sex_chromosomes:
        normal = _sex_remainder(case, keep)
        if normal:
            sexpart = normal + ","
    return prefix + sexpart + ",".join(items)


def _karyotype_count(case: Case) -> int:
    """Total chromosome number (4.2.1a).

    Sex chromosome aneuploidy is carried by the complement itself, not by
    +/- signs (5.3.1.1), so the sex contribution is read off
    ``sex_chromosomes``; only autosomal numerical changes and supernumerary
    structural chromosomes are added.
    """
    sx = case.sex_chromosomes or ""
    n_sex = 2 if sx in ("", "U") else len(sx.replace("c", ""))
    delta = sum(e.count_delta for e in case.events if e.chrom not in ("X", "Y"))
    delta += sum(1 for e in case.events
                 if e.chrom in ("X", "Y") and e.supernumerary)
    return 22 * case.ploidy + n_sex + delta


def _render_karyotype(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.2.1a: count, sex complement, abnormalities."""
    if case.sex_chromosomes is None:
        raise Unrenderable("karyotype", "кариотипная строка требует полового набора "
                                        "или обозначения U (5.2)")
    if any(e.mosaic_fraction is not None for e in case.events) and not case.cell_counts:
        raise Unrenderable("karyotype", "мозаичный результат без счёта клеток: "
                                        "кариотипная форма требует счёта метафаз "
                                        "(4.2.1c, 4.5.3b)")
    keep, _ = _karyotype_events(case)
    _require_mechanism(keep)
    count = _karyotype_count(case)
    sx = case.sex_chromosomes
    # 4.2.1b: a structurally abnormal sex chromosome is listed after the normal one
    parts: List[str] = []
    for e in sorted(keep, key=order_key_karyotype):
        if e.scale == "chromosome":
            if e.chrom in ("X", "Y"):
                continue          # carried by the complement (5.3.1.1)
            sign = "+" if e.direction == "gain" else "-"
            n = abs(e.copy_number - e.baseline)
            parts.extend([f"{sign}{e.chrom}"] * n)
        else:
            tok = _struct_token(e, bands)
            if e.supernumerary:
                tok = "+" + tok
            parts.append(tok)
    # the abnormal X/Y is written as a structural term, so it leaves the
    # complement and the normal homologue stays (4.2.1b, 5.5.11d)
    sx = _effective_sex(case, keep)
    body = ",".join([str(count), sx] + parts)
    return body


def _clone_order(clones: Sequence[Clone]) -> List[Clone]:
    """4.5.3d/e: largest abnormal line first, the normal line always last."""
    def key(c: Clone):
        size = c.cells if c.cells is not None else (c.fraction or 0)
        return (1 if c.normal else 0, -size)
    return sorted(clones, key=key)


def _clone_case(case: Case, clone: Clone) -> Case:
    return dataclasses.replace(
        case, events=clone.events, clones=[],
        sex_chromosomes=clone.sex_chromosomes or case.sex_chromosomes,
        ploidy=clone.ploidy, cell_counts={"present": 1})


def _render_clones(case: Case, bands: CytobandTable) -> str:
    parts = []
    for cl in _clone_order(case.clones):
        txt = _render_karyotype(_clone_case(case, cl), bands, FORMS["karyotype"])
        if cl.cells is not None:
            txt += f"[{cl.cells}]"
        parts.append(txt)
    return "mos " + "/".join(parts)                     # 4.4.1c, 4.5.3h


def _render_karyotype_mos(case: Case, bands: CytobandTable, form: Form) -> str:
    """4.5.3e/h: abnormal clones first, normal clone last, mos prefix."""
    if case.clones:
        if len(case.clones) < 2:
            raise Unrenderable("karyotype_mos", "задан один клон: форма клонов "
                                                "не применяется")
        return _render_clones(case, bands)
    if case.sex_chromosomes is None:
        raise Unrenderable("karyotype_mos", "требуется половой набор (5.2)")
    if any(e.chrom in ("X", "Y") and e.mosaic_fraction is not None
           for e in case.events):
        raise Unrenderable(
            "karyotype_mos", "мозаицизм по половой хромосоме: нормальный клон "
                             "по дозе не восстанавливается (45,X/46,XX и "
                             "45,X/46,XY дают одну дозу) — задайте клоны явно")
    frac = [e.mosaic_fraction for e in case.events if e.mosaic_fraction is not None]
    if not frac:
        raise Unrenderable("karyotype_mos", "немозаичный результат: форма клонов "
                                            "не применяется")
    abn = _render_karyotype(
        dataclasses.replace(
            case,
            events=[dataclasses.replace(e, mosaic_fraction=None) for e in case.events],
            cell_counts={"x": 1}),
        bands, FORMS["karyotype"])
    normal = f"{23 * case.ploidy},{case.sex_chromosomes}"
    counts = case.cell_counts or {}
    if counts:
        a = counts.get("abnormal")
        n = counts.get("normal")
        return f"mos {abn}[{a}]/{normal}[{n}]"
    return f"mos {abn}/{normal}"


_RENDERERS = {
    "abbrev": _render_abbrev,
    "abbrev_sex": _render_abbrev,
    "abbrev_range": _render_abbrev,
    "short": _render_short,
    "short_range": _render_short,
    "extended": _render_extended,
    "mixed": _render_mixed,
    "karyo_short": _render_karyo_bands,
    "karyo_detailed": _render_karyo_bands,
    "karyotype": _render_karyotype,
    "karyotype_mos": _render_karyotype_mos,
}


def _to_range_style(case: Case) -> Case:
    """8.1.1e alternative: express a mixed cell population as a copy-number
    range with a tilde instead of a proportion in square brackets."""
    evs = []
    for e in case.events:
        if (e.mosaic_fraction is not None or e.mosaic_fraction_range
                or e.mosaic_fraction_unknown) and not e.copy_number_range:
            lo, hi = sorted((e.baseline, e.copy_number))
            evs.append(dataclasses.replace(
                e, copy_number_range=(lo, hi), mosaic_fraction=None,
                mosaic_fraction_range=None, mosaic_fraction_unknown=False))
        else:
            evs.append(e)
    return dataclasses.replace(case, events=evs)


def _prepare(case: Case, form: str) -> Case:
    return _to_range_style(case) if form.endswith("_range") else case


def conformance_warnings(case: Case, form: str) -> List[str]:
    """Deviations that do not make the string ungrammatical but misreport how
    the result was obtained."""
    w: List[str] = []
    f = FORMS[form]
    if f.fmt == "karyotype" and f.sex_policy == "required" and not case.karyotype_available:
        w.append("кариотипный анализ не выполнялся: строка вида 46,XX относится к "
                 "бэндированному анализу метафаз, для результата секвенирования это "
                 "интерпретация, а не результат (11.2.2d)")
    if form == "karyotype_mos":
        w.append("mos с клонами предполагает счёт метафаз; для молекулярного "
                 "результата обозначение mos не применяется (4.5.2c, 4.5.3n)")
    if form.endswith("_range") and any(e.mosaic_fraction is not None
                                       for e in case.events):
        w.append("оценка доли известна численно, но диапазон копийности её не "
                 "передаёт (8.1.1e даёт обе записи как равнодопустимые)")
    if case.technique == "sseq" and any(
            e.start is not None and (e.end - e.start) < 5_000_000
            for e in case.events):
        w.append("для sseq стандарт указывает разрешение около 5–10 Мб (8.3b), "
                 "а заявлен сегмент меньше 5 Мб")
    return w


def render(case: Case, form: str, bands: CytobandTable) -> str:
    f = FORMS[form]
    case = _prepare(case, form)
    if case.clones and not case.events and f.fmt != "karyotype":
        raise Unrenderable(form, "задано клональное описание без суммарного "
                                 "профиля дозы: молекулярная форма требует "
                                 "измеренной копийности")
    if case.normal():
        return _render_normal(case, f)
    return _RENDERERS[form](case, bands, f)


def _normal_groups(case: Case) -> List[Tuple[List[str], int]]:
    """Chromosomes grouped by copy number, sex chromosomes first (8.1.1b/c).

    Reproduces the book's grouping: (X,1–22)×2 for a normal female,
    (X,Y)×1,(1–22)×2 for a male, (X)×2,(Y)×1,(1–22)×3 for 69,XXY.
    """
    sx = case.sex_chromosomes
    items: List[Tuple[str, int]] = []
    if sx not in (None, "U"):
        for letter in ("X", "Y"):
            n = sx.count(letter)
            if n:
                items.append((letter, n))
    items.append((f"1{NDASH}22", case.ploidy))
    groups: List[Tuple[List[str], int]] = []
    for tokname, n in items:
        for toks, gn in groups:
            if gn == n:
                toks.append(tokname)
                break
        else:
            groups.append(([tokname], n))
    return groups


def _render_normal(case: Case, form: Form) -> str:
    """11.3 / 8.2.1: a profile with no copy-number abnormality."""
    if form.fmt == "karyotype" and form.sex_policy == "required":
        sx = case.sex_chromosomes or "U"
        return f"{22 * case.ploidy + (2 if sx == 'U' else len(sx))},{sx}"
    if form.system in ("extended",) or form.needs_structure:
        raise Unrenderable(form.key, "нет аберраций: форма неприменима")
    tech = _tech_token(case, coords=False)
    cht = "cht" if case.chromatid else ""
    zyg = case.genome_zygosity or ""
    body = ",".join(f"({','.join(toks)}){TIMES}{n}" for toks, n in _normal_groups(case))
    return f"{tech}{cht}{body}{zyg}"


def render_all(case: Case, bands: CytobandTable) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for key in FORMS:
        try:
            out[key] = {"text": render(case, key, bands), "admissible": True,
                        "reason": ""}
        except Unrenderable as exc:
            out[key] = {"text": None, "admissible": False, "reason": exc.reason}
    return out


# ---------------------------------------------------------------------------
# Independent parser: a grammar, not a call back into the renderer.
# ---------------------------------------------------------------------------

GRAMMAR = r"""
start: line (PERIOD line)*

?line: micro | karyo

micro: TECH BUILD? SP? (MOSCHI SP)? body_items (SLASH body_items)*
body_items: (SEXCX COMMA)? body_item ((COMMA | SP OR SP) body_item)*
?body_item: group | short_item | extended_item | cx_item | karyo_item_alt

group: CHT? LPAR chromlist RPAR copyspec suffixes
chromlist: chromtok (COMMA chromtok)*
chromtok: CHR ARM?      -> single
        | CHR NDASH CHR -> crange

short_item: BANDSPAN LPAR NUM UNDER NUM RPAR copyspec suffixes
extended_item: BANDSPAN LPAR ext_part (COMMA ext_part)+ RPAR suffixes
ext_span: NUM UNDER NUM copyspec
ext_point: NUM copyspec
?ext_part: ext_span | ext_point
cx_item: (LPAR chromlist RPAR | BANDSPAN LPAR NUM UNDER NUM RPAR) CXSYM frac?

copyspec: TIMES GT? INT (TILDE INT)?
        | AMP
suffixes: zyg? inh? ctail? cxtail? frac?
ctail: SP? CONST
cxtail: CXSYM
zyg: ZYG
inh: SP? INH
frac: LBRACK FRACVAL RBRACK

karyo: clonepfx? clone (SLASH clone)*
clonepfx: MOSCHI SP
clone: COUNT (TILDE COUNT)? COMMA clone_body? cells?
clone_body: SEXCX (COMMA karyo_item_alt)*
          | karyo_item_alt (COMMA karyo_item_alt)*
cells: LBRACK INT RBRACK

karyo_item_alt: karyo_item (SP OR SP karyo_item)*
?karyo_item: numerical | struct_chain
?numerical: numchrom | numstruct
numchrom: SIGN CHR (SP? TAIL)?
numstruct: SIGN INT? struct_chain
struct_chain: structural+ (SP? TAIL)? copies?
copies: TIMES INT
structural: QM? (MODIF SP)? SYM LPAR chromgroup RPAR (LPAR bpgroup RPAR)?  -> plain
          | QM? (MODIF SP)? SYM LPAR chromgroup RPAR LPAR DETAIL RPAR      -> detailed
          | QM? (MODIF SP)? SYM LPAR DETAIL RPAR                           -> detailed
          | QM? (MODIF SP)? BARESYM                                        -> bare
chromgroup: chrtok? (SEMI chrtok?)*
?chrtok: CHR ARM? | QM
bpgroup: bandalt (SEMI bandalt)*
bandalt: bandtok (SP OR SP bandtok)*
bandtok: bandpart+
       | QM
?bandpart: BAND | TILDE BAND | TILDE INT | NDASH BAND | NDASH INT

TECH: "sseq" | "seq" | "arr" | "ogm"
BUILD: /\[[Gg][Rr][Cc][Hh]3[78]\]|\[hg(19|38)\]/
CHT: "cht"
GT: ">"
CONST: "c"
MOSCHI: "mos" | "chi"
OR: "or"
MODIF: "psu" | "neo"
CXSYM: "cth" | "cha" | "cpx" | "cx"
SYM: "idic" | "ider" | "der" | "dic" | "del" | "dup" | "inv" | "ins" | "rec" | "rob" | "trp" | "qdp" | "add" | "fis" | "trc" | "i" | "r" | "t"
BARESYM: "mar" | "ace"
INH: "mat pat" | "dmat" | "dpat" | "dinh" | "umat" | "upat" | "inh" | "mat" | "pat" | "dn"
ZYG: "hmz htz" | "hmz" | "htz"
AMP: "amp"
TAIL: ("c" | "dn" | "mat pat" | "mat" | "pat" | "inh" | "dmat" | "dpat" | "dinh")+
CHR: /X|Y|1[0-9]|2[0-2]|[1-9]/
ARM: /[pq]/
BAND: /[pq](?:[0-9?]*(?:\.[0-9?]+)?)/
BANDSPAN: /(chr)?(X|Y|1[0-9]|2[0-2]|[1-9])[pq][0-9]+(\.[0-9]+)?([~\u2013]?[pq]?[0-9]+(\.[0-9]+)?)?/
DETAIL: /[^()]*(::|\u2192)[^()]*/
NUM: /[0-9][0-9,]*/
INT: /[0-9]+/
COUNT: /[0-9]{2,3}/
FRACVAL: /(0?\.[0-9]+|1\.0)(~(0?\.[0-9]+|1\.0))?|\?/
SEXCX: /(X+Y*|Y+|U)\??c?/
SIGN: "+" | "-" | "\u2212"
TIMES: "\u00d7" | "x"
NDASH: "\u2013" | "-"
TILDE: "~"
UNDER: "_"
COMMA: ","
SEMI: ";"
SLASH: "/"
PERIOD: "."
LPAR: "("
RPAR: ")"
LBRACK: "["
RBRACK: "]"
SP: " "
QM: "?"
"""

_PARSER = None
_PARSER_AMBIG = None


def _get_parser(ambiguity: str = "resolve"):
    """Earley parser.  ``ambiguity='explicit'`` exposes every derivation, which
    is how grammatical ambiguity is measured rather than assumed."""
    global _PARSER, _PARSER_AMBIG
    from lark import Lark
    if ambiguity == "explicit":
        if _PARSER_AMBIG is None:
            _PARSER_AMBIG = Lark(GRAMMAR, start="start", parser="earley",
                                 ambiguity="explicit", keep_all_tokens=True)
        return _PARSER_AMBIG
    if _PARSER is None:
        _PARSER = Lark(GRAMMAR, start="start", parser="earley",
                       keep_all_tokens=True)
    return _PARSER


class ParseError(Exception):
    pass


class OutOfScope(Exception):
    """Valid ISCN that this parser deliberately does not cover."""


OUT_OF_SCOPE_MARKERS = ("ish", "nuc", "rsa", "//", "cp", "sl", "sdl", "idem",
                        "VAF", "dmin", "inc", "<", "fra", "hsr", "mar1", "mar2", "neo",
                        "mar3", "tas", "sep", "cha", "cth", "end", "nvt")


def _clean_num(tok: str) -> int:
    return int(tok.replace(",", ""))


def _split_bandspan(tok: str) -> Tuple[str, List[str]]:
    tok = tok.replace("chr", "").replace("~", "").replace("\u2013", "")
    m = re.fullmatch(r"(X|Y|1[0-9]|2[0-2]|[1-9])((?:[pq][0-9]+(?:\.[0-9]+)?)+)", tok)
    if m is None:
        raise ParseError(f"не разбирается обозначение полос {tok!r}")
    chrom, rest = m.group(1), m.group(2)
    bands = re.findall(r"[pq][0-9]+(?:\.[0-9]+)?", rest)
    return chrom, bands


def parse(text: str, ambiguity: str = "resolve") -> Dict[str, Any]:
    """Parse an ISCN string into facts.  Raises ParseError or OutOfScope."""
    from lark import Lark  # noqa: F401  (import error surfaces here, not later)
    from lark.exceptions import UnexpectedInput, VisitError
    stripped = text.strip()
    for m in OUT_OF_SCOPE_MARKERS:
        if re.search(rf"(?<![a-z]){re.escape(m)}", stripped):
            if m == "sl" and "slx" not in stripped and not stripped.startswith("sl "):
                continue
            raise OutOfScope(f"конструкция {m!r} вне области разборщика")
    try:
        tree = _get_parser(ambiguity).parse(stripped)
    except UnexpectedInput as exc:
        raise ParseError(str(exc).split("\n")[0]) from None
    return _facts_from_tree(tree, stripped)


def check_conventions(text: str, bands: Optional[CytobandTable] = None
                     ) -> List[Dict[str, str]]:
    """Audit a formula against the ISCN rules that a grammar cannot enforce.

    Returns one entry per violation with the rule number, so a rejected string
    is rejected for a stated reason rather than by opinion.  Used to screen
    records received from outside and as the rejection corpus check.
    """
    out: List[Dict[str, str]] = []

    STYLE = {"конвенция", "4.4.3a", "4.4.5d", "4.4.5e", "4.4.1a", "4.5.3h"}

    def flag(rule: str, msg: str) -> None:
        out.append({"rule": rule, "message": msg,
                    "severity": "style" if rule in STYLE else "error"})

    t = text.strip()
    parsed = parse(t)                    # raises for ungrammatical input
    meta = parsed["meta"]

    # --- typography ------------------------------------------------------
    if re.search(r"\)x[0-9]", t) or re.search(r"[0-9]\)x", t):
        flag("4.4.3a", "знак копийности записан латинской x вместо ×")
    if re.search(r"[0-9]{3}[,\u2013-][0-9]{3},[0-9]{3}\)", t) and "_" not in t:
        flag("4.4.5e", "границы участка разделены дефисом, стандарт требует "
                       "подчёркивания")
    if re.search(r"\u2212", t):
        # В перечне обозначений ISCN знак потери печатается дефисом; отдельного
        # пункта о недопустимости типографского минуса в тексте нет, поэтому
        # замечание выдаётся как внутренняя конвенция, а не как пункт стандарта.
        flag("конвенция", "использован типографский минус U+2212 вместо дефиса: "
                          "в машинной обработке это разные знаки")
    if "[" in t and re.search(r"\((?:[0-9]{1,3},)+[0-9]{3}\u00d7", t):
        flag("4.4.5d", "в расширенной системе присутствуют разделители тысяч")

    # --- build vs coordinates -------------------------------------------
    has_coords = any(f[0] == "cn" and f[3] is not None for f in parsed["facts"])
    if meta.get("technique"):
        if has_coords and not meta.get("build"):
            flag("4.4.5b", "приведены координаты, но не указана сборка генома")
        karyo_format = bool(re.search(r"[a-z]{1,4}\([0-9XY]", t.split(" ", 1)[-1]))
        if not has_coords and meta.get("build") and not karyo_format:
            flag("4.4.5c", "сборка указана, хотя координаты не приводятся")
        if meta["technique"] and " " not in t.split(meta["technique"], 1)[1][:1]:
            if not t.split(meta["technique"], 1)[1].startswith("["):
                flag("4.4.1a", "нет пробела между обозначением методики и результатом")

    # --- clones ----------------------------------------------------------
    clones = meta.get("clones") or []
    if len(clones) > 1 and meta.get("clone_prefix") is None:
        flag("4.5.3h", "несколько клонов без обозначения mos: стандарт допускает обе записи, для конституционального образца обозначение предпочтительно")
    if len(clones) == 1 and meta.get("clone_prefix"):
        flag("4.5.3h", "обозначение mos при единственном клоне")
    bodies = re.split(r"/", t.split(" ", 1)[-1] if meta.get("clone_prefix") else t)
    if len(bodies) > 1:
        def _is_normal_clone(b: str) -> bool:
            m = re.fullmatch(r"([0-9]{2,3}),(XX|XY|U)(?:\[[0-9]+\])?", b.strip())
            return bool(m) and int(m.group(1)) == 44 + (2 if m.group(2) == "U"
                                                        else len(m.group(2)))

        normal_at = [i for i, b in enumerate(bodies) if _is_normal_clone(b)]
        if normal_at and max(normal_at) != len(bodies) - 1:
            flag("4.5.3e", "нормальный клон записан не последним")
        counts = [int(m.group(1)) for b in bodies
                  if (m := re.search(r"\[([0-9]+)\]", b))]
        if len(counts) == len(bodies) and len(counts) > 1:
            abn = counts[:-1] if normal_at and max(normal_at) == len(bodies) - 1 \
                else counts
            if abn != sorted(abn, reverse=True):
                flag("4.5.3d", "клоны перечислены не по убыванию числа клеток")

    # --- ordering of abnormalities --------------------------------------
    for body in bodies:
        items = _split_top_level(body)
        keys = []
        for it in items:
            k = _order_probe(it)
            if k is not None:
                keys.append((k, it))
        seq = [k for k, _ in keys]
        if seq != sorted(seq):
            flag("4.3b", "порядок аберраций нарушен: ожидается X, Y, затем "
                         "аутосомы по возрастанию, внутри хромосомы — "
                         "числовые прежде структурных и от pter к qter "
                         f"(получено {[i for _, i in keys]})")

    # --- chromosome count ------------------------------------------------
    # Only the constitutional convention is checked: the sex contribution is
    # read off the complement plus one structural sex term per item (4.2.1b),
    # autosomal +/- change the count, and a whole-arm der/dic/trc replaces
    # several chromosomes with one.  Acquired notation (a c suffix, or a sex
    # +/- alongside a constitutional complement) is out of scope and skipped.
    for body in bodies:
        first_line = body.split(".")[0]
        m = re.match(r"^([0-9]{2,3})(?:~[0-9]{2,3})?,(X+Y*|Y+|U)(c?)", first_line.strip())
        if not m:
            continue
        if m.group(3) == "c" or re.search(r"[0-9]c(?![a-z])", first_line):
            continue                      # constitutional marker: neoplastic context
        if "cht" in first_line or re.search(r"\)\u00d7[0-9]", first_line):
            continue                      # chromatid counting (8.3) and copies of a
                                          # rearranged chromosome (5.6): the count
                                          # rule is not calibrated for these
        declared, sx = int(m.group(1)), m.group(2)
        rest = first_line[m.end():]
        if re.search(r"[+\-\u2212](X|Y)(?![0-9])", rest):
            continue                      # acquired sex change: out of scope
        n_sex = 2 if sx == "U" else len(sx)
        delta = 0
        for item in _split_top_level(rest):
            it = item.split(" or ")[0].strip()
            if not it:
                continue
            if re.fullmatch(r"\[[0-9]+\]", it):
                continue
            sign = it[0] if it[0] in "+-\u2212" else ""
            core = it[1:] if sign else it
            if re.fullmatch(r"[0-9]{1,2}", core):
                delta += 1 if sign == "+" else -1
                continue
            if re.match(r"^[a-z]{1,4}\((?:X|Y)[;)]", core) or \
               re.match(r"^[a-z]{1,4} ?[a-z]{0,4}\((?:X|Y)[;)]", core):
                n_sex += 1                # structural sex chromosome (4.2.1b)
            if sign == "+":
                delta += 1
            # a whole-arm rearrangement fuses several chromosomes into one
            for tm in re.finditer(r"(?:psu )?(der|dic|trc|rob)\(([^)]*)\)"
                                  r"(?:\(([^)]*)\))?", core):
                chroms = [c for c in tm.group(2).split(";") if c]
                bps = tm.group(3) or ""
                whole_arm = (tm.group(1) in ("dic", "trc", "rob")
                             or bool(re.search(r"[pq]10", bps)))
                if len(chroms) > 1 and whole_arm:
                    delta -= len(chroms) - 1
        expected = 22 * (n_sex and 2 or 2) + n_sex + delta
        options = {22 * n + n_sex + delta for n in (1, 2, 3, 4)}
        if declared not in options:
            flag("4.2.1a", f"число хромосом {declared} не согласуется с записью: "
                           f"при плоидности 1-4 ожидается одно из "
                           f"{sorted(options)} ({n_sex} половых, изменение "
                           f"{delta:+d})")

    # --- bands vs coordinates -------------------------------------------
    if bands is not None and has_coords:
        for c in check_band_consistency(t, bands):
            if c["relation"] == "inconsistent":
                flag("4.4.5a", f"полоса {c['declared']} не соответствует "
                               f"координатам {c['start']}_{c['end']}: по "
                               f"{c['build']} это {c['implied']}")
    return out


def _split_top_level(body: str) -> List[str]:
    """Split a karyotype body on commas that are not inside parentheses."""
    out, depth, buf = [], 0, ""
    for ch in body:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(buf)
            buf = ""
        else:
            buf += ch
    out.append(buf)
    return out


def _order_probe(item: str) -> Optional[Tuple]:
    """Sort key of one written abnormality, for the ordering check (4.3)."""
    it = item.strip()
    if not it or re.fullmatch(r"[0-9]{2,3}|X+Y*|Y+|Uc?|\[[0-9]+\]", it):
        return None
    m = re.match(r"^([+\-\u2212])?\s*([0-9]{1,2}|X|Y)$", it)
    if m:                                     # numerical
        return chrom_sort_key(m.group(2)) + (0, 0 if (m.group(1) or "+") == "+" else 1, "")
    m = re.match(r"^([+\-\u2212])?\s*[a-z]+\(([0-9]{1,2}|X|Y)", it)
    if m:                                     # structural
        band = re.search(r"\)\(([pq])", it)
        arm_rank = 0 if (band and band.group(1) == "p") else 1
        return chrom_sort_key(m.group(2)) + (1, arm_rank, it)
    return None


def _is_coarser(declared: str, b1: str, b2: str) -> bool:
    """A declared designation may name the parent band of the implied one."""
    parts = re.findall(r"[pq][0-9]+(?:\.[0-9]+)?", declared)
    imp = [b1] if b1 == b2 else [b1, b2]
    if len(parts) != len(imp):
        return False
    return all(i.startswith(p) for p, i in zip(parts, imp))


def check_band_consistency(text: str, bands: CytobandTable) -> List[Dict[str, Any]]:
    """Audit an existing string: are its band designations the ones its own
    coordinates imply (4.4.5a/f)?

    This is the invariant a generator satisfies by construction and a
    hand-written record does not, so it is the cheapest check to run against
    records received from outside.
    """
    out: List[Dict[str, Any]] = []
    parsed = parse(text)
    for f in parsed["facts"]:
        if f[0] != "cn" or f[3] is None:
            continue
        _, chrom, banddes, start, end = f[0], f[1], f[2], f[3], f[4]
        b1 = bands.band_at(chrom, start)
        b2 = bands.band_at(chrom, end)
        implied = b1 if b1 == b2 else f"{b1}{b2}"
        if banddes == implied:
            relation = "exact"
        elif banddes and _is_coarser(banddes, b1, b2):
            relation = "coarser"        # parent band: allowed low resolution
        else:
            relation = "inconsistent"
        out.append({"chromosome": chrom, "declared": banddes, "implied": implied,
                    "start": start, "end": end, "relation": relation,
                    "consistent": relation != "inconsistent",
                    "build": bands.build})
    return out


def count_derivations(text: str) -> int:
    """Number of distinct parse trees (grammatical ambiguity check)."""
    from lark import Tree
    tree = _get_parser("explicit").parse(text.strip())
    n = 1
    for sub in tree.iter_subtrees():
        if sub.data == "_ambig":
            n *= len(sub.children)
    return n


def _facts_from_tree(tree, text: str) -> Dict[str, Any]:
    """Walk the tree deliberately (not iter_subtrees) so that a term nested
    under another rule is counted once and in its correct role."""
    from lark import Tree, Token

    def tok(node, name):
        return next((c for c in node.children
                     if isinstance(c, Token) and c.type == name), None)

    def toks(node, name):
        return [c for c in node.children
                if isinstance(c, Token) and c.type == name]

    def kids(node, *names):
        return [c for c in node.children
                if isinstance(c, Tree) and (not names or c.data in names)]

    facts: List[Tuple] = []
    meta: Dict[str, Any] = {"technique": None, "build": None, "chromatid": False,
                            "lines": 0, "clones": [], "clone_prefix": None}

    def read_suffixes(node):
        zyg = inh = fr = None
        for s in kids(node, "suffixes"):
            for z in kids(s, "zyg"):
                zyg = str(tok(z, "ZYG"))
            for i in kids(s, "inh"):
                inh = str(tok(i, "INH"))
            for f in kids(s, "frac"):
                fr = str(tok(f, "FRACVAL"))
            for c in kids(s, "ctail"):
                meta["constitutional"] = True
            for c in kids(s, "cxtail"):
                facts.append(("complex", str(tok(c, "CXSYM"))))
        return zyg, inh, fr

    def read_cs(cs):
        ints = [int(str(c)) for c in cs.children
                if isinstance(c, Token) and c.type == "INT"]
        if not ints:
            return ("amp", None)
        return ("range", tuple(ints)) if len(ints) == 2 else ("cn", ints[0])

    def first_cs(node):
        cs = kids(node, "copyspec")
        return read_cs(cs[0]) if cs else (None, None)

    def read_term(node) -> Tuple:
        if node.data == "detailed":
            det = str(tok(node, "DETAIL"))
            return (str(tok(node, "SYM")), (),
                    tuple(s.strip() for s in det.split("::")))
        if node.data == "bare":
            return (str(tok(node, "BARESYM")), (), ())
        sym = str(tok(node, "SYM"))
        chroms = []
        for cg in kids(node, "chromgroup"):
            buf = ""
            for c in cg.scan_values(lambda v: True):
                if not isinstance(c, Token):
                    continue
                if c.type == "CHR":
                    if buf:
                        chroms.append(buf)
                    buf = str(c)
                elif c.type == "ARM":
                    buf += str(c)
                elif c.type == "QM":
                    if buf:
                        chroms.append(buf)
                    buf = "?"
            if buf:
                chroms.append(buf)
        chroms = tuple(chroms)
        groups = []
        for bg in kids(node, "bpgroup"):
            for ba in kids(bg, "bandalt"):
                for bt in kids(ba, "bandtok"):
                    groups.append(tuple(
                        str(c) for c in bt.children
                        if isinstance(c, Token) and c.type in ("BAND", "INT", "QM")))
        return (sym, chroms, tuple(groups))

    def read_chain(node) -> Tuple[Tuple, Optional[str]]:
        terms = tuple(read_term(t) for t in kids(node, "plain", "detailed", "bare"))
        tail = tok(node, "TAIL")
        return terms, (str(tail) if tail else None)

    def walk(node) -> None:
        d = node.data
        if d == "micro":
            meta["lines"] += 1
            meta["technique"] = str(tok(node, "TECH"))
            b = tok(node, "BUILD")
            meta["build"] = str(b)[1:-1] if b else None
            meta["chromatid"] = tok(node, "CHT") is not None
            for k in kids(node):
                walk(k)
        elif d == "group":
            kind, cn = first_cs(node)
            zyg, inh, fr = read_suffixes(node)
            for cl in kids(node, "chromlist"):
                for ct in kids(cl):
                    if ct.data == "single":
                        chrom = str(tok(ct, "CHR"))
                        arm = tok(ct, "ARM")
                        facts.append(("cn", chrom, str(arm) if arm else None,
                                      None, None, None,
                                      cn if kind == "cn" else None,
                                      cn if kind == "range" else None,
                                      fr, inh, zyg))
                    else:
                        chs = toks(ct, "CHR")
                        facts.append(("cnrange_group", str(chs[0]), str(chs[1]),
                                      cn, fr, zyg))
        elif d == "short_item":
            chrom, bnds = _split_bandspan(str(tok(node, "BANDSPAN")))
            nums = toks(node, "NUM")
            kind, cn = first_cs(node)
            zyg, inh, fr = read_suffixes(node)
            facts.append(("cn", chrom, "".join(bnds), _clean_num(str(nums[0])),
                          _clean_num(str(nums[1])), None,
                          cn if kind == "cn" else None,
                          cn if kind == "range" else None, fr, inh, zyg))
        elif d == "extended_item":
            chrom, bnds = _split_bandspan(str(tok(node, "BANDSPAN")))
            zyg, inh, fr = read_suffixes(node)
            parts = []
            for p in kids(node, "ext_span", "ext_point"):
                ns = [_clean_num(str(c)) for c in toks(p, "NUM")]
                cs = read_cs(kids(p, "copyspec")[0])
                parts.append((p.data, ns, cs))
            for i, (kind, ns, cs) in enumerate(parts):
                if kind != "ext_span":
                    continue
                left = next((q[1][0] for q in reversed(parts[:i])
                             if q[0] == "ext_point"), None)
                right = next((q[1][0] for q in parts[i + 1:]
                              if q[0] == "ext_point"), None)
                facts.append(("cn", chrom, "".join(bnds), ns[0], ns[1],
                              (left, right), cs[1] if cs[0] == "cn" else None,
                              cs[1] if cs[0] == "range" else None, fr, inh, zyg))
        elif d == "cx_item":
            facts.append(("complex", str(tok(node, "CXSYM"))))
        elif d == "struct_chain":
            terms, tail = read_chain(node)
            chrom = terms[0][1][0] if terms and terms[0][1] else None
            facts.append(("chain", chrom, terms, tail))
        elif d == "numchrom":
            facts.append(("num", str(tok(node, "CHR")), str(tok(node, "SIGN"))))
        elif d == "numstruct":
            sign = str(tok(node, "SIGN"))
            for sc in kids(node, "struct_chain"):
                terms, tail = read_chain(sc)
                facts.append(("num_chain", sign, terms, tail))
        elif d == "clone":
            counts = [str(c) for c in toks(node, "COUNT")]
            sx = tok(node, "SEXCX")
            cells = None
            for cc in kids(node, "cells"):
                cells = int(str(tok(cc, "INT")))
            meta["clones"].append({"count": counts, "cells": cells})
            facts.append(("count", "~".join(counts)))
            for k in kids(node, "clone_body"):
                walk(k)
        elif d in ("body_items", "clone_body"):
            sx = tok(node, "SEXCX")
            if sx is not None:
                facts.append(("sex", str(sx)))
            for k in kids(node):
                walk(k)
        elif d == "clonepfx":
            meta["clone_prefix"] = str(tok(node, "MOSCHI"))
        else:
            for k in kids(node):
                walk(k)

    walk(tree)
    return {"facts": _index_facts(facts), "meta": meta, "text": text}


def _index_facts(items: Sequence[Tuple]) -> frozenset:
    """Turn a list of facts into a set without losing multiplicity: a fact that
    occurs k times becomes k distinct entries (+21,+21 must not collapse)."""
    from collections import Counter
    seen: Counter = Counter()
    out = []
    for it in items:
        seen[it] += 1
        out.append(it + (seen[it],))
    return frozenset(out)


def _unused_facts_from_tree(tree, text: str) -> Dict[str, Any]:
    from lark import Tree, Token

    def tok(node, name):
        return next((c for c in node.children
                     if isinstance(c, Token) and c.type == name), None)

    def subs(node, name):
        return [c for c in node.children if isinstance(c, Tree) and c.data == name]

    facts: List[Tuple] = []
    meta: Dict[str, Any] = {"technique": None, "build": None, "chromatid": False,
                            "lines": 0, "clones": []}

    def read_suffixes(node) -> Tuple[Optional[str], Optional[str], Any]:
        zyg = inh = fr = None
        for s in subs(node, "suffixes"):
            for z in subs(s, "zyg"):
                zyg = str(tok(z, "ZYG"))
            for i in subs(s, "inh"):
                inh = str(tok(i, "INH"))
            for f in subs(s, "frac"):
                v = str(tok(f, "FRACVAL"))
                fr = v
        return zyg, inh, fr

    def read_copyspec(node):
        for cs in subs(node, "copyspec"):
            ints = [c for c in cs.children if isinstance(c, Token) and c.type == "INT"]
            if not ints:
                return ("amp", None)
            if len(ints) == 2:
                return ("range", (int(ints[0]), int(ints[1])))
            return ("cn", int(ints[0]))
        return (None, None)

    for line in tree.iter_subtrees():
        if line.data == "micro":
            meta["lines"] += 1
            meta["technique"] = str(tok(line, "TECH"))
            b = tok(line, "BUILD")
            meta["build"] = str(b)[1:-1] if b else None
            meta["chromatid"] = tok(line, "CHT") is not None

    for node in tree.iter_subtrees():
        if node.data == "group":
            kind, cn = read_copyspec(node)
            zyg, inh, fr = read_suffixes(node)
            for cl in subs(node, "chromlist"):
                for ct in cl.children:
                    if not isinstance(ct, Tree):
                        continue
                    if ct.data == "single":
                        chrom = str(tok(ct, "CHR"))
                        arm = tok(ct, "ARM")
                        facts.append(("cn", chrom, str(arm) if arm else None,
                                      None, None, None, cn if kind == "cn" else None,
                                      cn if kind == "range" else None, fr, inh, zyg))
                    elif ct.data == "crange":
                        chs = [str(c) for c in ct.children
                               if isinstance(c, Token) and c.type == "CHR"]
                        facts.append(("cnrange_group", chs[0], chs[1], cn, fr))
        elif node.data == "short_item":
            chrom, bnds = _split_bandspan(str(tok(node, "BANDSPAN")))
            nums = [c for c in node.children if isinstance(c, Token) and c.type == "NUM"]
            kind, cn = read_copyspec(node)
            zyg, inh, fr = read_suffixes(node)
            facts.append(("cn", chrom, "".join(bnds), _clean_num(nums[0]),
                          _clean_num(nums[1]), None, cn if kind == "cn" else None,
                          cn if kind == "range" else None, fr, inh, zyg))
        elif node.data == "extended_item":
            chrom, bnds = _split_bandspan(str(tok(node, "BANDSPAN")))
            nums = [_clean_num(str(c)) for c in node.children
                    if isinstance(c, Token) and c.type == "NUM"]
            specs = [read_copyspec_child(c) for c in node.children
                     if isinstance(c, Tree) and c.data == "copyspec"]
            zyg, inh, fr = read_suffixes(node)
            facts.append(("cn", chrom, "".join(bnds), nums[1], nums[2],
                          (nums[0], nums[3]), specs[1], None, fr, inh, zyg))
        elif node.data == "structural":
            sym = str(tok(node, "SYM"))
            chroms = [str(c) for cg in subs(node, "chromgroup") for c in cg.children
                      if isinstance(c, Token) and c.type == "CHR"]
            bps = []
            for bg in subs(node, "bpgroup"):
                for bt in subs(bg, "bandtok"):
                    bps.append("".join(str(c) for c in bt.children))
            tail = tok(node, "TAIL")
            facts.append(("struct", chroms[0] if chroms else None, sym,
                          tuple(chroms), tuple(bps),
                          str(tail) if tail else None))
        elif node.data == "detailed":
            sym = str(tok(node, "SYM"))
            det = str(tok(node, "DETAIL"))
            segs = tuple(s.strip() for s in det.split("::"))
            facts.append(("struct_detailed", sym, segs))
        elif node.data == "numerical":
            sign = str(tok(node, "SIGN"))
            ch = tok(node, "CHR")
            if ch is not None:
                facts.append(("num", str(ch), sign))
            else:
                facts.append(("num_struct", sign))
        elif node.data == "clone":
            cnt = [str(c) for c in node.children
                   if isinstance(c, Token) and c.type == "COUNT"]
            sx = tok(node, "SEXCX")
            cells = None
            for cc in subs(node, "cells"):
                cells = int(str(tok(cc, "INT")))
            meta["clones"].append({"count": cnt, "sex": str(sx) if sx else None,
                                   "cells": cells})
        elif node.data == "cx_item":
            cxs = str(tok(node, "CXSYM"))
            facts.append(("complex", cxs))

    return {"facts": frozenset(facts), "meta": meta, "text": text}


def read_copyspec_child(cs):
    from lark import Token
    ints = [int(str(c)) for c in cs.children
            if isinstance(c, Token) and c.type == "INT"]
    return ints[0] if len(ints) == 1 else tuple(ints)


# ---------------------------------------------------------------------------
# Projection: what a form is able to carry (executable loss matrix)
# ---------------------------------------------------------------------------


def _frac_text(e: Event) -> Optional[str]:
    if e.mosaic_fraction_unknown:
        return "?"
    if e.mosaic_fraction_range:
        lo, hi = e.mosaic_fraction_range
        return f"{lo:g}~{hi:g}"
    if e.mosaic_fraction is not None:
        return f"{e.mosaic_fraction:g}"
    return None


def _karyo_line_facts(case: Case, bands: CytobandTable, detailed: bool,
                      banded: bool, with_count: bool) -> List[Tuple]:
    """Facts one karyotype-format line carries, in the shape parse() returns."""
    out: List[Tuple] = []
    events, _ = _karyotype_events(case, banded=banded)
    for e in events:
        if e.scale == "chromosome":
            if e.chrom in ("X", "Y"):
                continue          # carried by the complement (5.3.1.1)
            n = abs(e.copy_number - e.baseline)
            sign = "+" if e.direction == "gain" else "-"
            out.extend([("num", e.chrom, sign)] * n)
        else:
            terms = []
            for t in _struct_chain(e, bands):
                if detailed and t.get("detailed"):
                    terms.append((t["symbol"], (),
                                  tuple(s.strip() for s in t["detailed"].split("::"))))
                else:
                    terms.append((t["symbol"], tuple(t.get("chroms") or ()),
                                  _bp_groups(t)))
            tail = _struct_tail(e) or None
            if e.supernumerary:
                out.append(("num_chain", "+", tuple(terms), tail))
            else:
                out.append(("chain", terms[0][1][0] if terms[0][1] else None,
                            tuple(terms), tail))
    if with_count:
        out.append(("count", str(_karyotype_count(case))))
        out.append(("sex", _effective_sex(case, events)))
    elif any(e.chrom in ("X", "Y") for e in events) and case.sex_chromosomes:
        rem = _sex_remainder(case, events)
        if rem:
            out.append(("sex", rem))
    return out


def project(case: Case, form: str, bands: CytobandTable) -> frozenset:
    """The fact set a form can express, in the same shape parse() returns."""
    f = FORMS[form]
    keep_attr = LOSS_MATRIX[form]
    case = _prepare(case, form)
    out: List[Tuple] = []

    if case.normal():
        if f.fmt == "karyotype" and f.sex_policy == "required":
            sx = case.sex_chromosomes or "U"
            out.append(("count", str(22 * case.ploidy +
                                     (2 if sx == "U" else len(sx)))))
            out.append(("sex", sx))
            return _index_facts(out)
        zyg = case.genome_zygosity
        for toks_, n in _normal_groups(case):
            for t in toks_:
                if NDASH in t or "-" in t:
                    a, b = re.split(r"[-\u2013]", t)
                    out.append(("cnrange_group", a, b, n, None, zyg))
                else:
                    out.append(("cn", t, None, None, None, None, n, None,
                                None, None, zyg))
        return _index_facts(out)

    if f.fmt == "microarray":
        for e in case.events:
            s, en = e.start, e.end
            if f.system == "abbreviated" or s is None:
                arm = None
                if e.scale == "arm" or (e.scale != "chromosome"
                                        and bands.covers_arm(e.chrom, s, en)):
                    arm = e.arm or bands.covers_arm(e.chrom, s, en)
                bandtxt, start, end = arm, None, None
            else:
                bandtxt = _band_designation(
                    dataclasses.replace(e, start=s, end=en), bands,
                    microarray=True)[len(e.chrom):]
                start, end = s, en
            flanks = ((e.flank_left, e.flank_right)
                      if keep_attr["flanking_normal"] else None)
            out.append(("cn", e.chrom, bandtxt, start, end, flanks,
                        e.copy_number if not e.copy_number_range else None,
                        e.copy_number_range,
                        _frac_text(e) if keep_attr["mosaic_fraction"] else None,
                        e.inheritance if keep_attr["parental_origin"] else None,
                        e.zygosity if keep_attr["zygosity"] else None))
        if keep_attr["sex_complement"] and case.sex_chromosomes not in (None, "U"):
            for letter in ("X", "Y"):
                n = case.sex_chromosomes.count(letter)
                if n:
                    out.append(("cn", letter, None, None, None, None, n, None,
                                None, None, None))

    elif f.fmt == "mixed":
        keep, _ = _karyotype_events(case, banded=False)
        for e in keep:
            if e.structure:
                terms = tuple((t["symbol"], tuple(t.get("chroms") or ()),
                               _bp_groups(t))
                              for t in _struct_chain(e, bands))
                tail = _struct_tail(e) or None
                if e.supernumerary:
                    out.append(("num_chain", "+", terms, tail))
                else:
                    out.append(("chain", terms[0][1][0] if terms[0][1] else None,
                                terms, tail))
            else:
                band = _band_designation(e, bands, microarray=True)[len(e.chrom):]
                out.append(("cn", e.chrom, band, e.start, e.end, None,
                            e.copy_number if not e.copy_number_range else None,
                            e.copy_number_range, _frac_text(e), e.inheritance,
                            e.zygosity))

    else:
        detailed = (f.system == "detailed")
        banded = (f.sex_policy == "required")
        if form == "karyotype_mos" and case.clones:
            for cl in _clone_order(case.clones):
                sub = _clone_case(case, cl)
                out.extend(_karyo_line_facts(sub, bands, detailed, banded, True))
        elif form == "karyotype_mos":
            stripped = dataclasses.replace(
                case,
                events=[dataclasses.replace(e, mosaic_fraction=None)
                        for e in case.events],
                cell_counts={"present": 1})
            out.extend(_karyo_line_facts(stripped, bands, detailed, banded, True))
            out.append(("count", str(23 * case.ploidy)))
            out.append(("sex", case.sex_chromosomes or "U"))
        else:
            out.extend(_karyo_line_facts(case, bands, detailed, banded, banded))
    return _index_facts(out)


def _sex_remainder(case: Case, events: Optional[Sequence[Event]] = None) -> str:
    """Normal sex chromosome named alongside a structurally abnormal one
    (4.2.1b, 11.4.2.1/ix)."""
    sx = case.sex_chromosomes or ""
    for e in (case.events if events is None else events):
        if e.chrom in ("X", "Y") and e.scale != "chromosome" and e.structure:
            sx = sx.replace(e.chrom, "", 1)
    return sx


def _effective_sex(case: Case, events: Optional[Sequence[Event]] = None) -> str:
    evs = case.events if events is None else events
    sx = case.sex_chromosomes or "U"
    if any(e.chrom in ("X", "Y") and e.scale != "chromosome" and e.structure
           for e in evs):
        return _sex_remainder(case, evs) or "U"
    return sx


def roundtrip_check(case: Case, form: str, bands: CytobandTable) -> Dict[str, Any]:
    """Render, parse back with the grammar, compare with the projection."""
    res: Dict[str, Any] = {"form": form, "text": None, "status": None,
                           "detail": "", "derivations": None}
    try:
        text = render(case, form, bands)
    except Unrenderable as exc:
        res["status"] = "not_admissible"
        res["detail"] = exc.reason
        return res
    res["text"] = text
    try:
        parsed = parse(text)
    except OutOfScope as exc:
        res["status"] = "out_of_scope"
        res["detail"] = str(exc)
        return res
    except ParseError as exc:
        res["status"] = "parse_failed"
        res["detail"] = str(exc)
        return res
    expected = project(case, form, bands)
    got = parsed["facts"]
    if got == expected:
        res["status"] = "ok"
    else:
        res["status"] = "mismatch"
        res["detail"] = json.dumps(
            {"missing": sorted(map(str, expected - got)),
             "extra": sorted(map(str, got - expected))}, ensure_ascii=False)
    try:
        res["derivations"] = count_derivations(text)
    except Exception as exc:                                   # pragma: no cover
        res["derivations"] = f"error: {exc}"
    return res


# ---------------------------------------------------------------------------
# Ambiguity
# ---------------------------------------------------------------------------


def ambiguity_report(case: Case, bands: CytobandTable) -> Dict[str, Any]:
    """Per-form: which attributes of THIS case the form drops, and whether the
    dropped attributes are clinically load-bearing."""
    present = {
        "coordinates": any(e.start is not None for e in case.events),
        "flanking_normal": any(e.flank_left is not None for e in case.events),
        "mosaic_fraction": any(e.mosaic_fraction is not None
                               or e.mosaic_fraction_range
                               or e.mosaic_fraction_unknown for e in case.events),
        "mechanism": any(e.structure for e in case.events),
        "parental_origin": any(e.inheritance for e in case.events),
        "zygosity": any(e.zygosity for e in case.events),
        "sex_complement": case.sex_chromosomes not in (None, "U"),
        "copy_number": True,
        "breakpoint_bands": any(e.start is not None for e in case.events),
        "segment_orientation": any((e.structure or {}).get("detailed")
                                   for e in case.events),
    }
    out: Dict[str, Any] = {}
    for key in FORMS:
        keep = LOSS_MATRIX[key]
        lost = [a for a, p in present.items() if p and not keep.get(a, False)]
        out[key] = {"lost_attributes": lost}
    return out


def collision_scan(cases: Sequence[Tuple[str, Case]], bands: CytobandTable
                   ) -> Dict[str, Any]:
    """Measure ambiguity: clinically different cases that a form prints alike.

    Only cases the form can actually express are counted, so the rate is the
    probability that two admissible results become indistinguishable.
    """
    report: Dict[str, Any] = {}
    for key in FORMS:
        rendered: Dict[str, List[str]] = {}
        for name, case in cases:
            try:
                txt = render(case, key, bands)
            except Unrenderable:
                continue
            rendered.setdefault(txt, []).append(name)
        n = sum(len(v) for v in rendered.values())
        coll = {k: v for k, v in rendered.items() if len(v) > 1}
        report[key] = {
            "admissible_cases": n,
            "distinct_strings": len(rendered),
            "colliding_groups": len(coll),
            "cases_in_collisions": sum(len(v) for v in coll.values()),
            "collision_rate": round(sum(len(v) for v in coll.values()) / n, 4)
            if n else 0.0,
            "examples": [{"text": k, "cases": v} for k, v in list(coll.items())[:6]],
        }
    return report


# ---------------------------------------------------------------------------
# JSON input/output
# ---------------------------------------------------------------------------


def _event_from_dict(ed: Dict[str, Any]) -> Event:
    ed = {k: v for k, v in ed.items() if v is not None or k in ("start", "end")}
    for k in ("copy_number_range", "mosaic_fraction_range"):
        if ed.get(k):
            ed[k] = tuple(ed[k])
    return Event(**ed)


def case_from_dict(d: Dict[str, Any]) -> Case:
    """Rebuild a Case from plain JSON, clones included.

    The inverse of ``case_to_dict``: the pair has to round-trip, otherwise a
    case handed over as JSON (command line, web stand) would silently lose its
    clone description and be rendered from dosage alone.
    """
    events = [_event_from_dict(ed) for ed in d.get("events", [])]
    clones = []
    for cd in d.get("clones", []):
        cd = dict(cd)
        cd["events"] = [_event_from_dict(ed) for ed in cd.get("events", [])]
        clones.append(Clone(**cd))
    kw = {k: v for k, v in d.items() if k not in ("events", "clones")}
    return Case(events=events, clones=clones, **kw)


def case_to_dict(case: Case) -> Dict[str, Any]:
    return dataclasses.asdict(case)


def analyse(case: Case, bands: CytobandTable) -> Dict[str, Any]:
    """Full result for one case: forms, reverse validation, ambiguity."""
    problems = validate_case(case, bands)
    out: Dict[str, Any] = {
        "sample": case.sample,
        "input_problems": problems,
        "forms": {},
        "ambiguity": ambiguity_report(case, bands),
        "rules_applied": sorted({f.citation for f in FORMS.values()}),
    }
    for key, form in FORMS.items():
        rt = roundtrip_check(case, key, bands)
        out["forms"][key] = {
            "name_ru": form.name_ru,
            "format": form.fmt,
            "system": form.system,
            "text": rt["text"],
            "admissible": rt["status"] not in ("not_admissible",),
            "roundtrip": rt["status"],
            "detail": rt["detail"],
            "derivations": rt["derivations"],
            "lost_attributes": out["ambiguity"][key]["lost_attributes"],
            "warnings": conformance_warnings(case, key),
        }
    canonical = "short" if any(e.start is not None for e in case.events) else "abbrev"
    if case.normal():
        canonical = "abbrev"
    out["canonical_form"] = canonical
    out["karyotype_omissions"] = karyotype_omissions(case)
    out["canonical_text"] = out["forms"][canonical]["text"]
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _default_bands() -> CytobandTable:
    for cand in ("cytoBand_hg38.txt.gz", "cytoBand_hg38.txt",
                 os.path.join(os.path.dirname(__file__), "cytoBand_hg38.txt.gz")):
        if os.path.exists(cand):
            return CytobandTable(cand)
    raise SystemExit("таблица цитобендов не найдена: укажите --cytobands")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Генератор, разборщик и анализ однозначности формул ISCN 2024")
    ap.add_argument("command", choices=["render", "parse", "check", "rules"])
    ap.add_argument("--input", help="JSON со случаем или файл со строками формул")
    ap.add_argument("--text", help="строка формулы для parse/check")
    ap.add_argument("--form", default="all")
    ap.add_argument("--cytobands")
    ap.add_argument("--out")
    args = ap.parse_args(argv)

    if args.command == "rules":
        for k in sorted(RULES):
            print(f"{k}\t{RULES[k]}")
        return 0

    bands = CytobandTable(args.cytobands) if args.cytobands else _default_bands()

    if args.command == "render":
        with open(args.input) as fh:
            data = json.load(fh)
        cases = data if isinstance(data, list) else [data]
        results = [analyse(case_from_dict(c), bands) for c in cases]
        text = json.dumps(results, ensure_ascii=False, indent=2)
        if args.out:
            with open(args.out, "w") as fh:
                fh.write(text)
        else:
            print(text)
        return 0

    texts: List[str] = []
    if args.text:
        texts = [args.text]
    elif args.input:
        with open(args.input) as fh:
            texts = [l.strip() for l in fh if l.strip()]
    for t in texts:
        try:
            p = parse(t)
            row: Dict[str, Any] = {"text": t, "status": "ok",
                                   "derivations": count_derivations(t)}
            if args.command == "check":
                row["violations"] = check_conventions(t, bands)
                row["bands"] = check_band_consistency(t, bands)
            else:
                row["facts"] = sorted(map(str, p["facts"]))
                row["meta"] = p["meta"]
            print(json.dumps(row, ensure_ascii=False))
        except OutOfScope as exc:
            print(json.dumps({"text": t, "status": "out_of_scope",
                              "detail": str(exc)}, ensure_ascii=False))
        except ParseError as exc:
            print(json.dumps({"text": t, "status": "invalid",
                              "detail": str(exc)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
