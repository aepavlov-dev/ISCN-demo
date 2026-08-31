"""Case battery for the ISCN formula generator.

Three groups:

``REAL``      the three worked samples of the run (low-pass WGS, ~12x);
``BOOK``      inputs whose expected output is a string published in ISCN 2024,
              so the generator is checked against the source verbatim;
``SYNTH``     PGT-A and PGT-SR situations covering the practice range.

Coordinates are 1-based inclusive.  ``from_bed`` is used where the value came
from the pipeline's bin table (0-based half-open).
"""

from iscn_formula import Case, Clone, Event, from_bed

# ---------------------------------------------------------------------------
# The three worked samples
# ---------------------------------------------------------------------------

REAL = [
    ("NA00496", "Мозаичная трисомия 8, доля 0,951 (мужской образец)",
     Case(sample="NA00496", technique="seq", sex_chromosomes="XY",
          events=[Event(chrom="8", copy_number=3, scale="chromosome",
                        mosaic_fraction=0.95,
                        label="chr8 500 000-145 000 000, CN=3, копийность "
                              "2,951 (2,946-2,957), доля клеток 0,951 "
                              "(0,946-0,957) — из pgt_report_NA00496_filled.md")])),

    ("NA02325-доза", "Три прироста дозы, механизм не установлен",
     Case(sample="NA02325", technique="seq", sex_chromosomes="XX",
          events=[from_bed("X", 154890327, 155335092, copy_number=3),
                  from_bed("16", 200000, 3700000, copy_number=3),
                  from_bed("22", 17100000, 30500000, copy_number=3)])),

    ("NA02325-производная", "Тот же образец с установленным механизмом: "
                            "сверхчисленная der(22)t(16;22) материнского происхождения",
     Case(sample="NA02325", technique="seq", sex_chromosomes="XX",
          karyotype_available=True,
          events=[from_bed("X", 154890327, 155335092, copy_number=3),
                  from_bed("16", 200000, 3700000, copy_number=3,
                           explained_by="22"),
                  from_bed("22", 17100000, 30500000, copy_number=3,
                           supernumerary=True,
                           structure={"chain": [
                               {"symbol": "der", "chroms": ["22"]},
                               {"symbol": "t", "chroms": ["16", "22"],
                                "breakpoints": ["p13.3", "q11.21"]}],
                               "tail": "mat"})])),

    ("NA13019-доза", "Одна X: потеря Xp и избыток Xq, доля 0,68",
     Case(sample="NA13019", technique="seq", sex_chromosomes="X",
          events=[from_bed("X", 2900000, 56900000, copy_number=1),
                  from_bed("X", 62800000, 155700000, copy_number=3,
                           mosaic_fraction=0.68)])),

    ("NA13019-клоны", "Тот же образец клонами: изохромосома Xq в 68 % и "
                      "моносомия X в 32 %",
     Case(sample="NA13019", technique="seq", sex_chromosomes="X",
          events=[from_bed("X", 2900000, 56900000, copy_number=1),
                  from_bed("X", 62800000, 155700000, copy_number=3,
                           mosaic_fraction=0.68)],
          clones=[Clone(events=[Event(chrom="X", copy_number=2, baseline=1,
                                      start=62800001, end=155700000,
                                      structure={"chain": [
                                          {"symbol": "i", "chroms": ["X"],
                                           "breakpoints": ["q10"]}]})],
                        sex_chromosomes="XX", fraction=0.68),
                  Clone(events=[Event(chrom="X", copy_number=1,
                                      scale="chromosome")],
                        sex_chromosomes="X", fraction=0.32)])),
]

# ---------------------------------------------------------------------------
# Inputs whose rendering must equal a string published in the book
# ---------------------------------------------------------------------------

BOOK = [
    ("8.2.1b", "abbrev", "arr (X,1\u201322)\u00d72",
     Case(sample="b1", technique="arr", sex_chromosomes="XX")),
    ("8.2.1b-male", "abbrev", "arr (X,Y)\u00d71,(1\u201322)\u00d72",
     Case(sample="b2", technique="arr", sex_chromosomes="XY")),
    ("8.2.6i", "abbrev", "arr (X,1\u201322)\u00d72hmz",
     Case(sample="b3", technique="arr", sex_chromosomes="XX",
          genome_zygosity="hmz")),
    ("8.2.6ii", "abbrev", "arr (X,1\u201322)\u00d72hmz htz",
     Case(sample="b4", technique="arr", sex_chromosomes="XX",
          genome_zygosity="hmz htz")),
    ("8.2.7iii", "abbrev", "arr (X)\u00d72,(Y)\u00d71,(1\u201322)\u00d73",
     Case(sample="b5", technique="arr", sex_chromosomes="XXY", ploidy=3)),
    ("8.2.7iv", "abbrev", "arr (X,1\u201322)\u00d73",
     Case(sample="b6", technique="arr", sex_chromosomes="XXX", ploidy=3)),
    ("8.2.2xxi", "short",
     "arr[GRCh38] 6q21q25.1(113,900,000_149,100,000)\u00d71,(21)\u00d73",
     Case(sample="b7", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="6", start=113900000, end=149100000,
                        copy_number=1),
                  Event(chrom="21", copy_number=3, scale="chromosome")])),
    ("8.2.2xxii", "short",
     "arr[GRCh38] 9p24.3p13.1(204,166_38,756,057)\u00d71,"
     "18q21.33q22.1(63,877,984_64,683,663)\u00d71,"
     "21q11.2q21.1(13,600,026_20,175,986)\u00d73",
     Case(sample="b8", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="9", start=204166, end=38756057, copy_number=1),
                  Event(chrom="18", start=63877984, end=64683663, copy_number=1),
                  Event(chrom="21", start=13600026, end=20175986,
                        copy_number=3)])),
    ("8.2.2xxiii", "short",
     "arr[GRCh38] 14q31.1(82,695,844_82,855,387)\u00d71,"
     "14q32.33(105,643,093_106,109,395)\u00d73",
     Case(sample="b9", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="14", start=82695844, end=82855387, copy_number=1),
                  Event(chrom="14", start=105643093, end=106109395,
                        copy_number=3)])),
    # Книга печатает 18q22.3, что соответствует hg19; при заявленной сборке
    # GRCh38 позиция 69 172 132 лежит в 18q22.2. Ожидаемое значение исправлено
    # по таблице цитобендов GRCh38, расхождение с книгой разобрано в отчёте.
    ("8.2.2xxiv", "short",
     "arr[GRCh38] 18p11.32p11.21(102,328_15,079,388)\u00d71,"
     "18q22.2q23(69,172,132_79,093,443)\u00d71",
     Case(sample="b10", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="18", start=102328, end=15079388, copy_number=1),
                  Event(chrom="18", start=69172132, end=79093443,
                        copy_number=1)])),
    ("8.2.2xix", "extended",
     "arr[GRCh38] 4q32.2q35.1(163002425\u00d72,163146681_183022312\u00d71,"
     "184322231\u00d72)",
     Case(sample="b11", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="4", start=163146681, end=183022312,
                        copy_number=1, flank_left=163002425,
                        flank_right=184322231)])),
    ("8.2.2xx", "extended",
     "arr[GRCh38] 11p12(37003221\u00d72,37741458_39209912\u00d73,39752007\u00d72)",
     Case(sample="b12", technique="arr", sex_chromosomes="XX",
          events=[Event(chrom="11", start=37741458, end=39209912,
                        copy_number=3, flank_left=37003221,
                        flank_right=39752007)])),
    ("8.2.3i", "short",
     "arr[GRCh38] Xq25(126,228,413_126,535,347)\u00d70mat",
     Case(sample="b13", technique="arr", sex_chromosomes="XY",
          events=[Event(chrom="X", start=126228413, end=126535347,
                        copy_number=0, baseline=1, inheritance="mat")])),
    ("11.4.1ix", "abbrev", "seq (X)\u00d71[0.6]",
     Case(sample="b14", technique="seq", sex_chromosomes="XX",
          events=[Event(chrom="X", copy_number=1, scale="chromosome",
                        mosaic_fraction=0.6)])),
    ("5.2iii", "karyotype", "46,U",
     Case(sample="b15", sex_chromosomes="U", karyotype_available=True)),
    ("5.3.1.1i", "karyotype", "45,X",
     Case(sample="b16", sex_chromosomes="X", karyotype_available=True,
          events=[Event(chrom="X", copy_number=1, scale="chromosome")])),
    ("5.4.1-iso", "karyotype", "46,X,i(X)(q10)",
     Case(sample="b17", sex_chromosomes="XX", karyotype_available=True,
          events=[Event(chrom="X", copy_number=2, baseline=1,
                        start=62800001, end=155700000,
                        structure={"chain": [{"symbol": "i", "chroms": ["X"],
                                              "breakpoints": ["q10"]}]})])),
    ("5.5.3-chain", "karyotype", "46,XY,der(1)t(1;3)(p32;q21)dup(1)(q25q42)",
     Case(sample="b18", sex_chromosomes="XY", karyotype_available=True,
          events=[Event(chrom="1", copy_number=3, start=1, end=248956422,
                        structure={"chain": [
                            {"symbol": "der", "chroms": ["1"]},
                            {"symbol": "t", "chroms": ["1", "3"],
                             "breakpoints": ["p32", "q21"]},
                            {"symbol": "dup", "chroms": ["1"],
                             "breakpoints": ["q25", "q42"]}]})])),
    ("8.2.4-rec", "karyotype", "46,XY,rec(18)dup(18q)inv(18)(p11.32q21)dpat",
     Case(sample="b19", sex_chromosomes="XY", karyotype_available=True,
          events=[Event(chrom="18", copy_number=3, start=1, end=80373285,
                        structure={"chain": [
                            {"symbol": "rec", "chroms": ["18"]},
                            {"symbol": "dup", "chroms": ["18q"],
                             "breakpoints": []},
                            {"symbol": "inv", "chroms": ["18"],
                             "breakpoints": ["p11.32", "q21"]}],
                            "tail": "dpat"})])),
]

# ---------------------------------------------------------------------------
# Synthetic PGT-A / PGT-SR battery
# ---------------------------------------------------------------------------


def _seg(chrom, start, end, cn, **kw):
    return Event(chrom=chrom, start=start, end=end, copy_number=cn, **kw)


def _whole(chrom, cn, **kw):
    return Event(chrom=chrom, copy_number=cn, scale="chromosome", **kw)


SYNTH = [
    # --- PGT-A: euploid and whole-chromosome aneuploidy -------------------
    ("A01", "Эуплоидный женский профиль",
     Case(sample="A01", sex_chromosomes="XX")),
    ("A02", "Эуплоидный мужской профиль",
     Case(sample="A02", sex_chromosomes="XY")),
    ("A03", "Эуплоидный профиль, пол не раскрывается",
     Case(sample="A03", sex_chromosomes="U", sex_disclosed=False)),
    ("A04", "Полная трисомия 21",
     Case(sample="A04", sex_chromosomes="XY", events=[_whole("21", 3)])),
    ("A05", "Полная моносомия 21",
     Case(sample="A05", sex_chromosomes="XX", events=[_whole("21", 1)])),
    ("A06", "Полная трисомия 16",
     Case(sample="A06", sex_chromosomes="XX", events=[_whole("16", 3)])),
    ("A07", "Двойная анеуплоидия: трисомия 21 и моносомия 22",
     Case(sample="A07", sex_chromosomes="XY",
          events=[_whole("21", 3), _whole("22", 1)])),
    ("A08", "Хаотический набор: пять полных событий",
     Case(sample="A08", sex_chromosomes="XX",
          events=[_whole("2", 3), _whole("7", 1), _whole("13", 3),
                  _whole("18", 1), _whole("21", 3)])),
    # --- PGT-A: mosaicism -------------------------------------------------
    ("A09", "Мозаичная трисомия 16, доля 0,35",
     Case(sample="A09", sex_chromosomes="XX",
          events=[_whole("16", 3, mosaic_fraction=0.35)])),
    ("A10", "Мозаичная трисомия 16 у нижнего предела, доля 0,20",
     Case(sample="A10", sex_chromosomes="XX",
          events=[_whole("16", 3, mosaic_fraction=0.20)])),
    ("A11", "Мозаичная моносомия 15, доля не определена",
     Case(sample="A11", sex_chromosomes="XY",
          events=[_whole("15", 1, mosaic_fraction_unknown=True)])),
    ("A12", "Мозаичная трисомия 8, доля дана интервалом 0,3~0,5",
     Case(sample="A12", sex_chromosomes="XY",
          events=[_whole("8", 3, mosaic_fraction_range=(0.3, 0.5))])),
    ("A13", "Мозаичная трисомия 21, счёт клеток известен (кариотип выполнен)",
     Case(sample="A13", sex_chromosomes="XY", karyotype_available=True,
          cell_counts={"abnormal": 12, "normal": 18},
          events=[_whole("21", 3, mosaic_fraction=0.4)])),
    # --- PGT-A: sex chromosomes ------------------------------------------
    ("A14", "Моносомия X",
     Case(sample="A14", sex_chromosomes="X", karyotype_available=True,
          events=[_whole("X", 1)])),
    ("A15", "47,XXY",
     Case(sample="A15", sex_chromosomes="XXY", karyotype_available=True,
          events=[_whole("X", 2, baseline=1)])),
    ("A16", "47,XXX",
     Case(sample="A16", sex_chromosomes="XXX", karyotype_available=True,
          events=[_whole("X", 3)])),
    ("A17", "47,XYY",
     Case(sample="A17", sex_chromosomes="XYY", karyotype_available=True,
          events=[_whole("Y", 2, baseline=1)])),
    ("A18", "Мозаичная моносомия X, доля 0,5",
     Case(sample="A18", sex_chromosomes="X",
          events=[_whole("X", 1, mosaic_fraction=0.5)])),
    # --- PGT-A: ploidy ----------------------------------------------------
    ("A19", "Триплоидия 69,XXX",
     Case(sample="A19", sex_chromosomes="XXX", ploidy=3,
          karyotype_available=True)),
    ("A20", "Триплоидия 69,XXY",
     Case(sample="A20", sex_chromosomes="XXY", ploidy=3,
          karyotype_available=True)),
    ("A21", "Гаплоидия 23,X",
     Case(sample="A21", sex_chromosomes="X", ploidy=1,
          karyotype_available=True)),
    ("A22", "Тетраплоидия 92,XXYY",
     Case(sample="A22", sex_chromosomes="XXYY", ploidy=4,
          karyotype_available=True)),
    ("A23", "Полный пузырный занос: диплоидия при полной гомозиготности",
     Case(sample="A23", sex_chromosomes="XX", genome_zygosity="hmz")),
    ("A24", "Двуспермный занос: гомозиготные и гетерозиготные участки",
     Case(sample="A24", sex_chromosomes="XX", genome_zygosity="hmz htz")),
    # --- PGT-A: segmental -------------------------------------------------
    ("A25", "Концевая делеция 4p16.3p15.32 (синдром Вольфа-Хиршхорна)",
     Case(sample="A25", sex_chromosomes="XX",
          events=[_seg("4", 68346, 17000000, 1)])),
    ("A26", "Внутренняя дупликация 15q11.2q13.1",
     Case(sample="A26", sex_chromosomes="XY",
          events=[_seg("15", 22600000, 28500000, 3)])),
    ("A27", "Мозаичная сегментная делеция 1p36, доля 0,45",
     Case(sample="A27", sex_chromosomes="XX",
          events=[_seg("1", 1000000, 10000000, 1, mosaic_fraction=0.45)])),
    ("A28", "Делеция с фланкирующими нормальными позициями (расширенная система)",
     Case(sample="A28", sex_chromosomes="XY",
          events=[_seg("22", 18730698, 21689521, 1,
                       flank_left=18650000, flank_right=21750000)])),
    ("A29", "Гомозиготный участок отцовского происхождения: "
            "однородительская дисомия 15",
     Case(sample="A29", sex_chromosomes="XY",
          events=[_seg("15", 22600000, 101991189, 2, zygosity="hmz",
                       inheritance="upat")])),
    ("A30", "Сегмент 3,5 Мб при малой глубине (sseq)",
     Case(sample="A30", technique="sseq", sex_chromosomes="XX",
          events=[_seg("7", 100000000, 103500000, 1)])),
    ("A31", "Первое полярное тельце: потеря хроматиды 21",
     Case(sample="A31", technique="arr", chromatid=True,
          sex_chromosomes=None,
          events=[_whole("21", 0, baseline=2)])),
    ("A32", "Комплексная перестройка одной хромосомы (хромотрипсис)",
     Case(sample="A32", sex_chromosomes="XY",
          events=[_seg("5", 1000000, 20000000, 1),
                  _seg("5", 20000001, 30000000, 3),
                  _seg("5", 40000000, 50000000, 1)])),
    # --- PGT-SR -----------------------------------------------------------
    ("S01", "Несбалансированный продукт реципрокной транслокации, "
            "сегрегация 2:2: потеря 5p, прирост 11q",
     Case(sample="S01", sex_chromosomes="XX", karyotype_available=True,
          events=[_seg("5", 1, 22000000, 1, explained_by="5"),
                  _seg("11", 100000000, 135086622, 3,
                       structure={"chain": [
                           {"symbol": "der", "chroms": ["5"]},
                           {"symbol": "t", "chroms": ["5", "11"],
                            "breakpoints": ["p14.1", "q22.3"]}],
                           "tail": "mat"})])),
    ("S02", "Тертиарная трисомия, сегрегация 3:1: сверхчисленная der(22)t(11;22)",
     Case(sample="S02", sex_chromosomes="XY", karyotype_available=True,
          events=[_seg("11", 116000000, 135086622, 3, explained_by="22"),
                  _seg("22", 1, 21000000, 3, supernumerary=True,
                       structure={"chain": [
                           {"symbol": "der", "chroms": ["22"]},
                           {"symbol": "t", "chroms": ["11", "22"],
                            "breakpoints": ["q23.3", "q11.21"]}],
                           "tail": "pat"})])),
    ("S03", "Несбалансированный продукт робертсоновской транслокации: "
            "трисомия 14 при der(13;14) у носителя",
     Case(sample="S03", sex_chromosomes="XX", karyotype_available=True,
          events=[_whole("14", 3),
                  _seg("13", 1, 114364328, 2, explained_by="14",
                       structure={"chain": [
                           {"symbol": "der", "chroms": ["13", "14"],
                            "breakpoints": [["q10"], ["q10"]]}]})])),
    ("S04", "Рекомбинантная хромосома от носителя инверсии",
     Case(sample="S04", sex_chromosomes="XY", karyotype_available=True,
          events=[_seg("18", 1, 80373285, 3,
                       structure={"chain": [
                           {"symbol": "rec", "chroms": ["18"]},
                           {"symbol": "dup", "chroms": ["18q"],
                            "breakpoints": []},
                           {"symbol": "inv", "chroms": ["18"],
                            "breakpoints": ["p11.32", "q21"]}],
                           "tail": "dpat"})])),
    ("S05", "Кольцевая хромосома 14",
     Case(sample="S05", sex_chromosomes="XY", karyotype_available=True,
          events=[_seg("14", 1, 107043718, 2,
                       structure={"chain": [
                           {"symbol": "r", "chroms": ["14"],
                            "breakpoints": ["p13", "q32"]}]})])),
    ("S06", "Изохромосома Xq у единственной X",
     Case(sample="S06", sex_chromosomes="XX", karyotype_available=True,
          events=[Event(chrom="X", copy_number=2, baseline=1,
                        start=62800001, end=155700000,
                        structure={"chain": [{"symbol": "i", "chroms": ["X"],
                                              "breakpoints": ["q10"]}]})])),
    ("S07", "Производная хромосома с тремя перестройками",
     Case(sample="S07", sex_chromosomes="XY", karyotype_available=True,
          events=[_seg("1", 1, 248956422, 3,
                       structure={"chain": [
                           {"symbol": "der", "chroms": ["1"]},
                           {"symbol": "t", "chroms": ["1", "3"],
                            "breakpoints": ["p32", "q21"]},
                           {"symbol": "dup", "chroms": ["1"],
                            "breakpoints": ["q25", "q42"]}]})])),
    ("S08", "Сегментное событие с установленной дозой, но неизвестным механизмом",
     Case(sample="S08", sex_chromosomes="XX",
          events=[_seg("7", 1, 25000000, 1)])),
]

ALL = ([(k, d, c, "образец") for k, d, c in REAL]
       + [(k, d, c, "синтетика") for k, d, c in SYNTH])
