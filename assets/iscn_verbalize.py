"""Словесное описание результата по формуле ISCN — для заключения врача.

Задача. Формула ISCN однозначна, но нечитаема без подготовки, поэтому в
заключении она сопровождается словесным описанием. Описание не пишется от руки:
оно порождается из того же набора событий, из которого построена формула, и
проверяется обратным извлечением фактов из готового текста. Тем самым
расхождение между формулой, таблицей находок и словами технически невозможно.

Четыре слоя описания (порядок соответствует практике европейских рекомендаций
по оформлению заключений, Claustres et al., EJHG 2014;22:160-170):

1. ``дословная расшифровка`` — что именно записано в формуле: хромосома, тип
   события, границы, размер, копийность, доля клеток. Слой механический:
   каждый его факт извлекаем из формулы и обратно сверяем с ней.
2. ``клиническая значимость`` — категория находки по пяти категориям
   упомянутых рекомендаций и, для полных анеуплоидий, общепринятое название
   состояния. Для сегментных событий вместо названия синдрома приводится
   перекрытие критического участка с долей перекрытия: это утверждение о
   координатах, а не диагноз.
3. ``ограничения`` — предел обнаружения, чем результат не является, оговорка о
   соответствии биопсии трофэктодермы внутренней клеточной массе.
4. ``что дальше`` — возможности, но не предписания: рекомендации прямо
   запрещают писать, что пренатальная диагностика «показана» или «необходима»;
   лаборатория вправе лишь сообщить, что она возможна.

Настройки. Класс ``Style`` задаёт стиль и состав описания; готовые профили в
``PROFILES`` собраны под роли получателей (репродуктолог, генетик, пациент,
досье). Профиль — это данные, а не код: врач меняет поля, не трогая шаблоны.
"""

from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import iscn_formula as IF

__all__ = ["Style", "PROFILES", "verbalize", "check_text", "SYNDROMES",
           "MOSAIC_LABELS", "MOSAIC_LABELS_PATIENT", "mosaic_bands", "load_regions", "describe_event"]


# ---------------------------------------------------------------------------
# Настройки стиля
# ---------------------------------------------------------------------------

@dataclass
class Style:
    """Настройки словесного описания.

    Поля разделены на три группы: состав (какие слои включены), терминология
    (какими словами называть события) и оформление (числа, лицо, длина).
    """

    name: str
    audience: str = "генетик"

    # --- состав -----------------------------------------------------------
    include_formula: bool = True
    formula_form: str = "short"          # какая форма ISCN приводится в тексте
    literal_layer: bool = True           # слой 1: дословная расшифровка
    significance_layer: bool = True      # слой 2: клиническая значимость
    limits_layer: bool = True            # слой 3: ограничения
    next_steps_layer: bool = True        # слой 4: что дальше
    counselling_note: bool = True        # напоминание о медико-генетическом
                                         # консультировании

    # --- терминология -----------------------------------------------------
    term_register: str = "профессиональный"   # профессиональный | смешанный | пациентский
    name_syndromes: bool = True          # называть состояние при полной анеуплоидии
    region_overlap: bool = True          # приводить перекрытие критического участка
    mosaic_policy: str = "класс_и_доля"  # класс_и_доля | только_доля | только_класс | не_указывать
    #: Границы полосы заявления мозаицизма: нижняя, разделяющая низкий и
    #: высокий уровень, верхняя. Значения по умолчанию (0,20; 0,50; 0,80) —
    #: внутренняя граница лаборатории, а не величина из первоисточника:
    #: приписывать их позиционным документам нельзя — они не получены. Опора —
    #: ESHRE hoaa017: решение о заявлении мозаицизма принимает сам центр по
    #: собственной валидации. Нижняя граница отвечает полосе
    #: шума платформы, верхняя — доле, выше которой находка считается полной.
    #: Лаборатория обязана подставить свои, полученные на своём наборе и своей
    #: платформе. Названия классов строятся из этих чисел, а не из подписей.
    mosaic_thresholds: Tuple[float, ...] = (0.20, 0.50, 0.80)
    #: Выделять вывод отдельным заголовком (ACGS, п. 2.1: «The overall result or
    #: conclusion must be clearly visible»).
    highlight_conclusion: bool = True
    #: Повторять клинический вопрос исследования (ACGS, п. 2.5.2: «reports must
    #: explicitly restate the clinical question being asked»). Текст вопроса
    #: подставляет вызывающая сторона: сам инструмент показания не знает.
    restate_question: bool = True
    #: Печатать указание о пригодности к переносу (ESHRE hoaa021 требует его в
    #: заключении цикла ПГТ). Указание НЕ выводится из находки: его формулирует
    #: врач, инструмент лишь переносит переданный текст.
    transfer_guidance: bool = False
    mosaic_uncertainty: bool = True      # приводить интервал доли
    sex_policy: str = "при_находке"      # не_раскрывать | при_находке | всегда
    size_units: str = "Мб"
    coordinates: bool = True             # приводить координаты и сборку

    # --- оформление -------------------------------------------------------
    person_form: str = "безличная"       # безличная | лабораторная
    decimal_comma: bool = True
    percent_for_fraction: bool = True    # доля клеток и в процентах
    max_sentences: Optional[int] = None  # ограничение длины слоёв 1-2
    bullet_layers: bool = False          # слои списком, а не абзацами
    heading_prefix: str = ""

    def with_(self, **kw: Any) -> "Style":
        """Копия профиля с изменёнными полями: врач правит настройку, не код."""
        return dataclasses.replace(self, **kw)


#: Готовые профили под роли получателей заключения.
PROFILES: Dict[str, Style] = {
    "репродуктолог": Style(
        name="репродуктолог", audience="репродуктолог",
        include_formula=True, formula_form="short",
        literal_layer=True, significance_layer=True, limits_layer=True,
        next_steps_layer=True, counselling_note=True,
        term_register="профессиональный", name_syndromes=True,
        region_overlap=True, mosaic_policy="класс_и_доля",
        sex_policy="не_раскрывать", coordinates=False,
        max_sentences=2, bullet_layers=True),
    "генетик": Style(
        name="генетик", audience="генетик",
        include_formula=True, formula_form="short",
        term_register="профессиональный", name_syndromes=True,
        region_overlap=True, mosaic_policy="класс_и_доля",
        mosaic_uncertainty=True, sex_policy="при_находке",
        coordinates=True, max_sentences=None, bullet_layers=False),
    "пациент": Style(
        name="пациент", audience="пациент",
        include_formula=False, literal_layer=True, significance_layer=True,
        limits_layer=True, next_steps_layer=True, counselling_note=True,
        term_register="пациентский", name_syndromes=True,
        region_overlap=False, mosaic_policy="только_класс",
        mosaic_uncertainty=False, sex_policy="не_раскрывать",
        coordinates=False, max_sentences=2, bullet_layers=True),
    "досье": Style(
        name="досье", audience="досье",
        include_formula=True, formula_form="short",
        term_register="профессиональный", name_syndromes=False,
        region_overlap=True, mosaic_policy="класс_и_доля",
        mosaic_uncertainty=True, sex_policy="всегда",
        coordinates=True, next_steps_layer=False, counselling_note=False,
        max_sentences=None, bullet_layers=False),
}


# ---------------------------------------------------------------------------
# Словари
# ---------------------------------------------------------------------------

#: Подписи классов мозаицизма. Числа в подписи подставляются из настройки
#: ``Style.mosaic_thresholds``: если границы заданы иные, а подпись осталась
#: прежней, текст заключения утверждал бы не то значение, по которому принято
#: решение.
MOSAIC_LABELS = ("доля ниже границы заявления мозаицизма {lo}",
                 "мозаицизм низкого уровня",
                 "мозаицизм высокого уровня",
                 "доля выше границы {hi}")

#: Те же классы для пациентского регистра: числовая граница читателю ничего не
#: говорит, поэтому называется наблюдаемое положение дел.
MOSAIC_LABELS_PATIENT = ("изменение затрагивает малую часть клеток",
                         "изменение затрагивает меньшую часть клеток",
                         "изменение затрагивает большую часть клеток",
                         "изменение затрагивает почти все исследованные клетки")
#: Границы 0,20 и 0,80 в практике ПГТ выбраны по полосе шума платформ на
#: микрочипах и малой глубине прочтения, а не по биологии. Настоящий метод даёт
#: доверительный интервал доли, поэтому решение «событие во всех клетках»
#: принимается по интервалу (см. iscn_report_bridge), а класс по этим границам
#: приводится только как принятый в практике ориентир и «во всех клетках» не
#: утверждает.

#: Общепринятые названия состояний при полной анеуплоидии.
ANEUPLOIDY_NAMES = {
    ("21", "gain"): "трисомия 21 (синдром Дауна)",
    ("18", "gain"): "трисомия 18 (синдром Эдвардса)",
    ("13", "gain"): "трисомия 13 (синдром Патау)",
    ("16", "gain"): "трисомия 16",
    ("22", "gain"): "трисомия 22",
    ("X", "loss"): "моносомия X (синдром Шерешевского-Тёрнера)",
}

#: Критические участки, перекрытие которых указывается для сегментных событий.
#: Координаты GRCh38. Утверждение о перекрытии — это утверждение о координатах;
#: диагноз по нему не ставится, что и оговаривается в тексте.
SYNDROMES = {
    ("22", 18900000, 21500000, "loss"): "критический участок делеции 22q11.2",
    ("22", 18900000, 21500000, "gain"): "критический участок дупликации 22q11.2",
    ("15", 22800000, 28400000, "loss"): "критический участок 15q11-q13 "
                                        "(области синдромов Прадера-Вилли и Ангельмана)",
    ("5", 10000, 12500000, "loss"): "критический участок 5p (область синдрома "
                                    "Кри-дю-Ша)",
    ("4", 10000, 2000000, "loss"): "критический участок 4p16.3 (область синдрома "
                                   "Вольфа-Хиршхорна)",
    ("7", 73300000, 74900000, "loss"): "критический участок 7q11.23 (область "
                                       "синдрома Вильямса)",
    ("17", 16800000, 20200000, "loss"): "критический участок 17p11.2 (область "
                                        "синдрома Смита-Магениса)",
    ("1", 10000, 3000000, "loss"): "критический участок 1p36",
}

#: Плечо хромосомы в трёх падежах: «трисомия по длинному плечу», «лишняя копия
#: длинного плеча», «затронуто длинное плечо». Без падежных форм шаблон даёт
#: несогласованную фразу.
ARM_FORMS = {
    "p": {"nom": "короткое плечо", "gen": "короткого плеча", "dat": "короткому плечу"},
    "q": {"nom": "длинное плечо", "gen": "длинного плеча", "dat": "длинному плечу"},
}
ARM_WORDS = {k: f["nom"] for k, f in ARM_FORMS.items()}

#: Термины по регистру: профессиональный, смешанный (термин плюс пояснение),
#: пациентский (без латинизмов).
TERMS = {
    "профессиональный": {
        "gain_whole": "трисомия по хромосоме {chrom}",
        "loss_whole": "моносомия по хромосоме {chrom}",
        "gain_arm": "трисомия по {arm_dat} хромосомы {chrom} ({chrom}{arm})",
        "loss_arm": "моносомия по {arm_dat} хромосомы {chrom} ({chrom}{arm})",
        "gain_seg": "дупликация участка {band}",
        "loss_seg": "делеция участка {band}",
        "nullisomy": "гомозиготная делеция участка {band}",
        "copies": "копийность {cn}",
    },
    "смешанный": {
        "gain_whole": "трисомия — три копии хромосомы {chrom} вместо двух",
        "loss_whole": "моносомия — одна копия хромосомы {chrom} вместо двух",
        "gain_arm": "лишняя копия {arm_gen} хромосомы {chrom} ({chrom}{arm})",
        "loss_arm": "утрата одной копии {arm_gen} хромосомы {chrom} ({chrom}{arm})",
        "gain_seg": "дупликация — лишняя копия участка {band}",
        "loss_seg": "делеция — утрата участка {band}",
        "nullisomy": "утрата обеих копий участка {band}",
        "copies": "число копий {cn}",
    },
    "пациентский": {
        "gain_whole": "лишняя, третья копия хромосомы {chrom}",
        "loss_whole": "нехватка одной из двух копий хромосомы {chrom}",
        "gain_arm": "лишняя копия части хромосомы {chrom}",
        "loss_arm": "нехватка части хромосомы {chrom}",
        "gain_seg": "лишняя копия участка хромосомы {chrom}",
        "loss_seg": "нехватка участка хромосомы {chrom}",
        "nullisomy": "отсутствие участка хромосомы {chrom}",
        "copies": "копий {cn}",
    },
}

#: Категории значимости по рекомендациям Claustres et al., EJHG 2014.
#: Слова «положительный» и «отрицательный» в заключении не употребляются:
#: рекомендации прямо называют их двусмысленными.
CATEGORIES = {
    "норма": "находки в пределах физиологической изменчивости",
    "патогномоничная": "находка, однозначно связанная с клинически значимым состоянием",
    "неопределённая": "находка неопределённой значимости",
    "сопутствующая": "сопутствующая находка, не относящаяся к поставленному вопросу",
    "неспецифическая": "неспецифическая находка без клинической значимости",
}


#: Пациентский регистр: те же утверждения без латинизмов и без слова
#: «находка». Регистр меняет слова, но не содержание: категория, доля клеток и
#: предел обнаружения остаются теми же величинами.
CATEGORIES_PATIENT = {
    "норма": "изменений, которые метод способен увидеть, не обнаружено",
    "патогномоничная": "обнаруженное изменение известно и связано с заболеванием",
    "неопределённая": "значение обнаруженного изменения на сегодня неизвестно",
    "сопутствующая": "изменение обнаружено попутно и к заданному вопросу "
                     "не относится",
    "неспецифическая": "изменение не связано с известными заболеваниями",
}

MOSAIC_PATIENT = ("часть клеток исследованного образца несёт изменение, "
                  "часть — нет")

LIMITS_PATIENT = {
    "предел": "изменения, затрагивающие менее {dl} клеток образца, метод мог "
              "не заметить",
    "трофэктодерма": "исследованы клетки наружного слоя, из которого "
                     "развивается плацента; клетки будущего плода исследовать "
                     "без вреда для эмбриона нельзя",
    "механизм": "метод показывает, сколько копий участка есть, но не показывает, "
                "как именно перестроена хромосома",
    "половые": "для X и Y применяется отдельная модель расчёта, её ограничения "
               "описаны в отчёте",
}


# ---------------------------------------------------------------------------
# Числа и согласование
# ---------------------------------------------------------------------------

def _num(x: float, digits: int = 2, comma: bool = True) -> str:
    """Число для текста заключения.

    Незначащие нули снимаются только после десятичной точки: без этой оговорки
    ``100`` превращается в ``1``, а ``80`` в ``8``.
    """
    s = f"{x:.{digits}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s.replace(".", ",") if comma else s


def _mb(bp: int, style: Style) -> str:
    return f"{_num(bp / 1e6, 1, style.decimal_comma)} {style.size_units}"


def _thousands(n: int) -> str:
    return f"{n:,}".replace(",", "\u00a0")


def _fraction_text(f: Optional[float], style: Style,
                   lo: Optional[float] = None,
                   hi: Optional[float] = None) -> str:
    if f is None:
        return ""
    parts = [_num(f, 2, style.decimal_comma)]
    if style.percent_for_fraction:
        parts.append(f"({_num(100 * f, 0, style.decimal_comma)} %)")
    txt = " ".join(parts)
    if style.mosaic_uncertainty and lo is not None and hi is not None:
        txt += (f", доверительный интервал {_num(lo, 2, style.decimal_comma)}"
                f"-{_num(hi, 2, style.decimal_comma)}")
    return txt


def mosaic_bands(style: Style) -> Tuple[Tuple[float, float, str], ...]:
    """Полосы долей мозаицизма с подписями по настроенным границам.

    Принимаются две границы (нижняя и верхняя) или три (с разделением низкого
    и высокого уровня). Границы обязаны возрастать и лежать в пределах от нуля
    до единицы: неупорядоченные границы дали бы полосы, в которые не попадает
    ни одно значение, и класс молча исчез бы из текста.
    """
    th = tuple(float(x) for x in style.mosaic_thresholds)
    if len(th) == 2:
        lo, hi = th
        mid = round((lo + hi) / 2, 4)
    elif len(th) == 3:
        lo, mid, hi = th
    else:
        raise ValueError("mosaic_thresholds: ожидаются две или три границы, "
                         f"получено {len(th)}")
    if not 0.0 < lo < mid < hi <= 1.0:
        raise ValueError("mosaic_thresholds: границы должны возрастать и лежать "
                         f"в пределах (0; 1], получено {th}")
    # у границы два знака после запятой всегда: «0,20», а не «0,2» —
    # значение порога принято записывать в том виде, в каком оно задано
    def fmt(x: float) -> str:
        s = f"{x:.2f}"
        return s.replace(".", ",") if style.decimal_comma else s
    if style.term_register == "пациентский":
        L = MOSAIC_LABELS_PATIENT
        return ((0.0, lo, L[0]), (lo, mid, L[1]), (mid, hi, L[2]), (hi, 1.01, L[3]))
    return ((0.0, lo, MOSAIC_LABELS[0].format(lo=fmt(lo))),
            (lo, mid, MOSAIC_LABELS[1]),
            (mid, hi, MOSAIC_LABELS[2]),
            (hi, 1.01, MOSAIC_LABELS[3].format(hi=fmt(hi))))


def mosaic_class(f: Optional[float], style: Style) -> Optional[str]:
    if f is None:
        return None
    for a, b, name in mosaic_bands(style):
        if a <= f < b:
            return name
    return None


def _overlap(ev_start: int, ev_end: int, a: int, b: int) -> float:
    """Доля критического участка, покрытая событием."""
    inter = max(0, min(ev_end, b) - max(ev_start, a))
    return inter / (b - a) if b > a else 0.0


def region_note(ev: IF.Event, style: Style, min_overlap: float = 0.9,
                regions: Optional[Dict[Tuple[str, int, int, str], str]] = None
                ) -> Optional[str]:
    """Перекрытие критического участка — утверждение о координатах.

    ``regions`` позволяет подставить справочник, построенный из кураторского
    перечня (см. ``load_regions``); без него берётся встроенный запас.
    """
    if not style.region_overlap or ev.start is None:
        return None
    best = None
    for (chrom, a, b, direction), label in (regions or SYNDROMES).items():
        if chrom != ev.chrom or direction != ev.direction:
            continue
        cov = _overlap(ev.start, ev.end, a, b)
        if cov >= min_overlap and (best is None or cov > best[0]):
            best = (cov, label)
    if best is None:
        return None
    cov, label = best
    return (f"событие перекрывает {label} на "
            f"{_num(100 * cov, 0, style.decimal_comma)} % его протяжённости")


#: Названия уровней плоидности. Плоидность выражена в формуле множителем
#: копийности, а не отдельным событием, поэтому словесный слой обязан читать её
#: самостоятельно: проверка на выборке (случаи A20, A21, A22) показала, что без
#: этого триплоидия, тетраплоидия и гаплоидия описывались как «аберраций не
#: выявлено».
PLOIDY_NAMES = {
    1: ("гаплоидия: по одной копии каждой хромосомы вместо двух",
        "один набор хромосом вместо двух"),
    3: ("триплоидия: по три копии каждой хромосомы вместо двух",
        "три набора хромосом вместо двух"),
    4: ("тетраплоидия: по четыре копии каждой хромосомы вместо двух",
        "четыре набора хромосом вместо двух"),
}

#: Ожидаемые половые наборы при двух наборах хромосом. Всё остальное при
#: плоидности 2 — числовая аномалия половых хромосом и потому находка.
NORMAL_SEX = ("XX", "XY")

SEX_ANEUPLOIDY_NAMES = {
    "X": "моносомия X (синдром Шерешевского-Тёрнера)",
    "XXY": "набор XXY (синдром Клайнфельтера)",
    "XXX": "набор XXX (трисомия X)",
    "XYY": "набор XYY",
    "XXYY": "набор XXYY",
}


def genome_findings(case: IF.Case, style: Style) -> List[str]:
    """Находки уровня генома: плоидность, половой набор, гомозиготность.

    Читаются из полей набора событий, а не из списка событий: в формуле они
    выражены множителем копийности и составом половой группы.
    """
    pat = style.term_register == "пациентский"
    out: List[str] = []
    if case.ploidy != 2:
        nm = PLOIDY_NAMES.get(case.ploidy)
        if nm:
            out.append(nm[1] if pat else nm[0])
        else:
            out.append(f"полиплоидия: наборов хромосом {case.ploidy}")
    sx = case.sex_chromosomes
    if case.ploidy == 2 and sx and sx != "U" and sx not in NORMAL_SEX:
        # находка на половых хромосомах — это результат, а не раскрытие пола для
        # выбора эмбриона, поэтому она сообщается при любой политике; скрывается
        # только эуплоидный набор
        nm = SEX_ANEUPLOIDY_NAMES.get(sx)
        out.append(f"число половых хромосом отличается от обычного: набор {sx}"
                   if pat else
                   f"числовая аномалия половых хромосом: {nm or 'набор ' + sx}")
    if case.genome_zygosity:
        if case.genome_zygosity == "hmz":
            out.append("все участки генома представлены в одном варианте, что "
                       "отвечает полному пузырному заносу" if pat else
                       "геном полностью гомозиготен, что отвечает полному "
                       "пузырному заносу")
        else:
            out.append(f"зиготность генома: {case.genome_zygosity}")
    if case.chromatid:
        out.append("исследовано полярное тельце" if pat else
                   "исследовано полярное тельце: копийность отнесена к "
                   "гаплоидному набору хроматид")
    return out


def has_genome_finding(case: IF.Case) -> bool:
    return bool(case.ploidy != 2 or case.genome_zygosity
                or (case.ploidy == 2 and case.sex_chromosomes
                    and case.sex_chromosomes not in NORMAL_SEX + ("U",)))


def load_regions(path: str, *, min_level: int = 3,
                 require_confirmed: bool = False) -> Dict[Tuple[str, int, int, str], str]:
    """Справочник критических участков из файла, а не из памяти.

    Файл — выдача построения справочника (`critical_regions_v1.csv`): координаты
    и уровень доказательности берутся из кураторского перечня дозовой
    чувствительности, обозначение — из названия болезни и русской затравки.
    Встроенный словарь ``SYNDROMES`` остаётся только как запас на случай
    отсутствия файла: его координаты записывались при разработке по памяти и
    расходятся с кураторскими (по участку 1p36 — почти вдвое по протяжённости).

    ``min_level`` отсекает по уровню доказательности (3 — достаточные
    доказательства, 2 — некоторые). ``require_confirmed`` оставляет только
    строки, подтверждённые врачом заказчика: до подтверждения справочник в
    выдаче заключений применять не следует.
    """
    import csv

    out: Dict[Tuple[str, int, int, str], str] = {}
    with open(path, encoding="utf-8") as fh:
        rows = csv.DictReader(l for l in fh if not l.startswith("#"))
        for r in rows:
            try:
                level = int(r["уровень_доказательности"])
            except (KeyError, TypeError, ValueError):
                continue
            if level < min_level:
                continue
            if require_confirmed and r.get("подтверждено_врачом", "нет").strip().lower() != "да":
                continue
            имя = (r.get("русское_обозначение_затравка") or "").strip()
            if not имя:
                болезнь = (r.get("название_болезни") or "").strip()
                имя = (f"критический участок {r['цитобенд']}"
                       + (f" ({болезнь})" if болезнь else ""))
            ключ = (str(r["хромосома"]), int(r["начало"]), int(r["конец"]),
                    r["направление"])
            out[ключ] = имя
    return out


# ---------------------------------------------------------------------------
# Слой 1: дословная расшифровка события
# ---------------------------------------------------------------------------

def describe_event(ev: IF.Event, bands: IF.CytobandTable, style: Style,
                   assembly: str = "GRCh38") -> str:
    """Одно событие словами. Каждый факт фразы есть в формуле, и наоборот."""
    T = TERMS[style.term_register]
    band = _band_label(ev, bands)
    arm = ev.arm or (band[len(ev.chrom):len(ev.chrom) + 1] if band else "")
    forms = ARM_FORMS.get(arm, {"nom": "плечо", "gen": "плеча", "dat": "плечу"})
    subs = {"chrom": ev.chrom, "band": band, "arm": arm,
            "arm_word": forms["nom"], "arm_gen": forms["gen"],
            "arm_dat": forms["dat"], "cn": ev.copy_number}

    if ev.scale == "chromosome":
        key = "gain_whole" if ev.direction == "gain" else "loss_whole"
    elif ev.scale == "arm":
        key = "gain_arm" if ev.direction == "gain" else "loss_arm"
    elif ev.copy_number == 0:
        key = "nullisomy"
    else:
        key = "gain_seg" if ev.direction == "gain" else "loss_seg"
    head = T[key].format(**subs)

    tail: List[str] = []
    if ev.scale != "chromosome" and ev.start is not None:
        size = _mb(ev.end - ev.start + 1, style)
        if style.coordinates:
            tail.append(f"размер {size}, координаты {_thousands(ev.start)}"
                        f"-{_thousands(ev.end)} по сборке {assembly}")
        else:
            tail.append(f"размер {size}")
    if ev.copy_number_range:
        a, b = ev.copy_number_range
        tail.append(f"{T['copies']} от {a} до {b}".replace("{cn}", ""))
    elif ev.scale != "chromosome":
        # для полной хромосомы копийность уже выражена словом «трисомия» или
        # «моносомия», и повторять её значило бы удваивать один факт
        tail.append(T["copies"].format(cn=ev.copy_number))

    if ev.mosaic_fraction is not None and style.mosaic_policy != "не_указывать":
        cls = mosaic_class(ev.mosaic_fraction, style)
        frac = _fraction_text(ev.mosaic_fraction, style)
        if style.mosaic_policy == "класс_и_доля":
            tail.append(f"в мозаичной форме, доля клеток с аберрацией {frac}"
                        + (f" — {cls}" if cls else ""))
        elif style.mosaic_policy == "только_доля":
            tail.append(f"в мозаичной форме, доля клеток с аберрацией {frac}")
        elif style.mosaic_policy == "только_класс":
            base = (MOSAIC_PATIENT if style.term_register == "пациентский"
                    else "в мозаичной форме")
            tail.append(base + (f", {cls}" if cls else ""))
    elif ev.mosaic_fraction_unknown:
        tail.append("в мозаичной форме, долю клеток определить не удалось")

    if ev.structure:
        tail.append(_mechanism_text(ev, style))
    if ev.inheritance in ("mat", "pat"):
        tail.append("унаследовано от матери" if ev.inheritance == "mat"
                    else "унаследовано от отца")
    elif ev.inheritance == "dn":
        tail.append("возникло впервые, у родителей не выявлено")

    note = region_note(ev, style)
    if note:
        tail.append(note)
    return head + (": " + "; ".join(t for t in tail if t) if tail else "")


def _band_label(ev: IF.Event, bands: IF.CytobandTable) -> str:
    if ev.scale == "chromosome":
        return ev.chrom
    if ev.scale == "arm" and ev.arm:
        return f"{ev.chrom}{ev.arm}"
    if ev.start is None:
        return ev.chrom
    b1, b2 = bands.band_span(ev.chrom, ev.start, ev.end)
    return f"{ev.chrom}{b1}" if b1 == b2 else f"{ev.chrom}{b1}{b2}"


MECHANISM_WORDS = {
    "i": "по строению — изохромосома",
    "der": "по строению — производная хромосома",
    "t": "по строению — транслокация",
    "del": "по строению — делеция",
    "dup": "по строению — дупликация",
    "inv": "по строению — инверсия",
    "ins": "по строению — вставка",
    "rec": "по строению — рекомбинантная хромосома",
    "rob": "по строению — робертсоновская транслокация",
    "idic": "по строению — изодицентрическая хромосома",
    "r": "по строению — кольцевая хромосома",
}


def _mechanism_text(ev: IF.Event, style: Style) -> str:
    """Механизм упоминается только если он задан во входных данных.

    Правило 5.5.15c ISCN 2024 допускает символ перестройки лишь для
    наблюдённого механизма; словесное описание повторяет это ограничение и
    прямо помечает механизм как заявленный, а не измеренный дозой.
    """
    st = ev.structure or {}
    chain = st.get("chain") or [st]
    syms = [c.get("symbol") for c in chain if c.get("symbol")]
    words = [MECHANISM_WORDS.get(s, f"по строению — {s}") for s in syms[:2]]
    txt = "; ".join(dict.fromkeys(words))
    if style.audience != "пациент":
        txt += " (механизм заявлен, а не выведен из копийности)"
    return txt


# ---------------------------------------------------------------------------
# Слой 2: клиническая значимость
# ---------------------------------------------------------------------------

def significance(case: IF.Case, bands: IF.CytobandTable, style: Style
                 ) -> Tuple[str, List[str]]:
    """Категория находки и пояснения к ней.

    Категория не зависит от профиля получателя: клиническая значимость есть
    свойство находки, а не читателя. От профиля зависит только то, какими
    словами она объясняется и приводится ли перекрытие критического участка.
    """
    gen = genome_findings(case, style)
    if not case.events and not case.clones and not gen:
        return "норма", []
    probe = style.with_(region_overlap=True)      # классификация вне стиля
    notes: List[str] = []
    category = "неопределённая"
    if gen:
        # плоидность и числовая аномалия половых хромосом — установленные
        # состояния, а не находки неопределённой значимости
        category = "патогномоничная"
        nm = SEX_ANEUPLOIDY_NAMES.get(case.sex_chromosomes or "")
        if nm and case.ploidy == 2 and style.name_syndromes:
            if style.term_register == "пациентский":
                # в пациентском регистре остаётся только принятое название
                # состояния, без слов «моносомия» и «трисомия»
                inner = re.search(r"\(([^)]+)\)", nm)
                if inner:
                    notes.append("это состояние известно как "
                                 + inner.group(1))
            else:
                notes.append(f"половые хромосомы: {nm}")
    for ev in case.events:
        if ev.scale == "chromosome":
            category = "патогномоничная"
            nm = ANEUPLOIDY_NAMES.get((ev.chrom, ev.direction))
            if nm and style.name_syndromes:
                notes.append(f"хромосома {ev.chrom}: {nm}")
        else:
            note = region_note(ev, probe)
            if note:
                category = "патогномоничная"
                if style.region_overlap:
                    notes.append(note + "; названием синдрома находка не "
                                        "подменяется: это утверждение о координатах")
    if any(ev.mosaic_fraction is not None for ev in case.events):
        notes.append("изменение есть не во всех исследованных клетках, поэтому "
                     "перенести вывод на весь эмбрион нельзя"
                     if style.term_register == "пациентский" else
                     "доля клеток измерена в биоптате и не переносится на весь "
                     "эмбрион без оговорки")
        # ESHRE hoaa017: «The clinical significance of transferring mosaic
        # embryos is currently unknown». Без этой оговорки мозаичная находка
        # подавалась как однозначно значимая, что первоисточнику противоречит.
        notes.append("последствия переноса эмбриона с мозаичной находкой на "
                     "сегодня неизвестны"
                     if style.term_register == "пациентский" else
                     "клиническая значимость переноса эмбриона с мозаичной "
                     "находкой на настоящий момент не установлена")
    if category == "неопределённая":
        # ACGS, п. 2.5.6.3: «If no clear diagnosis can be made from the evidence
        # available, this must be clear in the report». Названия категории для
        # этого недостаточно.
        notes.append("по имеющимся данным диагноз не устанавливается"
                     if style.term_register == "пациентский" else
                     "по имеющимся данным однозначный диагноз не ставится: "
                     "значимость находки не определена")
    return category, notes


# ---------------------------------------------------------------------------
# Слой 3: ограничения
# ---------------------------------------------------------------------------

def limits(case: IF.Case, style: Style, detection_limit: Optional[float] = None,
           extra: Sequence[str] = ()) -> List[str]:
    """Ограничения результата. Формулировки зависят от регистра, состав — нет.

    Состав опирается на прочитанные первоисточники: предел заявления назначает
    сама лаборатория внутренней валидацией (ESHRE hoaa017: «Each centre should
    decide whether or not to report mosaicism based on internal validation and
    recent literature»); точный уровень мозаицизма в биоптате определить нельзя
    (там же: «As the number of cells in a TE biopsy is unknown, the exact level
    of mosaicism in the sample cannot be determined»); при исследовании
    носительства перестройки норма и сбалансированное носительство по копийности
    неразличимы, и это требуется указывать прямо.

    Слово «отрицательный» о результате не употребляется: ACGS, п. 2.5.6.
    """
    pat = style.term_register == "пациентский"
    out: List[str] = []
    if detection_limit is not None:
        dl = _num(detection_limit, 3, style.decimal_comma)
        if pat:
            out.append(LIMITS_PATIENT["предел"].format(
                dl=f"{_num(100 * detection_limit, 1, style.decimal_comma)} %"))
        else:
            out.append(f"событие с долей затронутых клеток ниже {dl} методом не "
                       "выявляется; предел заявления установлен внутренней "
                       "валидацией лаборатории")
    if case.technique == "sseq":
        out.append("метод не различает изменения меньше 5-10 Мб" if pat else
                   "разрешение метода при малой глубине прочтения составляет "
                   "около 5-10 Мб: события меньшего размера метод не различает")
    if any(ev.chrom in ("X", "Y") for ev in case.events):
        out.append(LIMITS_PATIENT["половые"] if pat else
                   "для половых хромосом применяется отдельная модель, и её "
                   "ограничения указаны в разделе о канале половых хромосом")
    if any(ev.mosaic_fraction is not None for ev in case.events):
        out.append(LIMITS_PATIENT["трофэктодерма"] if pat else
                   "оценка доли клеток получена по биоптату трофэктодермы; "
                   "соответствие внутренней клеточной массе не гарантировано")
        out.append("доля клеток в биоптате точной величиной не является: число "
                   "клеток в биоптате неизвестно, поэтому истинный уровень "
                   "мозаицизма в образце установить нельзя"
                   if not pat else
                   "доля затронутых клеток — оценка, а не точное значение")
    if case.events and not any(ev.structure for ev in case.events):
        out.append(LIMITS_PATIENT["механизм"] if pat else
                   "механизм перестройки методом не наблюдается: по копийности "
                   "делеция, кольцевая хромосома и несбалансированный продукт "
                   "транслокации неразличимы")
    if case.karyotype_available or any(ev.structure for ev in case.events):
        out.append("сбалансированное носительство перестройки и нормальный "
                   "набор по копийности неразличимы: результат без находки не "
                   "отличает эмбрион-носитель от эмбриона без перестройки"
                   if not pat else
                   "исследование не отличает эмбрион, унаследовавший "
                   "перестройку в сбалансированном виде, от эмбриона без неё")
    out.extend(extra)
    return out


# ---------------------------------------------------------------------------
# Слой 4: что дальше
# ---------------------------------------------------------------------------

def next_steps(case: IF.Case, style: Style) -> List[str]:
    """Возможности, а не предписания.

    Рекомендации по оформлению заключений (Claustres et al., EJHG 2014) прямо
    запрещают формулировки «показана» и «необходима» в отношении пренатальной
    или досимптоматической диагностики: лаборатория сообщает, что исследование
    возможно, а решение принимает врач с пациентом.
    """
    out: List[str] = []
    if not case.events and not case.clones and not has_genome_finding(case):
        out.append("дополнительных исследований по результатам настоящего "
                   "анализа не требуется")
    else:
        out.append("результат можно подтвердить другим методом"
                   if style.term_register == "пациентский" else
                   "подтверждение находки независимым методом возможно "
                   "(кариотипирование, микрочиповое исследование, "
                   "количественная ПЦР по участку)")
        if any(ev.structure for ev in case.events):
            out.append("исследование кариотипа родителей возможно и позволяет "
                       "установить, унаследована ли перестройка")
        out.append("пренатальная диагностика при наступившей беременности "
                   "возможна независимо от результата настоящего исследования")
    if style.counselling_note:
        out.append("результат подлежит обсуждению на медико-генетическом "
                   "консультировании")
    return out


# ---------------------------------------------------------------------------
# Сборка описания
# ---------------------------------------------------------------------------

#: Формулировки, запрещённые в заключении. Первые две — по рекомендациям
#: Claustres et al. (EJHG 2014): слова «положительный» и «отрицательный»
#: двусмысленны, а предписывать пренатальную диагностику лаборатория не вправе.
FORBIDDEN = (
    (r"\bположительн\w*\s+результат", "слово «положительный» о результате "
                                      "двусмысленно"),
    # ACGS, п. 2.5.6: «The use of the terms ‘positive’ and ‘negative’ in
    # relation to a variant is not recommended but must be clearly defined if
    # used». Прежний образец искал только сочетания «означает» и «исключает»
    # и пропускал оборот «отрицательный результат заявляется до доли …», из-за
    # чего проверка не ловила нарушение того самого правила, которым
    # обосновано её существование.
    (r"\bотрицательн\w*\s+результат", "слово «отрицательный» о результате "
                                       "двусмысленно (ACGS, п. 2.5.6)"),
    (r"показан\w*\s+(пренатальн|инвазивн)", "лаборатория не вправе предписывать "
                                            "пренатальную диагностику"),
    # ESHRE hoaa021 прямо допускает указание на необходимость подтверждения
    # при наступившей беременности, поэтому запрет сужен: он касается
    # предписания диагностики как показания, а не подтверждения результата.
    (r"необходим\w*\s+(пренатальн|инвазивн)(?![^.]{0,80}подтвержд)",
     "предписывать пренатальную диагностику лаборатория не вправе; "
     "указание на необходимость подтверждения результата при наступившей "
     "беременности допустимо"),
    (r"\bздоров\w*\s+эмбрион", "«здоровый эмбрион» выходит за пределы "
                               "исследованного"),
    (r"\bнорм\w*\s+эмбрион", "то же"),
    (r"(?<!не )(?<!не\u00a0)гарантиру\w+", "результат гарантий не даёт, "
                                       "и обещать их в заключении нельзя"),
)

SEX_MARKERS = (r"\bXX\b", r"\bXY\b", r"\bмужск", r"\bженск", r"\bпол\s+эмбриона")


def _sex_sentence(case: IF.Case, style: Style) -> Optional[str]:
    if style.sex_policy == "не_раскрывать":
        return None
    sx = case.sex_chromosomes
    if not sx or sx == "U" or not case.sex_disclosed:
        return "половой набор не раскрывается"
    # при аномальной плоидности состав половой группы клинически значим: XXY
    # против XXX при триплоидии различает отцовское и материнское происхождение
    # лишнего набора, поэтому находка уровня генома тоже считается находкой
    has_sex_event = (any(ev.chrom in ("X", "Y") for ev in case.events)
                     or any(ev.chrom in ("X", "Y") for cl in case.clones
                            for ev in cl.events)
                     or has_genome_finding(case))
    if style.sex_policy == "при_находке" and not has_sex_event:
        return None
    return f"половой набор по данным исследования: {sx}"


def _cap(s: str) -> str:
    return s[0].upper() + s[1:] if s else s


def _cut(items: Sequence[str], style: Style) -> List[str]:
    """Ограничение длины для пояснений.

    К перечню находок не применяется никогда: сокращение пояснения — вопрос
    стиля, а сокращение перечня находок означает, что часть выявленного не
    попала в заключение. Проверка на выборке (случай A08 из пяти полных
    событий и образец NA02325 из трёх) показала, что прежняя сплошная обрезка
    именно это и делала.
    """
    if style.max_sentences is None:
        return list(items)
    return list(items[:style.max_sentences])


def verbalize(case: IF.Case, bands: IF.CytobandTable,
              style: Style | str = "генетик", *,
              detection_limit: Optional[float] = None,
              extra_limits: Sequence[str] = (),
              formula: Optional[str] = None,
              clinical_question: Optional[str] = None,
              transfer_note: Optional[str] = None,
              sex_nondisclosure_note: bool = True) -> Dict[str, Any]:
    """Словесное описание результата по набору событий.

    Возвращает разделы описания по отдельности и собранный текст, чтобы форма
    отчёта могла вставлять их поблочно, а бланк заключения — целиком.
    """
    st = PROFILES[style] if isinstance(style, str) else style
    out: Dict[str, Any] = {"профиль": st.name, "разделы": {}}

    # ACGS, п. 2.5.2: клинический вопрос повторяется в заключении. Текст
    # подставляет вызывающая сторона — показания инструмент не знает и
    # придумывать их не должен.
    if st.restate_question and clinical_question:
        out["разделы"]["вопрос"] = clinical_question.strip()

    if formula is None and st.include_formula:
        try:
            formula = IF.render(case, st.formula_form, bands)
        except IF.Unrenderable:
            formula = None
    if st.include_formula and formula:
        out["разделы"]["формула"] = formula

    gen = [_cap(s) + "." for s in genome_findings(case, st)]
    if case.normal() and not gen:
        lit = ["Аберраций числа и структуры хромосом в пределах разрешения "
               "метода не выявлено."]
    else:
        lit = [_cap(describe_event(ev, bands, st, case.assembly)) + "."
               for ev in case.events]
        for cl in case.clones:
            for ev in cl.events:
                lit.append(_cap(describe_event(ev, bands, st, case.assembly))
                           + " (отдельная клеточная линия).")
    lit = gen + lit
    sx = _sex_sentence(case, st)
    if sx:
        lit.append(sx[0].upper() + sx[1:] + ".")
    if st.literal_layer:
        # перечень находок не обрезается ни при каком значении max_sentences
        out["разделы"]["расшифровка"] = lit

    # Riggs и др.: политика неразглашения пола допустима, но должна быть
    # заявлена в заключении, иначе молчание неотличимо от отсутствия находки.
    if (sex_nondisclosure_note and st.sex_policy == "не_раскрывать"
            and st.literal_layer):
        out["разделы"].setdefault("расшифровка", lit).append(
            "Половой набор в заключении не приводится: это решение лаборатории, "
            "а не отсутствие данных.")

    if st.significance_layer:
        cat, notes = significance(case, bands, st)
        out["категория"] = cat
        words = (CATEGORIES_PATIENT if st.term_register == "пациентский"
                 else CATEGORIES)[cat]
        # строка категории остаётся всегда; ограничение длины действует только
        # на пояснения к ней
        block = [_cap(words) + "."]
        block += _cut([_cap(n) + "." for n in notes], st)
        out["разделы"]["значимость"] = block

    if st.limits_layer:
        out["разделы"]["ограничения"] = [
            s[0].upper() + s[1:] + "." for s in
            limits(case, st, detection_limit, extra_limits)]
    if st.next_steps_layer:
        out["разделы"]["что_дальше"] = [
            s[0].upper() + s[1:] + "." for s in next_steps(case, st)]

    # ESHRE hoaa021 требует в заключении цикла ПГТ указания о пригодности
    # эмбриона к переносу. Инструмент его НЕ выводит: пригодность решает врач,
    # и переносится только переданный текст. Отсутствие текста при включённой
    # настройке отмечается прямо, чтобы пропуск был видим.
    if st.transfer_guidance:
        out["разделы"]["пригодность"] = (
            transfer_note.strip() if transfer_note else
            "Указание о пригодности к переносу лабораторией не сформулировано: "
            "его вносит врач.")

    # ACGS, п. 2.1: вывод должен быть явно виден. Он собирается из категории и
    # перечня находок и ставится перед пояснениями.
    if st.highlight_conclusion and "значимость" in out["разделы"]:
        итог = out["разделы"]["значимость"][0]
        находки = out["разделы"].get("расшифровка") or []
        главное = " ".join(находки[:2]) if находки else ""
        out["разделы"]["вывод"] = (главное + " " + итог).strip()

    out["текст"] = _assemble(out["разделы"], st)
    out["проверка"] = check_text(out["текст"], case, bands, st)
    return out


HEADINGS = {"вопрос": "Клинический вопрос", "формула": "Запись по ISCN 2024",
            "расшифровка": "Что выявлено", "вывод": "Вывод",
            "значимость": "Клиническая значимость", "ограничения": "Ограничения",
            "пригодность": "Пригодность к переносу",
            "что_дальше": "Что можно сделать дальше"}

#: Порядок разделов в собранном тексте. Вывод стоит перед пояснениями: ACGS,
#: п. 2.1 требует, чтобы вывод был явно виден, а не отыскивался в конце.
SECTION_ORDER = ("вопрос", "вывод", "формула", "расшифровка", "значимость",
                 "ограничения", "пригодность", "что_дальше")


def _assemble(sections: Dict[str, Any], style: Style) -> str:
    parts: List[str] = []
    for key in SECTION_ORDER:
        if key not in sections:
            continue
        body = sections[key]
        head = style.heading_prefix + HEADINGS[key]
        if isinstance(body, str):
            parts.append(f"**{head}.** `{body}`")
        elif style.bullet_layers:
            parts.append(f"**{head}.**\n" + "\n".join(f"- {s}" for s in body))
        else:
            parts.append(f"**{head}.** " + " ".join(body))
    return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# Обратная проверка текста
# ---------------------------------------------------------------------------

def check_text(text: str, case: IF.Case, bands: IF.CytobandTable,
               style: Style) -> Dict[str, Any]:
    """Извлечь факты из готового текста и сверить с набором событий.

    Проверка обратная по построению: она разбирает текст, а не набор событий,
    поэтому пропущенная или лишняя хромосома, потерянная доля клеток и
    несогласованное направление события обнаруживаются механически. Отдельно
    проверяются запрещённые формулировки и политика раскрытия пола.
    """
    problems: List[str] = []

    # хромосомы, названные в тексте
    named = set(re.findall(r"хромосом[ыае]?\s+([0-9]{1,2}|X|Y)", text))
    # обозначение полос: перед номером хромосомы не должно стоять ни цифры, ни
    # точки, иначе «22q11.1q12.2» даст ложную хромосому 1
    named |= set(re.findall(r"(?<![0-9A-Za-z.])([0-9]{1,2}|X|Y)[pq][0-9]", text))
    # половой набор записан слитным обозначением: «набор XXY», «моносомия X».
    # Без этого правила проверка считала половые хромосомы потерянными, хотя
    # текст их называет (выявлено на случаях A20 и A22 выборки).
    for grp in re.findall(r"набор[а-я]*[^.;]{0,40}?\b([XY]{1,4})\b", text):
        named |= set(grp)
    for grp in re.findall(r"(?:моносомия|трисомия|дисомия|копия|копии)\s+([XY]{1,4})\b", text):
        named |= set(grp)
    required = {ev.chrom for ev in case.events}
    required |= {ev.chrom for cl in case.clones for ev in cl.events}
    # Заявленный половой набор упоминать не обязательно (профиль может его
    # скрывать), но и выдумкой его упоминание не является: при наборе XXY
    # названная в тексте Y — часть заявленного набора, а не лишняя хромосома.
    # Поэтому обязательное и допустимое разделены.
    allowed = set(required)
    if case.sex_chromosomes and case.sex_chromosomes != "U":
        allowed |= set(case.sex_chromosomes)
    if style.literal_layer:
        missing = required - named
        if missing:
            problems.append("в тексте не названы хромосомы: "
                            + ", ".join(sorted(missing)))
        extra = named - allowed
        if extra:
            problems.append("в тексте названы хромосомы, которых нет в наборе "
                            "событий: " + ", ".join(sorted(extra)))

    # направление события
    for ev in case.events:
        if ev.scale == "chromosome" and style.literal_layer:
            want = ("трисом", "лишн", "три копии") if ev.direction == "gain" \
                else ("моносом", "нехватк", "одна копия", "утрат")
            if not any(w in text.lower() for w in want):
                problems.append(f"направление события на хромосоме {ev.chrom} "
                                "в тексте не выражено")

    # доля клеток
    fr = [ev.mosaic_fraction for ev in case.events if ev.mosaic_fraction is not None]
    if fr and style.mosaic_policy in ("класс_и_доля", "только_доля"):
        for f in fr:
            if _num(f, 2, style.decimal_comma) not in text:
                problems.append(f"доля клеток {_num(f, 2, style.decimal_comma)} "
                                "в тексте отсутствует")
    if fr and style.mosaic_policy != "не_указывать":
        # в пациентском регистре слова «мозаичный» нет по построению, поэтому
        # проверка принимает и его описательную замену
        said = ("мозаич" in text.lower()
                or MOSAIC_PATIENT[:40] in text
                or "не во всех" in text.lower())
        if not said:
            problems.append("мозаичная форма события в тексте не упомянута")

    # запрещённые формулировки
    for pat, why in FORBIDDEN:
        if re.search(pat, text, re.I):
            problems.append(f"запрещённая формулировка ({why}): {pat}")

    # политика раскрытия пола
    if style.sex_policy == "не_раскрывать":
        for pat in SEX_MARKERS:
            if re.search(pat, text):
                problems.append("профиль запрещает раскрытие пола, но текст "
                                f"содержит {pat}")

    return {"итог": "ок" if not problems else "расхождение",
            "замечания": problems,
            "хромосом_в_тексте": sorted(named),
            "хромосом_в_наборе": sorted(required),
            "хромосом_допустимо": sorted(allowed)}


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    import json as _json
    ap = argparse.ArgumentParser(description="Словесное описание по формуле ISCN")
    ap.add_argument("--input", required=True, help="JSON набора событий")
    ap.add_argument("--cytobands", required=True)
    ap.add_argument("--profile", default="генетик", choices=sorted(PROFILES))
    ap.add_argument("--detection-limit", type=float)
    ap.add_argument("--set", action="append", default=[],
                    metavar="ПОЛЕ=ЗНАЧЕНИЕ", help="переопределить поле профиля")
    ap.add_argument("--out")
    args = ap.parse_args(argv)

    case = IF.case_from_dict(_json.load(open(args.input)))
    bands = IF.CytobandTable(args.cytobands)
    st = PROFILES[args.profile]
    for kv in args.set:
        k, _, v = kv.partition("=")
        cur = getattr(st, k)
        if isinstance(cur, bool):
            v = v.lower() in ("1", "true", "да", "yes")
        elif isinstance(cur, tuple):
            v = tuple(float(x) for x in re.split(r"[,;]", str(v)) if x.strip())
        elif isinstance(cur, float):
            v = float(str(v).replace(",", "."))
        elif isinstance(cur, int) and not isinstance(cur, bool):
            v = int(v)
        st = st.with_(**{k: v})
    res = verbalize(case, bands, st, detection_limit=args.detection_limit)
    text = res["текст"]
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    else:
        print(text)
    if res["проверка"]["итог"] != "ок":
        print("\n[проверка] " + "; ".join(res["проверка"]["замечания"]))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
