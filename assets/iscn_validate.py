"""Входной контроль записи ISCN: синтаксис и соответствие стандарту.

Зачем отдельный слой. Разборщик и проверки соглашений в ``iscn_formula``
построены для обратной проверки собственного вывода: своя строка всегда в
области применения и всегда синтаксически верна, а проверка либо молчит, либо
сообщает о дефекте генератора. Для записи, пришедшей извне — из другой
лаборатории, из архивного заключения, набранной врачом вручную, полученной
информационной системой, — этого недостаточно по трём причинам.

1.  Соматические обозначения (``dmin``, ``cp``, ``sl``, ``sdl``) выводят строку
    за область применения конституциональной номенклатуры, и проверки в
    ``iscn_formula`` бросают на них исключение. Внешней записи нужен вердикт
    «вне области применения», а не отказ работы.
2.  Верная по существу запись может не разобраться из-за оформления:
    типографский минус вместо дефиса, латинская ``x`` вместо знака умножения,
    неразрывные пробелы в разрядах. Такую запись следует привести к
    каноническому виду и сообщить об этом, а не отвергнуть.
3.  Не всякое отступление одинаково тяжело. Недопустимая доля мозаицизма — это
    ошибка, неоднозначный разбор — предупреждение, отступление от предпочтительной
    формы — замечание. Без уровней вердикт непригоден для принятия решения.

Модуль ничего не разбирает сам: он нормализует вход, вызывает разборщик и
проверки ``iscn_formula``, добавляет проверки полос кариотипной строки по
таблице цитобендов и сводит всё в вердикт с уровнями.
"""

from __future__ import annotations

import os
import re
import sys
import unicodedata
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

# При запуске из командной строки каталог модуля может не оказаться в путях
# поиска (окружение выставляет PYTHONSAFEPATH), поэтому добавляется явно.
if os.path.dirname(os.path.abspath(__file__)) not in sys.path:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import iscn_formula as IF

__all__ = ["normalize", "validate", "Finding", "VERDICTS", "band_findings"]

#: Итоговые вердикты. Порядок — от лучшего к худшему.
VERDICTS = ("годна", "годна с замечаниями", "требует решения", "негодна",
            "вне области применения")

SEVERITIES = ("ошибка", "предупреждение", "замечание")


@dataclass
class Finding:
    """Одно замечание к записи."""

    уровень: str
    правило: str
    сообщение: str
    место: str = ""

    def as_dict(self) -> Dict[str, str]:
        return {"уровень": self.уровень, "правило": self.правило,
                "сообщение": self.сообщение, "место": self.место}


# ---------------------------------------------------------------------------
# Нормализация оформления
# ---------------------------------------------------------------------------

#: Замены оформления. Каждая — приведение к тому виду, который печатает сам
#: стандарт; смысл записи ни одна из них не меняет.
SUBSTITUTIONS: Tuple[Tuple[str, str, str], ...] = (
    ("\u2212", "-", "типографский минус заменён дефисом"),
    ("\u2013", "-", "короткое тире заменено дефисом"),
    ("\u2014", "-", "длинное тире заменено дефисом"),
    ("\u00a0", " ", "неразрывный пробел заменён обычным"),
    ("\u2007", " ", "цифровой пробел заменён обычным"),
    ("\u202f", " ", "узкий неразрывный пробел заменён обычным"),
    ("\u00d7", "\u00d7", ""),          # знак умножения — канонический вид
    ("\u0445", "\u00d7", "русская «х» заменена знаком умножения"),
    ("\u2116", "N", "знак номера заменён латинской N"),
)


def normalize(text: str) -> Tuple[str, List[str]]:
    """Привести оформление к каноническому виду.

    Возвращает нормализованную строку и перечень выполненных замен. Замены
    касаются только оформления: знаков, пробелов и регистра обозначений
    техники. Символы номенклатуры не добавляются и не удаляются.
    """
    notes: List[str] = []
    out = unicodedata.normalize("NFC", text).strip()
    for src, dst, note in SUBSTITUTIONS:
        if src != dst and src in out:
            out = out.replace(src, dst)
            if note:
                notes.append(note)
    # латинская x между цифрой/скобкой и цифрой — это знак умножения
    fixed = re.sub(r"(?<=[\)\d])\s?[xX](?=\d)", "\u00d7", out)
    if fixed != out:
        notes.append("латинская x заменена знаком умножения")
        out = fixed
    # разряды координат: пробел внутри числа стандарт печатает запятой.
    # Замена повторяется до устойчивого состояния: в «46 000 000» пробелы идут
    # подряд, и однократная замена исправила бы только первый.
    before = out
    while True:
        fixed = re.sub(r"(?<=\d) (?=\d{3}(?:[ ,]\d{3})*(?!\d))", ",", out)
        if fixed == out:
            break
        out = fixed
    if out != before:
        notes.append("пробелы в разрядах координат заменены запятыми")
    fixed = re.sub(r"[ \t]{2,}", " ", out)
    if fixed != out:
        notes.append("повторяющиеся пробелы сжаты")
        out = fixed
    # обозначение техники стандарт печатает строчными
    m = re.match(r"^(ARR|SEQ|SSEQ|OGM|MOS|CHI)\b", out)
    if m:
        out = m.group(1).lower() + out[m.end():]
        notes.append("обозначение техники приведено к строчным буквам")
    return out, notes


# ---------------------------------------------------------------------------
# Проверки полос кариотипной строки по таблице цитобендов
# ---------------------------------------------------------------------------

#: Обозначение точки разрыва внутри кариотипной строки: хромосома в первых
#: скобках, полосы во вторых. Разбирается уже проверенная грамматикой строка,
#: поэтому достаточно поверхностного поиска.
_BP = re.compile(r"\(([^()]*?)\)\(([^()]*?)\)")
_BAND = re.compile(r"([pq])(\d+(?:\.\d+)?)")


def _band_names(bands: IF.CytobandTable, chrom: str) -> Dict[str, Tuple[int, int]]:
    try:
        rows = bands.bands[chrom]
    except (KeyError, AttributeError):
        return {}
    return {name: (start, end) for start, end, name in rows}


def band_findings(text: str, bands: IF.CytobandTable) -> List[Finding]:
    """Существование полосы и согласие плеча с таблицей цитобендов.

    В микрочиповом и молекулярном форматах координаты уже сверяются функцией
    ``check_band_consistency`` из ``iscn_formula``; кариотипная строка
    координат не содержит, и её полосы до сих пор не проверялись ничем.

    Порядок точек разрыва (проксимальная перед дистальной) здесь **не**
    проверяется: правило не удалось подтвердить выдержкой из первоисточника по
    доступным материалам, а вводить проверку по памяти нельзя — она отвергала
    бы верные записи. Место для неё оставлено осознанно.
    """
    out: List[Finding] = []
    if re.match(r"^\s*(arr|seq|sseq|ogm)\b", text):
        return out                      # не кариотипная строка
    for chroms_part, bands_part in _BP.findall(text):
        chroms = [c for c in re.split(r"[;,]", chroms_part) if re.fullmatch(r"X|Y|\d{1,2}", c.strip())]
        groups = [g for g in re.split(r"[;]", bands_part)]
        for i, grp in enumerate(groups):
            found = _BAND.findall(grp)
            if not found:
                continue
            chrom = chroms[i].strip() if i < len(chroms) else (chroms[0].strip() if chroms else None)
            if chrom is None:
                continue
            table = _band_names(bands, chrom)
            if not table:
                out.append(Finding("замечание", "4.4.5a",
                                   f"хромосомы {chrom} нет в таблице цитобендов: "
                                   "проверить полосы невозможно", grp))
                continue
            for arm, num in found:
                name = f"{arm}{num}"
                if name in table:
                    continue
                if num == "10":
                    # p10 и q10 обозначают центромеру в цельноплечевых
                    # перестройках (der(14;21)(q10;q10), i(X)(q10)); полосы с
                    # таким номером в таблице цитобендов нет по построению
                    continue
                # полоса могла быть указана до уровня подполосы, которого в
                # таблице нет: стандарт это допускает, если существует более
                # общая полоса
                # Указание до уровня полосы там, где таблица различает
                # подполосы (9q34 при наличии q34.11, q34.12, q34.13), —
                # законная запись, а не отступление: разрешение указания
                # определяется разрешением исследования. Замечания здесь нет,
                # иначе всякая настоящая кариотипная строка получала бы шум.
                prefix = [b for b in table if b.startswith(name)]
                broader = [b for b in table if name.startswith(b)]
                if prefix or broader:
                    pass
                else:
                    out.append(Finding("ошибка", "4.4.5a",
                                       f"полосы {chrom}{name} в сборке не "
                                       "существует", grp))
    return out


#: Пара координат в короткой или расширенной системе: «(начало_конец)».
_COORDS = re.compile(r"(X|Y|\d{1,2})([pq][0-9.]+(?:[pq][0-9.]+)?)?\s*\((\d[\d,]*)_(\d[\d,]*)\)")


def coordinate_findings(text: str, bands: IF.CytobandTable) -> List[Finding]:
    """Монотонность координат и их попадание в длину хромосомы.

    Ни разбор, ни проверка соответствия полос координатам этого не ловят:
    грамматика видит два числа, а проверка полос сравнивает обозначение с
    полосой, в которую попадает каждая координата, и при обратном порядке
    сообщает о согласии. Проверка арифметическая и в первоисточнике не
    нуждается.
    """
    out: List[Finding] = []
    for chrom, band, s_raw, e_raw in _COORDS.findall(text):
        start, end = int(s_raw.replace(",", "")), int(e_raw.replace(",", ""))
        место = f"{chrom}{band or ''}({s_raw}_{e_raw})"
        if start >= end:
            out.append(Finding("ошибка", "4.4.5c",
                               f"начало {s_raw} не меньше конца {e_raw}: "
                               "координаты приводятся от меньшей к большей",
                               место))
        try:
            length = max(e for _, e, _ in bands.bands[chrom])
        except (KeyError, AttributeError):
            continue
        for имя, знач in (("начало", start), ("конец", end)):
            if знач > length:
                out.append(Finding("ошибка", "4.4.5c",
                                   f"{имя} {знач} выходит за длину хромосомы "
                                   f"{chrom} ({length} по сборке)", место))
        if start == 0:
            out.append(Finding("предупреждение", "4.4.5c",
                               "координата 0: нумерация оснований в записи "
                               "ISCN начинается с единицы", место))
    return out


# ---------------------------------------------------------------------------
# Вердикт
# ---------------------------------------------------------------------------

#: Уровень, присваиваемый замечанию проверки соглашений, если она его не задала.
_DEFAULT_SEVERITY = {"error": "ошибка", "warning": "предупреждение",
                     "note": "замечание"}


def _as_findings(items: Sequence[Any], правило_по_умолчанию: str,
                 уровень_по_умолчанию: str) -> List[Finding]:
    out: List[Finding] = []
    for it in items or ():
        if isinstance(it, Finding):
            out.append(it)
        elif isinstance(it, dict):
            sev = _DEFAULT_SEVERITY.get(str(it.get("severity", "")).lower(),
                                        уровень_по_умолчанию)
            out.append(Finding(sev, str(it.get("rule", правило_по_умолчанию)),
                               str(it.get("message", it)), str(it.get("where", ""))))
        else:
            out.append(Finding(уровень_по_умолчанию, правило_по_умолчанию, str(it)))
    return out


def validate(text: str, bands: Optional[IF.CytobandTable] = None, *,
             normalise: bool = True, require_unambiguous: bool = True
             ) -> Dict[str, Any]:
    """Проверить внешнюю запись ISCN и выдать вердикт.

    Вердикты: ``годна``, ``годна с замечаниями``, ``требует решения``
    (запись разобрана, но неоднозначна либо содержит предупреждения),
    ``негодна`` (синтаксис или ошибка соответствия), ``вне области применения``
    (соматическая номенклатура и прочее, для чего настоящий комплект правил не
    предназначен).
    """
    исходная = text
    нормализована, замены = (normalize(text) if normalise else (text, []))
    findings: List[Finding] = [
        Finding("замечание", "оформление", note) for note in замены]

    res: Dict[str, Any] = {"строка": исходная, "нормализована": нормализована,
                           "замены_оформления": замены, "деревьев_разбора": None,
                           "факты": None, "вне_области": False}

    try:
        разбор = IF.parse(нормализована)
        res["факты"] = sorted(map(str, разбор["facts"]))
        res["мета"] = разбор.get("meta", {})
    except IF.OutOfScope as exc:
        findings.append(Finding("замечание", "область применения",
                                f"запись вне области конституциональной "
                                f"номенклатуры настоящего комплекта: {exc}"))
        res["вне_области"] = True
        res["итог"] = "вне области применения"
        res["замечания"] = [f.as_dict() for f in findings]
        res["по_уровням"] = _tally(findings)
        return res
    except IF.ParseError as exc:
        findings.append(Finding("ошибка", "синтаксис",
                                f"запись не разбирается грамматикой: {exc}"))
        res["итог"] = "негодна"
        res["замечания"] = [f.as_dict() for f in findings]
        res["по_уровням"] = _tally(findings)
        return res
    except Exception as exc:                              # pragma: no cover
        findings.append(Finding("ошибка", "синтаксис",
                                f"сбой разбора: {type(exc).__name__}: {exc}"))
        res["итог"] = "негодна"
        res["замечания"] = [f.as_dict() for f in findings]
        res["по_уровням"] = _tally(findings)
        return res

    # однозначность: несколько деревьев разбора означают, что запись читается
    # больше чем одним способом
    try:
        n = IF.count_derivations(нормализована)
        res["деревьев_разбора"] = n
        if n > 1 and require_unambiguous:
            findings.append(Finding("предупреждение", "однозначность",
                                    f"запись допускает {n} разборов: смысл "
                                    "определяется неоднозначно"))
    except Exception:
        pass

    if bands is not None:
        for fn, правило, уровень in ((IF.check_conventions, "соглашения", "предупреждение"),
                                     (IF.check_band_consistency, "полосы", "ошибка")):
            try:
                items = fn(нормализована, bands) or []
                if правило == "полосы":
                    # функция возвращает отчёт по каждому обозначению, а не
                    # перечень нарушений: строка с consistent=True означает,
                    # что обозначение согласуется с координатами
                    items = [_band_message(row) for row in items
                             if isinstance(row, dict) and row.get("consistent") is False]
                findings += _as_findings(items, правило, уровень)
            except IF.OutOfScope as exc:
                findings.append(Finding("замечание", "область применения",
                                        f"проверка «{правило}» неприменима: {exc}"))
            except ValueError as exc:
                # координаты вне хромосомы: причину сообщает проверка координат,
                # поэтому здесь остаётся одна пометка о невыполненных проверках
                findings.append(Finding("замечание", "проверки",
                                        "часть проверок не выполнена из-за "
                                        "недопустимых координат"))
            except Exception as exc:                      # pragma: no cover
                findings.append(Finding("замечание", правило,
                                        f"проверка не выполнена: "
                                        f"{type(exc).__name__}: {exc}"))
        findings += band_findings(нормализована, bands)
        findings += coordinate_findings(нормализована, bands)

    findings = _dedup(findings)
    tally = _tally(findings)
    if tally["ошибка"]:
        res["итог"] = "негодна"
    elif tally["предупреждение"]:
        res["итог"] = "требует решения"
    elif tally["замечание"]:
        res["итог"] = "годна с замечаниями"
    else:
        res["итог"] = "годна"
    res["замечания"] = [f.as_dict() for f in findings]
    res["по_уровням"] = tally
    return res


def _band_message(row: Dict[str, Any]) -> Dict[str, Any]:
    """Строка отчёта о полосах — в читаемое сообщение."""
    return {"rule": "4.4.5a", "severity": "error",
            "message": (f"обозначение {row.get('chromosome')}{row.get('declared')} "
                        f"не согласуется с координатами "
                        f"{row.get('start')}_{row.get('end')} по сборке "
                        f"{row.get('build')}: координаты отвечают "
                        f"{row.get('chromosome')}{row.get('implied')} "
                        f"(отношение: {row.get('relation')})"),
            "where": str(row.get("declared", ""))}


def _dedup(findings: Sequence[Finding]) -> List[Finding]:
    """Снять повтор одного и того же нарушения из двух проверок.

    Несоответствие полосы координатам сообщают и ``check_conventions``, и
    ``check_band_consistency``; в вердикте оно должно быть одно, с более
    подробной формулировкой.
    """
    подробные = {f.место for f in findings
                 if f.правило == "4.4.5a" and "не согласуется с координатами" in f.сообщение}
    out: List[Finding] = []
    видели: set = set()
    for f in findings:
        if ("не соответствует координатам" in f.сообщение
                and any(b and b in f.сообщение for b in подробные)):
            continue
        ключ = (f.уровень, f.правило, f.сообщение)
        if ключ in видели:
            continue
        видели.add(ключ)
        out.append(f)
    return out


def _tally(findings: Sequence[Finding]) -> Dict[str, int]:
    return {s: sum(1 for f in findings if f.уровень == s) for s in SEVERITIES}


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    import json
    import sys
    ap = argparse.ArgumentParser(description="Входной контроль записи ISCN")
    ap.add_argument("--text", action="append", default=[])
    ap.add_argument("--file", help="файл со строками, по одной на строку")
    ap.add_argument("--cytobands", required=True)
    ap.add_argument("--no-normalise", action="store_true")
    ap.add_argument("--json", action="store_true", help="выдача в JSON")
    args = ap.parse_args(argv)

    texts = list(args.text)
    if args.file:
        texts += [l.strip() for l in open(args.file) if l.strip()]
    if not texts:
        texts = [l.strip() for l in sys.stdin if l.strip()]
    bands = IF.CytobandTable(args.cytobands)

    worst = 0
    for t in texts:
        v = validate(t, bands, normalise=not args.no_normalise)
        worst = max(worst, VERDICTS.index(v["итог"]))
        if args.json:
            print(json.dumps(v, ensure_ascii=False))
        else:
            print(f"\n{t}\n  вердикт: {v['итог']}"
                  + (f"\n  нормализована: {v['нормализована']}"
                     if v["нормализована"] != t else ""))
            for f in v["замечания"]:
                print(f"  [{f['уровень']}] {f['правило']}: {f['сообщение']}")
    return 0 if worst <= VERDICTS.index("годна с замечаниями") else 1


if __name__ == "__main__":
    raise SystemExit(main())
