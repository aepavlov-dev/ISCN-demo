"""Программный слой демонстрационного стенда.

Единственная задача этого файла — принять запрос от страницы, вызвать
проверенный модуль `iscn_formula` и вернуть ответ в виде JSON. Ни одного
правила номенклатуры здесь нет: правила живут в модуле, который прошёл набор
проверок. Поэтому стенд не может показать результат, отличающийся от того,
что выдаёт командная строка.
"""

import json
import os
from typing import Dict

import iscn_formula as IF
import iscn_validate as IV
import iscn_verbalize as VB

BANDS = {}
CASES = []
EXERCISE = []
REGIONS = {}          # справочник критических участков, если файл рядом с сайтом

RT_RU = {
    "ok": "разобрано обратно, совпало",
    "not_admissible": "форма неприменима",
    "mismatch": "расхождение обратного разбора",
    "parse_failed": "строка не разобралась",
    "out_of_scope": "вне области разборщика",
}

LOSS_RU = {
    "chromosome": "хромосома", "direction": "направление изменения",
    "copy_number": "копийность", "breakpoint_bands": "полосы точек разрыва",
    "coordinates": "координаты", "flanking_normal": "фланкирующие позиции",
    "mosaic_fraction": "доля мозаицизма", "mechanism": "механизм перестройки",
    "parental_origin": "родительское происхождение", "zygosity": "зиготность",
    "sex_complement": "половой набор", "segment_orientation": "ориентация сегмента",
}


def boot(assets: str = ".") -> dict:
    """Загрузить таблицы цитобендов, примеры, упражнение и справочник участков."""
    global CASES, EXERCISE, REGIONS
    for build, fname in (("GRCh38", "cytoBand_hg38.txt"), ("hg19", "cytoBand_hg19.txt")):
        path = os.path.join(assets, fname)
        BANDS[build] = IF.CytobandTable(path, build=build)
    with open(os.path.join(assets, "cases.json")) as fh:
        CASES = json.load(fh)
    with open(os.path.join(assets, "exercise.json")) as fh:
        EXERCISE = json.load(fh)
    # справочник критических участков: если файла рядом нет, страница работает
    # на встроенном запасе модуля, о чём и сообщает
    path = os.path.join(assets, "critical_regions_v1.csv")
    if os.path.exists(path):
        REGIONS = VB.load_regions(path)
    return {"builds": sorted(BANDS), "cases": len(CASES), "forms": len(IF.FORMS),
            "rules": len(IF.RULES), "exercise": len(EXERCISE),
            "profiles": len(VB.PROFILES), "regions": len(REGIONS),
            "verdicts": list(IV.VERDICTS)}


# ---------------------------------------------------------------------------
# Чтение восстановленных фактов на русском языке
# ---------------------------------------------------------------------------

def _num(n):
    return f"{n:,}".replace(",", "\u202f") if isinstance(n, int) else n


def fact_ru(f) -> str:
    kind = f[0]
    if kind == "cn":
        _, chrom, band, start, end, flanks, cn, cnrange, frac, inh, zyg = f[:11]
        parts = [f"chr{chrom}" + (band or "")]
        if start is not None:
            parts.append(f"{_num(start)}–{_num(end)}")
        if cn is not None:
            parts.append(f"копий {cn}")
        if cnrange:
            parts.append(f"копий {cnrange[0]}~{cnrange[1]}")
        if frac:
            parts.append(f"доля {frac}")
        if zyg:
            parts.append(f"зиготность {zyg}")
        if inh:
            parts.append(f"происхождение {inh}")
        if flanks and any(flanks):
            parts.append("фланки " + "/".join(_num(x) for x in flanks if x))
        return ", ".join(parts)
    if kind == "cnrange_group":
        return f"хромосомы {f[1]}–{f[2]}: копий {f[3]}" + (f", доля {f[4]}" if f[4] else "")
    if kind == "num":
        return f"числовое изменение: {f[2]}{f[1]}"
    if kind in ("chain", "num_chain"):
        terms = f[2] if kind == "chain" else f[2]
        txt = "".join(
            t[0] + "(" + ";".join(t[1]) + ")" +
            ("(" + ";".join("".join(g) for g in t[2]) + ")" if t[2] else "")
            for t in terms)
        pre = "сверхчисленная " if kind == "num_chain" else ""
        tail = f", {f[3]}" if f[3] else ""
        return f"{pre}перестройка: {txt}{tail}"
    if kind == "count":
        return f"число хромосом {f[1]}"
    if kind == "sex":
        return f"половой набор {f[1]}"
    if kind == "complex":
        return f"комплексная перестройка {f[1]}"
    return " ".join(str(x) for x in f)


# ---------------------------------------------------------------------------
# Действия страницы
# ---------------------------------------------------------------------------

def _analyse(case_dict, build="GRCh38"):
    bands = BANDS[build]
    case = IF.case_from_dict(case_dict)
    res = IF.analyse(case, bands)
    forms = []
    for key, f in res["forms"].items():
        audit = []
        if f["text"]:
            try:
                audit = IF.check_conventions(f["text"], bands)
            except Exception as exc:                      # noqa: BLE001
                audit = [{"rule": "—", "severity": "error", "message": str(exc)}]
        forms.append({
            "key": key, "name": f["name_ru"], "format": f["format"],
            "system": f["system"], "text": f["text"],
            "admissible": f["admissible"],
            "roundtrip": RT_RU.get(f["roundtrip"], f["roundtrip"]),
            "roundtrip_code": f["roundtrip"],
            "detail": f["detail"], "derivations": f["derivations"],
            "lost": [LOSS_RU.get(x, x) for x in f["lost_attributes"]],
            "warnings": f["warnings"], "audit": audit,
            "citation": IF.FORMS[key].citation,
            "canonical": key == res["canonical_form"],
        })
    return {"sample": res["sample"], "canonical_form": res["canonical_form"],
            "canonical_text": res["canonical_text"],
            "input_problems": res["input_problems"],
            "omissions": res["karyotype_omissions"], "forms": forms}


def _audit(text, build="GRCh38"):
    bands = BANDS[build]
    out = {"text": text, "build": build}
    try:
        parsed = IF.parse(text)
    except IF.OutOfScope as exc:
        out["status"] = "вне области"
        out["detail"] = str(exc)
        return out
    except IF.ParseError as exc:
        out["status"] = "не разобрано"
        out["detail"] = str(exc)
        return out
    out["status"] = "разобрано"
    out["derivations"] = IF.count_derivations(text)
    out["facts"] = sorted(fact_ru(f) for f in parsed["facts"])
    out["meta"] = {k: v for k, v in parsed["meta"].items() if v}
    try:
        out["violations"] = IF.check_conventions(text, bands)
    except Exception as exc:                              # noqa: BLE001
        out["violations"] = [{"rule": "—", "severity": "error", "message": str(exc)}]
    try:
        out["bands"] = IF.check_band_consistency(text, bands)
    except Exception:                                     # noqa: BLE001
        out["bands"] = []
    return out


def _exercise_check(key, text, build="GRCh38"):
    """Сверка ответа участника со всеми допустимыми формами случая.

    Верным считается совпадение с любой из допустимых форм, а не только с
    канонической: формы равноправны, и участник вправе выбрать любую. Если
    ответ не совпал ни с одной, показывается разбор введённой строки и
    перечень форм, которые для этого случая неприменимы, — чаще всего
    участник пишет именно такую.
    """
    case = next(c for c in CASES if c["key"] == key)
    res = _analyse(case["case"], build)
    got = " ".join((text or "").split())
    # одна и та же строка может быть выдачей нескольких форм — тогда
    # показываются все: это и есть совпадение форм, измеряемое отдельно
    by_text: Dict[str, list] = {}
    for f in res["forms"]:
        if f["admissible"] and f["text"]:
            by_text.setdefault(f["text"], []).append(f)
    refused = [{"key": f["key"], "name": f["name"], "reason": f["detail"]}
               for f in res["forms"] if not f["admissible"]]
    out = {"key": key, "got": got, "canonical": res["canonical_text"],
           "admissible": [{"key": fs[0]["key"], "text": t,
                           "name": ", ".join(f["name"] for f in fs),
                           "keys": [f["key"] for f in fs]}
                          for t, fs in by_text.items()],
           "refused": refused}
    if got in by_text:
        fs = by_text[got]
        out["correct"] = True
        out["matched"] = {"key": fs[0]["key"], "keys": [f["key"] for f in fs],
                          "name": ", ".join(f["name"] for f in fs)}
        return out
    out["correct"] = False
    out["audit"] = _audit(got, build) if got else {"status": "пусто"}
    return out




# ---------------------------------------------------------------------------
# Словесное описание и входной контроль
# ---------------------------------------------------------------------------

STYLE_FIELDS_RU = {
    "include_formula": "приводить формулу", "formula_form": "форма формулы",
    "literal_layer": "слой расшифровки", "significance_layer": "слой значимости",
    "limits_layer": "слой ограничений", "next_steps_layer": "слой возможностей",
    "term_register": "терминологический регистр", "name_syndromes": "называть синдромы",
    "region_overlap": "перекрытие критических участков", "mosaic_policy": "запись мозаицизма",
    "mosaic_uncertainty": "неопределённость доли", "mosaic_thresholds": "границы мозаицизма",
    "sex_policy": "раскрытие полового набора", "coordinates": "координаты в тексте",
    "sizes": "размеры событий", "decimal_comma": "десятичная запятая",
    "bullet_layers": "разделы перечнем", "max_sentences": "предел пояснений",
    "heading_prefix": "приставка заголовка", "counselling_note": "напоминание о консультировании",
    "audience": "получатель", "highlight_conclusion": "выделять вывод",
    "restate_question": "повторять клинический вопрос",
    "transfer_guidance": "указание о пригодности к переносу",
}


def _style_from(req: dict) -> "VB.Style":
    """Профиль плюс точечные правки полей, пришедшие со страницы."""
    st = VB.PROFILES[req.get("profile", "генетик")]
    over = req.get("style") or {}
    fixed = {}
    for k, v in over.items():
        if not hasattr(st, k):
            continue
        cur = getattr(st, k)
        if isinstance(cur, bool):
            fixed[k] = bool(v)
        elif isinstance(cur, tuple):
            if isinstance(v, str):
                # при десятичной запятой разделителем служит точка с запятой
                # или пробел: делить по запятой нельзя — «0,10;0,30» распалось
                # бы на четыре величины
                import re as _re
                parts = [x for x in _re.split(r"[;\s]+", v) if x.strip()]
                v = parts if len(parts) > 1 else [x for x in v.split(",") if x.strip()]
            fixed[k] = tuple(float(str(x).replace(",", ".")) for x in v)
        elif isinstance(cur, int) and not isinstance(cur, bool):
            fixed[k] = None if v in ("", None) else int(v)
        elif isinstance(cur, float):
            fixed[k] = float(str(v).replace(",", "."))
        else:
            fixed[k] = v
    return st.with_(**fixed) if fixed else st


def _verbalize(req: dict) -> dict:
    """Описание по набору событий: разделы, текст, итог обратной проверки."""
    build = req.get("build", "GRCh38")
    bands = BANDS[build]
    case = (next(c for c in CASES if c["key"] == req["key"])["case"]
            if req.get("key") else req["case"])
    if isinstance(case, dict):
        case = IF.case_from_dict(case)
    st = _style_from(req)
    dl = req.get("detection_limit")
    res = VB.verbalize(case, bands, st,
                       detection_limit=(float(dl) if dl not in ("", None) else None),
                       clinical_question=req.get("question") or None,
                       transfer_note=req.get("transfer") or None)
    if REGIONS:
        res["участков_в_справочнике"] = len(REGIONS)
    return {
        "профиль": res["профиль"], "категория": res.get("категория"),
        "разделы": {k: (v if isinstance(v, str) else list(v))
                    for k, v in res["разделы"].items()},
        "текст": res["текст"], "проверка": res["проверка"],
        "запрещённые": [msg for pat, msg in VB.FORBIDDEN
                        if __import__("re").search(pat, res["текст"], __import__("re").I)],
    }


def _validate(req: dict) -> dict:
    """Входной контроль чужой записи: вердикт, замечания, нормализованный вид."""
    v = IV.validate(req["text"], BANDS[req.get("build", "GRCh38")])
    return {"итог": v["итог"], "нормализована": v.get("нормализована"),
            "замены_оформления": v.get("замены_оформления", []),
            "по_уровням": v.get("по_уровням", {}),
            "деревьев_разбора": v.get("деревьев_разбора"),
            "замечания": v.get("замечания", [])}


def _profiles() -> dict:
    """Профили и значения всех осей настройки — для таблицы на странице."""
    import dataclasses as dc
    names = [f.name for f in dc.fields(VB.Style) if f.name != "name"]
    return {"поля": [{"имя": n, "русское": STYLE_FIELDS_RU.get(n, n)} for n in names],
            "профили": {p: {n: (list(getattr(s, n)) if isinstance(getattr(s, n), tuple)
                               else getattr(s, n)) for n in names}
                        for p, s in VB.PROFILES.items()}}


def api(payload: str) -> str:
    """Единственная точка входа со страницы."""
    req = json.loads(payload) if isinstance(payload, str) else payload
    action = req.get("action")
    try:
        if action == "cases":
            data = [{"key": c["key"], "group": c["group"],
                     "description": c["description"]} for c in CASES]
        elif action == "case":
            data = next(c for c in CASES if c["key"] == req["key"])
        elif action == "analyse":
            data = _analyse(req["case"], req.get("build", "GRCh38"))
        elif action == "analyse_key":
            case = next(c for c in CASES if c["key"] == req["key"])
            data = _analyse(case["case"], req.get("build", "GRCh38"))
            data["description"] = case["description"]
        elif action == "audit":
            data = _audit(req["text"], req.get("build", "GRCh38"))
        elif action == "verbalize":
            data = _verbalize(req)
        elif action == "validate":
            data = _validate(req)
        elif action == "profiles":
            data = _profiles()
        elif action == "regions":
            data = {"участков": len(REGIONS),
                    "строки": [{"хромосома": k[0], "начало": k[1], "конец": k[2],
                                "направление": k[3], "обозначение": v}
                               for k, v in sorted(REGIONS.items())]}
        elif action == "exercise":
            data = [{"key": e["key"], "description": e["description"]}
                    for e in EXERCISE]
        elif action == "exercise_check":
            data = _exercise_check(req["key"], req.get("text", ""),
                                   req.get("build", "GRCh38"))
        elif action == "forms":
            data = {"forms": {k: {"name": f.name_ru, "format": f.fmt,
                                  "system": f.system, "citation": f.citation}
                              for k, f in IF.FORMS.items()},
                    "loss": {k: {LOSS_RU.get(a, a): v for a, v in row.items()}
                             for k, row in IF.LOSS_MATRIX.items()}}
        elif action == "rules":
            data = IF.RULES
        else:
            return json.dumps({"error": f"неизвестное действие {action!r}"},
                              ensure_ascii=False)
        return json.dumps({"ok": True, "data": data}, ensure_ascii=False,
                          default=str)
    except Exception as exc:                              # noqa: BLE001
        return json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}"},
                          ensure_ascii=False)
