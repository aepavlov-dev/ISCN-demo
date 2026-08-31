/* Демонстрационный стенд генератора формул ISCN 2024.
 *
 * Вся номенклатурная логика — в модуле iscn_formula.py, исполняемом здесь же
 * средой Pyodide. Этот файл только передаёт запросы и раскладывает ответы.
 *
 * Версия среды закреплена: обновление Pyodide может изменить поведение, а
 * стенд должен показывать ровно то, что показала проверка.
 */
const PYODIDE_VERSION = "314.0.6";   // сверено с перечнем версий jsDelivr 31.08.2026
const LARK_WHEEL = "assets/lark-1.3.1-py3-none-any.whl";
const PY_FILES = ["iscn_formula.py", "iscn_cases.py", "iscn_verbalize.py",
                  "iscn_validate.py", "iscn_web.py"];
const DATA_FILES = ["cytoBand_hg38.txt", "cytoBand_hg19.txt", "cases.json",
                    "exercise.json", "critical_regions_v1.csv"];

let pyApi = null;
let STATS = null;

const $ = (sel) => document.querySelector(sel);
const el = (tag, attrs = {}, ...kids) => {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") n.className = v;
    else if (k === "html") n.innerHTML = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined) n.setAttribute(k, v);
  }
  for (const c of kids.flat()) {
    if (c === null || c === undefined || c === false) continue;
    n.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  }
  return n;
};

function progress(pct, msg) {
  $("#bootfill").style.width = pct + "%";
  if (msg) $("#boot").textContent = msg;
}

/* ------------------------------------------------------------------ */
/* Загрузка среды                                                      */
/* ------------------------------------------------------------------ */

const PYODIDE_BASE = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;

/* Загрузчик подключается отсюда, а не из разметки: иначе номер версии
   пришлось бы держать в двух местах, и они разошлись бы. */
function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = src;
    s.onload = resolve;
    s.onerror = () => reject(new Error("не удалось загрузить " + src));
    document.head.appendChild(s);
  });
}

async function boot() {
  try {
    progress(5, "запуск среды Python…");
    await loadScript(PYODIDE_BASE + "pyodide.js");
    if (typeof loadPyodide !== "function") {
      throw new Error("файл загрузчика получен, но функция loadPyodide не " +
        "объявлена: вероятно, версия " + PYODIDE_VERSION + " несовместима. " +
        "Измените PYODIDE_VERSION в assets/stand.js");
    }
    const pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
    progress(35, "установка разборщика…");
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    await micropip.install(new URL(LARK_WHEEL, location.href).href);

    progress(60, "загрузка модуля и таблиц…");
    for (const f of [...PY_FILES, ...DATA_FILES]) {
      const r = await fetch("assets/" + f);
      if (!r.ok) throw new Error(`не найден файл assets/${f} (${r.status})`);
      pyodide.FS.writeFile("/home/pyodide/" + f, new Uint8Array(await r.arrayBuffer()));
    }
    STATS = await (await fetch("assets/stats.json")).json();

    progress(80, "разбор грамматики…");
    pyodide.runPython(`
import sys
sys.path.insert(0, "/home/pyodide")
import iscn_web
iscn_web.boot("/home/pyodide")
`);
    pyApi = pyodide.runPython("iscn_web.api");
    progress(100, "готово");
    $("#bootbar").hidden = true;
    $("#boot").textContent = `среда готова: ${STATS["батарея"]["форм"]} форм, ` +
      `${STATS["источник"]["пунктов_реестра"]} пунктов стандарта`;
    $("#app").hidden = false;
    initTabs();
    await initCases();
    initOwn();
    initAudit();
    await initExercise();
    await initForms();
  await initVerbalize();
  await initValidate();
  await initRegions();
  } catch (err) {
    $("#bootbar").hidden = true;
    $("#boot").textContent = "";
    document.body.appendChild(el("div", { class: "fatal" },
      el("h2", {}, "Среда не загрузилась"),
      el("p", {}, String(err.message || err)),
      el("p", {}, "Страница исполняет Python в браузере и требует доступа к " +
        "cdn.jsdelivr.net при первом открытии. Отчёт и слайды работают без этого: "),
      el("p", {}, el("a", { href: "index.html" }, "слайды"), " · ",
        el("a", { href: "report.html" }, "отчёт"))));
  }
}

function call(req) {
  const raw = pyApi(JSON.stringify(req));
  const res = JSON.parse(raw);
  if (!res.ok) throw new Error(res.error);
  return res.data;
}

/* ------------------------------------------------------------------ */
/* Вкладки                                                             */
/* ------------------------------------------------------------------ */

function initTabs() {
  document.querySelectorAll(".tab").forEach((t) => {
    t.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach((x) => x.classList.remove("active"));
      document.querySelectorAll(".panel").forEach((x) => x.classList.remove("active"));
      t.classList.add("active");
      $("#" + t.dataset.panel).classList.add("active");
    });
  });
}

/* ------------------------------------------------------------------ */
/* Общий вывод: таблица форм                                           */
/* ------------------------------------------------------------------ */

function formsTable(data) {
  const wrap = el("div", {});
  wrap.appendChild(el("div", { class: "canon" },
    el("span", { class: "lbl" }, "Каноническая форма"),
    el("code", { class: "big" }, data.canonical_text || "—")));

  if (data.input_problems && data.input_problems.length) {
    wrap.appendChild(el("div", { class: "warn" },
      el("b", {}, "Замечания к входным данным: "),
      data.input_problems.join("; ")));
  }
  if (data.omissions && data.omissions.length) {
    wrap.appendChild(el("div", { class: "note" },
      el("b", {}, "Выпало из кариотипной строки: "),
      el("ul", {}, data.omissions.map((o) => el("li", {}, o)))));
  }

  const rows = data.forms.map((f) => {
    const status = f.admissible
      ? el("span", { class: f.roundtrip_code === "ok" ? "ok" : "bad" }, f.roundtrip)
      : el("span", { class: "off" }, "неприменима");
    const cell = f.admissible
      ? el("div", {},
          el("code", { class: "res" }, f.text),
          el("button", { class: "copy", title: "скопировать",
            onclick: (e) => {
              navigator.clipboard.writeText(f.text);
              e.target.textContent = "✓";
              setTimeout(() => (e.target.textContent = "копировать"), 1200);
            } }, "копировать"))
      : el("span", { class: "reason" }, f.detail);
    const notes = el("div", { class: "notes" },
      f.derivations !== null && f.derivations !== undefined && f.admissible
        ? el("span", { class: "chip" }, `деревьев разбора: ${f.derivations}`) : null,
      (f.lost || []).length
        ? el("span", { class: "chip lost" }, "теряет: " + f.lost.join(", ")) : null,
      ...(f.warnings || []).map((w) => el("span", { class: "chip warnchip" }, w)),
      ...(f.audit || []).map((a) =>
        el("span", { class: "chip " + (a.severity === "error" ? "errchip" : "stylechip") },
          `${a.rule}: ${a.message}`)));
    return el("tr", { class: f.canonical ? "canonrow" : "" },
      el("td", {}, el("div", { class: "fname" }, (f.canonical ? "★ " : "") + f.name),
        el("div", { class: "fmeta" }, `${f.key} · ${f.citation}`)),
      el("td", {}, cell, notes),
      el("td", {}, status));
  });
  wrap.appendChild(el("table", { class: "forms" },
    el("thead", {}, el("tr", {}, el("th", {}, "Форма"), el("th", {}, "Запись"),
      el("th", {}, "Обратная проверка"))),
    el("tbody", {}, rows)));
  return wrap;
}

/* ------------------------------------------------------------------ */
/* Вкладка 1: примеры                                                  */
/* ------------------------------------------------------------------ */

async function initCases() {
  const cases = call({ action: "cases" });
  const sel = $("#caseSel");
  const groups = {};
  cases.forEach((c) => (groups[c.group] = groups[c.group] || []).push(c));
  for (const [g, items] of Object.entries(groups)) {
    const og = el("optgroup", { label: g });
    items.forEach((c) => og.appendChild(el("option", { value: c.key }, `${c.key} — ${c.description.slice(0, 58)}`)));
    sel.appendChild(og);
  }
  const show = () => {
    const key = sel.value;
    const c = cases.find((x) => x.key === key);
    $("#caseDesc").textContent = c ? c.description : "";
    const data = call({ action: "analyse_key", key });
    $("#caseOut").replaceChildren(formsTable(data));
  };
  sel.addEventListener("change", show);
  sel.value = "NA02325-производная";
  show();
  const hash = decodeURIComponent(location.hash.replace(/^#case=/, ""));
  if (hash && cases.some((c) => c.key === hash)) { sel.value = hash; show(); }
}

/* ------------------------------------------------------------------ */
/* Вкладка 2: свой случай                                              */
/* ------------------------------------------------------------------ */

const CHROMS = [...Array(22).keys()].map((i) => String(i + 1)).concat(["X", "Y"]);
const SYMS = ["", "del", "dup", "inv", "i", "idic", "der", "dic", "t", "ins",
              "rec", "r", "add", "trp", "qdp"];

function eventRow(n) {
  const row = el("div", { class: "event" },
    el("div", { class: "ehead" }, `Событие ${n}`,
      el("button", { class: "rm", onclick: (e) => e.target.closest(".event").remove() }, "убрать")),
    el("div", { class: "grid4" },
      el("label", {}, "Хромосома",
        el("select", { class: "e-chrom" }, CHROMS.map((c) => el("option", {}, c)))),
      el("label", {}, "Масштаб",
        el("select", { class: "e-scale" },
          el("option", { value: "segment" }, "сегмент"),
          el("option", { value: "chromosome" }, "хромосома"))),
      el("label", {}, "Начало", el("input", { class: "e-start", type: "number", placeholder: "1-based" })),
      el("label", {}, "Конец", el("input", { class: "e-end", type: "number" })),
      el("label", {}, "Копийность", el("input", { class: "e-cn", type: "number", value: "3", min: "0", max: "8" })),
      el("label", {}, "Ожидаемая", el("input", { class: "e-base", type: "number", value: "2", min: "1", max: "4" })),
      el("label", {}, "Доля", el("input", { class: "e-frac", placeholder: "0.35 или ?" })),
      el("label", {}, "Происхождение",
        el("select", { class: "e-inh" }, ["", "mat", "pat", "dn", "dmat", "dpat", "umat", "upat"]
          .map((v) => el("option", {}, v)))),
      el("label", {}, "Зиготность",
        el("select", { class: "e-zyg" }, ["", "hmz", "htz"].map((v) => el("option", {}, v)))),
      el("label", {}, "Механизм",
        el("select", { class: "e-sym" }, SYMS.map((v) => el("option", {}, v)))),
      el("label", {}, "Точки разрыва", el("input", { class: "e-bp", placeholder: "p13.3;q11.21" })),
      el("label", { class: "chk" }, el("input", { class: "e-super", type: "checkbox" }), " сверхчисленная")));
  return row;
}

function collectEvents() {
  return [...document.querySelectorAll("#events .event")].map((r) => {
    const g = (c) => r.querySelector(c);
    const num = (c) => (g(c).value === "" ? null : Number(g(c).value));
    const ev = {
      chrom: g(".e-chrom").value,
      scale: g(".e-scale").value,
      copy_number: Number(g(".e-cn").value),
      baseline: Number(g(".e-base").value),
      supernumerary: g(".e-super").checked,
    };
    if (ev.scale === "segment") { ev.start = num(".e-start"); ev.end = num(".e-end"); }
    const frac = g(".e-frac").value.trim();
    if (frac === "?") ev.mosaic_fraction_unknown = true;
    else if (frac) ev.mosaic_fraction = Number(frac.replace(",", "."));
    if (g(".e-inh").value) ev.inheritance = g(".e-inh").value;
    if (g(".e-zyg").value) ev.zygosity = g(".e-zyg").value;
    const sym = g(".e-sym").value;
    if (sym) {
      const bp = g(".e-bp").value.trim();
      const term = { symbol: sym, chroms: [ev.chrom] };
      if (bp) term.breakpoints = bp.split(";").map((s) => s.trim()).filter(Boolean);
      ev.structure = { chain: [term] };
    }
    return ev;
  });
}

function initOwn() {
  let n = 0;
  const add = () => { n += 1; $("#events").appendChild(eventRow(n)); };
  $("#addEvent").addEventListener("click", add);
  add();
  $("#runOwn").addEventListener("click", () => {
    const c = {
      sample: $("#oSample").value || "демонстрация",
      technique: $("#oTech").value,
      sex_chromosomes: $("#oSex").value,
      ploidy: Number($("#oPloidy").value),
      karyotype_available: $("#oKaryo").checked,
      events: collectEvents(),
    };
    if ($("#oZyg").value) c.genome_zygosity = $("#oZyg").value;
    try {
      $("#ownOut").replaceChildren(formsTable(call({ action: "analyse", case: c })));
    } catch (err) {
      $("#ownOut").replaceChildren(el("div", { class: "warn" }, String(err.message)));
    }
  });
}

/* ------------------------------------------------------------------ */
/* Вкладка 3: разбор чужой записи                                      */
/* ------------------------------------------------------------------ */

function auditView(d) {
  const kids = [el("div", { class: "canon" },
    el("span", { class: "lbl" }, "Итог разбора"),
    el("code", { class: "big " + (d.status === "разобрано" ? "ok" : "bad") }, d.status),
    d.derivations ? el("span", { class: "chip" }, `деревьев разбора: ${d.derivations}`) : null)];
  if (d.detail) kids.push(el("div", { class: "warn" }, d.detail));
  if (d.facts) {
    kids.push(el("h3", {}, "Восстановленные факты"),
      el("ul", { class: "facts" }, d.facts.map((f) => el("li", {}, f))));
  }
  if (d.violations) {
    kids.push(el("h3", {}, "Нарушенные пункты стандарта"));
    kids.push(d.violations.length
      ? el("table", { class: "viol" },
          el("thead", {}, el("tr", {}, el("th", {}, "Пункт"), el("th", {}, "Уровень"), el("th", {}, "Замечание"))),
          el("tbody", {}, d.violations.map((v) => el("tr", {},
            el("td", {}, el("code", {}, v.rule)),
            el("td", {}, el("span", { class: v.severity === "error" ? "errchip" : "stylechip" },
              v.severity === "error" ? "ошибка" : "типографика")),
            el("td", {}, v.message)))))
      : el("p", { class: "ok" }, "нарушений не найдено"));
  }
  if (d.bands && d.bands.length) {
    const RU = { exact: "совпадает", coarser: "названа крупнее (допустимо)", inconsistent: "расходится" };
    kids.push(el("h3", {}, "Полосы против координат"),
      el("table", { class: "viol" },
        el("thead", {}, el("tr", {}, el("th", {}, "Хромосома"), el("th", {}, "Заявлено"),
          el("th", {}, "Следует из координат"), el("th", {}, "Итог"))),
        el("tbody", {}, d.bands.map((b) => el("tr", { class: b.relation === "inconsistent" ? "badrow" : "" },
          el("td", {}, "chr" + b.chromosome), el("td", {}, el("code", {}, b.declared)),
          el("td", {}, el("code", {}, b.implied)), el("td", {}, RU[b.relation] || b.relation))))));
  }
  return el("div", {}, kids);
}

function initAudit() {
  const run = () => {
    const text = $("#auditText").value.trim();
    if (!text) return;
    $("#auditOut").replaceChildren(auditView(
      call({ action: "audit", text, build: $("#auditBuild").value })));
  };
  $("#runAudit").addEventListener("click", run);
  $("#auditText").addEventListener("keydown", (e) => {
    if ((e.ctrlKey || e.metaKey) && e.key === "Enter") run();
  });
  document.querySelectorAll(".samples .link").forEach((b) => {
    b.addEventListener("click", () => {
      $("#auditText").value = b.dataset.a;
      $("#auditBuild").value = b.dataset.b;
      run();
    });
  });
}

/* ------------------------------------------------------------------ */
/* Вкладка 4: упражнение                                               */
/* ------------------------------------------------------------------ */

const exState = {};

async function initExercise() {
  const items = call({ action: "exercise" });
  const list = $("#exList");
  items.forEach((it, i) => {
    const out = el("div", { class: "exout" });
    const inp = el("input", { class: "exinp", placeholder: "формула ISCN" });
    const check = () => {
      const d = call({ action: "exercise_check", key: it.key, text: inp.value });
      exState[it.key] = d.correct;
      out.replaceChildren(exResult(d));
      score(items.length);
    };
    inp.addEventListener("keydown", (e) => { if (e.key === "Enter") check(); });
    list.appendChild(el("div", { class: "excase" },
      el("div", { class: "exhead" }, el("b", {}, `${i + 1}. `), it.description,
        el("span", { class: "fmeta" }, " · " + it.key)),
      el("div", { class: "row" }, inp,
        el("button", { onclick: check }, "Проверить")),
      out));
  });
  score(items.length);
}

function exResult(d) {
  if (d.correct) {
    return el("div", { class: "ok box" }, "Верно — это ", el("b", {}, d.matched.name),
      ". Другие допустимые записи того же случая: ",
      el("ul", {}, d.admissible.filter((a) => a.key !== d.matched.key)
        .map((a) => el("li", {}, el("code", {}, a.text), " — " + a.name))));
  }
  const kids = [el("div", { class: "bad" }, "Не совпало ни с одной допустимой формой.")];
  if (d.audit && d.audit.status) {
    kids.push(el("p", {}, "Разбор введённой строки: ",
      el("b", { class: d.audit.status === "разобрано" ? "ok" : "bad" }, d.audit.status),
      d.audit.detail ? " — " + d.audit.detail : ""));
    (d.audit.violations || []).forEach((v) => kids.push(
      el("p", { class: "viol1" }, el("code", {}, v.rule), " ", v.message)));
  }
  kids.push(el("p", {}, "Допустимые записи этого случая:"),
    el("ul", {}, d.admissible.map((a) => el("li", {}, el("code", {}, a.text), " — " + a.name))));
  if (d.refused.length) {
    kids.push(el("details", {}, el("summary", {}, `Формы, неприменимые к этому случаю (${d.refused.length})`),
      el("ul", {}, d.refused.map((r) => el("li", {}, el("b", {}, r.name), " — " + r.reason)))));
  }
  return el("div", { class: "box" }, kids);
}

function score(total) {
  const done = Object.keys(exState).length;
  const right = Object.values(exState).filter(Boolean).length;
  $("#exScore").textContent = done
    ? `Проверено ${done} из ${total}, верно ${right}.`
    : "";
}

/* ------------------------------------------------------------------ */
/* Вкладка 5: формы и потери                                           */
/* ------------------------------------------------------------------ */

async function initForms() {
  const d = call({ action: "forms" });
  const keys = Object.keys(d.forms);
  const attrs = Object.keys(d.loss[keys[0]]);
  const cr = STATS["совпадения"];

  /* Заголовки — короткие обозначения форм горизонтально: повёрнутый текст
     длинных названий не вмещается в ячейку и выходит за пределы таблицы.
     Полные названия вынесены в перечень под таблицей. */
  const head = el("tr", {}, el("th", {}, "Признак события"),
    ...keys.map((k) => el("th", { class: "key" }, el("code", {}, k))));
  const body = attrs.map((a) => el("tr", {}, el("td", {}, a),
    ...keys.map((k) => el("td", { class: d.loss[k][a] ? "yes" : "no" },
      d.loss[k][a] ? "+" : "–"))));
  const adm = el("tr", { class: "rates" }, el("td", {}, "допустима в случаях"),
    ...keys.map((k) => el("td", {}, String(cr[k]["допустимых"]))));
  const rates = el("tr", { class: "rates" }, el("td", {}, "доля совпадений"),
    ...keys.map((k) => el("td", {}, String(cr[k]["доля"]).replace(".", ","))));

  const legend = el("table", { class: "legend" },
    el("thead", {}, el("tr", {}, el("th", {}, "Обозначение"), el("th", {}, "Форма записи"),
      el("th", {}, "Формат"), el("th", {}, "Система"), el("th", {}, "Пункты ISCN"))),
    el("tbody", {}, keys.map((k) => el("tr", {},
      el("td", {}, el("code", {}, k)), el("td", {}, d.forms[k].name),
      el("td", {}, d.forms[k].format), el("td", {}, d.forms[k].system),
      el("td", {}, el("code", {}, d.forms[k].citation))))));

  $("#formsOut").replaceChildren(
    el("div", { class: "scrollx" },
      el("table", { class: "loss" }, el("thead", {}, head),
        el("tbody", {}, [...body, adm, rates]))),
    el("h3", {}, "Обозначения форм"),
    el("div", { class: "scrollx" }, legend));
}

boot();

/* ------------------------------------------------------------------ */
/* Заключение: словесное описание в профилях                           */
/* ------------------------------------------------------------------ */

const VB_FLAGS = ["highlight_conclusion", "restate_question", "transfer_guidance",
                  "region_overlap", "name_syndromes", "mosaic_uncertainty",
                  "coordinates", "sizes", "bullet_layers", "counselling_note"];

async function initVerbalize() {
  const cases = await call({action: "cases"});
  const sel = el("vbCase");
  const groups = {};
  cases.forEach(c => { (groups[c.group] = groups[c.group] || []).push(c); });
  Object.keys(groups).forEach(g => {
    const og = document.createElement("optgroup");
    og.label = g;
    groups[g].forEach(c => {
      const o = document.createElement("option");
      o.value = c.key;
      o.textContent = c.key + " — " + c.description;
      og.appendChild(o);
    });
    sel.appendChild(og);
  });
  sel.selectedIndex = 0;

  const prof = await call({action: "profiles"});
  const ps = el("vbProfile");
  Object.keys(prof["профили"]).forEach(p => {
    const o = document.createElement("option");
    o.value = p; o.textContent = p;
    ps.appendChild(o);
  });
  ps.value = "генетик";

  const fs = el("vbFlags");
  VB_FLAGS.forEach(name => {
    const meta = prof["поля"].find(f => f["имя"] === name);
    const id = "vbf_" + name;
    const lab = document.createElement("label");
    lab.className = "check";
    lab.innerHTML = '<input type="checkbox" id="' + id + '"> ' +
                    (meta ? meta["русское"] : name);
    fs.appendChild(lab);
  });

  function syncFlags() {
    const vals = prof["профили"][ps.value] || {};
    VB_FLAGS.forEach(name => {
      const box = el("vbf_" + name);
      if (box) box.checked = !!vals[name];
    });
    const th = vals["mosaic_thresholds"];
    if (th) el("vbThresholds").value = th.map(x => String(x).replace(".", ",")).join("; ");
  }
  ps.addEventListener("change", syncFlags);
  syncFlags();

  async function run() {
    const style = {mosaic_thresholds: el("vbThresholds").value};
    VB_FLAGS.forEach(name => {
      const box = el("vbf_" + name);
      if (box) style[name] = box.checked;
    });
    const res = await call({
      action: "verbalize", key: sel.value, profile: ps.value, style: style,
      question: el("vbQuestion").value,
      detection_limit: el("vbLimit").value.replace(",", ".")
    });
    if (res && res.error) {
      el("vbVerdict").className = "verdict bad";
      el("vbVerdict").textContent = "ошибка настройки: " + res.error;
      el("vbOut").innerHTML = ""; el("vbCheck").innerHTML = "";
      return;
    }
    const cat = res["категория"] || "—";
    el("vbVerdict").className = "verdict " +
      (res["проверка"]["итог"] === "ок" ? "good" : "bad");
    el("vbVerdict").textContent = "категория значимости: " + cat +
      " · обратная проверка текста: " + res["проверка"]["итог"] +
      (res["запрещённые"].length
        ? " · запрещённых формулировок: " + res["запрещённые"].length : "");
    const order = ["вопрос", "вывод", "формула", "расшифровка", "значимость",
                   "ограничения", "пригодность", "что_дальше"];
    const names = {"вопрос": "Клинический вопрос", "вывод": "Вывод",
                   "формула": "Запись по ISCN 2024", "расшифровка": "Что выявлено",
                   "значимость": "Клиническая значимость", "ограничения": "Ограничения",
                   "пригодность": "Пригодность к переносу",
                   "что_дальше": "Что можно сделать дальше"};
    let html = "";
    order.forEach(k => {
      const body = res["разделы"][k];
      if (!body) return;
      html += "<h3>" + names[k] + "</h3>";
      if (typeof body === "string") html += "<pre class=\"formula\">" + body + "</pre>";
      else html += "<ul>" + body.map(s => "<li>" + s + "</li>").join("") + "</ul>";
    });
    el("vbOut").innerHTML = html;
    const chk = res["проверка"];
    const rows = Object.keys(chk).filter(k => k !== "итог").map(k =>
      [k, Array.isArray(chk[k]) ? (chk[k].join(", ") || "—") : String(chk[k])]);
    el("vbCheck").innerHTML = "<table class=\"grid\"><tbody>" +
      rows.map(r => "<tr><th>" + r[0] + "</th><td>" + r[1] + "</td></tr>").join("") +
      "</tbody></table>";
  }
  el("vbRun").addEventListener("click", run);
  sel.addEventListener("change", run);
  await run();
}

/* ------------------------------------------------------------------ */
/* Контроль записи                                                     */
/* ------------------------------------------------------------------ */

const IV_CLS = {"годна": "good", "годна с замечаниями": "warn",
                "требует решения": "warn", "негодна": "bad",
                "вне области применения": "info"};

async function initValidate() {
  async function run() {
    const res = await call({action: "validate", text: el("ivText").value});
    el("ivVerdict").className = "verdict " + (IV_CLS[res["итог"]] || "info");
    el("ivVerdict").textContent = "вердикт: " + res["итог"] +
      " · деревьев разбора: " + (res["деревьев_разбора"] === null ? "—" : res["деревьев_разбора"]);
    let html = "";
    if (res["нормализована"] && res["нормализована"] !== el("ivText").value)
      html += "<p>приведённый вид: <code>" + res["нормализована"] + "</code></p>";
    if (res["замены_оформления"] && res["замены_оформления"].length)
      html += "<p>замены оформления: " + res["замены_оформления"].join("; ") + "</p>";
    const f = res["замечания"] || [];
    if (f.length) {
      html += "<table class=\"grid\"><thead><tr><th>уровень</th><th>правило</th>" +
              "<th>замечание</th></tr></thead><tbody>" +
              f.map(x => "<tr><td>" + x["уровень"] + "</td><td>" + (x["правило"] || "") +
                         "</td><td>" + x["сообщение"] + "</td></tr>").join("") +
              "</tbody></table>";
    } else {
      html += "<p>замечаний нет.</p>";
    }
    el("ivOut").innerHTML = html;
  }
  el("ivRun").addEventListener("click", run);
  el("ivText").addEventListener("keydown", e => { if (e.key === "Enter") run(); });
  document.querySelectorAll(".ivEx").forEach(a => {
    a.addEventListener("click", e => {
      e.preventDefault();
      el("ivText").value = a.textContent.trim();
      run();
    });
  });
  await run();
}

/* ------------------------------------------------------------------ */
/* Справочник участков                                                 */
/* ------------------------------------------------------------------ */

async function initRegions() {
  const data = await call({action: "regions"});
  const rows = data["строки"] || [];
  const DIR = {loss: "утрата", gain: "прирост"};
  function draw(q) {
    const t = (q || "").trim().toLowerCase();
    const sel = t ? rows.filter(r =>
      (r["обозначение"] + " " + r["хромосома"]).toLowerCase().includes(t)) : rows;
    el("regCount").textContent = "показано " + sel.length + " из " + rows.length +
      " записей уровня «достаточные доказательства»";
    el("regOut").innerHTML = "<table class=\"grid\"><thead><tr><th>хромосома</th>" +
      "<th>границы</th><th>направление</th><th>обозначение</th></tr></thead><tbody>" +
      sel.slice(0, 200).map(r => "<tr><td>" + r["хромосома"] + "</td><td>" +
        r["начало"].toLocaleString("ru-RU") + "–" + r["конец"].toLocaleString("ru-RU") +
        "</td><td>" + (DIR[r["направление"]] || r["направление"]) + "</td><td>" +
        r["обозначение"] + "</td></tr>").join("") + "</tbody></table>";
  }
  el("regFilter").addEventListener("input", e => draw(e.target.value));
  draw("");
}
