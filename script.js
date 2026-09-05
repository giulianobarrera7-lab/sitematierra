"use strict";

/* ==========================================================================
   Datos de referencia
   ========================================================================== */

const SOIL_TYPES_SIMPLE = [
  { name: "Agua de mar", rho: 1 },
  { name: "Arcilla húmeda", rho: 30 },
  { name: "Arcilla limosa", rho: 60 },
  { name: "Tierra vegetal", rho: 100 },
  { name: "Arena húmeda", rho: 300 },
  { name: "Grava / piedra", rho: 1500 },
  { name: "Roca / granito", rho: 3000 },
];

// Tabla 771-C.VIII (AEA 90364, Anexo 771-C) — valor más probable y gamas por clima
const SOIL_TABLE_AEA = [
  { name: "Aluvial y arcillas livianas", a: 5, aRange: "—", bRange: "—", c: "1 a 5" },
  { name: "Arcillas (excl. aluvial)", a: 10, aRange: "5 a 20", bRange: "10 a 100", c: "3 a 10" },
  { name: "Greda", a: 20, aRange: "10 a 20", bRange: "50 a 300", c: "3 a 10" },
  { name: "Tierra calcárea porosa", a: 50, aRange: "30 a 100", bRange: "50 a 300", c: "3 a 10" },
  { name: "Arenisca porosa", a: 100, aRange: "30 a 300", bRange: "> 1000", c: "10 a 30" },
  { name: "Cuarzos y piedra caliza compacta y cristalina", a: 300, aRange: "100 a 1000", bRange: "> 1000", c: "30 a 100" },
  { name: "Pizarras arcillosas y esquistos pizarrosos", a: 1000, aRange: "300 a 3000", bRange: "> 1000", c: "30 a 100" },
  { name: "Granito", a: 1000, aRange: "300 a 3000", bRange: "> 1000", c: "30 a 100" },
  { name: "Pizarras rajadizas, rocas ígneas", a: 2000, aRange: "> 1000", bRange: "> 1000", c: "30 a 100" },
];

const SOIL_TABLE_AEA2 = [
  { name: "Terrenos pantanosos", rho: "1 a 30" },
  { name: "Limo", rho: "20 a 100" },
  { name: "Humus", rho: "10 a 150" },
  { name: "Turba húmeda", rho: "5 a 100" },
];

const RE_TABLE = [
  { name: '1/2"', d: 12.6, L: 1.5 },
  { name: '1/2"', d: 12.6, L: 3 },
  { name: '1/2"', d: 12.6, L: 6 },
  { name: '5/8"', d: 14.6, L: 1.5 },
  { name: '5/8"', d: 14.6, L: 3 },
  { name: '5/8"', d: 14.6, L: 6 },
  { name: '3/4"', d: 16.2, L: 1.5 },
  { name: '3/4"', d: 16.2, L: 3 },
  { name: '3/4"', d: 16.2, L: 6 },
];

/* ==========================================================================
   Funciones de cálculo (fórmulas AEA 90364 / Dwight / Sverak / Wenner)
   ========================================================================== */

// Wenner completa (considera la profundidad b de las picas)
function rhoWenner(a, b, R) {
  if (a <= 0 || R < 0) return NaN;
  const denom = 1 + (2 * a) / Math.sqrt(a * a + 4 * b * b) - a / Math.sqrt(a * a + b * b);
  if (denom === 0) return NaN;
  return (4 * Math.PI * a * R) / denom;
}

// Wenner simplificada (perfil de profundidad variable, a >> b)
function rhoWennerSimple(a, R) {
  return 2 * Math.PI * a * R;
}

// Jabalina vertical — fórmula de Dwight. d en mm, L y rho en m/Ω·m
function rodResistance(rho, L, d_mm) {
  const a = d_mm / 2 / 1000; // radio en m
  if (L <= 0 || a <= 0) return NaN;
  return (rho / (2 * Math.PI * L)) * (Math.log((4 * L) / a) - 1);
}

// Radio equivalente Re = L / (ln(4L/a) - 1)
function equivalentRadius(L, d_mm) {
  const a = d_mm / 2 / 1000;
  const denom = Math.log((4 * L) / a) - 1;
  if (denom === 0) return NaN;
  return L / denom;
}

// n jabalinas en paralelo con acoplamiento mutuo
function mutualResistance(rho, s) {
  return rho / (2 * Math.PI * s);
}

function parallelRods(R1, Rm, n) {
  if (n <= 0) return NaN;
  if (n === 1) return R1;
  return (R1 + (n - 1) * Rm) / n;
}

// Malla reticulada — Sverak (IEEE Std 80)
function meshResistanceSverak(rho, LT, A, h) {
  return rho * (1 / LT + (1 / Math.sqrt(20 * A)) * (1 + 1 / (1 + h * Math.sqrt(20 / A))));
}

// Malla reticulada — Laurent (comparación)
function meshResistanceLaurent(rho, LT, A) {
  return (rho / 4) * Math.sqrt(Math.PI / A) + rho / LT;
}

// Conductor horizontal enterrado — fórmula AEA 90364-7 (verificada contra el caso
// publicado en Aquino et al., "Diseño y Cálculo del Sistema de PAT para el
// Instituto de Biotecnología Misiones", UNaM, 2024: con ρ=198,76 Ω·m, L=60 m,
// h=0,9 m y d=6,45 mm da R≈6,72 Ω, igual al valor de esa referencia).
// R = ρ/(2πL) · [ln(2L/a) + ln(L/h) − 2 + 2h/L], con a = radio del conductor.
// Se trunca la serie de corrección en el término lineal (2h/L): los términos
// siguientes (h²/L², h⁴/L⁴, …) son despreciables cuando h ≪ L, como es el caso
// habitual de un conductor de PAT enterrado.
function buriedCableResistance(rho, L, a, h) {
  return (rho / (2 * Math.PI * L)) * (Math.log((2 * L) / a) + Math.log(L / h) - 2 + (2 * h) / L);
}

// Radio equivalente de un cable a partir de su sección (mm²) -> a en m
function cableEquivalentRadius(section_mm2) {
  return Math.sqrt(section_mm2 / Math.PI) / 1000;
}

// Combinación cable horizontal + jabalinas con acoplamiento
function combinedResistance(Rh, Rn, Rm) {
  const denom = Rh + Rn - 2 * Rm;
  if (denom === 0) return NaN;
  return (Rh * Rn - Rm * Rm) / denom;
}

// Sección mínima del PE según Tabla 771.18.III
function peSection(S) {
  if (S <= 16) return S;
  if (S <= 35) return 16;
  return S / 2;
}

function roundTo(x, decimals) {
  const f = 10 ** decimals;
  return Math.round(x * f) / f;
}

function mean(arr) {
  if (!arr.length) return NaN;
  return arr.reduce((a, b) => a + b, 0) / arr.length;
}

function stdev(arr) {
  if (arr.length < 2) return NaN;
  const m = mean(arr);
  const variance = arr.reduce((a, b) => a + (b - m) ** 2, 0) / (arr.length - 1);
  return Math.sqrt(variance);
}

function fmt(x, decimals = 2) {
  if (x === null || x === undefined || Number.isNaN(x) || !Number.isFinite(x)) return "—";
  return x.toLocaleString("es-AR", { minimumFractionDigits: decimals, maximumFractionDigits: decimals });
}

function escapeAttr(str) {
  return String(str).replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function okBadge(isOk, okText = "✅ Cumple", badText = "❌ No cumple") {
  const span = document.createElement("span");
  span.textContent = isOk ? okText : badText;
  span.className = isOk ? "ok" : "bad";
  return span;
}

function setResult(id, value, decimals = 2) {
  const el = document.getElementById(id);
  if (el) el.textContent = fmt(value, decimals);
}

function setBadge(id, isOk, okText, badText) {
  const el = document.getElementById(id);
  if (!el) return;
  el.textContent = "";
  el.appendChild(okBadge(isOk, okText, badText));
}

/* ==========================================================================
   Navegación por pestañas
   ========================================================================== */

function initTabs() {
  const buttons = document.querySelectorAll(".tab-btn");
  buttons.forEach((btn) => {
    btn.addEventListener("click", () => {
      buttons.forEach((b) => b.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach((p) => p.classList.remove("active"));
      btn.classList.add("active");
      document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
    });
  });
}

/* ==========================================================================
   Estado compartido entre pestañas
   ========================================================================== */

const state = {
  profileRows: [
    { a: 2, r: 1 },
    { a: 4, r: 22 },
    { a: 6, r: 19 },
    { a: 8, r: 16 },
    { a: 10, r: 14 },
  ],
  protocolRows: [{ id: "Punto 1", D: 20, dim: 1.5, p52: "", p57: "", p62: "", p67: "", p72: "" }],
};

/* ==========================================================================
   Pestaña 1 — Resistividad y jabalina
   ========================================================================== */

function getDesignRho() {
  const source = document.getElementById("rho_source").value;
  if (source === "manual") {
    return parseFloat(document.getElementById("rho_manual").value) || NaN;
  }
  if (source === "soil") {
    return parseFloat(document.getElementById("rho_soil").value) || NaN;
  }
  if (source === "mediciones") {
    return computeProfileStats().avg;
  }
  // wenner
  const a = parseFloat(document.getElementById("w_a").value);
  const b = parseFloat(document.getElementById("w_b").value);
  const R = parseFloat(document.getElementById("w_r").value);
  return rhoWenner(a, b, R);
}

function getRodDiameter() {
  const sel = document.getElementById("j_d").value;
  if (sel === "custom") {
    return parseFloat(document.getElementById("j_d_custom").value) || NaN;
  }
  return parseFloat(sel);
}

function updateTab1() {
  const a = parseFloat(document.getElementById("w_a").value);
  const b = parseFloat(document.getElementById("w_b").value);
  const R = parseFloat(document.getElementById("w_r").value);
  const rhoW = rhoWenner(a, b, R);
  setResult("w_rho", rhoW, 1);

  const rho = getDesignRho();
  setResult("rho_design", rho, 1);

  const L = parseFloat(document.getElementById("j_L").value);
  const d = getRodDiameter();
  const R1 = rodResistance(rho, L, d);
  setResult("j_r1", R1, 2);

  const s = parseFloat(document.getElementById("j_s").value);
  const n = parseInt(document.getElementById("j_n").value, 10);
  const Rm = mutualResistance(rho, s);
  const Rn = parallelRods(R1, Rm, n);
  setResult("j_rm", Rm, 2);
  setResult("j_rn", Rn, 2);

  const target = parseFloat(document.getElementById("j_target").value);
  let nNeeded;
  if (R1 <= target) {
    nNeeded = 1;
  } else if (target <= Rm) {
    nNeeded = NaN; // no alcanza solo con jabalinas a esa separación
  } else {
    nNeeded = Math.ceil((R1 - Rm) / (target - Rm));
  }
  const el_n = document.getElementById("j_n_needed");
  el_n.textContent = Number.isFinite(nNeeded) ? nNeeded : "— (usar mallado / más separación)";

  let RnNeeded = NaN;
  if (Number.isFinite(nNeeded)) {
    RnNeeded = parallelRods(R1, Rm, nNeeded);
  }
  setResult("j_rn_needed", RnNeeded, 2);
  setBadge("j_meets", Number.isFinite(Rn) && Rn <= target);

  const Re = equivalentRadius(L, d);
  setResult("j_re", Re, 3);
  setResult("j_dmin", 10 * Re, 2);
  setResult("j_smin", 2 * Re, 2);

  // Propaga a la pestaña 2
  const rhoDisplay = document.getElementById("m_rho_display");
  if (rhoDisplay) rhoDisplay.textContent = fmt(rho, 1);
  const rhoWarning = document.getElementById("m_rho_warning");
  if (rhoWarning) rhoWarning.classList.toggle("hidden", !(Number.isFinite(rho) && rho > 2000));
}

function initTab1() {
  const soilSelect = document.getElementById("rho_soil");
  SOIL_TYPES_SIMPLE.forEach((s) => {
    const opt = document.createElement("option");
    opt.value = s.rho;
    opt.textContent = `${s.name} (${s.rho} Ω·m)`;
    soilSelect.appendChild(opt);
  });
  soilSelect.value = 300;

  document.getElementById("rho_source").addEventListener("change", (e) => {
    document.getElementById("rho_manual_wrap").classList.toggle("hidden", e.target.value !== "manual");
    document.getElementById("rho_soil_wrap").classList.toggle("hidden", e.target.value !== "soil");
    updateAll();
  });

  document.getElementById("j_d").addEventListener("change", (e) => {
    document.getElementById("j_d_custom_wrap").classList.toggle("hidden", e.target.value !== "custom");
    updateAll();
  });

  ["w_a", "w_b", "w_r", "rho_manual", "rho_soil", "j_L", "j_d_custom", "j_s", "j_n", "j_target"].forEach((id) => {
    document.getElementById(id).addEventListener("input", updateAll);
  });
}

/* ==========================================================================
   Pestaña 2 — Malla / electrodo
   ========================================================================== */

function updateTab2() {
  const rho = getDesignRho();
  const mode = document.getElementById("m_mode").value;
  document.getElementById("m_modeA_fields").classList.toggle("hidden", mode !== "A");
  document.getElementById("m_modeB_fields").classList.toggle("hidden", mode !== "B");

  const h = parseFloat(document.getElementById("m_h").value);
  const section = parseFloat(document.getElementById("m_section").value);
  const a_cable = cableEquivalentRadius(section);
  const jd = parseFloat(document.getElementById("m_jd").value);

  let area = NaN,
    lhor = NaN,
    ljabtot = 0,
    ltotal = NaN,
    Rtotal = NaN,
    sverak = NaN,
    laurent = NaN,
    cableR = NaN,
    R1 = NaN,
    Rn = NaN;

  document.getElementById("m_row_sverak").classList.toggle("hidden", mode !== "A");
  document.getElementById("m_row_laurent").classList.toggle("hidden", mode !== "A");
  document.getElementById("m_row_cable").classList.toggle("hidden", mode !== "B");
  document.getElementById("m_row_r1").classList.toggle("hidden", mode !== "B");
  document.getElementById("m_row_rn").classList.toggle("hidden", mode !== "B");

  if (mode === "A") {
    const length = parseFloat(document.getElementById("m_length").value);
    const width = parseFloat(document.getElementById("m_width").value);
    const spacing = parseFloat(document.getElementById("m_spacing").value);
    area = length * width;
    const nLong = Math.floor(width / spacing) + 1;
    const nWide = Math.floor(length / spacing) + 1;
    lhor = nLong * length + nWide * width;
    ltotal = lhor;
    sverak = meshResistanceSverak(rho, ltotal, area, h);
    laurent = meshResistanceLaurent(rho, ltotal, area);
    Rtotal = sverak;
  } else {
    const cableLen = parseFloat(document.getElementById("m_cablelen").value);
    const nJab = parseInt(document.getElementById("m_njab").value, 10) || 0;
    const ljab = parseFloat(document.getElementById("m_ljab").value);
    lhor = cableLen;
    ljabtot = nJab * ljab;
    ltotal = lhor + ljabtot;

    cableR = buriedCableResistance(rho, cableLen, a_cable, h);

    if (nJab >= 1) {
      R1 = rodResistance(rho, ljab, jd);
      const s = nJab > 1 ? cableLen / (nJab - 1) : cableLen;
      const Rm_jab = mutualResistance(rho, s);
      Rn = parallelRods(R1, Rm_jab, nJab);
      const d = cableLen / 2 || ljab;
      const Rm_couple = mutualResistance(rho, d || 1);
      Rtotal = combinedResistance(cableR, Rn, Rm_couple);
    } else {
      Rtotal = cableR;
    }
  }

  setResult("m_area", area, 1);
  setResult("m_lhor", lhor, 1);
  setResult("m_ljabtot", ljabtot, 1);
  setResult("m_ltotal", ltotal, 1);
  setResult("m_sverak", sverak, 2);
  setResult("m_laurent", laurent, 2);
  setResult("m_cable", cableR, 2);
  setResult("m_r1", R1, 2);
  setResult("m_rn", Rn, 2);
  setResult("m_rtotal", Rtotal, 2);

  const rangeEl = document.getElementById("m_range");
  if (Number.isFinite(Rtotal)) {
    rangeEl.textContent = `Rango probable (±30%, suelo homogéneo — confirmar con telurímetro): ${fmt(0.7 * Rtotal, 1)} a ${fmt(1.3 * Rtotal, 1)} Ω`;
  } else {
    rangeEl.textContent = "Rango probable (±30%): —";
  }

  const target = parseFloat(document.getElementById("m_target").value);
  setBadge("m_meets", Number.isFinite(Rtotal) && Rtotal <= target);

  // Recomendaciones de protecciones
  const UL = parseFloat(document.getElementById("m_ul").value);
  const sFase = parseFloat(document.getElementById("m_sfase").value);
  const cond = document.getElementById("m_cond").value;

  [0.03, 0.3, 0.5].forEach((idn, i) => {
    const suffix = ["30", "300", "500"][i];
    const uc = Rtotal * idn;
    const ramax = UL / idn;
    setResult(`m_uc${suffix}`, uc, 2);
    setResult(`m_ramax${suffix}`, ramax, 2);
    setBadge(`m_apto${suffix}`, Number.isFinite(Rtotal) && Rtotal <= ramax, "✅ Apto", "❌ No");
  });

  let diffRec = "—";
  if (Number.isFinite(Rtotal)) {
    if (Rtotal <= UL / 0.03) diffRec = "30 mA";
    else if (Rtotal <= UL / 0.3) diffRec = "300 mA";
    else diffRec = "Reducir R (mejorar la puesta a tierra)";
  }
  document.getElementById("m_diff_rec").textContent = diffRec;

  const peSect = peSection(sFase);
  setResult("m_pe_section", peSect, 1);

  let groundSect;
  let groundLabel;
  if (cond === "desnudo") {
    groundSect = 25;
    groundLabel = "Cu desnudo enterrado sin protección contra corrosión: 25 mm² mínimo";
  } else if (cond === "corrosion") {
    groundSect = 16;
    groundLabel = "Con protección contra corrosión y sin protección mecánica: 16 mm² mínimo";
  } else {
    groundSect = Math.max(peSect, 2.5);
    groundLabel = "Protegido contra corrosión y mecánicamente: rige la sección del PE";
  }
  document.getElementById("m_ground_section").textContent = `${fmt(groundSect, 1)} — ${groundLabel}`;

  const sectionOk = section >= groundSect;
  const sectionOkEl = document.getElementById("m_section_ok");
  sectionOkEl.textContent = "";
  sectionOkEl.appendChild(
    okBadge(sectionOk, `✅ Sí — cable de ${section} mm² alcanza (mínimo ${fmt(groundSect, 1)} mm²)`, `❌ No — el mínimo exigido es ${fmt(groundSect, 1)} mm²`)
  );

  const bestGroundSect = Math.max(section, groundSect);
  const cableBuy =
    cond === "desnudo"
      ? `Cu desnudo de ${fmt(bestGroundSect, 1)} mm² para la zanja + cable verde-amarillo de ${fmt(peSect, 1)} mm² para el PE de la distribución`
      : `Cable verde-amarillo de ${fmt(bestGroundSect, 1)} mm² al electrodo + verde-amarillo de ${fmt(peSect, 1)} mm² para el PE de la distribución`;
  document.getElementById("m_cable_buy").textContent = cableBuy;

  // Colores
  document.getElementById("m_color_pe").textContent = `${fmt(peSect, 1)} mm²`;
  document.getElementById("m_color_ground_label").textContent = cond === "desnudo" ? "Cobre desnudo" : "Verde-Amarillo";
  document.getElementById("m_color_ground_section").textContent = `${fmt(bestGroundSect, 1)} mm²`;
  document.getElementById("m_color_n").textContent = `${fmt(sFase, 1)} mm²`;
  document.getElementById("m_color_f").textContent = `${fmt(sFase, 1)} mm²`;
}

function initTab2() {
  [
    "m_mode",
    "m_h",
    "m_section",
    "m_jd",
    "m_length",
    "m_width",
    "m_spacing",
    "m_cablelen",
    "m_njab",
    "m_ljab",
    "m_target",
    "m_ul",
    "m_sfase",
    "m_cond",
  ].forEach((id) => {
    document.getElementById(id).addEventListener("input", updateAll);
    document.getElementById(id).addEventListener("change", updateAll);
  });
}

/* ==========================================================================
   Pestaña 3 — Mediciones
   ========================================================================== */

function renderProfileTable() {
  const tbody = document.getElementById("profile-tbody");
  tbody.innerHTML = "";
  state.profileRows.forEach((row, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${i + 1}</td>
      <td><input type="number" step="0.1" class="prof-a" data-i="${i}" value="${row.a}"></td>
      <td><input type="number" step="0.1" class="prof-r" data-i="${i}" value="${row.r}"></td>
      <td class="prof-rho" data-i="${i}">—</td>
      <td><button class="btn-icon" data-remove="${i}" title="Eliminar">✕</button></td>
    `;
    tbody.appendChild(tr);
  });

  tbody.querySelectorAll(".prof-a, .prof-r").forEach((input) => {
    input.addEventListener("input", (e) => {
      const i = parseInt(e.target.dataset.i, 10);
      const field = e.target.classList.contains("prof-a") ? "a" : "r";
      state.profileRows[i][field] = parseFloat(e.target.value) || 0;
      updateAll();
    });
  });

  tbody.querySelectorAll("[data-remove]").forEach((btn) => {
    btn.addEventListener("click", (e) => {
      const i = parseInt(e.target.dataset.remove, 10);
      state.profileRows.splice(i, 1);
      renderProfileTable();
      updateAll();
    });
  });
}

function computeProfileStats() {
  const rhos = state.profileRows.map((r) => rhoWennerSimple(r.a, r.r)).filter((v) => Number.isFinite(v));
  return {
    avg: mean(rhos),
    max: Math.max(...rhos),
    min: Math.min(...rhos),
    std: stdev(rhos),
  };
}

function updateProfileResults() {
  const tbody = document.getElementById("profile-tbody");
  state.profileRows.forEach((row, i) => {
    const rho = rhoWennerSimple(row.a, row.r);
    const cell = tbody.querySelector(`.prof-rho[data-i="${i}"]`);
    if (cell) cell.textContent = fmt(rho, 1);
  });

  const stats = computeProfileStats();
  setResult("prof_avg", stats.avg, 1);
  document.getElementById("prof_minmax").textContent = `${fmt(stats.max, 1)} / ${fmt(stats.min, 1)}`;
  setResult("prof_std", stats.std, 1);

  const first = state.profileRows[0];
  const last = state.profileRows[state.profileRows.length - 1];
  let trend = "—";
  if (first && last && state.profileRows.length > 1) {
    const rFirst = rhoWennerSimple(first.a, first.r);
    const rLast = rhoWennerSimple(last.a, last.r);
    if (rLast < rFirst) trend = "Decrece (el suelo mejora en profundidad)";
    else if (rLast > rFirst) trend = "Crece (el suelo empeora en profundidad)";
    else trend = "Estable";
  }
  document.getElementById("prof_trend").textContent = trend;

  setResult("med_rho_final", stats.avg, 1);
}

// Dibuja el esquema topográfico E–S–H del método del 62% (picas en línea recta)
function renderMethod62Diagram(D, dmin, positions) {
  const svg = document.getElementById("pm_diagram");
  if (!svg) return;

  if (!Number.isFinite(D) || D <= 0) {
    svg.innerHTML = `<text x="500" y="115" text-anchor="middle" class="diagram-empty">Ingresá D para ver el esquema</text>`;
    return;
  }

  const W = 1000,
    marginL = 70,
    marginR = 70,
    usableW = W - marginL - marginR,
    axisY = 130;
  const xFor = (meters) => marginL + (meters / D) * usableW;

  const xE = xFor(0);
  const xH = xFor(D);
  const x62 = xFor(D * 0.62);
  const xDmin = Number.isFinite(dmin) ? xFor(Math.min(dmin, D)) : null;

  let svgParts = [];

  svgParts.push(`
    <defs>
      <marker id="arrowRed" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
        <path d="M0,0 L6,3 L0,6 Z" fill="#c0392b"></path>
      </marker>
      <marker id="arrowBlue" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
        <path d="M0,0 L6,3 L0,6 Z" fill="#1f5fa8"></path>
      </marker>
    </defs>
  `);

  // Línea de suelo (eje de las 3 picas)
  svgParts.push(`<line x1="${xE}" y1="${axisY}" x2="${xH}" y2="${axisY}" stroke="#9aa8a1" stroke-width="2"></line>`);

  // Distancia mínima recomendada (referencia)
  if (xDmin !== null && dmin < D) {
    svgParts.push(`<line x1="${xDmin}" y1="${axisY - 34}" x2="${xDmin}" y2="${axisY + 34}" stroke="#c9a227" stroke-width="1.5" stroke-dasharray="4,3"></line>`);
    svgParts.push(`<text x="${xDmin}" y="${axisY - 40}" text-anchor="middle" class="diagram-label diagram-label-muted">D mín. ${fmt(dmin, 1)} m</text>`);
  }

  // Flecha de corriente (I) — circuito E -> H, por debajo del eje
  svgParts.push(`<line x1="${xE}" y1="${axisY + 55}" x2="${xH - 4}" y2="${axisY + 55}" stroke="#c0392b" stroke-width="2.5" marker-end="url(#arrowRed)"></line>`);
  svgParts.push(`<text x="${(xE + xH) / 2}" y="${axisY + 72}" text-anchor="middle" class="diagram-label diagram-label-red">I — corriente de ensayo (circuito E↔H)</text>`);

  // Flecha de tensión (V) — E -> S al 62%, por encima del eje
  svgParts.push(`<line x1="${xE}" y1="${axisY - 55}" x2="${x62 - 4}" y2="${axisY - 55}" stroke="#1f5fa8" stroke-width="2.5" marker-end="url(#arrowBlue)"></line>`);
  svgParts.push(`<text x="${(xE + x62) / 2}" y="${axisY - 63}" text-anchor="middle" class="diagram-label diagram-label-blue">V — tensión medida (E↔S)</text>`);

  // Picas S (tensión) en cada posición del ensayo
  positions.forEach((p) => {
    const x = xFor(D * p.pct);
    const isMain = p.pct === 0.62;
    const r = isMain ? 7 : 5;
    svgParts.push(`<line x1="${x}" y1="${axisY - 10}" x2="${x}" y2="${axisY + 10}" stroke="#1f5fa8" stroke-width="1.5"></line>`);
    svgParts.push(`<circle cx="${x}" cy="${axisY}" r="${r}" fill="${isMain ? "#1f5fa8" : "#7fa8d9"}" stroke="#0d3a66" stroke-width="1"></circle>`);
    svgParts.push(`<text x="${x}" y="${axisY + 26}" text-anchor="middle" class="diagram-tick">${Math.round(p.pct * 100)}%</text>`);
    if (isMain) {
      svgParts.push(`<text x="${x}" y="${axisY - 20}" text-anchor="middle" class="diagram-label diagram-label-blue">S 62% ★</text>`);
    }
  });

  // Marcador E (electrodo bajo ensayo)
  svgParts.push(`<circle cx="${xE}" cy="${axisY}" r="9" fill="#1f6f4a" stroke="#164f34" stroke-width="1.5"></circle>`);
  svgParts.push(`<text x="${xE}" y="${axisY + 26}" text-anchor="middle" class="diagram-tick">0 m</text>`);
  svgParts.push(`<text x="${xE}" y="${axisY - 74}" text-anchor="middle" class="diagram-label diagram-label-green">E — bajo ensayo</text>`);

  // Marcador H (pica de corriente, fija, la más lejana)
  svgParts.push(`<circle cx="${xH}" cy="${axisY}" r="9" fill="#c0392b" stroke="#7a2016" stroke-width="1.5"></circle>`);
  svgParts.push(`<text x="${xH}" y="${axisY + 26}" text-anchor="middle" class="diagram-tick">${fmt(D, 1)} m</text>`);
  svgParts.push(`<text x="${xH}" y="${axisY - 74}" text-anchor="middle" class="diagram-label diagram-label-red">H — corriente</text>`);

  svg.innerHTML = svgParts.join("\n");
}

// Dibuja la curva R vs. % de D (la "curva del 62%"): banda de tolerancia +
// lecturas conectadas, para ver de un vistazo si el tramo central es plano.
function renderMethod62Chart(positions, readings, r62, tol, isFlat) {
  const svg = document.getElementById("pm_chart");
  if (!svg) return;

  const finiteReadings = readings.filter((v) => Number.isFinite(v));
  if (finiteReadings.length === 0) {
    svg.innerHTML = `<text x="500" y="140" text-anchor="middle" class="diagram-empty">Cargá las lecturas de R para ver la curva</text>`;
    return;
  }

  const W = 1000,
    H = 280,
    marginL = 70,
    marginR = 40,
    marginT = 30,
    marginB = 55;
  const usableW = W - marginL - marginR;
  const usableH = H - marginT - marginB;

  const pctMin = 0.52,
    pctMax = 0.72;
  const xFor = (pct) => marginL + ((pct - pctMin) / (pctMax - pctMin)) * usableW;

  const hasR62 = Number.isFinite(r62) && r62 !== 0;
  const bandLo = hasR62 ? r62 * (1 - tol) : null;
  const bandHi = hasR62 ? r62 * (1 + tol) : null;

  let yMin = Math.min(...finiteReadings);
  let yMax = Math.max(...finiteReadings);
  if (hasR62) {
    yMin = Math.min(yMin, bandLo);
    yMax = Math.max(yMax, bandHi);
  }
  if (yMin === yMax) {
    yMin -= 1;
    yMax += 1;
  }
  const pad = (yMax - yMin) * 0.2;
  yMin -= pad;
  yMax += pad;

  const yFor = (v) => marginT + (1 - (v - yMin) / (yMax - yMin)) * usableH;

  let parts = [];

  // Ejes
  parts.push(`<line x1="${marginL}" y1="${marginT}" x2="${marginL}" y2="${marginT + usableH}" stroke="#c7d2cb" stroke-width="1.5"></line>`);
  parts.push(`<line x1="${marginL}" y1="${marginT + usableH}" x2="${marginL + usableW}" y2="${marginT + usableH}" stroke="#c7d2cb" stroke-width="1.5"></line>`);

  // Banda de tolerancia alrededor de la lectura al 62%
  if (hasR62) {
    const yLo = yFor(bandHi);
    const yHi = yFor(bandLo);
    const bandColor = isFlat ? "#1f8a5c" : "#c9a227";
    parts.push(`<rect x="${marginL}" y="${yLo}" width="${usableW}" height="${Math.max(yHi - yLo, 1)}" fill="${bandColor}" fill-opacity="0.12"></rect>`);
    parts.push(`<line x1="${marginL}" y1="${yFor(r62)}" x2="${marginL + usableW}" y2="${yFor(r62)}" stroke="${bandColor}" stroke-width="1.5" stroke-dasharray="5,4"></line>`);
    parts.push(`<text x="${marginL + usableW - 4}" y="${yLo - 6}" text-anchor="end" class="diagram-tick">±${Math.round(tol * 100)}% de R(62%)</text>`);
  }

  // Línea que conecta las lecturas disponibles
  const pathPoints = positions
    .map((p, i) => (Number.isFinite(readings[i]) ? `${xFor(p.pct)},${yFor(readings[i])}` : null))
    .filter(Boolean);
  if (pathPoints.length >= 2) {
    parts.push(`<polyline points="${pathPoints.join(" ")}" fill="none" stroke="#1f5fa8" stroke-width="2.5"></polyline>`);
  }

  // Puntos de cada lectura
  positions.forEach((p, i) => {
    const v = readings[i];
    if (!Number.isFinite(v)) return;
    const x = xFor(p.pct);
    const y = yFor(v);
    const isMain = p.pct === 0.62;
    const withinBand = hasR62 ? v >= bandLo && v <= bandHi : true;
    const color = withinBand ? "#1f5fa8" : "#c0392b";
    parts.push(`<circle cx="${x}" cy="${y}" r="${isMain ? 8 : 6}" fill="${color}" stroke="#ffffff" stroke-width="1.5"></circle>`);
    parts.push(`<text x="${x}" y="${y - 14}" text-anchor="middle" class="diagram-tick">${fmt(v, 2)}</text>`);
    parts.push(`<text x="${x}" y="${marginT + usableH + 22}" text-anchor="middle" class="diagram-tick">${Math.round(p.pct * 100)}%</text>`);
  });

  // Etiquetas de eje
  parts.push(`<text x="${marginL - 10}" y="${marginT + 4}" text-anchor="end" class="diagram-tick">${fmt(yMax, 1)}</text>`);
  parts.push(`<text x="${marginL - 10}" y="${marginT + usableH}" text-anchor="end" class="diagram-tick">${fmt(yMin, 1)}</text>`);
  parts.push(`<text x="${marginL + usableW / 2}" y="${H - 8}" text-anchor="middle" class="diagram-tick">Posición de S (% de D)</text>`);
  parts.push(
    `<text x="${marginL}" y="${marginT - 10}" text-anchor="start" class="diagram-label ${isFlat ? "diagram-label-green" : "diagram-label-muted"}">${
      isFlat ? "✅ Curva plana en la zona de tolerancia" : "⚠ Curva todavía no está plana"
    }</text>`
  );

  svg.innerHTML = parts.join("\n");
}

function updateMethod62() {
  const dim = parseFloat(document.getElementById("pm_dim").value);
  const D = parseFloat(document.getElementById("pm_D").value);
  const factor = parseFloat(document.getElementById("pm_factor").value);
  const dmin = dim * factor;
  setResult("pm_dmin", dmin, 1);
  setBadge("pm_dok", D >= dmin, "✅ Sí — lectura confiable", "❌ No — alejar la pica H");

  const positions = [
    { pct: 0.52, posId: "pm_pos52", readId: "pm_r52", devId: "pm_dev52", moveId: "pm_move52" },
    { pct: 0.57, posId: "pm_pos57", readId: "pm_r57", devId: "pm_dev57", moveId: "pm_move57" },
    { pct: 0.62, posId: "pm_pos62", readId: "pm_r62", devId: "pm_dev62", moveId: "pm_move62" },
    { pct: 0.67, posId: "pm_pos67", readId: "pm_r67", devId: "pm_dev67", moveId: "pm_move67" },
    { pct: 0.72, posId: "pm_pos72", readId: "pm_r72", devId: "pm_dev72", moveId: "pm_move72" },
  ];

  renderMethod62Diagram(D, dmin, positions);

  const readings = [];
  let prevMeters = 0;
  positions.forEach((p) => {
    const meters = D * p.pct;
    document.getElementById(p.posId).textContent = fmt(meters, 2);
    document.getElementById(p.moveId).textContent = Number.isFinite(meters) ? fmt(meters - prevMeters, 2) : "—";
    prevMeters = Number.isFinite(meters) ? meters : prevMeters;
    const val = parseFloat(document.getElementById(p.readId).value);
    readings.push(val);
  });

  const r62 = readings[2];
  positions.forEach((p, i) => {
    const dev = Number.isFinite(readings[i]) && r62 ? (readings[i] - r62) / r62 : NaN;
    document.getElementById(p.devId).textContent = Number.isFinite(dev) ? `${(dev * 100).toFixed(1)}%` : "—";
  });

  const avg5 = mean(readings.filter((v) => Number.isFinite(v)));
  setResult("pm_avg5", avg5, 3);

  const devs = readings.filter((v) => Number.isFinite(v) && r62).map((v) => Math.abs((v - r62) / r62));
  const maxDev = devs.length ? Math.max(...devs) : NaN;
  document.getElementById("pm_maxdev").textContent = Number.isFinite(maxDev) ? `${(maxDev * 100).toFixed(1)}%` : "—";

  const tol = parseFloat(document.getElementById("pm_tol").value) / 100;
  const isFlat = Number.isFinite(maxDev) && roundTo(maxDev, 6) <= roundTo(tol, 6);
  setBadge("pm_flat", isFlat, "✅ Curva plana — lectura válida", "❌ Dispersa — alejar H y repetir");

  renderMethod62Chart(positions, readings, r62, tol, isFlat);

  const adoptMode = document.getElementById("pm_adopt").value;
  const adopted = adoptMode === "62" ? r62 : avg5;
  setResult("pm_radopted", adopted, 3);

  const limit = parseFloat(document.getElementById("pm_limit").value);
  setBadge("pm_complies", Number.isFinite(adopted) && adopted <= limit);
}

function initTab3() {
  renderProfileTable();

  document.getElementById("profile-add").addEventListener("click", () => {
    state.profileRows.push({ a: 0, r: 0 });
    renderProfileTable();
    updateAll();
  });

  [
    "pm_dim",
    "pm_D",
    "pm_factor",
    "pm_r52",
    "pm_r57",
    "pm_r62",
    "pm_r67",
    "pm_r72",
    "pm_tol",
    "pm_adopt",
    "pm_limit",
  ].forEach((id) => {
    document.getElementById(id).addEventListener("input", updateAll);
    document.getElementById(id).addEventListener("change", updateAll);
  });
}

/* ==========================================================================
   Pestaña 4 — Protocolo 62% multipunto
   ========================================================================== */

function computeProtocolRow(row) {
  const readings = [row.p52, row.p57, row.p62, row.p67, row.p72].map((v) => parseFloat(v)).filter((v) => Number.isFinite(v));
  const r62 = parseFloat(row.p62);
  const avg = mean(readings);
  let flat = null;
  let adopted = NaN;
  const tol = parseFloat(document.getElementById("pr_tol").value) / 100;

  if (readings.length === 5 && Number.isFinite(r62) && r62 !== 0) {
    const maxDev = Math.max(...readings.map((v) => Math.abs((v - r62) / r62)));
    flat = roundTo(maxDev, 6) <= roundTo(tol, 6);
    adopted = flat ? r62 : avg;
  } else if (readings.length > 0) {
    adopted = avg;
  }

  const dim = parseFloat(row.dim);
  const D = parseFloat(row.D);
  const factor = parseFloat(document.getElementById("pr_factor").value);
  const dOk = Number.isFinite(dim) && Number.isFinite(D) ? D >= dim * factor : null;

  const limit = parseFloat(document.getElementById("pr_limit").value);
  const complies = Number.isFinite(adopted) ? adopted <= limit : null;

  return { adopted, flat, complies, dOk };
}

function renderProtocolTable() {
  const tbody = document.getElementById("protocol-tbody");
  tbody.innerHTML = "";
  state.protocolRows.forEach((row, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${i + 1}</td>
      <td><input type="text" class="pr-id" data-i="${i}" value="${escapeAttr(row.id)}"></td>
      <td><input type="number" step="0.5" class="pr-D" data-i="${i}" value="${row.D}"></td>
      <td><input type="number" step="0.1" class="pr-dim" data-i="${i}" value="${row.dim}"></td>
      <td><input type="number" step="0.01" class="pr-p52" data-i="${i}" value="${row.p52}"></td>
      <td><input type="number" step="0.01" class="pr-p57" data-i="${i}" value="${row.p57}"></td>
      <td><input type="number" step="0.01" class="pr-p62" data-i="${i}" value="${row.p62}"></td>
      <td><input type="number" step="0.01" class="pr-p67" data-i="${i}" value="${row.p67}"></td>
      <td><input type="number" step="0.01" class="pr-p72" data-i="${i}" value="${row.p72}"></td>
      <td class="pr-adopted" data-i="${i}">—</td>
      <td class="pr-flat" data-i="${i}">—</td>
      <td class="pr-complies" data-i="${i}">—</td>
      <td><button class="btn-icon" data-remove="${i}" title="Eliminar">✕</button></td>
    `;
    tbody.appendChild(tr);
  });

  tbody.querySelectorAll("input").forEach((input) => {
    input.addEventListener("input", (e) => {
      const i = parseInt(e.target.dataset.i, 10);
      const cls = e.target.classList;
      if (cls.contains("pr-id")) state.protocolRows[i].id = e.target.value;
      else if (cls.contains("pr-D")) state.protocolRows[i].D = e.target.value;
      else if (cls.contains("pr-dim")) state.protocolRows[i].dim = e.target.value;
      else if (cls.contains("pr-p52")) state.protocolRows[i].p52 = e.target.value;
      else if (cls.contains("pr-p57")) state.protocolRows[i].p57 = e.target.value;
      else if (cls.contains("pr-p62")) state.protocolRows[i].p62 = e.target.value;
      else if (cls.contains("pr-p67")) state.protocolRows[i].p67 = e.target.value;
      else if (cls.contains("pr-p72")) state.protocolRows[i].p72 = e.target.value;
      updateAll();
    });
  });

  tbody.querySelectorAll("[data-remove]").forEach((btn) => {
    btn.addEventListener("click", (e) => {
      const i = parseInt(e.target.dataset.remove, 10);
      state.protocolRows.splice(i, 1);
      renderProtocolTable();
      updateAll();
    });
  });
}

function updateProtocolResults() {
  const tbody = document.getElementById("protocol-tbody");
  const adoptedValues = [];
  let compliesCount = 0;
  let notCompliesCount = 0;
  let dispersedCount = 0;

  state.protocolRows.forEach((row, i) => {
    const result = computeProtocolRow(row);
    const adoptedCell = tbody.querySelector(`.pr-adopted[data-i="${i}"]`);
    const flatCell = tbody.querySelector(`.pr-flat[data-i="${i}"]`);
    const compliesCell = tbody.querySelector(`.pr-complies[data-i="${i}"]`);

    if (adoptedCell) adoptedCell.textContent = fmt(result.adopted, 2);
    if (flatCell) flatCell.textContent = result.flat === null ? "—" : result.flat ? "✅ Sí" : "❌ Dispersa";
    if (compliesCell) compliesCell.textContent = result.complies === null ? "—" : result.complies ? "✅ Sí" : "❌ No";

    if (Number.isFinite(result.adopted)) adoptedValues.push(result.adopted);
    if (result.complies === true) compliesCount++;
    if (result.complies === false) notCompliesCount++;
    if (result.flat === false) dispersedCount++;
  });

  document.getElementById("pr_count").textContent = state.protocolRows.length;

  if (adoptedValues.length) {
    document.getElementById("pr_stats").textContent = `${fmt(Math.min(...adoptedValues), 2)} / ${fmt(mean(adoptedValues), 2)} / ${fmt(Math.max(...adoptedValues), 2)}`;
  } else {
    document.getElementById("pr_stats").textContent = "—";
  }

  document.getElementById("pr_compliance").textContent = `${compliesCount} cumplen / ${notCompliesCount} no cumplen`;

  const statusEl = document.getElementById("pr_status");
  if (adoptedValues.length === 0) {
    statusEl.textContent = "Sin lecturas cargadas todavía.";
  } else if (notCompliesCount === 0 && dispersedCount === 0) {
    statusEl.textContent = "✅ Todos los puntos medidos cumplen — protocolo apto";
  } else {
    statusEl.textContent = `❌ Revisar: ${notCompliesCount} fuera de límite, ${dispersedCount} curvas dispersas`;
  }
}

function initTab4() {
  renderProtocolTable();
  document.getElementById("protocol-add").addEventListener("click", () => {
    state.protocolRows.push({
      id: `Punto ${state.protocolRows.length + 1}`,
      D: 20,
      dim: 1.5,
      p52: "",
      p57: "",
      p62: "",
      p67: "",
      p72: "",
    });
    renderProtocolTable();
    updateAll();
  });

  ["pr_limit", "pr_tol", "pr_factor"].forEach((id) => {
    document.getElementById(id).addEventListener("input", updateAll);
  });
}

/* ==========================================================================
   Pestaña 5 — Referencias (tablas estáticas)
   ========================================================================== */

function initTab5() {
  const simpleBody = document.querySelector("#soil-table-simple tbody");
  SOIL_TYPES_SIMPLE.forEach((s) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${s.name}</td><td>${s.rho}</td>`;
    simpleBody.appendChild(tr);
  });

  const aeaBody = document.querySelector("#soil-table-aea tbody");
  SOIL_TABLE_AEA.forEach((s) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${s.name}</td><td>${s.a}</td><td>${s.aRange}</td><td>${s.bRange}</td><td>${s.c}</td>`;
    aeaBody.appendChild(tr);
  });

  const aea2Body = document.querySelector("#soil-table-aea2 tbody");
  SOIL_TABLE_AEA2.forEach((s) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${s.name}</td><td>${s.rho}</td>`;
    aea2Body.appendChild(tr);
  });

  const reBody = document.querySelector("#re-table tbody");
  RE_TABLE.forEach((r) => {
    const Re = equivalentRadius(r.L, r.d);
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${r.name}</td><td>${r.d}</td><td>${r.L}</td><td>${fmt(10 * Re, 1)}</td>`;
    reBody.appendChild(tr);
  });
}

/* ==========================================================================
   Orquestación general
   ========================================================================== */

function updateAll() {
  updateTab1();
  updateProfileResults();
  updateMethod62();
  updateTab2();
  updateProtocolResults();
}

document.addEventListener("DOMContentLoaded", () => {
  initTabs();
  initTab1();
  initTab2();
  initTab3();
  initTab4();
  initTab5();
  updateAll();
});
