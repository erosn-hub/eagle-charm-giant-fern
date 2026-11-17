import {
  ALKYNES,
  B3LYP_CHAIN_DH,
  B3LYP_DEFICIT,
  EXP_DG_ACETYLENE,
  EXP_DH_ACETYLENE,
  METHODS,
  METALS,
  R_KCAL,
  SOLVENTS,
  SUBSTITUENTS,
  type AlkyneId,
  type MetalId,
  type MethodId,
  type SolventId,
  type StandardState,
  type SubstituentId,
  type ThermalMode,
  type ThermoMode,
} from "./data";
import {
  buildGeometry,
  formingCC,
  minDistance,
  stationPieces,
  type GeomKind,
  type XYZ,
} from "./geometry";
import { collectionThermo, eyring, rtLnMolar, wignerKappa } from "./thermo";

export type JobInput = {
  alkynes: [AlkyneId, AlkyneId, AlkyneId];
  metal: MetalId;
  method: MethodId;
  solvent: SolventId;
  temperature: number;
  thermo: ThermoMode;
  standardState: StandardState;
  twoState: boolean;
  thermalMode: ThermalMode;
  wigner: boolean;
};

export type Anchor = "literatura" | "reconstrucción" | "modelo";

export type Station = {
  id: string;
  label: string;
  short: string;
  kind: "min" | "ts" | "mecp";
  dH298: number;
  dH: number;
  dS: number;
  dGatm: number;
  dG: number;
  uncertainty: number;
  anchor: Anchor;
  spin: string;
  imag: number;
  note: string;
  geom: GeomKind;
  nMolecules: number;
  inSpan: boolean;
};

export type Step = {
  id: string;
  name: string;
  dG: number;
  k: number;
  unit: string;
  imag: number;
  kappa: number;
};

export type SpanReport = {
  deltaE: number;
  tdts: string;
  tdi: string;
  tof: number;
  dGr: number;
};

export type Isomer = {
  name: string;
  pattern: string;
  share: number;
  dG: number;
  channel: string;
};

export type Profile = {
  input: JobInput;
  stations: Station[];
  steps: Step[];
  span: SpanReport | null;
  isomers: Isomer[];
  regime: string;
  headline: string;
  dHorganic: number;
  dGorganic: number;
  parent: boolean;
  symmetricReaction: boolean;
  warnings: string[];
};

const DISP: Record<string, number> = {
  pre: 0,
  bis: -2,
  "ts-oc": -2.4,
  mc: -1,
  adduct: -4,
  "ts-collapse": -4.5,
  arene: -7.5,
  "ts-release": -6,
  prod: -1,
  react: 0,
  "thermal-vdw": -3.2,
  "thermal-ts": -1.5,
  "thermal-birad": -1,
  "thermal-prod": -1,
};

const BSSE: Record<string, number> = {
  pre: 0,
  bis: 2.5,
  "ts-oc": 2.6,
  mc: 1,
  adduct: 3.5,
  "ts-collapse": 3.4,
  arene: 2,
  "ts-release": 2,
  prod: 0,
  react: 0,
  "thermal-vdw": 1.5,
  "thermal-ts": 1.2,
  "thermal-birad": 0.6,
  "thermal-prod": 0,
};

const POL: Record<string, number> = {
  pre: 0.35,
  bis: 0.55,
  "ts-oc": 1.25,
  mc: 0.85,
  adduct: 1.05,
  "ts-collapse": 1.15,
  arene: 0.6,
  "ts-release": 0.7,
  prod: 0.4,
  react: 0.2,
  "thermal-vdw": 0.3,
  "thermal-ts": 0.9,
  "thermal-birad": 1.1,
  "thermal-prod": 0.25,
};

function eOf(ids: AlkyneId[]) {
  return ids.reduce((s, id) => {
    const a = ALKYNES[id];
    return s + (SUBSTITUENTS[a.a].sigma + SUBSTITUENTS[a.b].sigma) / 2;
  }, 0);
}
function bOf(id: AlkyneId) {
  const a = ALKYNES[id];
  return (SUBSTITUENTS[a.a].bulk + SUBSTITUENTS[a.b].bulk) / 2;
}
function nPh(ids: AlkyneId[]) {
  return ids.reduce((s, id) => {
    const a = ALKYNES[id];
    return s + (a.a === "Ph" ? 1 : 0) + (a.b === "Ph" ? 1 : 0);
  }, 0);
}
function nEwg(ids: AlkyneId[]) {
  return ids.reduce((s, id) => {
    const a = ALKYNES[id];
    const add = (g: SubstituentId) => (g === "CO2Me" || g === "CN" || g === "CF3" ? 1 : 0);
    return s + add(a.a) + add(a.b);
  }, 0);
}

function orthoStrain(ids: [AlkyneId, AlkyneId, AlkyneId]): number {
  const bulks = ids.flatMap((id) => [SUBSTITUENTS[ALKYNES[id].a].bulk, SUBSTITUENTS[ALKYNES[id].b].bulk]);
  const heavy = bulks.filter((b) => b >= 1.6).length;
  return heavy >= 4 ? 1.6 * (heavy - 3) : 0;
}

export function organicDeltaH(ids: [AlkyneId, AlkyneId, AlkyneId], method: MethodId): number {
  const prod =
    -1.2 * nPh(ids) - 1.6 * nEwg(ids) + orthoStrain(ids) + 0.4 * ids.reduce((s, id) => s + bOf(id), 0);
  return B3LYP_CHAIN_DH + METHODS[method].recovery * B3LYP_DEFICIT + prod;
}

function kirkwood(eps: number) {
  if (eps <= 1.01) return 0;
  return (eps - 1) / (2 * eps + 1);
}

function substModes(g: SubstituentId, stiff: number): number[] {
  const s = stiff;
  switch (g) {
    case "Me":
      return [95 * s, 145 * s];
    case "tBu":
      return [40 * s, 62 * s, 88 * s];
    case "Ph":
      return [42 * s, 75 * s, 160 * s];
    case "CO2Me":
      return [50 * s, 85 * s, 130 * s];
    case "TMS":
      return [32 * s, 55 * s, 80 * s];
    case "CF3":
      return [48 * s, 95 * s];
    case "OMe":
      return [70 * s, 125 * s];
    case "CN":
      return [170 * s];
    default:
      return [];
  }
}

function skeleton(role: string, stiff: number): number[] {
  const m = (list: number[]) => list.map((v) => v * stiff);
  switch (role) {
    case "alkyne":
      return m([520, 640, 730, 780]);
    case "ligand":
      return m([992, 1120, 2320, 2327, 2335]);
    case "arene":
      return m([410, 610, 700, 850, 990]);
    case "complex":
      return m([55, 78, 110, 160, 230, 310, 420]);
    default:
      return m([80, 140, 220]);
  }
}

function modesFor(role: string, atoms: XYZ[], stiff: number, alkyneIds: AlkyneId[]): number[] {
  const base = skeleton(role, 1);
  const subs: SubstituentId[] = [];
  if (role === "alkyne" || role === "arene") {
    for (const id of alkyneIds) {
      subs.push(ALKYNES[id].a, ALKYNES[id].b);
    }
  }
  if (role === "complex") {
    // Los sustituyentes que viajan en el complejo dependen de la estación; se aproximan con los tres alquinos a media rigidez.
    for (const id of alkyneIds) subs.push(ALKYNES[id].a, ALKYNES[id].b);
  }
  const extra = subs.flatMap((g) => substModes(g, stiff));
  void atoms;
  return [...base.map((v) => v * (role === "complex" ? stiff : 1)), ...extra];
}

type Raw = {
  id: string;
  label: string;
  short: string;
  kind: "min" | "ts" | "mecp";
  dH: number;
  anchor: Anchor;
  spin: string;
  imag: number;
  note: string;
  geom: GeomKind;
  inSpan: boolean;
  stiff: number;
};

function solventShift(key: string, eps: number, polarityBoost: number, preKey: string) {
  const f = kirkwood(eps);
  const lam = 2.4 * polarityBoost;
  const abs = (k: string) => -lam * f * (POL[k] ?? 0.4);
  return abs(key) - abs(preKey);
}

function applyCorrections(
  key: string,
  dH: number,
  method: MethodId,
  solvent: SolventId,
  polBoost: number,
  ring: boolean,
  preKey: string,
): number {
  const spec = METHODS[method];
  let e = dH;
  e += (DISP[key] ?? 0) * spec.dispersionWeight;
  if (spec.tz) e += BSSE[key] ?? 0;
  e += solventShift(key, SOLVENTS[solvent].epsilon, polBoost, preKey);
  if (ring) e += spec.recovery * B3LYP_DEFICIT;
  return e;
}

export function compute(input: JobInput): Profile {
  const alkynes: [AlkyneId, AlkyneId, AlkyneId] = input.alkynes;
  const parent =
    alkynes.every((a) => a === "acetylene") && (input.metal === "CpCo" || input.metal === "none");
  const symmetricReaction = alkynes[0] === alkynes[1] && alkynes[1] === alkynes[2];
  const warnings: string[] = [];
  const method = METHODS[input.method];
  if (method.metalWarning && input.metal !== "none") warnings.push(method.metalWarning);
  if (input.metal === "Ni0")
    warnings.push(
      "Ni(PH₃)₂ es un truncamiento de Ni(PMe₃)ₙ. Cambia la energía de disociación de fosfina respecto del precatalizador real Ni(COD)₂.",
    );
  if (!symmetricReaction)
    warnings.push(
      "En la reacción cruzada el perfil dibuja el canal de menor acoplamiento oxidante. El reparto de isómeros está en Regioselectividad.",
    );

  const polBoost = 1 + 0.55 * Math.abs(eOf(alkynes));
  const stations =
    input.metal === "none"
      ? thermalStations(input, polBoost)
      : catalyticStations(input, polBoost);

  const thermoStations = attachThermo(stations, input);
  const steps = elementary(thermoStations, input);
  const span = input.metal === "none" ? null : energySpan(
    thermoStations.filter((s) => s.inSpan),
    input.temperature,
  );
  const isomers = regioselectivity(input, thermoStations);
  const dHorganic = organicDeltaH(alkynes, input.method);
  const dGorganic = freeOrganic(input);

  const rds = [...steps].sort((a, b) => b.dG - a.dG)[0];
  const headline =
    input.metal === "none"
      ? `Sin metal, la barrera concertada sigue siendo prohibitiva a temperatura ambiente. A ${input.temperature.toFixed(0)} K la constante del TS concertado es ${fmtRate(steps[0]?.k ?? 0)}.`
      : span
        ? `Span energético δE = ${span.deltaE.toFixed(1)} kcal/mol. TDTS: ${span.tdts}. TDI: ${span.tdi}. El paso elemental más alto es ${rds?.name ?? "—"} (${rds ? signed(rds.dG) : "—"} kcal/mol).`
        : "Perfil catalítico.";

  return {
    input,
    stations: thermoStations,
    steps,
    span,
    isomers,
    regime: input.metal === "none" ? "orgánica, sin catalizador" : "ciclo catalítico, referencia = precursor + 3 alquinos",
    headline,
    dHorganic,
    dGorganic,
    parent,
    symmetricReaction,
    warnings,
  };
}

function signed(x: number) {
  return `${x >= 0 ? "+" : ""}${x.toFixed(1)}`;
}

function catalyticStations(input: JobInput, polBoost: number): Raw[] {
  const metal = METALS[input.metal];
  const [ia, ib, ic] = input.alkynes;
  const e2 = (eOf([ia]) + eOf([ib])) / 2;
  const s2 = (bOf(ia) + bOf(ib)) / 2;
  const e3 = eOf([ic]);
  const s3 = bOf(ic);
  const prod =
    -1.2 * nPh(input.alkynes) -
    1.6 * nEwg(input.alkynes) +
    orthoStrain(input.alkynes);

  // Ancla Hardesty, referencia precursor + 3 acetilenos = 0.
  const dBis = -11.3 - 4 * e2 + 0.7 * s2;
  const barrier = 12.8 + metal.dBarrier - 8.5 * e2 + 1.6 * Math.max(0, s2 - 0.35);
  const dMc = dBis + (-13.1 + metal.dRxn - 9 * e2 + 1.3 * s2);
  const bind = -12.4 - 5 * e3 + 1.3 * s3 + 0.45 * s2;
  const dAd = dMc + bind;
  const collapse = Math.max(0.25, 0.5 - 0.8 * e3);
  const dTs2 = dAd + collapse;
  const dAr = dAd + (-81.4 + metal.dArene) + prod;
  const dProd = B3LYP_CHAIN_DH + prod;
  const dTs3 = dAr + Math.max(11, dProd - dAr + 8) ;

  const method = input.method;
  const solvent = input.solvent;
  const corr = (key: string, dH: number, ring: boolean) =>
    applyCorrections(key, dH, method, solvent, polBoost, ring, "pre");

  const rows: Raw[] = [
    {
      id: "pre",
      label: `${metal.fragment} + 3 alquinos`,
      short: "Prec",
      kind: "min",
      dH: 0,
      anchor: "literatura",
      spin: "singlete",
      imag: 0,
      note: "Referencia. Entrada al ciclo por sustitución de dos PH₃ (modelo de PPh₃ o CO).",
      geom: "pre",
      inSpan: true,
      stiff: 1,
    },
    {
      id: "bis",
      label: "Complejo bis-η² + alquino libre",
      short: "π₂",
      kind: "min",
      dH: corr("bis", dBis, false),
      anchor: input.metal === "CpCo" && e2 === 0 ? "literatura" : "modelo",
      spin: "singlete",
      imag: 0,
      note: "Hardesty: la doble sustitución PH₃ → acetileno libera 11.3 kcal/mol. Aquí se corrige por sustituyente, base y dispersión.",
      geom: "bis",
      inSpan: true,
      stiff: 1,
    },
    {
      id: "ts1",
      label: "TS de acoplamiento oxidante",
      short: "TS1",
      kind: "ts",
      dH: corr("ts-oc", dBis + barrier, false),
      anchor: input.metal === "CpCo" && e2 === 0 ? "literatura" : "modelo",
      spin: "singlete",
      imag: input.metal === "CpRh" ? 470 : 410,
      note: "Forma el metalaciclopentadieno. En el padre CpCo/acetileno, ΔH‡ = 12.8 kcal/mol desde el bis-complejo (Hardesty, B3LYP).",
      geom: "ts-oc",
      inSpan: true,
      stiff: 1.12,
    },
    {
      id: "mc",
      label: "Metalaciclopentadieno singlete + alquino",
      short: "MC",
      kind: "min",
      dH: corr("mc", dMc, false),
      anchor: input.metal === "CpCo" && e2 === 0 ? "literatura" : "modelo",
      spin: "singlete",
      imag: 0,
      note: "16 electrones si no entra otro ligando. En Co no es el fundamental: el triplete queda por debajo (Agenet, JACS 2007).",
      geom: "mc",
      inSpan: true,
      stiff: 1,
    },
  ];

  if (input.twoState && metal.allowsTwoState) {
    const dMecp = dMc + 10.2 + 0.4 * s2;
    const dTrip = dMc - 15;
    rows.push(
      {
        id: "mecp",
        label: "MECP singlete–triplete",
        short: "MECP",
        kind: "mecp",
        dH: corr("mc", dMecp, false),
        anchor: "modelo",
        spin: "cruce",
        imag: 0,
        note: "Ancla de orden: ~10 kcal/mol sobre el singlete en sistemas CpCo sustituidos (Rodríguez-García, Prieto, García, Chem. Eur. J. 2012). No es un punto estacionario de un solo determinante. Gaussian no optimiza MECP de forma nativa: EasyMECP / sobMECP.",
        geom: "mc",
        inSpan: false,
        stiff: 1,
      },
      {
        id: "mc3",
        label: "Metalaciclo triplete",
        short: "³MC",
        kind: "min",
        dH: corr("mc", dTrip, false),
        anchor: "modelo",
        spin: "triplete",
        imag: 0,
        note: "Pozo de trabajo −15 kcal/mol respecto del singlete (la literatura de sistemas sustituidos llega a −17). Puede ser el estado de reposo y ensanchar el span si hay que volver a cruzar al singlete antes del colapso.",
        geom: "mc",
        inSpan: false,
        stiff: 1,
      },
    );
  }

  rows.push(
    {
      id: "ad",
      label: "Aducto metalaciclo(η²-alquino)",
      short: "Ad",
      kind: "min",
      dH: corr("adduct", dAd, false),
      anchor: input.metal === "CpCo" && e3 === 0 ? "literatura" : "modelo",
      spin: "singlete",
      imag: 0,
      note: "Hardesty: el tercer acetileno se une por solo 12.4 kcal/mol, con el vector C≡C paralelo al eje Co–Cp. Es el origen de la coordinación débil.",
      geom: "adduct",
      inSpan: true,
      stiff: 1,
    },
    {
      id: "ts2",
      label: "TS de colapso [4+2]",
      short: "TS2",
      kind: "ts",
      dH: corr("ts-collapse", dTs2, false),
      anchor: input.metal === "CpCo" ? "literatura" : "modelo",
      spin: "singlete",
      imag: 190,
      note: "Barrera mínima en el padre (ΔH‡ = 0.5). El camino de inserción vía metalacicloheptatrieno está prohibido por simetría en el acetileno (Hardesty). Con alquenos, la inserción sí compite (Gandon, JACS 2006).",
      geom: "ts-collapse",
      inSpan: true,
      stiff: 1.08,
    },
    {
      id: "ar",
      label: "Complejo η⁴-areno",
      short: "η⁴",
      kind: "min",
      dH: corr("arene", dAr, true),
      anchor: input.metal === "CpCo" && prod === 0 ? "reconstrucción" : "modelo",
      spin: "singlete",
      imag: 0,
      note: "Pozo profundo. Hardesty da −81.4 kcal/mol desde el aducto. La recuperación del déficit de reacción de B3LYP se carga aquí porque el anillo ya está formado. Geometría de partida η⁴, no η⁶: el slip es parte de la liberación.",
      geom: "arene",
      inSpan: true,
      stiff: 0.92,
    },
    {
      id: "ts3",
      label: "TS de liberación del areno",
      short: "TS3",
      kind: "ts",
      dH: corr("ts-release", dTs3, true),
      anchor: "modelo",
      spin: "singlete",
      imag: 95,
      note: "No está en el abstract de Hardesty. Se incluye porque el η⁴ puede secuestrar el catalizador. Barrera de trabajo: al menos 11 kcal/mol sobre el complejo, y 8 sobre los productos si el pozo es poco profundo.",
      geom: "ts-release",
      inSpan: true,
      stiff: 1.05,
    },
    {
      id: "prod",
      label: `${metal.fragment} + areno`,
      short: "Prod",
      kind: "min",
      dH: corr("prod", dProd, true),
      anchor: "reconstrucción",
      spin: "singlete",
      imag: 0,
      note: `Cierre del ciclo B3LYP encadenado: ${B3LYP_CHAIN_DH.toFixed(1)} kcal/mol para el acetileno, frente a ΔH° experimental ${EXP_DH_ACETYLENE.toFixed(1)}. La diferencia es el déficit que la recuperación de método intenta cerrar. El catalizador vuelve al precursor.`,
      geom: "prod",
      inSpan: true,
      stiff: 1,
    },
  );
  return rows;
}

function thermalStations(input: JobInput, polBoost: number): Raw[] {
  const ids = input.alkynes;
  const e = eOf(ids) / 3;
  const s = ids.reduce((sum, id) => sum + bOf(id), 0) / 3;
  const ph = nPh(ids);
  const concerted =
    (input.method === "b3lyp" ? 58 : input.method === "m062x" ? 64 : input.method === "wb97xd" ? 63 : 61) -
    4 * e +
    2.4 * s -
    1.5 * ph;
  const birad = 49 - 3 * e + 1.2 * s - 4 * ph;
  const dProd = organicDeltaH(ids, input.method);
  const corr = (key: string, dH: number, ring: boolean) =>
    applyCorrections(key, dH, input.method, input.solvent, polBoost, ring, "react");

  // organicDeltaH ya incluye la recuperación. Evitar aplicarla dos veces.
  const rows: Raw[] = [
    {
      id: "react",
      label: "Tres alquinos separados",
      short: "3 A",
      kind: "min",
      dH: 0,
      anchor: "literatura",
      spin: "singlete",
      imag: 0,
      note: "Referencia de la reacción orgánica. ΔfH° de trabajo: acetileno 54.2, benceno 19.8 kcal/mol.",
      geom: "reactants",
      inSpan: false,
      stiff: 1,
    },
  ];
  if (input.thermalMode === "concerted") {
    rows.push(
      {
        id: "vdw",
        label: "Precomplejo de van der Waals",
        short: "vdW",
        kind: "min",
        dH: corr("thermal-vdw", -3.6 - 0.8 * Math.abs(e), false),
        anchor: "modelo",
        spin: "singlete",
        imag: 0,
        note: "Mínimo débil. En ΔG a 1 M suele desaparecer: la entropía traslacional no perdona un trímero laxo.",
        geom: "thermal-vdw",
        inSpan: false,
        stiff: 0.8,
      },
      {
        id: "ts-c",
        label: "TS concertado D₃h",
        short: "TS",
        kind: "ts",
        dH: corr("thermal-ts", concerted, false),
        anchor: "literatura",
        spin: "singlete",
        imag: 540,
        note: "Bach, Wolber y Schlegel, JACS 1985: ~62 kcal/mol con correlación (MP3); HF/6-31G* ~81. Woodward–Hoffmann la permite; la barrera es repulsión HOMO–HOMO. Santos, Polo y Andrés siguieron el IRC a B3LYP.",
        geom: "thermal-ts",
        inSpan: false,
        stiff: 1.15,
      },
    );
  } else {
    rows.push(
      {
        id: "ts-b",
        label: "TS birradical hacia ciclobutadieno",
        short: "TS-b",
        kind: "ts",
        dH: corr("thermal-birad", birad, false),
        anchor: "modelo",
        spin: "singlete abierto",
        imag: 480,
        note: "Jones y Krebs, Org. Biomol. Chem. 2017: la cascada birradical puede ganar al concertado. La cifra de 49 kcal/mol es un modelo de trabajo, no el valor del artículo. Exige estabilidad de onda (guess=mix) y, idealmente, un single-point multirreferencial.",
        geom: "thermal-birad",
        inSpan: false,
        stiff: 1.1,
      },
      {
        id: "cbd",
        label: "Ciclobutadieno + alquino",
        short: "CBD",
        kind: "min",
        dH: corr("thermal-birad", -6 + 2 * s, false),
        anchor: "modelo",
        spin: "singlete",
        imag: 0,
        note: "Estimación grosera con ΔfH(ciclobutadieno) ~102 kcal/mol: 2 acetileno → C₄H₄ queda cerca de −6 kcal/mol, con incertidumbre de ±12. El antiaromático dimeriza o sigue hacia Dewar benceno.",
        geom: "thermal-birad",
        inSpan: false,
        stiff: 1,
      },
    );
  }
  rows.push({
    id: "prod",
    label: "Areno",
    short: "Areno",
    kind: "min",
    dH: dProd,
    anchor: ids.every((a) => a === "acetylene") ? "literatura" : "modelo",
    spin: "singlete",
    imag: 0,
    note: `ΔH de la reacción orgánica. Experimental del padre: ${EXP_DH_ACETYLENE.toFixed(1)} kcal/mol. El método recupera una fracción del déficit B3LYP (${B3LYP_DEFICIT.toFixed(1)} kcal/mol).`,
    geom: "thermal-prod",
    inSpan: false,
    stiff: 1,
  });
  // Corrige doble conteo: organicDeltaH ya trae recovery. No pasar por applyCorrections.
  void polBoost;
  return rows;
}

function attachThermo(raw: Raw[], input: JobInput): Station[] {
  const method = METHODS[input.method];
  const T = input.temperature;
  const partsOf = (row: Raw) => {
    const pieces = stationPieces({
      kind: row.geom,
      metal: input.metal,
      alkynes: input.alkynes,
    });
    return pieces.map((p) => ({
      atoms: p.atoms,
      sigma: p.sigma,
      modes: modesFor(p.role, p.atoms, row.stiff, input.alkynes),
    }));
  };
  const refParts = partsOf(raw[0]);
  const ref = collectionThermo(refParts, T, input.thermo, method.scale);
  const ref298 = collectionThermo(refParts, 298.15, input.thermo, method.scale);
  const corrSS = input.standardState === "molar" ? rtLnMolar(T) : 0;

  return raw.map((row) => {
    const parts = partsOf(row);
    const now = collectionThermo(parts, T, input.thermo, method.scale);
    const at298 = collectionThermo(parts, 298.15, input.thermo, method.scale);
    const dS298 = at298.S - ref298.S;
    const dCp = at298.Cp - ref298.Cp;
    const dH = row.dH + (dCp * (T - 298.15)) / 1000;
    const dS = now.S - ref.S;
    const dGatm = dH - (T * dS) / 1000;
    const dn = now.n - ref.n;
    const dG = dGatm + dn * corrSS;
    const substU = input.alkynes.filter((a) => a !== "acetylene").length * 0.9;
    const baseU = row.anchor === "literatura" ? 1.5 : row.anchor === "reconstrucción" ? 4 : 6;
    const uncertainty = Math.sqrt(baseU ** 2 + method.uncertainty ** 2 + substU ** 2);
    return {
      id: row.id,
      label: row.label,
      short: row.short,
      kind: row.kind,
      dH298: row.dH,
      dH,
      dS,
      dGatm,
      dG,
      uncertainty,
      anchor: row.anchor,
      spin: row.spin,
      imag: row.imag,
      note: row.note,
      geom: row.geom,
      nMolecules: now.n,
      inSpan: row.inSpan,
    };
  });
}

function elementary(stations: Station[], input: JobInput): Step[] {
  const steps: Step[] = [];
  for (let i = 0; i < stations.length; i++) {
    const s = stations[i];
    if (s.kind !== "ts") continue;
    let prev = stations[0];
    for (let j = i - 1; j >= 0; j--) {
      if (stations[j].kind === "min") {
        prev = stations[j];
        break;
      }
    }
    const dG = s.dG - prev.dG;
    const kappa = input.wigner ? wignerKappa(s.imag, input.temperature) : 1;
    const k = eyring(dG, input.temperature, s.imag, input.wigner);
    steps.push({
      id: s.id,
      name: `${s.short} desde ${prev.short}`,
      dG,
      k,
      unit: "s⁻¹",
      imag: s.imag,
      kappa,
    });
  }
  return steps;
}

function energySpan(cycle: Station[], T: number): SpanReport | null {
  const mins = cycle.filter((s) => s.kind === "min");
  const tss = cycle.filter((s) => s.kind === "ts");
  if (mins.length < 2 || tss.length === 0) return null;
  const dGr = mins[mins.length - 1].dG - mins[0].dG;
  let best = -Infinity;
  let tdts = tss[0].label;
  let tdi = mins[0].label;
  const index = new Map(cycle.map((s, i) => [s.id, i]));
  for (const ts of tss) {
    for (const inter of mins) {
      const iT = index.get(ts.id) ?? 0;
      const iI = index.get(inter.id) ?? 0;
      const wrap = iT < iI;
      const span = ts.dG - inter.dG + (wrap ? dGr : 0);
      if (span > best) {
        best = span;
        tdts = ts.label;
        tdi = inter.label;
      }
    }
  }
  const tof = ((1.380649e-23 * T) / 6.62607015e-34) * Math.exp(-best / (R_KCAL * T));
  return { deltaE: best, tdts, tdi, tof, dGr };
}

function freeOrganic(input: JobInput): number {
  // ΔG del padre experimental se desplaza con el mismo déficit de entalpía; la entropía orgánica se aproxima constante.
  const dH = organicDeltaH(input.alkynes, input.method);
  const dS = -79.8; // cal/mol/K, padre. Los sustituyentes cambian esto; es una cota.
  const dGatm = dH - (input.temperature * dS) / 1000;
  if (input.standardState === "molar") {
    // 3 moléculas → 1, Δn = −2
    return dGatm + -2 * rtLnMolar(input.temperature);
  }
  return dGatm;
}

type End = SubstituentId;

function orientations(id: AlkyneId): { alpha: End; beta: End }[] {
  const d = ALKYNES[id];
  const base = [{ alpha: d.a, beta: d.b }];
  if (d.a !== d.b) base.push({ alpha: d.b, beta: d.a });
  return base;
}

function canon(subs: End[]): string {
  const s = subs.map((g) => SUBSTITUENTS[g].label);
  const forms: string[] = [];
  for (let k = 0; k < 6; k++) {
    const rot = [...s.slice(k), ...s.slice(0, k)];
    forms.push(rot.join("–"));
    forms.push([...rot].reverse().join("–"));
  }
  return forms.sort()[0];
}

function prettyPattern(key: string): string {
  const parts = key.split("–");
  const nonH = parts.map((p, i) => (p === "H" ? null : i + 1)).filter((x): x is number => x !== null);
  if (nonH.length === 0) return "benceno";
  if (nonH.length === 6) return "hexasustituido";
  // Renumerar para que el primer sustituyente sea la posición 1, si todos los no-H son iguales.
  return key;
}

export function regioselectivity(input: JobInput, stations: Station[]): Isomer[] {
  const ids = input.alkynes;
  const homo = ids[0] === ids[1] && ids[1] === ids[2];
  const terminal = ALKYNES[ids[0]].a !== ALKYNES[ids[0]].b;
  if (homo && !terminal) {
    return [
      {
        name: "Un solo regioisómero",
        pattern: `${ALKYNES[ids[0]].formula} × 3 → areno hexasustituido`,
        share: 1,
        dG: 0,
        channel: "Los dos extremos del alquino son iguales: no hay regioselectividad.",
      },
    ];
  }

  const ts = stations.find((s) => s.id === "ts1" || s.id === "ts-c" || s.id === "ts-b");
  const baseG = ts?.dG ?? 20;
  const channels: { key: string; dG: number; text: string }[] = [];

  const thirds = homo ? [0] : [0, 1, 2];
  for (const t of thirds) {
    const pair = [0, 1, 2].filter((i) => i !== t) as [number, number];
    const o1 = orientations(ids[pair[0]]);
    const o2 = orientations(ids[pair[1]]);
    const o3 = orientations(ids[t]);
    for (const a of o1) {
      for (const b of o2) {
        const costAlpha =
          0.85 * SUBSTITUENTS[a.alpha].bulk +
          0.85 * SUBSTITUENTS[b.alpha].bulk -
          1.05 * SUBSTITUENTS[a.alpha].sigma -
          1.05 * SUBSTITUENTS[b.alpha].sigma;
        for (const c of o3) {
          const approach = 0.55 * SUBSTITUENTS[c.alpha].bulk + 0.2 * SUBSTITUENTS[c.beta].bulk;
          const ring = [a.beta, a.alpha, c.alpha, c.beta, b.alpha, b.beta];
          const key = canon(ring);
          const dG = baseG + costAlpha + approach;
          const text = `Acopla primero ${ALKYNES[ids[pair[0]]].formula} + ${ALKYNES[ids[pair[1]]].formula}; el tercero es ${ALKYNES[ids[t]].formula}.`;
          channels.push({ key, dG, text });
        }
      }
    }
  }

  const bestBy = new Map<string, { dG: number; text: string }>();
  for (const ch of channels) {
    const prev = bestBy.get(ch.key);
    if (!prev || ch.dG < prev.dG) bestBy.set(ch.key, { dG: ch.dG, text: ch.text });
  }
  const list = [...bestBy.entries()].map(([key, v]) => ({ key, ...v }));
  const min = Math.min(...list.map((c) => c.dG));
  const weights = list.map((c) => Math.exp((-(c.dG - min)) / (R_KCAL * input.temperature)));
  const z = weights.reduce((s, w) => s + w, 0);
  return list
    .map((c, i) => ({
      name: homo && terminal ? (is124(c.key) ? "1,2,4-trisustituido" : is135(c.key) ? "1,3,5-trisustituido" : "otro patrón") : `Isómero ${i + 1}`,
      pattern: prettyPattern(c.key),
      share: weights[i] / z,
      dG: c.dG - min,
      channel: c.text,
    }))
    .sort((a, b) => b.share - a.share)
    .slice(0, 6);
}

function is124(key: string): boolean {
  const p = key.split("–");
  const idx = p.map((g, i) => (g === "H" ? -1 : i)).filter((i) => i >= 0);
  if (idx.length !== 3) return false;
  const gaps = [idx[1] - idx[0], idx[2] - idx[1], idx[0] + 6 - idx[2]];
  gaps.sort((a, b) => a - b);
  return gaps[0] === 1 && gaps[1] === 2 && gaps[2] === 3;
}
function is135(key: string): boolean {
  const p = key.split("–");
  const idx = p.map((g, i) => (g === "H" ? -1 : i)).filter((i) => i >= 0);
  if (idx.length !== 3) return false;
  const gaps = [idx[1] - idx[0], idx[2] - idx[1], idx[0] + 6 - idx[2]];
  return gaps.every((g) => g === 2);
}

export function methodScan(input: JobInput): { id: MethodId; label: string; dGts: number; dGr: number }[] {
  return (Object.keys(METHODS) as MethodId[]).map((id) => {
    const profile = compute({ ...input, method: id });
    const ts = profile.stations.find((s) => s.kind === "ts");
    const prod = profile.stations[profile.stations.length - 1];
    return {
      id,
      label: METHODS[id].short,
      dGts: ts ? ts.dG - profile.stations[0].dG : 0,
      dGr: prod.dG,
    };
  });
}

export function geometryReport(input: JobInput, geom: GeomKind): {
  atoms: XYZ[];
  minD: number;
  pair: string;
  forming: number | null;
  text: string;
} {
  const atoms = buildGeometry({ kind: geom, metal: input.metal === "none" ? "none" : input.metal, alkynes: input.alkynes });
  const { d, pair } = minDistance(atoms);
  return {
    atoms,
    minD: d,
    pair,
    forming: formingCC(atoms),
    text: atoms
      .map((a) => `${a.el.padEnd(2)} ${a.x.toFixed(6).padStart(11)} ${a.y.toFixed(6).padStart(11)} ${a.z.toFixed(6).padStart(11)}`)
      .join("\n"),
  };
}

export function strainCurve(input: JobInput): { xi: number; e: number; strain: number; interaction: number }[] {
  const metal = METALS[input.metal === "none" ? "CpCo" : input.metal];
  const dE = -13.1 + metal.dRxn;
  const dEa = Math.max(0.8, 12.8 + metal.dBarrier);
  // TS temprano si la reacción es exoérgica.
  const xiTs = dE < 0 ? 0.38 : 0.62;
  const pts: { xi: number; e: number; strain: number; interaction: number }[] = [];
  const strainAtTs = dEa + (input.metal === "CpRh" ? 8 : input.metal === "CpIr" ? 6 : 14);
  for (let i = 0; i <= 24; i++) {
    const xi = i / 24;
    const e = modelIrc(xi, xiTs, dEa, dE);
    const strain = strainAtTs * Math.pow(xi / Math.max(xiTs, 0.2), 1.15);
    pts.push({ xi, e, strain, interaction: e - strain });
  }
  return pts;
}

function modelIrc(xi: number, xiTs: number, dEa: number, dE: number): number {
  // Polinomio con E(0)=0, E(1)=dE, E(xiTs)=dEa, E'(xiTs)=0.
  // E = A xi^3 + B xi^2 + C xi
  const x = xiTs;
  // C + 2 B x + 3 A x^2 = 0
  // A + B + C = dE
  // A x^3 + B x^2 + C x = dEa
  // Solve linearly.
  // From derivative: C = -2 B x - 3 A x^2
  // Plug:
  // A + B + (-2 B x - 3 A x^2) = dE
  // A(1 - 3 x^2) + B(1 - 2 x) = dE
  // A x^3 + B x^2 + x(-2 B x - 3 A x^2) = dEa
  // A x^3 + B x^2 - 2 B x^2 - 3 A x^3 = dEa
  // A (x^3 - 3 x^3) + B (x^2 - 2 x^2) = dEa
  // A (-2 x^3) + B (-x^2) = dEa
  const a11 = 1 - 3 * x * x;
  const a12 = 1 - 2 * x;
  const a21 = -2 * x * x * x;
  const a22 = -(x * x);
  const det = a11 * a22 - a12 * a21 || 1e-8;
  const A = (dE * a22 - a12 * dEa) / det;
  const B = (a11 * dEa - dE * a21) / det;
  const C = -2 * B * x - 3 * A * x * x;
  return A * xi ** 3 + B * xi ** 2 + C * xi;
}

function fmtRate(k: number): string {
  if (!Number.isFinite(k)) return "—";
  if (k === 0) return "0";
  if (k < 1e-6 || k >= 1e6) return k.toExponential(2) + " s⁻¹";
  return k.toPrecision(3) + " s⁻¹";
}

