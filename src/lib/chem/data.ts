/** Constantes, catálogo y anclas bibliográficas. Unidades: kcal/mol, K, cm⁻¹. */

export type SubstituentId =
  | "H"
  | "Me"
  | "Ph"
  | "CO2Me"
  | "CN"
  | "TMS"
  | "CF3"
  | "OMe"
  | "tBu";

export type AlkyneId =
  | "acetylene"
  | "propyne"
  | "phac"
  | "butyne"
  | "dpa"
  | "dmad"
  | "tmsac"
  | "cyanoac"
  | "propiolate"
  | "tbuac";

export type MetalId = "none" | "CpCo" | "CpStarCo" | "CpRh" | "CpIr" | "CpRuCl" | "Ni0";

export type MethodId =
  | "b3lyp"
  | "b3lyp-d3"
  | "m06"
  | "m062x"
  | "wb97xd"
  | "pbe0-d3"
  | "composite";

export type SolventId = "gas" | "toluene" | "thf" | "meCN";
export type ThermoMode = "rrho" | "truhlar" | "grimme";
export type StandardState = "atm" | "molar";
export type ThermalMode = "concerted" | "biradical";

export type Substituent = {
  id: SubstituentId;
  label: string;
  /** Hammett σp aproximado. Positivo = aceptor. */
  sigma: number;
  /** Volumen estérico relativo (H = 0). */
  bulk: number;
  mass: number;
};

export const SUBSTITUENTS: Record<SubstituentId, Substituent> = {
  H: { id: "H", label: "H", sigma: 0, bulk: 0, mass: 1.00784 },
  Me: { id: "Me", label: "Me", sigma: -0.17, bulk: 1.0, mass: 15.0346 },
  Ph: { id: "Ph", label: "Ph", sigma: 0.06, bulk: 1.75, mass: 77.1039 },
  CO2Me: { id: "CO2Me", label: "CO₂Me", sigma: 0.45, bulk: 1.55, mass: 59.044 },
  CN: { id: "CN", label: "CN", sigma: 0.66, bulk: 0.55, mass: 26.0174 },
  TMS: { id: "TMS", label: "TMS", sigma: -0.04, bulk: 2.55, mass: 73.188 },
  CF3: { id: "CF3", label: "CF₃", sigma: 0.54, bulk: 1.35, mass: 69.006 },
  OMe: { id: "OMe", label: "OMe", sigma: -0.27, bulk: 1.15, mass: 31.0339 },
  tBu: { id: "tBu", label: "tBu", sigma: -0.2, bulk: 2.85, mass: 57.115 },
};

export type AlkyneDef = {
  id: AlkyneId;
  name: string;
  formula: string;
  a: SubstituentId;
  b: SubstituentId;
};

export const ALKYNES: Record<AlkyneId, AlkyneDef> = {
  acetylene: { id: "acetylene", name: "Acetileno", formula: "HC≡CH", a: "H", b: "H" },
  propyne: { id: "propyne", name: "Propino", formula: "MeC≡CH", a: "Me", b: "H" },
  phac: { id: "phac", name: "Fenilacetileno", formula: "PhC≡CH", a: "Ph", b: "H" },
  butyne: { id: "butyne", name: "2-Butino", formula: "MeC≡CMe", a: "Me", b: "Me" },
  dpa: { id: "dpa", name: "Difenilacetileno", formula: "PhC≡CPh", a: "Ph", b: "Ph" },
  dmad: { id: "dmad", name: "DMAD", formula: "MeO₂CC≡CCO₂Me", a: "CO2Me", b: "CO2Me" },
  tmsac: { id: "tmsac", name: "TMS-acetileno", formula: "TMSC≡CH", a: "TMS", b: "H" },
  cyanoac: { id: "cyanoac", name: "Cianoacetileno", formula: "NCC≡CH", a: "CN", b: "H" },
  propiolate: { id: "propiolate", name: "Propiolato de metilo", formula: "MeC≡CCO₂Me", a: "Me", b: "CO2Me" },
  tbuac: { id: "tbuac", name: "3,3-Dimetil-1-butino", formula: "tBuC≡CH", a: "tBu", b: "H" },
};

export const ALKYNE_LIST = Object.values(ALKYNES);

export type MetalDef = {
  id: MetalId;
  name: string;
  fragment: string;
  note: string;
  /** Desplazamiento de ΔH‡ del acoplamiento oxidante respecto de CpCo. */
  dBarrier: number;
  /** Desplazamiento de ΔH del acoplamiento (MC − bis). Positivo = menos exoérgico. */
  dRxn: number;
  /** Estabilización extra del complejo η⁴-areno (negativo = más unido). */
  dArene: number;
  /** Masa del fragmento metálico sin los dos PH₃ ni los alquinos, g/mol. */
  fragmentMass: number;
  element: string;
  allowsTwoState: boolean;
};

export const METALS: Record<MetalId, MetalDef> = {
  none: {
    id: "none",
    name: "Térmica, sin metal",
    fragment: "—",
    note: "Cicloadición [2+2+2] no catalizada. Permitida por Woodward–Hoffmann, pero con una barrera alta por repulsión de capas cerradas.",
    dBarrier: 0,
    dRxn: 0,
    dArene: 0,
    fragmentMass: 0,
    element: "",
    allowsTwoState: false,
  },
  CpCo: {
    id: "CpCo",
    name: "CpCo (ancla)",
    fragment: "CpCo(PH₃)₂",
    note: "Superficie ancla de Hardesty, Albright et al., JACS 1999. PH₃ trunca PPh₃/CO. Co(I), d⁸.",
    dBarrier: 0,
    dRxn: 0,
    dArene: 0,
    fragmentMass: 123.99,
    element: "Co",
    allowsTwoState: true,
  },
  CpStarCo: {
    id: "CpStarCo",
    name: "Cp*Co",
    fragment: "Cp*Co(PH₃)₂",
    note: "Cp* más dador y más voluminoso. Barrera de acoplamiento algo mayor; el areno se libera con más facilidad. Modelo, no un artículo único.",
    dBarrier: 2.0,
    dRxn: 3.0,
    dArene: 3.5,
    fragmentMass: 194.13,
    element: "Co",
    allowsTwoState: true,
  },
  CpRh: {
    id: "CpRh",
    name: "CpRh",
    fragment: "CpRh(PH₃)₂",
    note: "ΔE de reacción B3LYP/TZ2P de Dalla Tiezza, Bickelhaupt y Orian (ChemPhysChem 2018): Co −16.5, Rh +4.2 kcal/mol. La barrera se traslada con Evans–Polanyi (α = 0.45).",
    dBarrier: 9.3,
    dRxn: 20.68,
    dArene: -6,
    fragmentMass: 167.98,
    element: "Rh",
    allowsTwoState: false,
  },
  CpIr: {
    id: "CpIr",
    name: "CpIr",
    fragment: "CpIr(PH₃)₂",
    note: "Mismo análisis de activación-distorsión: Ir queda entre Co y Rh (ΔE ≈ +1.0 kcal/mol en B3LYP/TZ2P para el acoplamiento).",
    dBarrier: 7.9,
    dRxn: 17.53,
    dArene: -4,
    fragmentMass: 257.29,
    element: "Ir",
    allowsTwoState: false,
  },
  CpRuCl: {
    id: "CpRuCl",
    name: "CpRuCl",
    fragment: "CpRuCl(PH₃)₂",
    note: "Modelo de los sistemas CpRuCl usados en [2+2+2]. El desplazamiento frente a Co es una tendencia, no una barrera publicada para este truncamiento.",
    dBarrier: 3.5,
    dRxn: 4,
    dArene: -3,
    fragmentMass: 201.63,
    element: "Ru",
    allowsTwoState: false,
  },
  Ni0: {
    id: "Ni0",
    name: "Ni(0)",
    fragment: "Ni(PH₃)₂",
    note: "Truncamiento de Ni(PMe₃)ₙ / Ni(COD)₂. El Ni(0) acopla con barrera baja y suelta el areno antes que CpCo. Tendencia de modelo.",
    dBarrier: -4.5,
    dRxn: -8,
    dArene: 6,
    fragmentMass: 58.693,
    element: "Ni",
    allowsTwoState: false,
  },
};

export const METAL_LIST = Object.values(METALS);

export type MethodDef = {
  id: MethodId;
  label: string;
  short: string;
  functional: string;
  optBasis: string;
  spBasis: string;
  dispersionKeyword: string | null;
  /** Fracción del déficit de −28.5 kcal/mol del ciclo B3LYP reconstruido frente a ΔH° experimental. */
  recovery: number;
  /** true si la base del perfil es triple zeta (corrección BSSE). */
  tz: boolean;
  /** Fracción de la corrección de dispersión D3 de trabajo. */
  dispersionWeight: number;
  scale: number;
  uncertainty: number;
  metalWarning: string | null;
};

export const METHODS: Record<MethodId, MethodDef> = {
  b3lyp: {
    id: "b3lyp",
    label: "B3LYP/6-31G(d) · gas, sin dispersión",
    short: "B3LYP/DZ",
    functional: "B3LYP",
    optBasis: "6-31G(d)",
    spBasis: "6-311+G(d,p)",
    dispersionKeyword: null,
    recovery: 0,
    tz: false,
    dispersionWeight: 0,
    scale: 0.9614,
    uncertainty: 5,
    metalWarning: null,
  },
  "b3lyp-d3": {
    id: "b3lyp-d3",
    label: "B3LYP-D3(BJ)/def2-TZVP",
    short: "B3LYP-D3/TZ",
    functional: "B3LYP",
    optBasis: "def2SVP",
    spBasis: "def2TZVP",
    dispersionKeyword: "empiricaldispersion=GD3BJ",
    recovery: 0.25,
    tz: true,
    dispersionWeight: 1,
    scale: 0.96,
    uncertainty: 3.5,
    metalWarning: null,
  },
  m06: {
    id: "m06",
    label: "M06/def2-TZVP",
    short: "M06/TZ",
    functional: "M06",
    optBasis: "def2SVP",
    spBasis: "def2TZVP",
    dispersionKeyword: null,
    recovery: 0.65,
    tz: true,
    dispersionWeight: 0.55,
    scale: 0.97,
    uncertainty: 3,
    metalWarning: null,
  },
  m062x: {
    id: "m062x",
    label: "M06-2X/def2-TZVP",
    short: "M06-2X/TZ",
    functional: "M062X",
    optBasis: "def2SVP",
    spBasis: "def2TZVP",
    dispersionKeyword: null,
    recovery: 0.85,
    tz: true,
    dispersionWeight: 0.35,
    scale: 0.97,
    uncertainty: 2.4,
    metalWarning:
      "M06-2X (54 % de intercambio HF) está parametrizado para grupos principales. En Co/Rh/Ru/Ir sube el error: mejor M06, ωB97X-D o PBE0-D3 para el metal, y reservar M06-2X a la reacción térmica orgánica.",
  },
  wb97xd: {
    id: "wb97xd",
    label: "ωB97X-D/def2-TZVP",
    short: "ωB97X-D/TZ",
    functional: "wB97XD",
    optBasis: "def2SVP",
    spBasis: "def2TZVP",
    dispersionKeyword: null,
    recovery: 0.75,
    tz: true,
    dispersionWeight: 1,
    scale: 0.975,
    uncertainty: 2.8,
    metalWarning: null,
  },
  "pbe0-d3": {
    id: "pbe0-d3",
    label: "PBE0-D3(BJ)/def2-TZVP",
    short: "PBE0-D3/TZ",
    functional: "PBE1PBE",
    optBasis: "def2SVP",
    spBasis: "def2TZVP",
    dispersionKeyword: "empiricaldispersion=GD3BJ",
    recovery: 0.6,
    tz: true,
    dispersionWeight: 1,
    scale: 0.959,
    uncertainty: 3,
    metalWarning: null,
  },
  composite: {
    id: "composite",
    label: "Compuesto: ωB97X-D/def2-SVP opt+freq → SP def2-TZVPP + SMD",
    short: "Compuesto ωB97X-D",
    functional: "wB97XD",
    optBasis: "def2SVP",
    spBasis: "def2TZVPP",
    dispersionKeyword: null,
    recovery: 0.8,
    tz: true,
    dispersionWeight: 1,
    scale: 0.975,
    uncertainty: 2.5,
    metalWarning: null,
  },
};

export const METHOD_LIST = Object.values(METHODS);

export type SolventDef = { id: SolventId; name: string; epsilon: number; smd: string | null };

export const SOLVENTS: Record<SolventId, SolventDef> = {
  gas: { id: "gas", name: "Fase gas", epsilon: 1, smd: null },
  toluene: { id: "toluene", name: "Tolueno (SMD)", epsilon: 2.374, smd: "toluene" },
  thf: { id: "thf", name: "THF (SMD)", epsilon: 7.426, smd: "thf" },
  meCN: { id: "meCN", name: "Acetonitrilo (SMD)", epsilon: 35.688, smd: "acetonitrile" },
};

/** ΔH° experimental 3 acetileno → benceno, kcal/mol. */
export const EXP_DH_ACETYLENE = -142.8;
/** ΔG° experimental aproximado, kcal/mol, 298.15 K, 1 atm. */
export const EXP_DG_ACETYLENE = -119.6;

/**
 * Déficit del ciclo B3LYP reconstruido (Hardesty) frente al experimento:
 * el encadenado de ΔH publicados cierra en ≈ −114.3, no en −142.8.
 */
export const B3LYP_CHAIN_DH = -114.3;
export const B3LYP_DEFICIT = EXP_DH_ACETYLENE - B3LYP_CHAIN_DH; // −28.5

export const R_KCAL = 0.00198720425864083;
export const ANGSTROM_U_TO_KGM2 = 1.6605390666e-47;

export type Citation = {
  id: string;
  authors: string;
  year: number;
  title: string;
  source: string;
  doi: string;
  usedFor: string;
};

export const CITATIONS: Citation[] = [
  {
    id: "hardesty1999",
    authors: "Hardesty, J. H.; Koerner, J. B.; Albright, T. A.; Lee, G.-Y.",
    year: 1999,
    title: "Theoretical Study of the Acetylene Trimerization with CpCo",
    source: "J. Am. Chem. Soc. 121, 6055–6067",
    doi: "10.1021/ja983098e",
    usedFor:
      "Ancla de ΔH del ciclo CpCo(PH₃)₂: sustitución de ligandos −11.3; acoplamiento oxidante ΔH‡ 12.8 y ΔH −13.1; tercer acetileno unido por 12.4; colapso [4+2] ΔH‡ 0.5 y ΔH −81.4 kcal/mol.",
  },
  {
    id: "agenet2007",
    authors: "Agenet, N.; Gandon, V.; Vollhardt, K. P. C.; Malacria, M.; Aubert, C.",
    year: 2007,
    title: "Cobalt-Catalyzed Cyclotrimerization of Alkynes: The Answer to the Puzzle of Parallel Reaction Pathways",
    source: "J. Am. Chem. Soc. 129, 8860–8871",
    doi: "10.1021/ja072208r",
    usedFor:
      "El metalaciclopentadieno de cobalto es un triplete. Dos ciclos: [4+2] intramolecular asistido por el metal, o atrapamiento con PR₃ y Diels–Alder intermolecular vía cobaltanorborneno. Cruce singlete–triplete por MECP.",
  },
  {
    id: "gandon2006",
    authors: "Gandon, V.; Agenet, N.; Vollhardt, K. P. C.; Malacria, M.; Aubert, C.",
    year: 2006,
    title: "Cobalt-Mediated Cyclic and Linear 2:1 Cooligomerization of Alkynes with Alkenes: A DFT Study",
    source: "J. Am. Chem. Soc. 128, 8509–8520",
    doi: "10.1021/ja060756j",
    usedFor:
      "Contraste mecanístico: el alqueno se inserta en el enlace Co–C y da un metalaciclo de siete, a diferencia del colapso directo del alquino.",
  },
  {
    id: "roglans2021",
    authors: "Roglans, A.; Pla-Quintana, A.; Solà, M.",
    year: 2021,
    title: "Mechanistic Studies of Transition-Metal-Catalyzed [2+2+2] Cycloaddition Reactions",
    source: "Chem. Rev. 121, 1894–1979",
    doi: "10.1021/acs.chemrev.0c00062",
    usedFor:
      "Marco del campo: Co, Ni, Ru, Rh, Ir y otros; regioselectividad 1,2,4 frente a 1,3,5; el complejo de areno como posible estado de reposo.",
  },
  {
    id: "oriana2018",
    authors: "Dalla Tiezza, M.; Bickelhaupt, F. M.; Orian, L.",
    year: 2018,
    title: "Group 9 Metallacyclopentadienes as Key Intermediates in [2+2+2] Alkyne Cyclotrimerizations. Insight from Activation Strain Analyses",
    source: "ChemPhysChem 19, 1766–1773",
    doi: "10.1002/cphc.201800178",
    usedFor:
      "El acoplamiento oxidante es el paso de mayor barrera del ciclo CpM. ΔE B3LYP/TZ2P: Co −16.53, Rh +4.15, Ir +1.00 kcal/mol. El Co gana por una interacción más estabilizante.",
  },
  {
    id: "bach1985",
    authors: "Bach, R. D.; Wolber, G. J.; Schlegel, H. B.",
    year: 1985,
    title: "The Origin of the Barriers to Thermally Allowed, Six-Electron, Pericyclic Reactions: The Effect of HOMO–HOMO Interactions on the Trimerization of Acetylene",
    source: "J. Am. Chem. Soc. 107, 2837–2841",
    doi: "10.1021/ja00296a002",
    usedFor:
      "Barrera concertada ≈ 62 kcal/mol al incluir correlación (HF/6-31G* ≈ 81). ΔH° ≈ −143 kcal/mol. La barrera es repulsión de capas cerradas, no una prohibición de simetría.",
  },
  {
    id: "santos2005",
    authors: "Santos, J. C.; Polo, V.; Andrés, J.",
    year: 2005,
    title: "An electron localization function study of the trimerization of acetylene: Reaction mechanism and development of aromaticity",
    source: "Chem. Phys. Lett. 406, 393–397",
    doi: "10.1016/j.cplett.2005.02.102",
    usedFor:
      "IRC en D₃h a B3LYP/6-311+G(d) (Gaussian 98). La aromaticidad de plano σ y π aparece a lo largo del IRC, no de golpe en el TS.",
  },
  {
    id: "obc2017",
    authors: "Jones, G. O.; Krebs, Z. J.",
    year: 2017,
    title: "The trimerization of acetylenes involves a cascade of biradical and pericyclic processes",
    source: "Org. Biomol. Chem. 15, 8326–8333",
    doi: "10.1039/c7ob01885a",
    usedFor:
      "Para acetilenos neutros o ricos, la cascada birradical → ciclobutadieno → Dewar benceno → benceno puede ganar al concertado. La altura concreta de esa cascada en esta app es un modelo, no el número del artículo.",
  },
  {
    id: "prieto2012",
    authors: "Rodríguez-García, A. M.; Cebrián, C.; Prieto, P.; García, J. I.; de la Hoz, A.; Díaz-Ortiz, Á.",
    year: 2012,
    title:
      "DFT Studies on Cobalt-Catalyzed Cyclotrimerization Reactions: The Mechanism and Origin of Reaction Improvement under Microwave Irradiation",
    source: "Chem. Eur. J. 18, 6217–6224",
    doi: "10.1002/chem.201103560",
    usedFor:
      "Orden de magnitud del pozo triplete (~17 kcal/mol bajo el singlete) y del MECP (~10 kcal/mol sobre el singlete) en sistemas sustituidos. Aquí se usa como ancla aproximada, no como energía del acetileno padre.",
  },
  {
    id: "grimme2012",
    authors: "Grimme, S.",
    year: 2012,
    title: "Supramolecular Binding Thermodynamics by Dispersion-Corrected Density Functional Theory",
    source: "Chem. Eur. J. 18, 9955–9964",
    doi: "10.1002/chem.201200497",
    usedFor: "Cuasi-RRHO: interpolación vibración / rotor libre con ν₀ = 100 cm⁻¹.",
  },
  {
    id: "truhlar2011",
    authors: "Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.",
    year: 2011,
    title: "Use of Solution-Phase Vibrational Frequencies in Continuum Models for the Free Energy of Solvation",
    source: "J. Phys. Chem. B 115, 14556–14562",
    doi: "10.1021/jp205508z",
    usedFor: "Corte cuasi-armónico: modos por debajo de 100 cm⁻¹ se suben a 100 cm⁻¹ solo en la entropía.",
  },
  {
    id: "alecu2010",
    authors: "Alecu, I. M.; Zheng, J.; Zhao, Y.; Truhlar, D. G.",
    year: 2010,
    title: "Computational Thermochemistry: Scale Factor Databases and Scale Factors for Vibrational Frequencies Obtained from Electronic Model Chemistries",
    source: "J. Chem. Theory Comput. 6, 2872–2887",
    doi: "10.1021/ct100326h",
    usedFor: "Factores de escala. El de B3LYP/6-31G(d) usado aquí es 0.9614 (Scott y Radom). El resto son factores de trabajo: confirmar en CCCBDB para la base exacta.",
  },
  {
    id: "kozuch2011",
    authors: "Kozuch, S.; Shaik, S.",
    year: 2011,
    title: "How to Conceptualize Catalytic Cycles? The Energetic Span Model",
    source: "Acc. Chem. Res. 44, 101–110",
    doi: "10.1021/ar1000956",
    usedFor: "Span energético, TDTS y TDI. TOF ≈ (kBT/h) exp(−δE/RT).",
  },
  {
    id: "eyring",
    authors: "Eyring, H.; Wigner, E.",
    year: 1935,
    title: "Teoría del estado de transición y corrección de túnel de Wigner",
    source: "J. Chem. Phys. 1935, 3, 107; Z. Phys. Chem. 1932",
    doi: "10.1063/1.1749604",
    usedFor: "k = κ (kBT/h) exp(−ΔG‡/RT), con κ ≈ 1 + (1/24)(hν‡/kBT)².",
  },
];

export const MASSES: Record<string, number> = {
  H: 1.00784,
  C: 12.011,
  N: 14.007,
  O: 15.999,
  F: 18.998403,
  Si: 28.085,
  P: 30.973762,
  Cl: 35.45,
  Co: 58.933194,
  Ni: 58.6934,
  Ru: 101.07,
  Rh: 102.9055,
  Ir: 192.217,
};

export const PH3_MASS = MASSES.P + 3 * MASSES.H;
