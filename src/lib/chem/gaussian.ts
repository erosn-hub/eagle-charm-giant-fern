import { ALKYNES, METHODS, METALS, SOLVENTS, type MethodId } from "./data";
import { stationPieces, type GeomKind } from "./geometry";
import type { JobInput, Station } from "./surface";

const MULT: Record<string, number> = {
  singlete: 1,
  triplete: 3,
  cruce: 1,
  "singlete abierto": 1,
};

function xyzBlock(atoms: { el: string; x: number; y: number; z: number }[]): string {
  return atoms
    .map((a) => `${a.el.padEnd(2)} ${a.x.toFixed(6).padStart(12)} ${a.y.toFixed(6).padStart(12)} ${a.z.toFixed(6).padStart(12)}`)
    .join("\n");
}

function route(method: MethodId, job: "stable" | "opt" | "freq" | "sp" | "irc", kind: GeomKind, solventName: string | null): string {
  const m = METHODS[method];
  const disp = m.dispersionKeyword ? ` ${m.dispersionKeyword}` : "";
  const scrf = solventName ? ` scrf=(smd,solvent=${solventName})` : "";
  const scale = `scale=${m.scale.toFixed(4)}`;
  const baseOpt = m.optBasis;
  const func = m.functional;
  const ts = kind.startsWith("ts");
  if (job === "stable") return `#p stable=opt ${func}/${baseOpt}${disp}${scrf} guess=mix nosymm`;
  if (job === "opt") {
    const opt = ts ? "opt=(ts,calcfc,noeigentest,tight,maxstep=15)" : "opt=(calcfc,tight)";
    return `#p ${opt} ${func}/${baseOpt}${disp}${scrf} guess=read geom=allcheck nosymm int=ultrafine`;
  }
  if (job === "freq") return `#p freq=(noraman,${scale}) ${func}/${baseOpt}${disp}${scrf} guess=read geom=allcheck nosymm int=ultrafine`;
  if (job === "irc")
    return `#p irc=(calcfc,maxpoints=60,recorrect=never,stepsize=8,lqa) ${func}/${baseOpt}${disp}${scrf} guess=read geom=allcheck nosymm`;
  return `#p ${func}/${m.spBasis}${disp}${scrf} guess=read geom=allcheck nosymm int=ultrafine scf=(xqc,maxcycle=256)`;
}

export function gaussianDeck(input: JobInput, stations: Station[], nproc: number, memGb: number): string {
  const solvent = SOLVENTS[input.solvent].smd;
  const blocks: string[] = [];
  blocks.push(header(input, nproc));
  stations.forEach((st, idx) => {
    const pieces = stationPieces({ kind: st.geom, metal: input.metal, alkynes: input.alkynes });
    pieces.forEach((piece, p) => {
      const tag = `${st.id}_${p + 1}`.replace(/[^a-z0-9_]/gi, "");
      const mult = st.spin === "triplete" && piece.role === "complex" ? 3 : MULT[st.spin] ?? 1;
      const charge = 0;
      const title = `${st.label} · fragmento ${piece.label}`.slice(0, 180);
      const jobs: { route: string; read: boolean }[] = [
        { route: route(input.method, "stable", st.geom, solvent), read: false },
        { route: route(input.method, "opt", st.geom, solvent), read: true },
        { route: route(input.method, "freq", st.geom, solvent), read: true },
        { route: route(input.method, "sp", st.geom, solvent), read: true },
      ];
      if (st.kind === "ts" && piece.role === "complex") {
        jobs.push({ route: route(input.method, "irc", st.geom, solvent), read: true });
      }
      jobs.forEach((job, j) => {
        const lines = [
          `%chk=${tag}.chk`,
          `%nprocshared=${nproc}`,
          `%mem=${memGb}GB`,
          job.route,
          "",
          title,
          "",
        ];
        if (!job.read) {
          lines.push(`${charge} ${mult}`, xyzBlock(piece.atoms), "");
        } else {
          lines.push("");
        }
        const body = lines.join("\n");
        blocks.push((idx === 0 && p === 0 && j === 0 ? "" : "--Link1--\n") + body);
      });
    });
  });
  blocks.push(protocolNotes(input));
  return blocks.join("\n");
}

function header(input: JobInput, nproc: number): string {
  const m = METALS[input.metal];
  const method = METHODS[input.method];
  const alks = input.alkynes.map((id) => ALKYNES[id].formula).join(" + ");
  return [
    "! Trímera — baraja Gaussian 16 para la trimerización de alquinos",
    `! Sistema: ${alks}`,
    `! Catalizador modelo: ${m.fragment}`,
    `! Método de trabajo: ${method.label}`,
    `! Disolvente: ${SOLVENTS[input.solvent].name}`,
    "! Las coordenadas son geometrías INICIALES, no puntos estacionarios.",
    "! Orden de cada fragmento: stable=opt → opt → freq → SP de base mayor → IRC si es un TS del complejo.",
    "! No añadas empiricaldispersion sobre wB97XD ni sobre M06/M062X: ya llevan correlación de medio alcance o D2 propio.",
    `! Núcleos pedidos en la ruta: ${nproc}. Ajusta %nprocshared a la cola real.`,
    "! MECP: Gaussian no lo optimiza solo. Usa EasyMECP o sobMECP entre el chk singlete y el triplete del metalaciclo.",
    "! Para la reacción orgánica, un single-point DLPNO-CCSD(T) en ORCA sobre la geometría DFT es el control de oro.",
    "! Revisa <S**2> en el triplete (ideal ~2.0) y la estabilidad de la onda en el singlete del metalaciclo.",
    "",
  ].join("\n");
}

function protocolNotes(input: JobInput): string {
  return [
    "",
    "! --- Después de Gaussian ---",
    "! 1. Un mínimo tiene 0 frecuencias imaginarias; un TS tiene exactamente 1, y el modo debe corresponder al enlace que se forma.",
    "! 2. El IRC tiene que caer en los dos mínimos vecinos. Si no, el TS no pertenece a este paso.",
    "! 3. Termoquímica: el log trae ΔH y ΔG a 1 atm y 298.15 K. Para disolución resta o suma Δn·RT·ln(24.46) según el cambio de molecularidad,",
    "!    y trata los modos < 100 cm⁻¹ con cuasi-RRHO (Grimme, Chem. Eur. J. 2012) — GoodVibes lo hace sobre este log.",
    `! 4. Factor de escala escrito en freq=: ${METHODS[input.method].scale}. Confírmalo en CCCBDB para esta base.`,
    "! 5. counterpoise=2 en el aducto η² si la base de optimización es doble zeta: el tercer alquino está débilmente unido (Hardesty: 12.4 kcal/mol).",
    "! 6. No interpretes el perfil de esta aplicación como el resultado de este input: sustituye las energías de ancla por las del log.",
    "",
  ].join("\n");
}

export function routeCard(input: JobInput): { step: string; detail: string }[] {
  const m = METHODS[input.method];
  const solvent = SOLVENTS[input.solvent].smd;
  return [
    {
      step: "1 · Onda estable",
      detail: `stable=opt guess=mix ${m.functional}/${m.optBasis}. Obligatorio en Co: el metalaciclo singlete suele ser inestable hacia UHF.`,
    },
    {
      step: "2 · Optimización",
      detail: solvent
        ? `opt tight + scrf=(smd,solvent=${solvent}). Geometría en disolvente, no solo la energía.`
        : "opt=(calcfc,tight) en fase gas. Para un TS: opt=(ts,calcfc,noeigentest).",
    },
    {
      step: "3 · Frecuencias",
      detail: `freq scale=${m.scale}. Confirma ν‡ y usa las térmicas solo después del cuasi-armónico.`,
    },
    {
      step: "4 · Single point",
      detail: `${m.functional}/${m.spBasis}${m.dispersionKeyword ? " + " + m.dispersionKeyword : ""}${solvent ? " + SMD" : ""}. Es la energía electrónica que entra en ΔG.`,
    },
    {
      step: "5 · IRC",
      detail: "irc=(calcfc,maxpoints=60,lqa) desde cada TS. Sin IRC no hay estado de transición, solo un punto silla huérfano.",
    },
    {
      step: "6 · Espín",
      detail: input.metal === "CpCo" || input.metal === "CpStarCo"
        ? "Repite el metalaciclo como triplete (0 3). El MECP no sale de un Route estándar."
        : "Comprueba que el singlete del precursor no tiene inestabilidad RHF→UHF.",
    },
  ];
}
