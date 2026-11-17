import { ALKYNES, MASSES, METALS, type AlkyneId, type MetalId, type SubstituentId } from "./data";

export type XYZ = { el: string; x: number; y: number; z: number };
export type GeomKind =
  | "pre"
  | "bis"
  | "ts-oc"
  | "mc"
  | "adduct"
  | "ts-collapse"
  | "arene"
  | "ts-release"
  | "prod"
  | "reactants"
  | "thermal-vdw"
  | "thermal-ts"
  | "thermal-birad"
  | "thermal-prod";

export type Piece = { label: string; atoms: XYZ[]; sigma: number; role: "alkyne" | "ligand" | "arene" | "complex" };

type V = [number, number, number];
const add = (a: V, b: V): V => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a: V, b: V): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const scl = (a: V, s: number): V => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a: V, b: V) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: V, b: V): V => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = (a: V) => Math.hypot(a[0], a[1], a[2]);
const norm = (a: V): V => {
  const L = len(a) || 1;
  return scl(a, 1 / L);
};

function frame(u: V): { v: V; w: V } {
  const t: V = Math.abs(u[1]) < 0.85 ? [0, 1, 0] : [1, 0, 0];
  const v = norm(cross(u, t));
  return { v, w: norm(cross(u, v)) };
}

function push(atoms: XYZ[], el: string, p: V) {
  atoms.push({ el, x: round(p[0]), y: round(p[1]), z: round(p[2]) });
}
const round = (x: number) => Math.round(x * 1e6) / 1e6;

function addSubstituent(atoms: XYZ[], attach: V, outward: V, id: SubstituentId) {
  const u = norm(outward);
  if (id === "H") {
    push(atoms, "H", add(attach, scl(u, 1.08)));
    return;
  }
  const { v, w } = frame(u);
  const cap = (origin: V, axis: V, bond: number, el: string) => {
    for (let k = 0; k < 3; k++) {
      const ang = (k * 2 * Math.PI) / 3;
      const perp = add(scl(v, Math.cos(ang)), scl(w, Math.sin(ang)));
      const dir = norm(add(scl(axis, 1 / 3), scl(perp, Math.sqrt(8 / 9))));
      push(atoms, el, add(origin, scl(dir, bond)));
    }
  };
  if (id === "Me") {
    const c = add(attach, scl(u, 1.51));
    push(atoms, "C", c);
    cap(c, u, 1.09, "H");
    return;
  }
  if (id === "tBu") {
    const c = add(attach, scl(u, 1.54));
    push(atoms, "C", c);
    for (let k = 0; k < 3; k++) {
      const ang = (k * 2 * Math.PI) / 3;
      const perp = add(scl(v, Math.cos(ang)), scl(w, Math.sin(ang)));
      const dir = norm(add(scl(u, 1 / 3), scl(perp, Math.sqrt(8 / 9))));
      const mc = add(c, scl(dir, 1.53));
      push(atoms, "C", mc);
      cap(mc, dir, 1.09, "H");
    }
    return;
  }
  if (id === "CF3") {
    const c = add(attach, scl(u, 1.5));
    push(atoms, "C", c);
    cap(c, u, 1.34, "F");
    return;
  }
  if (id === "CN") {
    const c = add(attach, scl(u, 1.38));
    push(atoms, "C", c);
    push(atoms, "N", add(c, scl(u, 1.16)));
    return;
  }
  if (id === "OMe") {
    const o = add(attach, scl(u, 1.36));
    push(atoms, "O", o);
    const side = norm(add(scl(u, 0.3), scl(v, 0.95)));
    const c = add(o, scl(side, 1.43));
    push(atoms, "C", c);
    cap(c, side, 1.09, "H");
    return;
  }
  if (id === "CO2Me") {
    const c = add(attach, scl(u, 1.48));
    push(atoms, "C", c);
    push(atoms, "O", add(c, norm(add(scl(u, 0.35), scl(v, 1)))));
    const o2 = add(c, scl(norm(add(scl(u, 0.4), scl(v, -1))), 1.34));
    push(atoms, "O", o2);
    const mc = add(o2, scl(norm(add(scl(u, 1), scl(v, -0.45))), 1.42));
    push(atoms, "C", mc);
    cap(mc, norm(sub(mc, o2)), 1.09, "H");
    return;
  }
  if (id === "TMS") {
    const si = add(attach, scl(u, 1.88));
    push(atoms, "Si", si);
    for (let k = 0; k < 3; k++) {
      const ang = (k * 2 * Math.PI) / 3;
      const perp = add(scl(v, Math.cos(ang)), scl(w, Math.sin(ang)));
      const dir = norm(add(scl(u, 1 / 3), scl(perp, Math.sqrt(8 / 9))));
      const c = add(si, scl(dir, 1.88));
      push(atoms, "C", c);
      cap(c, dir, 1.09, "H");
    }
    return;
  }
  const center = add(attach, scl(u, 2.8));
  const ipso = add(center, scl(u, -1.4));
  push(atoms, "C", ipso);
  for (let k = 1; k < 6; k++) {
    const ang = Math.PI + (k * Math.PI) / 3;
    const p = add(center, add(scl(u, 1.4 * Math.cos(ang)), scl(v, 1.4 * Math.sin(ang))));
    push(atoms, "C", p);
    push(atoms, "H", add(p, scl(norm(sub(p, center)), 1.08)));
  }
}

export function alkyneMolecule(id: AlkyneId): XYZ[] {
  const d = ALKYNES[id];
  const atoms: XYZ[] = [];
  const c1: V = [0, 0, -0.62];
  const c2: V = [0, 0, 0.62];
  push(atoms, "C", c1);
  push(atoms, "C", c2);
  addSubstituent(atoms, c1, [0, 0, -1], d.a);
  addSubstituent(atoms, c2, [0, 0, 1], d.b);
  return atoms;
}

function ph3(): XYZ[] {
  const atoms: XYZ[] = [];
  push(atoms, "P", [0, 0, 0]);
  const u: V = [0, 0, 1];
  const { v, w } = frame(u);
  for (let k = 0; k < 3; k++) {
    const ang = (k * 2 * Math.PI) / 3;
    const perp = add(scl(v, Math.cos(ang)), scl(w, Math.sin(ang)));
    push(atoms, "H", add([0, 0, 0], scl(norm(add(scl(u, 0.12), scl(perp, 0.99))), 1.42)));
  }
  return atoms;
}

function cp(star: boolean): XYZ[] {
  const atoms: XYZ[] = [];
  const center: V = [0, 0, -1.76];
  const { v, w } = frame([0, 0, -1]);
  for (let k = 0; k < 5; k++) {
    const ang = (k * 2 * Math.PI) / 5;
    const radial = add(scl(v, Math.cos(ang)), scl(w, Math.sin(ang)));
    const c = add(center, scl(radial, 1.22));
    push(atoms, "C", c);
    if (star) {
      const me = add(c, scl(radial, 1.51));
      push(atoms, "C", me);
      const { v: mv, w: mw } = frame(radial);
      for (let h = 0; h < 3; h++) {
        const a = (h * 2 * Math.PI) / 3;
        const perp = add(scl(mv, Math.cos(a)), scl(mw, Math.sin(a)));
        const dir = norm(add(scl(radial, 1 / 3), scl(perp, Math.sqrt(8 / 9))));
        push(atoms, "H", add(me, scl(dir, 1.09)));
      }
    } else {
      push(atoms, "H", add(c, scl(radial, 1.08)));
    }
  }
  return atoms;
}

function metalAtoms(metal: MetalId): XYZ[] {
  if (metal === "none") return [];
  const atoms: XYZ[] = [];
  push(atoms, METALS[metal].element, [0, 0, 0]);
  if (metal !== "Ni0") atoms.push(...cp(metal === "CpStarCo"));
  if (metal === "CpRuCl") push(atoms, "Cl", [2.2, 0.4, 0.2]);
  return atoms;
}

function boundAlkyne(id: AlkyneId, center: V, along: V): XYZ[] {
  const d = ALKYNES[id];
  const u = norm(along);
  const c1 = add(center, scl(u, -0.63));
  const c2 = add(center, scl(u, 0.63));
  const atoms: XYZ[] = [];
  push(atoms, "C", c1);
  push(atoms, "C", c2);
  addSubstituent(atoms, c1, scl(u, -1), d.a);
  addSubstituent(atoms, c2, u, d.b);
  return atoms;
}

function twoPhosphines(close: boolean): XYZ[] {
  if (close) {
    return [...shift(ph3(), [1.7, 0.35, 1.15]), ...shift(ph3(), [-1.7, -0.35, 1.15])];
  }
  return [];
}

function shift(atoms: XYZ[], d: V): XYZ[] {
  return atoms.map((a) => ({ el: a.el, x: round(a.x + d[0]), y: round(a.y + d[1]), z: round(a.z + d[2]) }));
}

export function benzeneMolecule(ids: [AlkyneId, AlkyneId, AlkyneId]): XYZ[] {
  const atoms: XYZ[] = [];
  const ends: SubstituentId[] = ids.flatMap((id) => [ALKYNES[id].a, ALKYNES[id].b]);
  for (let k = 0; k < 6; k++) {
    const ang = (k * Math.PI) / 3 - Math.PI / 6;
    const radial: V = [Math.cos(ang), Math.sin(ang), 0];
    const c = scl(radial, 1.397);
    push(atoms, "C", c);
    addSubstituent(atoms, c, radial, ends[k]);
  }
  return atoms;
}

function precursor(metal: MetalId): XYZ[] {
  return [...metalAtoms(metal), ...twoPhosphines(true)];
}

function bis(metal: MetalId, a: AlkyneId, b: AlkyneId, approach: number): XYZ[] {
  const x = 1.45 - approach * 0.25;
  const y = 0.72 - approach * 0.35;
  return [
    ...metalAtoms(metal),
    ...boundAlkyne(a, [x, 0, 0.9], [0, y > 0.2 ? 1 : 0.4, 0]),
    ...boundAlkyne(b, [-x, 0, 0.9], [0, y > 0.2 ? 1 : 0.4, 0]),
  ];
}

function cycleCarbons(): V[] {
  return [
    [1.32, -0.9, 0.7],
    [1.05, -0.1, 1.82],
    [-0.1, 1.05, 1.82],
    [-0.9, 1.32, 0.7],
  ];
}

function metallacycle(metal: MetalId, a: AlkyneId, b: AlkyneId, third: AlkyneId | null, thirdZ: number): XYZ[] {
  const atoms = [...metalAtoms(metal)];
  const [a1, b1, b2, a2] = cycleCarbons();
  for (const p of [a1, b1, b2, a2]) push(atoms, "C", p);
  const A = ALKYNES[a];
  const B = ALKYNES[b];
  const out = (p: V) => norm([p[0], p[1], 0.3]);
  addSubstituent(atoms, a1, out(a1), A.a);
  addSubstituent(atoms, b1, out(b1), A.b);
  addSubstituent(atoms, b2, out(b2), B.a);
  addSubstituent(atoms, a2, out(a2), B.b);
  if (third) atoms.push(...boundAlkyne(third, [0, 0, thirdZ], [1, 0, 0]));
  return atoms;
}

function areneOnMetal(metal: MetalId, ids: [AlkyneId, AlkyneId, AlkyneId], z: number): XYZ[] {
  return [...metalAtoms(metal), ...shift(benzeneMolecule(ids), [0.3, 0, z])];
}

function threeFar(ids: [AlkyneId, AlkyneId, AlkyneId], R: number): XYZ[] {
  const atoms: XYZ[] = [];
  ids.forEach((id, i) => {
    const th = (i * 2 * Math.PI) / 3;
    const mid: V = [R * Math.cos(th), R * Math.sin(th), 0];
    const tan: V = [-Math.sin(th), Math.cos(th), 0];
    const d = ALKYNES[id];
    const c1 = add(mid, scl(tan, -0.62));
    const c2 = add(mid, scl(tan, 0.62));
    push(atoms, "C", c1);
    push(atoms, "C", c2);
    addSubstituent(atoms, c1, scl(tan, -1), d.a);
    addSubstituent(atoms, c2, tan, d.b);
  });
  return atoms;
}

export type GeomRequest = { kind: GeomKind; metal: MetalId; alkynes: [AlkyneId, AlkyneId, AlkyneId] };

function alkyneSigma(id: AlkyneId): number {
  const d = ALKYNES[id];
  if (d.a === "H" && d.b === "H") return 2;
  if (d.a === "Me" && d.b === "Me") return 6;
  if (d.a === d.b) return 2;
  if (d.a === "H" || d.b === "H") {
    const other = d.a === "H" ? d.b : d.a;
    if (other === "Me" || other === "CF3" || other === "tBu") return 3;
  }
  return 1;
}

function benzeneSigma(ids: [AlkyneId, AlkyneId, AlkyneId]): number {
  if (ids.every((id) => id === "acetylene")) return 12;
  if (ids.every((id) => id === "butyne")) return 6;
  return 1;
}

export function stationPieces(req: GeomRequest): Piece[] {
  const [ia, ib, ic] = req.alkynes;
  const free = (id: AlkyneId, label: string): Piece => ({
    label,
    atoms: alkyneMolecule(id),
    sigma: alkyneSigma(id),
    role: "alkyne",
  });
  const ligand = (n: number): Piece => ({ label: `PH₃ (${n})`, atoms: ph3(), sigma: 3, role: "ligand" });
  const benzene = (): Piece => ({
    label: "areno",
    atoms: benzeneMolecule(req.alkynes),
    sigma: benzeneSigma(req.alkynes),
    role: "arene",
  });
  const complex = (label: string, atoms: XYZ[]): Piece => ({ label, atoms, sigma: 1, role: "complex" });

  if (req.metal === "none") {
    if (req.kind === "thermal-prod" || req.kind === "prod") return [benzene()];
    if (req.kind === "thermal-vdw") return [complex("precomplejo", threeFar(req.alkynes, 3.2))];
    if (req.kind === "thermal-ts") return [complex("TS concertado", threeFar(req.alkynes, 1.72))];
    if (req.kind === "thermal-birad") return [complex("estructura birradical de partida", threeFar(req.alkynes, 2.15))];
    return [free(ia, "alquino A"), free(ib, "alquino B"), free(ic, "alquino C")];
  }

  if (req.kind === "pre") return [complex("precursor", precursor(req.metal)), free(ia, "A"), free(ib, "B"), free(ic, "C")];
  if (req.kind === "prod") return [complex("precursor", precursor(req.metal)), benzene()];
  if (req.kind === "bis") return [complex("bis-η²", bis(req.metal, ia, ib, 0)), ligand(1), ligand(2), free(ic, "C")];
  if (req.kind === "ts-oc") return [complex("TS acoplamiento", bis(req.metal, ia, ib, 0.85)), ligand(1), ligand(2), free(ic, "C")];
  if (req.kind === "mc") return [complex("metalaciclo", metallacycle(req.metal, ia, ib, null, 0)), ligand(1), ligand(2), free(ic, "C")];
  if (req.kind === "adduct") return [complex("aducto", metallacycle(req.metal, ia, ib, ic, 4.05)), ligand(1), ligand(2)];
  if (req.kind === "ts-collapse") return [complex("TS [4+2]", metallacycle(req.metal, ia, ib, ic, 3.55)), ligand(1), ligand(2)];
  if (req.kind === "arene") return [complex("η⁴-areno", areneOnMetal(req.metal, req.alkynes, 1.75)), ligand(1), ligand(2)];
  if (req.kind === "ts-release") return [complex("TS de liberación", areneOnMetal(req.metal, req.alkynes, 2.45)), ligand(1), ligand(2)];
  return [complex("especie", precursor(req.metal))];
}

export function buildGeometry(req: GeomRequest): XYZ[] {
  const pieces = stationPieces(req);
  const all: XYZ[] = [];
  pieces.forEach((p, i) => {
    all.push(...shift(p.atoms, [i * 14, 0, 0]));
  });
  return all;
}

export function minDistance(atoms: XYZ[]): { d: number; pair: string } {
  let d = Infinity;
  let pair = "";
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const r = Math.hypot(atoms[i].x - atoms[j].x, atoms[i].y - atoms[j].y, atoms[i].z - atoms[j].z);
      if (r < d) {
        d = r;
        pair = `${atoms[i].el}${i + 1}–${atoms[j].el}${j + 1}`;
      }
    }
  }
  return { d, pair };
}

export function formingCC(atoms: XYZ[]): number | null {
  const carbons = atoms.filter((a) => a.el === "C");
  let forming = Infinity;
  for (let i = 0; i < carbons.length; i++) {
    for (let j = i + 1; j < carbons.length; j++) {
      const r = Math.hypot(carbons[i].x - carbons[j].x, carbons[i].y - carbons[j].y, carbons[i].z - carbons[j].z);
      if (r > 1.45 && r < forming) forming = r;
    }
  }
  return Number.isFinite(forming) ? forming : null;
}

export function inertia(atoms: XYZ[], sigma: number): { linear: boolean; moments: [number, number, number]; sigma: number } {
  let mx = 0,
    my = 0,
    mz = 0,
    mass = 0;
  const pts = atoms.map((a) => {
    const m = MASSES[a.el] ?? 12;
    mass += m;
    mx += m * a.x;
    my += m * a.y;
    mz += m * a.z;
    return { m, x: a.x, y: a.y, z: a.z };
  });
  if (mass < 1e-6) return { linear: false, moments: [1e-46, 1e-46, 1e-46], sigma };
  mx /= mass;
  my /= mass;
  mz /= mass;
  let ixx = 0,
    iyy = 0,
    izz = 0,
    ixy = 0,
    ixz = 0,
    iyz = 0;
  for (const a of pts) {
    const x = a.x - mx,
      y = a.y - my,
      z = a.z - mz;
    ixx += a.m * (y * y + z * z);
    iyy += a.m * (x * x + z * z);
    izz += a.m * (x * x + y * y);
    ixy -= a.m * x * y;
    ixz -= a.m * x * z;
    iyz -= a.m * y * z;
  }
  const ev = jacobi(ixx, iyy, izz, ixy, ixz, iyz).sort((p, q) => p - q);
  const si = (v: number) => Math.max(v, 0) * 1.6605390666e-47;
  const moments: [number, number, number] = [si(ev[0]), si(ev[1]), si(ev[2])];
  const linear = moments[2] > 0 && moments[0] < moments[2] * 0.015;
  return { linear, moments, sigma };
}

function jacobi(ixx: number, iyy: number, izz: number, ixy: number, ixz: number, iyz: number) {
  const a = [
    [ixx, ixy, ixz],
    [ixy, iyy, iyz],
    [ixz, iyz, izz],
  ];
  for (let iter = 0; iter < 20; iter++) {
    let p = 0,
      q = 1,
      max = Math.abs(a[0][1]);
    if (Math.abs(a[0][2]) > max) {
      max = Math.abs(a[0][2]);
      q = 2;
    }
    if (Math.abs(a[1][2]) > max) {
      p = 1;
      q = 2;
    }
    if (max < 1e-10) break;
    const app = a[p][p],
      aqq = a[q][q],
      apq = a[p][q];
    const tau = (aqq - app) / (2 * apq || 1e-30);
    const t = Math.sign(tau || 1) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
    const c = 1 / Math.sqrt(1 + t * t);
    const s = t * c;
    a[p][p] = app - t * apq;
    a[q][q] = aqq + t * apq;
    a[p][q] = a[q][p] = 0;
    for (let r = 0; r < 3; r++) {
      if (r === p || r === q) continue;
      const arp = a[r][p],
        arq = a[r][q];
      a[r][p] = a[p][r] = c * arp - s * arq;
      a[r][q] = a[q][r] = s * arp + c * arq;
    }
  }
  return [a[0][0], a[1][1], a[2][2]];
}

void dot;
