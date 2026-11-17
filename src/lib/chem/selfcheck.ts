import { benzeneMolecule, alkyneMolecule, minDistance, stationPieces } from "./geometry";
import { compute, type JobInput } from "./surface";
import { moleculeThermo } from "./thermo";

const base: JobInput = {
  alkynes: ["acetylene", "acetylene", "acetylene"],
  metal: "CpCo",
  method: "b3lyp",
  solvent: "gas",
  temperature: 298.15,
  thermo: "rrho",
  standardState: "atm",
  twoState: false,
  thermalMode: "concerted",
  wigner: true,
};

const profile = compute(base);
console.log("--- CpCo acetileno B3LYP gas ---");
for (const s of profile.stations) {
  console.log(
    s.id.padEnd(6),
    "H",
    s.dH.toFixed(2).padStart(8),
    "G",
    s.dG.toFixed(2).padStart(8),
    "dS",
    s.dS.toFixed(2).padStart(7),
    "n",
    s.nMolecules,
  );
}
const acet = moleculeThermo(alkyneMolecule("acetylene"), 2, [612, 612, 730, 730, 1974, 3289, 3374], 298.15, "rrho", 1);
const benz = moleculeThermo(
  benzeneMolecule(["acetylene", "acetylene", "acetylene"]),
  12,
  [410, 606, 673, 849, 993, 1010, 1148, 1178, 1350, 1485],
  298.15,
  "rrho",
  1,
);
console.log("S acetileno", acet.S.toFixed(2), "cal (exp ~48.0)");
console.log("S benceno", benz.S.toFixed(2), "cal (exp ~64.3)");
console.log("masa acet", acet.mass.toFixed(2), "natoms", acet.nAtoms);

let clashes = 0;
for (const metal of ["CpCo", "CpStarCo", "Ni0", "CpRuCl"] as const) {
  for (const s of compute({ ...base, metal }).stations) {
    for (const pc of stationPieces({ kind: s.geom, metal, alkynes: base.alkynes })) {
      const md = minDistance(pc.atoms);
      if (pc.atoms.length < 3 || md.d < 0.72) {
        clashes++;
        console.log("PROBLEM", metal, s.id, pc.label, "d", md.d.toFixed(3), md.pair, "N", pc.atoms.length);
      }
    }
  }
}
console.log("clashes", clashes);

const propyne = compute({ ...base, alkynes: ["propyne", "propyne", "propyne"], method: "composite", solvent: "toluene", standardState: "molar", twoState: true });
console.log("--- propino regio ---");
for (const iso of propyne.isomers) console.log(iso.name, (iso.share * 100).toFixed(1) + "%", "ddG", iso.dG.toFixed(2), iso.pattern);

const rh = compute({ ...base, metal: "CpRh" });
const coTs = profile.stations.find((s) => s.id === "ts1")!;
const coBis = profile.stations.find((s) => s.id === "bis")!;
const rhTs = rh.stations.find((s) => s.id === "ts1")!;
const rhBis = rh.stations.find((s) => s.id === "bis")!;
console.log("barrera H Co", (coTs.dH - coBis.dH).toFixed(2), "Rh", (rhTs.dH - rhBis.dH).toFixed(2));

const dmad = compute({ ...base, alkynes: ["dmad", "dmad", "dmad"] });
const dmadTs = dmad.stations.find((s) => s.id === "ts1")!;
const dmadBis = dmad.stations.find((s) => s.id === "bis")!;
console.log("barrera H DMAD", (dmadTs.dH - dmadBis.dH).toFixed(2), "vs acet", (coTs.dH - coBis.dH).toFixed(2));

const thermal = compute({ ...base, metal: "none", temperature: 700 });
console.log("--- termica 700 K ---");
for (const s of thermal.stations) console.log(s.id, s.dH.toFixed(1), s.dG.toFixed(1));
for (const step of thermal.steps) console.log(step.name, "dG", step.dG.toFixed(1), "k", step.k.toExponential(2));
if (profile.span) console.log("span", profile.span.deltaE.toFixed(2), profile.span.tdts, "/", profile.span.tdi, "TOF", profile.span.tof.toExponential(2));
