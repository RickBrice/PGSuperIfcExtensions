"""Checks exports that update a copy of the IFC model a project was imported from (/IfcSourceModel, see
devdocs/RoundTripDesign.md).

check_ifc_model: after an IFC model's import, the project is exported as a copy of the model. The copy must keep the
model as it is: the model file isn't changed, every GlobalId of the model is in the copy, and no entity type has fewer
instances. The design-value IDS written with the copy (girders pinned by GlobalId) must pass on it.

check_added_reinforcement: the merge adds PGSuper's strands and bars to girders that have none in the model. The
round-trip export of a PGSuper project is stripped of its strands and bars, and the project imported from it is told
the model had none (its round trip registry, saved in the project, is edited). Its export into the stripped model must
put the strands and bars where an export of the project puts them (ifc_strands.py, ifc_rebar.py).

The copies are written to results/models/ (not tracked): the GlobalIds of what the merge adds differ on every run.
"""

import collections
import hashlib
import re
from pathlib import Path

import ifcopenshell
import ifcopenshell.api

import compare_bridge
import ifc_rebar
import ifc_strands
from check_standard import ifctester_ids

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "models"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def merge_summary(log):
    """what the merge did, from the export log"""
    if not log.exists():
        return "no export log"
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("Updated a copy of "):
            return line.split(": ", 1)[1] if ": " in line[20:] else line
        if "failed" in line:
            return line
    return "the export log has no merge summary"


def ids_failures(ifc, ids):
    """(failed requirements, specifications) of an IDS on a model"""
    spec = ifctester_ids().open(str(ids))
    spec.validate(ifc)
    failed = sum(len(getattr(r, "failures", [])) for s in spec.specifications for r in s.requirements)
    return failed, len(spec.specifications)


def compare_copy(source, copy):
    """problems of the copy: missing GlobalIds and entity types with fewer instances"""
    problems = []
    copy_ids = {e.GlobalId for e in copy.by_type("IfcRoot")}
    missing = [e for e in source.by_type("IfcRoot") if e.GlobalId not in copy_ids]
    if missing:
        problems.append(f"{len(missing)} GlobalIds of the model aren't in the copy (e.g. {missing[0].is_a()} {missing[0].GlobalId})")
    source_counts = collections.Counter(e.is_a() for e in source)
    copy_counts = collections.Counter(e.is_a() for e in copy)
    fewer = [t for t, n in source_counts.items() if copy_counts[t] < n]
    if fewer:
        problems.append("fewer " + ", ".join(f"{t} ({source_counts[t]} -> {copy_counts[t]})" for t in fewer[:5]))
    return problems


def export_copy(run_bridgelink, quote, exe, pgs, copy, ids, source=None):
    """exports the project pgs as a copy of its model (or of source). Returns (ok, message, log)"""
    log = copy.with_suffix(".export.log")
    for f in (copy, ids, log):
        f.unlink(missing_ok=True)
    flag = f"/IfcSourceModel={source}" if source else "/IfcSourceModel"
    ok, message = run_bridgelink(exe, [quote(f"/IfcExport={copy}"), quote(pgs), quote(flag), quote(f"/IfcIds={ids}"), quote(f"/IfcLog={log}")])
    if not ok or not copy.exists():
        return False, f"the export failed ({message}; {merge_summary(log)})", log
    return True, merge_summary(log), log


def check_ifc_model(run_bridgelink, quote, exe, run, model_ifc, pgs):
    """the round trip of an imported IFC model. Returns a message for the summary"""
    OUT.mkdir(parents=True, exist_ok=True)
    copy = OUT / f"{run}.roundtrip.ifc"
    ids = OUT / f"{run}.roundtrip.ids"
    before = sha256(model_ifc)
    ok, summary, _ = export_copy(run_bridgelink, quote, exe, pgs, copy, ids)
    if not ok:
        return f", round trip: {summary}"

    problems = []
    if sha256(model_ifc) != before:
        problems.append("the model file changed")
    copy_model = ifcopenshell.open(str(copy))
    problems += compare_copy(ifcopenshell.open(str(model_ifc)), copy_model)
    failed, total = ids_failures(copy_model, ids)
    result = f"design-value IDS: {failed} failures in {total} specifications"
    if problems:
        return f", round trip: {summary}; PROBLEMS: {'; '.join(problems)}; {result}"
    return f", round trip: {summary}; the model is kept (GlobalIds, entities); {result}"


def strip_reinforcement(source, stripped):
    """the model without its strands and bars (and their assemblies). Returns the number removed"""
    model = ifcopenshell.open(str(source))
    removed = 0
    for assembly in list(model.by_type("IfcElementAssembly")):
        parts = [p for rel in assembly.IsDecomposedBy for p in rel.RelatedObjects]
        if parts and all(p.is_a("IfcTendon") or p.is_a("IfcReinforcingBar") for p in parts):
            for p in parts:
                ifcopenshell.api.run("root.remove_product", model, product=p)
                removed += 1
            ifcopenshell.api.run("root.remove_product", model, product=assembly)
    model.write(str(stripped))
    return removed


def forget_reinforcement(pgs, edited):
    """a copy of the project whose round trip registry says its model has no strands or bars"""
    text = Path(pgs).read_text(encoding="utf-8")
    elements = re.findall(r'<Element version="1">.*?</Element>', text, re.S)
    if not elements:
        return False
    kept = [e for e in elements if not re.search(r"<Key>(tendons|bars):", e)]
    start = text.index('<Element version="1">')
    end = text.rindex("</Element>") + len("</Element>")
    text = text[:start] + "".join(kept) + text[end:]
    text = re.sub(r"<ElementCount>\d+</ElementCount>", f"<ElementCount>{len(kept)}</ElementCount>", text)
    Path(edited).write_text(text, encoding="utf-8")
    return True


def reinforcement(ifc):
    summary = ifc_strands.extract(str(ifc))
    summary["values"].update(ifc_rebar.extract(str(ifc))["values"])
    return summary


def check_added_reinforcement(run_bridgelink, quote, exe, run, export_ifc, pgs, reexport):
    """the strands and bars the merge adds, compared with an export of the project (reexport). Returns a message"""
    OUT.mkdir(parents=True, exist_ok=True)
    stripped = OUT / f"{run}.stripped.ifc"
    edited = OUT / f"{run}.stripped.pgs"
    copy = OUT / f"{run}.added.ifc"
    ids = OUT / f"{run}.added.ids"
    nRemoved = strip_reinforcement(export_ifc, stripped)
    if not forget_reinforcement(pgs, edited):
        return ", added reinforcement not checked: the project has no round trip registry"
    ok, summary, _ = export_copy(run_bridgelink, quote, exe, edited, copy, ids, source=stripped)
    if not ok:
        return f", added reinforcement: {summary}"

    comparison = compare_bridge.compare(reinforcement(copy), reinforcement(reexport))
    counts, compared = compare_bridge.score(comparison)
    mismatches = compared - counts["match"]
    return (f", added reinforcement: {nRemoved} strands and bars removed from the model, then {summary}; "
            f"{counts['match']} of {compared} strand and bar values as in the project's export" + ("" if mismatches == 0 else f" ({mismatches} differ)"))
