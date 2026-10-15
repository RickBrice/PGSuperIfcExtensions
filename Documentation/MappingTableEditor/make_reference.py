"""Writes src/reference.md, the list of element roles and targets, from IfcTargets.cpp.

The targets and element roles are defined in code; this keeps the documentation in step with them.
Run by BuildDocumentation.bat before Doxygen.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "IfcTargets.cpp"
OUTPUT = HERE / "src" / "reference.md"

KIND_TEXT = {
    "Stress": "stress",
    "Length": "length",
    "Angle": "angle",
    "Ratio": "ratio",
    "Count": "count",
    "Boolean": "true/false",
    "Text": "text",
    "TextList": "list of text",
    "Force": "force",
}


def main():
    code = SOURCE.read_text(encoding="utf-8", errors="replace")

    roles = re.findall(r'\{\s*ElementKind::(\w+),\s*"([a-z_]+)"\s*\}', code)
    role_of_kind = {kind: name for kind, name in roles}

    targets = []
    for m in re.finditer(r'\{\s*\.name\s*=\s*"([^"]+)"\s*,\s*\.element\s*=\s*ElementKind::(\w+)\s*,\s*\.kind\s*=\s*ValueKind::(\w+)(.*?)\}', code, re.S):
        name, element, kind, rest = m.groups()
        description = re.search(r'\.description\s*=\s*"([^"]*)"', rest)
        targets.append((name, role_of_kind.get(element, element), KIND_TEXT.get(kind, kind), description.group(1) if description else ""))

    if not roles or not targets:
        raise SystemExit(f"No element roles or targets found in {SOURCE}")

    lines = [
        "Element Roles and Targets {#reference}",
        "==========================",
        "",
        "This list is made from the IFC extension's source code when the documentation is built, so it matches this version of the extension.",
        "",
        "Element Roles",
        "-------------",
        "",
        "An element role is a kind of IFC element that PGSuper reads or writes. A table's selectors say which IFC elements play each role (see @ref element_roles).",
        "",
        "| Element role |",
        "|---|",
    ]
    lines += [f"| `{name}` |" for _, name in roles]
    lines += [
        "",
        "Targets",
        "-------",
        "",
        "A target is a PGSuper value with a fixed name. A table can't add targets; it says where each target is found in IFC (see @ref targets) and which exported properties hold it (see @ref property_sets).",
        "",
        "| Target | Element role | Kind | Description |",
        "|---|---|---|---|",
    ]
    lines += [f"| `{name}` | `{role}` | {kind} | {description} |" for name, role, kind, description in targets]
    lines.append("")

    OUTPUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"{OUTPUT}: {len(roles)} element roles, {len(targets)} targets")


if __name__ == "__main__":
    main()
