"""Compares the table-driven content of IFC exports: property sets, quantity sets, classifications,
and element names.

An IFC file is reduced to a canonical text listing, one line per property, grouped by element:

   IfcBeam "Span 1, Girder A" #0 | Pset_BeamCommon | Span = IFCPOSITIVELENGTHMEASURE(100.5) [FOOT]

Elements are identified by entity type, Name, and their order among elements with the same type and
Name, because GlobalIds differ from one export to the next. Property sets of an element are sorted by
name, so the order in which relationships are written doesn't matter. Properties keep their order.
Numbers are written with 10 significant digits.

Used to check that the table-driven exporter (devdocs/MappingTablesDesign.md, M3) writes the same
content as the exporter it replaces.

Usage:
   python compare_export.py dump model.ifc [listing.txt[.gz]]      write the listing (default: stdout)
   python compare_export.py compare old.ifc|old.txt[.gz] new.ifc   show the differences, exit code 1 if any
"""
import difflib
import gzip
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------------------------
# STEP (ISO 10303-21) parsing, enough for the DATA section of an IFC file

TOKEN = re.compile(r"""
    (?P<string>'(?:[^']|'')*')
  | (?P<ref>\#\d+)
  | (?P<enum>\.[A-Z0-9_]+\.)
  | (?P<number>[+-]?(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?))
  | (?P<keyword>[A-Z][A-Z0-9_]*)
  | (?P<symbol>[(),$*])
  | (?P<space>\s+)
""", re.VERBOSE)


class Ref(int):
    pass


class Enum(str):
    pass


class Typed:
    """a typed parameter, e.g. IFCLABEL('x')"""
    def __init__(self, name, args):
        self.name = name
        self.args = args


def tokenize(text):
    pos = 0
    while pos < len(text):
        m = TOKEN.match(text, pos)
        if not m:
            raise ValueError(f"can't parse STEP near: {text[pos:pos + 40]!r}")
        pos = m.end()
        kind = m.lastgroup
        if kind != "space":
            yield kind, m.group(kind)


def parse_args(tokens):
    """parses '(' ... ')' - the opening parenthesis has been read"""
    args = []
    for kind, value in tokens:
        if kind == "symbol":
            if value == ")":
                return args
            if value == ",":
                continue
            if value == "(":
                args.append(parse_args(tokens))
            elif value in "$*":
                args.append(None)
        elif kind == "string":
            args.append(value[1:-1].replace("''", "'"))
        elif kind == "ref":
            args.append(Ref(int(value[1:])))
        elif kind == "enum":
            args.append(Enum(value[1:-1]))
        elif kind == "number":
            args.append(float(value) if any(c in value for c in ".eE") else int(value))
        elif kind == "keyword":
            next(tokens)  # "("
            args.append(Typed(value, parse_args(tokens)))
    raise ValueError("unexpected end of arguments")


def read_step(path):
    """{id: (TYPE, args)} for the DATA section"""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    data = text[text.index("DATA;") + 5:text.rindex("ENDSEC;")]
    entities = {}
    # an instance ends with ");" outside of strings. Strings can contain ";" so split carefully
    for m in re.finditer(r"#(\d+)\s*=\s*([A-Z0-9_]+)\s*\(((?:[^;']|'(?:[^']|'')*')*)\)\s*;", data):
        tokens = tokenize(m.group(3) + ")")
        entities[int(m.group(1))] = (m.group(2), parse_args(tokens))
    return entities


# ---------------------------------------------------------------------------------------------
# canonical listing

def fmt_number(value):
    if isinstance(value, float):
        text = f"{value:.10g}"
        return text if any(c in text for c in ".e") else text + "."
    return str(value)


class Listing:
    def __init__(self, entities):
        self.e = entities

    def get(self, ref):
        return self.e[ref] if isinstance(ref, Ref) and ref in self.e else (None, [])

    def unit(self, ref):
        if ref is None:
            return ""
        kind, args = self.get(ref)
        if kind == "IFCCONVERSIONBASEDUNIT":
            return f" [{args[2]}]"
        if kind == "IFCSIUNIT":
            return f" [{(args[2] or '')}{args[3]}]"
        if kind == "IFCDERIVEDUNIT":
            return f" [{args[2] or args[1]}]"
        return f" [{kind}]" if kind else ""

    def value(self, v):
        if v is None:
            return "$"
        if isinstance(v, Typed):
            return f"{v.name}({', '.join(self.value(a) for a in v.args)})"
        if isinstance(v, list):
            return "(" + ", ".join(self.value(a) for a in v) + ")"
        if isinstance(v, Ref):
            kind, args = self.get(v)
            return kind or "#?"
        if isinstance(v, Enum):
            return f".{v}."
        if isinstance(v, str):
            return repr(v)
        return fmt_number(v)

    def spec(self, text):
        return f" {{{text}}}" if text else ""

    def prop(self, ref):
        kind, a = self.get(ref)
        if kind == "IFCPROPERTYSINGLEVALUE":
            return f"{a[0]} = {self.value(a[2])}{self.unit(a[3])}{self.spec(a[1])}"
        if kind == "IFCPROPERTYENUMERATEDVALUE":
            enum = ""
            if a[3] is not None:
                ek, ea = self.get(a[3])
                enum = f" of {ea[0]}{self.value(ea[1])}{self.unit(ea[2])}"
            return f"{a[0]} = enum {self.value(a[2])}{enum}{self.spec(a[1])}"
        if kind == "IFCPROPERTYLISTVALUE":
            return f"{a[0]} = list {self.value(a[2])}{self.unit(a[3])}{self.spec(a[1])}"
        if kind and kind.startswith("IFCQUANTITY"):
            # IfcPhysicalSimpleQuantity: Name, Description, Unit, value, Formula
            return f"{a[0]} = {kind}({fmt_number(a[3]) if a[3] is not None else '$'}){self.unit(a[2])}{self.spec(a[1])}"
        if kind == "IFCPHYSICALCOMPLEXQUANTITY":
            return f"{a[0]} = complex {self.value(a[2])}"
        return f"{kind}"

    def pset_lines(self, ref):
        """(title, [property lines]) of a property set, quantity set, or material properties"""
        kind, a = self.get(ref)
        if kind == "IFCPROPERTYSET":  # GlobalId, OwnerHistory, Name, Description, HasProperties
            return f"{a[2]}{self.spec(a[3])}", [self.prop(p) for p in a[4]]
        if kind == "IFCELEMENTQUANTITY":  # ..., Name, Description, MethodOfMeasurement, Quantities
            return f"{a[2]}{self.spec(a[3])} ({a[4]})", [self.prop(q) for q in a[5]]
        if kind == "IFCMATERIALPROPERTIES":  # Name, Description, Properties, Material
            return f"{a[0]}{self.spec(a[1])}", [self.prop(p) for p in a[2]]
        return kind or "?", []

    def dump(self):
        e = self.e
        objects = {}  # ref -> {"psets": [...], "classes": [...]}

        def slot(ref):
            return objects.setdefault(ref, {"psets": [], "classes": []})

        for ref, (kind, a) in e.items():
            if kind == "IFCRELDEFINESBYPROPERTIES":  # ..., RelatedObjects, RelatingPropertyDefinition
                for obj in a[4]:
                    slot(obj)["psets"].append(a[5])
            elif kind == "IFCRELASSOCIATESCLASSIFICATION":  # ..., RelatedObjects, RelatingClassification
                for obj in a[4]:
                    slot(obj)["classes"].append(a[5])
            elif kind == "IFCMATERIALPROPERTIES":
                slot(a[3])["psets"].append(Ref(ref))
            elif kind.endswith("TYPE") and len(a) > 5 and isinstance(a[5], list):
                # IfcTypeObject: GlobalId, OwnerHistory, Name, Description, ApplicableOccurrence, HasPropertySets
                for p in a[5]:
                    slot(Ref(ref))["psets"].append(p)

        # element keys: type, name, and order among elements with the same type and name
        counts = {}
        keys = {}
        for ref in sorted(objects):
            kind, a = self.get(ref)
            name = a[0] if kind in ("IFCMATERIAL",) else (a[2] if len(a) > 2 else None)
            base = f'{kind} "{name}"' if name is not None else kind
            n = counts.get(base, 0)
            counts[base] = n + 1
            keys[ref] = f"{base} #{n}"

        lines = []
        for ref in sorted(objects, key=lambda r: keys[r]):
            key = keys[ref]
            kind, a = self.get(ref)
            if kind not in ("IFCMATERIAL",) and len(a) > 4 and not kind.endswith("TYPE"):
                lines.append(f"{key} | ObjectType = {self.value(a[4])}")
            for c in objects[ref]["classes"]:
                ck, ca = self.get(c)
                if ck == "IFCCLASSIFICATIONREFERENCE":  # Location, Identification, Name, ReferencedSource, ...
                    src = self.get(ca[3])[1]  # IfcClassification: Source, Edition, EditionDate, Name, ...
                    source = src[3] if src and len(src) > 3 else ""
                    lines.append(f"{key} | classification {source} {ca[1]} {self.value(ca[2])} {self.value(ca[0])}")
                else:
                    lines.append(f"{key} | classification {ck}")
            psets = sorted((self.pset_lines(p) for p in objects[ref]["psets"]), key=lambda t: (t[0], t[1]))
            for title, props in psets:
                if not props:
                    lines.append(f"{key} | {title} | (empty)")
                for p in props:
                    lines.append(f"{key} | {title} | {p}")
        return lines


def listing(path):
    """the listing of an IFC file, or a listing saved as .txt or .txt.gz"""
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".txt.gz"):
        return gzip.decompress(path.read_bytes()).decode("utf-8").splitlines()
    if name.endswith(".txt"):
        return path.read_text(encoding="utf-8").splitlines()
    return Listing(read_step(path)).dump()


def save(lines, path):
    """saves a listing as .txt, or .txt.gz"""
    path = Path(path)
    data = ("\n".join(lines) + "\n").encode("utf-8")
    if path.name.lower().endswith(".gz"):
        path.write_bytes(gzip.compress(data, mtime=0))
    else:
        path.write_bytes(data)


def compare(old, new, context=0, limit=200):
    """unified diff lines between two listings"""
    diff = list(difflib.unified_diff(old, new, "old", "new", n=context, lineterm=""))
    return diff[:limit] + ([f"... {len(diff) - limit} more lines"] if len(diff) > limit else [])


def count_differences(old, new):
    """number of lines only in one of the listings"""
    return sum(1 for line in difflib.unified_diff(old, new, n=0, lineterm="")
               if line[:1] in "+-" and not line.startswith(("+++", "---")))


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "dump":
        lines = listing(sys.argv[2])
        if len(sys.argv) > 3:
            save(lines, sys.argv[3])
        else:
            sys.stdout.write("\n".join(lines) + "\n")
    elif len(sys.argv) == 4 and sys.argv[1] == "compare":
        diff = compare(listing(sys.argv[2]), listing(sys.argv[3]))
        print("\n".join(diff) if diff else "same")
        sys.exit(1 if diff else 0)
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
