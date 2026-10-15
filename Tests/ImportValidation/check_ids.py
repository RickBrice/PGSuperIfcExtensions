"""Checks an IFC file against an IDS with ifctester and lists the failures.

Usage: python check_ids.py <file.ifc> <file.ids> [--spec NAME]

--spec checks only the specifications with that name. For each specification the applicable entity
classes are listed, then each failing requirement with the element and the reason.
"""
import argparse
import sys

import ifcopenshell

from check_standard import ifctester_ids


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ifc")
    parser.add_argument("ids")
    parser.add_argument("--spec", help="check only the specifications with this name")
    args = parser.parse_args()

    ids = ifctester_ids()
    ifc = ifcopenshell.open(args.ifc)
    spec = ids.open(args.ids)
    if args.spec:
        spec.specifications = [s for s in spec.specifications if s.name == args.spec]
        if not spec.specifications:
            sys.exit(f'No specification named "{args.spec}"')
    spec.validate(ifc)

    failures = 0
    for s in spec.specifications:
        classes = sorted({e.is_a() for e in s.applicable_entities})
        failed = [(r, f) for r in s.requirements for f in getattr(r, "failures", [])]
        failures += len(failed)
        if not failed and not args.spec:
            continue
        print(f'{s.name}: applies to {", ".join(classes) or "nothing"}')
        for r, f in failed:
            element = f["element"]
            what = getattr(r, "baseName", None) or type(r).__name__
            print(f'   FAIL {what}: {element.is_a()} "{getattr(element, "Name", "")}" - {f.get("reason")}')
    print(f"{len(spec.specifications)} specifications, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
