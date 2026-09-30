"""Compares design-value IDS files (buildingSMART Information Delivery Specification).

An IDS file is reduced to a canonical listing, one line per facet of each specification:

   spec "Girder design values - Span 1, Girder A" | requirement property Pset_PrecastConcreteElementGeneral.ReleaseStrength ...

The facets of a specification are sorted, so their order doesn't matter, and the date in <ids:info>
is left out (it's the day the IDS was written). Used to check that the design-value IDS exporter
writes the same specifications when its property locations come from the mapping table
(devdocs/MappingTablesDesign.md, M4).

Usage:
   python compare_ids.py dump design.ids [listing.txt]
   python compare_ids.py compare old.ids|old.txt new.ids     exit code 1 if they differ
"""
import difflib
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def local(tag):
    return tag.split('}', 1)[-1]


def text_of(element):
    """a facet as text: its tag, attributes, and children, recursively, with namespaces removed"""
    attributes = ' '.join(f'{local(k)}="{v}"' for k, v in sorted(element.attrib.items()))
    children = [text_of(child) for child in element]
    text = (element.text or '').strip()
    inner = ' '.join(filter(None, [text] + children))
    return f"{local(element.tag)}{('[' + attributes + ']') if attributes else ''}({inner})"


def listing(path):
    path = Path(path)
    if path.suffix.lower() == '.txt':
        return path.read_text(encoding='utf-8').splitlines()

    root = ET.parse(path).getroot()
    lines = []
    for child in root:
        if local(child.tag) == 'info':
            for item in child:
                if local(item.tag) != 'date':
                    lines.append(f"info | {text_of(item)}")
        elif local(child.tag) == 'specifications':
            for spec in child:
                name = spec.attrib.get('name', '')
                attributes = ' '.join(f'{local(k)}="{v}"' for k, v in sorted(spec.attrib.items()) if local(k) != 'name')
                facets = []
                for section in spec:
                    kind = local(section.tag)  # applicability, requirements
                    if kind in ('applicability', 'requirements'):
                        section_attributes = ' '.join(f'{local(k)}="{v}"' for k, v in sorted(section.attrib.items()))
                        if section_attributes:
                            facets.append(f"{kind} [{section_attributes}]")
                        facets += [f"{kind} {text_of(facet)}" for facet in section]
                    else:
                        facets.append(text_of(section))
                lines.append(f'spec "{name}" | [{attributes}]')
                lines += [f'spec "{name}" | {facet}' for facet in sorted(facets)]
    return lines


def save(lines, path):
    Path(path).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def compare(old, new, context=0, limit=200):
    diff = list(difflib.unified_diff(old, new, 'old', 'new', n=context, lineterm=''))
    return diff[:limit] + ([f'... {len(diff) - limit} more lines'] if len(diff) > limit else [])


def count_differences(old, new):
    return sum(1 for line in difflib.unified_diff(old, new, n=0, lineterm='')
               if line[:1] in '+-' and not line.startswith(('+++', '---')))


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == 'dump':
        lines = listing(sys.argv[2])
        if len(sys.argv) > 3:
            save(lines, sys.argv[3])
        else:
            sys.stdout.write('\n'.join(lines) + '\n')
    elif len(sys.argv) == 4 and sys.argv[1] == 'compare':
        diff = compare(listing(sys.argv[2]), listing(sys.argv[3]))
        print('\n'.join(diff) if diff else 'same')
        sys.exit(1 if diff else 0)
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == '__main__':
    main()
