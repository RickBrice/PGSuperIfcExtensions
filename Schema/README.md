# IDS schema binding

buildingSMART **IDS (Information Delivery Specification) 1.0.0** schema and the
CodeSynthesis xsd-cxx C++/Tree binding generated from it. Used by `CIdsExporter`
(`IdsExporter.cpp`) to write a project-specific `.ids` for a bridge's girder concrete
strength design values.

| File | What it is |
|---|---|
| `ids.xsd` | The **pristine, official** IDS 1.0.0 schema (`targetNamespace http://standards.buildingsmart.org/IDS`), from <https://github.com/buildingSMART/IDS>. Reference / validation gate only — **not** compiled. |
| `ids-binding.xsd` | A locally-adapted copy of `ids.xsd` that `xsd.exe` can process. This is what the build compiles (see the `CXX_Tree_Mapping_Rule` item in `PGSuperIfcExtensions.vcxproj`). |
| `ids-binding.hxx` / `ids-binding.cxx` | **Generated** by the `xsd cxx-tree` build step from `ids-binding.xsd`. Checked in for IntelliSense and first-build convenience; regenerated automatically when `ids-binding.xsd` or the project file changes. |

## Why a locally-adapted schema

`xsd.exe` cannot consume the pristine `ids.xsd` directly. Two constructs pull in the
XML-Schema-for-schemas (`<xs:import namespace="http://www.w3.org/2001/XMLSchema">`):

1. `idsValue` offers `<xs:element ref="xs:restriction"/>` as an alternative to
   `<ids:simpleValue>` (this is how an IDS embeds `xs:enumeration` / `xs:pattern` /
   `xs:minInclusive` value constraints).
2. `applicabilityType` uses `<xs:attributeGroup ref="xs:occurs"/>` for its
   `minOccurs` / `maxOccurs` attributes.

`ids-binding.xsd` removes the three `<xs:import>`s and replaces those two spots with:

1. `<xs:any namespace="##other" processContents="skip"/>` — xsd-cxx maps this to a
   `xercesc::DOMElement` wildcard (`--generate-wildcard`). `CIdsExporter` builds the
   `xs:restriction` subtree by hand as DOM (only needed for the non-default `Minimum`
   strength-match mode).
2. Plain optional `minOccurs` / `maxOccurs` string attributes.

Every other line is identical to `ids.xsd`. Output written by `CIdsExporter` is still
validated against the pristine `ids.xsd` — that is the real conformance gate.

## Regenerating by hand

```
cd Schema
"%XSDDIR%\bin\xsd.exe" cxx-tree --generate-serialization --generate-wildcard ^
   --root-element ids --char-type char --output-dir . ids-binding.xsd
```
