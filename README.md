# Nearby PAM Guide Selector

A production-oriented Streamlit app for finding approved guides from a consolidated Excel workbook in a supplied reference DNA sequence, comparing exact 1-nt PAM-anchor neighbors, and ranking each nearby pair/connected group by a selected Indel% metric.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Upload `Sum_samples_guides_hJAG1_dKO_consolidated.xlsx` in the sidebar. If a file with that name is placed in `fixtures/`, the app auto-loads it as the default fixture; the workbook path is never hardcoded as a required input.

Run tests:

```bash
pytest -q
```

## Expected workbook schema

Required sheet: `All Guides`. It should contain an Ionis#/Ion_ID column and/or a Spacer Sequence column. The app preserves the requested whitelist fields when present: `Source Sheet`, `Source Row`, `gRNA Name`, `Ionis#`, `Description`, `135 Indel%`, `45 Indel%`, `11 Indel%`, `0 Indel%`, `Spacer Sequence`, `PAM`, and `Notes`.

Optional sheets:

- `all_spacer_data.registered_JAG1`
- `all_spacer_data.registered_hJAG`

Spacer data is joined first by normalized Ionis#/Ion_ID and falls back to normalized spacer sequence. Excel numeric IDs such as `2000471.0` normalize to `2000471`. Duplicate All Guides rows are aggregated with semicolon-delimited distinct values and their source rows are retained. Spacer records are deduplicated by spacer, nuclease, strand, pam_site, and cut_site while retaining source-row provenance.

## Coordinates and PAM convention

Reference coordinates are 1-based and inclusive. The app searches every forward spacer and reverse-complement hit. A forward guide is reported as `+`; a reverse-complement hit is `-`.

The single `pam_anchor` coordinate is the PAM base immediately adjacent to the protospacer on the reference sequence. For a forward guide with a 3′ PAM, this is the PAM start; for a reverse guide with a 3′ PAM, this is the PAM end. The analogous adjacent base is used for 5′ PAMs. `pam_site` from the workbook is retained as source metadata and is not used as a reference coordinate.

Nearby pairs require all of the following: same source sheet, same nuclease/enzyme, same detected strand, and `abs(pam_anchor_a - pam_anchor_b) == 1`. PAMs 2 nt apart are deliberately excluded. Connected components of pair edges form nearby groups.

## Edge cases and warnings

The UI and Data checks tab report missing sheets/columns, unmatched whitelist rows, missing spacer sequences, missing or invalid Indel% values, multiple hits, ambiguous N bases, missing/out-of-range PAMs, and PAM-validation failures. Missing or nonnumeric selected Indel% values do not win; ties are explicitly labeled.

## Project layout

```text
app.py
src/
  io_utils.py
  sequence_utils.py
  guide_matching.py
  nearby_pam.py
  ranking.py
  visualization.py
tests/
  test_sequence_utils.py
  test_nearby_pam.py
  test_ranking.py
  test_io_utils.py
requirements.txt
README.md
```
