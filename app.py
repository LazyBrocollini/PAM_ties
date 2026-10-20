from __future__ import annotations

# Generated single-file deployment; rebuild with scripts/build_standalone.py


# --- src/data_utils.py ---
"""Spreadsheet normalization only; no reference sequence mapping."""
import math
import re
from decimal import Decimal, InvalidOperation
import pandas as pd

METRICS = ('135 Indel%', '45 Indel%', '11 Indel%', '0 Indel%')

def text(value):
    if value is None:
        return ''
    try:
        if pd.isna(value):
            return ''
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def identifier(value):
    value = text(value)
    if re.fullmatch(r'[+-]?\d+\.0+', value):
        return value.split('.')[0]
    return value


def spacer_text(value):
    # Used only as a workbook-join fallback, never for genomic searching.
    return re.sub(r'\s+', '', text(value)).upper()


def distinct(values):
    result = []
    for value in values:
        for token in text(value).split(';'):
            token = token.strip()
            if token and token not in result:
                result.append(token)
    return '; '.join(result)


def numeric_values(value):
    out = []
    for token in text(value).split(';'):
        try:
            number = float(token.strip().removesuffix('%').strip())
        except (ValueError, TypeError):
            continue
        if math.isfinite(number):
            out.append(number)
    return out


def numeric_indel(value):
    numbers = numeric_values(value)
    return max(numbers) if numbers else None


def integer_site(value):
    try:
        number = Decimal(text(value))
        if number.is_finite() and number >= 0 and number == number.to_integral_value():
            return int(number)
    except InvalidOperation:
        pass
    return None


def normalize_strand(value):
    value = text(value).lower()
    if value in ('+', '+1', '1', 'forward', 'plus', 'fwd'):
        return '+'
    if value in ('-', '-1', 'reverse', 'minus', 'rev'):
        return '-'
    return ''

# --- src/io_utils.py ---
from collections import defaultdict
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import BadZipFile
import pandas as pd

SPACER_SHEETS = ('all_spacer_data.registered_JAG1', 'all_spacer_data.registered_hJAG')
ALIASES = {
    'Ionis#': ('Ionis#', 'Ion_ID', 'Ion ID', 'guide_id'),
    'gRNA Name': ('gRNA Name', 'guide_name'),
    'spacer': ('spacer', 'Spacer Sequence'),
    'pam_site': ('pam_site', 'PAM site'),
    'nuclease': ('nuclease', 'enzyme', 'nuclease/enzyme'),
    'strand': ('strand', 'orientation'),
    'source_sheet': ('source_sheet', 'source dataset', 'dataset'),
    'reference_id': ('reference_id', 'assembly', 'reference'),
    'chromosome': ('chromosome', 'chrom', 'chr', 'contig'),
    'pam': ('pam', 'PAM sequence'),
    'cut_site': ('cut_site',),
    'Description': ('Description',),
    'Notes': ('Notes', 'Note'),
    **{metric: (metric, metric.replace('%', '')) for metric in METRICS},
}
REQUIRED = ('Ionis#', 'pam_site', 'nuclease', 'strand', 'source_sheet')


def _key(value):
    return ''.join(ch.lower() for ch in str(value) if ch.isalnum())


def _columns(frame):
    lookup = {_key(col): col for col in frame.columns}
    return {field: next((lookup[_key(alias)] for alias in aliases if _key(alias) in lookup), None)
            for field, aliases in ALIASES.items()}


def _get(row, columns, field):
    col = columns.get(field)
    return text(row.get(col, '')) if col is not None else ''


def _check(checks, category, detail, source_sheet='', source_row='', severity='warning'):
    checks.append(dict(category=category, severity=severity, source_sheet=source_sheet, source_row=source_row, detail=detail))


def _raw_records(frame, sheet, checks, flat):
    cols = _columns(frame)
    required = REQUIRED if flat else ('pam_site', 'nuclease', 'strand')
    missing = [field for field in required if cols.get(field) is None]
    if not flat and cols['Ionis#'] is None and cols['spacer'] is None:
        missing.append('Ionis# or spacer')
    if flat and not any(cols[metric] is not None for metric in METRICS):
        missing.append('at least one Indel% metric column')
    if missing:
        _check(checks, 'Missing columns', f'{sheet}: missing {", ".join(missing)}.', sheet, severity='error')
        return []
    records = []
    for index, row in frame.iterrows():
        if not any(text(v) for v in row):
            continue
        rec = {f'input__{col}': text(row[col]) for col in frame.columns}
        rec.update({field: _get(row, cols, field) for field in ALIASES})
        rec['source_sheet'] = rec['source_sheet'] or (sheet if not flat else '')
        rec['input_sheet'] = sheet
        rec['source_row'] = str(index + 2)
        rec['Ionis#'] = identifier(rec['Ionis#'])
        rec['spacer'] = spacer_text(rec['spacer'])
        records.append(rec)
    return records


def _consolidated(frames, checks):
    whitelist = frames['All Guides']
    columns = _columns(whitelist)
    if columns['Ionis#'] is None and columns['spacer'] is None:
        _check(checks, 'Missing columns', 'All Guides needs Ionis# or Spacer Sequence.', severity='error')
        return []
    by_id, by_spacer = defaultdict(list), defaultdict(list)
    whitelist_records = []
    for index, row in whitelist.iterrows():
        if not any(text(v) for v in row):
            continue
        rec = {field: _get(row, columns, field) for field in ALIASES}
        rec.update({f'all_guides__{col}': text(row[col]) for col in whitelist.columns})
        rec['Ionis#'] = identifier(rec['Ionis#'])
        rec['spacer'] = spacer_text(rec['spacer'])
        rec['all_guides_excel_rows'] = str(index + 2)
        whitelist_records.append(rec)
        if rec['Ionis#']:
            by_id[rec['Ionis#']].append(rec)
        if rec['spacer']:
            by_spacer[rec['spacer']].append(rec)
    used, records = set(), []
    for sheet in SPACER_SHEETS:
        if sheet not in frames:
            _check(checks, 'Missing optional sheet', f'{sheet} is not present.')
            continue
        for rec in _raw_records(frames[sheet], sheet, checks, flat=False):
            def in_scope(candidate):
                scope = candidate['source_sheet'].casefold()
                return not scope or scope == rec['source_sheet'].casefold()
            candidates = [candidate for candidate in by_id.get(rec['Ionis#'], []) if in_scope(candidate)]
            method = 'Ionis#'
            if not candidates and rec['spacer']:
                candidates = [candidate for candidate in by_spacer.get(rec['spacer'], []) if in_scope(candidate)]
                method = 'Spacer Sequence fallback'
            if not candidates:
                continue  # Never admit spacer-data records outside the whitelist.
            used.update(candidate['all_guides_excel_rows'] for candidate in candidates)
            # A spacer fallback can match distinct IDs; keep their indel values separate.
            candidate_groups = defaultdict(list)
            for candidate in candidates:
                candidate_groups[candidate['Ionis#'] or rec['Ionis#']].append(candidate)
            for guide_id, associated in candidate_groups.items():
                joined = dict(rec)
                for field in ('gRNA Name', 'Description', *METRICS, 'Notes'):
                    joined[field] = distinct(candidate[field] for candidate in associated)
                joined['Ionis#'] = guide_id
                joined['spacer'] = joined['spacer'] or distinct(candidate['spacer'] for candidate in associated)
                joined['pam'] = joined['pam'] or distinct(candidate['pam'] for candidate in associated)
                joined['all_guides_excel_rows'] = distinct(candidate['all_guides_excel_rows'] for candidate in associated)
                joined['join_method'] = method
                for key in associated[0]:
                    if key.startswith('all_guides__'):
                        joined[key] = distinct(candidate.get(key) for candidate in associated)
                records.append(joined)
    for rec in whitelist_records:
        if rec['all_guides_excel_rows'] not in used:
            _check(checks, 'Unmatched All Guides row', f"No spacer-data join for Ionis# {rec['Ionis#'] or '(blank)'}. Check ID, spacer, and source dataset.", 'All Guides', rec['all_guides_excel_rows'])
    return records


def _finalize(records, checks):
    valid, excluded = [], []
    for rec in records:
        rec = dict(rec)
        problems = []
        rec['pam_site_raw'] = rec['pam_site']
        site = integer_site(rec['pam_site'])
        if site is None:
            problems.append('pam_site must be a nonnegative integer, not a range or free text')
        for field in ('Ionis#', 'nuclease', 'source_sheet'):
            if not text(rec[field]):
                problems.append(f'missing {field}')
        strand = normalize_strand(rec['strand'])
        if not strand:
            problems.append('missing or unrecognized strand')
        rec['strand_original'] = rec['strand']
        rec['strand'] = strand
        if problems:
            rec['exclusion_reason'] = '; '.join(problems)
            excluded.append(rec)
            _check(checks, 'Excluded row', rec['exclusion_reason'], rec['source_sheet'], rec['source_row'])
            continue
        rec['pam_site'] = site
        if not rec['spacer']:
            _check(checks, 'Missing spacer (informational)', 'Spacer is optional; this record remains eligible.', rec['source_sheet'], rec['source_row'], 'info')
        valid.append(rec)
    # Include ID, dataset and coordinate context: never collapse distinct guide IDs.
    buckets = defaultdict(list)
    for rec in valid:
        key = tuple(text(rec.get(field)).casefold() for field in
                    ('source_sheet', 'reference_id', 'chromosome', 'Ionis#', 'spacer', 'nuclease', 'strand', 'pam_site', 'cut_site'))
        buckets[key].append(rec)
    merged = []
    for key, rows in buckets.items():
        result = dict(rows[0])
        for col in set().union(*(row.keys() for row in rows)):
            if col != 'pam_site':
                result[col] = distinct(row.get(col) for row in rows)
        result['pam_site'] = rows[0]['pam_site']
        result['record_id'] = 'R-' + sha256(repr(key).encode()).hexdigest()[:20]
        if len(rows) > 1:
            _check(checks, 'Duplicate records merged', f'{len(rows)} records merged; all distinct values and row numbers retained.', result['source_sheet'], result['source_row'], 'info')
        for metric in METRICS:
            if numeric_indel(result.get(metric)) is None:
                _check(checks, 'Missing numeric Indel%', metric, result['source_sheet'], result['source_row'])
            elif len(set(numeric_values(result.get(metric)))) > 1:
                _check(checks, 'Multiple Indel% values', f'{metric}: retained all values; ranking uses their maximum for this same record.', result['source_sheet'], result['source_row'], 'info')
        merged.append(result)
    guides = pd.DataFrame(merged)
    if not guides.empty:
        guides = guides.sort_values(['source_sheet', 'nuclease', 'strand', 'pam_site', 'record_id'], kind='stable').reset_index(drop=True)
        for _, rows in guides.groupby(['source_sheet', 'Ionis#'], dropna=False):
            if len(rows) > 1:
                _check(checks, 'Repeated guide identifier', f"Ionis# {rows.iloc[0]['Ionis#']} occurs at multiple sites or annotations; distinct record IDs retained.", rows.iloc[0]['source_sheet'])
    return guides, pd.DataFrame(excluded)


def sheet_names(content):
    with pd.ExcelFile(BytesIO(content), engine='openpyxl') as workbook:
        return workbook.sheet_names


def load_input(content, filename, sheet=None):
    checks = []
    mode, records = 'flat table', []
    try:
        if Path(filename).suffix.lower() == '.csv':
            frame = pd.read_csv(BytesIO(content), dtype=object, keep_default_na=False)
            records = _raw_records(frame, 'CSV', checks, flat=True)
            sheets = ['CSV']
        else:
            with pd.ExcelFile(BytesIO(content), engine='openpyxl') as workbook:
                sheets = workbook.sheet_names
                if 'All Guides' in sheets and any(name in sheets for name in SPACER_SHEETS):
                    mode = 'consolidated whitelist workbook'
                    frames = {name: pd.read_excel(workbook, sheet_name=name, dtype=object, keep_default_na=False)
                              for name in ('All Guides', *SPACER_SHEETS) if name in sheets}
                    records = _consolidated(frames, checks)
                else:
                    selected = sheet or ('Guides' if 'Guides' in sheets else sheets[0])
                    frame = pd.read_excel(workbook, sheet_name=selected, dtype=object, keep_default_na=False)
                    records = _raw_records(frame, selected, checks, flat=True)
        guides, excluded = _finalize(records, checks)
    except (ValueError, OSError, KeyError, pd.errors.ParserError, EOFError, BadZipFile, UnicodeDecodeError) as error:
        _check(checks, 'Input error', str(error), severity='error')
        guides, excluded, sheets = pd.DataFrame(), pd.DataFrame(), []
    return dict(guides=guides, excluded=excluded, checks=pd.DataFrame(checks, columns=['category', 'severity', 'source_sheet', 'source_row', 'detail']), mode=mode, sheets=sheets)

# --- src/nearby_pam.py ---
from collections import defaultdict
from itertools import product
import pandas as pd

SCOPE = ('source_sheet', 'reference_id', 'chromosome', 'nuclease', 'strand')


def find_nearby(guides):
    """Compare spreadsheet pam_site integers directly; never map sequence."""
    if guides.empty:
        return pd.DataFrame(), pd.DataFrame(), {}, {}
    indexed = guides.set_index('record_id', drop=False)
    parent = {record: record for record in indexed.index}
    def find(record):
        while parent[record] != record:
            parent[record] = parent[parent[record]]
            record = parent[record]
        return record
    buckets = defaultdict(lambda: defaultdict(list))
    for record, row in indexed.iterrows():
        key = tuple(str(row.get(field, '')).strip().casefold() for field in SCOPE)
        buckets[key][int(row['pam_site'])].append(record)
    pairs = []
    for key in sorted(buckets):
        sites = buckets[key]
        for site in sorted(sites):
            for a, b in product(sorted(sites[site]), sorted(sites.get(site + 1, []))):
                parent[find(b)] = find(a)
                pair = dict(pair_id=f'P{len(pairs)+1:04}', record_id_a=a, record_id_b=b, pam_distance_nt=1)
                for side, record in (('a', a), ('b', b)):
                    for field, value in indexed.loc[record].items():
                        pair[f'{field}_{side}'] = value
                pairs.append(pair)
    components = defaultdict(list)
    involved = {pair[f'record_id_{side}'] for pair in pairs for side in ('a', 'b')}
    for record in sorted(involved):
        components[find(record)].append(record)
    groups, record_to_group, record_to_pairs = [], {}, defaultdict(list)
    for members in sorted(components.values(), key=lambda members: sorted(members)):
        group_id = f'G{len(groups)+1:04}'
        rows = indexed.loc[members]
        groups.append(dict(group_id=group_id, record_ids='; '.join(members), group_size=len(members),
                           pam_site_min=int(rows['pam_site'].min()), pam_site_max=int(rows['pam_site'].max()),
                           **{field: rows.iloc[0].get(field, '') for field in SCOPE}))
        record_to_group.update({record: group_id for record in members})
    for pair in pairs:
        pair['group_id'] = record_to_group[pair['record_id_a']]
        for side in ('a', 'b'):
            record_to_pairs[pair[f'record_id_{side}']].append(pair['pair_id'])
    return pd.DataFrame(pairs), pd.DataFrame(groups), record_to_group, dict(record_to_pairs)

# --- src/ranking.py ---
import pandas as pd


def rank_results(guides, pairs, groups, record_to_group, record_to_pairs, metric):
    if metric not in METRICS:
        raise ValueError('Unrecognized Indel% metric')
    guides, pairs, groups = guides.copy(), pairs.copy(), groups.copy()
    if guides.empty:
        return guides, pairs, groups, pd.DataFrame()
    guides['selected_indel_metric'] = metric
    guides['selected_indel_value'] = [numeric_indel(value) for value in guides[metric]]
    guides['group_id'] = guides['record_id'].map(record_to_group).fillna('')
    guides['pair_id'] = guides['record_id'].map(lambda record: '; '.join(record_to_pairs.get(record, [])))
    guides['pam_distance_nt'] = guides['record_id'].map(lambda record: 1 if record_to_pairs.get(record) else None)
    guides['is_pair_or_group_winner'] = False
    guides['winner_reason'] = ''
    records = guides.set_index('record_id', drop=False)
    best, reasons = [], {record: [] for record in records.index}
    for entity_type, frame in (('pair', pairs), ('group', groups)):
        for index, entity in frame.iterrows():
            entity_id = entity[f'{entity_type}_id']
            members = [entity['record_id_a'], entity['record_id_b']] if entity_type == 'pair' else entity['record_ids'].split('; ')
            values = {record: numeric_indel(records.loc[record, metric]) for record in members}
            usable = {record: value for record, value in values.items() if value is not None}
            highest = max(usable.values()) if usable else None
            winners = [record for record, value in usable.items() if value == highest]
            reason = ('No numeric Indel% available' if not winners else
                      'Tie for highest selected Indel%' if len(winners) > 1 else 'Highest selected Indel%')
            frame.at[index, 'selected_indel_metric'] = metric
            frame.at[index, 'winner_record_ids'] = '; '.join(winners)
            frame.at[index, 'winner_reason'] = reason
            frame.at[index, 'Winner'] = '; '.join(str(records.loc[record, 'gRNA Name'] or records.loc[record, 'Ionis#']) for record in winners) if winners else reason
            frame.at[index, 'winning_indel_value'] = highest
            if entity_type == 'pair':
                for side in ('a', 'b'):
                    frame.at[index, f'selected_indel_value_{side}'] = values[entity[f'record_id_{side}']]
            for record in winners or [None]:
                row = records.loc[record].to_dict() if record else {}
                row.update(entity_type=entity_type, entity_id=entity_id, winner_record_id=record or '',
                           selected_indel_metric=metric, selected_indel_value=highest, winner_reason=reason,
                           is_pair_or_group_winner=bool(record))
                best.append(row)
                if record:
                    reasons[record].append(f'{entity_id}: {reason}')
            if not winners:
                for record in members:
                    reasons[record].append(f'{entity_id}: {reason}')
    winning = {row['winner_record_id'] for row in best if row['winner_record_id']}
    guides['is_pair_or_group_winner'] = guides['record_id'].isin(winning)
    guides['winner_reason'] = guides['record_id'].map(lambda record: '; '.join(reasons[record]))
    return guides, pairs, groups, pd.DataFrame(best)

# --- src/templates.py ---
"""Synthetic example data: not experimental results."""
import pandas as pd

REQUIRED_COLUMNS = ['Ionis#', 'pam_site', 'nuclease', 'strand', 'source_sheet', '135 Indel%']


def example_table():
    rows = [
        ('DEMO001', 'Group A low', 100, 'SpCas9', '+', 'DEMO_A', 'chr1', 10, 5, 2, 0, '100 and 102 are NOT a pair; connected through 101.'),
        ('DEMO002', 'Group A best', 101, 'SpCas9', '+', 'DEMO_A', 'chr1', 40, 20, 9, 1, 'Highest 135 Indel% in the three-guide connected group.'),
        ('DEMO003', 'Group A mid', 102, 'SpCas9', '+', 'DEMO_A', 'chr1', 25, 12, 6, 0, 'Neighbor of 101 only.'),
        ('DEMO004', 'Tie reverse A', 200, 'SpCas9', '-', 'DEMO_A', 'chr1', 55, 15, 7, 1, 'Tie at 135 Indel%.'),
        ('DEMO005', 'Tie reverse B', 201, 'SpCas9', '-', 'DEMO_A', 'chr1', 55, 30, 8, 0, 'Tie at 135 Indel%; winner changes at 45 Indel%.'),
        ('DEMO006', 'Missing A', 300, 'Cas12a', '+', 'DEMO_A', 'chr1', '', '', '', '', 'Pair has no numeric 135 Indel%.'),
        ('DEMO007', 'Missing B', 301, 'Cas12a', '+', 'DEMO_A', 'chr1', 'not measured', '', '', '', 'Invalid metric is treated as missing.'),
        ('DEMO008', 'Isolated', 400, 'SpCas9', '+', 'DEMO_A', 'chr1', 30, 10, 5, 0, '401 has a different nuclease, so no pair.'),
        ('DEMO009', 'Different nuclease', 401, 'Cas12a', '+', 'DEMO_A', 'chr1', 90, 40, 10, 1, 'Different nuclease and source restrictions apply.'),
        ('DEMO010', 'Different dataset', 402, 'Cas12a', '+', 'DEMO_B', 'chr1', 80, 30, 9, 1, 'Not paired with 401: different source dataset.'),
        ('DEMO011', 'Different chromosome', 101, 'SpCas9', '+', 'DEMO_A', 'chr2', 99, 50, 20, 1, 'Never compared with chr1 coordinates.'),
    ]
    columns = ['Ionis#', 'gRNA Name', 'pam_site', 'nuclease', 'strand', 'source_sheet', 'chromosome', '135 Indel%', '45 Indel%', '11 Indel%', '0 Indel%', 'Notes']
    frame = pd.DataFrame(rows, columns=columns)
    frame['reference_id'] = 'DEMO_ASSEMBLY'
    frame['spacer'] = ''
    frame['pam'] = ''
    return frame


def schema_table():
    return pd.DataFrame([
        ('Ionis#', 'Required', 'Guide identifier. Ion_ID or guide_id are accepted aliases.'),
        ('pam_site', 'Required', 'Nonnegative integer from your spreadsheet. Used directly, without sequence validation or coordinate conversion.'),
        ('nuclease', 'Required', 'Enzyme label. Only identical labels (case-insensitive) are compared.'),
        ('strand', 'Required', '+ or -. Forward/reverse and 1/-1 are accepted aliases.'),
        ('source_sheet', 'Required in flat input', 'Dataset identifier; guides from different datasets are never paired.'),
        ('135 Indel%', 'At least one metric column', 'Numeric values such as 40 represent 40 percentage points. Blank/text values are missing.'),
        ('45 Indel%', 'Optional metric', 'Use the same numeric scale as the other metric columns.'),
        ('11 Indel%', 'Optional metric', 'Use the same numeric scale as the other metric columns.'),
        ('0 Indel%', 'Optional metric', 'Use the same numeric scale as the other metric columns.'),
        ('gRNA Name', 'Optional', 'Human-readable guide name.'),
        ('chromosome', 'Recommended if multiple chromosomes', 'Additional comparison boundary; different values are never paired.'),
        ('reference_id', 'Recommended if multiple assemblies/loci', 'Additional comparison boundary; use one consistent coordinate system per value.'),
        ('spacer', 'Optional', 'Metadata only, except as the original consolidated-workbook join fallback. No reference matching.'),
        ('pam', 'Optional', 'Metadata only; not extracted or sequence-validated.'),
        ('Notes', 'Optional', 'Free-text comments preserved in output.'),
    ], columns=['Column', 'Requirement', 'Meaning'])

# --- src/exports.py ---
from io import BytesIO
import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment


def safe_table(frame):
    out = frame.copy()
    def safe(value):
        if isinstance(value, str) and len(value) > 1 and value.lstrip().startswith(('=', '+', '-', '@')):
            return "'" + value
        return value
    for col in out:
        out[col] = out[col].map(safe)
    return out


def csv_bytes(frame):
    return safe_table(frame).to_csv(index=False).encode('utf-8-sig')


def excel_bytes(tables):
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        for name, frame in tables.items():
            safe_table(frame).to_excel(writer, sheet_name=name[:31], index=False)
            ws = writer.sheets[name[:31]]
            ws.freeze_panes = 'A2'
            for cell in ws[1]:
                cell.font = Font(name='Arial', bold=True, color='FFFFFF')
                cell.fill = PatternFill('solid', fgColor='145A48')
                cell.alignment = Alignment(wrap_text=True, vertical='center')
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    cell.font = Font(name='Arial', size=10)
                    cell.alignment = Alignment(vertical='top')
            for col in ws.columns:
                width = min(52, max(13, max(len(str(cell.value or '')) for cell in list(col)[:100]) + 2))
                ws.column_dimensions[col[0].column_letter].width = width
            ws.row_dimensions[1].height = 32
    return buffer.getvalue()

# --- Streamlit interface ---
import hashlib
import sys
from pathlib import Path
import pandas as pd
import streamlit as st

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


st.set_page_config(page_title='Nearby PAM Guide Selector', page_icon='🧬', layout='wide')


@st.cache_data(show_spinner=False)
def cached_load(content, filename, sheet):
    return load_input(content, filename, sheet)


@st.cache_data(show_spinner=False)
def cached_sheets(content):
    return sheet_names(content)


@st.cache_data(show_spinner=False)
def cached_analyze(guides, metric):
    pairs, groups, record_to_group, record_to_pairs = find_nearby(guides)
    return rank_results(guides, pairs, groups, record_to_group, record_to_pairs, metric)


def display_table(frame, best=False):
    if frame.empty:
        st.info('No rows to display.')
        return
    leading = ['entity_type', 'entity_id', 'group_id', 'pair_id', 'gRNA Name', 'Ionis#', 'pam_site',
               'nuclease', 'strand', 'selected_indel_metric', 'selected_indel_value', 'source_sheet', 'source_row',
               'Winner', 'winner_reason', 'pam_distance_nt']
    ordered = [col for col in leading if col in frame] + [col for col in frame if col not in leading]
    frame = frame[ordered]
    if best and 'selected_indel_value' in frame:
        try:
            frame = frame.style.map(lambda value: 'background-color:#D8F3E6;color:#145A48;font-weight:bold' if pd.notna(value) else '', subset=['selected_indel_value'])
        except (ImportError, AttributeError):
            st.success('Winner results below retain all ties. Cell color formatting is unavailable in this environment.')
    st.dataframe(frame, use_container_width=True, hide_index=True)


def main():
    st.title('Nearby PAM Guide Selector')
    st.caption('Spreadsheet-only analysis: compare pam_site values exactly 1 nucleotide apart and rank nearby pairs and connected groups by Indel%. No reference sequence or spacer mapping is used.')
    with st.expander('Input columns and example files'):
        st.write('Required flat-table columns: Ionis#, pam_site, nuclease, strand, source_sheet, and at least one Indel% metric column. Spacer sequence is optional.')
        st.warning('pam_site values must already use a consistent coordinate convention within each dataset/reference/chromosome. The app does not validate or convert genomic coordinates.')
        st.download_button('Download example Excel', excel_bytes({'Guides': example_table(), 'Column Guide': schema_table()}), 'example_pam_site_guides.xlsx', key='template_excel')
        st.download_button('Download example CSV', csv_bytes(example_table()), 'example_pam_site_guides.csv', 'text/csv', key='template_csv')
        display_table(schema_table())

    with st.sidebar:
        st.header('Spreadsheet input')
        upload = st.file_uploader('Guide data (.xlsx or .csv)', type=['xlsx', 'csv'])
        use_example = st.checkbox('Use synthetic example data', value=False)
    if upload is not None:
        content, filename = upload.getvalue(), upload.name
        if use_example:
            st.info('The uploaded file takes precedence over synthetic example data.')
    elif use_example:
        content, filename = csv_bytes(example_table()), 'example_pam_site_guides.csv'
    else:
        default = PROJECT_ROOT / 'fixtures' / 'Sum_samples_guides_hJAG1_dKO_consolidated.xlsx'
        if default.exists():
            content, filename = default.read_bytes(), default.name
        else:
            st.info('Upload a workbook/CSV or enable synthetic example data in the sidebar.')
            return
    selected_sheet = None
    if filename.lower().endswith('.xlsx'):
        try:
            names = cached_sheets(content)
        except Exception as error:
            st.error(f'Cannot open Excel file: {error}')
            return
        if not ('All Guides' in names and any(name in names for name in SPACER_SHEETS)):
            selected_sheet = st.sidebar.selectbox('Worksheet', names, index=names.index('Guides') if 'Guides' in names else 0)
    loaded = cached_load(content, filename, selected_sheet)
    raw_guides = loaded['guides']
    st.caption(f"File: {filename} · Input mode: {loaded['mode']} · Worksheets: {', '.join(loaded['sheets'])}")
    if raw_guides.empty:
        st.error('No eligible guide records were loaded. Review the data checks below.')
        display_table(loaded['checks'])
        display_table(loaded['excluded'])
        return
    with st.sidebar:
        metric = st.selectbox('Indel% metric', METRICS)
        filters = {}
        for col, label in [('nuclease', 'Nuclease/enzyme'), ('strand', 'Strand'), ('source_sheet', 'Source dataset')]:
            filters[col] = st.multiselect(label, sorted(raw_guides[col].unique()))
        analyze = st.button('Analyze nearby PAM guides', type='primary', use_container_width=True)
    guides = raw_guides.copy()
    for col, selected in filters.items():
        if selected:
            guides = guides[guides[col].isin(selected)]
    signature = hashlib.sha256(content + repr((filename, selected_sheet, metric, filters)).encode()).hexdigest()
    if analyze:
        st.session_state['pam_site_analysis'] = (signature, cached_analyze(guides, metric))
    stored = st.session_state.get('pam_site_analysis')
    if stored is None or stored[0] != signature:
        st.info(f'{len(guides)} eligible guide records after filters. Select Analyze nearby PAM guides.')
        with st.expander('Input data checks'):
            display_table(loaded['checks'])
        return
    ranked, pairs, groups, best = stored[1]
    metrics = st.columns(4)
    for column, label, value in zip(metrics, ['Eligible guides', 'Guides in nearby groups', 'Exact 1-nt pairs', 'Connected groups'],
                                    [len(ranked), ranked['group_id'].ne('').sum() if not ranked.empty else 0, len(pairs), len(groups)]):
        column.metric(label, int(value))
    st.info('Rule: abs(pam_site A − pam_site B) == 1. Same dataset, nuclease, strand, reference_id and chromosome required. Distances of 0 or 2 nt are excluded.')
    st.caption('Indel values are used on their supplied numeric scale. Multiple values for the same deduplicated guide record are preserved; its ranking value is the highest numeric value. Ties are retained.')
    if pairs.empty:
        st.warning('No exact 1-nt nearby-PAM pairs were found after filtering.')
    best_tab, pairs_tab, groups_tab, guides_tab, checks_tab = st.tabs(['Best guides', 'Nearby pairs', 'Connected groups', 'All guides', 'Data checks'])
    with best_tab:
        kind = st.radio('Show winners for', ['Groups and pairs', 'Groups only', 'Pairs only'], horizontal=True)
        shown = best if kind == 'Groups and pairs' or best.empty else best[best['entity_type'].eq('group' if kind == 'Groups only' else 'pair')]
        display_table(shown, best=True)
    with pairs_tab:
        display_table(pairs)
    with groups_tab:
        display_table(groups)
    with guides_tab:
        search = st.text_input('Filter by guide name, ID, spacer or note')
        shown = ranked
        if search and not ranked.empty:
            mask = pd.Series(False, index=ranked.index)
            for col in ['gRNA Name', 'Ionis#', 'spacer', 'Notes']:
                mask |= ranked[col].astype(str).str.contains(search, case=False, regex=False, na=False)
            shown = ranked[mask]
        sort_by = st.selectbox('Sort guides by', ['pam_site', 'selected_indel_value'])
        display_table(shown.sort_values(sort_by, ascending=sort_by == 'pam_site', na_position='last') if not shown.empty else shown)
    with checks_tab:
        st.warning('No reference matching, PAM extraction, strand inference or sequence-based validation is performed. All coordinates and strand annotations come from the input file.')
        display_table(loaded['checks'])
        if not loaded['excluded'].empty:
            st.subheader('Excluded records')
            display_table(loaded['excluded'])
    st.subheader('Download analysis')
    st.download_button('Download guide results CSV', csv_bytes(ranked), 'pam_site_guide_results.csv', 'text/csv')
    st.download_button('Download nearby pairs CSV', csv_bytes(pairs), 'pam_site_pairs.csv', 'text/csv')
    st.download_button('Download Excel results', excel_bytes({'Guides': ranked, 'Nearby Pairs': pairs, 'Connected Groups': groups, 'Best Guides': best, 'Data Checks': loaded['checks'], 'Excluded Records': loaded['excluded']}),
                       'pam_site_results.xlsx', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')


if __name__ == '__main__':
    main()
