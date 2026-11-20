#!/usr/bin/env python3
"""A1 design table: one row per representative strain (K order) with keys, analysis flags, covariates, CrAg platform and grade. Standard library only.

Inputs (relative to the project root): admissions and isolates tables, the errata table, taxa.tsv and K.tsv of the G1 output. Output: design.tsv next to this script.
Usage (from the project root): python3 association/build_design.py   (optional: --root, --out-dir)
"""
import argparse
import collections
import csv
import datetime
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

ADMISSIONS = 'data/admissions.csv'
ISOLATES = 'data/isolates.csv'
G1_OUT = 'results/g1'   # G1 output directory, relative to the project root
TAXA = G1_OUT + '/03_l1/taxa.tsv'
KMAT = G1_OUT + '/03_l1/K.tsv'
ERRATA = 'data/errata_clinical.tsv'

# Cn13 is excluded from the A1 main analysis (sequenced DNA is a two-genotype mixture)
A1_EXCLUDED = {'Cn13': 'D1b: sequenced DNA is a ~6:4 two-genotype mixture (D1a)'}
# Lab-confirmed CrAg platforms, keyed by (admission_id, 'csf' | 'blood'); empty until confirmations are recorded
CONFIRMED_PLATFORM = {}
# provisional transition window for presumed platforms (exploratory operating rule)
WINDOW = (datetime.date(2016, 8, 1), datetime.date(2016, 8, 31))
SERIES_2 = {2 ** k for k in range(1, 14)}
SERIES_5 = {5 * 2 ** k for k in range(0, 11)}
# Preset grades (exploratory; not validated cross-platform equivalents)
GRADE_CUTS = {'1:2': (64, 1024), '1:5': (80, 1280)}  # <= low_max -> low; >= high_min -> high; otherwise mid
SPECIMEN_GROUP = {'Cerebrospinal fluid': 'csf', 'Whole blood': 'blood'}


def p(rel):
    return os.path.join(ROOT, rel)


def read_csv(rel, delimiter=','):
    with open(p(rel), encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def num(v):
    return float(v) if v not in ('', None) else None


def tsv_text(rows, cols):
    out = ['\t'.join(cols)]
    out += ['\t'.join('' if r[c] is None else str(r[c]) for c in cols) for r in rows]
    return '\n'.join(out) + '\n'


def titer_series(t):
    t = int(t)
    if t in SERIES_2:
        return '1:2'
    if t in SERIES_5:
        return '1:5'
    return None


def date_side(d):
    if d < WINDOW[0]:
        return 'pre'
    if d > WINDOW[1]:
        return 'post'
    return 'window'


# CrAg platform is presumed (no laboratory confirmation): from the titre series (1:2 = powers of 2, 1:5 = 5 x powers of 2; a titre of 1 belongs to neither), with the date as auxiliary information and August 2016 as a provisional transition window.
# A series that disagrees with the date outside the window is marked insufficient and gets no grade. The CrAg test date is not recorded, so the admission date of the representative's admission stands in for it.
def crag(adm_id, site, result, titer, test_date):
    """Return platform, status, basis, grade, positive record and log2 titre for one CrAg result."""
    pos = 1 if result == 'positive' else 0
    log2t = round(math.log2(titer), 4) if titer else None
    if result == '':
        return {'platform': '', 'platform_status': '', 'platform_basis': '', 'grade': '', 'pos_record': pos, 'titer_log2': log2t}
    side = date_side(test_date)
    by_date = {'pre': '1:2', 'post': '1:5'}.get(side)
    if (adm_id, site) in CONFIRMED_PLATFORM:
        plat, status, basis = CONFIRMED_PLATFORM[(adm_id, site)], 'confirmed', 'lab_record'
    elif titer:
        s = titer_series(titer)
        if s is None:
            plat, status, basis = (by_date, 'presumed', 'date') if by_date else ('', 'insufficient', 'titre_outside_series_in_window')
        elif side == 'window':
            plat, status, basis = s, 'presumed', 'titre_series_in_window'
        elif s == by_date:
            plat, status, basis = s, 'presumed', 'titre_series_and_date'
        else:
            plat, status, basis = '', 'insufficient', 'series_date_conflict'
    else:
        plat, status, basis = (by_date, 'presumed', 'date') if by_date else ('', 'insufficient', 'qualitative_in_window')
    if result == 'negative':
        grade = 'negative'
    elif not titer:
        grade = 'unknown'
    elif status == 'insufficient':
        grade = 'NA_platform'
    else:
        low_max, high_min = GRADE_CUTS[plat]
        grade = 'low' if titer <= low_max else ('high' if titer >= high_min else 'mid')
    return {'platform': plat, 'platform_status': status, 'platform_basis': basis, 'grade': grade, 'pos_record': pos, 'titer_log2': log2t}


def main():
    global ROOT, OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--root', default=ROOT, help='project root (default: inferred from this script)')
    ap.add_argument('--out-dir', default=None, help='output directory (default: next to this script)')
    a = ap.parse_args()
    ROOT = os.path.abspath(a.root)
    OUT_DIR = os.path.abspath(a.out_dir) if a.out_dir else OUT_DIR

    errors = []

    adm = {r['admission_id']: r for r in read_csv(ADMISSIONS)}
    iso = {r['sample']: r for r in read_csv(ISOLATES)}
    taxa = [line.split('\t')[0] for line in open(p(TAXA), encoding='utf-8').read().splitlines() if line.strip()]
    k_header = open(p(KMAT), encoding='utf-8').readline().rstrip('\n').split('\t')[1:]

    # Representative set and order
    reps = {s for s, r in iso.items() if r['wgs_role'] == 'representative'}
    if set(taxa) != reps:
        errors.append(f'taxa.tsv ({len(taxa)}) differs from the cleaned representatives by {sorted(set(taxa) ^ reps)}')
    if k_header != taxa:
        errors.append('K.tsv column order differs from taxa.tsv')
    if len({iso[s]['admission_id'] for s in taxa if s in iso}) != len(taxa):
        errors.append('representatives do not map to distinct admissions')

    # Errata: applied only if the archived value matches
    errata = read_csv(ERRATA, delimiter='\t')
    applied = collections.defaultdict(dict)
    for e in errata:
        if e['table'] != 'admissions':   # errata of the isolates table are applied by the figure and table scripts
            continue
        if e['key_field'] != 'admission_id':
            errors.append(f'{e["errata_id"]}: only admissions keyed by admission_id are supported')
            continue
        row = adm.get(e['key_value'])
        if row is None or e['field'] not in row:
            errors.append(f'{e["errata_id"]}: admission {e["key_value"]} or field {e["field"]} not found')
            continue
        if row[e['field']] != e['archived_value']:
            errors.append(f'{e["errata_id"]}: archived value {e["archived_value"]!r} does not match the table ({row[e["field"]]!r})')
            continue
        applied[e['key_value']][e['field']] = e
    if errors:
        sys.exit('\n'.join(['FAILED, nothing written:'] + errors))

    rows = []
    for idx, s in enumerate(taxa):
        i = iso[s]
        a = adm[i['admission_id']]
        fixes = applied.get(i['admission_id'], {})
        def val(field):
            return fixes[field]['corrected_value'] if field in fixes else a[field]
        for field, allowed in [('csf_crag_result', {'', 'positive', 'negative'}), ('blood_crag_result', {'', 'positive', 'negative'}),
                               ('alive_10w', {'', '0', '1'}), ('alive_6m', {'', '0', '1'}), ('cryptococcal_meningitis', {'0', '1'}),
                               ('hiv_positive', {'0', '1'}), ('altered_consciousness', {'0', '1'}), ('opening_pressure_censored', {'', '0', '1'}),
                               ('art_at_admission', {'', '0', '1'}),
                               ('discharge_outcome', {'improved', 'died', 'discharged_against_advice', 'transferred'})]:
            v = val(field)
            if v.endswith('.0'):
                v = v[:-2]
            if v not in allowed:
                errors.append(f'{s}: {field} = {v!r} outside {sorted(allowed)}')
        sdate = datetime.date.fromisoformat(i['sampling_date'])
        adate = datetime.date.fromisoformat(a['admission_date'])
        cm = int(a['cryptococcal_meningitis'])
        main = 0 if s in A1_EXCLUDED else 1
        csf3 = all(a[c] != '' for c in ('csf_wbc', 'csf_protein', 'csf_glucose'))
        op = num(a['opening_pressure'])
        cd4 = num(a['cd4_count'])
        r = {'sample': s, 'k_index': idx, 'admission_id': i['admission_id'], 'patient_id': i['patient_id'],
             'a1_main': main, 'a1_excluded_reason': A1_EXCLUDED.get(s, ''), 'cm': cm, 'hiv_positive': int(a['hiv_positive']),
             'age_years': a['age_years'], 'sampling_date': i['sampling_date'], 'sampling_year': sdate.year, 'admission_date': a['admission_date'],
             'specimen_type': i['specimen_type'], 'specimen_group': SPECIMEN_GROUP.get(i['specimen_type'], 'other'),
             'cd4_count': a['cd4_count'], 'cd4_available': int(cd4 is not None), 'art_at_admission': a['art_at_admission'].replace('.0', ''),
             'csf_wbc': a['csf_wbc'], 'csf_protein': a['csf_protein'], 'csf_glucose': a['csf_glucose'],
             'opening_pressure': a['opening_pressure'], 'opening_pressure_censored': a['opening_pressure_censored'].replace('.0', ''),
             'altered_consciousness': a['altered_consciousness'], 'alive_10w': a['alive_10w'].replace('.0', ''), 'alive_6m': a['alive_6m'].replace('.0', ''),
             'discharge_outcome': a['discharge_outcome'],
             'pop_csf': int(main and cm and csf3), 'pop_op': int(main and cm and op is not None),
             'pop_cd4_pheno': int(main and a['hiv_positive'] == '1' and cd4 is not None),
             'crag_date_proxy': a['admission_date']}
        for site in ('csf', 'blood'):
            f_res, f_tit = f'{site}_crag_result', f'{site}_crag_titer'
            res = val(f_res)
            tit = num(a[f_tit])
            if tit is not None and (tit != int(tit) or tit < 1):
                errors.append(f'{s}: {f_tit} = {a[f_tit]!r} is not a positive integer')
                continue
            tit = int(tit) if tit is not None else None
            if res == 'negative' and tit:
                errors.append(f'{s}: {site} CrAg negative with a titre')
            c = crag(i['admission_id'], site, res, tit, adate)
            r.update({f'{site}_crag_result_archived': a[f_res], f'{site}_crag_result': res, f'{site}_crag_errata': fixes[f_res]['errata_id'] if f_res in fixes else '',
                      f'{site}_crag_titer': tit, f'{site}_crag_titer_log2': c['titer_log2'], f'{site}_crag_pos_record': c['pos_record'],
                      f'{site}_crag_platform': c['platform'], f'{site}_crag_platform_status': c['platform_status'],
                      f'{site}_crag_platform_basis': c['platform_basis'], f'{site}_crag_grade': c['grade']})
        rows.append(r)
    if errors:
        sys.exit('\n'.join(['FAILED, nothing written:'] + errors))

    cols = ['sample', 'k_index', 'admission_id', 'patient_id', 'a1_main', 'a1_excluded_reason', 'cm', 'hiv_positive', 'age_years',
            'sampling_date', 'sampling_year', 'admission_date', 'specimen_type', 'specimen_group', 'cd4_count', 'cd4_available',
            'art_at_admission', 'csf_wbc', 'csf_protein', 'csf_glucose', 'pop_csf', 'opening_pressure', 'opening_pressure_censored', 'pop_op',
            'altered_consciousness', 'alive_10w', 'alive_6m', 'discharge_outcome', 'pop_cd4_pheno', 'crag_date_proxy']
    for site in ('csf', 'blood'):
        cols += [f'{site}_crag_{x}' for x in ('result_archived', 'result', 'errata', 'titer', 'titer_log2', 'pos_record', 'platform',
                                              'platform_status', 'platform_basis', 'grade')]

    M = [r for r in rows if r['a1_main']]

    with open(os.path.join(OUT_DIR, 'design.tsv'), 'w', encoding='utf-8', newline='') as fh:
        fh.write(tsv_text(rows, cols))
    print(f'design table: {len(rows)} rows x {len(cols)} columns; A1 main {len(M)}; errata applied: {sum(len(v) for v in applied.values())}')


if __name__ == '__main__':
    main()
