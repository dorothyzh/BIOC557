"""BIOC557 teaching pipeline. Shared by Python, the API and Pyodide.
Samples are rows and features are columns. Missing values stay missing for
Welch tests and MOFA; median imputation is used only for PCA and distances.
"""
from __future__ import annotations
import contextlib
import io
import json
import platform
import warnings
from dataclasses import asdict, dataclass
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist

@dataclass(frozen=True)
class Parameters:
    transform: str = 'auto'
    normalization: str = 'none'
    max_missing: float = 0.30
    model_features: int = 100
    standardize: bool = True
    feature_scaling: str = 'unit_variance'
    mofa_center: bool = True
    mofa_scale_views: bool = False
    factors: int = 3
    iterations: int = 100
    seed: int = 557
    challenge: str = 'none'
    mask_fraction: float = 0.10
    method: str = 'pca'
    feature_alignment: str = 'matched_genes'


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    good = np.flatnonzero(np.isfinite(p))
    if not len(good):
        return q
    order = good[np.argsort(p[good], kind='stable')]
    adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
    q[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return q


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (np.ndarray, list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def validate_dataset(data):
    ids = [str(x) for x in data['sample_ids']]
    if not 6 <= len(ids) <= 300:
        raise ValueError('Use 6–300 samples for this classroom lab.')
    if len(set(ids)) != len(ids) or any(not s.strip() for s in ids):
        raise ValueError('Sample IDs must be nonempty and unique.')
    groups = [str(x) for x in data['groups']]
    if len(groups) != len(ids):
        raise ValueError('The metadata must provide one group for every sample.')
    if data.get('metadata'):
        meta = data['metadata']
        mids = [str(m['sample_id']) for m in meta]
        if len(mids) != len(set(mids)) or set(mids) != set(ids):
            raise ValueError('Metadata sample IDs must match the views exactly, without duplicates.')
        mg = {str(m['sample_id']): str(m.get('group', groups[ids.index(str(m['sample_id']))])) for m in meta}
        if any(mg[s] != g for s, g in zip(ids, groups)):
            raise ValueError('Metadata group labels disagree with the study group vector.')
    levels = list(data.get('group_order', list(dict.fromkeys(groups))))
    if set(levels) != set(groups) or len(levels) != len(set(levels)):
        raise ValueError('Group order must name each observed group exactly once.')
    if len(levels) != 2 or min(groups.count(g) for g in levels) < 3:
        raise ValueError('Use exactly two groups with at least three samples each.')
    if not 2 <= len(data['views']) <= 4:
        raise ValueError('Provide 2–4 molecular views.')
    names = [str(v['name']) for v in data['views']]
    if len(set(names)) != len(names):
        raise ValueError('Each molecular view needs a unique name.')
    views = []
    for v in data['views']:
        vi = [str(x) for x in v['sample_ids']]
        if len(vi) != len(set(vi)) or set(vi) != set(ids):
            raise ValueError(f"{v['name']}: sample IDs must match the metadata exactly, with no duplicates.")
        features = [str(x) for x in v['features']]
        if not 3 <= len(features) <= 5000 or len(set(features)) != len(features):
            raise ValueError(f"{v['name']}: use 3–5000 uniquely named features.")
        x = np.asarray(v['values'], dtype=float)
        if x.shape != (len(vi), len(features)) or np.isinf(x).any():
            raise ValueError(f"{v['name']}: values must be a numeric samples-by-features matrix; blanks represent missing data.")
        positions = {s: i for i, s in enumerate(vi)}
        x = x[[positions[s] for s in ids]].copy()
        views.append({'name': str(v['name']), 'features': features, 'values': x,
                      'unit': str(v.get('unit', 'abundance')),
                      'input_scale': str(v.get('input_scale', 'linear')),
                      'rows_reordered': vi != ids})
    return ids, np.asarray(groups), levels, views


def _pca(x, n_components=3):
    x = x - x.mean(axis=0)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    k = min(n_components, len(s), x.shape[0] - 1)
    scores = u[:, :k] * s[:k]
    # Stable signs improve comparisons; sign itself has no biological meaning.
    for j in range(k):
        idx = np.argmax(np.abs(vt[j]))
        if vt[j, idx] < 0:
            scores[:, j] *= -1
            vt[j] *= -1
    total = np.sum(s ** 2)
    explained = s[:k] ** 2 / total if total > 0 else np.zeros(k)
    return {'scores': scores, 'explained': explained, 'loadings': vt[:k].T}


def _differential(x, groups, levels, features):
    a = x[groups == levels[0]]
    b = x[groups == levels[1]]
    n_a = np.sum(np.isfinite(a), axis=0)
    n_b = np.sum(np.isfinite(b), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tests = stats.ttest_ind(b, a, axis=0, equal_var=False, nan_policy='omit')
        delta = np.nanmean(b, axis=0) - np.nanmean(a, axis=0)
    p = np.asarray(tests.pvalue)
    p[(n_a < 3) | (n_b < 3)] = np.nan
    q = bh_adjust(p)
    rows = [dict(feature=f, mean_difference=delta[j], p_value=p[j], q_value=q[j],
                 observed_group_1=n_a[j], observed_group_2=n_b[j])
            for j, f in enumerate(features)]
    return rows


def _fit_mofa(matrices, feature_names, ids, names, params):
    from mofapy2.run.entry_point import entry_point
    with contextlib.redirect_stdout(io.StringIO()):
        entry = entry_point()
        entry.set_data_options(scale_views=params.mofa_scale_views, scale_groups=False, center_groups=params.mofa_center)
        entry.set_data_matrix([[a] for a in matrices], likelihoods=['gaussian'] * len(matrices),
                              views_names=names, groups_names=['cohort'],
                              samples_names=[ids], features_names=[[name + '::' + f for f in fs] for name, fs in zip(names, feature_names)])
        # build() replaces NaNs by zeros internally; preserve the observed-entry mask first.
        inputs = [np.asarray(x).copy() for x in entry.data]
        entry.set_model_options(factors=params.factors, spikeslab_weights=True, ard_weights=True)
        entry.set_train_options(iter=params.iterations, convergence_mode='fast',
                                startELBO=1, freqELBO=5, dropR2=None,
                                gpu_mode=False, verbose=False, quiet=True, seed=params.seed)
        entry.build()
        entry.run()
    ex = entry.model.getExpectations()
    scores = ex['Z']['E'].copy()
    weights = [v['E'].copy() for v in ex['W']]
    for j in range(scores.shape[1]):
        merged = np.concatenate([w[:, j] for w in weights])
        if merged[np.argmax(np.abs(merged))] < 0:
            scores[:, j] *= -1
            for w in weights:
                w[:, j] *= -1
    # Use the actual post-centering/post-view-scaling matrices used by MOFA.
    r2, total_r2 = [], []
    for x, w in zip(inputs, weights):
        observed = np.isfinite(x)
        denom = np.nansum(x ** 2)
        def reconstruction(prediction):
            return 1 - np.nansum(np.where(observed, (x - prediction) ** 2, np.nan)) / denom
        r2.append([reconstruction(np.outer(scores[:, j], w[:, j])) for j in range(scores.shape[1])])
        total_r2.append(reconstruction(scores @ w.T))
    # Present the strongest fitted sources of variation first, using no group labels.
    order = np.argsort(-np.mean(r2, axis=0), kind='stable')
    return {'scores': scores[:, order], 'weights': [w[:, order] for w in weights], 'variance_explained': np.asarray(r2)[:, order],
            'total_r2': total_r2, 'model_input': inputs,
            'factors': scores.shape[1], 'model': 'MOFA+ / mofapy2 0.7.5',
            'iterations_requested': params.iterations,
            'center': params.mofa_center, 'scale_views': params.mofa_scale_views,
            'baseline': 'centered variation' if params.mofa_center else 'uncentered input energy, including mean levels',
            'note': 'Factor-wise reconstruction uses actual MOFA input and a zero reference; contributions overlap. Total reconstruction is computed from all factors together, not by adding factor-wise scores.'}


def gene_pairs(data):
    """Only explicit, one-to-one identifiers; never match features by position."""
    views = {v['name']: v for v in data['views']}
    if not {'Proteomics', 'RNA'} <= views.keys():
        return []
    protein, rna = (set(views[n]['features']) for n in ['Proteomics', 'RNA'])
    pairs = data.get('shared_gene_pairs')
    if pairs is None:
        pairs = [dict(protein_feature_id=f, rna_feature_id=f, gene_symbol=f, gene_id=f)
                 for f in sorted(protein & rna)]
    seen_p, seen_r, seen_g, valid = set(), set(), set(), []
    for pair in pairs:
        p, r = pair['protein_feature_id'], pair['rna_feature_id']
        g = pair.get('gene_id', p if p == r else pair['gene_symbol'])
        if p not in protein or r not in rna:
            raise ValueError('Gene crosswalk references a feature absent from its view.')
        if p in seen_p or r in seen_r or g in seen_g:
            raise ValueError('Gene crosswalk must be one-to-one; resolve duplicate or ambiguous mappings first.')
        seen_p.add(p); seen_r.add(r); seen_g.add(g)
        valid.append(dict(pair, gene_id=g))
    return valid


def analyze(data, parameters=None):
    options = dict(parameters or {})
    if 'feature_scaling' not in options and options.get('standardize') is False:
        options['feature_scaling'] = 'none'
    params = Parameters(**options)
    if params.feature_alignment not in ['matched_genes', 'independent']:
        raise ValueError('Choose matched genes or independent feature selection.')
    if params.feature_scaling not in ['none', 'unit_variance', 'pareto']:
        raise ValueError('Choose no feature scaling, unit variance or Pareto scaling.')
    if not isinstance(params.mofa_center, bool) or not isinstance(params.mofa_scale_views, bool):
        raise ValueError('MOFA center and view scaling must be true or false.')
    if params.transform not in ['auto', 'log2', 'none'] or params.normalization not in ['none', 'sample_median']:
        raise ValueError('Choose a supported transform and normalization.')
    if not 0 <= params.max_missing <= 0.8 or not 10 <= params.model_features <= 500:
        raise ValueError('Set maximum missingness to 0–0.8 and model features to 10–500 per view.')
    if not 2 <= params.factors <= 6 or not 25 <= params.iterations <= 300:
        raise ValueError('Choose 2–6 factors and 25–300 MOFA iterations.')
    if params.method not in ['pca', 'mofa']:
        raise ValueError('Choose PCA or MOFA.')
    if params.challenge not in ['none', 'shuffle_rows', 'break_pairing', 'missing']:
        raise ValueError('Unknown challenge.')
    if not 0 <= params.mask_fraction <= 0.5:
        raise ValueError('The missing-value challenge can mask at most 50%.')
    ids, groups, levels, source_views = validate_dataset(data)
    rng = np.random.default_rng(params.seed)
    result = {'dataset': data.get('title', 'Uploaded teaching data'), 'sample_ids': ids,
              'annotation_scope': 'bee' if data.get('teaching', {}).get('annotation_mode') == 'bee_genes' else 'human',
              'groups': groups.tolist(), 'group_order': levels, 'parameters': asdict(params),
              'views': [], 'challenge': {'kind': params.challenge, 'synthetic': params.challenge in ['break_pairing', 'missing']}}
    if params.challenge == 'break_pairing':
        perm = rng.permutation(len(ids))
        source_views[-1]['values'] = source_views[-1]['values'][perm]
        result['challenge']['description'] = 'Teaching perturbation: values in the final view were deliberately assigned to the wrong sample IDs.'
    elif params.challenge == 'missing':
        for v in source_views:
            mask = rng.random(v['values'].shape) < params.mask_fraction
            v['values'][mask] = np.nan
        result['challenge']['description'] = 'Teaching perturbation: extra values were randomly hidden, independently of group.'
    elif params.challenge == 'shuffle_rows':
        # Re-run the actual ID alignment after a harmless row-order permutation.
        permuted = dict(data, views=[])
        for v in data['views']:
            perm = rng.permutation(len(v['sample_ids']))
            columns = rng.permutation(len(v['features']))
            permuted['views'].append(dict(v, sample_ids=[v['sample_ids'][i] for i in perm],
                                         features=[v['features'][j] for j in columns],
                                         values=[[v['values'][i][j] for j in columns] for i in perm]))
        _, _, _, source_views = validate_dataset(permuted)
        result['challenge']['description'] = 'Rows and feature columns were shuffled with their IDs intact. Alignment by sample ID and gene ID restores the pairing.'
    matrices, model_features, retained_views, prepared_views = [], [], [], []
    pairs = gene_pairs(data)
    for v in source_views:
        original = v['values']
        x = original.copy()
        finite = np.isfinite(x)
        apply_log = params.transform == 'log2' or (params.transform == 'auto' and v['input_scale'] == 'linear')
        if apply_log:
            if v['input_scale'] != 'linear':
                raise ValueError(f"{v['name']} is already transformed. Choose No log transform.")
            if np.any(x[finite] <= 0):
                raise ValueError(f"{v['name']} has nonpositive values. Supply positive linear abundances, or use No log transform for logged data.")
            x[finite] = np.log2(x[finite])
        if params.normalization == 'sample_median':
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                med = np.nanmedian(x, axis=1)
            if not np.isfinite(med).all():
                raise ValueError(f"{v['name']}: a sample has no observed values.")
            x = x - med[:, None] + np.nanmedian(med)
        missing = np.mean(~np.isfinite(x), axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            variance = np.nanvar(x, axis=0, ddof=1)
        keep = (missing <= params.max_missing) & np.isfinite(variance) & (variance > 1e-12)
        if np.sum(keep) < max(3, params.factors + 1):
            raise ValueError(f"{v['name']}: too few variable features remain. Relax the missingness filter.")
        x = x[:, keep]
        features = np.asarray(v['features'])[keep].tolist()
        retained_views.append({'name': v['name'], 'features': features, 'values': x.copy()})
        retained_var = variance[keep]
        prepared_views.append(dict(view=v, original=original, x=x, features=features,
                                   variance=retained_var, apply_log=apply_log))
    paired_selection, aligned_pairs = {}, []
    if params.feature_alignment == 'matched_genes' and {'Proteomics', 'RNA'} <= {v['name'] for v in source_views}:
        by_name = {p['view']['name']: p for p in prepared_views}
        indices = {n: {f: j for j, f in enumerate(by_name[n]['features'])} for n in ['Proteomics', 'RNA']}
        # Percentile ranks give equal weight to both assays despite different units.
        ranks = {n: stats.rankdata(by_name[n]['variance'], method='average') / len(by_name[n]['variance'])
                 for n in ['Proteomics', 'RNA']}
        eligible = [p for p in pairs if p['protein_feature_id'] in indices['Proteomics'] and p['rna_feature_id'] in indices['RNA']]
        eligible.sort(key=lambda p: (-sum(ranks[n][indices[n][p[k]]] for n, k in [('Proteomics', 'protein_feature_id'), ('RNA', 'rna_feature_id')]), p['gene_id']))
        aligned_pairs = eligible[:params.model_features]
        if len(aligned_pairs) < max(4, params.factors + 1):
            raise ValueError('Too few verified shared genes pass both views’ QC. Use independent selection or provide a gene crosswalk.')
        for n, k in [('Proteomics', 'protein_feature_id'), ('RNA', 'rna_feature_id')]:
            paired_selection[n] = [indices[n][p[k]] for p in aligned_pairs]
        result['gene_alignment'] = dict(mode='matched_genes', verified_input_pairs=len(pairs),
            eligible_pairs=len(eligible), model_pairs=aligned_pairs,
            note='The same genes, in the same order, in protein and RNA model matrices. Sample rows align independently by sample ID. Other views retain their own features. QC and differential tests retain all eligible source features.')
    else:
        result['gene_alignment'] = dict(mode='independent', verified_input_pairs=len(pairs), model_pairs=[],
            note='Features selected within each view; verified same-gene pairs are retained for cross-layer comparisons.')
    for prepared in prepared_views:
        v, original, x, features, retained_var, apply_log = (prepared[k] for k in ['view', 'original', 'x', 'features', 'variance', 'apply_log'])
        select = paired_selection.get(v['name'], np.argsort(-retained_var, kind='stable')[:params.model_features])
        model_x = x[:, select].copy()
        feat = [features[i] for i in select]
        # Keep centering separate from SD division so the no-center option retains offsets.
        if params.mofa_center:
            model_x = model_x - np.nanmean(model_x, axis=0)
        sd = np.nanstd(model_x, axis=0, ddof=1)
        if params.feature_scaling == 'unit_variance':
            model_x = model_x / sd
        elif params.feature_scaling == 'pareto':
            model_x = model_x / np.sqrt(sd)
        # PCA needs complete data. MOFA receives model_x with its original missing values.
        imputed = model_x.copy()
        medians = np.nanmedian(imputed, axis=0)
        miss = np.where(~np.isfinite(imputed))
        imputed[miss] = medians[miss[1]]
        pca = _pca(imputed)
        matrices.append(model_x)
        model_features.append(feat)
        differential = _differential(x, groups, levels, features)
        result['views'].append({'name': v['name'], 'unit': v['unit'],
            'input_features': original.shape[1], 'input_scale': v['input_scale'],
            'analysis_scale': 'log2 abundance' if apply_log else v['unit'],
            'transformation_applied': 'log2' if apply_log else 'none', 'retained_features': len(features),
            'selected_features': feat, 'retained_feature_ids': features, 'model_feature_count': len(feat),
            'missing_per_sample': np.mean(~np.isfinite(original), axis=1),
            'overall_missing': np.mean(~np.isfinite(original)),
            'missing_after_filter': np.mean(~np.isfinite(x)),
            'imputed_for_pca': len(miss[0]), 'rows_reordered': v['rows_reordered'],
            'pca': pca, 'differential': differential,
            'processed_matrix': model_x, 'pca_matrix': imputed,
            'distribution_raw': [row[np.isfinite(row)].tolist() for row in original],
            'distribution_processed': [row[np.isfinite(row)].tolist() for row in x]})
    # Equal block weighting prevents the larger view from dominating solely by feature count.
    joint = np.concatenate([v['pca_matrix'] / np.sqrt(v['model_feature_count']) for v in result['views']], axis=1)
    result['joint_pca'] = _pca(joint)
    result['joint_pca']['weighting'] = 'Feature-standardized when selected, then each view divided by sqrt(number of selected features).'
    distances = [pdist(v['pca_matrix']) for v in result['views']]
    agreement = np.eye(len(distances))
    for i in range(len(distances)):
        for j in range(i + 1, len(distances)):
            agreement[i, j] = agreement[j, i] = stats.spearmanr(distances[i], distances[j]).statistic
    result['cross_view_distance_spearman'] = agreement[0, 1]
    result['cross_view_distance_matrix'] = agreement
    result['mean_cross_view_distance_spearman'] = np.nanmean(agreement[np.triu_indices(len(distances), 1)])
    result['cross_view_metric_note'] = 'Descriptive correlation between within-view sample distances. Pairwise distances are dependent; no inferential p value is reported.'
    # Matched gene examples use retained, prepared values before model z-scoring.
    examples = []
    retained_by_name = {v['name']: v for v in retained_views}
    if 'Proteomics' in retained_by_name and 'RNA' in retained_by_name:
        protein, rna = retained_by_name['Proteomics'], retained_by_name['RNA']
        pi, ri = ({f: j for j, f in enumerate(v['features'])} for v in [protein, rna])
        for pair in pairs:
            pf, rf = pair['protein_feature_id'], pair['rna_feature_id']
            if pf not in pi or rf not in ri:
                continue
            y, x = protein['values'][:, pi[pf]], rna['values'][:, ri[rf]]
            observed = np.isfinite(x) & np.isfinite(y)
            if observed.sum() < 6 or np.nanvar(x) <= 0 or np.nanvar(y) <= 0:
                continue
            examples.append({'gene': pair['gene_symbol'], 'gene_id': pair['gene_id'],
                'protein_feature_id': pf, 'rna_feature_id': rf, 'sample_ids': ids,
                'rna': x, 'protein': y,
                'observed_n': int(observed.sum()), 'spearman': stats.spearmanr(x[observed], y[observed]).statistic,
                'selection_variance_product': np.nanvar(x, ddof=1) * np.nanvar(y, ddof=1)})
        examples.sort(key=lambda g: -g['selection_variance_product'])
    teaching_markers = [] if result['annotation_scope']=='bee' else ['EPCAM', 'VIM', 'EGFR', 'ERBB2', 'TYR', 'AXL']
    examples.sort(key=lambda g: (teaching_markers.index(g['gene']) if g['gene'] in teaching_markers else len(teaching_markers), -g['selection_variance_product']))
    result['shared_gene_examples'] = examples[:40]
    result['shared_gene_example_note'] = 'Matched gene symbols present in both retained views. Preselected teaching markers appear first, followed by genes ranked by the product of within-view variances without group labels. Values precede model z-scoring. Correlations pool the recorded study groups.'
    if params.method == 'mofa':
        result['mofa'] = _fit_mofa(matrices, model_features, ids, [v['name'] for v in source_views], params)
    result['effect_label'] = ('Difference in mean log₂ abundance' if params.transform == 'log2' else 'Difference in mean input-scale abundance')
    result['contrast'] = f'{levels[1]} − {levels[0]}'
    result['test'] = 'Unadjusted two-group Welch tests on observed values, BH FDR within each view; at least 3 observed samples per group. No imputed values enter the tests.'
    result['versions'] = {'pipeline': '3.0.0', 'python': platform.python_version(), 'numpy': np.__version__}
    return clean_json(result)


def run_json(dataset_json, parameters_json='{}'):
    return json.dumps(analyze(json.loads(dataset_json), json.loads(parameters_json)), allow_nan=False)
