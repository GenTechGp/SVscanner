"""
RM element selection methods for filter_rm().

Three strategies for choosing a representative subset of RepeatMasker hits
within a repeat class when multiple partially-overlapping hits cover the SV.

select_strict  — current behaviour: sort by sv_coverage desc, strict no-overlap
select_greedy  — accept if unique contribution / sv_length >= min_sv_coverage
select_dp      — minimum interval cover DP (optimal: max coverage, min count)
"""


def select_strict(elements):
    """Sort by sv_coverage descending; accept only non-overlapping elements."""
    selected = []
    for e in sorted(elements, key=lambda e: e['sv_coverage'], reverse=True):
        if all(not _overlaps(e, s) for s in selected):
            selected.append(e)
    return selected


def select_greedy(elements, sv_length, min_sv_coverage):
    """
    Sort by sv_coverage descending; accept an element if its unique contribution
    (bp not already covered by selected elements) / sv_length >= min_sv_coverage.

    Selected elements may overlap each other.
    Order-dependent when sv_coverage values are tied.
    """
    selected = []
    covered = []  # (start, end) intervals already covered

    for e in sorted(elements, key=lambda e: e['sv_coverage'], reverse=True):
        if _unique_bp(e['repeat_start'], e['repeat_end'], covered) / sv_length >= min_sv_coverage:
            selected.append(e)
            covered.append((e['repeat_start'], e['repeat_end']))

    return selected


def select_dp(elements):
    """
    Minimum interval cover DP.

    Finds the minimum-cardinality subset with maximum union coverage, globally
    optimal regardless of element ordering.  Tie-break: maximum sum of
    sv_coverage weights.

    Phase 1 — decompose into connected components (gaps separate independent
               subproblems).
    Phase 2 — DP per component, O(n^2); adequate since n per SV is small.
    """
    if not elements:
        return []

    # Phase 1: connected components (sort by start, merge overlapping/touching)
    comp, components, max_end = [], [], -1
    for e in sorted(elements, key=lambda e: e['repeat_start']):
        if not comp or e['repeat_start'] <= max_end:
            comp.append(e)
            max_end = max(max_end, e['repeat_end'])
        else:
            components.append(comp)
            comp = [e]
            max_end = e['repeat_end']
    if comp:
        components.append(comp)

    # Phase 2: DP per component
    selected = []
    for comp in components:
        selected.extend(_dp_component(comp))
    return selected


def _dp_component(elements):
    """DP for one connected component; returns the selected element subset."""
    by_right = sorted(elements, key=lambda e: e['repeat_end'])
    a = min(e['repeat_start'] for e in elements)
    b = max(e['repeat_end'] for e in elements)

    INF = float('inf')
    # dp[x] = (min_count, max_weight, prev_x, index_into_by_right)
    dp = {a: (0, 0.0, None, None)}

    for i, e in enumerate(by_right):
        l, r, w = e['repeat_start'], e['repeat_end'], e['sv_coverage']

        best_x, best_count, best_weight = None, INF, -INF
        for x, (count, weight, _, _) in dp.items():
            if l <= x < r:
                if count < best_count or (count == best_count and weight > best_weight):
                    best_x, best_count, best_weight = x, count, weight

        if best_x is None:
            continue

        new_count, new_weight = best_count + 1, best_weight + w
        if r not in dp or new_count < dp[r][0] or (new_count == dp[r][0] and new_weight > dp[r][1]):
            dp[r] = (new_count, new_weight, best_x, i)

    if b not in dp or dp[b][3] is None:
        return elements  # fallback: return all elements

    path, x = [], b
    while dp[x][3] is not None:
        path.append(by_right[dp[x][3]])
        x = dp[x][2]
    path.reverse()
    return path


# ── helpers ────────────────────────────────────────────────────────────────────

def _overlaps(e1, e2):
    return not (e1['repeat_end'] <= e2['repeat_start'] or e2['repeat_end'] <= e1['repeat_start'])


def _unique_bp(start, end, covered):
    """Count bp in [start, end) not covered by any interval in covered."""
    if not covered:
        return end - start
    clipped = [(max(s, start), min(e, end)) for s, e in covered if s < end and e > start]
    if not clipped:
        return end - start
    clipped.sort()
    cs, ce = clipped[0]
    union = 0
    for s, e in clipped[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            union += ce - cs
            cs, ce = s, e
    union += ce - cs
    return (end - start) - union
