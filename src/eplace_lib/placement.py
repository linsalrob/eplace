"""
Placement planning for ePLACE query sequences.

This module separates two concepts that were previously entangled:

1. BLAST/MMseqs confidence for direct taxonomic classification.
2. Whether a query sequence has enough evidence to enter phylogenetic placement.

The core design principle is that low pairwise search identity should reduce
classification confidence, not silently remove a query from tree-building.

Typical routes:

- classification_hit: retained/high-confidence hits are available.
- placement_hit: no retained hits, but moderate raw hits can anchor placement.
- rescue_hit: only weak raw hits are available, useful for broad context.
- backbone_only: no usable search hit, but a backbone route may still place it.
- no_evidence: no usable hit and no backbone route was requested/available.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

from .blast_analysis import BlastHit


CLASSIFICATION_HIT = "classification_hit"
PLACEMENT_HIT = "placement_hit"
RESCUE_HIT = "rescue_hit"
BACKBONE_ONLY = "backbone_only"
NO_EVIDENCE = "no_evidence"

PLACEMENT_ROUTES = {
    CLASSIFICATION_HIT,
    PLACEMENT_HIT,
    RESCUE_HIT,
    BACKBONE_ONLY,
    NO_EVIDENCE,
}

DEFAULT_PLACEMENT_MIN_IDENTITY = 70.0
DEFAULT_PLACEMENT_MIN_COVERAGE = 50.0
DEFAULT_RESCUE_MIN_IDENTITY = 50.0
DEFAULT_RESCUE_MIN_COVERAGE = 40.0


@dataclass(frozen=True)
class PlacementThresholds:
    """
    Thresholds used to route queries through the placement workflow.

    classification_min_identity / classification_min_coverage are normally the
    existing ePLACE retained-hit thresholds (currently 90/80 by CLI default).
    placement_* and rescue_* are lower thresholds used only to select reference
    context for phylogenetic placement, not to claim confident taxonomy.
    """

    classification_min_identity: float = 90.0
    classification_min_coverage: float = 80.0
    placement_min_identity: float = DEFAULT_PLACEMENT_MIN_IDENTITY
    placement_min_coverage: float = DEFAULT_PLACEMENT_MIN_COVERAGE
    rescue_min_identity: float = DEFAULT_RESCUE_MIN_IDENTITY
    rescue_min_coverage: float = DEFAULT_RESCUE_MIN_COVERAGE


@dataclass
class QueryPlacementPlan:
    """
    Placement decision for a single query sequence.

    Attributes:
        query_id: Query/ZOTU/ASV identifier from the input FASTA.
        route: One of PLACEMENT_ROUTES.
        retained_hits: Hits that passed the strict classification filter.
        raw_hits: All search hits parsed from raw BLAST/MMseqs output.
        tree_hits: Hits selected for phylogenetic tree construction.
        best_raw_hit: Best raw hit by bit score, identity, then coverage.
        best_tree_hit: Best hit used for tree construction.
        reason: Human-readable explanation of the route.
    """

    query_id: str
    route: str
    retained_hits: List[BlastHit] = field(default_factory=list)
    raw_hits: List[BlastHit] = field(default_factory=list)
    tree_hits: List[BlastHit] = field(default_factory=list)
    best_raw_hit: Optional[BlastHit] = None
    best_tree_hit: Optional[BlastHit] = None
    reason: str = ""

    @property
    def phylogenetic_placement_attempted(self) -> str:
        """Return Yes/No/Backbone pending status for TSV-style output."""
        if self.route in {CLASSIFICATION_HIT, PLACEMENT_HIT, RESCUE_HIT}:
            return "Yes"
        if self.route == BACKBONE_ONLY:
            return "Backbone route required"
        return "No"

    def to_summary_fields(self) -> Dict[str, object]:
        """
        Return flat fields suitable for classification.tsv output.

        The values are intentionally string/simple numeric friendly so callers
        can merge them directly into the existing classification dict.
        """
        return {
            "placement_route": self.route,
            "placement_hit_count": len(self.tree_hits),
            "placement_best_percent_identity": _format_hit_float(
                self.best_tree_hit, "percent_identity"
            ),
            "placement_best_query_coverage": _format_hit_float(
                self.best_tree_hit, "query_coverage"
            ),
            "placement_best_bit_score": _format_hit_float(
                self.best_tree_hit, "bit_score"
            ),
            "placement_best_subject_id": (
                self.best_tree_hit.subject_id if self.best_tree_hit else "N/A"
            ),
            "placement_best_taxon_name": _best_taxon_name(self.best_tree_hit),
            "phylogenetic_placement_attempted": self.phylogenetic_placement_attempted,
            "phylogenetic_placement_reason": self.reason,
        }


def _format_hit_float(hit: Optional[BlastHit], attr: str, digits: int = 3) -> str:
    """Format a numeric BlastHit attribute while preserving N/A."""
    if hit is None:
        return "N/A"
    try:
        return f"{float(getattr(hit, attr)):.{digits}f}"
    except Exception:
        return "N/A"

def _best_taxon_name(hit: Optional[BlastHit]) -> str:
    """
    Return the most informative readable taxon name for a placement hit.

    Preference order:
    species > genus > family > order > class > phylum > domain > subject_id
    """
    if hit is None:
        return "N/A"

    if getattr(hit, "subject_taxonomy", None):
        for rank in ["species", "genus", "family", "order", "class", "phylum", "domain"]:
            if rank in hit.subject_taxonomy:
                _taxid, name = hit.subject_taxonomy[rank]
                if name:
                    return name

    return getattr(hit, "subject_id", "N/A")
    
def _group_hits_by_query(hits: Iterable[BlastHit]) -> Dict[str, List[BlastHit]]:
    """Group BlastHit objects by query_id."""
    grouped: Dict[str, List[BlastHit]] = defaultdict(list)
    for hit in hits:
        grouped[hit.query_id].append(hit)
    return dict(grouped)


def _sort_hits_best_first(hits: Sequence[BlastHit]) -> List[BlastHit]:
    """Sort hits by bit score, percent identity, coverage, then alignment length."""
    return sorted(
        hits,
        key=lambda h: (
            float(getattr(h, "bit_score", 0.0) or 0.0),
            float(getattr(h, "percent_identity", 0.0) or 0.0),
            float(getattr(h, "query_coverage", 0.0) or 0.0),
            int(getattr(h, "alignment_length", 0) or 0),
        ),
        reverse=True,
    )


def _passes_thresholds(
    hit: BlastHit,
    min_identity: float,
    min_coverage: float,
) -> bool:
    """Return True if a hit passes identity and query coverage thresholds."""
    try:
        identity = float(hit.percent_identity)
        coverage = float(hit.query_coverage)
    except Exception:
        return False

    return identity >= min_identity and coverage >= min_coverage


def _select_hits_passing(
    hits: Sequence[BlastHit],
    min_identity: float,
    min_coverage: float,
    max_hits: Optional[int] = None,
) -> List[BlastHit]:
    """Return best-first hits passing the supplied thresholds."""
    selected = [
        hit
        for hit in _sort_hits_best_first(hits)
        if _passes_thresholds(hit, min_identity, min_coverage)
    ]
    if max_hits is not None:
        return selected[:max_hits]
    return selected


def _best_hit(hits: Sequence[BlastHit]) -> Optional[BlastHit]:
    """Return the best hit using the same sorting convention as placement."""
    sorted_hits = _sort_hits_best_first(hits)
    return sorted_hits[0] if sorted_hits else None


def _deduplicate_hits_by_subject(hits: Sequence[BlastHit]) -> List[BlastHit]:
    """
    Keep the best hit per subject_id/accession-like identifier.

    This avoids adding identical reference sequences multiple times when raw hit
    rescue uses a broad candidate set.
    """
    best_by_subject: Dict[str, BlastHit] = {}

    for hit in _sort_hits_best_first(hits):
        key = hit.get_accession() if hasattr(hit, "get_accession") else hit.subject_id
        if key not in best_by_subject:
            best_by_subject[key] = hit

    return list(best_by_subject.values())


def build_query_placement_plan(
    sequences: Mapping[str, str],
    retained_hits: Sequence[BlastHit],
    raw_hits: Sequence[BlastHit],
    thresholds: Optional[PlacementThresholds] = None,
    *,
    max_placement_hits_per_query: Optional[int] = None,
    max_rescue_hits_per_query: Optional[int] = 25,
    backbone_available: bool = False,
) -> Dict[str, QueryPlacementPlan]:
    """
    Build placement routes for every input query sequence.

    Args:
        sequences: Mapping of all query IDs to sequence strings, usually from
            FastaReader.read_fasta(). Every key receives a placement plan.
        retained_hits: Hits that passed the strict ePLACE search filter.
        raw_hits: Raw/unfiltered hits parsed from search output.
        thresholds: Optional PlacementThresholds object. Defaults are used when
            omitted.
        max_placement_hits_per_query: Optional cap on moderate placement hits.
            None means keep all passing placement hits.
        max_rescue_hits_per_query: Optional cap on weak rescue hits. The default
            caps rescue hits at 25 to prevent very noisy trees.
        backbone_available: Whether the caller has a backbone alignment/tree
            route for queries with no usable search hit.

    Returns:
        Dictionary mapping query_id to QueryPlacementPlan.
    """
    thresholds = thresholds or PlacementThresholds()

    retained_by_query = _group_hits_by_query(retained_hits)
    raw_by_query = _group_hits_by_query(raw_hits)

    plan: Dict[str, QueryPlacementPlan] = {}

    for query_id in sequences.keys():
        query_retained_hits = _sort_hits_best_first(retained_by_query.get(query_id, []))
        query_raw_hits = _sort_hits_best_first(raw_by_query.get(query_id, []))

        # Route 1: strict retained/classification hits.
        if query_retained_hits:
            tree_hits = _deduplicate_hits_by_subject(query_retained_hits)
            plan[query_id] = QueryPlacementPlan(
                query_id=query_id,
                route=CLASSIFICATION_HIT,
                retained_hits=query_retained_hits,
                raw_hits=query_raw_hits,
                tree_hits=tree_hits,
                best_raw_hit=_best_hit(query_raw_hits),
                best_tree_hit=_best_hit(tree_hits),
                reason=(
                    "One or more hits passed the strict retained-hit thresholds "
                    f"({thresholds.classification_min_identity:.1f}% identity, "
                    f"{thresholds.classification_min_coverage:.1f}% coverage). "
                    "Use these hits for confident classification and local phylogenetic placement."
                ),
            )
            continue

        # Route 2: moderate raw hits that can anchor phylogenetic placement.
        placement_hits = _select_hits_passing(
            query_raw_hits,
            thresholds.placement_min_identity,
            thresholds.placement_min_coverage,
            max_hits=max_placement_hits_per_query,
        )
        placement_hits = _deduplicate_hits_by_subject(placement_hits)

        if placement_hits:
            plan[query_id] = QueryPlacementPlan(
                query_id=query_id,
                route=PLACEMENT_HIT,
                retained_hits=query_retained_hits,
                raw_hits=query_raw_hits,
                tree_hits=placement_hits,
                best_raw_hit=_best_hit(query_raw_hits),
                best_tree_hit=_best_hit(placement_hits),
                reason=(
                    "No hits passed the strict retained-hit thresholds, but one or more raw hits "
                    f"passed the placement thresholds ({thresholds.placement_min_identity:.1f}% identity, "
                    f"{thresholds.placement_min_coverage:.1f}% coverage). Use these hits to select "
                    "reference context for phylogenetic placement, but report taxonomy cautiously."
                ),
            )
            continue

        # Route 3: weak raw hits used only for broad rescue placement.
        rescue_hits = _select_hits_passing(
            query_raw_hits,
            thresholds.rescue_min_identity,
            thresholds.rescue_min_coverage,
            max_hits=max_rescue_hits_per_query,
        )
        rescue_hits = _deduplicate_hits_by_subject(rescue_hits)

        if rescue_hits:
            plan[query_id] = QueryPlacementPlan(
                query_id=query_id,
                route=RESCUE_HIT,
                retained_hits=query_retained_hits,
                raw_hits=query_raw_hits,
                tree_hits=rescue_hits,
                best_raw_hit=_best_hit(query_raw_hits),
                best_tree_hit=_best_hit(rescue_hits),
                reason=(
                    "No retained or moderate placement hits were found, but weak raw hits passed "
                    f"the rescue thresholds ({thresholds.rescue_min_identity:.1f}% identity, "
                    f"{thresholds.rescue_min_coverage:.1f}% coverage). Use only for broad, "
                    "low-confidence phylogenetic context."
                ),
            )
            continue

        # Route 4/5: no usable search evidence. A backbone can still carry the query.
        if backbone_available:
            route = BACKBONE_ONLY
            reason = (
                "No raw search hits passed retained, placement, or rescue thresholds. "
                "A backbone placement route is available, so this query should be appended "
                "to the backbone alignment/tree as unresolved database evidence."
            )
        else:
            route = NO_EVIDENCE
            reason = (
                "No raw search hits passed retained, placement, or rescue thresholds, "
                "and no backbone placement route was available."
            )

        plan[query_id] = QueryPlacementPlan(
            query_id=query_id,
            route=route,
            retained_hits=query_retained_hits,
            raw_hits=query_raw_hits,
            tree_hits=[],
            best_raw_hit=_best_hit(query_raw_hits),
            best_tree_hit=None,
            reason=reason,
        )

    return plan


def collect_tree_candidate_hits(
    placement_plan: Mapping[str, QueryPlacementPlan],
    routes: Optional[Iterable[str]] = None,
) -> List[BlastHit]:
    """
    Collect hits that should be used for representative extraction/tree building.

    Args:
        placement_plan: Output from build_query_placement_plan().
        routes: Optional subset of routes to include. By default, includes
            classification_hit, placement_hit, and rescue_hit.

    Returns:
        De-duplicated list of BlastHit objects selected for tree construction.
    """
    include_routes = set(routes or {CLASSIFICATION_HIT, PLACEMENT_HIT, RESCUE_HIT})
    hits: List[BlastHit] = []

    for query_plan in placement_plan.values():
        if query_plan.route in include_routes:
            hits.extend(query_plan.tree_hits)

    return hits


def group_tree_hits_by_query(
    placement_plan: Mapping[str, QueryPlacementPlan],
    routes: Optional[Iterable[str]] = None,
) -> Dict[str, List[BlastHit]]:
    """
    Return query_id -> tree_hits for routes that should enter tree construction.

    This is a drop-in replacement for building hits_by_query_map from
    filtered_hits in cli.py.
    """
    include_routes = set(routes or {CLASSIFICATION_HIT, PLACEMENT_HIT, RESCUE_HIT})

    return {
        query_id: query_plan.tree_hits
        for query_id, query_plan in placement_plan.items()
        if query_plan.route in include_routes and query_plan.tree_hits
    }


def placement_plan_summary_counts(
    placement_plan: Mapping[str, QueryPlacementPlan],
) -> Dict[str, int]:
    """Count queries assigned to each placement route."""
    counts = {route: 0 for route in sorted(PLACEMENT_ROUTES)}
    for query_plan in placement_plan.values():
        counts[query_plan.route] = counts.get(query_plan.route, 0) + 1
    return counts
