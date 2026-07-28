"""Importable toolkit for pulling NCBI Protein (eutils) records as a historical
corpus source.

Built incrementally, one step at a time, per FINDINGS_NCBI.md. Every function
here is meant to be called from a notebook or another script, e.g.:

    from scripts.ncbi.pull_data import list_searchable_fields, search_count

There is no `if __name__ == "__main__":`-only logic and no hidden global
state (beyond the rate-limiter's internal clock, which callers never touch
directly) -- import and call.

Rate limiting: NCBI eutils allows ~3 requests/sec unauthenticated. All eutils
calls in this module go through `_rate_limited_get`, so every function here
automatically respects that limit -- callers do not need to add their own
sleeps.
"""

from __future__ import annotations

import http.client
import random
import threading
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from typing import Any

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Month abbreviations as they appear in INSDSeq date fields, e.g. "24-JUL-2026".
_INSD_MONTHS = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}

# Unauthenticated eutils rate limit is 3 requests/sec -> minimum 1/3 s between
# requests. Shared across every function in this module via a module-level
# lock + timestamp, so concurrent callers (e.g. two functions called from the
# same notebook cell) never jointly exceed the limit.
_MIN_INTERVAL_SECONDS = 1.0 / 3.0
_rate_lock = threading.Lock()
_last_request_monotonic = 0.0


_MAX_RETRIES = 3


def _rate_limited_get(
    endpoint: str, params: dict[str, Any], timeout: float = 60.0, max_retries: int = _MAX_RETRIES
) -> bytes:
    """Issue a rate-limited GET request to an eutils endpoint and return the raw body.

    Internal helper used by every public function in this module -- callers
    should not need to call this directly, but it is not name-mangled so it
    can be reused if a future step needs a raw call this module doesn't yet
    wrap.

    Retries on transient connection failures (e.g. `socket.timeout`,
    `IncompleteRead` -- both observed in practice during long batch runs of
    hundreds of sequential calls) with exponential backoff, since a single
    dropped connection out of hundreds of calls shouldn't kill an entire
    sampling run. Non-transient errors (e.g. a malformed query returning
    HTTP 400) are not retried -- they propagate immediately.

    Args:
        endpoint: eutils CGI script name, e.g. "einfo.fcgi", "esearch.fcgi".
        params: query parameters to URL-encode (db, term, retmax, etc.).
        timeout: socket timeout in seconds for the HTTP request.
        max_retries: number of retry attempts after an initial failure, with
            backoff of 1s, 2s, 4s, ... between attempts.

    Returns:
        The raw response body as bytes (caller parses XML/JSON as needed).
    """
    global _last_request_monotonic
    last_error: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            with _rate_lock:
                elapsed = time.monotonic() - _last_request_monotonic
                if elapsed < _MIN_INTERVAL_SECONDS:
                    time.sleep(_MIN_INTERVAL_SECONDS - elapsed)
                url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
                with urllib.request.urlopen(url, timeout=timeout) as resp:
                    body = resp.read()
                _last_request_monotonic = time.monotonic()
            return body
        except (TimeoutError, OSError, http.client.HTTPException) as exc:
            last_error = exc
            if attempt < max_retries:
                time.sleep(2**attempt)
    assert last_error is not None
    raise last_error


def list_searchable_fields(db: str = "protein") -> list[dict[str, Any]]:
    """Return the full set of esearch-searchable fields for an NCBI database via einfo.

    "Searchable field" here means a field that can be used in an esearch
    `term` query with `[FieldAbbreviation]` syntax (e.g. `spike[Protein
    Name]`), as opposed to a field that only appears inside a *fetched*
    record (efetch) but isn't independently query-able. einfo is the eutils
    endpoint that enumerates these; it does not take a query, just a `db`.

    Args:
        db: NCBI database name, e.g. "protein", "nuccore".

    Returns:
        One dict per field, in einfo's own order, with keys:
            name: short field abbreviation used in queries, e.g. "PROT".
            full_name: human-readable name, e.g. "Protein Name".
            description: einfo's description of what the field indexes.
            term_count: number of distinct indexed terms for this field.
            is_date: True if the field is a date field usable in [Date]-style
                range queries.
            is_numerical: True if the field holds numeric values.
            single_token: True if values are indexed as a single token
                (affects whether phrase quoting matters in queries).
            hierarchy: True if the field supports hierarchical/exploded
                matching (e.g. taxonomy fields with the `[Organism:exp]`
                modifier).
            is_hidden: True if einfo marks the field hidden from normal UI
                field pickers (still usable in queries).

    Example:
        >>> from scripts.ncbi.pull_data import list_searchable_fields
        >>> fields = list_searchable_fields("protein")
        >>> [f["name"] for f in fields if f["is_date"]]
        ['PDAT', 'MDAT']
    """
    body = _rate_limited_get("einfo.fcgi", {"db": db})
    root = ET.fromstring(body)
    fields: list[dict[str, Any]] = []
    for field_el in root.iter("Field"):
        fields.append(
            {
                "name": field_el.findtext("Name"),
                "full_name": field_el.findtext("FullName"),
                "description": field_el.findtext("Description"),
                "term_count": int(field_el.findtext("TermCount") or 0),
                "is_date": field_el.findtext("IsDate") == "Y",
                "is_numerical": field_el.findtext("IsNumerical") == "Y",
                "single_token": field_el.findtext("SingleToken") == "Y",
                "hierarchy": field_el.findtext("Hierarchy") == "Y",
                "is_hidden": field_el.findtext("IsHidden") == "Y",
            }
        )
    return fields


def _esearch(db: str, term: str, **params: Any) -> ET.Element:
    """Run esearch for `term` against `db` and return the parsed XML root.

    Internal helper shared by every function below that needs to run a
    query (as opposed to `list_searchable_fields`, which never queries).

    Args:
        db: NCBI database name, e.g. "protein".
        term: esearch query string, e.g. "txid11118[Organism:exp] AND spike[Title]".
        **params: additional esearch parameters (retmax, retstart, etc.).

    Returns:
        The parsed `<eSearchResult>` XML root element.
    """
    query_params: dict[str, Any] = {"db": db, "term": term, **params}
    body = _rate_limited_get("esearch.fcgi", query_params)
    return ET.fromstring(body)


def search_count(query: str, db: str = "protein") -> int:
    """Return the total number of records matching an esearch query.

    Reads esearch's `<Count>` field directly (`retmax=0`, so no record IDs
    are actually fetched) -- this is the count of *all* matching records,
    not the size of any sample. Generalized over both `db` and the full
    query string, so it works for any field combination (`[Protein Name]`,
    `[Title]`, unqualified free text, taxid filters, `[PDAT]` date ranges,
    etc.), not just the family-wide spike query it was first written for.

    Args:
        query: full esearch term string.
        db: NCBI database name, e.g. "protein".

    Returns:
        Total matching record count.

    Example:
        >>> from scripts.ncbi.pull_data import search_count
        >>> search_count("txid11118[Organism:exp] AND spike[Title]")
        34205
    """
    root = _esearch(db, query, retmax=0)
    return int(root.findtext("Count") or 0)


def sample_titles(query: str, n: int = 20, db: str = "protein", retstart: int = 0) -> list[dict[str, str]]:
    """Fetch `n` record titles matching `query`, via esearch + esummary.

    Grabs the first `n` UIDs esearch returns starting at `retstart` (esearch's
    default result order for db=protein is reverse-chronological by UID, so
    with the default retstart=0 this samples the most recently added records
    -- fine for the eyeballing/title-inspection purpose this is for, not a
    statistically representative sample; use `retstart` to look at a
    different slice if needed). No full records are fetched -- esummary
    returns just the docsum (title, accession, etc.), not sequence data.

    Args:
        query: full esearch term string.
        n: number of titles to fetch.
        db: NCBI database name, e.g. "protein".
        retstart: esearch result offset to start from.

    Returns:
        List of dicts with keys "uid" and "title", in esearch's result order.
        Empty list if the query matches nothing.

    Example:
        >>> from scripts.ncbi.pull_data import sample_titles
        >>> sample_titles("txid11118[Organism:exp] AND peplomer[Title]", n=5)
        [{'uid': '...', 'title': '...'}, ...]
    """
    root = _esearch(db, query, retmax=n, retstart=retstart)
    id_list = root.find("IdList")
    uids = [id_el.text for id_el in id_list.findall("Id")] if id_list is not None else []
    if not uids:
        return []
    body = _rate_limited_get("esummary.fcgi", {"db": db, "id": ",".join(uids), "retmode": "xml"})
    summary_root = ET.fromstring(body)
    results: list[dict[str, str]] = []
    for docsum in summary_root.findall("DocSum"):
        uid = docsum.findtext("Id") or ""
        title = ""
        for item in docsum.findall("Item"):
            if item.get("Name") == "Title":
                title = item.text or ""
                break
        results.append({"uid": uid, "title": title})
    return results


def _random_sample_uids(query: str, sample_size: int, db: str = "protein", chunk_size: int = 50) -> list[str]:
    """Draw a uniform-random sample of UIDs across a query's full result set.

    esearch's default ordering for db=protein is reverse-chronological, so a
    single contiguous page (e.g. `retmax=250, retstart=0`) is a near-single
    point in time, not a representative sample. This instead pulls chunks of
    `chunk_size` contiguous UIDs starting at uniformly-random offsets across
    the whole result set ("cluster random sampling"), pooling chunks until
    `sample_size` unique UIDs are collected.

    Args:
        query: full esearch term string.
        sample_size: target number of unique UIDs to return.
        db: NCBI database name, e.g. "protein".
        chunk_size: UIDs fetched per esearch call.

    Returns:
        Up to `sample_size` unique UIDs (fewer if the query has fewer total
        matches than requested).
    """
    total = search_count(query, db=db)
    if total == 0:
        return []
    max_offset = max(total - chunk_size, 0)
    uids: list[str] = []
    seen: set[str] = set()
    attempts = 0
    max_attempts = sample_size * 4  # generous cap so a pathological query can't loop forever
    while len(uids) < sample_size and attempts < max_attempts:
        attempts += 1
        offset = random.randint(0, max_offset) if max_offset > 0 else 0
        root = _esearch(db, query, retmax=chunk_size, retstart=offset)
        id_list = root.find("IdList")
        if id_list is None:
            continue
        for id_el in id_list.findall("Id"):
            uid = id_el.text
            if uid and uid not in seen:
                seen.add(uid)
                uids.append(uid)
                if len(uids) >= sample_size:
                    break
    return uids


def _parse_insd_date(value: str) -> date:
    """Parse an INSDSeq-style date string, e.g. "24-JUL-2026", into a `date`."""
    day_s, mon_s, year_s = value.split("-")
    return date(int(year_s), _INSD_MONTHS[mon_s.upper()], int(day_s))


def _parse_pdat_date(value: str) -> date:
    """Parse a PDAT-style date string, e.g. "2020/01/01", into a `date`."""
    year_s, mon_s, day_s = value.split("/")
    return date(int(year_s), int(mon_s), int(day_s))


def fetch_create_dates(uids: list[str], db: str = "protein", batch_size: int = 200) -> dict[str, date]:
    """Fetch full records for `uids` and extract each one's `INSDSeq_create-date`.

    Uses efetch with `rettype=gbc&retmode=xml` (the INSDSeq XML format used
    throughout this project for create-date extraction), batched at
    `batch_size` UIDs per request. Relies on efetch returning records in the
    same order as the requested UID list to map dates back to UIDs (standard
    eutils behavior, not independently re-verified here).

    Args:
        uids: list of protein UIDs (GI-style numeric IDs from esearch).
        db: NCBI database name, e.g. "protein".
        batch_size: UIDs fetched per efetch call.

    Returns:
        Dict mapping UID -> parsed `INSDSeq_create-date` as a `date`. A UID
        is omitted if its record has no create-date (not expected, but not
        assumed impossible).
    """
    dates: dict[str, date] = {}
    for i in range(0, len(uids), batch_size):
        batch = uids[i : i + batch_size]
        body = _rate_limited_get(
            "efetch.fcgi", {"db": db, "id": ",".join(batch), "rettype": "gbc", "retmode": "xml"}
        )
        root = ET.fromstring(body)
        seqs = root.findall("INSDSeq")
        for uid, seq in zip(batch, seqs):
            raw = seq.findtext("INSDSeq_create-date")
            if raw:
                dates[uid] = _parse_insd_date(raw)
    return dates


def validate_pdat(
    query: str,
    mindate: str,
    maxdate: str,
    sample_size: int = 250,
    db: str = "protein",
) -> dict[str, Any]:
    """Check whether esearch's server-side `PDAT` date filter agrees with the
    true `INSDSeq_create-date` parsed client-side from full records.

    Method:
      1. Draw a uniform-random sample of `sample_size` UIDs from `query`
         (via `_random_sample_uids`, not a single contiguous page).
      2. Fetch full records for those UIDs and parse each one's true
         `INSDSeq_create-date` (ground truth).
      3. Separately, run `query AND mindate:maxdate[PDAT]` through esearch
         to get its server-side total count for the same range.
      4. For each sampled UID, check membership in the PDAT-filtered query
         directly (`query AND {uid}[UID] AND mindate:maxdate[PDAT]`, count
         0 or 1) and compare that against whether the UID's true create-date
         actually falls in `[mindate, maxdate]`. This is a per-record check,
         not just an aggregate count comparison -- two populations can have
         matching totals while disagreeing on which specific records they
         contain, and this catches that case.

    Args:
        query: full esearch term string (no date filter).
        mindate: range start, PDAT format "YYYY/MM/DD".
        maxdate: range end, PDAT format "YYYY/MM/DD".
        sample_size: number of full records to fetch and check.
        db: NCBI database name, e.g. "protein".

    Returns:
        Dict with:
            total_query_count: total matches for `query` alone.
            esearch_pdat_filtered_count: total matches for `query` with the
                PDAT range filter added (server-side count, whole population).
            sample_size_used: number of sampled UIDs actually checked
                (<= sample_size; fewer only if the query has fewer matches).
            actual_in_range_count: of the sample, how many have a true
                create-date inside [mindate, maxdate].
            agreement_count / agreement_rate: of the sample, how many/what
                fraction have PDAT-filter membership matching their true
                in-range status.
            disagreements: list of dicts (uid, create_date, actual_in_range,
                pdat_says_in_range) for every sampled record where they
                didn't match.
    """
    total = search_count(query, db=db)
    mind = _parse_pdat_date(mindate)
    maxd = _parse_pdat_date(maxdate)
    uids = _random_sample_uids(query, sample_size, db=db)
    create_dates = fetch_create_dates(uids, db=db)
    esearch_pdat_filtered_count = search_count(f"{query} AND {mindate}:{maxdate}[PDAT]", db=db)

    agreement_count = 0
    actual_in_range_count = 0
    disagreements: list[dict[str, Any]] = []
    for uid in uids:
        create_dt = create_dates.get(uid)
        if create_dt is None:
            continue
        actual_in_range = mind <= create_dt <= maxd
        if actual_in_range:
            actual_in_range_count += 1
        member_count = search_count(f"{query} AND {uid}[UID] AND {mindate}:{maxdate}[PDAT]", db=db)
        pdat_says_in_range = member_count >= 1
        if pdat_says_in_range == actual_in_range:
            agreement_count += 1
        else:
            disagreements.append(
                {
                    "uid": uid,
                    "create_date": create_dt.isoformat(),
                    "actual_in_range": actual_in_range,
                    "pdat_says_in_range": pdat_says_in_range,
                }
            )

    sample_size_used = len(create_dates)
    return {
        "query": query,
        "mindate": mindate,
        "maxdate": maxdate,
        "total_query_count": total,
        "esearch_pdat_filtered_count": esearch_pdat_filtered_count,
        "sample_size_used": sample_size_used,
        "actual_in_range_count": actual_in_range_count,
        "agreement_count": agreement_count,
        "agreement_rate": agreement_count / sample_size_used if sample_size_used else None,
        "disagreements": disagreements,
    }
