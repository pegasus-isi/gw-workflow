#!/usr/bin/env python3
"""Fetch one detector's public GWOSC strain segment for one GWTC-1 event.

SPEC.md D1-D3, D8. Uses the documented GWOSC event API to resolve the event's
GPS time and its strain-file URLs -- GPS times are deliberately NOT taken from
references/gwtc1_reference.json, which omits them rather than asserting derived
values (SPEC.md D8).

Selects the 4096 s / 4096 Hz HDF5 file (SPEC.md D2-D3): the 4 kHz data is valid
to ~1700 Hz, comfortably above this analysis's 512 Hz ceiling, and the 4096 s
file is required because a 32 s snippet cannot support 16 s PSD blocks with a
16 s inverse-PSD truncation.
"""

import argparse
import json
import logging
import os
import sys
import time

import requests

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# The catalog segment is REQUIRED: /eventapi/json/{event}/ returns 404.
# Verified against the live service; the root lists the available catalogs.
ROOT_API = "https://gwosc.org/eventapi/json/"
EVENT_API = "https://gwosc.org/eventapi/json/{catalog}/{event}/"
DEFAULT_CATALOG = "GWTC-1-confident"     # this workflow's scope (SPEC.md N12)
TIMEOUT = 60          # SPEC: explicit timeout on every request
RETRIES = 4


def _get(url, allow_missing=False, **kw):
    """GET with an explicit timeout, status check, and bounded backoff.

    CLAUDE.md §4: GWOSC is a shared public service. Concurrency is capped by the
    workflow's DAGMan category; this adds per-request politeness.

    Client errors other than 429 are NOT retried -- a 404 is an answer, and
    retrying it four times with backoff wastes 15 s and pointlessly re-hits the
    service. `allow_missing` returns None on 404 so a caller can probe.
    """
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.get(url, timeout=TIMEOUT, **kw)
            if r.status_code == 404 and allow_missing:
                return None
            if 400 <= r.status_code < 500 and r.status_code != 429:
                r.raise_for_status()          # definitive; do not retry
            r.raise_for_status()
            return r
        except requests.HTTPError as exc:
            code = exc.response.status_code if exc.response is not None else None
            if code and 400 <= code < 500 and code != 429:
                raise
            last = exc
        except Exception as exc:
            last = exc
        wait = 2 ** attempt
        logger.warning("GET %s failed (%s); retry in %ds", url, last, wait)
        time.sleep(wait)
    raise RuntimeError(f"giving up on {url}: {last}")


def resolve_event(event, catalog=None):
    """Return (gps, strain_file_list, api_name, catalog).

    Tries the workflow's own catalog first, then falls back to enumerating the
    catalogs the API advertises -- so a re-catalogued event still resolves
    instead of failing on a hard-coded path.
    """
    tried = []
    for cat in ([catalog] if catalog else [DEFAULT_CATALOG]):
        tried.append(cat)
        r = _get(EVENT_API.format(catalog=cat, event=event), allow_missing=True)
        if r is not None:
            return _parse_event(r.json(), event, cat)

    if catalog is None:
        logger.warning("%s not in %s; enumerating catalogs", event, DEFAULT_CATALOG)
        for cat in (_get(ROOT_API).json() or {}):
            if cat in tried:
                continue
            tried.append(cat)
            r = _get(EVENT_API.format(catalog=cat, event=event), allow_missing=True)
            if r is not None:
                return _parse_event(r.json(), event, cat)

    raise RuntimeError(
        f"{event!r} not found in any GWOSC catalog (tried: {', '.join(tried)})")


def _parse_event(data, event, catalog):
    events = data.get("events") or {}
    if not events:
        raise RuntimeError(f"GWOSC event API returned no events for {event!r}")
    # Prefer the highest catalogue version (keys look like "GW150914-v3").
    def version(key):
        tail = key.rsplit("-v", 1)
        try:
            return int(tail[1])
        except (IndexError, ValueError):
            return -1
    name = sorted(events, key=version)[-1]
    rec = events[name]
    gps = rec.get("GPS")
    if gps is None:
        raise RuntimeError(f"no GPS time in event API record for {name}")
    return float(gps), rec.get("strain") or [], name, catalog


def pick_strain_file(strain, detector, duration, sample_rate):
    """Select the exact file SPEC.md D2/D3 calls for."""
    cands = [
        s for s in strain
        if s.get("detector") == detector
        and str(s.get("format", "")).lower() in ("hdf5", "hdf")
        and int(s.get("duration", 0)) == duration
        and int(s.get("sampling_rate", 0)) == sample_rate
    ]
    if not cands:
        have = sorted({
            (s.get("detector"), s.get("duration"), s.get("sampling_rate"),
             s.get("format")) for s in strain
        })
        raise RuntimeError(
            f"no {duration}s/{sample_rate}Hz hdf5 file for {detector}; "
            f"API offers {have}"
        )
    return cands[0]


def main():
    p = argparse.ArgumentParser(description="Fetch GWOSC strain for one event/detector")
    p.add_argument("--event", required=True, help="Event name, e.g. GW150914")
    p.add_argument("--detector", required=True, help="H1, L1 or V1")
    p.add_argument("--duration", type=int, default=4096, help="Segment length (s)")
    p.add_argument("--sample-rate", type=int, default=4096, help="Sample rate (Hz)")
    p.add_argument("--catalog", default=None,
                   help=f"GWOSC catalog (default: {DEFAULT_CATALOG}, with "
                        f"fallback to enumerating all catalogs)")
    p.add_argument("--output", required=True, help="Output HDF5 path")
    p.add_argument("--meta-output", required=True, help="Output provenance JSON")
    args = p.parse_args()

    for path in (args.output, args.meta_output):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    logger.info("Resolving %s from the GWOSC event API", args.event)
    gps, strain, api_name, catalog = resolve_event(args.event, args.catalog)
    logger.info("%s -> GPS %.4f (record %s in catalog %s)",
                args.event, gps, api_name, catalog)

    entry = pick_strain_file(strain, args.detector, args.duration, args.sample_rate)
    url = entry["url"]
    logger.info("Downloading %s", url)

    r = _get(url, stream=True)
    total = 0
    with open(args.output, "wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 20):
            if chunk:
                fh.write(chunk)
                total += len(chunk)
    logger.info("Wrote %s (%.1f MB)", args.output, total / 1e6)

    if total == 0:
        logger.error("downloaded 0 bytes from %s", url)
        sys.exit(1)

    # SPEC.md V5c: record provenance for every analysed segment.
    meta = {
        "event": args.event,
        "gwosc_api_record": api_name,
        "gwosc_catalog": catalog,
        "detector": args.detector,
        "event_gps": gps,
        "url": url,
        "gps_start": entry.get("GPSstart"),
        "duration": entry.get("duration"),
        "sample_rate": entry.get("sampling_rate"),
        "bytes": total,
        "gps_source": "GWOSC event API (SPEC.md D8)",
    }
    gps_start = entry.get("GPSstart")
    if gps_start is not None:
        before = gps - float(gps_start)
        after = float(gps_start) + float(entry.get("duration", 0)) - gps
        meta["event_offset_in_file_s"] = round(before, 3)
        meta["margin_after_event_s"] = round(after, 3)
        logger.info("Event sits %.0fs into the file, %.0fs from its end",
                    before, after)
        if min(before, after) < 200:
            logger.warning("only %.0fs of margin on one side; a long "
                           "--analysis-length may not fit (see condition_strain)",
                           min(before, after))
    with open(args.meta_output, "w") as fh:
        json.dump(meta, fh, indent=2)

    if not os.path.exists(args.output):
        logger.error("expected output missing: %s", args.output)
        sys.exit(1)
    logger.info("Output: %s, %s", args.output, args.meta_output)


if __name__ == "__main__":
    main()
