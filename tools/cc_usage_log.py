#!/usr/bin/env python3
"""Generate a token-usage / prompt / list-price-cost report for a Claude Code session.

Reads a session transcript (~/.claude/projects/<slug>/<session-id>.jsonl) and
writes a Markdown report. Every figure is computed from the transcript; nothing
is hand-entered.

Usage
-----
  tools/cc_usage_log.py                        # latest session in this project
  tools/cc_usage_log.py --list                 # list sessions, newest first
  tools/cc_usage_log.py --session f6ef630b      # id or unique prefix
  tools/cc_usage_log.py --session /path/to.jsonl
  tools/cc_usage_log.py --project ~/work/other -o /tmp/other.md
  tools/cc_usage_log.py --stdout               # print instead of writing

Why dedupe by requestId
-----------------------
One API response is split across several transcript rows (thinking, text,
tool_use), and each row repeats the *same* usage block. Summing rows inflates
every total by roughly two thirds. Totals here are summed over unique
requestId.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, date
from pathlib import Path

# ---------------------------------------------------------------------------
# Pricing — USD per million tokens, list rates.
#
# Keep this table in sync with the model catalog; a stale rate here produces a
# confidently wrong number. Cost is computed per call from the model the call
# actually used, so a session that switched models is priced correctly.
# ---------------------------------------------------------------------------
PRICES = {
    "claude-fable-5":    (10.00, 50.00),
    "claude-mythos-5":   (10.00, 50.00),
    "claude-opus-5":     ( 5.00, 25.00),
    "claude-opus-4-8":   ( 5.00, 25.00),
    "claude-opus-4-7":   ( 5.00, 25.00),
    "claude-opus-4-6":   ( 5.00, 25.00),
    "claude-opus-4-5":   ( 5.00, 25.00),
    "claude-sonnet-5":   ( 3.00, 15.00),
    "claude-sonnet-4-6": ( 3.00, 15.00),
    "claude-sonnet-4-5": ( 3.00, 15.00),
    "claude-haiku-4-5":  ( 1.00,  5.00),
}

# Claude Sonnet 5 introductory pricing, applies to calls on or before this date.
SONNET5_INTRO = (2.00, 10.00)
SONNET5_INTRO_UNTIL = date(2026, 8, 31)

# Fast mode (research preview) is premium-priced on the Opus tier.
FAST_PRICES = {
    "claude-opus-5":   (10.00, 50.00),
    "claude-opus-4-8": (10.00, 50.00),
}

# Cache multipliers applied to a model's input rate.
CACHE_WRITE_1H, CACHE_WRITE_5M, CACHE_READ = 2.00, 1.25, 0.10


def rates(model: str, when: date, fast: bool) -> tuple[float, float] | None:
    """Input/output rate for one call, or None if the model is unpriced."""
    if fast and model in FAST_PRICES:
        return FAST_PRICES[model]
    if model == "claude-sonnet-5" and when <= SONNET5_INTRO_UNTIL:
        return SONNET5_INTRO
    return PRICES.get(model)


# ---------------------------------------------------------------------------
# Transcript location
# ---------------------------------------------------------------------------
def transcript_dir(project: Path) -> Path:
    """Claude Code slugs a project path by replacing '/' with '-'."""
    slug = str(project.resolve()).replace("/", "-")
    return Path.home() / ".claude" / "projects" / slug


def find_sessions(project: Path) -> list[Path]:
    d = transcript_dir(project)
    if not d.is_dir():
        return []
    return sorted(d.glob("*.jsonl"), key=lambda p: p.stat().st_mtime, reverse=True)


def resolve_session(project: Path, want: str | None) -> Path:
    if want and ("/" in want or want.endswith(".jsonl")):
        p = Path(want).expanduser()
        if not p.is_file():
            sys.exit(f"error: no such transcript: {p}")
        return p
    sessions = find_sessions(project)
    if not sessions:
        sys.exit(f"error: no transcripts under {transcript_dir(project)}\n"
                 f"       (is {project} a Claude Code project?)")
    if not want:
        return sessions[0]
    hits = [p for p in sessions if p.stem.startswith(want)]
    if not hits:
        sys.exit(f"error: no session matching {want!r}; try --list")
    if len(hits) > 1:
        sys.exit("error: ambiguous prefix {!r} matches:\n  {}".format(
            want, "\n  ".join(p.stem for p in hits)))
    return hits[0]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
ZERO = dict(inp=0, w1h=0, w5m=0, cr=0, out=0, think=0, n=0, cost=0.0, unpriced=0)


def add(a: dict, b: dict) -> None:
    for k, v in b.items():
        a[k] = a.get(k, 0) + v


def load(path: Path) -> list[dict]:
    rows = []
    for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            print(f"warning: skipping malformed line {n}", file=sys.stderr)
    return rows


def user_prompts(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        if r.get("type") != "user" or r.get("isMeta"):
            continue
        c = r.get("message", {}).get("content")
        txt = (" ".join(b.get("text", "") for b in c
                        if isinstance(b, dict) and b.get("type") == "text")
               if isinstance(c, list) else (c or ""))
        txt = re.sub(r"<system-reminder>.*?</system-reminder>", "", txt, flags=re.S).strip()
        if not txt or txt.startswith("<local-command-stdout>"):
            continue
        out.append({"ts": r.get("timestamp", ""), "text": txt, "u": dict(ZERO)})
    return out


def api_calls(rows: list[dict]) -> tuple[list[dict], int]:
    """Unique API responses, newest-wins-never (first row per requestId)."""
    seen: dict[str, dict] = {}
    n_rows = 0
    for r in rows:
        if r.get("type") != "assistant":
            continue
        n_rows += 1
        # A response with no requestId (rare) is kept on its own uuid.
        key = r.get("requestId") or f"uuid:{r.get('uuid')}"
        seen.setdefault(key, r)
    return sorted(seen.values(), key=lambda r: r.get("timestamp", "")), n_rows


def call_usage(r: dict, unknown: set) -> dict:
    u = r.get("message", {}).get("usage") or {}
    cc = u.get("cache_creation") or {}
    w1h = cc.get("ephemeral_1h_input_tokens", 0)
    w5m = cc.get("ephemeral_5m_input_tokens", 0)
    tot = u.get("cache_creation_input_tokens", 0)
    if w1h + w5m == 0 and tot:
        # No TTL breakdown recorded; attribute to 1h (the pessimistic rate)
        # rather than guessing low.
        w1h = tot
    d = dict(
        inp=u.get("input_tokens", 0), w1h=w1h, w5m=w5m,
        cr=u.get("cache_read_input_tokens", 0),
        out=u.get("output_tokens", 0),
        think=(u.get("output_tokens_details") or {}).get("thinking_tokens", 0),
        n=1, unpriced=0,
    )
    model = r.get("message", {}).get("model", "")
    ts = r.get("timestamp", "")
    when = datetime.fromisoformat(ts.replace("Z", "+00:00")).date() if ts else date.today()
    rt = rates(model, when, (u.get("speed") == "fast"))
    if rt is None:
        unknown.add(model or "<unknown>")
        d["unpriced"] = 1
        d["cost"] = 0.0
    else:
        ir, orr = rt
        d["cost"] = (d["inp"] / 1e6 * ir
                     + d["w1h"] / 1e6 * ir * CACHE_WRITE_1H
                     + d["w5m"] / 1e6 * ir * CACHE_WRITE_5M
                     + d["cr"] / 1e6 * ir * CACHE_READ
                     + d["out"] / 1e6 * orr)
    return d


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def fm(n) -> str:
    return f"{n:,}"


def build(path: Path, project: Path) -> str:
    rows = load(path)
    prompts = user_prompts(rows)
    calls, n_rows = api_calls(rows)
    if not calls:
        sys.exit(f"error: no assistant responses in {path.name}")

    unknown: set[str] = set()
    per_call = [(r, call_usage(r, unknown)) for r in calls]

    T = dict(ZERO)
    for _, d in per_call:
        add(T, d)

    # Attribute each call to the most recent preceding prompt.
    for r, d in per_call:
        idx = 0
        for i, p in enumerate(prompts):
            if p["ts"] and p["ts"] <= r.get("timestamp", ""):
                idx = i
        if prompts:
            add(prompts[idx]["u"], d)

    models = sorted({r.get("message", {}).get("model", "?") for r in calls})
    speeds = sorted({(r.get("message", {}).get("usage") or {}).get("speed", "standard")
                     for r in calls})
    prompt_tok = T["inp"] + T["w1h"] + T["w5m"] + T["cr"]

    # No-cache counterfactual, priced per call so mixed-model sessions are right.
    nocache = 0.0
    for r, d in per_call:
        model = r.get("message", {}).get("model", "")
        ts = r.get("timestamp", "")
        when = datetime.fromisoformat(ts.replace("Z", "+00:00")).date() if ts else date.today()
        rt = rates(model, when, ((r.get("message", {}).get("usage") or {}).get("speed") == "fast"))
        if rt:
            nocache += ((d["inp"] + d["w1h"] + d["w5m"] + d["cr"]) / 1e6 * rt[0]
                        + d["out"] / 1e6 * rt[1])

    peak = max(d["inp"] + d["cr"] + d["w1h"] + d["w5m"] for _, d in per_call)
    first = prompts[0]["ts"] if prompts else calls[0].get("timestamp", "")
    last = calls[-1].get("timestamp", "")
    mins = 0.0
    if first and last:
        mins = (datetime.fromisoformat(last.replace("Z", "+00:00"))
                - datetime.fromisoformat(first.replace("Z", "+00:00"))).total_seconds() / 60

    L: list[str] = []
    w = L.append
    w("# Claude Code session usage log")
    w("")
    w("Token usage, prompts, and list-price cost for one Claude Code session,")
    w("computed from its transcript. Regenerate with `tools/cc_usage_log.py`.")
    w("")
    w("| | |")
    w("|---|---|")
    w(f"| Session ID | `{path.stem}` |")
    w(f"| Project | `{project}` |")
    w(f"| Model(s) | {', '.join(f'`{m}`' for m in models)} |")
    if speeds != ["standard"]:
        w(f"| Speed | {', '.join(speeds)} |")
    w(f"| Session start | {first[:19]}Z |")
    w(f"| Snapshot taken | {last[:19]}Z |")
    w(f"| Wall clock | {mins:.0f} min |")
    w(f"| User prompts | {len(prompts)} |")
    w(f"| API calls | {len(calls)} (from {n_rows} transcript rows) |")
    w("")
    w("> **This is a list-price estimate, not a bill.** A Claude Code")
    w("> subscription is not metered per token — the figures below are what the")
    w("> same traffic would cost at published API rates. Read them as a relative")
    w("> measure of where tokens went, not as an invoice.")
    w("")
    if unknown:
        w(f"> ⚠️ **Unpriced model(s):** {', '.join(sorted(unknown))} — "
          f"{T['unpriced']} call(s) counted in the token totals but contributing")
        w("> $0 to cost. Add the rate to `PRICES` in the generator to fix.")
        w("")
    w("---")
    w("")
    w("## 1. Token usage")
    w("")
    w("| Category | Tokens | Multiplier | Cost |")
    w("|---|---:|---:|---:|")
    w(f"| Uncached input | {fm(T['inp'])} | 1.00x | — |")
    w(f"| Cache writes (1 h TTL) | {fm(T['w1h'])} | {CACHE_WRITE_1H:.2f}x | — |")
    w(f"| Cache writes (5 min TTL) | {fm(T['w5m'])} | {CACHE_WRITE_5M:.2f}x | — |")
    w(f"| Cache reads | {fm(T['cr'])} | {CACHE_READ:.2f}x | — |")
    w(f"| Output | {fm(T['out'])} | output rate | — |")
    w(f"| **Total** | **{fm(prompt_tok + T['out'])}** | | **${T['cost']:.2f}** |")
    w("")
    w("Per-category cost is left blank because rates are applied per call from")
    w("each call's own model; only the total is meaningful when a session spans")
    w("more than one model. For a single-model session the split is the token")
    w("count times that model's rate and multiplier.")
    w("")
    if T["out"]:
        w(f"Of {fm(T['out'])} output tokens, {fm(T['think'])} "
          f"({T['think'] / T['out'] * 100:.0f}%) were thinking tokens.")
        w("")
    w("### What caching saved")
    w("")
    if prompt_tok:
        w(f"- Prompt tokens processed: **{fm(prompt_tok)}**")
        w(f"- Served from cache: **{fm(T['cr'])}** — a "
          f"**{T['cr'] / prompt_tok * 100:.1f}%** hit rate")
        w(f"- Same traffic with no caching: **${nocache:.2f}**")
        if nocache > 0:
            w(f"- Actual: **${T['cost']:.2f}** — caching avoided "
              f"**${nocache - T['cost']:.2f}** "
              f"({(1 - T['cost'] / nocache) * 100:.0f}%)")
        w("")
    w(f"Largest single-request context: **{fm(peak)}** tokens.")
    w("")
    w("---")
    w("")
    w("## 2. Prompts")
    w("")
    w("Verbatim user prompts in order. Cost is attributed to the turn each")
    w("prompt opened.")
    w("")
    for i, p in enumerate(prompts, 1):
        u = p["u"]
        w(f"### {i}. {p['ts'][11:19]}Z")
        w("")
        w("```")
        w(p["text"].strip())
        w("```")
        w("")
        if u["n"] == 0:
            w("*No API call (slash command or local output).*")
        else:
            w(f"{u['n']} API calls · {fm(u['w1h'] + u['w5m'])} cache write · "
              f"{fm(u['cr'])} cache read · {fm(u['out'])} output "
              f"({fm(u['think'])} thinking) · **${u['cost']:.2f}**")
        w("")
    w("---")
    w("")
    w("## 3. Cost by turn")
    w("")
    w("| # | Prompt | Calls | Cache write | Cache read | Output | Cost |")
    w("|---|---|---:|---:|---:|---:|---:|")
    for i, p in enumerate(prompts, 1):
        u = p["u"]
        lbl = re.sub(r"\s+", " ", p["text"])[:44]
        m = re.match(r"<command-name>([^<]+)</command-name>", p["text"])
        if m:
            lbl = f"{m.group(1)} (slash command)"
        lbl = lbl.replace("|", "\\|")
        w(f"| {i} | {lbl} | {u['n']} | {fm(u['w1h'] + u['w5m'])} | "
          f"{fm(u['cr'])} | {fm(u['out'])} | ${u['cost']:.2f} |")
    w(f"| | **Total** | **{T['n']}** | **{fm(T['w1h'] + T['w5m'])}** | "
      f"**{fm(T['cr'])}** | **{fm(T['out'])}** | **${T['cost']:.2f}** |")
    w("")
    w("---")
    w("")
    w("## 4. Method and caveats")
    w("")
    w(f"**Source.** `{path}`")
    w("")
    w(f"**Deduplication.** {n_rows} assistant rows collapse to {len(calls)} "
      "distinct `requestId`s: one API response is split across thinking, text, "
      "and tool_use rows, each repeating the *same* usage block. Totals are "
      "summed over unique `requestId` — summing rows inflates every figure by "
      "roughly two thirds.")
    w("")
    w("**Pricing.** List rates per million tokens, selected per call from that "
      "call's own model (and `usage.speed`, so fast-mode calls price at the "
      "premium rate). Cache multipliers on the input rate: writes "
      f"{CACHE_WRITE_1H:.2f}x at 1-hour TTL and {CACHE_WRITE_5M:.2f}x at "
      f"5-minute TTL, reads {CACHE_READ:.2f}x. A call whose model is absent "
      "from the table contributes tokens but $0, and is flagged at the top "
      "rather than silently priced.")
    w("")
    w("**Caveats.**")
    w("")
    w("- Not a bill (see the note above).")
    w("- The snapshot ends at the last API call in the transcript, so the calls "
      "that generate and save this file are not counted. The log cannot measure "
      "its own writing; re-run later to catch up.")
    w("- A call is attributed to the most recent preceding prompt, so a prompt "
      "sent mid-turn takes the work from that point on. Turn boundaries are "
      "approximate.")
    w("- Thinking tokens are billed as output and are included in the output row.")
    w("- Where no cache TTL breakdown is recorded, writes are charged at the "
      "1-hour multiplier — the pessimistic choice, so the total is never "
      "understated.")
    w("")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Report token usage, prompts, and list-price cost for a Claude Code session.")
    ap.add_argument("--project", type=Path, default=Path.cwd(),
                    help="project directory whose transcripts to read (default: cwd)")
    ap.add_argument("--session", metavar="ID|PATH",
                    help="session id, unique id prefix, or path to a .jsonl (default: newest)")
    ap.add_argument("-o", "--output", type=Path,
                    help="output file (default: <project>/cc-usage-log.md)")
    ap.add_argument("--stdout", action="store_true", help="print instead of writing")
    ap.add_argument("--list", action="store_true", dest="do_list",
                    help="list this project's sessions, newest first, and exit")
    args = ap.parse_args()

    if args.do_list:
        sessions = find_sessions(args.project)
        if not sessions:
            sys.exit(f"no transcripts under {transcript_dir(args.project)}")
        print(f"{'session id':38}  {'modified':19}  {'size':>9}")
        for p in sessions:
            st = p.stat()
            print(f"{p.stem:38}  "
                  f"{datetime.fromtimestamp(st.st_mtime):%Y-%m-%d %H:%M:%S}  "
                  f"{st.st_size / 1024:8.0f}K")
        return

    path = resolve_session(args.project, args.session)
    report = build(path, args.project.resolve())

    if args.stdout:
        sys.stdout.write(report)
        return
    out = args.output or (args.project / "cc-usage-log.md")
    out.write_text(report)
    print(f"wrote {out}  (session {path.stem[:8]})")


if __name__ == "__main__":
    main()
