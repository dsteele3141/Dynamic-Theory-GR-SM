#!/usr/bin/env python3
# ============================================================================
# linter.py  --  Dynamic Theory: Sage convention linter
#
# The Sage analogue of AxiomCheck (the Lean axiom audit): a meta-check over
# ALL Sage scripts in the corpus.
#
# Scans every  *.sage  in
#     Layer1/Sage , Layer2/Sage , Layer3/Sage
# and the matching  <name>.out.txt  in each layer's  sage-results/ .
# Figure-generator scripts (plot_*.sage) emit PNGs rather than a verification
# result and are skipped.
#
# Checks
#   * Signature (+,-,-,-,-): a 5D mass component g[4,4] (or a [1,-1,-1,-1,x]
#     signature list) must be SPACELIKE (g44 < 0); g44 = +1 is FAIL.
#   * Dirac: where the 5D Clifford algebra is built, the corrected COMPLEX
#     A-matrix [[0,iI],[-iI,0]] must be present (the real [[0,I],[I,0]] fails
#     anticommutation).
#   * Self-energy (T2.6): the corrected cube-root R0 = (...)^(1/3) (Eq 103)
#     and/or the corrected coefficient 3 e^2/(16 pi eps0 R0) (factor of 3,
#     vs the manuscript's printed 1/16) must be present.
#   * Verdict honesty: a literal [PASS]/[OK]/[VERIFIED]/[CONFIRMED] printed next
#     to an interpolated value/residual, NOT derived from a computed boolean,
#     is flagged WARN (the hard-coded-verdict anti-pattern; this is the textual
#     signal of a print-diagnostic [D] script vs an assertion-driven [A] one).
#   * Sage-API regression guards: .symmetric_tensor_field(...),
#     block_matrix(...).matrix(), connection().christoffel_symbols(),
#     .automorphism_field(0,2,...), <comp>.simplify_full() without .expr().
#   * Output-path convention: scripts must rely on stdout capture into the
#     layer's sage-results/ -- NO self-managed output path. Any residual
#     os.makedirs(...), Tee stdout redirection, 'Stage1_Postulates',
#     '02_sage-results' or '04_codex_code(_results)' path is flagged.
#   * Runtime: each script's matching sage-results/<name>.out.txt must exist
#     and contain no traceback / Error / [FAIL].
#
# Pure Python (no Sage needed):  python3 linter.py
# Writes linter.out.txt (repo root) and prints to stdout.
# Exit 0 if clean (no WARN/FAIL), 1 otherwise.
# ============================================================================

import os
import re
import glob
import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))   # repo root
ROOT = SCRIPT_DIR
LAYERS = ["Layer1", "Layer2", "Layer3"]
CODE_SUBDIR = "Sage"
RESULTS_SUBDIR = "sage-results"
REPORT_PATH = os.path.join(ROOT, "linter.out.txt")

SEV_OK, SEV_INFO, SEV_WARN, SEV_FAIL = "OK", "INFO", "WARN", "FAIL"

findings = []          # (layer, sev, fname, line_or_None, msg)
g44_records = []       # (layer, fname, line, value, spacelike)


def add(layer, sev, fname, line, msg):
    findings.append((layer, sev, fname, line, msg))


def lineno(text, idx):
    return text.count("\n", 0, idx) + 1


# ----------------------------------------------------------------------------
# Static checks (operate on one .sage file's text)
# ----------------------------------------------------------------------------

def check_signature(layer, fname, text):
    for m in re.finditer(r'(?m)^[ \t]*([A-Za-z_]\w*)\s*\[\s*4\s*,\s*4\s*\]\s*=\s*([^#\n]+)', text):
        ln = lineno(text, m.start())
        val = m.group(2).strip().rstrip(';').strip()
        spacelike = val.startswith('-')
        g44_records.append((layer, fname, ln, val, spacelike))
        if val in ('1', '+1'):
            add(layer, SEV_FAIL, fname, ln,
                "g44 = +1 (timelike mass) violates adjudicated (+,-,-,-,-) signature; expected g44 < 0.")
        elif spacelike:
            add(layer, SEV_OK, fname, ln, "signature: g44 = %s (spacelike mass)." % val)
        else:
            add(layer, SEV_WARN, fname, ln, "g44 = %s : confirm spacelike (g44 < 0)." % val)
    for m in re.finditer(r'\[\s*1\s*,\s*-1\s*,\s*-1\s*,\s*-1\s*,\s*(\+?1|-1)\s*\]', text):
        ln = lineno(text, m.start())
        last = m.group(1)
        spacelike = (last == '-1')
        g44_records.append((layer, fname, ln, "signature-list[%s]" % last, spacelike))
        if spacelike:
            add(layer, SEV_OK, fname, ln, "signature list [1,-1,-1,-1,-1]: spacelike mass.")
        else:
            add(layer, SEV_FAIL, fname, ln,
                "signature list ends in +1 (timelike mass); expected -1 (spacelike).")


def check_dirac_A(layer, fname, text):
    triggers = bool(re.search(r'A_(?:source|correct)\s*=\s*block', text)) or ('check_clifford' in text) \
        or bool(re.search(r'\bA\s*=\s*(?:block_matrix|matrix)', text)) and 'dirac' in fname.lower()
    if not triggers:
        return
    has_complex = bool(re.search(r'-\s*[iI]\s*\*\s*I2', text)) and bool(re.search(r'(?<!-)\b[iI]\s*\*\s*I2', text))
    has_real_source = bool(re.search(r'A_source\s*=\s*block', text))
    if has_complex:
        add(layer, SEV_OK, fname, None, "dirac: corrected complex A-matrix [[0,iI],[-iI,0]] present.")
    elif has_real_source:
        add(layer, SEV_FAIL, fname, None,
            "dirac: real A-matrix used without the corrected imaginary A; fails Clifford anticommutation.")
    else:
        add(layer, SEV_WARN, fname, None,
            "dirac: Clifford construction detected but corrected imaginary A-matrix not confirmed.")


def check_self_energy(layer, fname, text):
    if 'self_energy' not in fname.lower() and 'T2.6' not in text and 'T2_6' not in text:
        return
    despaced = re.sub(r'\s+', '', text)
    has_cube = bool(re.search(r'\^\(\s*1\s*/\s*3\s*\)', text)) or 'cbrt' in text or '(1/3)' in despaced
    has_coeff = bool(re.search(r'3\*?e\^?2?/\(16', despaced))
    if has_cube:
        add(layer, SEV_OK, fname, None, "self-energy: corrected cube-root R0 = (...)^(1/3) present (Eq 103).")
    if has_coeff:
        add(layer, SEV_OK, fname, None, "self-energy: corrected coefficient 3 e^2/(16 pi eps0 R0) present (factor of 3).")
    if not has_cube and not has_coeff:
        add(layer, SEV_WARN, fname, None,
            "self-energy: neither corrected cube-root R0 nor 3/16 coefficient detected.")


def check_verdict_honesty(layer, fname, text):
    """Merge of both upstream heuristics: a literal verdict glued to an
    interpolated value/residual, not derived via 'if ... else', is hard-coded."""
    for i, line in enumerate(text.splitlines(), start=1):
        derived = (' if ' in line and ' else ' in line)
        if derived:
            continue
        # (A) sympy-style:  {value} [PASS]   or   %s ... [PASS]
        if re.search(r'\}\s*\[(PASS|FAIL|OK|VERIFIED|CONFIRMED)\]', line):
            add(layer, SEV_WARN, fname, i,
                "verdict hard-coded next to an interpolated value ('{..} [PASS]'); "
                "derive it: [{'PASS' if cond else 'FAIL'}].")
            continue
        # (B) sage-style: literal verdict on a line that also shows a residual/comparison
        if re.search(r'\[(PASS|OK|VERIFIED|CONFIRMED)\]', line):
            shows_result = ('%s' in line or re.search(r'\{[^}]+\}', line)) and \
                           (re.search(r'==\s*0', line) or '.simplify_full()' in line or
                            re.search(r'[Dd]ifference|[Rr]esidual', line))
            if shows_result:
                add(layer, SEV_WARN, fname, i,
                    "verdict hard-coded next to a computed residual/comparison; "
                    "derive it: %% ('VERIFIED' if cond else 'FAIL').")


def check_api_patterns(layer, fname, text):
    lines = text.splitlines()
    conn_vars = set(re.findall(r'(\w+)\s*=\s*[\w.]+\.connection\(\)', text))
    for i, line in enumerate(lines, start=1):
        if re.search(r'\.symmetric_tensor_field\s*\(', line):
            add(layer, SEV_FAIL, fname, i,
                "API: .symmetric_tensor_field(...) does not exist; use .tensor_field(k,l,name,sym=[(0,1)]).")
        if re.search(r'block_matrix\s*\(.*\)\s*\.\s*matrix\s*\(\s*\)', line):
            add(layer, SEV_FAIL, fname, i,
                "API: block_matrix(...).matrix() is invalid; use block_matrix(..., subdivide=False).")
        if re.search(r'\.connection\(\)\s*\.\s*christoffel_symbols\s*\(', line):
            add(layer, SEV_FAIL, fname, i,
                "API: connection().christoffel_symbols() is invalid; use connection().coef().")
        for v in conn_vars:
            if re.search(r'\b' + re.escape(v) + r'\s*\.\s*christoffel_symbols\s*\(', line):
                add(layer, SEV_FAIL, fname, i,
                    "API: '%s' is a connection; .christoffel_symbols() invalid -> use .coef()." % v)
        if re.search(r'\.automorphism_field\s*\(\s*\d', line):
            add(layer, SEV_FAIL, fname, i,
                "API: automorphism_field takes no (k,l) tensor type; use .tensor_field(...).")
        if re.search(r'\b(Ric|Riem|g|gd|gh|gs|gg|gm|g_metric|gd_sigma)\s*\[[^\]]*\]\s*\.\s*simplify_full\s*\(', line):
            add(layer, SEV_WARN, fname, i,
                "API: simplify_full() on a metric/curvature component (chart function); call .expr() first.")


def check_path_convention(layer, fname, text):
    """Scripts must rely on stdout capture into the layer's sage-results/,
    with NO self-managed output path."""
    flagged = False
    if 'Stage1_Postulates' in text:
        ln = lineno(text, text.index('Stage1_Postulates'))
        add(layer, SEV_FAIL, fname, ln, "path: writes to 'Stage1_Postulates' tree (wrong repo).")
        flagged = True
    if re.search(r"makedirs\(\s*['\"]/", text):
        add(layer, SEV_FAIL, fname, None, "path: hardcoded absolute output path in makedirs().")
        flagged = True
    if re.search(r'(?m)^[ \t]*os\.makedirs\s*\(', text):
        add(layer, SEV_WARN, fname, None,
            "path: residual os.makedirs(...) self-write; should rely on stdout capture into sage-results/.")
        flagged = True
    if re.search(r"['\"]04_codex_code(_results)?['\"]|['\"]02_sage-results['\"]", text):
        add(layer, SEV_WARN, fname, None,
            "path: residual old audit-corpus path ('04_codex_code*' / '02_sage-results') referenced.")
        flagged = True
    if re.search(r'(?m)^[ \t]*sys\.stdout\s*=\s*Tee', text):
        add(layer, SEV_WARN, fname, None,
            "path: residual 'sys.stdout = Tee(...)' redirection; should rely on plain stdout capture.")
        flagged = True
    if not flagged:
        add(layer, SEV_OK, fname, None, "path: stdout-capture convention (no self-managed output path).")


# ----------------------------------------------------------------------------
# Runtime check (per-script matching result file)
# ----------------------------------------------------------------------------

def check_runtime(layer, fname, res_dir):
    base = fname[:-5] if fname.endswith(".sage") else fname
    path = os.path.join(res_dir, base + ".out.txt")
    if not os.path.isfile(path):
        add(layer, SEV_WARN, base + ".out.txt", None, "runtime: no matching result file found.")
        return
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except OSError as exc:
        add(layer, SEV_WARN, base + ".out.txt", None, "could not read result file: %s" % exc)
        return
    if "Traceback (most recent call last)" in content or re.search(r'(?m)^\w*Error:', content):
        add(layer, SEV_FAIL, base + ".out.txt", None, "runtime: result contains a traceback / Error (crashed).")
    elif "[FAIL]" in content:
        add(layer, SEV_FAIL, base + ".out.txt", None, "runtime: result contains a [FAIL] verdict.")
    else:
        add(layer, SEV_OK, base + ".out.txt", None, "runtime: no traceback/Error/FAIL in result output.")


# ----------------------------------------------------------------------------
# Cadabra cross-checks (independent second engine; Layer*/Cadabra/*.cdb)
# ----------------------------------------------------------------------------

def check_cadabra(layer, fname, res_dir):
    base = fname[:-4]
    path = os.path.join(res_dir, base + ".out.txt")
    if not os.path.isfile(path):
        add(layer, SEV_WARN, base + ".out.txt", None, "cadabra: no matching result file found.")
        return
    with open(path, encoding="utf-8", errors="replace") as fh:
        content = fh.read()
    n_ok = len(re.findall(r"(?m)^CHECK .*: OK$", content))
    if "Traceback (most recent call last)" in content or "Script ended execution" in content:
        add(layer, SEV_FAIL, base + ".out.txt", None, "cadabra: result contains a traceback (assertion failed).")
    elif n_ok == 0 or "SUMMARY" not in content:
        add(layer, SEV_FAIL, base + ".out.txt", None, "cadabra: no CHECK ... OK lines / no SUMMARY (incomplete run).")
    else:
        add(layer, SEV_OK, base + ".out.txt", None, "cadabra: %d assertion check(s) OK, run complete." % n_ok)


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def main():
    n_files = 0
    n_skipped = 0
    n_cdb = 0
    for layer in LAYERS:
        code_dir = os.path.join(ROOT, layer, CODE_SUBDIR)
        res_dir = os.path.join(ROOT, layer, RESULTS_SUBDIR)
        for path in sorted(glob.glob(os.path.join(code_dir, "*.sage"))):
            fname = os.path.basename(path)
            if fname.startswith("plot_"):      # figure generators emit PNGs, not a result
                n_skipped += 1
                continue
            n_files += 1
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            check_signature(layer, fname, text)
            check_dirac_A(layer, fname, text)
            check_self_energy(layer, fname, text)
            check_verdict_honesty(layer, fname, text)
            check_api_patterns(layer, fname, text)
            check_path_convention(layer, fname, text)
            check_runtime(layer, fname, res_dir)
        cdb_dir = os.path.join(ROOT, layer, "Cadabra")
        cdb_res = os.path.join(ROOT, layer, "cadabra-results")
        for path in sorted(glob.glob(os.path.join(cdb_dir, "*.cdb"))):
            n_cdb += 1
            check_cadabra(layer, os.path.basename(path), cdb_res)

    # cross-file signature consistency
    if g44_records:
        timelike = [r for r in g44_records if not r[4]]
        if timelike:
            add("(repository)", SEV_FAIL, "(cross-file)", None,
                "signature DRIFT: %d g44 assignment(s) not spacelike." % len(timelike))
        else:
            add("(repository)", SEV_OK, "(cross-file)", None,
                "signature consistent: all %d g44 definition(s) spacelike (g44 < 0)." % len(g44_records))

    n_ok = sum(1 for f in findings if f[1] == SEV_OK)
    n_info = sum(1 for f in findings if f[1] == SEV_INFO)
    n_warn = sum(1 for f in findings if f[1] == SEV_WARN)
    n_fail = sum(1 for f in findings if f[1] == SEV_FAIL)

    sev_order = {SEV_FAIL: 0, SEV_WARN: 1, SEV_INFO: 2, SEV_OK: 3}
    out = []
    out.append("=" * 72)
    out.append("DYNAMIC THEORY -- SAGE LINTER  (linter.py)")
    out.append("Run: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    present = [l for l in LAYERS if os.path.isdir(os.path.join(ROOT, l, CODE_SUBDIR))]
    out.append("Scanned: %s (repository root)" % os.path.basename(ROOT))
    out.append("Sage scripts: %d  (across %s; %d plot_* skipped)" % (n_files, ", ".join(present), n_skipped))
    out.append("Cadabra cross-checks: %d" % n_cdb)
    out.append("=" * 72)

    for layer in LAYERS + ["(repository)"]:
        layer_findings = [f for f in findings if f[0] == layer]
        if not layer_findings:
            continue
        out.append("")
        out.append("### %s" % layer)
        by_file = {}
        for _, sev, fname, line, msg in layer_findings:
            by_file.setdefault(fname, []).append((sev, line, msg))
        for fname in sorted(by_file):
            out.append(fname)
            for sev, line, msg in sorted(by_file[fname], key=lambda t: (sev_order[t[0]], t[1] or 0)):
                loc = (":%d" % line) if line else ""
                out.append("  [%-4s]%s %s" % (sev, loc, msg))

    out.append("")
    out.append("=" * 72)
    out.append("SUMMARY    OK: %d    INFO: %d    WARN: %d    FAIL: %d" % (n_ok, n_info, n_warn, n_fail))
    clean = (n_warn == 0 and n_fail == 0)
    out.append("RESULT: %s" % ("CLEAN -- no inconsistencies detected."
                               if clean else
                               "INCONSISTENCIES FOUND (see WARN/FAIL above)."))
    out.append("=" * 72)

    report = "\n".join(out) + "\n"
    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print("Report written to: %s" % os.path.normpath(REPORT_PATH))
    return 0 if clean else 1


if __name__ == "__main__":
    raise SystemExit(main())
