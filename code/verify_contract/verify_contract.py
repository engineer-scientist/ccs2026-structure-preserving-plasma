#
#   VERIFY{} contract: generated PETSc programs print one machine-readable line
#       VERIFY {"code": ..., "invariants": {"<name>": <value>, ...}, "finite": true}
#   and the pipeline (compile-run / analysis agent) accepts the code only if every
#   invariant in the contract is within tolerance. Verification, not LLM confidence, decides.
#
import json, re

VERIFY_RE = re.compile(r'^VERIFY\s+(\{.*\})\s*$', re.M)

# default tolerances for structure-preserving plasma codes (overridable per problem by NA agent)
DEFAULT_CONTRACT = {
  'energy_rel_drift_max': 1e-6,     # discrete-gradient / exact-energy methods: ~1e-12; symplectic: bounded O(dt^2)
  'p_phi_rel_drift_max':  1e-2,     # canonical toroidal momentum (axisymmetric fields)
  'divB_max':             1e-12,    # constrained transport / FEEC: machine precision
  'charge_rel_drift_max': 1e-12,    # PIC charge conservation
  'phase_volume_det_minus_1': 1e-8, # symplectic / volume-preserving maps
}

def parse_verify(stdout: str) -> list:
  out = []
  for m in VERIFY_RE.finditer(stdout or ''):
    try:
      out.append(json.loads(m.group(1)))
    except json.JSONDecodeError:
      out.append({'parse_error': m.group(1)})
  return out

def check(stdout: str, contract: dict = None) -> dict:
  '''Return {'pass': bool, 'failures': [...], 'records': [...]}. Missing VERIFY line is a failure.'''
  contract = dict(DEFAULT_CONTRACT if contract is None else contract)
  recs = parse_verify(stdout)
  failures = []
  if not recs:
    failures.append('no VERIFY{} line in program output')
  for r in recs:
    if 'parse_error' in r:
      failures.append('unparseable VERIFY line: ' + r['parse_error'][:200]); continue
    if r.get('finite') is False:
      failures.append(f"{r.get('code','?')}: non-finite / blown-up solution")
    for name, val in r.get('invariants', {}).items():
      tol = contract.get(name)
      if tol is None: continue
      if not (isinstance(val, (int, float)) and abs(val) <= tol):
        failures.append(f"{r.get('code','?')}[{r.get('method','')}]: {name} = {val} exceeds tolerance {tol}")
  return {'pass': not failures, 'failures': failures, 'records': recs}

VERIFY_INSTRUCTIONS = ('The program MUST end by printing exactly one line of the form '
  'VERIFY {"code": "<name>", "method": "<integrator>", "invariants": {"<invariant>": <max relative drift>, ...}, "finite": true|false} '
  'reporting the maximum relative drift over the whole run of every conserved quantity of the model '
  '(energy, momentum, canonical toroidal momentum, charge, max |div B|, ...). Use these names when applicable: '
  + ', '.join(DEFAULT_CONTRACT) + '.')
