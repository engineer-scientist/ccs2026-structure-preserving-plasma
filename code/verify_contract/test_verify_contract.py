"""Offline test of the VERIFY{} contract on real structure-preserving pusher runs (no LLM needed)."""
import subprocess, json, verify_contract as vc
P = '/workspace/research/sp/pusher/pusher_v'
cases = [('boris', 'solovev', ''), ('dg', 'solovev', ''), ('rk4', 'solovev', ''), ('bsi', 'uniform', '-ts_basicsymplectic_type 2')]
contract = {'energy_rel_drift_max': 1e-6, 'p_phi_rel_drift_max': 1e-2}
rows = []
for m, f, extra in cases:
    cmd = f'{P} -field {f} -method {m} -dt 0.2 -steps 50000 -nout 50000 -u0 1.1,0,0,0,0.004,0.003 -o /tmp/vc.dat -snes_rtol 1e-14 -snes_atol 1e-15 {extra}'
    out = subprocess.run(cmd.split(), capture_output=True, text=True).stdout
    r = vc.check(out, contract)
    rows.append({'method': m, 'field': f, 'pass': r['pass'], 'failures': r['failures'], 'invariants': r['records'][0].get('invariants') if r['records'] else None})
    print(m, f, 'PASS' if r['pass'] else 'FAIL', r['failures'])
assert vc.check('no verify line here')['pass'] is False
json.dump(rows, open('/workspace/research/results/verify_contract_test.json', 'w'), indent=1)
