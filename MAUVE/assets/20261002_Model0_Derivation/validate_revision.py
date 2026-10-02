"""Targeted falsification of new algebra and final artifact structure."""
from pathlib import Path
import json
import re
import hashlib
import numpy as np
from scipy.integrate import solve_ivp

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent.parent
MD=ROOT/'20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md'
text=MD.read_text()
errors=[]
for gamma,tau,tauphi in [(.3,.003,1.),(2.,.5,.3)]:
    t=np.linspace(0,1,501)
    expected=(1+t/tauphi)*np.exp(-gamma*t)
    sol=solve_ivp(lambda u,y: [((1+u/tauphi)*np.exp(-gamma*u)-y[0])/tau],
        (0,1),[1.],t_eval=t,rtol=1e-11,atol=1e-12,max_step=.001)
    difference=1/tau-gamma
    if abs(difference)<1e-12:
        analytic=np.exp(-t/tau)*(1+t/tau+t*t/(2*tau*tauphi))
    else:
        # Algebraic stable form of the displayed two elementary integrals.
        analytic=(np.exp(-t/tau)+(np.exp(-gamma*t)-np.exp(-t/tau))/(tau*difference)
            +(np.exp(-gamma*t)*(difference*t-1)+np.exp(-t/tau))/(tau*tauphi*difference**2))
    error=float(np.max(np.abs(analytic-sol.y[0])))
    assert sol.success and error<1e-9
    errors.append(error)
numbers=list(map(int,re.findall(r'\\tag\{(\d+)\}',text)))
assert numbers==list(range(1,len(numbers)+1))
assert '[eq:' not in text and '(t)(t' not in text
assert r'\widetilde\Sigma_{\mathrm{HI}}(t)' not in text
main=text.split('# Appendix A.')[0]
inverse=main[main.index('## 8.6'):main.index('## 8.7')]
assert r'\mathcal L_{B,' not in inverse
assert 'two rates are equal' not in main
assert '1.5' in main and 'negative temporal derivative' in main
anchors=set(re.findall(r'<span id="([^"]+)"',text))
assert set(re.findall(r'\]\(#(ref-[^)]+)\)',text))<=anchors
for path in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',text):
    assert (ROOT/path).exists()
audit=json.loads((OUT/'numerical_audit.json').read_text())
assert audit['line_budget'][2]['r_holmes_required_at_fiducial']<0
assert audit['minimum_decline_time_gyr']>3.9
assert 0<audit['response_relative_correction_at_1gyr']<.001
result=dict(equal_rate_kernel_max_abs=errors,numbered_equations=len(numbers),
    reference_anchors=len(anchors),md_sha256=hashlib.sha256(MD.read_bytes()).hexdigest(),
    old_md_sha256=hashlib.sha256((ROOT/'20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md').read_bytes()).hexdigest(),
    old_pdf_sha256=hashlib.sha256((ROOT/'20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.pdf').read_bytes()).hexdigest(),
    structural_checks='passed',classification='checks of stated equations; no new data fit')
(OUT/'revision_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
