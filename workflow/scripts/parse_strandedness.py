#!/usr/bin/env python3
"""
Parse RSeQC infer_experiment.py output and print the featureCounts strandedness code.

Output (stdout): single integer
  0 = unstranded
  1 = forward-stranded (e.g. Takara SMARTer)
  2 = reverse-stranded (e.g. TruSeq Stranded mRNA, NEBNext Directional)

Decision rule: >60% of assignable reads in one direction → stranded.
Usage: python parse_strandedness.py infer_experiment_output.txt
"""
import sys

fwd = rev = 0.0

with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.strip()
        # Paired-end patterns
        if "1++,1--,2+-,2-+" in line or '"++,--"' in line:
            fwd = float(line.split(":")[-1].strip())
        elif "1+-,1-+,2++,2--" in line or '"+-,-+"' in line:
            rev = float(line.split(":")[-1].strip())

THRESHOLD = 0.60
if rev > THRESHOLD:
    code = 2
elif fwd > THRESHOLD:
    code = 1
else:
    code = 0

sys.stderr.write(
    f"infer_experiment: fwd={fwd:.3f}  rev={rev:.3f}  threshold={THRESHOLD}  "
    f"→ strandedness={code} "
    f"({'reverse' if code == 2 else 'forward' if code == 1 else 'unstranded'})\n"
)
print(code)
