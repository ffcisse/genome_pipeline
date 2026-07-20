#!/usr/bin/env bash
#
# submit_phase5.sh -- submit run_phase5.sbatch with account/QOS/partition/
# mail-user read from config/config.yaml's `slurm:` block (the same block
# submit_phase1.sh/submit_phase2a.sh/submit_phase2b.sh/submit_phase3.sh/
# submit_phase4.sh read), so switching clusters or allocations is a config
# edit, not a script edit.
#
# Usage: workflow/scripts/submit_phase5.sh

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO_ROOT}"

read -r ACCOUNT QOS PARTITION MAIL_USER <<EOF
$(python3 -c "
import yaml
with open('config/config.yaml') as fh:
    cfg = yaml.safe_load(fh)
slurm = cfg.get('slurm', {})
def field(v):
    return v if v else '-'
print(
    field(slurm.get('account')),
    field(slurm.get('qos')),
    field(slurm.get('partition')),
    field(slurm.get('mail_user')),
)
")
EOF

SBATCH_ARGS=()
[ "${ACCOUNT}" != "-" ] && SBATCH_ARGS+=(--account="${ACCOUNT}")
[ "${QOS}" != "-" ] && SBATCH_ARGS+=(--qos="${QOS}")
[ "${PARTITION}" != "-" ] && SBATCH_ARGS+=(--partition="${PARTITION}")
[ "${MAIL_USER}" != "-" ] && SBATCH_ARGS+=(--mail-user="${MAIL_USER}")

echo "Submitting with: sbatch ${SBATCH_ARGS[*]} workflow/scripts/run_phase5.sbatch"
sbatch "${SBATCH_ARGS[@]}" workflow/scripts/run_phase5.sbatch
