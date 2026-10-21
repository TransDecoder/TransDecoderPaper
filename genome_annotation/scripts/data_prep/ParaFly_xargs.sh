#!/bin/bash
# Drop-in replacement for ParaFly (the bioconda 2.15.2 build runs commands serially regardless of -CPU).
# Supports -c <cmds> -CPU <n> [-shuffle] [-v] [-failed_cmds <file>] [-max_retry <n>].
# Successful commands are appended to <cmds>.completed and skipped when rerun, as ParaFly does.
set -uo pipefail
cmds= ; cpu=1 ; shuffle=0 ; failed= ; retry=1
while [ $# -gt 0 ]; do
  case "$1" in
    -c) cmds=$2; shift 2;;
    -CPU) cpu=$2; shift 2;;
    -shuffle) shuffle=1; shift;;
    -v|-vv) shift;;
    -failed_cmds) failed=$2; shift 2;;
    -max_retry) retry=$2; shift 2;;
    *) echo "ParaFly wrapper: unknown option $1" >&2; exit 2;;
  esac
done
[ -n "$cmds" ] && [ -f "$cmds" ] || { echo "ParaFly wrapper: need -c <cmds file>" >&2; exit 2; }
completed="$cmds.completed"; touch "$completed"
failed=${failed:-FailedCommands}; : > "$failed"
todo=$(mktemp "$cmds.todo.XXXXXX")
grep -v '^[[:space:]]*$' "$cmds" | grep -Fxv -f "$completed" > "$todo.all" || true
if [ $shuffle = 1 ]; then shuf "$todo.all" > "$todo"; else mv "$todo.all" "$todo"; fi
rm -f "$todo.all"
echo "Number of Commands: $(grep -c . "$cmds"); to run: $(grep -c . "$todo")" >&2
export completed failed retry
run_one() {
  local i
  for ((i=1; i<=retry; i++)); do
    if bash -c "$1"; then
      { flock 9; printf '%s\n' "$1" >&9; } 9>>"$completed"; return 0
    fi
  done
  { flock 9; printf '%s\n' "$1" >&9; } 9>>"$failed"; return 0
}
export -f run_one
tr '\n' '\0' < "$todo" | xargs -0 -P "$cpu" -I CMD bash -c 'run_one "$1"' _ CMD
rm -f "$todo"
nfail=$(grep -c . "$failed")
if [ "$nfail" -gt 0 ]; then
  echo "ParaFly wrapper: $nfail commands failed; see $failed" >&2; exit 1
fi
echo "All commands completed successfully. :-)" >&2
