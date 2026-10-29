#!/usr/bin/env bash
# Start qaview on a server and reach it from this laptop through an ssh tunnel.
#
#   qaview-connect.sh <host> [--port 5010] [--dir ~/hifi_assembly_QA]
#                            [--project path/to/project.toml] [--python PATH] [--open]
#
# Runs on any laptop with bash, ssh and curl (macOS, Linux, WSL or Git Bash on
# Windows); nothing of qaview needs to be installed here. Get it alone with
#   curl -O https://raw.githubusercontent.com/AntonZamyatin/hifi_assembly_QA/master/scripts/qaview-connect.sh
#
# <host> is anything ssh accepts: a Host from ~/.ssh/config, user@name or an IP.
# One ssh session forwards localhost:<port> here to localhost:<port> on the
# server and runs qaview there; the link is printed as soon as qaview answers.
# Ctrl-C stops qaview on the server, background recruitment and scoring
# included; so does losing the connection.
#
#   --port N      port on both ends (default 5010); users sharing a server need
#                 one each
#   --dir PATH    your clone of the repository on the server (default
#                 ~/hifi_assembly_QA; ~ is the server's home)
#   --project F   open this project file (a path on the server, relative to
#                 --dir or absolute); without it the page asks for one
#   --python P    the server's Python with qaview installed; by default the
#                 clone's venv/bin/python, else the conda env "assemblyqa"
#   --open        also open the link in the default browser
#
# Defaults can be kept in the environment instead: QAVIEW_HOST, QAVIEW_PORT,
# QAVIEW_DIR, QAVIEW_PYTHON (e.g. in ~/.bashrc or ~/.zshrc).

set -euo pipefail

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit "${1:-0}"; }

HOST="${QAVIEW_HOST:-}"
PORT="${QAVIEW_PORT:-5010}"
DIR="${QAVIEW_DIR:-~/hifi_assembly_QA}"
PROJECT=""
PYTHON="${QAVIEW_PYTHON:-}"
OPEN=0
HOST_GIVEN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage 0 ;;
        --port)    PORT="${2:?--port needs a number}"; shift 2 ;;
        --dir)     DIR="${2:?--dir needs a path}"; shift 2 ;;
        --project) PROJECT="${2:?--project needs a path}"; shift 2 ;;
        --python)  PYTHON="${2:?--python needs a path}"; shift 2 ;;
        --open)    OPEN=1; shift ;;
        -*)        echo "unknown option: $1" >&2; usage 2 ;;
        *)         if [[ "$HOST_GIVEN" == 0 ]]; then HOST="$1"; HOST_GIVEN=1; shift; else echo "one host only" >&2; usage 2; fi ;;
    esac
done
[[ -n "$HOST" ]] || usage 2
[[ "$PORT" =~ ^[0-9]+$ ]] || { echo "--port must be a number: $PORT" >&2; exit 2; }

URL="http://localhost:${PORT}/"

# The port here must be free, or ssh would forward it somewhere else silently.
if (exec 3<>"/dev/tcp/127.0.0.1/${PORT}") 2>/dev/null; then
    echo "port ${PORT} is already in use on this laptop; pick another with --port" >&2
    exit 1
fi

# What runs on the server. qaview runs in a session of its own (setsid), so the
# tty's Ctrl-C reaches only this wrapper, which passes exactly one TERM on --
# for Ctrl-C (INT), a dropped connection (HUP) or TERM alike. One signal is a
# clean shutdown that also stops recruitment and scoring; a second one would
# make uvicorn force-exit and skip that, leaving the workers running.
q() { printf '%q' "$1"; }                    # quoting for bash, inside REMOTE
# The same, but a leading ~ stays the server's home.
qh() { if [[ "$1" == "~" || "$1" == "~/"* ]]; then printf '~/%s' "$(q "${1#\~/}")"; else q "$1"; fi; }
# POSIX single quotes for the outer layer: ssh hands the command to the login
# shell on the server, which may be sh (dash), and dash cannot read bash's $'...'.
sq() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }
# The Python: --python, else the clone's venv, else the conda env "assemblyqa"
# (environment.yml), looked up through conda when it is on PATH and in the usual
# install folders when it is not, as in a non-interactive ssh shell.
REMOTE="cd $(qh "$DIR") || { echo 'no clone of the repository at ${DIR//\'/} on the server; give --dir' >&2; exit 1; }
PY=${PYTHON:+$(qh "$PYTHON")}
[ -z \"\$PY\" ] && [ -x venv/bin/python ] && PY=venv/bin/python
if [ -z \"\$PY\" ]; then
    for c in \"\${CONDA_EXE:-}\" \$(command -v mamba conda 2>/dev/null); do
        b=\$(\"\$c\" info --base 2>/dev/null) && [ -x \"\$b/envs/assemblyqa/bin/python\" ] && { PY=\$b/envs/assemblyqa/bin/python; break; }
    done
fi
if [ -z \"\$PY\" ]; then
    for p in ~/*/envs/assemblyqa/bin/python /opt/*/envs/assemblyqa/bin/python; do [ -x \"\$p\" ] && { PY=\$p; break; }; done
fi
[ -n \"\$PY\" ] || { echo 'no Python with qaview found: make the env (mamba env create -f environment.yml) or give --python' >&2; exit 1; }
echo \"qaview from \$PY\"
S=; command -v setsid >/dev/null && S=setsid
\$S \"\$PY\" -m qaview --port $(q "$PORT")${PROJECT:+ $(qh "$PROJECT")} &
pid=\$!
trap 'kill -TERM \$pid 2>/dev/null' INT HUP TERM
wait \$pid; wait \$pid"

# Announce the link once qaview answers through the tunnel.
announce() {
    for _ in $(seq 1 120); do
        if curl -s -o /dev/null --max-time 2 "${URL}api/project"; then
            printf '\n  qaview is ready:  %s\n  (Ctrl-C here stops it)\n\n' "$URL"
            if [[ "$OPEN" == 1 ]]; then
                if command -v open >/dev/null; then open "$URL"
                elif command -v xdg-open >/dev/null; then xdg-open "$URL" >/dev/null 2>&1
                fi
            fi
            return
        fi
        sleep 1
    done
    echo "qaview did not answer on ${URL} within 2 minutes; see the messages above" >&2
}
announce &
WATCH=$!
trap 'kill $WATCH 2>/dev/null || true' EXIT

echo "connecting to ${HOST}: starting qaview in ${DIR} on port ${PORT} ..."
# -t: a tty, so Ctrl-C reaches the server; ExitOnForwardFailure: stop at once
# if the tunnel cannot be made rather than run qaview unreachable.
ssh -t -o ExitOnForwardFailure=yes -o ServerAliveInterval=30 \
    -L "${PORT}:localhost:${PORT}" "$HOST" "bash -c $(sq "$REMOTE")"
