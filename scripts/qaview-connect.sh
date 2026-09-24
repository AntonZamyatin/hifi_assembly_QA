#!/usr/bin/env bash
# Start qaview on a server and reach it from this laptop through an ssh tunnel.
#
#   qaview-connect.sh <host> [--port 5010] [--dir /scratch/projects/assemblyQA]
#                            [--project path/to/project.toml] [--open]
#
# <host> is anything ssh accepts: a Host from ~/.ssh/config, user@name or an IP.
# One ssh session forwards localhost:<port> here to localhost:<port> on the
# server and runs qaview there; the link is printed as soon as qaview answers.
# Ctrl-C stops qaview on the server, background recruitment and scoring
# included; so does losing the connection.
#
#   --port N      port on both ends (default 5010)
#   --dir PATH    the assemblyQA checkout on the server, with its venv
#                 (default /scratch/projects/assemblyQA)
#   --project F   open this project file (a path on the server); without it the
#                 page starts on a form to choose one
#   --open        also open the link in the default browser

set -euo pipefail

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit "${1:-0}"; }

HOST=""
PORT=5010
DIR=/scratch/projects/assemblyQA
PROJECT=""
OPEN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage 0 ;;
        --port)    PORT="${2:?--port needs a number}"; shift 2 ;;
        --dir)     DIR="${2:?--dir needs a path}"; shift 2 ;;
        --project) PROJECT="${2:?--project needs a path}"; shift 2 ;;
        --open)    OPEN=1; shift ;;
        -*)        echo "unknown option: $1" >&2; usage 2 ;;
        *)         if [[ -z "$HOST" ]]; then HOST="$1"; shift; else echo "one host only" >&2; usage 2; fi ;;
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
# POSIX single quotes for the outer layer: ssh hands the command to the login
# shell on the server, which may be sh (dash), and dash cannot read bash's $'...'.
sq() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }
REMOTE="cd $(q "$DIR") || exit 1
[ -x venv/bin/python ] || { echo \"no venv/bin/python in $(q "$DIR")\" >&2; exit 1; }
S=; command -v setsid >/dev/null && S=setsid
\$S venv/bin/python -m qaview --port $(q "$PORT")${PROJECT:+ $(q "$PROJECT")} &
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
