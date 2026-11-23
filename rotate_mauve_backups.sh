#!/bin/bash
set -euo pipefail

readonly root=/scratch/pawsey1308
readonly source_dir="$root/mauve"
readonly backup="$root/mauve_backup"
readonly older_backup="$root/mauve_backupp"
readonly index_lock='/._dbindex_/.lock'

fail() {
    printf '%s\n' "$*" >&2
    exit 1
}

dry_run=false
case "${1-}" in
    "") [[ $# -eq 0 ]] || fail "Usage: $0 [--dry-run]" ;;
    --dry-run) [[ $# -eq 1 ]] || fail "Usage: $0 [--dry-run]"
               dry_run=true ;;
    *) fail "Usage: $0 [--dry-run]" ;;
esac

[[ -d "$root" && ! -L "$root" && -w "$root" && -x "$root" ]] ||
    fail "Cannot write to project directory: $root"

if ! "$dry_run"; then
    exec 9>>"$root/.mauve-backup.lock"
    flock -n 9 || fail "Could not acquire the MAUVE backup lock."
fi

[[ -d "$source_dir" && ! -L "$source_dir" && -r "$source_dir" && -x "$source_dir" ]] ||
    fail "Source is not a readable directory: $source_dir"
for directory in "$backup" "$older_backup"; do
    [[ ! -L "$directory" && ( ! -e "$directory" || -d "$directory" ) ]] ||
        fail "Unexpected backup path: $directory"
done

if "$dry_run"; then
    printf 'Source: %s\nLatest: %s\nPrevious: %s\n' "$source_dir" "$backup" "$older_backup"
    printf 'Excluded: %s%s\n' "$source_dir" "$index_lock"
    if [[ -d "$backup" ]]; then
        printf 'Would stage the source, then replace the previous backup with the latest and publish the new copy.\n'
    else
        printf 'Would preserve any previous backup and copy the source to recreate latest.\n'
    fi
    exit 0
fi

printf '%s Starting MAUVE backup rotation\n' "$(date -Is)"
staging=$(mktemp -d -- "$root/.mauve-backup.partial.XXXXXXXX")
trap 'rm -rf --one-file-system -- "$staging"' EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

# Finish copying before changing either backup generation.
rsync -rlD --exclude="$index_lock" -- "$source_dir/" "$staging/mauve/"

# If a previous copy was interrupted, keep the last good generation.
if [[ -d "$backup" ]]; then
    rm -rf --one-file-system -- "$older_backup"
    mv -T -- "$backup" "$older_backup"
fi
mv -T -- "$staging/mauve" "$backup"
printf '%s Finished MAUVE backup rotation\n' "$(date -Is)"
