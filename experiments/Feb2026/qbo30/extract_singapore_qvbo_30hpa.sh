#!/usr/bin/env bash
# Extract daily Singapore/Changi QVBO proxy: zonal wind interpolated to 30 hPa.
#
# Example:
#   ./extract_singapore_qvbo_30hpa.sh --start 2006-01-01 --end 2020-12-31

set -euo pipefail

START_DATE="2006-01-01"
END_DATE="$(date -u +%F)"
OUTPUT="singapore_qvbo_30hpa_daily.csv"
URL="https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/SNM00048698-data.txt.zip"

usage() {
    cat <<'EOF'
Usage: extract_singapore_qvbo_3hpa.sh [options]

Download the NOAA IGRA Singapore/Changi (WMO 48698) radiosonde archive and
write daily mean zonal wind at 30 hPa. Wind is interpolated linearly in log
pressure within each individual sounding before daily averaging.

Options:
  --start YYYY-MM-DD   First UTC day to retain (default: 2006-01-01)
  --end YYYY-MM-DD     Last UTC day to retain (default: today)
  --output FILE        CSV output path (default: singapore_qvbo_30hpa_daily.csv)
  --help               Show this message
EOF
}

while (($#)); do
    case "$1" in
        --start)
            START_DATE="$2"
            shift 2
            ;;
        --end)
            END_DATE="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! "$START_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ||
      ! "$END_DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    printf 'Dates must use YYYY-MM-DD.\n' >&2
    exit 2
fi
if [[ "$START_DATE" > "$END_DATE" ]]; then
    printf '--start must not be later than --end.\n' >&2
    exit 2
fi

temporary_directory="$(mktemp -d)"
trap 'rm -rf "$temporary_directory"' EXIT
archive="$temporary_directory/SNM00048698-data.txt.zip"

curl --fail --location --retry 3 --retry-delay 2 --output "$archive" "$URL"

ARCHIVE="$archive" START_DATE="$START_DATE" END_DATE="$END_DATE" OUTPUT="$OUTPUT" \
python3 - <<'PY'
import csv
import math
import os
import zipfile
from collections import defaultdict
from datetime import date

target_pressure_pa = 3000  # 30 hPa
start_date = date.fromisoformat(os.environ["START_DATE"])
end_date = date.fromisoformat(os.environ["END_DATE"])


def integer_field(line, start, stop):
    """Decode an IGRA fixed-width integer field and its missing-value marker."""
    value = line[start:stop].strip().rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    try:
        parsed = int(value)
    except ValueError:
        return None
    return None if parsed <= -8888 else parsed


def zonal_wind_at_30hpa(level_lines):
    """Interpolate u linearly against log pressure for one sounding."""
    levels = []
    for line in level_lines:
        pressure_pa = integer_field(line, 9, 15)
        direction_degrees = integer_field(line, 40, 45)
        speed_tenths_ms = integer_field(line, 46, 51)
        if (
            pressure_pa is None
            or direction_degrees is None
            or speed_tenths_ms is None
            or pressure_pa <= 0
        ):
            continue
        # Meteorological direction is where wind comes from; u is eastward.
        zonal_wind = -speed_tenths_ms / 10.0 * math.sin(
            math.radians(direction_degrees)
        )
        levels.append((pressure_pa, zonal_wind))

    levels.sort(reverse=True)
    for (upper_pressure, upper_u), (lower_pressure, lower_u) in zip(
        levels, levels[1:]
    ):
        if lower_pressure <= target_pressure_pa <= upper_pressure:
            log_range = math.log(upper_pressure / lower_pressure)
            # Reject brackets wider than approximately a factor of 2.7 in
            # pressure rather than extrapolating a sparse profile.
            if log_range > 1.0:
                continue
            fraction = math.log(upper_pressure / target_pressure_pa) / log_range
            return upper_u + fraction * (lower_u - upper_u)
    return None


daily_values = defaultdict(list)
with zipfile.ZipFile(os.environ["ARCHIVE"]) as archive:
    with archive.open(archive.namelist()[0]) as source:
        header_date = None
        level_lines = []

        def finish_sounding():
            if header_date is None:
                return
            value = zonal_wind_at_30hpa(level_lines)
            if value is not None:
                daily_values[header_date].append(value)

        for raw_line in source:
            line = raw_line.decode("ascii", "replace")
            if line.startswith("#"):
                finish_sounding()
                level_lines = []
                try:
                    header_date = date(
                        int(line[13:17]), int(line[18:20]), int(line[21:23])
                    )
                except ValueError:
                    header_date = None
            elif header_date is not None:
                level_lines.append(line)
        finish_sounding()

with open(os.environ["OUTPUT"], "w", newline="", encoding="utf-8") as destination:
    writer = csv.writer(destination)
    writer.writerow(("date", "u_30hpa_ms", "soundings"))
    for day in sorted(daily_values):
        if start_date <= day <= end_date:
            values = daily_values[day]
            writer.writerow((day.isoformat(), f"{sum(values) / len(values):.6f}", len(values)))
PY

printf 'Wrote %s\n' "$OUTPUT"
