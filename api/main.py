from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
from inputs.vents import table  # noqa: E402  (project import after sys.path patch)

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

GRIB_DIR = ROOT / "data" / "grib_vent"
_df_cache: dict = {}


def _load(filename: str):
    if filename not in _df_cache:
        df = table(str(GRIB_DIR / filename))
        # Normalise longitude 0-360 → -180/180
        df["longitude"] = df["longitude"].apply(lambda x: x - 360 if x > 180 else x)
        _df_cache[filename] = df
    return _df_cache[filename]


@app.get("/files")
def list_files():
    return sorted(f.name for f in GRIB_DIR.glob("*.grb2"))


@app.get("/wind/{filename}/meta")
def get_meta(filename: str):
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    df = _load(filename)
    step_df = (
        df[["temps", "valid_time"]]
        .drop_duplicates("temps")
        .sort_values("temps")
    )
    times = step_df["temps"].tolist()
    # Ensure UTC ISO strings for every step
    def _iso(ts):
        t = pd.Timestamp(ts)
        return (t if t.tzinfo else t.tz_localize("UTC")).isoformat()
    valid_times = [_iso(vt) for vt in step_df["valid_time"]]
    return {
        "times": times,
        "valid_times": valid_times,
        "bbox": [
            float(df["longitude"].min()),
            float(df["latitude"].min()),
            float(df["longitude"].max()),
            float(df["latitude"].max()),
        ],
    }


@app.get("/wind/{filename}/step/{time_idx}")
def get_step(filename: str, time_idx: int, stride: int = 2):
    if not (GRIB_DIR / filename).exists():
        raise HTTPException(404, "File not found")
    df = _load(filename)
    times = sorted(df["temps"].unique())
    if not (0 <= time_idx < len(times)):
        raise HTTPException(400, "Invalid time index")

    t = times[time_idx]
    df_t = df[df["temps"] == t].dropna(subset=["u10", "v10"])

    lats = sorted(df_t["latitude"].unique())
    lons = sorted(df_t["longitude"].unique())
    lats_k = set(lats[::stride])
    lons_k = set(lons[::stride])
    df_t = df_t[df_t["latitude"].isin(lats_k) & df_t["longitude"].isin(lons_k)]

    cols = df_t[["latitude", "longitude", "force", "direction"]].copy()
    cols.columns = ["lat", "lon", "speed", "dir"]
    cols = cols.round({"lat": 4, "lon": 4, "speed": 2, "dir": 1})

    return {"time_h": float(t), "data": cols.to_dict(orient="records")}
