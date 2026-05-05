import re

def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.#=-]+", "_", name)
