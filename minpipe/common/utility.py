def disk_usage(path: str = "/") -> dict:
    import shutil
    total, used, free = [i / 1073741824 for i in shutil.disk_usage(path)]
    return {"total": total, "used": used, "free": free,
            "used_p": round((used / total) * 100, 2), "free_p": round((free / total) * 100, 2)}