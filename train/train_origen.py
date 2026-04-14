from pathlib import Path
import runpy


if __name__ == "__main__":
    root_script = Path(__file__).resolve().parents[1] / "10-train_origen.py"
    runpy.run_path(str(root_script), run_name="__main__")
