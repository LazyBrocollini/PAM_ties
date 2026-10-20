"""Import smoke tests only: Streamlit itself is stubbed, not UI-tested."""
import runpy
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def _check_import(path):
    import tempfile
    code = """
import runpy, sys, types
st = types.ModuleType('streamlit')
st.set_page_config = lambda **kwargs: None
st.cache_data = lambda **kwargs: (lambda function: function)
sys.modules['streamlit'] = st
namespace = runpy.run_path(sys.argv[1], run_name='deployment_import_check')
assert namespace['normalize_reference']('>ref\\nACGT')[0] == 'ACGT'
print('Import smoke check passed')
"""
    with tempfile.TemporaryDirectory() as directory:
        result = subprocess.run([sys.executable, '-c', code, str(path)], cwd=directory, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr


def test_modular_import_from_another_directory():
    _check_import(ROOT / 'app.py')


def test_standalone_has_no_local_package_dependencies():
    import tempfile
    text = (ROOT / 'standalone' / 'app.py').read_text()
    assert 'from src.' not in text
    assert 'from .sequence_utils' not in text
    with tempfile.TemporaryDirectory() as directory:
        standalone = Path(directory) / 'app.py'
        standalone.write_text(text)
        _check_import(standalone)
