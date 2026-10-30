from pytest import mark
from typer.testing import CliRunner

from gramep.cli import app

runner = CliRunner()


@mark.parametrize(
    'objective', ['get-mutations', 'get-intersection', 'classify', 'predict']
)
def test_gramep_from_configfile(objective, caplog):
    result = runner.invoke(
        app,
        [
            'from-configfile',
            '--objective',
            objective,
            '--config-file',
            'data/parameters.ini',
        ],
    )
    assert result.exit_code == 0
    assert 'Done!' in caplog.text
