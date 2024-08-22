from gramep.cli import app
from pytest import mark
from typer.testing import CliRunner

runner = CliRunner()
SAVE_PATH = 'data/output/mutations/'
DIR_PATH = 'data/VOCs/'


def test_cli_classify(caplog):
    result = runner.invoke(
        app,
        [
            'classify',
            '-w',
            '10',
            '-s',
            '1',
            '--save-path',
            SAVE_PATH,
            '--dir-path',
            DIR_PATH,
            '--dictonary',
            'ACTG',
        ],
    )
    assert result.exit_code == 0
    assert 'Done!' in caplog.text


@mark.parametrize(
    'variant',
    [
        'Alpha.fasta',
        'Beta.fasta',
        'Delta.fasta',
        'Gamma.fasta',
    ],
)
def test_cli_predict(variant, caplog):
    result = runner.invoke(
        app,
        [
            'predict',
            '-w',
            '10',
            '-s',
            '1',
            '--save-path',
            'data/output/mutations/',
            '-pseqpath',
            DIR_PATH + variant,
            '-dpath',
            'data/output/classify/',
            '-d',
            'ACTG',
            '-lrpath',
            'data/output/classify/' + 'model/ranges.sav',
            '-lmpath',
            'data/output/classify/' + 'model/model.sav',
        ],
    )
    assert result.exit_code == 0
    assert 'Done!' in caplog.text
