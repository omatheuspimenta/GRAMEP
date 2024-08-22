from gramep.cli import app
from pytest import mark
from typer.testing import CliRunner

runner = CliRunner()

REF_PATH = 'data/reference/SARS-CoV2_wuhan_refseq.fasta'
SEQ_PATH = 'data/VOCs/'


@mark.parametrize(
    'variant',
    [
        'Alpha.fasta',
        'Beta.fasta',
        'Delta.fasta',
        'Gamma.fasta',
    ],
)
def test_grid_search(variant, caplog):
    result = runner.invoke(
        app,
        [
            'grid-search',
            '--rpath',
            REF_PATH,
            '--spath',
            SEQ_PATH + variant,
            '-minw',
            '5',
            '-maxw',
            '10',
            '-mins',
            '1',
            '-maxs',
            '3',
        ],
    )
    assert result.exit_code == 0
    assert 'Done!' in caplog.text
