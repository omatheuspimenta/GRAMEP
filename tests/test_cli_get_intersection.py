# from pytest import mark
# from typer.testing import CliRunner

# from gramep.cli import app

# runner = CliRunner()
# SAVE_PATH = 'data/output/mutations/'


# def test_get_intersection_default(caplog):
#     result = runner.invoke(
#         app,
#         ['get-intersection', '--save-path', SAVE_PATH],
#     )
#     assert result.exit_code == 0
#     assert 'Done!' in caplog.text


# @mark.parametrize(
#     'intersection_seletion', ['Beta-Delta', 'Beta-Alpha', 'Beta-Delta-Alpha']
# )
# def test_get_intersection_with_selection(intersection_seletion, caplog):
#     result = runner.invoke(
#         app,
#         [
#             'get-intersection',
#             '--save-path',
#             SAVE_PATH,
#             '--intersection-seletion',
#             intersection_seletion,
#         ],
#     )
#     assert result.exit_code == 0
#     assert 'Done!' in caplog.text
