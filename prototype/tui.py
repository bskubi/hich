from textual.app import App, ComposeResult
from textual.widgets import DirectoryTree, Button, Input, Static
from textual.validation import ValidationResult, Validator, Function
from textual.containers import ScrollableContainer
import pathvalidate

def is_valid_filename(filename):
    try:
        pathvalidate.validate_filename(filename)
        return True
    except pathvalidate.ValidationError:
        return False

class FileSaver(Static):
    def compose(self) -> ComposeResult:
        yield Input(placeholder = "nextflow.config", validators = [Function(is_valid_filename)], id = "filename")
        yield Button("Save New File", id = "save")
        yield Button("Overwrite", id = "overwrite")

class FilePicker(App):
    CSS_PATH = "tui.tcss"

    def compose(self) -> ComposeResult:
        yield DirectoryTree("./")
        yield ScrollableContainer(FileSaver())



if __name__ == "__main__":
    app = FilePicker()
    app.run()