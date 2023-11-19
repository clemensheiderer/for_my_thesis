from ete3 import Tree
from app_helper import modify_tree_branch_length, subtree_to_fasta, fasta_copy
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QWidget, QApplication, QLabel, QTextEdit, QPushButton
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import QFileInfo
import os


class NewickFileWindow(QWidget):
    def __init__(self, label_text):
        super().__init__()
        self.setWindowTitle("Subtree Creator - Newick File Selection")
        self.setGeometry(100, 100, 400, 200)

        self.label = QLabel(label_text)
        self.text_edit = QTextEdit()
        self.text_edit.setAcceptDrops(True)

        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_file)

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.text_edit)
        layout.addWidget(self.browse_button)

        self.setLayout(layout)

        self.text_edit.dragEnterEvent = self.drag_enter_event
        self.text_edit.dropEvent = self.drop_event

    def drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and url.toString().lower().endswith(
                    (".txt", ".newick", ".nwk", ".fasta.treefile", "phy.treefile", ".node.treefile", ".nex.treefile"))
                   for url in urls):
                event.acceptProposedAction()

    def drop_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and url.toString().lower().endswith(
                    (".txt", ".newick", ".nwk", ".fasta.treefile", "phy.treefile", ".node.treefile", ".nex.treefile"))
                   for url in urls):
                file_path = urls[0].toLocalFile()
                self.text_edit.setText(file_path)
                event.acceptProposedAction()

    def browse_file(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            filter="All Files (*.txt *.newick *.nwk *.fasta.treefile *phy.treefile *.treefile *node.treefile)")
        if file_path:
            self.text_edit.setText(file_path)


class FastaFileWindow(QWidget):
    def __init__(self, label_text):
        super().__init__()
        self.setWindowTitle("Subtree Creator - Fasta File Selection")
        self.setGeometry(100, 100, 400, 200)

        self.label = QLabel(label_text)
        self.text_edit = QTextEdit()
        self.text_edit.setAcceptDrops(True)

        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_file)

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.text_edit)
        layout.addWidget(self.browse_button)

        self.setLayout(layout)

        self.text_edit.dragEnterEvent = self.drag_enter_event
        self.text_edit.dropEvent = self.drop_event

    def drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and url.toString().lower().endswith((".fa", ".fasta")) for url in urls):
                event.acceptProposedAction()

    def drop_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and url.toString().lower().endswith((".fa", ".fasta")) for url in urls):
                file_path = urls[0].toLocalFile()
                self.text_edit.setText(file_path)
                event.acceptProposedAction()

    def browse_file(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(filter="FASTA Files (*.fa *.fasta)")
        if file_path:
            self.text_edit.setText(file_path)


class DirectoryWindow(QWidget):
    def __init__(self, label_text):
        super().__init__()
        self.setWindowTitle("Subtree Creator - Directory Selection")
        self.setGeometry(100, 100, 400, 200)

        self.label = QLabel(label_text)
        self.text_edit = QTextEdit()
        self.text_edit.setAcceptDrops(True)

        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_directory)

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.text_edit)
        layout.addWidget(self.browse_button)

        self.setLayout(layout)

        self.text_edit.dragEnterEvent = self.drag_enter_event
        self.text_edit.dropEvent = self.drop_event

    def drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and QFileInfo(url.toLocalFile()).isDir() for url in urls):
                event.acceptProposedAction()

    def drop_event(self, event):
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.isLocalFile() and QFileInfo(url.toLocalFile()).isDir() for url in urls):
                directory_path = urls[0].toLocalFile()
                self.text_edit.setText(directory_path)
                event.acceptProposedAction()

    def browse_directory(self):
        directory_dialog = QFileDialog()
        directory_dialog.setFileMode(QFileDialog.Directory)
        directory_dialog.setOption(QFileDialog.ShowDirsOnly)
        directory_path = directory_dialog.getExistingDirectory()
        self.text_edit.setText(directory_path)


class NodeWindow(QWidget):
    def __init__(self, node_name):
        super().__init__()
        self.setWindowTitle(f"Node: {node_name}")
        self.setGeometry(100, 100, 400, 200)

        self.label = QLabel(f"{node_name}:")
        self.text_edit = QTextEdit()

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.text_edit)

        self.setLayout(layout)


class SubtreeCreator(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Subtree Creator")
        self.setGeometry(100, 100, 400, 200)
        self.subtree_widgets = []
        self.node_values = []
        self.original_height = self.height()

        self.file_window_newick = NewickFileWindow("Drag and Drop a tree file here:")
        self.file_window_fasta = FastaFileWindow("Drag and Drop an alignment in fasta format here:")
        self.directory_window_output = DirectoryWindow("Drag and Drop an Output Directory here:")

        self.button_create_another_subtree = QPushButton("Create another subtree")
        self.button_create_another_subtree.clicked.connect(self.create_another_subtree)

        self.button_create_newick_subtrees = QPushButton("Create subtrees")
        self.button_create_newick_subtrees.clicked.connect(self.create_newick_subtrees)

        self.button_create_nwk_fasta_subtrees = QPushButton("Create subtrees and alignment in fasta format")
        self.button_create_nwk_fasta_subtrees.clicked.connect(self.create_nwk_fasta_subtrees)

        self.label_node_a = QLabel("NodeA:")
        self.text_edit_node_a = QTextEdit()

        self.label_node_b = QLabel("NodeB:")
        self.text_edit_node_b = QTextEdit()

        labels_layout = QHBoxLayout()
        labels_layout.addWidget(self.label_node_a)
        labels_layout.addWidget(self.label_node_b)

        nodes_layout = QHBoxLayout()
        nodes_layout.addWidget(self.text_edit_node_a)
        nodes_layout.addWidget(self.text_edit_node_b)

        buttons_layout = QVBoxLayout()
        buttons_layout.addWidget(self.button_create_another_subtree)
        buttons_layout.addWidget(self.button_create_newick_subtrees)
        buttons_layout.addWidget(self.button_create_nwk_fasta_subtrees)

        layout = QVBoxLayout()

        image_label = QLabel(self)
        pixmap = QPixmap('/home/clemens/Downloads/cut-mistle.png')
        image_label.setPixmap(pixmap)
        image_label.setAlignment(Qt.AlignCenter)

        layout.addWidget(image_label)
        layout.addSpacing(10)  # Add spacing between the image and subsequent widgets

        layout.addWidget(self.file_window_newick)
        layout.addWidget(self.file_window_fasta)
        layout.addWidget(self.directory_window_output)
        layout.addLayout(buttons_layout)
        layout.addLayout(labels_layout)
        layout.addLayout(nodes_layout)

        self.setLayout(layout)

        # self.setWindowFlags(self.windowFlags() | Qt.WindowMinimizeButtonHint | Qt.WindowResizable)
        # self.setWindowFlags(self.windowFlags() | Qt.WindowMinimizeButtonHint | Qt.WindowFlag(Qt.WindowResizable))
        self.show()

    def showEvent(self, event):
        super().showEvent(event)
        self.setMinimumHeight(self.original_height)  # Set the minimum height for the window after it's shown

    def create_another_subtree(self):
        node_a = self.text_edit_node_a.toPlainText()
        node_b = self.text_edit_node_b.toPlainText()

        widget = QWidget()
        widget.setWindowTitle("New Subtree")
        widget.setGeometry(100, 100, 400, 200)

        node_window_a = NodeWindow("NodeA")
        node_window_b = NodeWindow("NodeB")

        nodes_layout = QHBoxLayout()
        nodes_layout.addWidget(node_window_a)
        nodes_layout.addWidget(node_window_b)

        widget.setLayout(nodes_layout)
        widget.show()
        self.subtree_widgets.append(widget)
        self.node_values.append((node_a, node_b))  # Store NodeA and NodeB values

        widget.node_window_a = node_window_a  # Assign NodeA window to the widget
        widget.node_window_b = node_window_b  # Assign NodeB window to the widget

    def read_newick_file(self, file_path):
        try:
            with open(file_path, 'r') as file:
                tree_string = file.read()
            return tree_string
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error reading Newick file: {str(e)}")
            return None

    def read_fasta_file(self, file_path):
        try:
            with open(file_path, 'r') as file:
                sequence = file.read()

            return sequence
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error reading FASTA file: {str(e)}")
            return None

    def node_a_node_b_collector(self):
        all_node_values = []

        # Handle the NodeA and NodeB windows of the main app
        main_node_a = self.text_edit_node_a.toPlainText()
        main_node_b = self.text_edit_node_b.toPlainText()

        all_node_values.append((main_node_a, main_node_b))

        # Handle the NodeA and NodeB windows of the new widgets
        for widget in self.subtree_widgets:
            node_window_a = widget.node_window_a
            node_window_b = widget.node_window_b
            node_a = node_window_a.text_edit.toPlainText()
            node_b = node_window_b.text_edit.toPlainText()

            all_node_values.append((node_a, node_b))

        return all_node_values

    def save_newick_file(self, tree_string, output_directory):
        if output_directory:
            output_file_path = os.path.join(output_directory, "output.nwk")
            try:
                with open(output_file_path, 'w') as file:
                    file.write(tree_string)

                QMessageBox.information(self, "Success", "Newick file saved successfully.")
            except Exception as e:

                QMessageBox.critical(self, "Error", f"Error saving Newick file: {str(e)}")
        else:

            QMessageBox.critical(self, "Error", "No output directory selected.")

    def create_newick_subtrees(self):
        file_path = self.file_window_newick.text_edit.toPlainText()
        tree_file_prefix = os.path.splitext(os.path.basename(file_path))[0]

        if file_path:
            tree_string = self.read_newick_file(file_path)
            if tree_string:
                t = Tree(tree_string, format=1)

                subtrees = []
                current_NodeA = set()
                remaining = t.write(format=1)

                all_node_values = self.node_a_node_b_collector()

                for i, node_value in enumerate(all_node_values, start=1):
                    NodeA, NodeB = node_value
                    current_NodeA.add(NodeA)
                    t, subtree = modify_tree_branch_length(t, NodeA, NodeB)
                    subtrees.append((NodeA, NodeB, subtree))

                    # Update the remaining tree after each modification
                    remaining = t.write(format=1)

                output_directory = self.directory_window_output.text_edit.toPlainText()
                if output_directory:
                    if not os.path.exists(output_directory):
                        os.makedirs(output_directory)

                    for i, (NodeA, NodeB, subtree) in enumerate(subtrees, start=1):
                        filename = os.path.join(output_directory,
                                                f"{tree_file_prefix}_subtree{i}_{NodeA}_{NodeB}.nwk")
                        with open(filename, 'w') as file:
                            file.write(subtree.write(format=1))

                    remaining_newick_filename = os.path.join(output_directory, f"{tree_file_prefix}_remaining.nwk")
                    with open(remaining_newick_filename, 'w') as file:
                        file.write(remaining)

                    QMessageBox.information(self, "Success", "Subtree files saved successfully.")
                else:
                    QMessageBox.critical(self, "Error", "No output directory selected.")
            else:
                QMessageBox.critical(self, "Error", "Error reading Newick file.")
        else:
            QMessageBox.critical(self, "Error", "No Newick file selected.")

    def create_nwk_fasta_subtrees(self):
        fasta_file_path = self.file_window_fasta.text_edit.toPlainText()

        if fasta_file_path:
            fasta_sequence = self.read_fasta_file(fasta_file_path)
            if fasta_sequence:
                file_path = self.file_window_newick.text_edit.toPlainText()
                tree_file_prefix = os.path.splitext(os.path.basename(file_path))[0]

                if file_path:
                    tree_string = self.read_newick_file(file_path)
                    if tree_string:
                        t = Tree(tree_string, format=1)
                        subtrees = []
                        current_NodeA = set()
                        remaining = t.write(format=1)

                        all_node_values = self.node_a_node_b_collector()

                        for i, node_value in enumerate(all_node_values, start=1):
                            NodeA, NodeB = node_value
                            current_NodeA.add(NodeA)
                            t, subtree = modify_tree_branch_length(t, NodeA, NodeB)
                            subtrees.append((NodeA, NodeB, subtree))

                            # Update the remaining tree after each modification
                            remaining = t.write(format=1)

                        output_directory = self.directory_window_output.text_edit.toPlainText()

                        if output_directory:
                            if not os.path.exists(output_directory):
                                os.makedirs(output_directory)

                            subtrees_fasta_directory = os.path.join(output_directory, 'subtrees_fasta')
                            if not os.path.exists(subtrees_fasta_directory):
                                os.makedirs(subtrees_fasta_directory)

                            for i, (NodeA, NodeB, subtree) in enumerate(subtrees, start=1):
                                filename = os.path.join(output_directory,
                                                        f"{tree_file_prefix}_subtree{i}_{NodeA}_{NodeB}.nwk")
                                with open(filename, 'w') as file:
                                    file.write(subtree.write(format=1))

                                subtree = subtree.write(format=1)
                                taxas = subtree_to_fasta(subtree)

                                sequences = fasta_copy(fasta_sequence, taxas)

                                output_prefix = os.path.splitext(os.path.basename(fasta_file_path))[0]
                                output_filename = f"{output_prefix}_subtree{i}_{NodeA}_{NodeB}.fasta"
                                output_file = os.path.join(subtrees_fasta_directory, output_filename)
                                with open(output_file, 'w') as file:
                                    file.write('\n'.join(sequences))

                            taxa_remaining = subtree_to_fasta(remaining)
                            # print(f"taxa_remaining: {taxa_remaining}")
                            sequences_remaining = fasta_copy(fasta_sequence, taxa_remaining)
                            remaining_output_filename = f"{output_prefix}_remaining.fasta"
                            remaining_output_file = os.path.join(subtrees_fasta_directory, remaining_output_filename)
                            with open(remaining_output_file, 'w') as file:
                                file.write('\n'.join(sequences_remaining))

                            remaining_newick_filename = os.path.join(output_directory,
                                                                     f"{tree_file_prefix}_remaining.nwk")
                            with open(remaining_newick_filename, 'w') as file:
                                file.write(remaining)

                            QMessageBox.information(self, "Success", "Subtree files saved successfully.")
                        else:
                            QMessageBox.critical(self, "Error", "No output directory selected.")
                    else:
                        QMessageBox.critical(self, "Error", "Error reading Newick file.")
                else:
                    QMessageBox.critical(self, "Error", "No Newick file selected.")
            else:
                QMessageBox.critical(self, "Error", "Error reading FASTA file.")
        else:
            QMessageBox.critical(self, "Error", "No FASTA file selected.")


app = QApplication([])
icon = QIcon('/home/clemens/PycharmProjects/app/subtree_icon.ico')
app.setWindowIcon(icon)
window = SubtreeCreator()
window.show()
app.exec_()


