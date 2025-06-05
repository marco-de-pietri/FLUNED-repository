# extract_names.py
import os
import ast

def collect_names_from_file(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as f:
        try:
            tree = ast.parse(f.read(), filename=path)
        except SyntaxError:
            return []
    names: list[str] = []

    class FuncVisitor(ast.NodeVisitor):
        def __init__(self):
            self.current_class: str | None = None

        def visit_ClassDef(self, node: ast.ClassDef):
            prev = self.current_class
            self.current_class = node.name
            self.generic_visit(node)
            self.current_class = prev

        def visit_FunctionDef(self, node: ast.FunctionDef):
            if self.current_class:
                names.append(f"{self.current_class}.{node.name}")
            else:
                names.append(node.name)
            # don’t recurse into inner functions/methods further
            # but you can call generic_visit if you want nested defs
            # self.generic_visit(node)

    visitor = FuncVisitor()
    visitor.visit(tree)
    return names

def main():
    all_names: list[str] = []
    for root, _, files in os.walk("."):
        for fn in files:
            if fn.endswith(".py"):
                full = os.path.join(root, fn)
                all_names += collect_names_from_file(full)

    # write to output file
    with open("all_function_names.txt", "w", encoding="utf-8") as out:
        for name in sorted(all_names):
            out.write(name + "\n")

if __name__ == "__main__":
    main()

