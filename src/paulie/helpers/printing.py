from paulie.common.pauli_string import PauliString

def print_vertix(debug, vertix: PauliString, title=""):
    if debug:
        print(f"{title} {vertix}")

def print_vertices(debug, vertices: list[PauliString], title = ""):
    if debug is False:
        return

    print(f"----{title}--lenght = {len(vertices)}")
    for v in vertices:
        print_vertix(debug, v)
    print("-------------------")


def print_lit_vertices(debug, vertices: list[PauliString], lits, title = ""):
    if debug is False:
        return

    print(f"----{title}--lenght = {len(vertices)}")
    for v in vertices:
        title = ""
        if v in lits:
            title = "*"
        print_vertix(debug, v, title)
    print("-------------------")


class Debug:
      def __init__(self, debug):
          self.debug = debug
          self.save_debug = debug

      def get_debug(self):
           return self.debug
     
      def set_debug(self, debug):
          self.debug = debug

      def debuging(self):
          self.debug = True

      def restore(self):
          self.debug = self.save_debug

      def print_vertix(self, vertix: PauliString, title=""):
          print_vertix(self.debug, vertix, title)


      def print_vertices(self, vertices: list[PauliString], title=""):
          print_vertices(self.debug, vertices, title)

      def print_title(self, title):
          if self.debug:
              if title != "":
                  print(f"{title}")

      def print_lit_vertices(self, vertices:list[PauliString], lits:list[PauliString], title = ""):
          print_lit_vertices(self.debug, vertices, lits, title)

      def is_pauli_string(self, vertix: PauliString, paulistring:PauliString):
          return paulistring == vertix
